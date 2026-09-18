#ifndef AGC_TEXT_ORACLE_HPP
#define AGC_TEXT_ORACLE_HPP

// sxgc Bit 5: AGC random-access text oracle.
//
// Mirrors suffixient::uncompressed_text_oracle's interface, but text bytes are
// served from an AGC archive via the ragc-ffi cdylib (dlopen — no build-time
// Rust dependency), using the shard sidecar (<text>.names.tsv:
// name \t offset \t len) for flat-offset -> stored-contig-name mapping and
// <text>.agc-ref (one line: absolute archive path; or $SXGC_AGC) for the
// archive. Bytes at contig boundaries are the '$' separators, exactly as the
// flat shard text produced by agc2flat.
//
// Sequence cache: fixed-size windows over the flat shard text (byte-budgeted LRU,
// default 64 KiB windows / 256 MiB) — LCP walks are cache-hits and misses cost a
// partial-range fetch, not a full contig decompression.

#include <common.hpp>

#include <dlfcn.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <list>
#include <string>
#include <unordered_map>
#include <vector>

namespace suffixient {

class agc_text_oracle
{
private:

    // ----- ragc-ffi signatures -----
    typedef void* agc_handle;
    typedef agc_handle (*fn_open)(const char*);
    typedef uint64_t (*fn_len)(agc_handle, const char*);
    typedef int (*fn_range)(agc_handle, const char*, uint64_t, uint64_t, unsigned char*);
    typedef void (*fn_close)(agc_handle);

    fn_open ffi_open = nullptr;
    fn_len ffi_len = nullptr;
    fn_range ffi_range = nullptr;
    fn_close ffi_close = nullptr;
    void* ffi_lib = nullptr;
    agc_handle handle = nullptr;

    // ----- sidecar (flat offsets) -----
    struct Seg { std::string cname; usafe_t offset; usafe_t len; };
    std::vector<Seg> segs;       // sorted by offset
    usafe_t N = 0;

    // ----- window cache (the AGC page cache) -----
    // The flat shard text is treated as one virtual byte array; it is cached in
    // fixed-size windows rather than whole contigs. Miss cost = one partial-range
    // fetch (64 KiB), not a full contig decompression; budget is in bytes so it
    // scales to collections with tens of thousands of contigs.
    usafe_t win_w = usafe_t(1) << 16;                 // 64 KiB windows
    usafe_t win_budget = usafe_t(256) << 20;          // 256 MiB byte budget
    usafe_t win_bytes = 0;
    // stats (env SXGC_ORACLE_STATS=1 prints at exit)
    usafe_t stat_fills = 0, stat_bytes = 0, stat_bytes_at = 0;
    std::list<usafe_t> win_order;                     // LRU: front = oldest
    std::unordered_map<usafe_t,
        std::pair<std::vector<unsigned char>, std::list<usafe_t>::iterator>> win_cache;

    static std::string sidecar_path(std::string textPath)
    {
        if(textPath.length() > 4 && textPath.substr(textPath.length()-4) == ".txt")
            return textPath.substr(0, textPath.length()-4) + ".names.tsv";
        return textPath + ".names.tsv";
    }

    static std::string archive_ref(std::string textPath)
    {
        {
            std::ifstream f(textPath + ".agc-ref");
            std::string l;
            if(f.is_open() && std::getline(f, l)) { if(!l.empty() && l.back()=='\r') l.pop_back(); if(!l.empty()) return l; }
        }
        {
            const char* e = std::getenv("SXGC_AGC");
            if(e) return std::string(e);
        }
        std::cerr << "agc_text_oracle: no archive path — write <text>.agc-ref (one line: archive path) or set $SXGC_AGC" << std::endl;
        exit(1);
    }

    void load_sidecar(std::string textPath)
    {
        std::ifstream f(sidecar_path(textPath));
        if(!f.is_open())
        {
            std::cerr << "agc_text_oracle: cannot open sidecar " << sidecar_path(textPath) << std::endl;
            exit(1);
        }
        std::string line;
        while(std::getline(f, line))
        {
            if(line.empty()) continue;
            if(!line.empty() && line.back()=='\r') line.pop_back();
            auto t1 = line.find('\t'); if(t1 == std::string::npos) continue;
            auto t2 = line.find('\t', t1+1); if(t2 == std::string::npos) continue;
            Seg s;
            s.cname = line.substr(0, t1);
            s.offset = std::stoull(line.substr(t1+1, t2-t1-1));
            s.len = std::stoull(line.substr(t2+1));
            segs.push_back(s);
        }
        std::sort(segs.begin(), segs.end(), [](const Seg& a, const Seg& b){ return a.offset < b.offset; });
        N = 0;
        for(auto& s : segs) N = s.offset + s.len + 1;   // incl. trailing '$'
    }

    void ffi_load(std::string archivePath)
    {
        ffi_lib = dlopen("libragc_ffi.so", RTLD_NOW);
        if(!ffi_lib)
        {
            std::cerr << "agc_text_oracle: cannot dlopen libragc_ffi.so: " << dlerror()
                      << " — build sxgc/ragc-ffi and set LD_LIBRARY_PATH" << std::endl;
            exit(1);
        }
        ffi_open = (fn_open)dlsym(ffi_lib, "sxgc_agc_open");
        ffi_len = (fn_len)dlsym(ffi_lib, "sxgc_agc_len");
        ffi_range = (fn_range)dlsym(ffi_lib, "sxgc_agc_range");
        ffi_close = (fn_close)dlsym(ffi_lib, "sxgc_agc_close");
        if(!ffi_open || !ffi_len || !ffi_range || !ffi_close)
        {
            std::cerr << "agc_text_oracle: missing sxgc_agc_* symbols" << std::endl;
            exit(1);
        }
        handle = ffi_open(archivePath.c_str());
        if(!handle)
        {
            std::cerr << "agc_text_oracle: sxgc_agc_open failed for " << archivePath << std::endl;
            exit(1);
        }
    }

    // stored contig containing flat offset i (last seg with offset <= i)
    usafe_t shard_of(usafe_t i)
    {
        usafe_t lo = 0, hi = segs.size();
        while(lo + 1 < hi)
        {
            usafe_t mid = (lo + hi) / 2;
            if(segs[mid].offset <= i) lo = mid; else hi = mid;
        }
        return lo;
    }

    void fill_window(usafe_t wid)
    {
        usafe_t a = wid * win_w, b = a + win_w;
        std::vector<unsigned char> buf(win_w, 0);
        // shards intersecting [a, b); layout: shard bytes [s0, s0+len), '$' at s0+len
        for(usafe_t z = shard_of(a); z < segs.size() && segs[z].offset < b; ++z)
        {
            Seg& s = segs[z];
            usafe_t s0 = s.offset, s1 = s0 + s.len;
            usafe_t x0 = std::max(a, s0), x1 = std::min(b, s1);
            if(x1 > x0)
            {
                int got = ffi_range(handle, s.cname.c_str(), x0 - s0, x1 - s0,
                                    buf.data() + (x0 - a));
                if(got < 0 || (usafe_t)got != x1 - x0)
                {
                    std::cerr << "agc_text_oracle: read failed for " << s.cname << std::endl;
                    exit(1);
                }
                stat_bytes += (x1 - x0);
            }
            if(s1 >= a && s1 < b) buf[s1 - a] = '$';
        }
        win_cache[wid] = {std::move(buf), win_order.insert(win_order.end(), wid)};
        win_bytes += win_w;
        stat_fills++;
        while(win_bytes > win_budget)
        {
            usafe_t old = win_order.front();
            win_bytes -= win_w;
            win_order.pop_front();
            win_cache.erase(old);
        }
    }

    // shard byte at flat offset i ('$' at boundaries)
    unsigned char byte_at(usafe_t i)
    {
        stat_bytes_at++;
        if(i >= N) return '$';
        usafe_t wid = i / win_w;
        auto it = win_cache.find(wid);
        if(it == win_cache.end())
        {
            fill_window(wid);
            it = win_cache.find(wid);
        }
        else
        {
            // O(1) LRU touch
            win_order.splice(win_order.end(), win_order, it->second.second);
        }
        return it->second.first[i - wid * win_w];
    }

public:

    agc_text_oracle(){}

    ~agc_text_oracle()
    {
        if(std::getenv("SXGC_ORACLE_STATS"))
            std::cerr << "agc_oracle_stats: fills=" << stat_fills
                      << " ffi_bytes=" << stat_bytes
                      << " byte_at=" << stat_bytes_at << std::endl;
    }

    void build(std::string input_file_path)
    {
        load_sidecar(input_file_path);
        ffi_load(archive_ref(input_file_path));
    }

    // nothing to persist beyond the .agc-ref link (written here for symmetry)
    usafe_t text_length(){ return this->N; }

    usafe_t store(std::string output_file_path)
    {
        std::string a = archive_ref(output_file_path);
        std::ofstream f(output_file_path + ".agc-ref");
        f << a << std::endl;
        return N;
    }

    void load(std::string input_file_path)
    {
        build(input_file_path);
    }

    usafe_t size(){ return N; }

    unsigned char extract(usafe_t i){ return byte_at(i); }

    std::string display(usafe_t i, usafe_t j)
    {
        std::string res(j - i + 1, 0);
        for(usafe_t y = i; y <= j; ++y) res[y - i] = byte_at(y);
        return res;
    }

    // backward matcher used by the binary search (mirrors uncompressed_text_oracle)
    std::pair<usafe_t,char_t> LCS_char(std::string& pattern, usafe_t p, usafe_t t)
    {
        usafe_t matched_chars = 0;
        usafe_t available_chars = std::min(p+1,t+1);
        while(available_chars > 0)
        {
            if(pattern[p-matched_chars] != byte_at(t-matched_chars))
                return std::make_pair(matched_chars, byte_at(t-matched_chars));
            matched_chars++;
            available_chars--;
        }
        return std::make_pair(matched_chars, static_cast<char_t>(-1));
    }

    usafe_t LCP(std::string& pattern, usafe_t p, usafe_t t)
    {
        usafe_t matched_chars = 0;
        usafe_t available_chars = std::min((pattern.size() - p), (this->N - t));
        while(available_chars > 0)
        {
            if(pattern[p + matched_chars] != byte_at(t + matched_chars))
                return matched_chars;
            matched_chars++;
            available_chars--;
        }
        return matched_chars;
    }
};

}

#endif
