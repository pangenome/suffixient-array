// rindex_query.cpp — sxgc r-index toehold (Bit 7): backward search + locate.
//
// Index = run-length BWT (of reverse(input) + 0-sentinel) built by rindex_build.
// Query = standard backward search over the RLBWT. Patterns are given in
// FORWARD orientation; the tool reverses them internally (the index is over the
// reversed text, so the backward search reads P[0], P[1], ... in order).
// Occurrences are reported in the same ">name\npos len" format as locate, with
// the position mapping calibrated empirically (see SXGC_RI_OFFSET).
//
// Zero text access: every operation is index arithmetic.

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cstring>
#include <cstdlib>
#include <cstdint>
#include <fstream>
#include <sdsl/int_vector.hpp>

struct RIndex
{
    uint64_t n = 0;          // reversed+sentinel length
    uint64_t R = 0;          // number of runs
    std::vector<uint64_t> C;         // 256 prefix counts
    std::vector<unsigned char> run_char;
    std::vector<uint32_t> run_len;
    std::vector<uint64_t> sa_sample;   // v2: 64-bit

    // derived
    std::vector<uint64_t> run_start_blk;   // every BLK runs: cumulative start
    static const uint64_t BLK = 64;
    // per-char: sorted run indices and prefix sums of lengths
    std::vector<std::vector<uint32_t>> cruns;
    std::vector<std::vector<uint64_t>> csum;
    std::vector<uint64_t> total;           // count of each char

    void load(const std::string& path)
    {
        std::ifstream f(path, std::ios::binary);
        if(!f){ std::cerr << "cannot open " << path << std::endl; exit(1); }
        uint32_t magic, version, sigma;
        f.read((char*)&magic, 4);
        if(!f || magic != 0x5258'5349){ std::cerr << "bad magic\n"; exit(1); }
        f.read((char*)&version, 4); f.read((char*)&n, 8); f.read((char*)&sigma, 4); f.read((char*)&R, 8);
        C.resize(256); f.read((char*)C.data(), 256*8);
        run_char.resize(R); f.read((char*)run_char.data(), R);
        run_len.resize(R); f.read((char*)run_len.data(), 4*R);
        sa_sample.resize(R);
        if(version == 3)
        {
            // packed SA: sdsl int_vector of width bits(n) (its header carries size+width)
            sdsl::int_vector<> sa_iv;
            sa_iv.load(f);
            for(uint64_t r = 0; r < R; ++r) sa_sample[r] = sa_iv[r];
        }
        else if(version == 2)
        {
            f.read((char*)sa_sample.data(), 8*R);
        }
        else
        {
            std::vector<uint32_t> tmp(R); f.read((char*)tmp.data(), 4*R);
            for(uint64_t r = 0; r < R; ++r) sa_sample[r] = tmp[r];
        }
        f.close();

        // block index of run starts
        run_start_blk.reserve(R/BLK + 2);
        uint64_t acc = 0;
        for(uint64_t r = 0; r < R; ++r)
        {
            if(r % BLK == 0) run_start_blk.push_back(acc);
            acc += run_len[r];
        }
        run_start_blk.push_back(acc);   // n

        // per-char run index + length prefix sums
        cruns.resize(256); csum.resize(256); total.assign(256, 0);
        std::vector<uint64_t> acc256(256, 0);
        for(uint64_t r = 0; r < R; ++r)
        {
            unsigned char c = run_char[r];
            cruns[c].push_back((uint32_t)r);
            csum[c].push_back(acc256[c]);
            acc256[c] += run_len[r];
        }
        for(int c = 0; c < 256; ++c) total[c] = acc256[c];
        // sentinel: csum[c][size] = total count of c (rank past the last c-run)
        for(int c = 0; c < 256; ++c) csum[c].push_back(total[c]);
    }

    uint64_t run_start(uint64_t r) const
    {
        uint64_t s = run_start_blk[r / BLK];
        for(uint64_t k = (r / BLK) * BLK; k < r; ++k) s += run_len[k];
        return s;
    }

    // run containing BWT position i (i in [0, n))
    uint64_t run_of(uint64_t i) const
    {
        // binary search over block bases, then linear
        uint64_t lo = 0, hi = run_start_blk.size() - 1;
        while(lo + 1 < hi)
        {
            uint64_t mid = (lo + hi) / 2;
            if(run_start_blk[mid] <= i) lo = mid; else hi = mid;
        }
        uint64_t r = lo * BLK;
        uint64_t s = run_start_blk[lo];
        while(r + 1 < R && s + run_len[r] <= i){ s += run_len[r]; ++r; }
        return r;
    }

    // count of char c in BWT[0..i)   (i may be n)
    uint64_t rank(unsigned char c, uint64_t i) const
    {
        if(i == 0) return 0;
        if(i >= n) return total[c];
        uint64_t r = run_of(i);
        uint64_t s = run_start(r);
        // count of c before run r
        uint64_t before = 0;
        if(run_char[r] == c)
        {
            // index of r among c-runs
            const auto& v = cruns[c];
            uint64_t j = std::lower_bound(v.begin(), v.end(), (uint32_t)r) - v.begin();
            before = csum[c][j] + (i - s);
        }
        else
        {
            const auto& v = cruns[c];
            uint64_t j = std::lower_bound(v.begin(), v.end(), (uint32_t)r) - v.begin();
            before = csum[c][j];
        }
        return before;
    }

    // LF mapping
    uint64_t lf(uint64_t i) const
    {
        uint64_t r = run_of(i);
        unsigned char c = run_char[r];
        return C[c] + rank(c, i);
    }

    // backward search; returns half-open [l, r) of suffixes prefixed by pattern
    std::pair<uint64_t,uint64_t> search(const std::string& p, bool trace = false) const
    {
        uint64_t l = 0, r = n;
        for(size_t k = 0; k < p.size() && l < r; ++k)
        {
            unsigned char c = (unsigned char)p[k];
            uint64_t rl = rank(c, l), rr = rank(c, r);
            if(trace)
                std::cerr << "step " << k << " c=" << (int)c << " [l,r)=[" << l << "," << r
                          << ") rank(c,l)=" << rl << " rank(c,r)=" << rr
                          << " -> [" << C[c]+rl << "," << C[c]+rr << ")" << std::endl;
            l = C[c] + rl;
            r = C[c] + rr;
        }
        return {l, r};
    }

    // SA value at BWT position j (walk LF to nearest run end)
    uint64_t sa_at(uint64_t j) const
    {
        uint64_t steps = 0;
        uint64_t pos = j;
        while(true)
        {
            uint64_t r = run_of(pos);
            uint64_t e = run_start(r) + run_len[r];   // exclusive run end
            if(pos == e - 1) return (sa_sample[r] + steps) % n;
            pos = lf(pos);
            steps++;
        }
    }
};

int main(int argc, char** argv)
{
    if(argc < 3)
    {
        std::cerr << "usage: " << argv[0] << " <index.ri> <patterns.fa> [-o out.occs] [-N flatlen]\n";
        return 1;
    }
    std::string idx_path = argv[1], pat_path = argv[2];
    std::string out_path;
    uint64_t flatlen = 0;
    for(int i = 3; i < argc; ++i)
    {
        if(!strcmp(argv[i], "-o") && i+1 < argc) out_path = argv[++i];
        else if(!strcmp(argv[i], "-N") && i+1 < argc) flatlen = strtoull(argv[++i], nullptr, 10);
    }

    RIndex ri;
    ri.load(idx_path);
    std::cerr << "loaded: n=" << ri.n << " R=" << ri.R << std::endl;

    std::ofstream of;
    std::ostream* out = &std::cout;
    if(!out_path.empty()){ of.open(out_path); out = &of; }

    // read patterns
    std::vector<std::pair<std::string,std::string>> pats; // name, seq
    {
        std::ifstream pf(pat_path);
        std::string line, name, seq;
        while(std::getline(pf, line))
        {
            if(line.empty()) continue;
            if(!line.empty() && line.back()=='\r') line.pop_back();
            if(line[0] == '>'){ if(!name.empty()) pats.push_back({name, seq}); name = line.substr(1); seq.clear(); }
            else seq += line;
        }
        if(!name.empty()) pats.push_back({name, seq});
    }

    // optional empirical position offset calibration:
    // forward_start = sa - SXGC_RI_OFFSET   (or plain sa if env unset)
    long long ri_offset = -1;
    { const char* e = std::getenv("SXGC_RI_OFFSET"); if(e) ri_offset = atoll(e); }

    for(auto& pr : pats)
    {
        static const bool trace = std::getenv("SXGC_RI_TRACE");
        auto [l, r] = ri.search(pr.second, trace);
        (*out) << ">" << pr.first << std::endl;
        if(l >= r){ (*out) << "-1 " << pr.second.size() << std::endl; continue; }
        for(uint64_t j = l; j < r; ++j)
        {
            uint64_t sa = ri.sa_at(j);
            long long pos = (ri_offset >= 0) ? (long long)sa - ri_offset
                                             : (long long)(flatlen - pr.second.size() - sa);
            (*out) << pos << " " << pr.second.size() << std::endl;
        }
    }
    if(!out_path.empty()) of.close();
    return 0;
}
