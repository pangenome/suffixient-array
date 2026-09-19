// disk_vector.hpp — sxgc rung 3: out-of-core mmap-backed arrays.
//
// Minimal file-backed vector for the PFP machinery at HPRC-v2 scale
// (n = 1.4 Tbp): the dictionary SA/LCP, the parse SA, pos_T, s_lcp_T and
// ilist together exceed available RAM, so they live in scratch files on the
// work filesystem and are mmap'd. Element addresses are stable for the
// mapping's lifetime, so code taking raw pointers into these arrays (e.g.
// the pfp priority queue over ilist) keeps working unchanged; pages flow
// through the page cache, so resident RAM stays bounded by the kernel.
//
// Scratch files are unlinked on close (anonymous-then-deleted semantics);
// persistent arrays simply stay mapped for the process lifetime.

#ifndef _SXGC_DISK_VECTOR_HH
#define _SXGC_DISK_VECTOR_HH

#include <cerrno>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

template <typename T>
class disk_vector {
public:
    typedef T value_type;
    typedef T* iterator;
    typedef const T* const_iterator;
    typedef size_t size_type;

    disk_vector() {}
    ~disk_vector() { close(); }

    disk_vector(const disk_vector&) = delete;
    disk_vector& operator=(const disk_vector&) = delete;

    // Create a fresh sparse writable backing file of n elements.
    void create(const std::string& path, size_t n, bool scratch = true)
    {
        close();
        path_ = path; scratch_ = scratch; n_ = n;
        if (n_ == 0) n_ = 1; // keep a valid (never dereferenced) mapping
        fd_ = ::open(path.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0644);
        if (fd_ < 0) fail("create " + path);
        if (ftruncate(fd_, (off_t)n_ * sizeof(T)) != 0) fail("ftruncate " + path);
        map_ = (T*)mmap(nullptr, n_ * sizeof(T), PROT_READ | PROT_WRITE, MAP_SHARED, fd_, 0);
        if (map_ == MAP_FAILED) fail("mmap " + path);
    }

    // Attach an existing file; size from stat when n == 0. Writable so callers
    // can patch in-place (e.g. sacak terminators).
    void attach(const std::string& path, size_t n = 0, bool scratch = false)
    {
        close();
        path_ = path; scratch_ = scratch;
        fd_ = ::open(path.c_str(), O_RDWR);
        if (fd_ < 0) fail("open " + path);
        struct stat st;
        if (fstat(fd_, &st) != 0) fail("fstat " + path);
        n_ = n ? n : (size_t)st.st_size / sizeof(T);
        if (n_ > 0)
        {
            map_ = (T*)mmap(nullptr, n_ * sizeof(T), PROT_READ | PROT_WRITE, MAP_SHARED, fd_, 0);
            if (map_ == MAP_FAILED) fail("mmap " + path);
        }
    }

    // Stream-copy n elements from open fd src at byte offset src_off into
    // this mapping at element offset dst_off (chunked; never pins the whole
    // source in anonymous RAM).
    void copy_from_fd(int src_fd, off_t src_off, size_t n, size_t dst_off)
    {
        const size_t CHUNK = 16 << 20; // bytes
        std::string buf(CHUNK, '\0');
        size_t done = 0;
        while (done < n)
        {
            size_t want = std::min(CHUNK, n - done);
            ssize_t got = pread(src_fd, &buf[0], want, src_off + (off_t)done);
            if (got <= 0) fail("pread (short read)");
            memcpy((char*)(map_ + dst_off) + done, &buf[0], (size_t)got);
            done += (size_t)got;
        }
    }

    T& operator[](size_t i) { return map_[i]; }
    const T& operator[](size_t i) const { return map_[i]; }
    T* data() { return map_; }
    const T* data() const { return map_; }
    T* begin() { return map_; }
    T* end() { return map_ + n_; }
    size_t size() const { return n_; }
    bool empty() const { return n_ == 0; }

    // Unmap; unlink the backing file when scratch.
    void close()
    {
        if (map_) { munmap(map_, n_ * sizeof(T)); map_ = nullptr; }
        if (fd_ >= 0)
        {
            ::close(fd_); fd_ = -1;
            if (scratch_) ::unlink(path_.c_str());
        }
        n_ = 0;
    }

private:
    void fail(const std::string& what) const
    {
        throw std::runtime_error("disk_vector: " + what + ": " +
                                  std::string(strerror(errno)));
    }
    int fd_ = -1;
    T* map_ = nullptr;
    size_t n_ = 0;
    bool scratch_ = false;
    std::string path_;
};

// Scratch file naming: <pfp-base>.dv.<what>.<pid> — pid suffix so concurrent
// builds over the same base never collide; stale files after a crash are
// reclaimed by the work-dir sweep.
inline std::string dv_path(const std::string& base, const char* what)
{
    return base + ".dv." + what + "." + std::to_string((long)getpid());
}

#endif /* _SXGC_DISK_VECTOR_HH */
