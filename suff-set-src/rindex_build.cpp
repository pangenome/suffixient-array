// rindex_build.cpp — sxgc r-index toehold (Bit 7).
//
// Consumes the same PFP stream the -A scan uses (BWT of reverse(input) +
// 0-sentinel, in SA order) and emits the run-length BWT:
//   run_char[R], run_len[R], sa_sample[R]  (SA value at each run's LAST position)
// plus the C prefix-count table. Query side needs zero text access.
//
// Output file layout (little-endian):
//   "SXRI" u32, version u32=1
//   n u64 (reversed+sentinel length), sigma u32 (max byte + 1), R u64
//   C[256] u64
//   run_char[R] u8, run_len[R] u32, sa_sample[R] u32

#include <iostream>
#include <vector>
#include <cstring>

#include <common.hpp>
#include <pfp.hpp>
#include <pfp_iterator.hpp>
#include <malloc_count.h>

void help(){
  std::cout << "rindex_build [options]" << std::endl <<
  "Options:" << std::endl <<
  "-i <arg>    Basepath for the PFP data structures." << std::endl <<
  "-o <arg>    Output .ri file." << std::endl <<
  "-w <arg>    PFP trigger string size." << std::endl <<
  "-n <arg>    Text length (incl. sentinel)." << std::endl;
}

int main(int argc, char** argv)
{
  std::string input_path, output_path;
  int w = 10; uint64_t n = 0;
  int opt;
  while((opt = getopt(argc, argv, "i:o:w:n:h")) != -1)
  {
    switch(opt){
      case 'i': input_path = optarg; break;
      case 'o': output_path = optarg; break;
      case 'w': w = atoi(optarg); break;
      case 'n': n = strtoull(optarg, nullptr, 10); break;
      default: help(); return 1;
    }
  }
  if(input_path.empty() || output_path.empty()){ help(); return 1; }

  pf_parsing pf(input_path, w);
  pfp_iterator iter(pf, input_path);

  std::vector<unsigned char> run_char;
  std::vector<uint32_t> run_len;
  std::vector<uint32_t> sa_sample;

  ++iter; // first stream element
  unsigned char b = (unsigned char)iter.get_bwt();
  uint64_t prev_sa = iter.get_sa();
  uint64_t cnt = 1;
  uint64_t scanned = 1;
  while(++iter)
  {
    unsigned char c = (unsigned char)iter.get_bwt();
    uint64_t cur_sa = iter.get_sa();
    if(c != b)
    {
      run_char.push_back(b); run_len.push_back((uint32_t)cnt); sa_sample.push_back((uint32_t)prev_sa);
      b = c; cnt = 1;
    }
    else cnt++;
    prev_sa = cur_sa;
    scanned++;
  }
  run_char.push_back(b); run_len.push_back((uint32_t)cnt); sa_sample.push_back((uint32_t)prev_sa);

  uint64_t R = run_char.size();
  std::cout << "Stream length = " << scanned << std::endl;
  std::cout << "BWT equal-letter runs = " << R << std::endl;

  // C table (count of chars < c) over byte values
  std::vector<uint64_t> occ(256, 0);
  for(size_t r = 0; r < R; ++r) occ[run_char[r]] += run_len[r];
  uint64_t check = 0; for(int c = 0; c < 256; ++c) check += occ[c];
  std::vector<uint64_t> C(256, 0);
  for(int c = 1; c < 256; ++c) C[c] = C[c-1] + occ[c-1];

  FILE* f = fopen(output_path.c_str(), "wb");
  if(!f){ std::cerr << "cannot open " << output_path << std::endl; return 1; }
  uint32_t magic = 0x5258'5349; // "IXSR"-ish tag
  uint32_t version = 1;
  uint32_t sigma = 256;
  fwrite(&magic, 4, 1, f); fwrite(&version, 4, 1, f);
  uint64_t n_out = (n ? n : scanned);
  fwrite(&n_out, 8, 1, f); fwrite(&sigma, 4, 1, f); fwrite(&R, 8, 1, f);
  fwrite(C.data(), 8, 256, f);
  fwrite(run_char.data(), 1, R, f);
  fwrite(run_len.data(), 4, R, f);
  fwrite(sa_sample.data(), 4, R, f);
  fclose(f);

  uint64_t bytes = 8 + 8 + 4 + 8 + 256*8 + R*9;
  std::cout << "r-index written: " << output_path << " (" << R << " runs, "
            << bytes/1048576 << " MiB raw layout)" << std::endl;
  return 0;
}
