// rindex_build.cpp — sxgc r-index (rung 1): run-length BWT from the streamed PFP.
//
// Consumes the same PFP stream the -A scan uses (BWT of reverse(input) +
// 0-sentinel, in SA order) and emits the run-length BWT:
//   run_char[R], run_len[R], sa_sample[R]  (SA value at each run's LAST position)
// plus the C prefix-count table. Query side needs zero text access.
//
// Output file layout (little-endian):
//   "SXRI" u32, version u32
//     v1: sa_sample[R] u32                       (overflows at n > 4.29 Gbp)
//     v2: sa_sample[R] u64                       (13 B/run)
//     v3: header + C + run_char u8 + run_len u32 raw (FILE-style prefix,
//         written via stream), then sdsl::int_vector sa packed to bits(n)
//         (HPRC v2: 41 bits; 6.125 B/run ~= 15.5 GB projected, R=2.53B)

#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>

#include <sdsl/int_vector.hpp>
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
  std::vector<uint64_t> sa_sample;

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
      run_char.push_back(b); run_len.push_back((uint32_t)cnt); sa_sample.push_back(prev_sa);
      b = c; cnt = 1;
    }
    else cnt++;
    if(cnt > 0xFFFFFFFFULL){ std::cerr << "run length exceeds u32 (n too large?)" << std::endl; return 1; }
    prev_sa = cur_sa;
    scanned++;
  }
  run_char.push_back(b); run_len.push_back((uint32_t)cnt); sa_sample.push_back(prev_sa);

  uint64_t R = run_char.size();
  std::cout << "Stream length = " << scanned << std::endl;
  std::cout << "BWT equal-letter runs = " << R << std::endl;

  // C table (count of chars < c) over byte values
  std::vector<uint64_t> occ(256, 0);
  for(size_t r = 0; r < R; ++r) occ[run_char[r]] += run_len[r];
  uint64_t check = 0; for(int c = 0; c < 256; ++c) check += occ[c];
  std::vector<uint64_t> C(256, 0);
  for(int c = 1; c < 256; ++c) C[c] = C[c-1] + occ[c-1];

  uint64_t n_final = (n ? n : scanned);
  // SA sample width: smallest w with n_final-1 < 2^w (SA values are in [0, n))
  uint8_t sa_w = 1;
  while(sa_w < 64 && (uint64_t(1) << sa_w) <= n_final - 1) ++sa_w;

  std::ofstream f(output_path.c_str(), std::ios::binary);
  if(!f){ std::cerr << "cannot open " << output_path << std::endl; return 1; }
  uint32_t magic = 0x5258'5349; // "IXSR"-ish tag
  uint32_t version = 3;
  uint32_t sigma = 256;
  f.write((char*)&magic, 4);   f.write((char*)&version, 4);
  f.write((char*)&n_final, 8);  f.write((char*)&sigma, 4); f.write((char*)&R, 8);
  f.write((char*)C.data(), 256*8);
  f.write((char*)run_char.data(), R);
  f.write((char*)run_len.data(), 4*R);
  sdsl::int_vector<> sa_iv(R, 0, sa_w);
  for(uint64_t r = 0; r < R; ++r) sa_iv[r] = sa_sample[r];
  sa_iv.serialize(f);
  f.close();

  uint64_t bytes = 8 + 8 + 4 + 8 + 256*8 + R*5 + ((uint64_t)R * sa_w + 7) / 8;
  std::cout << "r-index written: " << output_path << " (" << R << " runs, sa_w="
            << (int)sa_w << ", " << bytes/1048576 << " MiB raw layout)" << std::endl;
  return 0;
}
