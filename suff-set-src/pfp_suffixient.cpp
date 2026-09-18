// Copyright (c) 2024, REGINDEX.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

#include <iostream>
#include <limits>

#include <common.hpp>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>

#include <pfp.hpp>
#include <pfp_iterator.hpp>

#include <malloc_count.h>

struct lcp_maxima
{
  int64_t len;
  uint64_t pos;
  bool active;
};

constexpr int sigma = 128; 
constexpr int SIGMA = 128;

void help(){

  std::cout << "suffixient [options]" << std::endl <<
  "Input: Path to PFP data structures. Output: smallest suffixient set." << std::endl <<
  "Options:" << std::endl <<
  "-h          Print usage info." << std::endl << 
  "-i <arg>    Basepath for the PFP data structures." << std::endl << 
  "-o <arg>    Store output to file using 64-bits unsigned integers. If not specified, output is streamed to standard output in human-readable format." << std::endl <<
  "-w <arg>    PFP trigger string size." << std::endl << 
  "-n <arg>    Text length." << std::endl <<
  "-p          Print to standard output size of suffixient set. Default: false." << std::endl <<
  "-r          Print to standard output number of equal-letter runs in the BWT of reverse text. Default: false." << std::endl;
  exit(0);
}

inline void eval(int64_t l, uint64_t& size, std::vector<lcp_maxima>& R, 
                 std::string output_file, FILE *suffixient_file)
{
  for(uint8_t c = 1; c < sigma; ++c)
    if(l < R[c].len)
    {
      // process an active candidate
      if(R[c].active)
      {
        size++;
        if(output_file.length() == 0){
          std::cout << R[c].pos << " ";
        }
        else
          if (fwrite(&R[c].pos, SSABYTES, 1, suffixient_file) != 1)
            error("S write error 1");
      }
      // update to inactive state
      R[c] = {l,0,false};
    }
}

int main(int argc, char* const argv[])
{
  if(argc<2) help();

  std::string output_file, input_path;

  bool sort=false, chi=false, runs=false, dump=false, emitA=false;

  FILE *suffixient_file;

  int w, N;

  int opt;
  while ((opt = getopt(argc, argv, "Adprsho:w:n:i:")) != -1){
    switch (opt){
      case 'h':
        help();
      break;
      case 'o':
        output_file = std::string(optarg);
      break;
      case 'i':
        input_path = std::string(optarg);
      break;
      case 'w':
        w = atoi(optarg);
      break;
      case 'n':
        N = atoi(optarg);
      break;
      case 'p':
        chi = true;
      break;
      case 'r':
        runs = true;
      break;
      case 'd':
        dump = true;
      break;
      case 'A':
        emitA = true;
      break;
      default:
        help();
      return -1;
    }
  }

  FILE* dump_file = nullptr;
  std::string dump_path = input_path + ".triples";
  // compute PFP data structures
  pf_parsing pf(input_path, w);

  // compute PFP iterator
  pfp_iterator iter(pf, input_path);

  // opening output files
  if(output_file.length() != 0)
  {
    if ((suffixient_file = fopen(output_file.c_str(), "w")) == nullptr)
        error("open() file " + output_file + " failed");
  }
  else
    std::cout << "\nSmallest suffixient set: ";

  /*
  * algorithm: compute suffixient-nexessary set by streaming SA, LCP, and BWT using the PFP data structures.
  */

  if(emitA)
  {
    // ===== Bit 4: sA components scan (mirror one_pass_build_index -t sA) =====
    struct lcp_max4 { int64_t len; uint64_t pos; int64_t lcs; bool active; };
    int64_t max_byte = 1;
    std::vector<lcp_max4> R4(SIGMA, {-1, 0, -1, false});
    std::vector<uint64_t> S0; std::vector<int64_t> L; std::vector<uint64_t> A(SIGMA, 0);
    int64_t m4 = std::numeric_limits<int64_t>::max();
    uint64_t i4 = 1; ++iter;
    char p4 = iter.get_bwt(); uint64_t p4_sa = iter.get_sa();
    max_byte = std::max<int64_t>(max_byte, (int64_t)(unsigned char)p4);
    auto eval4 = [&](int64_t l)
    {
      for(int c = 1; c < (int)SIGMA; ++c)
        if(l < R4[c].len)
        {
          if(R4[c].active)
          {
            S0.push_back(R4[c].pos - 1);          // 0-based (build_index convention)
            L.push_back(R4[c].lcs + 1);
            A[c]++;
          }
          R4[c] = {l, 0, l, false};               // lcs carried = run-min l
        }
    };
    while(++iter)
    {
      m4 = std::min(m4, int64_t(iter.get_lcp()));
      char c4 = iter.get_bwt(); uint64_t c4_sa = iter.get_sa();
      max_byte = std::max<int64_t>(max_byte, (int64_t)(unsigned char)c4);
      if(c4 != p4)
      {
        eval4(m4);
        for(uint64_t ip = i4 - 1; ip < i4 + 1; ++ip)
        {
          char bc = (ip == i4 - 1) ? p4 : c4;                 // BWT(ip): p4 = get_bwt at i4-1 (prev), c4 at i4
          uint64_t sap = (ip == i4 - 1) ? p4_sa : c4_sa;
          (void)sap;
          // NOTE: build_index uses BWT(ip) and SA[ip] of its own stream:
          // at the boundary, BWT(i-1) = p4 with SA(i-1) = p4_sa; BWT(i) = c4 with SA(i) = c4_sa
          auto &Rc = R4[(int)(unsigned char)bc];
          if(int64_t(iter.get_lcp()) > Rc.len)
            Rc = {int64_t(iter.get_lcp()), N - sap, Rc.lcs, true};
        }
        m4 = std::numeric_limits<int64_t>::max();
      }
      p4 = c4; p4_sa = c4_sa; i4++;
    }
    eval4(-1);
    A.resize(max_byte + 1);
    if(output_file.length() != 0)
    {
      std::string b = output_file;
      { std::ofstream o(b + ".suff", std::ios::binary);
        for (const auto& x : S0) { o.write(reinterpret_cast<const char*>(&x), 5); } }
      { std::ofstream o(b + ".lcs", std::ios::binary);
        for (const auto& x : L) { int64_t y = x; o.write(reinterpret_cast<const char*>(&y), 5); } }
      { std::ofstream o(b + ".mult", std::ios::binary);
        for (const auto& x : A) { o.write(reinterpret_cast<const char*>(&x), 5); } }
      std::cout << "Size of smallest suffixient set: " << S0.size() << std::endl;
    }
    if(dump_file) fclose(dump_file);
    return 0;
  }

  // move forward pfp iterator to first position
  uint64_t i=1; ++iter; 
  char p = iter.get_bwt(), c;
  //std::cout << p;
  uint64_t p_sa = iter.get_sa(), c_sa;
  if(dump && output_file.length() != 0)
  {
    if ((dump_file = fopen(dump_path.c_str(), "w")) == nullptr)
      error("open() dump file failed");
    fprintf(dump_file, "%d %lld %llu\n", (int)(unsigned char)p, (long long)iter.get_lcp(), (unsigned long long)p_sa);
  }
  //std::cout << p_sa << " ";
  
  uint64_t bwtruns=1, suffixient_size=0; //tot_size = 1;
  int64_t m = std::numeric_limits<int64_t>::max();
  // vector STORING candidate suffixient right-extensions
  std::vector<lcp_maxima> r_ext(sigma,{-1,0,false});
  
  // iterate until all values have been streamed
  while( ++iter )
  {
    // read current values from the stream
    m = std::min(m,int64_t(iter.get_lcp()));
    c = iter.get_bwt();
    c_sa = iter.get_sa();
    if(dump_file)
      fprintf(dump_file, "%d %lld %llu\n", (int)(unsigned char)c, (long long)iter.get_lcp(), (unsigned long long)c_sa);
    //std::cout << c_sa << " ";
    //tot_size++;

    if(c != p)
    {
      // evaluate sigma candidates
      eval(m,suffixient_size,r_ext,output_file,suffixient_file);
      // update p and c candidates
      if(iter.get_lcp() > r_ext[p].len) 
        r_ext[p] = {iter.get_lcp(),N - p_sa,true};
      if(iter.get_lcp() > r_ext[c].len) 
        r_ext[c] = {iter.get_lcp(),N - c_sa,true};  
      // reset LCP value
      m = std::numeric_limits<int64_t>::max();
      // increment number of runs
      bwtruns++;
    }
    
    // update the previous BWT character and SA entry
    p = c; p_sa = c_sa;
  }
  // evaluate last active candidates
  eval(-1,suffixient_size,r_ext,output_file,suffixient_file);

  if(output_file.length() == 0)
      std::cout << std::endl;
  else
      fclose(suffixient_file);

  if(chi)
    std::cout << "Size of smallest suffixient set: " << suffixient_size << std::endl;
  if(runs)
    std::cout << "Number of equal-letter runs: " << bwtruns << std::endl;
  
  return 0;
}