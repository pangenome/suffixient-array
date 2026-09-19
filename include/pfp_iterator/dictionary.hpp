/* pfp-dictionary - prefix free parsing dictionary
    Copyright (C) 2020 Massimiliano Rossi

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*!
   \file dictionary.hpp
   \brief dictionary.hpp define and build the prefix free parsing dictionary data structure.
   \author Massimiliano Rossi
   \date 25/06/2020
   \note This is a short version of the dictionary in https://github.com/maxrossi91/pfp-data-structures

   sxgc rung 3: the big arrays (d, saD, isaD, lcpD) are out-of-core
   disk_vector mappings (scratch files, unlinked on close) so the structure
   scales to HPRC-v2-sized dictionaries that exceed RAM. The succinct
   structures (b_d, rmq_lcp_D) stay in RAM (~2-3 bits per element).
*/

#ifndef _PFP_DICTIONARY_HH
#define _PFP_DICTIONARY_HH

#include <queue>

#include "common.hpp"
#include "disk_vector.hpp"

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
extern "C" {
    #include<gsacak.h>
}


// TODO: Extend it to integer alphabets
class dictionary{
public:
  disk_vector<uint8_t> d;
  disk_vector<uint_t> saD;
  disk_vector<uint_t> isaD;
  disk_vector<int_t> lcpD;
  sdsl::rmq_succinct_sct<> rmq_lcp_D;
  sdsl::bit_vector b_d; // Starting position of each phrase in D
  sdsl::bit_vector::rank_1_type rank_b_d;
  sdsl::bit_vector::select_1_type select_b_d;

  std::vector<uint8_t> alphabet;

  std::string dv_base; // pfp basepath for scratch files

  typedef size_t size_type;

  // default constructor for load.
  dictionary() {}

  dictionary(std::string filename,
             size_t w)
  {
    // Building dictionary from file (out-of-core: streamed into scratch,
    // prepending w dollars as before)
    dv_base = filename;
    std::string tmp_filename = filename + std::string(".dict");
    int fdd = ::open(tmp_filename.c_str(), O_RDONLY);
    if (fdd < 0) { std::cerr << "cannot open " << tmp_filename << std::endl; exit(1); }
    struct stat st;
    if (fstat(fdd, &st) != 0) { std::cerr << "fstat failed" << std::endl; exit(1); }
    size_t file_bytes = (size_t)st.st_size;

    // Count the leading dollars so we can prepend (w - n_dollars) of them.
    size_t n_dollars = 0;
    {
        std::string head(std::min<size_t>(file_bytes, 1 << 20), '\0');
        ssize_t got = pread(fdd, &head[0], head.size(), 0);
        while (n_dollars < (size_t)got && (uint8_t)head[n_dollars] == Dollar)
            ++n_dollars;
    }

    d.create(dv_path(dv_base, "d"), file_bytes + (w - n_dollars), true);
    // dollars prefix
    memset(d.data(), Dollar, w - n_dollars);
    if (file_bytes > 0)
        d.copy_from_fd(fdd, 0, file_bytes, w - n_dollars);
    ::close(fdd);

    assert(d[0] == Dollar);

    build();

  }

  inline size_t length_of_phrase(size_t id){
    assert(id > 0);
    return select_b_d(id+1)-select_b_d(id) - 1; // to remove the EndOfWord
  }

  inline size_t n_phrases(){
    return rank_b_d(d.size()-1);
  }

  size_t longest_common_phrase_prefix(size_t a, size_t b)
  {
    if(a == 0 || b == 0)
      return 0;
    // Compute the lcp between phrases a and b
    auto a_in_sa = isaD[select_b_d(a)]; // position of the phrase a in saD
    auto b_in_sa = isaD[select_b_d(b)]; // position of the phrase b in saD

    auto lcp_left = std::min(a_in_sa, b_in_sa) + 1;
    auto lcp_right = std::max(a_in_sa, b_in_sa);

    size_t lcp_a_b_i = rmq_lcp_D(lcp_left, lcp_right);
    return lcpD[lcp_a_b_i];

  }

  void build(){

    // Constructing the alphabet
    std::vector<bool> visit(256,false);
    for(auto elem: d)
    {
      if(!visit[elem])
      {
        visit[elem] = true;
        alphabet.push_back(elem);
      }
    }

    // Building the bitvector with a 1 in each starting position of each phrase in D
    b_d.resize(d.size());
    for(size_t i = 0; i < b_d.size(); ++i) b_d[i] = false; // bug in resize
    b_d[0] = true; // Mark the first phrase
    for(size_t i = 1; i < d.size(); ++i )
      b_d[i] = (d[i-1]==EndOfWord);
    b_d[d.size()-1] = true; // This is necessary to get the length of the last phrase

    rank_b_d = sdsl::bit_vector::rank_1_type(&b_d);
    select_b_d = sdsl::bit_vector::select_1_type(&b_d);

    // out-of-core SA/LCP of the dictionary (gsacak writes through the mapping)
    saD.create(dv_path(dv_base, "saD"), d.size(), true);
    lcpD.create(dv_path(dv_base, "lcpD"), d.size(), true);
    verbose("Computing SA, LCP, and DA of dictionary");
    _elapsed_time(
      gsacak(&d[0], &saD[0], &lcpD[0], nullptr, d.size())
    );

    // inverse suffix array of the dictionary (scratch; consumed by
    // compute_s_lcp_T, then unmapped and unlinked by the loader)
    verbose("Computing ISA of dictionary");
    _elapsed_time(
      {
        isaD.create(dv_path(dv_base, "isaD"), d.size(), true);
        for(size_t i = 0; i < saD.size(); ++i){
          isaD[saD[i]] = i;
        }
      }
    );

    verbose("Computing RMQ over LCP of dictionary");
    // Compute the LCP rank of D (sequential read over the mapping)
    _elapsed_time(
      rmq_lcp_D = sdsl::rmq_succinct_sct<>(&lcpD)
    );

  }

};





#endif /* end of include guard: _PFP_DICTIONARY_HH */
