/* pfp-parse - prefix free parsing parse
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
   \file parse.hpp
   \brief parse.hpp define and build the prefix-free parse data structure.
   \author Massimiliano Rossi
   \date 25/06/2020
   \note This is a short version of the parse in https://github.com/maxrossi91/pfp-data-structures
*/

#ifndef _PFP_PARSE_HH
#define _PFP_PARSE_HH

#include "common.hpp"
#include "disk_vector.hpp"

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
extern "C" {
    #include<gsacak.h>
}

// TODO: Extend it to non-integer alphabets
class parse{
public:
  disk_vector<uint32_t> p;   // sxgc rung 3: out-of-core (scratch copy of .parse + terminator)
  disk_vector<uint_t> saP;   // out-of-core scratch (freed after load)
  disk_vector<uint_t> isaP;  // out-of-core scratch (freed after load)

  disk_vector<int_t> ilist; // Inverted list of phrases of P in BWT_P (kept mapped; raw pointers into it are held by the pfp priority queue — mmap addresses are stable, so this keeps working)
  sdsl::bit_vector ilist_s; // The ith 1 is in correspondence of the first occurrence of the ith phrase
  sdsl::bit_vector::select_1_type select_ilist_s;

  size_t alphabet_size;

  std::string dv_base; // pfp basepath for scratch files

  typedef size_t size_type;

  // Default constructor for load
  parse() {}

  parse(  std::string filename,
          size_t alphabet_size_):
          alphabet_size(alphabet_size_)
  {
    // sxgc rung 3: stream the .parse file into a scratch mapping and append
    // the 0 terminator required by sacak (the original .parse stays intact)
    dv_base = filename;
    std::string tmp_filename = filename + std::string(".parse");
    struct stat st;
    if (stat(tmp_filename.c_str(), &st) != 0)
    { std::cerr << "cannot stat " << tmp_filename << std::endl; exit(1); }
    size_t entries = (size_t)st.st_size / sizeof(uint32_t);
    int fds = ::open(tmp_filename.c_str(), O_RDONLY);
    if (fds < 0) { std::cerr << "cannot open " << tmp_filename << std::endl; exit(1); }
    p.create(dv_path(dv_base, "p"), entries + 1, true);
    if (entries > 0)
      p.copy_from_fd(fds, 0, entries * sizeof(uint32_t), 0);
    ::close(fds);
    p[entries] = 0; // this is the terminator for the sacak algorithm

    compute_freq();

    build();

  }

  void build(){

    // suffix array of the parsing (out-of-core scratch; sacak writes
    // through the mapping)
    saP.create(dv_path(dv_base, "saP"), p.size(), true);
    verbose("Computing SA of the parsing");
    _elapsed_time(
      sacak_int((int_text*)&p[0], &saP[0], p.size(), alphabet_size);
    );



    // inverted list of the parsing.
    verbose("Computing ilist");
    _elapsed_time(
      compute_ilist()
    );


    // inverse suffix array of the parsing (out-of-core scratch).
    verbose("Computing ISA of the parsing");
    _elapsed_time(
      {
        isaP.create(dv_path(dv_base, "isaP"), p.size(), true);
        for(size_t i = 0; i < saP.size(); ++i){
          isaP[saP[i]] = i;
        }
      }
    );


  }


  void compute_ilist()
  {

    ilist.create(dv_path(dv_base, "ilist"), p.size(), true);
    ilist_s = sdsl::bit_vector(p.size() + 1, 0);

    // computing the bucket boundaries
    size_t j = 0;
    ilist_s[j++] = 1; // this is equivalent to set pf.freq[0]=1;
    ilist_s[j] = 1;
    for (size_t i = 1; i < freq.size(); ++i)
    {
      j += freq[i];
      ilist_s[j] = 1;
      freq[i] = 0;
    }

    select_ilist_s = sdsl::bit_vector::select_1_type(&ilist_s);


    for(size_t i = 0; i < saP.size(); ++i)
    {
      size_t prec_phrase_index = (saP[i] == 0 ? p.size() : saP[i]) - 1;
      uint_t prec_phrase = p[prec_phrase_index];

      size_t ilist_p = select_ilist_s(prec_phrase + 1) + freq[prec_phrase]++;
      ilist[ilist_p] = i;
    }

    freq.clear();
    freq.shrink_to_fit();
  }

  void compute_freq()
  {
    uint32_t max_p = p[0];
    for(size_t i = 1; i < p.size(); ++i)
      max_p = std::max(max_p, p[i]);

    freq.resize(max_p+1,0);
    for(size_t i = 0; i < p.size(); ++i)
      freq[p[i]]++;

  }

  // Serialize to a stream.
  // sxgc rung 3: removed — out-of-core members cannot be my_serialize'd and
  // no tool in this tree serializes the PFP structures.

private:
  std::vector<uint_t> freq;
};

#endif /* end of include guard: _PFP_PARSE_HH */
