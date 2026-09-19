/* pfp - prefix free parsing 
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
   \file pfp.hpp
   \brief pfp.hpp define and build the prefix-free parsing data structures.
   \author Massimiliano Rossi
   \date 25/06/2020
   \note This is a short version of the prefix free parsing in https://github.com/maxrossi91/pfp-data-structures
*/

#ifndef _PFP_HH
#define _PFP_HH

#include "common.hpp"

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/io.hpp>
extern "C" {
    #include<gsacak.h>
}

#include "dictionary.hpp"
#include "parse.hpp"

class pf_parsing{
public:
  
  dictionary dict;
  parse pars;
  // std::vector<uint_t> freq;
  size_t n; // Size of the text
  size_t w; // Size of the window

  disk_vector<size_t> s_lcp_T; // sxgc rung 3: out-of-core — LCP array of T sampled at the beginning of each phrase
  sdsl::rmq_succinct_sct<> rmq_s_lcp_T;
  
  disk_vector<size_t> pos_T; // sxgc rung 3: out-of-core — for each suffix of P we store the starting position of that suffix in T.
  // std::vector<int_t>  ilist;            // Inverted list of phrases of P in BWT_P
  // sdsl::bit_vector ilist_s; // The ith 1 is in correspondence of the first occurrence of the ith phrase
  // sdsl::bit_vector::select_1_type select_ilist_s;

  // sdsl::bit_vector b_p;
  // sdsl::bit_vector::rank_1_type rank_b_p;
  // sdsl::bit_vector::select_1_type select_b_p;

  typedef size_t size_type;

  std::string dv_base; // pfp basepath for scratch files

  // Default constructor for load
  pf_parsing() {}

  // sxgc rung 3: the RAM-vector constructor is gone — the arrays are out of
  // core and only the file-based path exists in this tree.

  pf_parsing( std::string filename, size_t w_):
              dict(filename, w_),
              pars(filename,dict.n_phrases()+1),
              // freq(),
              // ilist(pars.p.size()),
              // ilist_s(pars.p.size()+ 1, 0),
              w(w_)
  {
    dv_base = filename;

    // Compute the length of the string;
    compute_n();

    verbose("Computing pos_T");
    _elapsed_time(compute_pos_T());

    verbose("Computing s_lcp_T");
    _elapsed_time(compute_s_lcp_T());

    print_sizes();

    print_stats();

    // Clear unnecessary elements
    clear_unnecessary_elements();
  }

  void print_sizes()
  {

    verbose("Parse");
    verbose("Size of pars.p: ", pars.p.size() * sizeof(pars.p[0]));
    verbose("Size of pars.saP: ", pars.saP.size() * sizeof(pars.saP[0]));
    verbose("Size of pars.isaP: ", pars.isaP.size() * sizeof(pars.isaP[0]));

    verbose("Size of pars.ilist: ", pars.ilist.size() * sizeof(pars.ilist[0]));
    verbose("Size of pars.ilist_s: ", sdsl::size_in_bytes(pars.ilist_s));
    verbose("Size of pars.select_ilist_s: ", sdsl::size_in_bytes(pars.select_ilist_s));

    verbose("Dictionary");
    verbose("Size of dict.d: ", dict.d.size() * sizeof(dict.d[0]));
    verbose("Size of dict.saD: ", dict.saD.size() * sizeof(dict.saD[0]));
    verbose("Size of dict.isaD: ", dict.isaD.size() * sizeof(dict.isaD[0]));
    verbose("Size of dict.lcpD: ", dict.lcpD.size() * sizeof(dict.lcpD[0]));
    verbose("Size of dict.rmq_lcp_D: ", sdsl::size_in_bytes(dict.rmq_lcp_D));

    verbose("Size of dict.b_d: ", sdsl::size_in_bytes(dict.b_d));
    verbose("Size of dict.rank_b_d: ", sdsl::size_in_bytes(dict.rank_b_d));
    verbose("Size of dict.select_b_d: ", sdsl::size_in_bytes(dict.select_b_d));

    verbose("PFP");
    verbose("Size of s_lcp_T: ", s_lcp_T.size() * sizeof(s_lcp_T[0]));
    verbose("Size of rmq_s_lcp_T: ", sdsl::size_in_bytes(rmq_s_lcp_T));
    verbose("Size of pos_T: ", pos_T.size() * sizeof(pos_T[0]));

  }

  void print_stats()
  {

    verbose("Text length: ", n);
    verbose("Parse length: ", pars.p.size());
    verbose("Dictionary length: ", dict.d.size());
    verbose("Number of phrases: ", dict.n_phrases());

  }

  void compute_n(){
    // Compute the length of the string;
    n = 0;
    for (size_t j = 0; j < pars.p.size() - 1; ++j)
    {
      // parse.p[j]: phrase_id
      assert(pars.p[j] != 0);
      n += dict.length_of_phrase(pars.p[j]) - w;
    }
    //n += w; // + w because n is the length including the last w markers
    //n += w - 1; // Changed after changind b_d in dict // -1 is for the first dollar + w because n is the length including the last w markers
  }

  void compute_pos_T(){
    // sxgc rung 3: out-of-core — poss is a scratch mapping (sequential write,
    // random single reads), unlinked at scope exit; pos_T stays mapped for
    // the iterator pass.
    disk_vector<size_t> poss;
    poss.create(dv_path(dv_base, "poss"), pars.p.size(), true);
    for (size_t j = 1; j < pars.p.size(); ++j)
    {
      poss[j] = poss[j-1] + dict.length_of_phrase(pars.p[j-1]) - w;
    }

    pos_T.create(dv_path(dv_base, "pos_T"), pars.p.size(), true);
    for (size_t j = 0; j < pars.p.size(); ++j)
    {
      size_t pos = poss[pars.saP[j]];
      if(pos == 0)
        pos = n; // - w;
      pos_T[j] = pos;
    }

  }

  // Return the frequency of the phrase
  size_t get_freq(size_t phrase)
  {
    return pars.select_ilist_s(phrase + 2) - pars.select_ilist_s(phrase + 1);
    // return select_ilist_s(phrase + 2) - select_ilist_s(phrase + 1);
  }

  // Customized Kasai et al.
  void compute_s_lcp_T()
  {
    size_t n = pars.saP.size();
    // sxgc rung 3: out-of-core scratch mapping (random single writes through
    // the mapping; writeback is sequential). Stays mapped for the iterator
    // pass (rmq_s_lcp_T + point reads).
    s_lcp_T.create(dv_path(dv_base, "s_lcp_T"), n, true);

    size_t l = 0;
    size_t lt = 0;
    for (size_t i = 0; i < n; ++i)
    {
      // if i is the last character LCP is not defined
      size_t k = pars.isaP[i];
      if (k > 0)
      {
        size_t j = pars.saP[k - 1];
        // I find the longest common prefix of the i-th suffix and the j-th suffix.
        while (pars.p[i + l] == pars.p[j + l])
        {
          lt += dict.length_of_phrase(pars.p[i + l]) - w; // I remove the last w overlapping characters
          l++;
        }
        size_t lcpp = dict.longest_common_phrase_prefix(pars.p[i + l], pars.p[j + l]);

        // l stores the length of the longest common prefix between the i-th suffix and the j-th suffix
        s_lcp_T[k] = lt + lcpp;
        if (l > 0)
        {
          l--;
          lt -= dict.length_of_phrase(pars.p[i]) - w; // I have to remove the length of the first matching phrase
        }

      }
    }

    rmq_s_lcp_T = sdsl::rmq_succinct_sct<>(&s_lcp_T);
  }

  void clear_unnecessary_elements(){
    // sxgc rung 3: out-of-core — unmap and unlink the scratch arrays that are
    // no longer needed. This frees mapping space AND disk (at HPRC v2 scale
    // these are ~1 TB of scratch).
    pars.isaP.close();   // sacak ISA of the parse — done
    pars.saP.close();    // sacak SA of the parse — done
    pars.p.close();     // scratch copy of .parse (the .parse file itself stays)

    // isaD is only consumed by compute_s_lcp_T (via
    // longest_common_phrase_prefix), which has already run. The pfp_iterator
    // never touches it.
    dict.isaD.close();
    
  }

  // sxgc rung 3: serialize/load removed — out-of-core members; no tool in
  // this tree serializes the PFP structures.

  std::string filesuffix() const
  {
    return ".pf.ds.thr";
  }


};

#endif /* end of include guard: _PFP_HH */
