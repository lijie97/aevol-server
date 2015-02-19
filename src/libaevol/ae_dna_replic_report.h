// Aevol - An in silico experimental evolution platform
//
// Copyright: See the AUTHORS file provided with the package or <www.aevol.fr>
// Web: http://www.aevol.fr/
// E-mail: See <http://www.aevol.fr/contact/>
// Original Authors : Guillaume Beslon, Carole Knibbe, David Parsons
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef AEVOL_DNA_REPLIC_REPORT_H
#define AEVOL_DNA_REPLIC_REPORT_H

#include <cinttypes>
#include <list>
#include <cstdlib>
#include <cassert>

#include "ae_mutation.h"

namespace aevol {

class ae_dna;

class ae_dna_replic_report {
  friend class ae_dna;
  
 public :
  // Constructor, copy contructor and destructor are synthesized by compiler.

  // Accessors
  const std::list<ae_mutation> get_mutations()      const { return _mutations; };
  const std::list<ae_mutation> get_rearrangements() const { return _rearrangements; };
  const std::list<ae_mutation> get_HT()             const { return _HT; };

  int32_t get_nb(enum ae_mutation_type t)  const { return _nb_mut[t]; };
  // The following functions should be integrated in get_nb. But enums
  // are not extensible, so there is a semantics challenge to define
  // this properly.
  inline int32_t get_nb_small_mutations()  const;
  inline int32_t get_nb_rearrangements()   const;
  inline int32_t get_nb_HT()               const;
  inline int32_t get_nb_small_insertions() const;
  inline int32_t get_nb_indels()           const;

  // Public Methods
  void compute_stats();  // useful when we inspect a tree file
  inline void add_rear(const ae_mutation& rear);
  inline void add_mut(const ae_mutation& mut);
  inline void add_HT(const ae_mutation& HT);

 protected :
  std::list<ae_mutation> _mutations;       // Lists of mutations, rearrangements and undergone
  std::list<ae_mutation> _rearrangements;  // by the genetic unit at last replication
  std::list<ae_mutation> _HT;
  int32_t _nb_mut[10]; // Number of mutations/rearrangements/HT of each type undergone
};

inline int32_t ae_dna_replic_report::get_nb_small_mutations() const {
  assert(_mutations.size() == static_cast<size_t>(_nb_mut[SWITCH] + _nb_mut[S_INS] + _nb_mut[S_DEL]));
  return _mutations.size();
}

inline int32_t ae_dna_replic_report::get_nb_rearrangements() const {
  assert(_rearrangements.size() ==
         static_cast<size_t>(_nb_mut[DUPL] + _nb_mut[DEL] + _nb_mut[TRANS] + _nb_mut[INV]));
  return _rearrangements.size();
}

inline int32_t ae_dna_replic_report::get_nb_HT() const {
  assert(_HT.size() == static_cast<size_t>(_nb_mut[INS_HT] + _nb_mut[REPL_HT]));
  return _HT.size();
}

inline int32_t ae_dna_replic_report::get_nb_indels() const {
  return _nb_mut[S_INS] + _nb_mut[S_DEL];
}

inline void ae_dna_replic_report::add_mut(const ae_mutation& mut) {
  _mutations.push_back(mut);
  _nb_mut[mut.get_mut_type()]++;
}

inline void ae_dna_replic_report::add_rear(const ae_mutation& rear) {
  _rearrangements.push_back(rear);
  _nb_mut[rear.get_mut_type()]++;
}

inline void ae_dna_replic_report::add_HT(const ae_mutation& HT) {
  _HT.push_back(HT);
  _nb_mut[HT.get_mut_type()]++;
}

} // namespace aevol
#endif // AEVOL_DNA_REPLIC_REPORT_H
