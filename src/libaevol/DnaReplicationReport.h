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

#ifndef AEVOL_DNA_REPLICATION_REPORT_H__
#define AEVOL_DNA_REPLICATION_REPORT_H__

#include <list>

#include "Mutation.h"

namespace aevol {

class Dna;

class DnaReplicReport {
  friend class Dna;
  
 public :
  // Constructor, copy contructor and destructor are synthesized by compiler.

  // Accessors
  const std::list<Mutation> get_mutations()      const { return mutations_; };
  const std::list<Mutation> get_rearrangements() const { return rearrangements_; };
  const std::list<Mutation> get_HT()             const { return ht_; };
  int32_t get_nb(MutationType t) const;

  // Public Methods
  void compute_stats();  // useful when we inspect a tree file
  void add_rear(const Mutation & rear);
  void add_rear(Mutation && rear);
  void add_mut(const Mutation & mut);
  void add_mut(Mutation && mut);
  void add_HT(const Mutation & HT);
  void add_HT(Mutation && HT);


 protected :
  std::list<Mutation> mutations_;       // Lists of mutations, rearrangements and undergone
  std::list<Mutation> rearrangements_;  // by the genetic unit at last replication
  std::list<Mutation> ht_;
  int32_t _nb_mut[10]; // Number of mutations/rearrangements/HT of each (simple) type undergone
};

} // namespace aevol
#endif // AEVOL_DNA_REPLICATION_REPORT_H__
