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

#include <list>

#include "ae_mutation.h"

namespace aevol {

class Dna;

class DnaReplicReport {
  friend class Dna;
  
 public :
  // Constructor, copy contructor and destructor are synthesized by compiler.

  // Accessors
  const std::list<ae_mutation> get_mutations()      const { return mutations_; };
  const std::list<ae_mutation> get_rearrangements() const { return rearrangements_; };
  const std::list<ae_mutation> get_HT()             const { return ht_; };
  int32_t get_nb(MutationType t) const;

  // Public Methods
  void compute_stats();  // useful when we inspect a tree file
  void add_rear(const ae_mutation& rear);
  void add_rear(ae_mutation&& rear);
  void add_mut(const ae_mutation& mut);
  void add_mut(ae_mutation&& mut);
  void add_HT(const ae_mutation& HT);
  void add_HT(ae_mutation&& HT);


 protected :
  std::list<ae_mutation> mutations_;       // Lists of mutations, rearrangements and undergone
  std::list<ae_mutation> rearrangements_;  // by the genetic unit at last replication
  std::list<ae_mutation> ht_;
  int32_t _nb_mut[10]; // Number of mutations/rearrangements/HT of each (simple) type undergone
};

} // namespace aevol
#endif // AEVOL_DNA_REPLIC_REPORT_H
