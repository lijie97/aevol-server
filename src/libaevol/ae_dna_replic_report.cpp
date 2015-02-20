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

#include <cinttypes>
#include <list>
#include <cstdlib>
#include <cassert>

#include "ae_dna_replic_report.h"
#include "ae_mutation.h"

namespace aevol {

int32_t ae_dna_replic_report::get_nb(MutationType t)  const {
  switch (t) {
    case S_MUT:
      assert(_mutations.size() ==
             static_cast<size_t>(_nb_mut[SWITCH] +
                                 _nb_mut[S_INS] +
                                 _nb_mut[S_DEL]));
      return _mutations.size();
    case REARR:
      assert(_rearrangements.size() ==
             static_cast<size_t>(_nb_mut[DUPL] +
                                 _nb_mut[DEL] +
                                 _nb_mut[TRANS] +
                                 _nb_mut[INV]));
      return _rearrangements.size();
    case HT:
      assert(_HT.size() ==
             static_cast<size_t>(_nb_mut[INS_HT] +
                                 _nb_mut[REPL_HT]));
      return _HT.size();
    case INDEL:
      return _nb_mut[S_INS] + _nb_mut[S_DEL];
    default: // Simple mutation type.
      return _nb_mut[t];
  };
}

void ae_dna_replic_report::add_mut(const ae_mutation& mut) {
  _mutations.push_back(mut);
  _nb_mut[mut.get_mut_type()]++;
}

void ae_dna_replic_report::add_rear(const ae_mutation& rear) {
  _rearrangements.push_back(rear);
  _nb_mut[rear.get_mut_type()]++;
}

void ae_dna_replic_report::add_HT(const ae_mutation& HT) {
  _HT.push_back(HT);
  _nb_mut[HT.get_mut_type()]++;
}


/// Useful when we inspect a tree file
/// because stats are not saved in the file.
void ae_dna_replic_report::compute_stats( void )
{
  _nb_mut[SWITCH] = 0;
  _nb_mut[S_INS]  = 0;
  _nb_mut[S_DEL]  = 0;
  _nb_mut[DUPL]   = 0;
  _nb_mut[DEL]    = 0;
  _nb_mut[TRANS]  = 0;
  _nb_mut[INV]    = 0;
  _nb_mut[INS_HT] = 0;
  _nb_mut[REPL_HT]= 0;
  
  for (const auto& event: _HT) {
    assert(event.get_mut_type() == INS_HT or
           event.get_mut_type() == REPL_HT);
    _nb_mut[event.get_mut_type()]++;
  }
  
  for (const auto& event: _rearrangements) {
    assert(event.get_mut_type() == DUPL or
           event.get_mut_type() == DEL or
           event.get_mut_type() == TRANS or
           event.get_mut_type() == INV);
    _nb_mut[event.get_mut_type()]++;
  }
  
  for (const auto& event: _mutations) {
    assert(event.get_mut_type() == SWITCH or
           event.get_mut_type() == S_INS or
           event.get_mut_type() == S_DEL);
    _nb_mut[event.get_mut_type()]++;
  }
}
} // namespace aevol
