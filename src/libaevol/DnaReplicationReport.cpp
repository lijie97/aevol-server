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

#if __cplusplus == 201103L
#include "make_unique.h"
#endif

#include "DnaReplicationReport.h"
#include "InsertionHT.h"
#include "ReplacementHT.h"

#include "Mutation.h"
#include "Duplication.h"
#include "Translocation.h"
#include "Inversion.h"
#include "Deletion.h"
#include "SmallDeletion.h"
#include "PointMutation.h"
#include "SmallInsertion.h"

namespace aevol {

DnaReplicationReport::DnaReplicationReport(const DnaReplicationReport& other) {
  for (auto& ht : other.ht_)
    add_HT(ht->Clone());
  for (auto& rear : other.rearrangements_)
    add_rear(rear->Clone());
  for (auto& mut : other.mutations_)
    add_local_mut(mut->Clone());
}

int32_t DnaReplicationReport::get_nb(MutationType t)  const {
  switch (t) {
    case S_MUT:
      assert(mutations_.size() ==
             static_cast<size_t>(_nb_mut[SWITCH] +
                                 _nb_mut[S_INS] +
                                 _nb_mut[S_DEL]));
      return mutations_.size();
    case REARR:
      assert(rearrangements_.size() ==
             static_cast<size_t>(_nb_mut[DUPL] +
                                 _nb_mut[DEL] +
                                 _nb_mut[TRANS] +
                                 _nb_mut[INV]));
      return rearrangements_.size();
    case H_T:
      assert(ht_.size() ==
             static_cast<size_t>(_nb_mut[INS_HT] +
                                 _nb_mut[REPL_HT]));
      return ht_.size();
    case INDEL:
      return _nb_mut[S_INS] + _nb_mut[S_DEL];
    default: // Simple mutation type.
      return _nb_mut[t];
  };
}

void DnaReplicationReport::add_mut(Mutation* mut) {
  if (mut->is_local_mut()) {
    add_local_mut(mut);
  }
  else if (mut->is_rear()) {
    add_rear(mut);
  }
  else if (mut->is_ht()) {
    add_HT(mut);
  }
}

void DnaReplicationReport::add_local_mut(Mutation* mut) {
  assert(mut->is_local_mut());
  std::unique_ptr<const LocalMutation> cmut = nullptr;
  switch(mut->get_mut_type()) {
    case SWITCH:
#if __cplusplus == 201103L
      cmut = make_unique<const PointMutation>(static_cast<PointMutation&>(*mut));
#else
      cmut = std::make_unique<const PointMutation>(static_cast<PointMutation&>(*mut));
#endif
      break;
    case S_DEL:
#if __cplusplus == 201103L
      cmut = make_unique<const SmallDeletion>(static_cast<SmallDeletion&>(*mut));
#else
      cmut = std::make_unique<const SmallDeletion>(static_cast<SmallDeletion&>(*mut));
#endif
      break;
    case S_INS:
#if __cplusplus == 201103L
      cmut = make_unique<const SmallInsertion>(static_cast<SmallInsertion&>(*mut));
#else
      cmut = std::make_unique<const SmallInsertion>(static_cast<SmallInsertion&>(*mut));
#endif
      break;
    default:
#if __cplusplus == 201103L
      cmut = make_unique<const SmallInsertion>(static_cast<SmallInsertion&>(*mut));
#else
      cmut = std::make_unique<const SmallInsertion>(static_cast<SmallInsertion&>(*mut));
#endif
      break;
  }
  mutations_.push_back(std::move(cmut));
  _nb_mut[mut->get_mut_type()]++;
}

void DnaReplicationReport::add_rear(Mutation* mut) {
  assert(mut->is_rear());

  std::unique_ptr<const Rearrangement> cmut = nullptr;
  switch(mut->get_mut_type()) {
    case DUPL:
#if __cplusplus == 201103L
      cmut = make_unique<const Duplication>(static_cast<Duplication&>(*mut));
#else
      cmut = std::make_unique<const Duplication>(static_cast<Duplication&>(*mut));
#endif
      break;
    case DEL:
#if __cplusplus == 201103L
      cmut = make_unique<const Deletion>(static_cast<Deletion&>(*mut));
#else
      cmut = std::make_unique<const Deletion>(static_cast<Deletion&>(*mut));
#endif
      break;
    case TRANS:
#if __cplusplus == 201103L
      cmut = make_unique<const Translocation>(static_cast<Translocation&>(*mut));
#else
      cmut = std::make_unique<const Translocation>(static_cast<Translocation&>(*mut));
#endif
      break;
    case INV:
#if __cplusplus == 201103L
      cmut = make_unique<const Inversion>(static_cast<Inversion&>(*mut));
#else
      cmut = std::make_unique<const Inversion>(static_cast<Inversion&>(*mut));
#endif
      break;
    default:
#if __cplusplus == 201103L
      cmut = make_unique<const Inversion>(static_cast<Inversion&>(*mut));
#else
      cmut = std::make_unique<const Inversion>(static_cast<Inversion&>(*mut));
#endif
      break;
  }
  rearrangements_.push_back(std::move(cmut));
  _nb_mut[mut->get_mut_type()]++;
}

void DnaReplicationReport::add_HT(Mutation* mut) {
  assert(mut->is_ht());

  std::unique_ptr<const HT> cmut = nullptr;
  switch(mut->get_mut_type()) {
    case INS_HT:
#if __cplusplus == 201103L
      cmut = make_unique<const InsertionHT>(static_cast<InsertionHT&>(*mut));
#else
      cmut = std::make_unique<const InsertionHT>(static_cast<InsertionHT&>(*mut));
#endif
      break;
    case REPL_HT:
#if __cplusplus == 201103L
      cmut = make_unique<const ReplacementHT>(static_cast<ReplacementHT&>(*mut));
#else
      cmut = std::make_unique<const ReplacementHT>(static_cast<ReplacementHT&>(*mut));
#endif
      break;
    default:
#if __cplusplus == 201103L
      cmut = make_unique<const ReplacementHT>(static_cast<ReplacementHT&>(*mut));
#else
      cmut = std::make_unique<const ReplacementHT>(static_cast<ReplacementHT&>(*mut));
#endif
      break;
  }
  ht_.push_back(std::move(cmut));
  _nb_mut[mut->get_mut_type()]++;
}


/// Useful when we inspect a tree file
/// because stats are not saved in the file.
void DnaReplicationReport::compute_stats( void )
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
  
  for (const auto& ht : ht_) {
    assert(ht->get_mut_type() == INS_HT or
           ht->get_mut_type() == REPL_HT);
    _nb_mut[ht->get_mut_type()]++;
  }

  for (const auto& rear : rearrangements_) {
    assert(rear->get_mut_type() == DUPL or
           rear->get_mut_type() == DEL or
           rear->get_mut_type() == TRANS or
           rear->get_mut_type() == INV);
    _nb_mut[rear->get_mut_type()]++;
  }

  for (const auto& mut : mutations_) {
    assert(mut->get_mut_type() == SWITCH or
           mut->get_mut_type() == S_INS or
           mut->get_mut_type() == S_DEL);
    _nb_mut[mut->get_mut_type()]++;
  }
}

void DnaReplicationReport::write_to_tree_file(gzFile tree_file) const {
  // Write the mutations and rearrangements undergone during replication
  // Store HT
  int32_t nb_HT = get_nb(H_T);
  gzwrite(tree_file, &nb_HT, sizeof(nb_HT));
  for (const auto& ht : ht_) {
    switch(ht->get_mut_type()) {
      case INS_HT:
        ht->save(tree_file);
        break;
      case REPL_HT:
        ht->save(tree_file);
        break;
      default:
        ht->save(tree_file);
        break;
    }
  }


  // Store rearrangements
  int32_t nb_rears = get_nb(REARR);
  gzwrite(tree_file, &nb_rears, sizeof(nb_rears));
  for (const auto& rear : rearrangements_) {
    switch(rear->get_mut_type()) {
      case DUPL:
        rear->save(tree_file);
        break;
      case DEL:
        rear->save(tree_file);
        break;
      case TRANS:
        rear->save(tree_file);
        break;
      case INV:
        rear->save(tree_file);
        break;
      default:
        rear->save(tree_file);
        break;
    }
  }

  // Store mutations
  int32_t nb_muts = get_nb(S_MUT);
  gzwrite(tree_file, &nb_muts, sizeof(nb_muts));
  for (const auto& mut : mutations_)
    switch(mut->get_mut_type()) {
      case SWITCH:
        mut->save(tree_file);
        break;
      case S_DEL:
        mut->save(tree_file);
        break;
      case S_INS:
        mut->save(tree_file);
        break;
      default:
        mut->save(tree_file);
        break;
    }
}

void DnaReplicationReport::read_from_tree_file(gzFile tree_file) {
  int32_t nb_rears, nb_muts, nb_HT;

  gzread(tree_file, &nb_HT, sizeof(nb_HT));
  for (int i = 0 ; i < nb_HT ; i++)
    add_HT(Mutation::Load(tree_file));

  gzread(tree_file, &nb_rears, sizeof(nb_rears));
  for (int i = 0 ; i < nb_rears ; i++)
    add_rear(Mutation::Load(tree_file));

  gzread(tree_file, &nb_muts, sizeof(nb_muts));
  for(int i = 0 ; i < nb_muts ; i++)
    add_mut(Mutation::Load(tree_file));
}
} // namespace aevol
