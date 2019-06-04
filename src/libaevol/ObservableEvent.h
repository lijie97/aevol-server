// ****************************************************************************
//
//          Aevol - An in silico experimental evolution platform
//
// ****************************************************************************
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
//
// ****************************************************************************

#ifndef AEVOL_OBSERVABLEEVENT_H_
#define AEVOL_OBSERVABLEEVENT_H_

/**
 *
 */


enum ObservableEvent {
  NEW_INDIV,
  MUTATION,
  END_REPLICATION,
  END_GENERATION
};

namespace  aevol {

class Individual;
class Internal_SIMD_Struct;

class NewIndivEvent {
 public:
    NewIndivEvent(Individual* childx, Individual* parentx, int xx, int yx, int indiv_id, int parent_id) {
      child = childx;
      parent = parentx;
      simd_parent = nullptr;
      simd_child = nullptr;
        indiv_id_ = indiv_id;
        parent_id_ = parent_id;
      x = xx;
      y = yx;
    }

    NewIndivEvent(Internal_SIMD_Struct* childx, Internal_SIMD_Struct* parentx, int xx, int yx, int indiv_id, int parent_id) {
        simd_child = childx;
        simd_parent = parentx;
        child = nullptr;
        parent = nullptr;

        indiv_id_ = indiv_id;
        parent_id_ = parent_id;
        x = xx;
        y = yx;
    }


    Individual* child;
    Individual* parent;
    Internal_SIMD_Struct* simd_child;
    Internal_SIMD_Struct* simd_parent;

    int x;
    int y;

    int indiv_id_;
    int parent_id_;
};

class EndReplicationEvent {
 public:
    EndReplicationEvent(Individual* childx, int xx, int yx) {
      child = childx;
      x = xx;
      y = yx;
    }


    EndReplicationEvent(Internal_SIMD_Struct* childx, int xx, int yx) {
        simd_child = childx;
        x = xx;
        y = yx;
    }

    Individual* child;
    Internal_SIMD_Struct* simd_child;

    int x;
    int y;
};
}

#endif //AEVOL_OBSERVABLEEVENT_H_
