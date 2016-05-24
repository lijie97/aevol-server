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

class NewIndivEvent {
 public:
    NewIndivEvent(Individual* childx, Individual* parentx, int xx, int yx) {
      child = childx;
      parent = parentx;
      x = xx;
      y = yx;
    }

    Individual* child;
    Individual* parent;

    int x;
    int y;
};

class EndReplicationEvent {
 public:
    EndReplicationEvent(Individual* childx, int xx, int yx) {
      child = childx;
      x = xx;
      y = yx;
    }

    Individual* child;

    int x;
    int y;
};
}

#endif //AEVOL_OBSERVABLEEVENT_H_
