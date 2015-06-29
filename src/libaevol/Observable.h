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

#ifndef OBSERVER_OBSERVABLE_H__
#define OBSERVER_OBSERVABLE_H__


// ============================================================================
//                                   Includes
// ============================================================================
#include "ObservableEvent.h"
#include "Observer.h"

#include <list>
#include <map>

using std::list;
using std::map;




/**
 *
 */
  class Observable {
   public :
  // ==========================================================================
  //                               Constructors
  // ==========================================================================
  Observable(void) = default; //< Default ctor
  Observable(const Observable&) = delete; //< Copy ctor
  Observable(Observable&&) = delete; //< Move ctor

  // ==========================================================================
  //                                Destructor
  // ==========================================================================
  virtual ~Observable(void) = default; //< Destructor

  // ==========================================================================
  //                                Operators
  // ==========================================================================
  Observable& operator=(const Observable& other); //< Copy assignment
  Observable& operator=(const Observable&& other); //< Move assignment

  // ==========================================================================
  //                              Public Methods
  // ==========================================================================
  void addObserver(Observer* o, ObservableEvent e) {
    observers_[e].emplace_back(o);
  };
  void deleteObserver(Observer* o, ObservableEvent e) {
    observers_[e].remove(o);
  };
  void notifyObservers(ObservableEvent e, void* arg = nullptr);

  // ==========================================================================
  //                                 Getters
  // ==========================================================================

  // ==========================================================================
  //                                 Setters
  // ==========================================================================

 protected :
  // ==========================================================================
  //                            Protected Methods
  // ==========================================================================

  // ==========================================================================
  //                               Attributes
  // ==========================================================================
  map<ObservableEvent, list<Observer*> > observers_;
};


#endif //OBSERVER_OBSERVABLE_H__
