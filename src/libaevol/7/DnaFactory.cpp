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

#include "DnaFactory.h"

#include <list>
namespace aevol {
    void DnaFactory::init(int init_size) {
        for (int i = 0; i < pool_size_; i++) {
            list_unused_dna_.push_back(new Dna_7(init_size,this));
        }
    }

    Dna_7 *DnaFactory::get_dna(int request_size) {
        request_size++; // Count the end \0
        int req_block = Dna_7::nb_blocks(request_size);
        if (policy_ == DnaFactory_Policy::FIRST) {
          Dna_7 *pop = nullptr;
            //printf("DNA Factory -- Length %d -- %d\n",request_size,list_unused_dna_.size());

            #pragma omp critical(pop_dna)
            {
                if (list_unused_dna_.empty()) {
                    //printf("DnaFactory -- Empty Factory !!");
                    pop = new Dna_7(request_size,this);
                } else {
                    pop = list_unused_dna_.front();
                    list_unused_dna_.pop_front();
                }
            }
            return pop;
        } else if (policy_ == DnaFactory_Policy::FIRSTFIT) {
          Dna_7 *pop = nullptr;

#pragma omp critical(pop_dna)
            {
                if (list_unused_dna_.empty()) {
                    //printf("DnaFactory -- Empty Factory !!");
                    pop = new Dna_7(request_size,this);
                } else {
                    std::list<Dna_7 *>::iterator found_it;
                    for (auto it = list_unused_dna_.begin(); it != list_unused_dna_.end(); it++) {
                        if ((*it)->nb_block() >= req_block) {
                            found_it = it;
                            pop = (*it);
                            break;
                        }
                    }

                    if (pop == nullptr) {
                        pop = list_unused_dna_.front();
                        list_unused_dna_.pop_front();
                    } else {
                        list_unused_dna_.erase(found_it);
                    }
                }
            }
            return pop;
        }
    }

    void DnaFactory::give_back(Dna_7 *dna) {
        #pragma omp critical(pop_dna)
        {
            list_unused_dna_.push_back(dna);
        }
    }
}