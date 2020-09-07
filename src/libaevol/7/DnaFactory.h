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

#ifndef AEVOL_DNAFACTORY_H
#define AEVOL_DNAFACTORY_H

#include "Dna_7.h"
namespace aevol {
    enum DnaFactory_Policy {
        FIRST = 0,
        FIRSTFIT = 1,
        BESTFIT = 2
    };

    class DnaFactory {
    public:
     DnaFactory(DnaFactory_Policy policy, int pool_size, int init_size) {
            policy_ = policy;
            pool_size_ = pool_size;
            init(init_size);
        }

        ~DnaFactory() {
            for (std::list<Dna_7 *>::iterator it_dna = list_unused_dna_.begin();
                 it_dna != list_unused_dna_.end(); it_dna++) {
                delete (*(it_dna));
            }
        }

        void init(int init_size);

        void stats() {
            int total_length_ = 0;
            for (std::list<Dna_7 *>::iterator it_dna = list_unused_dna_.begin();
                 it_dna != list_unused_dna_.end(); it_dna++) {
                total_length_ += (*it_dna)->nb_block()*BLOCK_SIZE* sizeof(char);
            }
            total_length_/=1000000;
            printf("DNA_FACTORY_STATS -- Number of DNAs %ld - Combined size %d Mb\n",list_unused_dna_.size(),total_length_);
        }

        Dna_7 *get_dna(int request_size);

        void give_back(Dna_7 *dna);


    private:
        std::list<Dna_7 *> list_unused_dna_;
        DnaFactory_Policy policy_;
        int pool_size_;
    };

}
#endif //AEVOL_DNAFACTORY_H
