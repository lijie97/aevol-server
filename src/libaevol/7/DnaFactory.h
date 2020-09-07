//
// Created by arrouan on 20/03/20.
//

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
