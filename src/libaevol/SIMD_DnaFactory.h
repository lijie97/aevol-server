//
// Created by arrouan on 20/03/20.
//

#ifndef AEVOL_SIMD_DNAFACTORY_H
#define AEVOL_SIMD_DNAFACTORY_H

#include "Dna_SIMD.h"
namespace aevol {
    enum DnaFactory_Policy {
        FIRSTFIT = 1,
        BESTFIT = 2
    };

    class SIMD_DnaFactory {
    public:
        SIMD_DnaFactory(DnaFactory_Policy policy, int pool_size, int init_size) {
            policy_ = policy;
            pool_size_ = pool_size;
            init(init_size);
        }

        void init(int init_size);

        Dna_SIMD *get_dna(int request_size);

        void give_back(Dna_SIMD *dna);


    private:
        std::list<Dna_SIMD *> list_unused_dna_;
        DnaFactory_Policy policy_;
        int pool_size_;
    };

}
#endif //AEVOL_SIMD_DNAFACTORY_H
