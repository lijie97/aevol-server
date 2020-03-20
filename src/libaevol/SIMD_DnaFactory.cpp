//
// Created by arrouan on 20/03/20.
//
#include <list>
#include "SIMD_DnaFactory.h"
namespace aevol {
    void SIMD_DnaFactory::init(int init_size) {
        for (int i = 0; i < pool_size_; i++) {
            list_unused_dna_.push_back(new Dna_SIMD(init_size,this));
        }
    }

    Dna_SIMD *SIMD_DnaFactory::get_dna(int request_size) {
        request_size++; // Count the end \0
        if (policy_ == DnaFactory_Policy::FIRSTFIT) {
            Dna_SIMD *pop = nullptr;
            //printf("DNA Factory -- Length %d -- %d\n",request_size,list_unused_dna_.size());

            #pragma omp critical(pop_dna)
            {
                if (list_unused_dna_.empty()) {
                    //printf("DnaFactory -- Empty Factory !!");
                    pop = new Dna_SIMD(request_size,this);
                } else {
                    pop = list_unused_dna_.front();
                    list_unused_dna_.pop_front();
                }
            }
            return pop;
        }
    }

    void SIMD_DnaFactory::give_back(Dna_SIMD *dna) {
        #pragma omp critical(pop_dna)
        {
            list_unused_dna_.push_back(dna);
        }
    }
}