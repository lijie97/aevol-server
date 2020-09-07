//
// Created by arrouan on 27/07/17.
//

#include "Individual_7.h"

#include "DnaMutator.h"
#include "Dna_7.h"
#include "ExpManager.h"
#include "HybridFuzzy.h"
#include "Stats_7.h"
//#include "Abstract_Metadata.h"
#include "7/DynTab_Metadata.h"
#include "7/List_Metadata.h"
#include "7/Map_Metadata.h"
#include "7/Vector_Fuzzy.h"
#include "7/Protein_7.h"

#include "Fuzzy.h"
#include "stdint.h"
#include "stdio.h"

#include <algorithm>
#include <chrono>
#include <err.h>
#include <omp.h>
#include <sys/stat.h>

namespace aevol {

/** Individual_7 Constructor and Destructor **/
Individual_7::Individual_7(ExpManager* exp_m, double w_max,
                           DnaFactory* dna_factory) {
        exp_m_ = exp_m;
        w_max_ = w_max;

        if (exp_m_->exp_s()->get_simd_metadata_flavor() == SIMDMetadataFlavor::STD_MAP)
            metadata_ = new Map_Metadata(this);
        else if (exp_m_->exp_s()->get_simd_metadata_flavor() == SIMDMetadataFlavor::DYN_TAB)
            metadata_ = new DynTab_Metadata(this);
        else if (exp_m_->exp_s()->get_simd_metadata_flavor() == SIMDMetadataFlavor::STD_LIST)
            metadata_ = new List_Metadata(this);

        dna_factory_ = dna_factory;
    }

    Individual_7::Individual_7(ExpManager* exp_m,
                               Individual_7* clone,
                               DnaFactory* dna_factory) {
    w_max_ = clone->w_max_;

  exp_m_ = exp_m;

  usage_count_ = 1;
  dna_ = dna_factory->get_dna(clone->dna_->length());
  //printf("DNA Factory -- %p %p\n",dna_,dna_->data_);
  dna_->set_indiv(clone->dna_,this);

  dna_factory_ = dna_factory;


  //promoters.resize(clone->promoters.size());
  if (exp_m_->exp_s()->get_simd_metadata_flavor() == SIMDMetadataFlavor::STD_MAP)
      metadata_ = new Map_Metadata(this,dynamic_cast<Map_Metadata*>(clone->metadata_));
  else if (exp_m_->exp_s()->get_simd_metadata_flavor() == SIMDMetadataFlavor::DYN_TAB)
      metadata_ = new DynTab_Metadata(this,dynamic_cast<DynTab_Metadata*>(clone->metadata_));
  else if (exp_m_->exp_s()->get_simd_metadata_flavor() == SIMDMetadataFlavor::STD_LIST)
      metadata_ = new List_Metadata(this,dynamic_cast<List_Metadata*>(clone->metadata_));

  fitness = clone->fitness;
  metaerror = clone->metaerror;
  //leading_prom_pos = clone->leading_prom_pos;
  //lagging_prom_pos = clone->lagging_prom_pos;

}

Individual_7::~Individual_7() {
  //printf("DESTRUCTOR SIMD INDIV -- Metadata %p\n",metadata_);

  dna_factory_->give_back(dna_);

#ifndef PHENOTYPE_VECTOR
  delete phenotype;
#endif

  delete metadata_;

  clearAllObserver();
}

/**
 * We need some index for the promoter optimization
 */
void Individual_7::rebuild_index() {
        if (exp_m_->exp_s()->get_simd_metadata_flavor() == SIMDMetadataFlavor::STD_MAP)
            dynamic_cast<Map_Metadata*>(metadata_)->rebuild_index();
}

    bool Protein_7::operator<(const Protein_7& other){
        return (h <  other.h)
               || (h == other.h && m < other.m)
               || (h == other.h && m == other.m && w < other.w);
    }
}
