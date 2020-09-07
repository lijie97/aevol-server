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


#include "Individual_7.h"
#include "DnaMutator.h"
#include "ExpManager.h"
#include "Stats_7.h"
#include "7/DynTab_Metadata.h"
#include "7/List_Metadata.h"
#include "7/Map_Metadata.h"
#include "7/Protein_7.h"

#include "Fuzzy.h"

#include <algorithm>

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


  if (exp_m_->exp_s()->get_simd_metadata_flavor() == SIMDMetadataFlavor::STD_MAP)
      metadata_ = new Map_Metadata(this,dynamic_cast<Map_Metadata*>(clone->metadata_));
  else if (exp_m_->exp_s()->get_simd_metadata_flavor() == SIMDMetadataFlavor::DYN_TAB)
      metadata_ = new DynTab_Metadata(this,dynamic_cast<DynTab_Metadata*>(clone->metadata_));
  else if (exp_m_->exp_s()->get_simd_metadata_flavor() == SIMDMetadataFlavor::STD_LIST)
      metadata_ = new List_Metadata(this,dynamic_cast<List_Metadata*>(clone->metadata_));

  fitness = clone->fitness;
  metaerror = clone->metaerror;

}

Individual_7::~Individual_7() {
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
