//
// Created by arrouan on 01/09/16.
//

//#ifdef __OPENMP_GPU
#ifndef RAEVOL_CUDA_ALGORITHM_CUDA_H
#define RAEVOL_CUDA_ALGORITHM_CUDA_H

#include <vector>
#include <list>
#include  <functional>

#include "raevol/Protein_R.h"
#include "raevol/Rna_R.h"
#include "ae_enums.h"

using std::vector;
using std::list;

namespace aevol {
  class algorithm_cuda {
   public:
      static std::list<Protein_R>::iterator find_if_protein(
          std::list<Protein_R>::iterator first,
          std::list<Protein_R>::iterator last,
          int32_t shine_dal_pos);

      static std::list<Rna_R>::iterator find_if_rna_1(
          std::list<Rna_R>::iterator first,
          std::list<Rna_R>::iterator last,
          int32_t pos, aevol::Strand);

      static std::list<Rna_R>::iterator find_if_rna_2(
          std::list<Rna_R>::iterator first,
          std::list<Rna_R>::iterator last,
          int32_t pos);

      static std::list<Rna_R>::iterator find_if_rna_3(
          std::list<Rna_R>::iterator first,
          std::list<Rna_R>::iterator last,
          int32_t pos);

      static std::list<Rna_R>::iterator find_if_rna_4(
          std::list<Rna_R>::iterator first,
          std::list<Rna_R>::iterator last,
          int32_t pos);
  };
}


#endif //RAEVOL_CUDA_ALGORITHM_CUDA_H
//#endif
