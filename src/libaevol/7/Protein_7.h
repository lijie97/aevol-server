//
// Created by arrouan on 07/09/2020.
//

#ifndef AEVOL_PROTEIN_7_H
#define AEVOL_PROTEIN_7_H


#include <cstdint>

namespace aevol {
class Protein_7 {
 public:
  Protein_7(){};
  Protein_7(int32_t t_protein_start,
            int32_t t_protein_end,
            int32_t t_protein_length,
            int8_t t_leading_lagging,
            double t_e) {
    protein_start   = t_protein_start;
    protein_end     = t_protein_end;
    protein_length  = t_protein_length;
    leading_lagging = t_leading_lagging;
    e               = t_e;
    is_init_        = true;
  }

  bool operator<(const Protein_7& other);

  int32_t protein_start;
  int32_t protein_end;
  int32_t protein_length;
  int8_t leading_lagging; // 0 = leading, 1 = lagging
  double m;
  double w;
  double h;
  double e;
  bool is_functional;

  bool is_init_ = false;
};
}

#endif //AEVOL_PROTEIN_7_H
