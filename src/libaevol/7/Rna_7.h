//
// Created by arrouan on 07/09/2020.
//

#ifndef AEVOL_RNA_7_H
#define AEVOL_RNA_7_H


#include <list>
#include <cstdint>

namespace aevol {
class Rna_7 {
 public:
  Rna_7(){};
  Rna_7(int32_t t_begin,
        int32_t t_end,
        int8_t t_leading_lagging,
        double t_e,
        int32_t t_length) {
    begin             = t_begin;
    end               = t_end;
    leading_lagging   = t_leading_lagging;
    e                 = t_e;
    length            = t_length;
    is_coding_        = false;
    is_init_          = true;
    start_prot_count_ = 0;
  }

  ~Rna_7() {}

  int32_t begin;
  int32_t end;
  int8_t leading_lagging; // 0 = leading, 1 = lagging
  double e;
  std::list<int32_t> start_prot;
  int32_t start_prot_count_;
  int32_t length;
  bool is_coding_;

  bool is_init_ = false;
};

}

#endif //AEVOL_RNA_7_H
