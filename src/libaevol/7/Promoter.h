//
// Created by arrouan on 07/09/2020.
//

#ifndef AEVOL_PROMOTER_H
#define AEVOL_PROMOTER_H

#include <map>

namespace aevol {

class PromoterStruct {
 public:
  PromoterStruct(int32_t t_pos, int8_t t_error, bool lead) {
    pos                = t_pos;
    error              = t_error;
    leading_or_lagging = lead;
  }

  PromoterStruct(const PromoterStruct& clone) {
    pos                = clone.pos;
    error              = clone.error;
    leading_or_lagging = clone.leading_or_lagging;
  }

  PromoterStruct(PromoterStruct* clone) {
    pos                = clone->pos;
    error              = clone->error;
    leading_or_lagging = clone->leading_or_lagging;
  }

  int32_t pos  = -1;
  int8_t error = -1;
  bool leading_or_lagging; // TRUE = leading / FALSE = lagging
};

class PromoterList {
  PromoterList() = default;

  std::map<int32_t, int> leading_prom_pos;
  std::map<int32_t, int> lagging_prom_pos;
};
}

#endif //AEVOL_PROMOTER_H
