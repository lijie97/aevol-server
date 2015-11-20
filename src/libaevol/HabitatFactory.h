//
// Created by arrouan on 19/11/15.
//

#ifndef AEVOL_REGUL_HABITATFACTORY_H
#define AEVOL_REGUL_HABITATFACTORY_H

#include "Habitat.h"
#include "raevol/Habitat_R.h"

namespace aevol {


class HabitatFactory {
 public:
  static std::unique_ptr<Habitat> create_unique_habitat(Habitat& habitat,
                                                        bool share_phenotypic_target);

};

}
#endif //AEVOL_REGUL_HABITATFACTORY_H
