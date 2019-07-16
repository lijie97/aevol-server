//
// Created by arrouan on 19/11/15.
//

#ifndef AEVOL_REGUL_HABITATFACTORY_H
#define AEVOL_REGUL_HABITATFACTORY_H

#include "Habitat.h"
#include "raevol/Habitat_R.h"
#if __cplusplus == 201103L
#include "make_unique.h"
#endif

namespace aevol {


class HabitatFactory {
 public:
  static std::unique_ptr<Habitat> create_unique_habitat(Habitat& habitat,
                                                        bool share_phenotypic_target);
  static std::unique_ptr<Habitat_R> create_unique_habitat(Habitat_R& habitat,
														  bool share_phenotypic_target);
  static
#ifdef __REGUL
    std::unique_ptr<Habitat_R>
#else
    std::unique_ptr<Habitat>
#endif
    create_unique_habitat(gzFile backup_file,
                 					  PhenotypicTargetHandler* phenotypic_target_handler);

};

}
#endif //AEVOL_REGUL_HABITATFACTORY_H
