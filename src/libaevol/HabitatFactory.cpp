//
// Created by arrouan on 19/11/15.
//

#include "HabitatFactory.h"

namespace aevol {


std::unique_ptr<Habitat> HabitatFactory::create_unique_habitat(Habitat& habitat,
                                                      bool share_phenotypic_target) {

#ifndef __REGUL
#if __cplusplus == 201103L
  return make_unique<Habitat> (habitat, share_phenotypic_target);
#else
  return std::make_unique<Habitat> (habitat, share_phenotypic_target);
#endif
#else
#if __cplusplus == 201103L
  return make_unique<Habitat_R> (dynamic_cast<Habitat_R&>(habitat), share_phenotypic_target);
#else
  return std::make_unique<Habitat_R> (dynamic_cast<Habitat_R&>(habitat), share_phenotypic_target);
#endif
#endif
}

std::unique_ptr<Habitat> 
HabitatFactory::create_unique_habitat(gzFile backup_file,
                 					  PhenotypicTargetHandler* phenotypic_target_handler) {
#ifndef __REGUL
#if __cplusplus == 201103L
  return make_unique<Habitat>(backup_file, phenotypic_target_handler);
#else
  return std::make_unique<Habitat>(backup_file, phenotypic_target_handler);
#endif
#else
#if __cplusplus == 201103L
  return make_unique<Habitat_R>(backup_file, dynamic_cast<PhenotypicTargetHandler_R*>(phenotypic_target_handler));
#else
  return std::make_unique<Habitat_R>(backup_file, dynamic_cast<PhenotypicTargetHandler_R*>(phenotypic_target_handler));
#endif

 #endif

}




}
