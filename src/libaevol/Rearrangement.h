//
// Created by dparsons on 02/07/15.
//

#ifndef AEVOL_REARRANGEMENT_H__
#define AEVOL_REARRANGEMENT_H__


// ============================================================================
//                                   Includes
// ============================================================================
#include "Mutation.h"

namespace aevol {

/**
 * Abstract base class for local mutations
 */
class Rearrangement : public Mutation {

 public :
  // ==========================================================================
  //                               Constructors
  // ==========================================================================
  virtual Mutation* Clone() const override = 0;

  // ==========================================================================
  //                                Destructor
  // ==========================================================================

  // ==========================================================================
  //                                Operators
  // ==========================================================================

  // ==========================================================================
  //                              Public Methods
  // ==========================================================================
  virtual void save(gzFile backup_file) const override = 0;

  // ==========================================================================
  //                                 Getters
  // ==========================================================================
  virtual MutationType get_mut_type() const override = 0;
  bool is_rear() const override { return true; };

  // ==========================================================================
  //                                 Setters
  // ==========================================================================

 protected :
  // ==========================================================================
  //                            Protected Methods
  // ==========================================================================

  // ==========================================================================
  //                               Attributes
  // ==========================================================================
};

} // namespace aevol
#endif //AEVOL_REARRANGEMENT_H__
