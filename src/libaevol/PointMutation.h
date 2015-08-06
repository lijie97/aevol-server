//
// Created by dparsons on 01/07/15.
//

#ifndef AEVOL_POINTMUTATION_H__
#define AEVOL_POINTMUTATION_H__


// ============================================================================
//                                   Includes
// ============================================================================
#include <cstdint>

#include "LocalMutation.h"

namespace aevol {

/**
 *
 */
class PointMutation : public LocalMutation {

 public :
  // ==========================================================================
  //                               Constructors
  // ==========================================================================
  PointMutation() = default; //< Default ctor
  PointMutation(const PointMutation&) = default; //< Copy ctor
  PointMutation(PointMutation&&) = delete; //< Move ctor
  PointMutation(int32_t pos);

  virtual Mutation* Clone() const override { return new PointMutation(*this); };

  // ==========================================================================
  //                                Destructor
  // ==========================================================================
  virtual ~PointMutation() = default; //< Destructor

  // ==========================================================================
  //                                Operators
  // ==========================================================================
  /// Copy assignment
  PointMutation& operator=(const PointMutation& other) = default;
  /// Move assignment
  PointMutation& operator=(PointMutation&& other) = delete;

  // ==========================================================================
  //                              Public Methods
  // ==========================================================================
  virtual void save(gzFile backup_file) const override;
  virtual void load(gzFile backup_file) override;
  void get_generic_description_string(char* str) const override;

  // ==========================================================================
  //                                 Getters
  // ==========================================================================
  virtual MutationType get_mut_type() const override {
    return SWITCH;
  };

  int32_t pos() const {
    return pos_;
  }

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
  int32_t pos_;
};

} // namespace aevol
#endif //AEVOL_POINTMUTATION_H__
