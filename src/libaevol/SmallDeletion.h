//
// Created by dparsons on 02/07/15.
//

#ifndef AEVOL_SMALLDELETION_H__
#define AEVOL_SMALLDELETION_H__


// ============================================================================
//                                   Includes
// ============================================================================
#include <cstdint>

#include "LocalMutation.h"

namespace aevol {

/**
 *
 */
class SmallDeletion : public LocalMutation {

 public :
  // ==========================================================================
  //                               Constructors
  // ==========================================================================
  SmallDeletion() = default; //< Default ctor
  SmallDeletion(const SmallDeletion&) = default; //< Copy ctor
  SmallDeletion(SmallDeletion&&) = default; //< Move ctor
  SmallDeletion(int32_t pos, int16_t length);

  virtual Mutation* Clone() const override { return new SmallDeletion(*this); };

  // ==========================================================================
  //                                Destructor
  // ==========================================================================
  virtual ~SmallDeletion() = default; //< Destructor

  // ==========================================================================
  //                                Operators
  // ==========================================================================
  /// Copy assignment
  SmallDeletion& operator=(const SmallDeletion& other) = default;

  /// Move assignment
  SmallDeletion& operator=(SmallDeletion&& other) = delete;

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
    return S_DEL;
  };

  int32_t pos() const {
    return pos_;
  };

  int16_t length() const {
    return length_;
  };

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
  int16_t length_;
};

} // namespace aevol
#endif //AEVOL_SMALLDELETION_H__
