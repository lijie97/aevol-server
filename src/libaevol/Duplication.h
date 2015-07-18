//
// Created by dparsons on 02/07/15.
//

#ifndef AEVOL_DUPLICATION_H__
#define AEVOL_DUPLICATION_H__


// ============================================================================
//                                   Includes
// ============================================================================

#include "Rearrangement.h"
#include "sys/types.h"

namespace aevol {

/**
 *
 */
class Duplication : public Rearrangement {

 public :
  // ==========================================================================
  //                               Constructors
  // ==========================================================================
  Duplication() = default; //< Default ctor
  Duplication(const Duplication&) = default; //< Copy ctor
  Duplication(Duplication&&) = delete; //< Move ctor
  Duplication(int32_t pos1, int32_t pos2, int32_t pos3,
              int32_t length, int16_t align_score = -1);

  virtual Mutation* Clone() const override { return new Duplication(*this); };

  // ==========================================================================
  //                                Destructor
  // ==========================================================================
  virtual ~Duplication() = default; //< Destructor

  // ==========================================================================
  //                                Operators
  // ==========================================================================
  /// Copy assignment
  Duplication& operator=(const Duplication& other) = default;

  /// Move assignment
  Duplication& operator=(Duplication&& other) = delete;

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
    return DUPL;
  };

  int32_t pos1() const {
    return pos1_;
  }

  int32_t pos2() const {
    return pos2_;
  }

  int32_t pos3() const {
    return pos3_;
  }

  int32_t length() const {
    return length_;
  }

  int16_t align_score() const {
    return align_score_;
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
  int32_t pos1_, pos2_, pos3_;
  int32_t length_;
  int16_t align_score_ = -1;
};

} // namespace aevol
#endif //AEVOL_DUPLICATION_H__
