//
// Created by dparsons on 02/07/15.
//

#ifndef AEVOL_INVERSION_H__
#define AEVOL_INVERSION_H__


// ============================================================================
//                                   Includes
// ============================================================================

#include "Rearrangement.h"

namespace aevol {

/**
 *
 */
class Inversion : public Rearrangement{

 public :
  // ==========================================================================
  //                               Constructors
  // ==========================================================================
  Inversion() = default; //< Default ctor
  Inversion(const Inversion&) = default; //< Copy ctor
  Inversion(Inversion&&) = delete; //< Move ctor
  Inversion(int32_t pos1, int32_t pos2,
            int32_t length, int16_t align_score = -1);

  virtual Mutation* Clone() const override { return new Inversion(*this); };

  // ==========================================================================
  //                                Destructor
  // ==========================================================================
  virtual ~Inversion() = default; //< Destructor

  // ==========================================================================
  //                                Operators
  // ==========================================================================
  /// Copy assignment
  Inversion& operator=(const Inversion& other) = default;

  /// Move assignment
  Inversion& operator=(Inversion&& other) = delete;

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
    return INV;
  };

  int32_t pos1() const {
    return pos1_;
  }

  int32_t pos2() const {
    return pos2_;
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
  int32_t pos1_, pos2_;
  int32_t length_;
  int16_t align_score_ = -1;
};

} // namespace aevol
#endif //AEVOL_INVERSION_H__
