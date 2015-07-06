//
// Created by dparsons on 02/07/15.
//

#ifndef AEVOL_TRANSLOCATION_H__
#define AEVOL_TRANSLOCATION_H__


// ============================================================================
//                                   Includes
// ============================================================================

#include "Rearrangement.h"

namespace aevol {

/**
 *
 */
class Translocation : public Rearrangement {

 public :
  // ==========================================================================
  //                               Constructors
  // ==========================================================================
  Translocation() = default; //< Default ctor
  Translocation(const Translocation&) = default; //< Copy ctor
  Translocation(Translocation&&) = delete; //< Move ctor

  Translocation(int32_t pos1, int32_t pos2, int32_t pos3, int32_t pos4,
                int32_t length,
                bool invert,
                int16_t align_score_1 = -1, int16_t align_score_2 = -1);

  virtual Mutation* Clone() const override { return new Translocation(*this); };

  // ==========================================================================
  //                                Destructor
  // ==========================================================================
  virtual ~Translocation() = default; //< Destructor

  // ==========================================================================
  //                                Operators
  // ==========================================================================
  /// Copy assignment
  Translocation& operator=(const Translocation& other) = default;

  /// Move assignment
  Translocation& operator=(Translocation&& other) = delete;

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
    return TRANS;
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

  int32_t pos4() const {
    return pos4_;
  }

  int32_t length() const {
    return length_;
  }

  bool invert() const {
    return invert_;
  }

  int16_t align_score_1() const {
    return align_score_1_;
  }

  int16_t align_score_2() const {
    return align_score_2_;
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
  int32_t pos1_, pos2_, pos3_, pos4_;
  int32_t length_;
  bool invert_;
  int16_t align_score_1_ = -1;
  int16_t align_score_2_ = -1;
};

} // namespace aevol
#endif //AEVOL_TRANSLOCATION_H__
