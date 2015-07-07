//
// Created by dparsons on 01/07/15.
//

#ifndef AEVOL_SMALLINSERTION_H__
#define AEVOL_SMALLINSERTION_H__


// ============================================================================
//                                   Includes
// ============================================================================
#include <cstdint>

#include "LocalMutation.h"

namespace aevol {

/**
 *
 */
class SmallInsertion : public LocalMutation {

 public :
  // ==========================================================================
  //                               Constructors
  // ==========================================================================
  SmallInsertion() = default; //< Default ctor
  SmallInsertion(const SmallInsertion&); //< Copy ctor
  SmallInsertion(SmallInsertion&&) = delete; //< Move ctor
  SmallInsertion(int32_t pos, int32_t length, const char* seq);

  virtual Mutation* Clone() const override { return new SmallInsertion(*this); };

  // ==========================================================================
  //                                Destructor
  // ==========================================================================
  virtual ~SmallInsertion(); //< Destructor

  // ==========================================================================
  //                                Operators
  // ==========================================================================
  /// Copy assignment
  SmallInsertion& operator=(const SmallInsertion& other);

  /// Move assignment
  SmallInsertion& operator=(SmallInsertion&& other) = delete;

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
    return S_INS;
  };

  int32_t pos() const {
    return pos_;
  };

  int32_t length() const {
    return length_;
  };

  char* seq() const {
    return seq_;
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
  int32_t length_;
  char* seq_;
};

} // namespace aevol
#endif //AEVOL_SMALLINSERTION_H__
