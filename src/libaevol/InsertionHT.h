//
// Created by dparsons on 03/07/15.
//

#ifndef AEVOL_INSERTION_HT_H__
#define AEVOL_INSERTION_HT_H__


// ============================================================================
//                                   Includes
// ============================================================================

#include "HT.h"
#include "VisAVis.h"

namespace aevol {

/**
 *
 */
class InsertionHT : public HT {

 public :
  // ==========================================================================
  //                               Constructors
  // ==========================================================================
  InsertionHT() = default; //< Default ctor
  InsertionHT(const InsertionHT&) = default; //< Copy ctor
  InsertionHT(InsertionHT&&) = delete; //< Move ctor
  InsertionHT(VisAVis& donor_donor_align, VisAVis& exo_recv_align,
              int32_t length, char* seq, int32_t donor_id);

  virtual Mutation* Clone() const override { return new InsertionHT(*this); };

  // ==========================================================================
  //                                Destructor
  // ==========================================================================
  virtual ~InsertionHT(); //< Destructor

  // ==========================================================================
  //                                Operators
  // ==========================================================================
  /// Copy assignment
  InsertionHT& operator=(const InsertionHT& other) = default;

  /// Move assignment
  InsertionHT& operator=(InsertionHT&& other) = delete;

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
    return INS_HT;
  };

  int32_t donor_pos1() const {
    return donor_donor_align_.get_i_1();
  }

  int32_t donor_pos2() const {
    return donor_donor_align_.get_i_2();
  }

  int32_t exogenote_pos() const {
    return exo_recv_align_.get_i_1();
  }

  int32_t receiver_pos() const {
    return exo_recv_align_.get_i_2();
  }

  AlignmentSense sense() const {
    return exo_recv_align_.get_sense();
  }

  int32_t length() const {
    return length_;
  }

  char* seq() const {
    return seq_;
  }

  int32_t donor_id() const {
    return donor_id_;
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
  VisAVis donor_donor_align_, exo_recv_align_;
  int32_t length_;
  char* seq_ = nullptr;
  int32_t donor_id_ = -1;
};

} // namespace aevol
#endif //AEVOL_INSERTION_HT_H__
