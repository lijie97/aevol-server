//
// Created by dparsons on 03/07/15.
//

#ifndef AEVOL_REPLACEMENT_HT_H__
#define AEVOL_REPLACEMENT_HT_H__


// ============================================================================
//                                   Includes
// ============================================================================

#include "HT.h"
#include "VisAVis.h"

namespace aevol {

/**
 *
 */
class ReplacementHT : public HT {

 public :
  // ==========================================================================
  //                               Constructors
  // ==========================================================================
  ReplacementHT() = default; //< Default ctor
  ReplacementHT(const ReplacementHT&) = default; //< Copy ctor
  ReplacementHT(ReplacementHT&&) = delete; //< Move ctor
  ReplacementHT(const VisAVis& align1, const VisAVis& align2,
                int32_t length, int32_t replaced_seq_length,
                char* seq, int32_t donor_id);

  virtual Mutation* Clone() const override { return new ReplacementHT(*this); };

  // ==========================================================================
  //                                Destructor
  // ==========================================================================
  virtual ~ReplacementHT(); //< Destructor

  // ==========================================================================
  //                                Operators
  // ==========================================================================
  /// Copy assignment
  ReplacementHT& operator=(const ReplacementHT& other) = default;

  /// Move assignment
  ReplacementHT& operator=(ReplacementHT&& other) = delete;

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
    return REPL_HT;
  };

  int32_t donor_pos1() const {
    return align1_.get_i_2();
  }

  int32_t donor_pos2() const {
    return align2_.get_i_2();
  }

  int32_t receiver_pos1() const {
    return align1_.get_i_1();
  }

  int32_t receiver_pos2() const {
    return align2_.get_i_1();
  }

  AlignmentSense sense() const {
    return align2_.get_sense();
  }

  char* seq() const {
    return seq_;
  }

  int32_t length() const {
    return length_;
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
  /**
   * Alignments between the receiver and the donor
   *
   * The first point of each alignment corresponds to the receiver, the sequence
   * that will be replaced lies between align1_->i1 and align2_->i1. It will be
   * replaced by the seq. btw align1_->i2 and align2_->i2 on the donor.
   */
  VisAVis align1_, align2_;
  int32_t length_;
  int32_t replaced_seq_length_;
  char* seq_ = nullptr;
  int32_t donor_id_ = -1;
};

} // namespace aevol
#endif //AEVOL_REPLACEMENT_HT_H__
