//
// Created by dparsons on 03/07/15.
//

// ============================================================================
//                                   Includes
// ============================================================================
#include "InsertionHT.h"

namespace aevol {

// ============================================================================
//                       Definition of static attributes
// ============================================================================

// ============================================================================
//                                Constructors
// ============================================================================
InsertionHT::InsertionHT(VisAVis& donor_donor_align, VisAVis& exo_recv_align,
                         int32_t length, char* seq, int32_t donor_id)
    : donor_donor_align_(donor_donor_align), exo_recv_align_(exo_recv_align),
      length_(length), donor_id_(donor_id) {
  assert(donor_donor_align_.get_sense() == DIRECT);
  seq_ = new char[length_ + 1];
  memcpy(seq_, seq, length_ + 1);
}

// ============================================================================
//                                 Destructor
// ============================================================================
InsertionHT::~InsertionHT() {
  delete seq_;
}

// ============================================================================
//                                   Methods
// ============================================================================
void InsertionHT::save(gzFile backup_file) const {
  Utils::ExitWithDevMsg("Not implemented yet", __FILE__, __LINE__);
}

void InsertionHT::load(gzFile backup_file) {
  Utils::ExitWithDevMsg("Not implemented yet", __FILE__, __LINE__);
}

void InsertionHT::get_generic_description_string(char* str) const {
  sprintf(str, "%" PRId8 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32
      " %" PRId8 " %" PRId16 " %" PRId16 " %" PRId32 " %" PRId32,
          get_mut_type(),
          donor_pos1(), donor_pos2(), exogenote_pos(), receiver_pos(),
          sense(),
          donor_donor_align_.get_score(), exo_recv_align_.get_score(),
          length_, (int32_t) -1);
}

// ============================================================================
//                            Non inline accessors
// ============================================================================

} // namespace aevol
