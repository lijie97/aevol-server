//
// Created by dparsons on 03/07/15.
//

// ============================================================================
//                                   Includes
// ============================================================================
#include "ReplacementHT.h"

namespace aevol {

// ============================================================================
//                       Definition of static attributes
// ============================================================================

// ============================================================================
//                                Constructors
// ============================================================================
ReplacementHT::ReplacementHT(const VisAVis& align1, const VisAVis& align2,
                             int32_t length, int32_t replaced_seq_length,
                             char* seq, int32_t donor_id)
    : align1_(align1), align2_(align2),
      length_(length), replaced_seq_length_(replaced_seq_length),
      donor_id_(donor_id) {
  seq_ = new char[length_ + 1];
  memcpy(seq_, seq, length_ + 1);
}

// ============================================================================
//                                 Destructor
// ============================================================================
ReplacementHT::~ReplacementHT() {
  delete seq_;
}

// ============================================================================
//                                   Methods
// ============================================================================
void ReplacementHT::save(gzFile backup_file) const {
  Utils::ExitWithDevMsg("Not implemented yet", __FILE__, __LINE__);
}

void ReplacementHT::load(gzFile backup_file) {
  Utils::ExitWithDevMsg("Not implemented yet", __FILE__, __LINE__);
}

void ReplacementHT::get_generic_description_string(char* str) const {
  sprintf(str, "%" PRId8 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32
      " %" PRId8 " %" PRId16 " %" PRId16 " %" PRId32 " %" PRId32 " ",
          get_mut_type(),
          donor_pos1(), donor_pos2(),
          receiver_pos1(), receiver_pos2(),
          sense(),
          align1_.get_score(), align2_.get_score(),
          length_, replaced_seq_length_);
}

// ============================================================================
//                            Non inline accessors
// ============================================================================

} // namespace aevol
