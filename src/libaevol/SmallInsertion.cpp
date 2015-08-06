//
// Created by dparsons on 01/07/15.
//

// ============================================================================
//                                   Includes
// ============================================================================
#include "SmallInsertion.h"

#include <cstring>


namespace aevol {

// ============================================================================
//                       Definition of static attributes
// ============================================================================

// ============================================================================
//                                Constructors
// ============================================================================
SmallInsertion::SmallInsertion(const SmallInsertion& other) :
    pos_(other.pos_),
    length_(other.length_) {
  seq_ = new char[length_ + 1];
  memcpy(seq_, other.seq_, length_ + 1);
}

SmallInsertion::SmallInsertion(int32_t pos, int32_t length, const char* seq) :
    pos_(pos),
    length_(length) {
  seq_ = new char[length_ + 1];
  memcpy(seq_, seq, length_);
  seq_[length_] = '\0';
}

// ============================================================================
//                                 Destructor
// ============================================================================
SmallInsertion::~SmallInsertion() {
  delete [] seq_;
}

// ============================================================================
//                                 Operators
// ============================================================================
/// Copy assignment
SmallInsertion& SmallInsertion::operator=(const SmallInsertion& other) {
  pos_ = other.pos_;
  length_ = other.length_;
  seq_ = new char[length_ + 1];
  memcpy(seq_, other.seq_, length_ + 1);

  return *this;
}

// ============================================================================
//                                   Methods
// ============================================================================
void SmallInsertion::save(gzFile backup_file) const {
  int8_t tmp_mut_type = S_INS;
  gzwrite(backup_file, &tmp_mut_type, sizeof(tmp_mut_type));
  gzwrite(backup_file, &pos_, sizeof(pos_));
  gzwrite(backup_file, &length_, sizeof(length_));
  gzwrite(backup_file, seq_, length_ * sizeof(seq_[0]));
}

void SmallInsertion::load(gzFile backup_file) {
  gzread(backup_file, &pos_, sizeof(pos_));
  gzread(backup_file, &length_, sizeof(length_));
  seq_ = new char[length_ + 1];
  gzread(backup_file, seq_, length_ * sizeof(seq_[0]));
}

void SmallInsertion::get_generic_description_string(char* str) const {
  sprintf(str, "%" PRId8 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32
      " %" PRId8 " %" PRId16 " %" PRId16 " %" PRId32 " %" PRId32,
          get_mut_type(), pos_, -1, -1, -1, -1, -1, -1, length_, -1);
}

// ============================================================================
//                            Non inline accessors
// ============================================================================

} // namespace aevol
