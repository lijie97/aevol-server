//
// Created by dparsons on 02/07/15.
//

// ============================================================================
//                                   Includes
// ============================================================================
#include "SmallDeletion.h"

namespace aevol {

// ============================================================================
//                       Definition of static attributes
// ============================================================================

// ============================================================================
//                                Constructors
// ============================================================================
SmallDeletion::SmallDeletion(int32_t pos, int16_t length) :
    pos_(pos),
    length_(length) {
}

// ============================================================================
//                                 Destructor
// ============================================================================

// ============================================================================
//                                   Methods
// ============================================================================
void SmallDeletion::save(gzFile backup_file) const {
  int8_t tmp_mut_type = S_DEL;
  gzwrite(backup_file, &tmp_mut_type, sizeof(tmp_mut_type));
  gzwrite(backup_file, &pos_, sizeof(pos_));
  gzwrite(backup_file, &length_, sizeof(length_));
}

void SmallDeletion::load(gzFile backup_file) {
  gzread(backup_file, &pos_, sizeof(pos_));
  gzread(backup_file, &length_, sizeof(length_));
}

void SmallDeletion::get_generic_description_string(char* str) const {
  sprintf(str, "%" PRId8 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32
      " %" PRId8 " %" PRId16 " %" PRId16 " %" PRId32 " %" PRId32,
          get_mut_type(), pos_, -1, -1, -1, -1, -1, -1, length_, -1);
}

// ============================================================================
//                            Non inline accessors
// ============================================================================

} // namespace aevol
