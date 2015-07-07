//
// Created by dparsons on 01/07/15.
//

// ============================================================================
//                                   Includes
// ============================================================================
#include "PointMutation.h"


namespace aevol {

// ============================================================================
//                       Definition of static attributes
// ============================================================================

// ============================================================================
//                                Constructors
// ============================================================================
PointMutation::PointMutation(int32_t pos) : pos_(pos) {}

// ============================================================================
//                                 Destructor
// ============================================================================

// ============================================================================
//                                   Methods
// ============================================================================
void PointMutation::save(gzFile backup_file) const {
  int8_t tmp_mut_type = SWITCH;
  gzwrite(backup_file, &tmp_mut_type, sizeof(tmp_mut_type));
  gzwrite(backup_file, &pos_, sizeof(pos_));
}

void PointMutation::load(gzFile backup_file) {
  gzread(backup_file, &pos_, sizeof(pos_));
}

void PointMutation::get_generic_description_string(char* str) const {
  sprintf(str, "%" PRId8 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32
      " %" PRId8 " %" PRId16 " %" PRId16 " %" PRId32 " %" PRId32,
          get_mut_type(), pos_, -1, -1, -1, -1, -1, -1, -1, -1);
}

// ============================================================================
//                            Non inline accessors
// ============================================================================

} // namespace aevol
