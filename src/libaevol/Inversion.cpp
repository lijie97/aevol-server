//
// Created by dparsons on 02/07/15.
//

// ============================================================================
//                                   Includes
// ============================================================================
#include "Inversion.h"

namespace aevol {

// ============================================================================
//                       Definition of static attributes
// ============================================================================

// ============================================================================
//                                Constructors
// ============================================================================
Inversion::Inversion(int32_t pos1, int32_t pos2, int32_t length,
                     int16_t align_score) :
    pos1_(pos1), pos2_(pos2),
    length_(length), align_score_(align_score) {
}

// ============================================================================
//                                 Destructor
// ============================================================================

// ============================================================================
//                                   Methods
// ============================================================================
void Inversion::save(gzFile backup_file) const {
  int8_t tmp_mut_type = INV;
  gzwrite(backup_file, &tmp_mut_type, sizeof(tmp_mut_type));
  gzwrite(backup_file, &pos1_, sizeof(pos1_));
  gzwrite(backup_file, &pos2_, sizeof(pos2_));
  gzwrite(backup_file, &length_, sizeof(length_));
  gzwrite(backup_file, &align_score_, sizeof(align_score_));
}

void Inversion::load(gzFile backup_file) {
  gzread(backup_file, &pos1_, sizeof(pos1_));
  gzread(backup_file, &pos2_, sizeof(pos2_));
  gzread(backup_file, &length_, sizeof(length_));
  gzread(backup_file, &align_score_, sizeof(align_score_));
}

void Inversion::get_generic_description_string(char* str) const {
  sprintf(str, "%" PRId8 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32
      " %" PRId8 " %" PRId16 " %" PRId16 " %" PRId32 " %" PRId32,
          get_mut_type(), pos1(), pos2(), -1, -1,
          -1, align_score_, -1, length_, -1);
}

// ============================================================================
//                            Non inline accessors
// ============================================================================

} // namespace aevol
