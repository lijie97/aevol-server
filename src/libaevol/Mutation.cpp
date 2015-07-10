// ****************************************************************************
//
//          Aevol - An in silico experimental evolution platform
//
// ****************************************************************************
//
// Copyright: See the AUTHORS file provided with the package or <www.aevol.fr>
// Web: http://www.aevol.fr/
// E-mail: See <http://www.aevol.fr/contact/>
// Original Authors : Guillaume Beslon, Carole Knibbe, David Parsons
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//*****************************************************************************




// =================================================================
//                              Includes
// =================================================================
#include <cstdlib>

#include "Mutation.h"
#include "PointMutation.h"
#include "SmallInsertion.h"
#include "SmallDeletion.h"
#include "Duplication.h"
#include "Deletion.h"
#include "Translocation.h"
#include "Inversion.h"
#include "InsertionHT.h"
#include "ReplacementHT.h"

namespace aevol {


//#############################################################################
//
//                                Class Mutation
//
//#############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
Mutation* Mutation::Load(gzFile backup_file) {
  // Retrieve mutation type
  int8_t tmp_mut_type;
  gzread(backup_file, &tmp_mut_type,  sizeof(tmp_mut_type));
  MutationType mut_type = (MutationType) tmp_mut_type;

  // Call the appropriate constructor accordingly
  Mutation* mut;
  switch (mut_type) {
    case SWITCH :
      mut = new PointMutation();
      break;
    case S_INS :
      mut = new SmallInsertion();
      break;
    case S_DEL :
      mut = new SmallDeletion();
      break;
    case DUPL :
      mut = new Duplication();
      break;
    case DEL :
      mut = new Deletion();
      break;
    case TRANS :
      mut = new Translocation();
      break;
    case INV :
      mut = new Inversion();
      break;
    case INS_HT :
      mut = new InsertionHT();
      break;
    case REPL_HT :
      mut = new ReplacementHT();
      break;
    default :
      Utils::ExitWithDevMsg("invalid mutation type ", __FILE__, __LINE__);
  }

  // Load from backup file
  mut->load(backup_file);
  return mut;
}




//Mutation::Mutation(gzFile backup_file)
//{
//  _pos = NULL;
//  _length = NULL;
//  _seq = NULL;
//  _align_score = NULL;
//
//  int8_t tmp_mut_type;
//  gzread(backup_file, &tmp_mut_type,  sizeof(tmp_mut_type));
//  _mut_type = (MutationType) tmp_mut_type;
//  //~ printf("mut type %d\n", _mut_type);
//
//  switch (_mut_type)
//  {
//    case SWITCH :
//    {
//      _pos = new int32_t;
//      gzread(backup_file, _pos,  sizeof(*_pos));
//      break;
//    }
//    case S_INS :
//    {
//      _pos = new int32_t;
//      gzread(backup_file, _pos,      sizeof(*_pos));
//      _length = new int32_t;
//      gzread(backup_file, _length,  sizeof(*_length));
//
//      _seq = new char[_length[0] + 1];
//      gzread(backup_file, _seq,  _length[0] * sizeof(_seq[0]));
//      _seq[_length[0]] = '\0';
//      break;
//    }
//    case S_DEL :
//    {
//      _pos = new int32_t;
//      gzread(backup_file, _pos,      sizeof(*_pos));
//      _length = new int32_t;
//      gzread(backup_file, _length,  sizeof(*_length));
//      break;
//    }
//    case DUPL :
//    {
//      _pos = new int32_t[3];
//      gzread(backup_file, _pos,  3 * sizeof(_pos[0]));
//      _length = new int32_t;
//      gzread(backup_file, _length,  sizeof(*_length));
//      _align_score = new int16_t;
//      gzread(backup_file, _align_score, sizeof(*_align_score));
//
//      break;
//    }
//    case DEL :
//    {
//      _pos = new int32_t[2];
//      gzread(backup_file, _pos,  2 * sizeof(_pos[0]));
//      _length = new int32_t;
//      gzread(backup_file, _length,  sizeof(*_length));
//      _align_score = new int16_t;
//      gzread(backup_file, _align_score, sizeof(*_align_score));
//
//      break;
//    }
//    case TRANS :
//    {
//      _pos = new int32_t[4];
//      gzread(backup_file, _pos,  4 * sizeof(_pos[0]));
//      int8_t tmp_invert;
//      gzread(backup_file, &tmp_invert,  sizeof(tmp_invert));
//      _invert = (tmp_invert != 0);
//      _length = new int32_t;
//      gzread(backup_file, _length,  sizeof(*_length));
//      _align_score = new int16_t[2];
//      gzread(backup_file, _align_score, 2 * sizeof(_align_score[0]));
//
//      break;
//    }
//    case INV :
//    {
//      _pos = new int32_t[2];
//      gzread(backup_file, _pos,  2 * sizeof(_pos[0]));
//      _length = new int32_t;
//      gzread(backup_file, _length,  sizeof(*_length));
//      _align_score = new int16_t;
//      gzread(backup_file, _align_score, sizeof(*_align_score));
//
//      break;
//    }
//    case INSERT:
//    {
//      _pos = new int32_t;
//      gzread(backup_file, _pos,  sizeof(*_pos));
//     _length = new int32_t;
//      gzread(backup_file, _length,  sizeof(*_length));
//      _seq = new char[_length[0] + 1];
//      gzread(backup_file, _seq,  _length[0] * sizeof(_seq[0]));
//      _seq[_length[0]] = '\0';
//    }
//    case INS_HT:
//    {
//      _pos = new int32_t[4];
//      gzread(backup_file, _pos,  4 * sizeof(_pos[0]));
//      _length = new int32_t;
//      gzread(backup_file, _length,  sizeof(*_length));
//      _seq = new char[_length[0] + 1];
//      gzread(backup_file, _seq,  _length[0] * sizeof(_seq[0]));
//      _seq[_length[0]] = '\0';
//      _align_score = new int16_t[2];
//      gzread(backup_file, _align_score, 2 * sizeof(_align_score[0]));
//      gzread(backup_file, &_donor_id,  sizeof(_donor_id));
//      gzread(backup_file, &_sense,  sizeof(_sense));
//      break;
//    }
//    case REPL_HT:
//    {
//      _pos = new int32_t[4];
//      gzread(backup_file, _pos,  4 * sizeof(_pos[0]));
//      _length = new int32_t[2];
//      gzread(backup_file, _length, 2 *sizeof(_length[0]));
//      _seq = new char[_length[1] + 1];
//      gzread(backup_file, _seq,  _length[1] * sizeof(_seq[0]));
//      _seq[_length[1]] = '\0';
//      _align_score = new int16_t[2];
//      gzread(backup_file, _align_score, 2 * sizeof(_align_score[0]));
//      gzread(backup_file, &_donor_id,  sizeof(_donor_id));
//      gzread(backup_file, &_sense,  sizeof(_sense));
//      break;
//    }
//    default :
//    {
//      fprintf(stderr, "ERROR, invalid mutation type \"%d\" in file %s:%d.\n", _mut_type, __FILE__, __LINE__);
//      exit(EXIT_FAILURE);
//      break;
//    }
//  }
//
//
//}





// =================================================================
//                             Destructor
// =================================================================

// =================================================================
//                            Public Methods
// =================================================================

// =================================================================
//                           Protected Methods
// =================================================================
} // namespace aevol
