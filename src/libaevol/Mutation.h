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
// ****************************************************************************


#ifndef AEVOL_MUTATION_H__
#define AEVOL_MUTATION_H__


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>




// =================================================================
//                            Project Files
// =================================================================
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "ae_enums.h"
#include <zlib.h>

namespace aevol {

// =================================================================
//                          Class declarations
// =================================================================


enum MutationType
{
  // Simple mutation types.
  SWITCH  = 0,
  S_INS,
  S_DEL,
  DUPL,
  DEL,
  TRANS,
  INV,
  INSERT,
  INS_HT,
  REPL_HT,

  // Composite mutation types follow. They represent categories of
  // several simple mutation types. Therefore, they should not be used
  // as array index for counters.
  //
  // The composite mutations should extend the ae_mutation_type. But
  // C++ enums can't be inherited directly.
  S_MUT, // SWITCH or S_INS or S_DEL
  REARR, // DUPL or DEL or TRANS or INV
  HT,    // INS_HT or REPL_HT
  INDEL  // S_INS or S_DEL
};

// enum composite_mutation_type {
// };

class Mutation
{
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    Mutation( void );
    Mutation( const Mutation &model );
    Mutation( gzFile backup_file );

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~Mutation( void );

    // =================================================================
    //                        Accessors: Getters
    // =================================================================
    inline MutationType get_mut_type( void ) const;
    inline int32_t          get_length( void ) const;

    void get_infos_point_mutation( int32_t* pos ) const;
    void get_infos_small_insertion( int32_t* pos, int32_t* length ) const; // everything except the sequence
    void get_sequence_small_insertion( char* seq ) const; // seq must be a char array, large enough to contain _length+1 characters
    void get_infos_small_deletion( int32_t* pos, int32_t* length ) const;
    void get_infos_duplication( int32_t* pos1, int32_t* pos2, int32_t* pos3, int16_t* align_score = NULL ) const;
    void get_infos_deletion( int32_t* pos1, int32_t* pos2, int16_t* align_score = NULL ) const;
    void get_infos_translocation( int32_t* pos1, int32_t* pos2, int32_t* pos3, int32_t* pos4, bool* invert,
                                  int16_t* align_score_1 = NULL, int16_t* align_score_2 = NULL ) const;
    void get_infos_inversion( int32_t* pos1, int32_t* pos2, int16_t* align_score = NULL ) const;
    void get_infos_insertion( int32_t* pos, int32_t* length ) const;
    void get_sequence_insertion( char* seq ) const;
    void get_infos_ins_HT( int32_t* pos1, int32_t* pos2, int32_t* pos3, int32_t* pos4, ae_sense* sense, int32_t* length ) const;
    void get_sequence_ins_HT( char* seq ) const;
    void get_infos_repl_HT( int32_t* pos1, int32_t* pos2, int32_t* pos3, int32_t* pos4, ae_sense* sense, int32_t* length ) const;
    void get_sequence_repl_HT( char* seq ) const;

    // =================================================================
    //                        Accessors: Setters
    // =================================================================


    // =================================================================
    //                            Public Methods
    // =================================================================
    void report_point_mutation( int32_t pos );
    void report_small_insertion( int32_t pos, int32_t length, const char* seq );
    void report_small_deletion( int32_t pos, int32_t length );
    void report_duplication( int32_t pos_1, int32_t pos_2, int32_t pos_3, int32_t length, int16_t align_score = -1 );
    void report_deletion( int32_t pos_1, int32_t pos_2, int32_t length, int16_t align_score = -1 );
    void report_translocation( int32_t pos_1, int32_t pos_2, int32_t pos_3, int32_t pos_4, int32_t length,
                                bool invert, int16_t align_score_1 = -1, int16_t align_score_2 = -1 );
    void report_inversion( int32_t pos_1, int32_t pos_2, int32_t length, int16_t align_score = -1 );
    void report_insertion( int32_t pos, int32_t length, const char* seq );
    void report_ins_HT(int32_t alignment_1_donor_pos_1, int32_t alignment_1_donor_pos_2, int32_t alignment_2_ind_pos, int32_t alignment_2_donor_pos, int32_t length, int16_t alignment_1_score, int16_t alignment_2_score, int32_t donor_id, ae_sense sense, const char* seq);
    void report_repl_HT(int32_t alignment_1_donor_pos_1, int32_t alignment_1_donor_pos_2, int32_t alignment_2_ind_pos, int32_t alignment_2_donor_pos, int32_t repl_seq_length, int32_t donor_seq_length, int16_t alignment_1_score, int16_t alignment_2_score, int32_t donor_id, ae_sense sense, const char* seq);

    void get_generic_description_string( char * str ) const;

    /* DEPRECATED, use get_length instead */
    int32_t segment_length( int32_t gen_unit_len );

    void save( gzFile backup_file ) const;

    // =================================================================
    //                           Public Attributes
    // =================================================================





  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    //~ Mutation( void )
    //~ {
      //~ printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      //~ exit( EXIT_FAILURE );
    //~ };


    // =================================================================
    //                           Protected Methods
    // =================================================================

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    MutationType  _mut_type;
    int32_t*          _pos;
    int32_t*           _length;
    char*             _seq;
    bool              _invert;
    int16_t*          _align_score;
    int32_t           _donor_id;   // CK: There should also be another int to indicate which genetic unit in the donor
    ae_sense          _sense;
};


// =====================================================================
//                          Accessors' definitions
// =====================================================================
MutationType Mutation::get_mut_type( void ) const
{
  return _mut_type;
}

inline int32_t Mutation::get_length( void ) const
{
  return _length[0];
}



// =====================================================================
//                       Inline functions' definition
// =====================================================================

} // namespace aevol
#endif // AEVOL_MUTATION_H__
