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


#ifndef __AE_VIS_A_VIS_H__
#define __AE_VIS_A_VIS_H__


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
#include "ae_dna.h"
#include "ae_utils.h"


namespace aevol {


// =================================================================
//                          Class declarations
// =================================================================







class ae_vis_a_vis
{
  friend class ae_align;
  
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    ae_vis_a_vis( const ae_dna* chrom_1, const ae_dna* chrom_2, int32_t i_1, int32_t i_2, ae_sense sense = DIRECT );
    ae_vis_a_vis( const ae_vis_a_vis & orig );
    //~ ae_vis_a_vis( const ae_vis_a_vis * const orig );

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_vis_a_vis( void );

    // =================================================================
    //                              Accessors
    // =================================================================
    inline const ae_dna*  get_chrom_1( void ) const;
    inline const ae_dna*  get_chrom_2( void ) const;
    inline int32_t        get_i_1( void ) const;
    inline int32_t        get_i_2( void ) const;
    inline int16_t        get_score( void ) const;
    inline ae_sense       get_sense( void ) const;

    // =================================================================
    //                              Operators
    // =================================================================
    inline bool operator <  ( ae_vis_a_vis &cmp );
    inline bool operator <= ( ae_vis_a_vis &cmp );    
    inline bool operator >  ( ae_vis_a_vis &cmp );
    inline bool operator >= ( ae_vis_a_vis &cmp );

    // =================================================================
    //                            Public Methods
    // =================================================================
    inline bool match( void );
    inline void step_fwd( void );
    inline void step_back( void );
    inline void add( int common_inc );
    inline void add( int inc_1, int inc_2 );
    inline void sub( int common_inc );
    inline void sub( int inc_1, int inc_2 );
    inline void swap( void );
    
    inline void copy( ae_vis_a_vis * source );
    inline void check_indices( void );

    // =================================================================
    //                           Public Attributes
    // =================================================================





  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
  
    //~ ae_vis_a_vis( void )
    //~ {
      //~ _chrom_1 = _chrom_2 = new ae_dna(); // const members must be initialized
      //~ printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      //~ exit( EXIT_FAILURE );
    //~ };
    /*ae_vis_a_vis( const ae_vis_a_vis &model )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
      };*/


    // =================================================================
    //                           Protected Methods
    // =================================================================

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    const ae_dna *  _chrom_1;
    const ae_dna *  _chrom_2;
    int32_t         _i_1; // Index on chrom_1
    int32_t         _i_2; // Index on chrom_2
    int16_t         _score;
    ae_sense        _sense; // Sense (DIRECT or INDIRECT) of the vis_a_vis (alignement)
    // Say we have the following sequences :
    //    0 1 2 3 4 5 6 7 8 9             0 1 2 3 4 5 6 7 8 9
    //    |a|b|c|d|e|f|g|h|i|j|           |a|b|c|d|e|f|g|h|i|j|
    //
    // The DIRECT vis_a_vis between _i_1 = 3 and _i_2 = 7 is 'd' with 'h' (caracters at indices 3 and 7 resp.).
    //   The corresponding alignment would be "defgh" with "hijkl"
    //
    // WARNING! The INDIRECT vis_a_vis between the same _i_1 = 3 and i_2 = 7 is 'd' with 'g' (and not 'h'!).
    // This is because we are reading backwards (towards the left). Directly left to index 7 is 'g' which corresponds to index 6.
    //   The corresponding alignment would hence be "defgh" with "!g!f!e!d!c" ("!x" means "complementary of x")
};


// =====================================================================
//                          Accessors' definitions
// =====================================================================
inline const ae_dna* ae_vis_a_vis::get_chrom_1( void ) const
{
  return _chrom_1;
}

inline const ae_dna* ae_vis_a_vis::get_chrom_2( void ) const
{
  return _chrom_2;
}

inline int32_t ae_vis_a_vis::get_i_1( void ) const
{
  return _i_1;
}

inline int32_t ae_vis_a_vis::get_i_2( void ) const
{
  return _i_2;
}

inline int16_t ae_vis_a_vis::get_score( void ) const
{
  return _score;
}

inline ae_sense ae_vis_a_vis::get_sense( void ) const
{
  return _sense;
}



// =====================================================================
//                          Operators' definitions
// =====================================================================
inline bool ae_vis_a_vis::operator < ( ae_vis_a_vis &cmp )
{
  return ( _i_1 < cmp._i_1 );
}

inline bool ae_vis_a_vis::operator <= ( ae_vis_a_vis &cmp )
{
  return ( _i_1 <= cmp._i_1 );
}

inline bool ae_vis_a_vis::operator > ( ae_vis_a_vis &cmp )
{
  return ( _i_1 > cmp._i_1 );
}

inline bool ae_vis_a_vis::operator >= ( ae_vis_a_vis &cmp )
{
  return ( _i_1 >= cmp._i_1 );
}


// =====================================================================
//                       Inline functions' definition
// =====================================================================
inline bool ae_vis_a_vis::match( void )
{
  if ( _sense == DIRECT )
  {
    return ( _chrom_1->get_data()[ae_utils::mod(_i_1, _chrom_1->get_length())] == _chrom_2->get_data()[ae_utils::mod(_i_2, _chrom_2->get_length())] );
  }
  else // ( _sense == INDIRECT )
  {
    // Note that we are reading the sequence backwards, The nucleotide corresponding to a breakpoint at point <i>
    // is hence stored at index <i-1>
    //    a b c d e f g h i j
    //    |_|_|_|_|_|_|_|_|_|_|
    //    | | | | | | | | | | |
    //      9 8 7 6 5 4 3 2 1 0
    //
    // The breakpoint F-5 puts into a vis_a_vis the nucleotide at index F on seq1 and that at index 4 (not 5!!!) on seq2
    return ( _chrom_1->get_data()[ae_utils::mod(_i_1, _chrom_1->get_length())] != _chrom_2->get_data()[ae_utils::mod(_i_2-1, _chrom_2->get_length())] );
  }
}

inline void ae_vis_a_vis::step_fwd( void )
{
  if ( _sense == DIRECT )
  {
    _i_1++;
    _i_2++;
  }
  else // ( _sense == INDIRECT )
  {
    _i_1++;
    _i_2--;
  }
}

inline void ae_vis_a_vis::step_back( void )
{
  if ( _sense == DIRECT )
  {
    _i_1--;
    _i_2--;
  }
  else // ( _sense == INDIRECT )
  {
    _i_1--;
    _i_2++;
  }
}

inline void ae_vis_a_vis::add( int common_inc )
{
  if ( _sense == DIRECT )
  {
    _i_1 += common_inc;
    _i_2 += common_inc;
  }
  else // ( _sense == INDIRECT )
  {
    _i_1 += common_inc;
    _i_2 -= common_inc;
  }
}

inline void ae_vis_a_vis::add( int inc_1, int inc_2 )
{
  if ( _sense == DIRECT )
  {
    _i_1 += inc_1;
    _i_2 += inc_2;
  }
  else // ( _sense == INDIRECT )
  {
    _i_1 += inc_1;
    _i_2 -= inc_2;
  }
}

inline void ae_vis_a_vis::sub( int common_inc )
{
  if ( _sense == DIRECT )
  {
    _i_1 -= common_inc;
    _i_2 -= common_inc;
  }
  else // ( _sense == INDIRECT )
  {
    _i_1 -= common_inc;
    _i_2 += common_inc;
  }
}

inline void ae_vis_a_vis::sub( int inc_1, int inc_2 )
{
  if ( _sense == DIRECT )
  {
    _i_1 -= inc_1;
    _i_2 -= inc_2;
  }
  else // ( _sense == INDIRECT )
  {
    _i_1 -= inc_1;
    _i_2 += inc_2;
  }
}

inline void ae_vis_a_vis::swap( void )
{
  const ae_dna *  tmp_chrom = _chrom_1;
  int32_t         tmp_i     = _i_1;
  
  _chrom_1  = _chrom_2;
  _i_1      = _i_2;
  
  _chrom_2  = tmp_chrom;
  _i_2      = tmp_i;
}

inline void ae_vis_a_vis::copy( ae_vis_a_vis * source )
{
  _i_1 = source->_i_1;
  _i_2 = source->_i_2;
  _chrom_1 = source->_chrom_1;
  _chrom_2 = source->_chrom_2;
  _sense = source->_sense;
  _score = source->_score;
}

inline void ae_vis_a_vis::check_indices( void )
{
  _i_1 = ae_utils::mod( _i_1, _chrom_1->get_length() );
  _i_2 = ae_utils::mod( _i_2, _chrom_2->get_length() );
}
} // namespace aevol

#endif // __AE_VIS_A_VIS_H__
