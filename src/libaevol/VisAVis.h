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


#ifndef AEVOL_VIS_A_VIS_H__
#define AEVOL_VIS_A_VIS_H__


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
#include "Dna.h"
#include "Utils.h"


namespace aevol {


// =================================================================
//                          Class declarations
// =================================================================







class VisAVis
{
  friend class Alignment;
  
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    VisAVis() = default;
    VisAVis(const Dna * chrom_1, const Dna * chrom_2,
            int32_t i_1, int32_t i_2, AlignmentSense sense = DIRECT);
    VisAVis(const VisAVis & orig);

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~VisAVis(void);

    // =================================================================
    //                              Accessors
    // =================================================================
    inline const Dna*  get_chrom_1( void ) const;
    inline const Dna*  get_chrom_2( void ) const;
    inline int32_t        get_i_1( void ) const;
    inline int32_t        get_i_2( void ) const;
    inline int16_t        get_score( void ) const;
    inline AlignmentSense get_sense( void ) const;

    // =================================================================
    //                              Operators
    // =================================================================
    inline bool operator <  ( VisAVis &cmp );
    inline bool operator <= ( VisAVis &cmp );
    inline bool operator >  ( VisAVis &cmp );
    inline bool operator >= ( VisAVis &cmp );

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
    
    inline void copy( VisAVis * source );
    inline void check_indices( void );

    // =================================================================
    //                           Public Attributes
    // =================================================================





  protected :
    // =================================================================
    //                           Protected Methods
    // =================================================================

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    const Dna* _chrom_1 = nullptr;
    const Dna* _chrom_2 = nullptr;
    int32_t _i_1 = -1; //< Index on chrom_1
    int32_t _i_2 = -1; //< Index on chrom_2
    int16_t _score = -1;
    // Sense (DIRECT or INDIRECT) of the vis_a_vis (alignement)
    AlignmentSense _sense = DIRECT;
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
inline const Dna *VisAVis::get_chrom_1( void ) const
{
  return _chrom_1;
}

inline const Dna *VisAVis::get_chrom_2( void ) const
{
  return _chrom_2;
}

inline int32_t VisAVis::get_i_1( void ) const
{
  return _i_1;
}

inline int32_t VisAVis::get_i_2( void ) const
{
  return _i_2;
}

inline int16_t VisAVis::get_score( void ) const
{
  return _score;
}

inline AlignmentSense VisAVis::get_sense( void ) const
{
  return _sense;
}



// =====================================================================
//                          Operators' definitions
// =====================================================================
inline bool VisAVis::operator < ( VisAVis &cmp )
{
  return ( _i_1 < cmp._i_1 );
}

inline bool VisAVis::operator <= ( VisAVis &cmp )
{
  return ( _i_1 <= cmp._i_1 );
}

inline bool VisAVis::operator > ( VisAVis &cmp )
{
  return ( _i_1 > cmp._i_1 );
}

inline bool VisAVis::operator >= ( VisAVis &cmp )
{
  return ( _i_1 >= cmp._i_1 );
}


// =====================================================================
//                       Inline functions' definition
// =====================================================================
inline bool VisAVis::match( void )
{
  if ( _sense == DIRECT )
  {
    return (_chrom_1->data()[Utils::mod(_i_1, _chrom_1->length())] ==
        _chrom_2->data()[Utils::mod(_i_2, _chrom_2->length())] );
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
    return (_chrom_1->data()[Utils::mod(_i_1, _chrom_1->length())] !=
        _chrom_2->data()[Utils::mod(_i_2-1, _chrom_2->length())] );
  }
}

inline void VisAVis::step_fwd( void )
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

inline void VisAVis::step_back( void )
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

inline void VisAVis::add( int common_inc )
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

inline void VisAVis::add( int inc_1, int inc_2 )
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

inline void VisAVis::sub( int common_inc )
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

inline void VisAVis::sub( int inc_1, int inc_2 )
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

inline void VisAVis::swap( void )
{
  const Dna *  tmp_chrom = _chrom_1;
  int32_t         tmp_i     = _i_1;
  
  _chrom_1  = _chrom_2;
  _i_1      = _i_2;
  
  _chrom_2  = tmp_chrom;
  _i_2      = tmp_i;
}

inline void VisAVis::copy( VisAVis * source )
{
  _i_1 = source->_i_1;
  _i_2 = source->_i_2;
  _chrom_1 = source->_chrom_1;
  _chrom_2 = source->_chrom_2;
  _sense = source->_sense;
  _score = source->_score;
}

inline void VisAVis::check_indices( void )
{
  _i_1 = Utils::mod( _i_1, _chrom_1->length() );
  _i_2 = Utils::mod( _i_2, _chrom_2->length() );
}
} // namespace aevol

#endif // AEVOL_VIS_A_VIS_H__
