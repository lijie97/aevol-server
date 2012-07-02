//*****************************************************************************
//
//                         aevol - Artificial Evolution
//
// Copyright (C) 2004  LIRIS.
// Web: https://liris.cnrs.fr/
// E-mail: carole.knibbe@liris.cnrs.fr
// Original Authors : Guillaume Beslon, Carole Knibbe, Virginie Lefort
//                    David Parsons
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//*****************************************************************************


/** \class
 *  \brief
 */


#ifndef __AE_ALIGN_H__
#define __AE_ALIGN_H__


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_object.h>
#include <ae_vis_a_vis.h>




// =================================================================
//                          Class declarations
// =================================================================
class ae_dna;






class ae_align : public ae_object
{
  public :

    // =================================================================
    //                             Constructors
    // =================================================================

    // =================================================================
    //                             Destructors
    // =================================================================

    // =================================================================
    //                              Accessors
    // =================================================================

    // =================================================================
    //                            Public Methods
    // =================================================================
    static ae_vis_a_vis* search_alignment_direct( const ae_dna* chrom1, const int32_t seed1,
                                                  const ae_dna* chrom2, const int32_t seed2, int16_t needed_score );
    static ae_vis_a_vis* search_alignment_indirect( const ae_dna* chrom1, const int32_t seed1,
                                                    const ae_dna* chrom2, const int32_t seed2, int16_t needed_score );

    // =================================================================
    //                           Public Attributes
    // =================================================================

    static bool with_alignments;

    static ae_align_fun_shape align_fun_shape;
    
    static double  align_sigm_lambda;
    static int16_t align_sigm_mean;
    static int16_t align_lin_min;
    static int16_t align_lin_max;
    
    // Maximum shift of one seq on the other
    static int16_t align_max_shift;
    // Work zone half length
    static int16_t align_w_zone_h_len;
    // Corresponding residues match bonus
    static int16_t align_match_bonus;
    // Corresponding residues mismatch cost
    static int16_t align_mismatch_cost;





  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    ae_align( void )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
    ae_align( const ae_align &model )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };


    // =================================================================
    //                           Protected Methods
    // =================================================================

    // =================================================================
    //                          Protected Attributes
    // =================================================================
};


// =====================================================================
//                          Accessors' definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================


#endif // __AE_ALIGN_H__
