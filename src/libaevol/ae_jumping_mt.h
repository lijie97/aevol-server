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


/*! \class
    \brief
*/


#ifndef __AE_JUMPING_MT_H__
#define __AE_JUMPING_MT_H__


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>
#include <zlib.h>
#include <math.h>
#include <SFMT-src-1.4/SFMT.h>
#include <SFMT-src-1.4/jump/SFMT-jump.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_object.h>




// =================================================================
//                          Class declarations
// =================================================================




// MT_RAND_MAX = 2^32-1
#define MT_RAND_MAX         4294967295.0
#define MT_RAND_MAX_PLUS_1  4294967296.0





class ae_jumping_mt : public ae_object
{
  public :
    // =================================================================
    //                             Constructors
    // =================================================================
    ae_jumping_mt( const uint32_t& simple_seed );   // Initialize with a simple uint32_t
    ae_jumping_mt( const ae_jumping_mt &model );    // Create a copy of an existing generator
    ae_jumping_mt( gzFile* backup_file );           // Load from a gz backup file

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_jumping_mt( void );

    // =================================================================
    //                        Accessors: getters
    // =================================================================

    // =================================================================
    //                        Accessors: setters
    // =================================================================

    // =================================================================
    //                              Operators
    // =================================================================

    // =================================================================
    //                            Public Methods
    // =================================================================
    inline double   random( void );         // Double in [0, 1[ (uniform distribution)
    inline int8_t   random( int8_t max );   // ~
    inline int16_t  random( int16_t max );  // ~
    inline int32_t  random( int32_t max );  // ~ > Integer in [0, max[ (uniform distribution)
    inline int64_t  random( int64_t max );  // ~
    int32_t         binomial_random( int32_t nb, double prob ); // Binomial drawing of parameters (nb, prob)
    double          gaussian_random( void );                    // Double following a Standard Normal distribution
    void            multinomial_drawing ( int32_t* destination, double* source, int32_t nb_drawings, int32_t colors );
    // Multinomial drawing of parameters ( nb, {source[0], source[1], ... source[colors-1]} )
    
    void jump( void );
    
    void save( gzFile* backup_file ) const;

    // =================================================================
    //                           Public Attributes
    // =================================================================
    static int32_t nb_jumps;
    static double  jump_time;





  protected :
    
    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    ae_jumping_mt( void )
    {
      printf( "%s:%d: error: call to forbidden constructor.\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
    /*ae_jumping_mt( const ae_jumping_mt &model )
    {
      printf( "%s:%d: error: call to forbidden constructor.\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };*/


    // =================================================================
    //                           Protected Methods
    // =================================================================
    static double gammln( double X );

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    sfmt_t* _sfmt;
};


// =====================================================================
//                           Getters' definitions
// =====================================================================

// =====================================================================
//                           Setters' definitions
// =====================================================================

// =====================================================================
//                          Operators' definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================
/*!
  Draw a double precision real-number in [0, 1) with a uniform distribution
 */
inline double ae_jumping_mt::random( void )
{
  return sfmt_genrand_real2( _sfmt );
}

/*!
  Draw an 8-bit integer in [0, max[ with a uniform distribution
 */
inline int8_t ae_jumping_mt::random( int8_t max )
{
  return (int8_t) floor( ((double)max) * sfmt_genrand_real2(_sfmt) );
}

/*!
  Draw an 16-bit integer in [0, max[ with a uniform distribution
 */
inline int16_t ae_jumping_mt::random( int16_t max )
{
  return (int16_t) floor( ((double)max) * sfmt_genrand_real2(_sfmt) );
}

/*!
  Draw an 32-bit integer in [0, max[ with a uniform distribution
 */
inline int32_t ae_jumping_mt::random( int32_t max )
{
  return (int32_t) floor( ((double)max) * sfmt_genrand_real2(_sfmt) );
}

/*!
  Draw an 64-bit integer in [0, max[ with a uniform distribution
 */
inline int64_t ae_jumping_mt::random( int64_t max )
{
  return (int64_t) floor( ((double)max) * sfmt_genrand_real2(_sfmt) );
}


#endif // __AE_JUMPING_MT_H__
