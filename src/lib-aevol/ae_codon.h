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


 #ifndef __AE_CODON_H__
#define  __AE_CODON_H__


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>





// =================================================================
//                            Project Files
// =================================================================
#include <ae_object.h>
#include <ae_dna.h>




// =================================================================
//                          Class declarations
// =================================================================






class ae_codon : public ae_object
{
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    ae_codon( void );
    ae_codon( const ae_codon &model );
    ae_codon( int8_t value );
    ae_codon( ae_dna* genome, ae_strand strand, int32_t index );
    ae_codon( gzFile backup_file );

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_codon( void );

    // =================================================================
    //                              Accessors
    // =================================================================
    inline int8_t get_value( void );

    // =================================================================
    //                            Public Methods
    // =================================================================
    inline bool is_start( void );
    inline bool is_stop( void );
    inline ae_codon* copy( void );
    void   write_to_backup( gzFile backup_file );

    // =================================================================
    //                           Public Attributes
    // =================================================================





  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    //~ ae_codon( void )
    //~ {
      //~ printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      //~ exit( EXIT_FAILURE );
    //~ };
    //~ ae_codon( const ae_codon &model )
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
    int8_t _value;
};


// =====================================================================
//                          Accessors' definitions
// =====================================================================
int8_t ae_codon::get_value( void )
{
  return _value;
}

bool ae_codon::is_start( void )
{
  return _value == CODON_START;
}

bool ae_codon::is_stop( void )
{
  return _value == CODON_STOP;
}

// =====================================================================
//                       Inline functions' definition
// =====================================================================
ae_codon* ae_codon::copy( void )
{
  return new ae_codon( _value );
}

#endif // __AE_CODON_H__
