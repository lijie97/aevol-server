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


 #ifndef __AE_PROTEIN_R_H__
#define  __AE_PROTEIN_R_H__


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>

// =================================================================
//                            Project Files
// =================================================================
#include <ae_protein.h>
#include <ae_rna_R.h>

// =================================================================
//                          Class declarations
// =================================================================
class ae_genetic_unit;


class ae_protein_R : public ae_protein
{
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    ae_protein_R( ae_genetic_unit* gen_unit, const ae_protein_R &model );
    ae_protein_R( ae_genetic_unit* gen_unit, ae_list* codon_list, ae_strand strand, int32_t shine_dal_pos, ae_rna* rna ); // TODO ae_rna_R?
	ae_protein_R( gzFile* backup_file );
	
    // =================================================================
    //                             Destructors
    // =================================================================
    ~ae_protein_R( void );

    // =================================================================
    //                              Accessors
    // =================================================================
    inline ae_list* get_influence_list( void );
    inline void     set_inherited( bool is_inherited );
    inline bool     is_inherited( void );

    // =================================================================
    //                            Public Methods
    // =================================================================
    //inline ae_protein_R* copy( void );
    inline void multiply_concentration( double factor );
    inline void update_concentration( void );
    void        compute_delta_concentration( void );
    int8_t      get_codon( int32_t index );
    void        add_influence( ae_influence_R* influence );
    void        save( gzFile* backup_file );

    // =================================================================
    //                           Public Attributes
    // =================================================================

  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    ae_protein_R( const ae_protein_R &model )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };

    // =================================================================
    //                           Protected Methods
    // =================================================================
    void remove_influences( void );

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    ae_list*  _influence_list;
    double    _delta_concentration;
    bool      _inherited;
};

// =====================================================================
//                          Accessors definitions
// =====================================================================
ae_list* ae_protein_R::get_influence_list( void )
{
  return _influence_list;
}

// =====================================================================
//                       Inline functions' definition
// =====================================================================
inline void ae_protein_R::update_concentration( void )
{
  _concentration += _delta_concentration;
}

inline void ae_protein_R::set_inherited( bool inherited )
{
  _inherited = inherited;
}

inline bool ae_protein_R::is_inherited( void )
{
  return _inherited;
}

/*
ae_protein_R* ae_protein_R::copy( void )
{
  ae_protein_R* new_prot = new ae_protein_R( this );
  new_prot->_shine_dal_pos = -1;

  return new_prot;
}
*/

void ae_protein_R::multiply_concentration( double factor )
{
  _concentration *= factor;
}

#endif // __AE_PROTEIN_R_H__
