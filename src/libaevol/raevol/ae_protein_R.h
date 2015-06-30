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

#ifndef AE_PROTEIN_R_H
#define AE_PROTEIN_R_H

// =================================================================
//                              Libraries
// =================================================================

// =================================================================
//                            Project Files
// =================================================================
#include "ae_protein.h"
#include "ae_codon.h"

namespace aevol {

//##############################################################################
//                                                                             #
//                           Class ae_protein_R                                #
//                                                                             #
//##############################################################################
class ae_genetic_unit;

class ae_protein_R : public Protein
{
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    ae_protein_R( GeneticUnit* gen_unit, const ae_protein_R &model );
    ae_protein_R( GeneticUnit* gen_unit,
    		const std::list<ae_codon*> codon_list,
    		ae_strand strand,
    		int32_t shine_dal_pos,
    		ae_rna* rna ); // TODO ae_rna_R?
    ae_protein_R( const std::list<ae_codon*> codon_list, double concentration);
	ae_protein_R( gzFile backup_file );

    // =================================================================
    //                             Destructors
    // =================================================================
    ~ae_protein_R( void );

    // =================================================================
    //                              Accessors
    // =================================================================
//    inline std::vector<ae_influence_R*> get_influence_list( void );
    inline void     set_inherited( bool is_inherited );
    inline void     set_signal( bool is_signal);
    inline bool     is_inherited( void );
    inline bool     is_signal( void );

    // =================================================================
    //                            Public Methods
    // =================================================================
    //inline ae_protein_R* copy( void );
    inline void    multiply_concentration( double factor );
    inline void    set_concentration ( double concentration);
    inline void    update_concentration( void );
    inline void    reset_concentration( void );
    inline void    set_initial_concentration( void );
           void    compute_delta_concentration( void );
           int8_t  get_codon( int32_t index );
//           void    add_influence( ae_influence_R* influence );
	         void    save( gzFile backup_file );
//	         void    remove_influence( ae_influence_R* influence );


    // =================================================================
    //                           Public Attributes
    // =================================================================
	bool not_pure_TF;

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
    std::vector<ae_rna_R*>  _rna_R_list;
    double    _delta_concentration;
    bool      _inherited;
    bool      _signal;
    double    _initial_concentration; // concentration at cell birth
};

// =====================================================================
//                          Accessors definitions
// =====================================================================
//std::vector<ae_influence_R*> ae_protein_R::get_influence_list( void )
//{
//  return _influence_list;
//}

// =====================================================================
//                       Inline functions' definition
// =====================================================================
inline void ae_protein_R::update_concentration( void )
{
  _concentration += _delta_concentration;
}

inline void ae_protein_R::set_inherited( bool is_inherited )
{
  _inherited = is_inherited;
}

inline void ae_protein_R::set_signal( bool is_signal )
{
  _signal = is_signal;
}

inline void ae_protein_R::reset_concentration( void )
{
  _concentration = _initial_concentration;
}

inline void ae_protein_R::set_initial_concentration( void )
{
  _initial_concentration = _concentration;
}

inline bool ae_protein_R::is_inherited( void )
{
  return _inherited;
}

inline bool ae_protein_R::is_signal( void )
{
  return _signal;
}

/*
ae_protein_R* ae_protein_R::copy( void )
{
  ae_protein_R* new_prot = new ae_protein_R( this );
  new_prot->_shine_dal_pos = -1;

  return new_prot;
}
*/

inline void ae_protein_R::multiply_concentration( double factor )
{
  _concentration *= factor;
}

inline void ae_protein_R::set_concentration( double concentration )
{
  _concentration = concentration;
}

} // namespace aevol

#endif // AE_PROTEIN_R_H
