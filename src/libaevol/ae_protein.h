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


#ifndef __AE_PROTEIN_H__
#define __AE_PROTEIN_H__


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>

#include <list>

// =================================================================
//                            Project Files
// =================================================================
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "macros.h"
#include "ae_dna.h"
#include "ae_codon.h"

namespace aevol {



// =================================================================
//                          Class declarations
// =================================================================
class Individual;
class ae_rna;





class ae_protein
{
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    ae_protein() = delete;
    ae_protein(const ae_protein &model) = delete;
    ae_protein( GeneticUnit* gen_unit, const ae_protein &model );
    ae_protein(GeneticUnit* gen_unit,
               const std::list<ae_codon*> codon_list,
               ae_strand strand,
               int32_t shine_dal,
               ae_rna* rna );
    //ae_protein( ae_protein* parent );
    ae_protein( gzFile backup_file );
    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_protein( void );

    // =================================================================
    //                              Accessors
    // =================================================================
    inline ae_strand          get_strand( void )                const;
    inline const std::list<ae_rna*> get_rna_list()          const;
    inline int32_t            get_shine_dal_pos( void )         const;
    inline int32_t            get_first_translated_pos( void )  const;
    inline int32_t            get_last_translated_pos( void )   const;
           int32_t            get_last_STOP_base_pos( void )    const;
    inline double             get_mean( void )                  const;
    inline double             get_width( void )                 const; // returns the half-width
    inline double             get_height( void )                const;
    inline int32_t            get_length( void )                const; // Number of Amino-Acids (not including START and STOP)
    inline double             get_concentration( void )         const;
    inline  bool              get_is_functional( void )         const;
    
    Individual * get_indiv( void ) const;

    // =================================================================
    //                            Public Methods
    // =================================================================
            void  add_RNA( ae_rna* rna );
    char* get_AA_sequence(char separator = ' ') const; // WARNING : creates a new char[...] (up to you to delete it!)
    virtual void  save( gzFile backup_file );

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
    GeneticUnit*  _gen_unit;
    ae_strand         _strand;
    std::list<ae_rna*> rna_list;              // RNAs transcribing this protein
    int32_t           _shine_dal_pos;         // Index of the corresponding shine dalgarno sequence in the genome
    int32_t           _first_translated_pos;  // Index of the first base following the START codon
    int32_t           _last_translated_pos;   // Index of the last base before the STOP codon
    int32_t           _length;                // Number of Amino-Acids (START and STOP codon do NOT produce AAs)
    double            _concentration;
    bool              _is_functional;
    
    std::list<ae_codon*> _AA_list;

    // Phenotypic contribution (triangle) parameters
    double _mean;
    double _width;   // in fact, half-width
    double _height;
};


// =====================================================================
//                          Accessors definitions
// =====================================================================
inline ae_strand ae_protein::get_strand( void ) const
{
  return _strand;
}

inline const std::list<ae_rna*> ae_protein::get_rna_list() const {
  return rna_list;
}

int32_t ae_protein::get_shine_dal_pos( void ) const
{
  return _shine_dal_pos;
}

int32_t ae_protein::get_first_translated_pos( void ) const
{
  return _first_translated_pos;
}

int32_t ae_protein::get_last_translated_pos( void ) const
{
  return _last_translated_pos;
}

double ae_protein::get_mean( void ) const
{
  return _mean;
}

double ae_protein::get_width( void ) const
{
  return _width;
}

double ae_protein::get_height( void ) const
{
  return _height;
}

int32_t ae_protein::get_length( void ) const
{
  return _length;
}

double ae_protein::get_concentration( void ) const
{
  return _concentration;
}

bool ae_protein::get_is_functional( void ) const
{
  return _is_functional;
}


// =====================================================================
//                       Inline functions' definition
// =====================================================================
} // namespace aevol
#endif // __AE_PROTEIN_H__
