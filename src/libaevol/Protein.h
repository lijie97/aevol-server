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


#ifndef AEVOL_PROTEIN_H__
#define AEVOL_PROTEIN_H__


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
#include "Dna.h"
#include "Codon.h"

namespace aevol {



// =================================================================
//                          Class declarations
// =================================================================
class Individual;
class Rna;





class Protein
{
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    Protein() = delete;
    Protein(const Protein &model) = delete;
    Protein(GeneticUnit* gen_unit, const Protein &model);
    Protein(GeneticUnit* gen_unit,
               const std::list<Codon*>& codon_list,
               Strand strand,
               int32_t shine_dal_pos,
               Rna * rna,
               double w_max);
    //Protein( Protein* parent );
    Protein( gzFile backup_file );
    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~Protein( void );

    // =================================================================
    //                              Accessors
    // =================================================================
    inline Strand get_strand( void )                const;
    inline const std::list<Rna *> get_rna_list()          const;
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
    void  add_RNA( Rna * rna );
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
    Strand _strand;
    std::list<Rna *> rna_list;              // RNAs transcribing this protein
    int32_t           _shine_dal_pos;         // Index of the corresponding shine dalgarno sequence in the genome
    int32_t           _first_translated_pos;  // Index of the first base following the START codon
    int32_t           _last_translated_pos;   // Index of the last base before the STOP codon
    int32_t           _length;                // Number of Amino-Acids (START and STOP codon do NOT produce AAs)
    double            _concentration;
    bool              _is_functional;
    
    std::list<Codon *> _AA_list;

    // Phenotypic contribution (triangle) parameters
    double _mean;
    double _width;   // in fact, half-width
    double _height;
};


// =====================================================================
//                          Accessors definitions
// =====================================================================
inline Strand Protein::get_strand( void ) const
{
  return _strand;
}

inline const std::list<Rna *> Protein::get_rna_list() const {
  return rna_list;
}

int32_t Protein::get_shine_dal_pos( void ) const
{
  return _shine_dal_pos;
}

int32_t Protein::get_first_translated_pos( void ) const
{
  return _first_translated_pos;
}

int32_t Protein::get_last_translated_pos( void ) const
{
  return _last_translated_pos;
}

double Protein::get_mean( void ) const
{
  return _mean;
}

double Protein::get_width( void ) const
{
  return _width;
}

double Protein::get_height( void ) const
{
  return _height;
}

int32_t Protein::get_length( void ) const
{
  return _length;
}

double Protein::get_concentration( void ) const
{
  return _concentration;
}

bool Protein::get_is_functional( void ) const
{
  return _is_functional;
}


// =====================================================================
//                       Inline functions' definition
// =====================================================================
} // namespace aevol
#endif // AEVOL_PROTEIN_H__
