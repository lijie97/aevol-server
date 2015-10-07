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


#ifndef AEVOL_RNA_H__
#define AEVOL_RNA_H__


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
#include "Protein.h"

namespace aevol {



// =================================================================
//                          Class declarations
// =================================================================
class Individual;
class GeneticUnit;

class Rna
{
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    Rna(void) = delete;
    Rna(const GeneticUnit&) = delete;
    Rna( GeneticUnit* gen_unit, const Rna &model );
    Rna( GeneticUnit* gen_unit, Strand strand, int32_t index, int8_t diff );
    //Rna( Rna* parent );

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~Rna( void );

    // =================================================================
    //                              Accessors
    // =================================================================

    // <DEBUG>
    void check( GeneticUnit* gen_unit ) { assert( gen_unit == _gen_unit ); };
    //~ void* get_indiv( void ) const { return (void*)indiv_; };
    // </DEBUG>

    inline const GeneticUnit * get_genetic_unit(void) const;
    inline void set_genetic_unit(const GeneticUnit*  gen_unit);
    inline Strand get_strand( void ) const;
    inline void       set_strand( Strand strand );
    inline int32_t    get_promoter_pos( void ) const;
    inline void       set_promoter_pos( int32_t pos );
    inline double     get_basal_level( void ) const;
    inline int32_t    get_transcript_length( void ) const; // The promoter is NOT transcribed.
    inline void       set_transcript_length( int32_t length );
    inline bool       is_coding( void ) const;

    inline const std::list<Protein *>& get_transcribed_proteins() const;
    inline void clear_transcribed_proteins() { transcribed_proteins.clear(); };

    // =================================================================
    //                            Public Methods
    // =================================================================
    int32_t get_first_transcribed_pos( void ) const;   // The promoter is NOT transcribed.
    int32_t get_last_transcribed_pos( void ) const;    // The terminator is transcribed.
    inline void add_transcribed_protein( Protein * prot );
    inline void shift_position( int32_t delta_pos, int32_t genome_length );

    // =================================================================
    //                           Public Attributes
    // =================================================================





  protected :
    // =================================================================
    //                         Forbidden Constructors
    // =================================================================

    // =================================================================
    //                           Protected Methods
    // =================================================================

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    const GeneticUnit*  _gen_unit;
    Strand _strand;
    int32_t _pos; // Index of the promoter on the genome.
                  // The promoter itself is NOT transcribed
                  // The terminator is transcribed.
    int32_t _transcript_length;
    double _basal_level;

    // Access list to the proteins transcribed by this rna
    std::list<Protein*> transcribed_proteins;
};


// =====================================================================
//                          Accessors' definitions
// =====================================================================
inline const GeneticUnit*Rna::get_genetic_unit( void ) const
{
  return _gen_unit;
}

inline void Rna::set_genetic_unit(const GeneticUnit*  gen_unit)
{
  _gen_unit = gen_unit;
}

inline Strand Rna::get_strand( void ) const
{
  return _strand;
}

inline void Rna::set_strand( Strand strand )
{
  _strand = strand;
}

void Rna::set_promoter_pos( int32_t pos )
{
  _pos = pos;
}

inline int32_t Rna::get_promoter_pos( void ) const
{
  return _pos;
}

inline double Rna::get_basal_level( void ) const
{
  return _basal_level;
}

inline int32_t Rna::get_transcript_length( void ) const
{
  return _transcript_length;
}

inline void Rna::set_transcript_length( int32_t transcript_length )
{
  _transcript_length = transcript_length;
}

inline const std::list<Protein *>&Rna::get_transcribed_proteins() const {
  return transcribed_proteins;
}

inline bool Rna::is_coding( void ) const
{
  return (not transcribed_proteins.empty());
}

// =====================================================================
//                       Inline functions' definition
// =====================================================================
void Rna::add_transcribed_protein( Protein * prot )
{
  transcribed_proteins.push_back(prot);
}

void Rna::shift_position( int32_t delta_pos, int32_t genome_length )
{
  _pos = Utils::mod( _pos + delta_pos, genome_length );
}

} // namespace aevol
#endif // AEVOL_RNA_H__
