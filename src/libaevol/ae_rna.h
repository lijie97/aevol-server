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


#ifndef __AE_RNA_H__
#define __AE_RNA_H__


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
#include "ae_protein.h"

#include "ae_list.h" // TODO vld: to be removed with utility function `distinct`

namespace aevol {



// =================================================================
//                          Class declarations
// =================================================================
class ae_individual;
class ae_genetic_unit;


/// Tell if there are no two consecutive nodes pointing to the same object.
/// Please read me vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
/// ***!#'$       TO REMOVE ALONG WITH INCLUDE AE_LIST !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
/// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
inline bool distinct(ae_list<ae_rna*>* l) {
  for (ae_list_node<ae_rna*>* n = l->get_first();
       n != nullptr and n->get_next() != nullptr;
       n = n->get_next())
    if (n->get_obj() == n->get_next()->get_obj())
      return false;
  return true;
}

inline bool distinct(std::list<ae_rna*>& l) {
  for (auto it = l.begin(); it != l.end() and next(it) != l.end(); ++it)
    if (*it == *next(it))
      return false;
  return true;
}

class ae_rna
{
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    ae_rna(void) = delete;
    ae_rna(const ae_genetic_unit&) = delete;
    ae_rna( ae_genetic_unit* gen_unit, const ae_rna &model );
    ae_rna( ae_genetic_unit* gen_unit, ae_strand strand, int32_t index, int8_t diff );
    //ae_rna( ae_rna* parent );

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_rna( void );

    // =================================================================
    //                              Accessors
    // =================================================================

    // <DEBUG>
    void check( ae_genetic_unit* gen_unit ) { assert( gen_unit == _gen_unit ); };
    //~ void* get_indiv( void ) const { return (void*)_indiv; };
    // </DEBUG>

    inline const ae_genetic_unit * get_genetic_unit(void) const;
    inline void set_genetic_unit(const ae_genetic_unit*  gen_unit);
    inline ae_strand  get_strand( void ) const;
    inline void       set_strand( ae_strand strand );
    inline int32_t    get_promoter_pos( void ) const;
    inline void       set_promoter_pos( int32_t pos );
    inline double     get_basal_level( void ) const;
    inline int32_t    get_transcript_length( void ) const; // The promoter is NOT transcribed.
    inline void       set_transcript_length( int32_t length );
    inline bool       is_coding( void ) const;

    inline const std::list<ae_protein*>& get_transcribed_proteins() const;
    inline void clear_transcribed_proteins() { transcribed_proteins.clear(); };

    // =================================================================
    //                            Public Methods
    // =================================================================
    int32_t get_first_transcribed_pos( void ) const;   // The promoter is NOT transcribed.
    int32_t get_last_transcribed_pos( void ) const;    // The terminator is transcribed.
    inline void add_transcribed_protein( ae_protein* prot );
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
    const ae_genetic_unit*  _gen_unit;
    ae_strand         _strand;
    int32_t           _pos; // Index of the promoter on the genome.
                            // The promoter itself is NOT transcribed
                            // The terminator is transcribed.
    int32_t           _transcript_length;
    double            _basal_level;

    std::list<ae_protein*> transcribed_proteins;
};


// =====================================================================
//                          Accessors' definitions
// =====================================================================
inline const ae_genetic_unit* ae_rna::get_genetic_unit( void ) const
{
  return _gen_unit;
}

inline void ae_rna::set_genetic_unit(const ae_genetic_unit*  gen_unit)
{
  _gen_unit = gen_unit;
}

inline ae_strand ae_rna::get_strand( void ) const
{
  return _strand;
}

inline void ae_rna::set_strand( ae_strand strand )
{
  _strand = strand;
}

inline void ae_rna::set_promoter_pos( int32_t pos )
{
  _pos = pos;
}

inline int32_t ae_rna::get_promoter_pos( void ) const
{
  return _pos;
}

inline double ae_rna::get_basal_level( void ) const
{
  return _basal_level;
}

inline int32_t ae_rna::get_transcript_length( void ) const
{
  return _transcript_length;
}

inline void ae_rna::set_transcript_length( int32_t transcript_length )
{
  _transcript_length = transcript_length;
}

inline const std::list<ae_protein*>& ae_rna::get_transcribed_proteins() const {
  return transcribed_proteins;
}

inline bool ae_rna::is_coding( void ) const
{
  return (not transcribed_proteins.empty());
}

// =====================================================================
//                       Inline functions' definition
// =====================================================================
void ae_rna::add_transcribed_protein( ae_protein* prot )
{
  transcribed_proteins.push_back(prot);
}

void ae_rna::shift_position( int32_t delta_pos, int32_t genome_length )
{
  _pos = ae_utils::mod( _pos + delta_pos, genome_length );
}

} // namespace aevol
#endif // __AE_RNA_H__
