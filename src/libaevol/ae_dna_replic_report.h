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


#ifndef __AE_DNA_REPLIC_REPORT_H__
#define __AE_DNA_REPLIC_REPORT_H__


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_object.h>
#include <ae_list.h>
#include <ae_mutation.h>




// =================================================================
//                          Class declarations
// =================================================================
class ae_dna;





class ae_dna_replic_report : public ae_object
{
  friend class ae_dna;
  
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    ae_dna_replic_report( void );
    ae_dna_replic_report( const ae_dna_replic_report &model );

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_dna_replic_report( void );

    // =================================================================
    //                              Accessors
    // =================================================================
    inline ae_list<ae_mutation*>* get_mutations( void )      const;
    inline ae_list<ae_mutation*>* get_rearrangements( void ) const;

    inline int32_t get_nb_small_mutations( void )  const;
    inline int32_t get_nb_rearrangements( void )   const;
    inline int32_t get_nb_switch( void )           const;
    inline int32_t get_nb_small_insertions( void ) const;
    inline int32_t get_nb_small_deletions( void )  const;
    inline int32_t get_nb_indels( void )           const;
    inline int32_t get_nb_duplications( void )     const;
    inline int32_t get_nb_deletions( void )        const;
    inline int32_t get_nb_translocations( void )   const;
    inline int32_t get_nb_inversions( void )       const;

    // =================================================================
    //                            Public Methods
    // =================================================================
    void compute_stats( void );  // useful when we inspect a tree file
    inline void add_rear( ae_mutation* rear );
    inline void add_mut( ae_mutation* mut );

    // =================================================================
    //                           Public Attributes
    // =================================================================





  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    /*    ae_dna_replic_report( void )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
    ae_dna_replic_report( const ae_dna_replic_report &model )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };*/


    // =================================================================
    //                           Protected Methods
    // =================================================================

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    ae_list<ae_mutation*>* _mutations;       // Lists of mutations and rearrangements undergone
    ae_list<ae_mutation*>* _rearrangements;  // by the genetic unit at last replication

    int32_t _nb_mut[7]; // Number of mutations/rearrangements of each type undergone
};


// =====================================================================
//                          Accessors' definitions
// =====================================================================
inline ae_list<ae_mutation*>* ae_dna_replic_report::get_mutations( void ) const
{
  return _mutations;
}

inline ae_list<ae_mutation*>* ae_dna_replic_report::get_rearrangements( void ) const
{
  return _rearrangements;
}

inline int32_t ae_dna_replic_report::get_nb_small_mutations( void ) const
{
  assert( _mutations->get_nb_elts() == _nb_mut[SWITCH] + _nb_mut[S_INS] + _nb_mut[S_DEL] );
  return _mutations->get_nb_elts();
}

inline int32_t ae_dna_replic_report::get_nb_rearrangements( void ) const
{
  assert( _rearrangements->get_nb_elts() == _nb_mut[DUPL] + _nb_mut[DEL] + _nb_mut[TRANS] + _nb_mut[INV] );
  return _rearrangements->get_nb_elts();
}

inline int32_t ae_dna_replic_report::get_nb_switch( void ) const
{
  return _nb_mut[SWITCH];
}

inline int32_t ae_dna_replic_report::get_nb_small_insertions( void ) const
{
  return _nb_mut[S_INS];
}

inline int32_t ae_dna_replic_report::get_nb_small_deletions( void ) const
{
  return _nb_mut[S_DEL];
}

inline int32_t ae_dna_replic_report::get_nb_indels( void ) const
{
  return _nb_mut[S_INS] + _nb_mut[S_DEL];
}

inline int32_t ae_dna_replic_report::get_nb_duplications( void ) const
{
  return _nb_mut[DUPL];
}

inline int32_t ae_dna_replic_report::get_nb_deletions( void ) const
{
  return _nb_mut[DEL];
}

inline int32_t ae_dna_replic_report::get_nb_translocations( void ) const
{
  return _nb_mut[TRANS];
}

inline int32_t ae_dna_replic_report::get_nb_inversions( void ) const
{
  return _nb_mut[INV];
}

// =====================================================================
//                       Inline functions' definition
// =====================================================================
inline void ae_dna_replic_report::add_mut( ae_mutation* mut )
{
  assert( mut != NULL );
  
  _mutations->add( mut );
  _nb_mut[ mut->get_mut_type() ]++;
}

inline void ae_dna_replic_report::add_rear( ae_mutation* rear )
{
  assert( rear != NULL );
  
  _rearrangements->add( rear );
  _nb_mut[ rear->get_mut_type() ]++;
}


#endif // __AE_DNA_REPLIC_REPORT_H__
