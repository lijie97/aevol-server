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


#ifndef __AE_DNA_REPLIC_REPORT_H__
#define __AE_DNA_REPLIC_REPORT_H__


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
#include <ae_list.h>
#include <ae_mutation.h>
namespace aevol {




// =================================================================
//                          Class declarations
// =================================================================
class ae_dna;





class ae_dna_replic_report
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
    inline ae_list<ae_mutation*>* get_HT( void ) const;

    inline int32_t get_nb_small_mutations( void )  const;
    inline int32_t get_nb_rearrangements( void )   const;
    inline int32_t get_nb_HT( void )   const;
    inline int32_t get_nb_switch( void )           const;
    inline int32_t get_nb_small_insertions( void ) const;
    inline int32_t get_nb_small_deletions( void )  const;
    inline int32_t get_nb_indels( void )           const;
    inline int32_t get_nb_duplications( void )     const;
    inline int32_t get_nb_deletions( void )        const;
    inline int32_t get_nb_translocations( void )   const;
    inline int32_t get_nb_inversions( void )       const;
    inline int32_t get_nb_ins_HT( void )           const;
    inline int32_t get_nb_repl_HT( void )          const;

    // =================================================================
    //                            Public Methods
    // =================================================================
    void compute_stats( void );  // useful when we inspect a tree file
    inline void add_rear( ae_mutation* rear );
    inline void add_mut( ae_mutation* mut );
    inline void add_HT( ae_mutation* HT);

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
    ae_list<ae_mutation*>* _mutations;       // Lists of mutations, rearrangements and undergone
    ae_list<ae_mutation*>* _rearrangements;  // by the genetic unit at last replication
    ae_list<ae_mutation*>* _HT;              

    int32_t _nb_mut[10]; // Number of mutations/rearrangements/HT of each type undergone
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

inline ae_list<ae_mutation*>* ae_dna_replic_report::get_HT( void ) const
{
  return _HT;
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

inline int32_t ae_dna_replic_report::get_nb_HT( void ) const
{
  assert( _HT->get_nb_elts() == _nb_mut[INS_HT] + _nb_mut[REPL_HT] );
  return _HT->get_nb_elts();
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

inline int32_t ae_dna_replic_report::get_nb_ins_HT( void ) const
{
  return _nb_mut[INS_HT];
}

inline int32_t ae_dna_replic_report::get_nb_repl_HT( void ) const
{
  return _nb_mut[REPL_HT];
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

inline void ae_dna_replic_report::add_HT( ae_mutation* HT )
{
  assert( HT != NULL );
  _HT->add( HT );
  _nb_mut[ HT->get_mut_type() ]++;
}

} // namespace aevol
#endif // __AE_DNA_REPLIC_REPORT_H__
