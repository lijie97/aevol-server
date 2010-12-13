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
 
 
 #ifndef __AE_DNA_H__
#define  __AE_DNA_H__
 
 
// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>
#include <stdlib.h>
#include <zlib.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_string.h>
#include <ae_list.h>
#include <ae_mutation.h>
#include <ae_dna_replic_report.h>




// =================================================================
//                          Class declarations
// =================================================================
class ae_genetic_unit;





enum ae_strand
{
  LEADING = 0,
  LAGGING = 1
};




 
class ae_dna : public ae_string
{  
  public :
  
    // =================================================================
    //                             Constructors
    // =================================================================
    ae_dna( ae_genetic_unit* gen_unit, int32_t length );
    ae_dna( ae_genetic_unit* gen_unit, const ae_dna &model );
    ae_dna( ae_genetic_unit* gen_unit, ae_dna* const parent_dna );
    ae_dna( ae_genetic_unit* gen_unit, char* seq, int32_t length );
    ae_dna( ae_genetic_unit* gen_unit, gzFile* backup_file );
  
    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_dna( void );
  
    // =================================================================
    //                              Accessors
    // =================================================================
    // From ae_string
    //   inline const char*   get_data( void ) const;
    //   inline       void    set_data( char* data, int32_t length = -1 );
    //   inline       int32_t get_length( void ) const;
    inline ae_dna_replic_report*  get_replic_report( void ) const;
    inline void                   set_replic_report( ae_dna_replic_report * rep ); // for post-treatment only
    
    char* get_subsequence( int32_t from, int32_t to, ae_strand strand ) const; // WARNING : creates a new char[...] (up to you to delete it!)


    // =================================================================
    //                            Public Methods
    // =================================================================
    // Perform all the mutations (local mutations and rearrangements)
    void perform_mutations( void );
    
    // Perform all the local mutations (point mutations and indels) of the replication
    void do_small_mutations( void );
    
    // Perform all the chromosomic rearrangements (duplications, deletions, translocations and inversions)
    // of the replication
    void do_rearrangements( void );
    void do_rearrangements_with_align( void );
    
    // Perform a single local mutation at a random position
    ae_mutation* do_switch( void );
    ae_mutation* do_small_insertion( void );
    ae_mutation* do_small_deletion( void );

    // Perform a single local mutation at a specified position (useful to replay the evolution)
    bool do_switch( int32_t pos );
    bool do_small_insertion( int32_t pos, int16_t nb_insert, char * seq );
    bool do_small_deletion( int32_t pos, int16_t nb_del );
    
    // Perform a single rearrangement at random positions
    ae_mutation* do_duplication( void );
    ae_mutation* do_deletion( void );
    ae_mutation* do_translocation( void );
    ae_mutation* do_inter_GU_translocation( void );
    ae_mutation* do_inversion( void );
    ae_mutation* do_insertion( const char* seq_to_insert, int32_t seq_length = -1 );
  
    // Perform a single rearrangement at specified positions
    bool do_duplication( int32_t pos_1, int32_t pos_2, int32_t pos_3 );
    bool do_deletion( int32_t pos_1, int32_t pos_2 );
    bool do_translocation( int32_t pos_1, int32_t pos_2, int32_t pos_3, int32_t pos_4, bool invert );
    bool do_inter_GU_translocation( int32_t pos_1, int32_t pos_2, int32_t pos_3, int32_t pos_4, bool invert );
    bool do_inversion( int32_t pos_1, int32_t pos_2 );
    bool do_insertion( int32_t pos, const char* seq_to_insert, int32_t seq_length );
    
    void undergo_this_mutation( ae_mutation * mut ); // useful when we replay the evolution

    void compute_statistical_data( void );
    
    
    static void set_GU( ae_list** rna_list, ae_genetic_unit* GU );
  
  
    // =================================================================
    //                           Public Attributes
    // =================================================================
  
  
  
  
  
  protected :
  
    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    ae_dna( void )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
    ae_dna( const ae_dna &model )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
  
    // =================================================================
    //                           Protected Methods
    // =================================================================
    ae_genetic_unit*  extract_into_new_GU( int32_t pos_1, int32_t pos_2 );
    void insert_GU( ae_genetic_unit* GU_to_insert, int32_t pos_B, int32_t pos_D, bool invert );
    
    void ABCDE_to_ADCBE(   int32_t pos_B, int32_t pos_C, int32_t pos_D, int32_t pos_E );
    void ABCDE_to_ADBpCpE( int32_t pos_B, int32_t pos_C, int32_t pos_D, int32_t pos_E );
    void ABCDE_to_ACpDpBE( int32_t pos_B, int32_t pos_C, int32_t pos_D, int32_t pos_E );
    void inter_GU_ABCDE_to_ACDBE( int32_t pos_B, int32_t pos_C, int32_t pos_E );
    void inter_GU_ABCDE_to_BDCAE( int32_t pos_B, int32_t pos_C, int32_t pos_E );
    
    
    // =================================================================
    //                          Protected Attributes
    // =================================================================
    ae_genetic_unit*      _gen_unit; // Genetic unit which the genetic unit belongs to
    ae_dna_replic_report* _replic_report;
};


// =====================================================================
//                          Accessors definitions
// =====================================================================
inline ae_dna_replic_report* ae_dna::get_replic_report( void ) const
{
  return _replic_report;
}

 // for post-treatment only
inline void ae_dna::set_replic_report( ae_dna_replic_report * rep )
{
  _replic_report = rep;
}


// =====================================================================
//                       Inline functions' definition
// =====================================================================


#endif // __AE_DNA_H__
