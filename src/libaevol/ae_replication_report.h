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


 #ifndef __AE_REPLICATION_REPORT_H__
#define  __AE_REPLICATION_REPORT_H__


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_object.h>
#include <ae_list.h>
#include <zlib.h>




// =================================================================
//                          Class declarations
// =================================================================
class ae_individual;





class ae_replication_report : public ae_object
{
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    ae_replication_report( ae_individual * indiv, ae_individual * parent, ae_individual * donor = NULL );

    // Creates a completely independent copy of the original report
    ae_replication_report( const ae_replication_report &model );

    ae_replication_report( gzFile * tree_file, ae_individual * indiv );

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_replication_report( void );

    // =================================================================
    //                              Accessors
    // =================================================================
    inline ae_individual * get_indiv( void ) const;
    int32_t         get_id( void ) const;
    void            set_id( int32_t id );
    int32_t         get_rank( void ) const;
    void            set_rank( int32_t rank );
    inline void     set_parent_id( int32_t parent_id );
    inline int32_t  get_parent_id( void ) const;
    inline void     set_donor_id( int32_t donor_id );
    inline int32_t  get_donor_id( void ) const;
    inline void     set_parent_metabolic_error( double parent_metabolic_error );
    inline double   get_parent_metabolic_error( void ) const;
    inline void     set_parent_secretion_error( double parent_secretion_error );
    inline double   get_parent_secretion_error( void ) const;
    inline void     set_donor_metabolic_error( double donor_metabolic_error );
    inline double   get_donor_metabolic_error( void ) const;
    inline void     set_parent_genome_size( int32_t parent_genome_size );
    inline int32_t  get_parent_genome_size( void ) const;
    inline void     set_donor_genome_size( int32_t donor_genome_size );
    inline int32_t  get_donor_genome_size( void ) const;
    inline double   get_mean_align_score( void ) const;
    
    inline ae_list* get_dna_replic_reports( void ) const;


    // =================================================================
    //                            Public Methods
    // =================================================================
    //~ inline void add_mutation_list( ae_list* mut_list );
    //~ inline void add_rearrangement_list( ae_list* rear_list );
    
    void signal_end_of_replication( void );
    void write_to_tree_file( gzFile* tree_file ) const;


    // =================================================================
    //                           Public Attributes
    // =================================================================





  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    ae_replication_report( void )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };



    // =================================================================
    //                           Protected Methods
    // =================================================================

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    ae_individual*  _indiv;
    
    // ********** Data recorded in the tree **********
    int32_t         _id;
    int32_t         _rank;
    int32_t         _parent_id;
    int32_t         _donor_id; // Horizontal transfer donor
    
    int32_t         _genome_size;
    double          _metabolic_error;
    int16_t         _nb_genes_activ;
    int16_t         _nb_genes_inhib;
    int16_t         _nb_non_fun_genes;
    int16_t         _nb_coding_RNAs;
    int16_t         _nb_non_coding_RNAs;
    
    ae_list*  _dna_replic_reports; // List of each genetic unit's replication report
    // ********** END Data recorded in the tree **********
    
    double          _parent_metabolic_error;
    double          _parent_secretion_error;
    double          _donor_metabolic_error;
    int32_t         _parent_genome_size;
    int32_t         _donor_genome_size;
    
    double          _mean_align_score;
    
    //~ ae_list*  _mutations;       // These are lists of lists of mutations (resp rearrangements)
    //~ ae_list*  _rearrangements;  // undergone during the replication.
                                //~ //  The first node of the list will correspond to the list of mutations
                                //~ // (resp rears) of the first genetic unit of the individual and so forth.
                                //~ // Example : to access the first rearrangement that occured on the first plasmid
                                //~ // (i.e the second genetic unit), one will have to type : 
                                //~ //    ae_list* rear_list = (ae_list*) _rearrangements->get_first()->get_next()->get_obj();
                                //~ //    ae_mutation* my_mut = (ae_mutation*) rear_list->get_first()->get_obj();
};


// =====================================================================
//                          Accessors' definitions
// =====================================================================
inline ae_individual * ae_replication_report::get_indiv( void ) const
{
  return _indiv;
}

inline int32_t ae_replication_report::get_id( void ) const
{
  return _id;
}

inline void ae_replication_report::set_id( int32_t id )
{
  _id = id;
}

inline int32_t ae_replication_report::get_rank( void ) const
{
  return _rank;
}

inline void ae_replication_report::set_rank( int32_t rank )
{
  _rank = rank;
}


int32_t ae_replication_report::get_parent_id( void ) const
{
  return _parent_id;
}

void ae_replication_report::set_parent_id( int32_t parent_id )
{
  _parent_id = parent_id;
}

int32_t ae_replication_report::get_donor_id( void ) const
{
  return _donor_id;
}

void ae_replication_report::set_donor_id( int32_t donor_id )
{
  _donor_id = donor_id;
}

void ae_replication_report::set_parent_metabolic_error( double parent_metabolic_error )
{
  _parent_metabolic_error = parent_metabolic_error;
}

double ae_replication_report::get_parent_metabolic_error( void ) const
{
  return _parent_metabolic_error;
}

void ae_replication_report::set_parent_secretion_error( double parent_secretion_error )
{
  _parent_secretion_error = parent_secretion_error;
}

double ae_replication_report::get_parent_secretion_error( void ) const
{
  return _parent_secretion_error;
}

void ae_replication_report::set_donor_metabolic_error( double donor_metabolic_error )
{
  _donor_metabolic_error = donor_metabolic_error;
}

double ae_replication_report::get_donor_metabolic_error( void ) const
{
  return _donor_metabolic_error;
}


void ae_replication_report::set_parent_genome_size( int32_t parent_genome_size )
{
  _parent_genome_size = parent_genome_size;
}

int32_t ae_replication_report::get_parent_genome_size( void ) const
{
  return _parent_genome_size;
}

void ae_replication_report::set_donor_genome_size( int32_t donor_genome_size )
{
  _donor_genome_size = donor_genome_size;
}

int32_t ae_replication_report::get_donor_genome_size( void ) const
{
  return _donor_genome_size;
}

inline double ae_replication_report::get_mean_align_score( void ) const
{
  return _mean_align_score;
}

inline ae_list* ae_replication_report::get_dna_replic_reports( void ) const
{
  return _dna_replic_reports;
}



// =====================================================================
//                       Inline functions' definition
// =====================================================================
//~ void ae_replication_report::add_mutation_list( ae_list* mut_list )
//~ {
  //~ _mutations->add( mut_list );
//~ }

//~ void ae_replication_report::add_rearrangement_list( ae_list* rear_list )
//~ {
  //~ _rearrangements->add( rear_list );
//~ }


#endif // __AE_REPLICATION_REPORT_H__
