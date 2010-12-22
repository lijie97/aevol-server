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
    ae_replication_report( ae_individual * indiv );

    // Creates a completely independent copy of the original report
    ae_replication_report( const ae_replication_report &model );

    ae_replication_report( gzFile * backup_file );

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_replication_report( void );

    // =================================================================
    //                              Accessors
    // =================================================================
    inline ae_individual * get_indiv( void ) const;
    int32_t         get_index( void ) const;
    void            set_index( int32_t index );
    int32_t         get_rank( void ) const;
    void            set_rank( int32_t rank );
    inline void     set_parent_index( int32_t parent_index );
    inline int32_t  get_parent_index( void ) const;
    inline void     set_donnor_index( int32_t donnor_index );
    inline int32_t  get_donnor_index( void ) const;
    inline void     set_parent_metabolic_error( double parent_metabolic_error );
    inline double   get_parent_metabolic_error( void ) const;
    inline void     set_parent_secretion_error( double parent_secretion_error );
    inline double   get_parent_secretion_error( void ) const;
    inline void     set_donnor_metabolic_error( double donnor_metabolic_error );
    inline double   get_donnor_metabolic_error( void ) const;
    inline void     set_parent_genome_size( int32_t parent_genome_size );
    inline int32_t  get_parent_genome_size( void ) const;
    inline void     set_donnor_genome_size( int32_t donnor_genome_size );
    inline int32_t  get_donnor_genome_size( void ) const;
    
    inline ae_list*  get_dna_replic_reports( void ) const;


    // =================================================================
    //                            Public Methods
    // =================================================================
    //~ inline void add_mutation_list( ae_list* mut_list );
    //~ inline void add_rearrangement_list( ae_list* rear_list );
    
    void write_to_backup( gzFile* backup_file );


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
    int32_t         _index;
    int32_t         _rank;
    int32_t         _parent_index;
    int32_t         _donnor_index; // Horizontal transfer donnor
    double          _parent_metabolic_error;
    double          _parent_secretion_error;
    double          _donnor_metabolic_error;
    int32_t         _parent_genome_size;
    int32_t         _donnor_genome_size;
    
    ae_list*  _dna_replic_reports; // List of each genetic unit's replicatino report
    
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

inline int32_t ae_replication_report::get_index( void ) const
{
  return _index;
}

inline void ae_replication_report::set_index( int32_t index )
{
  _index = index;
}

inline int32_t ae_replication_report::get_rank( void ) const
{
  return _rank;
}

inline void ae_replication_report::set_rank( int32_t rank )
{
  _rank = rank;
}


int32_t ae_replication_report::get_parent_index( void ) const
{
  return _parent_index;
}

void ae_replication_report::set_parent_index( int32_t parent_index )
{
  _parent_index = parent_index;
}

int32_t ae_replication_report::get_donnor_index( void ) const
{
  return _donnor_index;
}

void ae_replication_report::set_donnor_index( int32_t donnor_index )
{
  _donnor_index = donnor_index;
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

void ae_replication_report::set_donnor_metabolic_error( double donnor_metabolic_error )
{
  _donnor_metabolic_error = donnor_metabolic_error;
}

double ae_replication_report::get_donnor_metabolic_error( void ) const
{
  return _donnor_metabolic_error;
}


void ae_replication_report::set_parent_genome_size( int32_t parent_genome_size )
{
  _parent_genome_size = parent_genome_size;
}

int32_t ae_replication_report::get_parent_genome_size( void ) const
{
  return _parent_genome_size;
}

void ae_replication_report::set_donnor_genome_size( int32_t donnor_genome_size )
{
  _donnor_genome_size = donnor_genome_size;
}

int32_t ae_replication_report::get_donnor_genome_size( void ) const
{
  return _donnor_genome_size;
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
