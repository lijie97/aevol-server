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
 
 
// =================================================================
//                              Libraries
// =================================================================



// =================================================================
//                            Project Files
// =================================================================
#include <ae_replication_report.h>
#include <ae_dna_replic_report.h>
#include <ae_mutation.h>
#include <ae_individual.h>




//##############################################################################
//                                                                             #
//                         Class ae_replication_report                         #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
ae_replication_report::ae_replication_report( ae_individual * indiv )
{
  _indiv = indiv;
  
  _index = indiv->get_index_in_population();
  _rank  = indiv->get_rank_in_population();
  
  _parent_index = -1;
  _donnor_index = -1;
  
  _parent_metabolic_error = 0.0;
  _parent_secretion_error = 0.0;
  _donnor_metabolic_error = 0.0;
  _parent_genome_size     = 0;
  _donnor_genome_size     = 0;
  
  _dna_replic_reports = new ae_list();
}


// Creates an independent copy of the original report
ae_replication_report::ae_replication_report( const ae_replication_report &model )
{
  _parent_index = model._parent_index;
  _donnor_index = model._donnor_index;
  
  _index = model._index;
  _rank  = model._rank;

  _parent_metabolic_error = model._parent_metabolic_error;
  _parent_secretion_error = model._parent_secretion_error;
  _donnor_metabolic_error = model._donnor_metabolic_error;
  _parent_genome_size     = model._parent_genome_size;
  _donnor_genome_size     = model._donnor_genome_size;
  
  _dna_replic_reports = new ae_list();

  ae_list_node * node = (model._dna_replic_reports)->get_first();
  ae_dna_replic_report * dnarep = NULL;
  while ( node != NULL )
  {
    dnarep = (ae_dna_replic_report *) node->get_obj();
    _dna_replic_reports->add( new ae_dna_replic_report(*dnarep) );
    node = node->get_next();
  }
}


ae_replication_report::ae_replication_report( gzFile* backup_file )
{
  gzread( backup_file, &_index,         sizeof(_index)        );
  gzread( backup_file, &_rank,          sizeof(_rank)         );
  gzread( backup_file, &_parent_index,  sizeof(_parent_index) );
  gzread( backup_file, &_donnor_index,  sizeof(_donnor_index) );
  //~ printf( "  _index : %"PRId32"\n",         _index );
  //~ printf( "  _rank : %"PRId32"\n",          _rank );
  //~ printf( "  _parent_index : %"PRId32"\n", _parent_index );
  //~ printf( "  _donnor_index : %"PRId32"\n", _donnor_index );

  gzread( backup_file, &_parent_metabolic_error,  sizeof(_parent_metabolic_error) );
  gzread( backup_file, &_parent_secretion_error,  sizeof(_parent_secretion_error) );
  gzread( backup_file, &_donnor_metabolic_error,  sizeof(_donnor_metabolic_error) );
  gzread( backup_file, &_parent_genome_size,      sizeof(_parent_genome_size) );
  gzread( backup_file, &_donnor_genome_size,      sizeof(_donnor_genome_size) );
  //~ printf( "  _parent_metabolic_error : %lf\n", _parent_metabolic_error );
  //~ printf( "  _parent_secretion_error : %lf\n", _parent_secretion_error );
  //~ printf( "  _donnor_metabolic_error : %lf\n", _donnor_metabolic_error );
  //~ printf( "  _parent_genome_size : %"PRId32"\n", _parent_genome_size );
  //~ printf( "  _donnor_genome_size : %"PRId32"\n", _donnor_genome_size );
  
  int32_t nb_dna_replic_reports;
  gzread( backup_file, &nb_dna_replic_reports, sizeof(nb_dna_replic_reports) );
  //~ printf( "  nb_dna_replic_reports : %"PRId32"\n", nb_dna_replic_reports );
  
  
  _dna_replic_reports = new ae_list();

  int32_t mydnareport, myevent;
  int32_t nb_rears, nb_muts;
  ae_dna_replic_report * dnareport = NULL;
  ae_mutation * event = NULL;

  for ( mydnareport = 0 ; mydnareport < nb_dna_replic_reports ; mydnareport++ )
  {
    dnareport = new ae_dna_replic_report();
    
    gzread( backup_file, &nb_rears, sizeof(nb_rears) );
    //~ printf( "  nb_rears : %"PRId32"\n", nb_rears );
    for ( myevent  = 0 ; myevent < nb_rears ; myevent++ )
    {
      event = new ae_mutation( backup_file );
      dnareport->get_rearrangements()->add( event );
    }
    
    gzread( backup_file, &nb_muts, sizeof(nb_muts) );
    //~ printf( "  nb_rears : %"PRId32"\n", nb_rears );
    for(myevent  = 0 ; myevent < nb_muts ; myevent++ )
    {
      event = new ae_mutation( backup_file );
      dnareport->get_mutations()->add( event );
    }

    dnareport->compute_stats();
    _dna_replic_reports->add( dnareport );
  }
}








// =================================================================
//                             Destructors
// =================================================================
ae_replication_report::~ae_replication_report( void )
{
  _dna_replic_reports->erase( DELETE_OBJ );
  delete _dna_replic_reports;
}

// =================================================================
//                            Public Methods
// =================================================================

void ae_replication_report::write_to_backup( gzFile* backup_file )
{
  // Store individual identifiers and rank
  gzwrite( backup_file, &_index,        sizeof(_index)        );
  gzwrite( backup_file, &_rank,         sizeof(_rank)         );
  gzwrite( backup_file, &_parent_index, sizeof(_parent_index) );
  gzwrite( backup_file, &_donnor_index, sizeof(_donnor_index) );
  //~ printf( "  _index : %"PRId32"\n",         _index );
  //~ printf( "  _rank : %"PRId32"\n",          _rank );
  //~ printf( "  _parent_index : %"PRId32"\n",  _parent_index );
  //~ printf( "  _donnor_index : %"PRId32"\n",  _donnor_index );

  // Store parent and donnor's basic info
  gzwrite( backup_file, &_parent_metabolic_error, sizeof(_parent_metabolic_error) );
  gzwrite( backup_file, &_parent_secretion_error, sizeof(_parent_secretion_error) );
  gzwrite( backup_file, &_donnor_metabolic_error, sizeof(_donnor_metabolic_error) );
  gzwrite( backup_file, &_parent_genome_size,     sizeof(_parent_genome_size) );
  gzwrite( backup_file, &_donnor_genome_size,     sizeof(_donnor_genome_size) );
  //~ printf( "  _parent_metabolic_error : %lf\n", _parent_metabolic_error );
  //~ printf( "  _parent_secretion_error : %lf\n", _parent_secretion_error );
  //~ printf( "  _donnor_metabolic_error : %lf\n", _donnor_metabolic_error );
  //~ printf( "  _parent_genome_size : %"PRId32"\n", _parent_genome_size );
  //~ printf( "  _donnor_genome_size : %"PRId32"\n", _donnor_genome_size );
  
  // For each genetic unit, write the mutations and rearrangements undergone during replication
  int32_t nb_dna_replic_reports = _dna_replic_reports->get_nb_elts();
  gzwrite( backup_file, &nb_dna_replic_reports, sizeof(nb_dna_replic_reports) );
  //~ printf( "  nb_dna_replic_reports : %"PRId32"\n", nb_dna_replic_reports );

  ae_list_node*         report_node = _dna_replic_reports->get_first();
  ae_dna_replic_report* report      = NULL;
  
  while ( report_node != NULL )
  {
    report = (ae_dna_replic_report*) report_node->get_obj();
    
    // Store rearrangements
    int32_t nb_rears = report->get_nb_rearrangements();
    gzwrite( backup_file, &nb_rears, sizeof(nb_rears) );
    //~ printf( "  nb_rears : %"PRId32"\n", nb_rears );

    ae_list_node* rear_node  = report->get_rearrangements()->get_first();
    while ( rear_node != NULL )
    {
      ((ae_mutation*)rear_node->get_obj())->write_to_backup( backup_file );
      rear_node = rear_node->get_next();
    }
    
    // Store mutations
    int32_t nb_muts = report->get_nb_small_mutations();
    gzwrite( backup_file, &nb_muts, sizeof(nb_muts) );
    //~ printf( "  nb_muts : %"PRId32"\n", nb_muts );

    ae_list_node* mut_node  = report->get_mutations()->get_first();
    while ( mut_node != NULL )
    {
      ((ae_mutation*)mut_node->get_obj())->write_to_backup( backup_file );
      mut_node = mut_node->get_next();
    }
    
    report_node = report_node->get_next();
  }
}






// =================================================================
//                           Protected Methods
// =================================================================





// =================================================================
//                          Non inline accessors
// =================================================================

