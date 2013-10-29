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
ae_replication_report::ae_replication_report( ae_individual * indiv, ae_individual * parent, ae_individual * donor /* = NULL */ )
{
  _indiv = indiv;
  
  _id   = indiv->get_id();
  _rank = indiv->get_rank();
  
  _parent_id = parent->get_id();
  // _donor_id is set further down
    
  _genome_size        = 0;
  _metabolic_error    = 0.0;
  _nb_genes_activ     = 0;
  _nb_genes_inhib     = 0;
  _nb_non_fun_genes   = 0;
  _nb_coding_RNAs     = 0;
  _nb_non_coding_RNAs = 0;
  
  _parent_metabolic_error = parent->get_dist_to_target_by_feature( METABOLISM );
  _parent_secretion_error = parent->get_dist_to_target_by_feature( SECRETION );
  _parent_genome_size     = parent->get_total_genome_size();
  _mean_align_score       = 0.0;
  
  if ( donor == NULL )
  {
    _donor_id               = -1;
    _donor_metabolic_error  = 0.0;
    _donor_secretion_error	= 0.0;
    _donor_genome_size      = 0;
  }
  else
  {
    _donor_id              = donor->get_id();
    _donor_metabolic_error = donor->get_dist_to_target_by_feature( METABOLISM );
    _donor_secretion_error = donor->get_dist_to_target_by_feature( SECRETION );
    _donor_genome_size     = donor->get_total_genome_size();
  }
    
  _HT_ins 							            = false;
  _HT_ins_sense 					          = DIRECT;
  _HT_ins_donor_id 					        = 0;
  _HT_ins_alignment_1_donor_pos_1 	= 0;
  _HT_ins_alignment_1_donor_pos_2 	= 0; 
  _HT_ins_alignment_1_score 		    = 0;
  _HT_ins_alignment_2_ind_pos 		  = 0;
  _HT_ins_alignment_2_donor_pos 	  = 0; 
  _HT_ins_alignment_2_score 		    = 0;   
  _HT_repl 							            = false;
  _HT_repl_sense 					          = DIRECT;
  _HT_repl_donor_id 				        = 0;
  _HT_repl_alignment_1_ind_pos 		  = 0;
  _HT_repl_alignment_1_donor_pos 	= 0; 
  _HT_repl_alignment_1_score 		    = 0;
  _HT_repl_alignment_2_ind_pos 		  = 0;
  _HT_repl_alignment_2_donor_pos 	= 0; 
  _HT_repl_alignment_2_score 		    = 0;  
  
  _dna_replic_reports = new ae_list<ae_dna_replic_report*>();
}


// Creates an independent copy of the original report
ae_replication_report::ae_replication_report( const ae_replication_report &model )
{
  _parent_id  = model._parent_id;
  _donor_id   = model._donor_id;
  
  _id   = model._id;
  _rank = model._rank;
    
  _genome_size        = model._genome_size;
  _metabolic_error    = model._metabolic_error;
  _nb_genes_activ     = model._nb_genes_activ;
  _nb_genes_inhib     = model._nb_genes_inhib;
  _nb_non_fun_genes   = model._nb_non_fun_genes;
  _nb_coding_RNAs     = model._nb_coding_RNAs;
  _nb_non_coding_RNAs = model._nb_non_coding_RNAs;

  _parent_metabolic_error = model._parent_metabolic_error;
  _parent_secretion_error = model._parent_secretion_error;
  _donor_metabolic_error  = model._donor_metabolic_error;
  _donor_secretion_error  = model._donor_secretion_error;
  _parent_genome_size     = model._parent_genome_size;
  _donor_genome_size      = model._donor_genome_size;
  _mean_align_score       = model._mean_align_score;
  
  _HT_ins 							            = model._HT_ins;
  _HT_ins_sense 					          = model._HT_ins_sense;
  _HT_ins_donor_id 					        = model._HT_ins_donor_id;
  _HT_ins_alignment_1_donor_pos_1 	= model._HT_ins_alignment_1_donor_pos_1;
  _HT_ins_alignment_1_donor_pos_2 	= model._HT_ins_alignment_1_donor_pos_2; 
  _HT_ins_alignment_1_score 	    	= model._HT_ins_alignment_1_score;
  _HT_ins_alignment_2_ind_pos 	  	= model._HT_ins_alignment_2_ind_pos;
  _HT_ins_alignment_2_donor_pos  	  = model._HT_ins_alignment_2_donor_pos; 
  _HT_ins_alignment_2_score 		    = model._HT_ins_alignment_2_score;   
  _HT_repl 						            	= model._HT_repl;
  _HT_repl_sense 				          	= model._HT_repl_sense;
  _HT_repl_donor_id 	        			= model._HT_repl_donor_id;
  _HT_repl_alignment_1_ind_pos 		  = model._HT_repl_alignment_1_ind_pos;
  _HT_repl_alignment_1_donor_pos 	  = model._HT_repl_alignment_1_donor_pos; 
  _HT_repl_alignment_1_score 		    = model._HT_repl_alignment_1_score;
  _HT_repl_alignment_2_ind_pos 	  	= model._HT_repl_alignment_2_ind_pos;
  _HT_repl_alignment_2_donor_pos 	  = model._HT_repl_alignment_2_donor_pos; 
  _HT_repl_alignment_2_score 	    	= model._HT_repl_alignment_2_score;  
  
  _dna_replic_reports = new ae_list<ae_dna_replic_report*>();

  ae_list_node<ae_dna_replic_report*>* node = (model._dna_replic_reports)->get_first();
  ae_dna_replic_report * dnarep = NULL;
  while ( node != NULL )
  {
    dnarep = node->get_obj();
    _dna_replic_reports->add( new ae_dna_replic_report(*dnarep) );
    node = node->get_next();
  }
}


ae_replication_report::ae_replication_report( gzFile tree_file, ae_individual * indiv )
{
  _indiv = indiv;
    
  gzread( tree_file, &_id,        sizeof(_id)         );
  gzread( tree_file, &_rank,      sizeof(_rank)       );
  gzread( tree_file, &_parent_id, sizeof(_parent_id)  );
  gzread( tree_file, &_donor_id,  sizeof(_donor_id)   );
  
  gzread( tree_file, &_genome_size,         sizeof(_genome_size) );
  gzread( tree_file, &_metabolic_error,     sizeof(_metabolic_error) );
  gzread( tree_file, &_nb_genes_activ,      sizeof(_nb_genes_activ) );
  gzread( tree_file, &_nb_genes_inhib,      sizeof(_nb_genes_inhib) );
  gzread( tree_file, &_nb_non_fun_genes,    sizeof(_nb_non_fun_genes) );
  gzread( tree_file, &_nb_coding_RNAs,      sizeof(_nb_coding_RNAs) );
  gzread( tree_file, &_nb_non_coding_RNAs,  sizeof(_nb_non_coding_RNAs) );
  
  gzread( tree_file, &_HT_ins,         					        sizeof(_HT_ins) );
  gzread( tree_file, &_HT_ins_sense,         			      sizeof(_HT_ins_sense) );
  gzread( tree_file, &_HT_ins_donor_id,         		    sizeof(_HT_ins_donor_id) );
  gzread( tree_file, &_HT_ins_alignment_1_donor_pos_1,  sizeof(_HT_ins_alignment_1_donor_pos_1) );
  gzread( tree_file, &_HT_ins_alignment_1_donor_pos_2,  sizeof(_HT_ins_alignment_1_donor_pos_2) ); 
  gzread( tree_file, &_HT_ins_alignment_1_score,        sizeof(_HT_ins_alignment_1_score) );
  gzread( tree_file, &_HT_ins_alignment_2_ind_pos,      sizeof(_HT_ins_alignment_2_ind_pos) );
  gzread( tree_file, &_HT_ins_alignment_2_donor_pos,    sizeof(_HT_ins_alignment_2_donor_pos) ); 
  gzread( tree_file, &_HT_ins_alignment_2_score,        sizeof(_HT_ins_alignment_2_score) );   
  gzread( tree_file, &_HT_repl,         				        sizeof(_HT_repl) );
  gzread( tree_file, &_HT_repl_sense,         			    sizeof(_HT_repl_sense) );
  gzread( tree_file, &_HT_repl_donor_id,         		    sizeof(_HT_repl_donor_id) );
  gzread( tree_file, &_HT_repl_alignment_1_ind_pos,     sizeof(_HT_repl_alignment_1_ind_pos) );
  gzread( tree_file, &_HT_repl_alignment_1_donor_pos,   sizeof(_HT_repl_alignment_1_donor_pos) ); 
  gzread( tree_file, &_HT_repl_alignment_1_score,       sizeof(_HT_repl_alignment_1_score) );
  gzread( tree_file, &_HT_repl_alignment_2_ind_pos,     sizeof(_HT_repl_alignment_2_ind_pos) );
  gzread( tree_file, &_HT_repl_alignment_2_donor_pos,   sizeof(_HT_repl_alignment_2_donor_pos) ); 
  gzread( tree_file, &_HT_repl_alignment_2_score,       sizeof(_HT_repl_alignment_2_score) );    
  
  int32_t nb_dna_replic_reports;
  gzread( tree_file, &nb_dna_replic_reports, sizeof(nb_dna_replic_reports) );
  //~ printf( "  nb_dna_replic_reports : %"PRId32"\n", nb_dna_replic_reports );
  
  _dna_replic_reports = new ae_list<ae_dna_replic_report*>();

  int32_t mydnareport, myevent;
  int32_t nb_rears, nb_muts;
  ae_dna_replic_report * dnareport = NULL;
  ae_mutation * event = NULL;

  for ( mydnareport = 0 ; mydnareport < nb_dna_replic_reports ; mydnareport++ )
  {
    dnareport = new ae_dna_replic_report();
    
    gzread( tree_file, &nb_rears, sizeof(nb_rears) );
    //~ printf( "  nb_rears : %"PRId32"\n", nb_rears );
    for ( myevent  = 0 ; myevent < nb_rears ; myevent++ )
    {
      event = new ae_mutation( tree_file );
      dnareport->get_rearrangements()->add( event );
    }
    
    gzread( tree_file, &nb_muts, sizeof(nb_muts) );
    //~ printf( "  nb_rears : %"PRId32"\n", nb_rears );
    for(myevent  = 0 ; myevent < nb_muts ; myevent++ )
    {
      event = new ae_mutation( tree_file );
      dnareport->get_mutations()->add( event );
    }

    dnareport->compute_stats();
    _dna_replic_reports->add( dnareport );
  }
  
  _parent_metabolic_error = -1;
  _parent_secretion_error = -1;
  _donor_metabolic_error  = -1;
  _parent_genome_size     = -1;
  _donor_genome_size      = -1;
  _mean_align_score       = 0.0;
}


// =================================================================
//                             Destructors
// =================================================================
ae_replication_report::~ae_replication_report( void )
{
  _dna_replic_reports->erase( true );
  delete _dna_replic_reports;
}

// =================================================================
//                            Public Methods
// =================================================================
/**
 * Method called at the end of the replication. Actions such as finalize the calculation of average values can be done here.
 */
void ae_replication_report::signal_end_of_replication( void )
{
  // Retreive data from the individual
  _genome_size        = _indiv->get_total_genome_size();
  _metabolic_error    = _indiv->get_dist_to_target_by_feature( METABOLISM );
  _nb_genes_activ     = _indiv->get_nb_genes_activ();
  _nb_genes_inhib     = _indiv->get_nb_genes_inhib();
  _nb_non_fun_genes   = _indiv->get_nb_functional_genes();
  _nb_coding_RNAs     = _indiv->get_nb_coding_RNAs();
  _nb_non_coding_RNAs = _indiv->get_nb_non_coding_RNAs();
  
  
  // Compute the mean alignment score
  int32_t nb_align = 0;
  ae_list_node<ae_dna_replic_report*>* dna_rep_node  = _dna_replic_reports->get_first();
  ae_dna_replic_report* dna_rep = NULL;
  ae_list_node<ae_mutation*>* rear_node     = NULL;
  ae_mutation* rear = NULL;
  ae_mutation_type        mut_type;
  int32_t*  int32_trash = new int32_t;
  bool*     bool_trash  = new bool;
  int16_t*  align_scores = new int16_t[2];
  while ( dna_rep_node != NULL )
  {
    dna_rep = dna_rep_node->get_obj();
    
    nb_align += dna_rep->get_nb_duplications();
    nb_align += dna_rep->get_nb_deletions();
    nb_align += 2 * dna_rep->get_nb_translocations();
    nb_align += dna_rep->get_nb_inversions();
    
    rear_node = dna_rep->get_rearrangements()->get_first();
    while ( rear_node != NULL )
    {
      rear = rear_node->get_obj();
      
      mut_type = rear->get_mut_type();
      
      switch( mut_type )
      {
        case DUPL:
        {
          rear->get_infos_duplication( int32_trash, int32_trash, int32_trash, &align_scores[0] );
          _mean_align_score += align_scores[0];
          break;
        }
        case DEL:
        {
          rear->get_infos_deletion( int32_trash, int32_trash, &align_scores[0] );
          _mean_align_score += align_scores[0];
          break;
        }
        case TRANS:
        {
          rear->get_infos_translocation( int32_trash, int32_trash, int32_trash, int32_trash, bool_trash, &align_scores[0], &align_scores[1] );
          _mean_align_score += align_scores[0] + align_scores[1];
          break;
        }
        case INV:
        {
          rear->get_infos_inversion( int32_trash, int32_trash, &align_scores[0] );
          _mean_align_score += align_scores[0];
          break;
        }
        default:
        {
          fprintf( stderr, "ERROR, invalid mutation type \"%d\" in file %s:%d.\n", mut_type, __FILE__, __LINE__ );
          exit( EXIT_FAILURE );
        }
      }
      
      rear_node = rear_node->get_next();
    }
    
    dna_rep_node = dna_rep_node->get_next();
  }
  
  if ( nb_align != 0 )
  {
    _mean_align_score /= nb_align;
  }
  else
  {
    assert( _mean_align_score == 0.0 );
  }
  
  delete int32_trash;
  delete bool_trash;
  delete [] align_scores;
}

void ae_replication_report::write_to_tree_file( gzFile tree_file ) const
{
  // Store individual identifiers and rank
  gzwrite( tree_file, &_id,         sizeof(_id)         );
  gzwrite( tree_file, &_rank,       sizeof(_rank)       );
  gzwrite( tree_file, &_parent_id,  sizeof(_parent_id)  );
  gzwrite( tree_file, &_donor_id,   sizeof(_donor_id)   );
  
  gzwrite( tree_file, &_genome_size,         sizeof(_genome_size) );
  gzwrite( tree_file, &_metabolic_error,     sizeof(_metabolic_error) );
  gzwrite( tree_file, &_nb_genes_activ,      sizeof(_nb_genes_activ) );
  gzwrite( tree_file, &_nb_genes_inhib,      sizeof(_nb_genes_inhib) );
  gzwrite( tree_file, &_nb_non_fun_genes,    sizeof(_nb_non_fun_genes) );
  gzwrite( tree_file, &_nb_coding_RNAs,      sizeof(_nb_coding_RNAs) );
  gzwrite( tree_file, &_nb_non_coding_RNAs,  sizeof(_nb_non_coding_RNAs) );
  
  gzwrite( tree_file, &_HT_ins,         					        sizeof(_HT_ins) );
  gzwrite( tree_file, &_HT_ins_sense,         				    sizeof(_HT_ins_sense) );
  gzwrite( tree_file, &_HT_ins_donor_id,         			    sizeof(_HT_ins_donor_id) );
  gzwrite( tree_file, &_HT_ins_alignment_1_donor_pos_1,	  sizeof(_HT_ins_alignment_1_donor_pos_1) );
  gzwrite( tree_file, &_HT_ins_alignment_1_donor_pos_2, 	sizeof(_HT_ins_alignment_1_donor_pos_2) ); 
  gzwrite( tree_file, &_HT_ins_alignment_1_score,        	sizeof(_HT_ins_alignment_1_score) );
  gzwrite( tree_file, &_HT_ins_alignment_2_ind_pos,      	sizeof(_HT_ins_alignment_2_ind_pos) );
  gzwrite( tree_file, &_HT_ins_alignment_2_donor_pos,   	sizeof(_HT_ins_alignment_2_donor_pos) ); 
  gzwrite( tree_file, &_HT_ins_alignment_2_score,        	sizeof(_HT_ins_alignment_2_score) );   
  gzwrite( tree_file, &_HT_repl,         					        sizeof(_HT_repl) );
  gzwrite( tree_file, &_HT_repl_sense,         				    sizeof(_HT_repl_sense) );
  gzwrite( tree_file, &_HT_repl_donor_id,         			  sizeof(_HT_repl_donor_id) );
  gzwrite( tree_file, &_HT_repl_alignment_1_ind_pos,     	sizeof(_HT_repl_alignment_1_ind_pos) );
  gzwrite( tree_file, &_HT_repl_alignment_1_donor_pos,  	sizeof(_HT_repl_alignment_1_donor_pos) ); 
  gzwrite( tree_file, &_HT_repl_alignment_1_score,       	sizeof(_HT_repl_alignment_1_score) );
  gzwrite( tree_file, &_HT_repl_alignment_2_ind_pos,     	sizeof(_HT_repl_alignment_2_ind_pos) );
  gzwrite( tree_file, &_HT_repl_alignment_2_donor_pos,  	sizeof(_HT_repl_alignment_2_donor_pos) ); 
  gzwrite( tree_file, &_HT_repl_alignment_2_score,       	sizeof(_HT_repl_alignment_2_score) );    
  
  
  // For each genetic unit, write the mutations and rearrangements undergone during replication
  int32_t nb_dna_replic_reports = _dna_replic_reports->get_nb_elts();
  gzwrite( tree_file, &nb_dna_replic_reports, sizeof(nb_dna_replic_reports) );
  //~ printf( "  nb_dna_replic_reports : %"PRId32"\n", nb_dna_replic_reports );

  ae_list_node<ae_dna_replic_report*>* report_node = _dna_replic_reports->get_first();
  ae_dna_replic_report* report = NULL;
  
  while ( report_node != NULL )
  {
    report = report_node->get_obj();
    
    // Store rearrangements
    int32_t nb_rears = report->get_nb_rearrangements();
    gzwrite( tree_file, &nb_rears, sizeof(nb_rears) );
    //~ printf( "  nb_rears : %"PRId32"\n", nb_rears );

    ae_list_node<ae_mutation*>* rear_node = report->get_rearrangements()->get_first();
    while ( rear_node != NULL )
    {
      rear_node->get_obj()->save( tree_file );
      rear_node = rear_node->get_next();
    }
    
    // Store mutations
    int32_t nb_muts = report->get_nb_small_mutations();
    gzwrite( tree_file, &nb_muts, sizeof(nb_muts) );
    //~ printf( "  nb_muts : %"PRId32"\n", nb_muts );

    ae_list_node<ae_mutation*>* mut_node = report->get_mutations()->get_first();
    while ( mut_node != NULL )
    {
      mut_node->get_obj()->save( tree_file );
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

