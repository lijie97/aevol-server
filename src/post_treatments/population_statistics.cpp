/** \class
 *  \brief
 */

// =================================================================
//                              Libraries
// =================================================================
#include <stdlib.h>
#include <math.h>
#include <algorithm>

// =================================================================
//                            Project Files
// =================================================================
#include "population_statistics.h"


//##############################################################################
//                                                                             #
//                  Class population_statistics                                #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================


population_statistics::population_statistics( void )
{
    _output_file = fopen("robustness.out","w");
    _replication_output = NULL;
    
    _nb_children = 1000;
    _wanted_rank = -1;
    _wanted_index = -1;
    _details = false;

    _f_nu_pop = NULL;
    _reprod_proba = NULL;
    _fitness = NULL;
    
    _type = ROBUSTNESS;	
}


population_statistics::population_statistics( analysis_type type, int nb_children, char* output_dir, int wanted_rank, int wanted_index,  bool details)
{   
    _nb_children = nb_children;
    _wanted_rank = wanted_rank;
    _wanted_index = wanted_index;
    _details = details;
    _type = type;
    
    _f_nu_pop = NULL;
    _reprod_proba = NULL;
    _fitness = NULL;
    
    _output_file = NULL;
    _replication_output = NULL;
     
    char* file_name = NULL;
     
    if (output_dir==NULL) // default = .
    {
       if (_type == ROBUSTNESS) 
       {
        file_name = new char[ strlen("robustness.out") + 1 ];
        sprintf(file_name,"%s","robustness.out");
       }
       else 
       {
        file_name = new char[ strlen("evolvability.out") + 1];
        sprintf(file_name,"%s","evolvability.out");
       }
       
       _output_file = fopen(file_name,"w");
       if (_output_file==NULL) { fprintf(stderr, "Warning: could not open file %s.\n", file_name); }
       if (_details == true )
       {
        _replication_output=fopen("replications.out","w");
        if (_replication_output==NULL) { fprintf(stderr, "Warning: could not open file replications.out.\n"); }
       }
     }
     else
     {       
       if (_type == ROBUSTNESS) 
       {
        file_name = new char[ strlen(output_dir) + strlen("robustness.out") + 2];
        sprintf(file_name,"%s/%s",output_dir,"robustness.out");
       }
       else 
       {
        file_name = new char[ strlen(output_dir) + strlen("evolvability.out") + 2];
        sprintf(file_name,"%s/%s",output_dir,"evolvability.out");
       }
       
       char* replication_file_name = NULL;
       replication_file_name = new char[ strlen(output_dir) + strlen("replication.out") + 2];
       sprintf(replication_file_name,"%s/%s",output_dir,"replication.out");
       
       _output_file = fopen(file_name,"w");
       if ( _output_file == NULL ) { fprintf(stderr, "Warning: could not open file %s.\n", file_name); }
       if ( details == true )
       { 
        _replication_output=fopen(replication_file_name,"w");
        if (_replication_output==NULL) { fprintf(stderr, "Warning: could not open file %s.\n", replication_file_name); }
       }
       
       delete [] replication_file_name;
     }
     delete [] file_name;
     fflush( stderr );    
}


// =================================================================
//                             Destructors
// =================================================================

population_statistics::~population_statistics( void )
{
  if (_output_file != NULL) { fclose(_output_file); }
  if (_replication_output != NULL) { fclose(_replication_output); }
  
  delete [] _f_nu_pop;
  delete [] _reprod_proba;
  delete [] _fitness;
  
  //ae_common::clean();
}


// =================================================================
//                            Public Methods
// =================================================================

void population_statistics::compute_population_f_nu(ae_exp_manager* exp_manager)
{
  // ----------------------------------------
  //              Initializing....
  // ----------------------------------------
  
  // Load the simulation
  /*#ifndef __NO_X
    ae_exp_manager* exp_manager = new ae_exp_manager_X11();
  #else
    ae_exp_manager* exp_manager = new ae_exp_manager();
  #endif
  exp_manager->load_experiment( _exp_setup_file_name, _out_prof_file_name, environment_file_name, genomes_file_name, _sp_struct_file_name, true );*/

  
  // TODO: update this
  // load from lineage
  /* int last_backup_time = (generation / param::back_step) * param::back_step;
     char * command = new char[100];
     sprintf(command, "./lineage --begin %d --end %d --ind %d > /dev/null",last_backup_time,generation,indiv_id);
     printf("running %s...",command);
     fflush(NULL);
     system(command);
     printf(" done\n");
     delete [] command;
      
     ae_individual * initial_indiv = get_final_individual_using_dstory();
      
     // remove dstory output
     remove("dstory.bak.gz"); */

 
  // --------------------------------------------------------------
  //       Get genome(s) of interest and compute Fv
  // --------------------------------------------------------------
  
  _pop_size = exp_manager->get_nb_indivs();
  
  ae_individual* initial_indiv = NULL;
  ae_list_node<ae_individual*>* node = exp_manager->get_pop()->get_indivs()->get_last();
  int current_rank = 1;
  int current_index = -1;
  
  double neutral_or_better = -1;
  double exp_fv = -1;
  double th_fv;
  
  if (_f_nu_pop != NULL) { delete [] _f_nu_pop;}
  if (_fitness != NULL) { delete [] _fitness;}
  if (_reprod_proba != NULL) { delete [] _reprod_proba;}
  
  _f_nu_pop = new double[_pop_size];
  _fitness = new double[_pop_size];  
  _reprod_proba = new double [_pop_size];
  
  exp_manager->get_exp_s()->get_sel()->compute_prob_reprod();
  double* tmp_reprod = exp_manager->get_exp_s()->get_sel()->get_prob_reprod();
  
  // simply parse through individuals and keep those that are wanted
  while (node != NULL)
  {
    ae_individual* tmpind = (ae_individual*) node->get_obj();
    node = node->get_prev();

    if ( (_wanted_rank==-1) || (_wanted_rank == current_rank) )
    {
      current_index = tmpind->get_id();
      if ( (_wanted_index == -1) || (_wanted_index == current_index) )
      {
	
	// ------------------------------------
	//         Get initial individual
	// ------------------------------------
	      
      #ifdef __NO_X
	#ifndef __REGUL
	initial_indiv = new ae_individual(tmpind, current_index);
	#else
	initial_indiv = new ae_individual_R( (dynamic_cast<ae_individual_R*>(tmpind)), current_index );
	#endif
      #elif defined __X11
	#ifndef __REGUL
	initial_indiv = new ae_individual_X11( (dynamic_cast<ae_individual_X11*>(tmpind)), current_index );
	#else
	initial_indiv = new ae_individual_R_X11( (dynamic_cast<ae_individual_R_X11*>(tmpind)), current_index );
	#endif
      #endif

	initial_indiv->evaluate(exp_manager->get_env());
	//printf("found indiv %"PRId32" with fitness %le (rank %d)\n", tmpind->get_index_in_population(), tmpind->get_fitness(), current_rank);
	//printf("%d ", current_rank);
	fflush( stdout);
	
	_fitness[current_rank-1] = initial_indiv->get_fitness();
	_reprod_proba[current_rank-1] = tmp_reprod[_pop_size-current_rank];
	
	// ------------------------------------
	//              Compute Fv
	// ------------------------------------
	
	
		
	if ( _details == true && _type == ROBUSTNESS)
	{ 
	  // 'overloaded' function
	  //exp_fv = compute_experimental_fv( initial_indiv, &neutral_or_better );// TO DO
	  exp_fv = initial_indiv->compute_experimental_f_nu( _nb_children, &neutral_or_better );
	}
	else
	{
	  // aevol function
	  exp_fv = initial_indiv->compute_experimental_f_nu( _nb_children, &neutral_or_better );
	}
	
	_f_nu_pop[current_rank-1] = exp_fv;
	  
	if ( _type == ROBUSTNESS)
	{
	  th_fv = initial_indiv->compute_theoritical_f_nu();
	
	  // ------------------------------------
	  //            Write to file
	  // ------------------------------------ 
	  if ( _output_file != NULL )
	  {
	    fprintf( _output_file, "%d %d %le %le %le %le %le %"PRId32"\n",
		   current_rank, current_index, initial_indiv->get_fitness(),
		   _reprod_proba[current_rank-1],
		   exp_fv, neutral_or_better, th_fv, _pop_size );
	  // fflush(fv_output);
	  
	  }
	}
	      
	delete initial_indiv;
      }
    }
    current_rank++;
  }
  //ae_common::clean();
  //sdelete exp_manager;
}


void population_statistics::compute_evolvability_stats(int32_t num_gener)
{
  printf("\n\nCompute evolvability stats\n");
  double variation = 0;
  double population_variability = 0; 
  
  double best_variability = 1. - _f_nu_pop[0];
  
  int fitness_number = 1;
  int quasi_species_size = 1;
  
  char* file_name = NULL;
  FILE* f_nu_file;
  file_name = new char[ strlen("F_nu/gen_.out") + 10];
  sprintf(file_name,"%s%d%s","F_nu/gen_",num_gener,".out");
  
  f_nu_file = fopen(file_name,"w");
  delete [] file_name;
  
  population_variability = _f_nu_pop[0] * _reprod_proba[0];
  fprintf( f_nu_file, "%le %le %le\n",_f_nu_pop[0],_reprod_proba[0]*_pop_size,_fitness[0]);
  for (int i = 1; i < _pop_size; i++)
  {
    if ( fabs(_fitness[i] - _fitness[i-1]) > 1e-10*std::max(_fitness[i],_fitness[i-1]))
    {
      fitness_number += 1;
    }
    else if ( fitness_number == 1)
    {
      quasi_species_size += 1;
    }
    
    fprintf( f_nu_file, "%le %le %le\n",_f_nu_pop[i],_reprod_proba[i]*_pop_size,_fitness[i]);
    population_variability += _f_nu_pop[i] * _reprod_proba[i];
  }
  fflush( f_nu_file );
  
  population_variability = 1. - population_variability;
  variation = (static_cast<double>(fitness_number)) / (static_cast<double>(_pop_size));
  
  printf("\n Evolvability stats : \n \tgeneration : %"PRId32"\n \tpopulation size : %"PRId32"\n \tvariation : %le\n \tpopulation variability : %le\n \tbest f nu : %le\n \tbest variability : %le\n \tbest quasi species size : %d\n \tfitness number : %d\n", num_gener, 
	   _pop_size, variation, population_variability, _f_nu_pop[0], best_variability, 
	   quasi_species_size, fitness_number);
  
  fprintf( _output_file, "%"PRId32" %"PRId32" %le %le %le %d\n", num_gener, 
	   _pop_size, variation, population_variability, best_variability, 
	   quasi_species_size);
  fflush( _output_file );
}

// =================================================================
//                            Protected Methods
// =================================================================

// function copied from ae_individual's computation of experimental fv.
// In addition, it provides and prints information about replications
// TO DO : Change this function for adaptation to the new version
double population_statistics::compute_experimental_fv( ae_individual* indiv, double* neutral_or_better)
{
  double initial_fitness = indiv->get_fitness();
  int32_t index = indiv->get_id();
  double Fv = 0;
  if ( neutral_or_better != NULL) *neutral_or_better = 0;
  
  // ------------------------------------------
  //      Simulate fitness degradation
  // ------------------------------------------
  double childfit = 0;
  double * fitnesses = new double [_nb_children];
  int * nb_aff_genes = new int [_nb_children];
	      
  // reproduce this individual to create 'nb_children' children 
  ae_individual * child = NULL;
  for (int i = 0; i < _nb_children; i++)
  {
    //child = ae_common::sim->get_pop()->do_replication(indiv, index); 
    childfit = child->get_fitness(); // child is automatically evaluated
    
    if (_replication_output != NULL) // only interesting if we put that into a file...
    {
      fitnesses[i]    = childfit;
      nb_aff_genes[i] = count_affected_genes(indiv, child);
    }
    
    //count neutral offspring
    if ( fabs(initial_fitness - childfit) < 1e-10*std::max(initial_fitness,childfit) )
    { 
      Fv += 1;
      if ( neutral_or_better != NULL ) *neutral_or_better += 1;
    }
    else if ( (neutral_or_better != NULL) && (childfit > initial_fitness) )
    {
       *neutral_or_better += 1;
    }
    delete child;
  }
  
  //compute Fv
  Fv /= _nb_children;
  if ( neutral_or_better != NULL ) *neutral_or_better /= _nb_children;
  	      
  // ------------------------------------------
  //      print statistics to file
  // ------------------------------------------
  if ( _replication_output != NULL )
  {
    //general statistics of the individual
    fprintf( _replication_output, "#%"PRId32" %"PRId32" %"PRId32" %d %le",
	    indiv->get_id(),
	    indiv->get_total_genome_size(),
	    indiv->get_nb_functional_genes(),
	    (int) indiv->get_overall_size_coding_RNAs(),
	    indiv->get_fitness());
    //statistics of replication
    print_replication_stats( indiv, fitnesses, nb_aff_genes );    
  }  

  if ( fitnesses != NULL ) delete [] fitnesses;
  if ( nb_aff_genes != NULL ) delete [] nb_aff_genes;

  return Fv;
}


// count how many proteins were modified after replication
int population_statistics::count_affected_genes( ae_individual* parent, ae_individual* child )
{
  // ------------------------------------------------------------
  //       list all functional proteins of the child
  // ------------------------------------------------------------
  
  // Copy protein list
  ae_list<ae_protein*>* child_protein_list = new ae_list<ae_protein*>( *(child->get_protein_list()) );
  // Delete all nodes containing degenerated proteins
  ae_list_node<ae_protein*>* current_prot_node = child_protein_list->get_first();
  ae_list_node<ae_protein*>* next_prot_node = NULL;
  ae_protein* current_prot = NULL;
  while ( current_prot_node != NULL )
  {
    current_prot = current_prot_node->get_obj();
    next_prot_node = current_prot_node->get_next();
    if ( current_prot->get_is_functional() == false ) 
    {
      //delete the node but not the protein
      child_protein_list->remove( current_prot_node, true, false );
    }
    current_prot_node = next_prot_node;
  }
  
  // ----------------------------------------------------------------
  //    for each functional protein of the parent, search for
  //    the same one in the child's proteins and pop it out of the  
  //    list if found
  // ----------------------------------------------------------------
  ae_list_node<ae_protein*>* parent_prot_node = parent->get_protein_list()->get_first();
  ae_protein* parent_prot = NULL;
  ae_list_node<ae_protein*>* child_prot_node = NULL;
  ae_protein* child_prot = NULL;
  bool found = false;

  int nb_proteins_left = 0;

  // parse parent protein list
  while ( parent_prot_node != NULL )
  {
    parent_prot = parent_prot_node->get_obj();
    if ( parent_prot->get_is_functional() == true )
    {
      found = false;
      // parse remaining protein list of the child
      child_prot_node = child_protein_list->get_first();
      while ( child_prot_node != NULL )
      {
	child_prot = (ae_protein*) child_prot_node->get_obj();
	// compare the two proteins
	if ( (fabs(child_prot->get_mean() - parent_prot->get_mean()) < 1e-20)
	     && (fabs(child_prot->get_width() - parent_prot->get_width()) < 1e-20)
	     && (fabs(child_prot->get_height() - parent_prot->get_height()) < 1e-20) )
	{
	  found = true;
	  // delete the node but not the protein
	  child_protein_list->remove( child_prot_node, true, false);
	  child_prot_node=NULL;
	}
	else
	{
	  child_prot_node=child_prot_node->get_next();
	}
      }
      // ------------------------------------------------------------------------------------
      //  if the protein was not found, it means that one of the genes were modified in child
      // ------------------------------------------------------------------------------------
      if (found==false) nb_proteins_left++;
    }
    parent_prot_node = parent_prot_node->get_next();
  }
  
  // --------------------------------------------------------------
  //   with low probability, some genes [dis]appeared
  //   therefore, the number of proteins of the parent not
  //   retreived in the child might not be equal to the number
  //   of proteins left in child. We keep the biggest element,
  //   neglecting possible compensations
  // --------------------------------------------------------------
  //~ printf("%d %d \n",child_protein_list->get_nb_elts(), nb_proteins_left);
  int result = (child_protein_list->get_nb_elts() > nb_proteins_left) ? child_protein_list->get_nb_elts() : nb_proteins_left;

  delete child_protein_list;

  return result;
}


// print stats about replications of the individual
void population_statistics::print_replication_stats( ae_individual* initial_indiv, double* fitnesses, int* nb_aff_genes )
{
  if ( _replication_output == NULL ) return;
  
  double fitness_initial_indiv = initial_indiv->get_fitness();

  // --------------------------------------------------------------
  //              compute statistics on fitness variations
  // --------------------------------------------------------------
  double average_fit = 0.0;
  double sd = 0.0;
  double prop_better_than_parent = 0;
  double prop_same_as_parent = 0;
  double prop_no_genes_affected = 0;
  double avg_nb_aff_genes = 0.0;

  for (int i = 0; i < _nb_children; i++)
  {
    average_fit += fitnesses[i];
    avg_nb_aff_genes += nb_aff_genes[i];
    sd += fitnesses[i]*fitnesses[i];
    if ( fabs(fitnesses[i] - fitness_initial_indiv) < 1e-10*std::max(fitnesses[i],fitness_initial_indiv) ) { prop_same_as_parent += 1; }
    else if ( fitnesses[i] > fitness_initial_indiv ) { prop_better_than_parent += 1; }
    if ( nb_aff_genes[i] == 0 ) { prop_no_genes_affected += 1; }
  }
  average_fit /= _nb_children;
  avg_nb_aff_genes /= _nb_children;
  prop_no_genes_affected /= _nb_children;
  sd /= _nb_children;
  sd -= (average_fit)*(average_fit);
  sd = sqrt(sd);
  prop_better_than_parent /= _nb_children;
  prop_same_as_parent /= _nb_children;

  // --------------------------------------------------------------
  //         print statistics on current line
  // --------------------------------------------------------------
  fprintf(_replication_output, " %le %le %le %le %le %le ", average_fit, sd, prop_same_as_parent, prop_better_than_parent, prop_no_genes_affected, avg_nb_aff_genes);
  fprintf(_replication_output, "\n");
  
  // --------------------------------------------------------------
  //         print fitness statistics for each child
  // --------------------------------------------------------------
  for (int i = 0; i < _nb_children; i++)
  {
    fprintf(_replication_output, "%le %d \n", fitnesses[i], nb_aff_genes[i]);
  }
}


/*
  ae_individual * get_final_individual_using_dstory()
  {
  // warning: delete the genome when you don't need it anymore

  gzFile dstory_file = gzopen("dstory.bak.gz","r");

  int begin_gener, end_gener, final_index;
  gzread(dstory_file, &begin_gener, sizeof(int));
  gzread(dstory_file, &end_gener, sizeof(int));
  gzread(dstory_file, &final_index, sizeof(int));

  // ----------------------
  // getting initial genome
  // ----------------------

  int size;
  gzread(dstory_file, &size, sizeof(int));
  char * gen_str = new char[size+1];
  gzread(dstory_file, gen_str, sizeof(char)*(size+1));
  ae_individual * indiv = new ae_individual(gen_str, -1, -1, -1); // phenotype computed
  delete [] gen_str;
  

  // ----------------------------------------------------
  //  rebuild successive genomes and write in 
  //  stat_file at each time step
  // ----------------------------------------------------

  int nb_muts;
  f_list mut_list;
  f_list_node * m_node;
  int type = -1;
  int at = -1;
  int begin_seg = -1;
  int end_seg = -1;
  int seq_length = 0;
  char * seq = NULL;
  ae_mutation * mut = NULL;

  for(int t = begin_gener + 1; t <= end_gener; t++) // t = generation to be done
  {

  // get mutation information

  gzread(dstory_file, &nb_muts, sizeof(int));
  for (int m = 0; m < nb_muts; m++)
  {
  gzread(dstory_file, &type, sizeof(int));
  gzread(dstory_file, &at, sizeof(int));
  gzread(dstory_file, &begin_seg, sizeof(int));
  gzread(dstory_file, &end_seg, sizeof(int));
  gzread(dstory_file, &seq_length, sizeof(int));
  if (seq_length != 0)
  {
  seq = new char[seq_length + 1];
  gzread(dstory_file, seq, sizeof(char)*(seq_length+1));
  }
  else
  {
  seq = NULL;
  }
  mut = new ae_mutation(type, at, begin_seg, end_seg, seq);
  if (seq!=NULL) delete[] seq;
  mut_list.add(mut);
  }  

  // execute mutations

  m_node = mut_list._first;
  while(m_node != NULL)
  {
  mut = (ae_mutation *) m_node->_obj;
  mut->execute(indiv->get_genome());  
  m_node = m_node->_next;
  }

  // prepare next time step

  mut_list.erase(true);
  } 

  gzclose(dstory_file);
  
  indiv->compute_phenotype(true);
  return indiv;
  }
*/
