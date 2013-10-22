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
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <err.h>
#include <errno.h>
#include <sys/stat.h>

// =================================================================
//                            Project Files
// =================================================================
#include "population_statistics.h"

#define STATS_DIR                   "stats"
#define ROBUSTNESS_FILE             STATS_DIR"/robustness_%06"PRId32".out"
#define REPLICATION_FILE            STATS_DIR"/replication_%06"PRId32".out"
#define POP_STAT_FILE               STATS_DIR"/global_pop_stats.out"
#define POP_STAT_DIR                STATS_DIR"/pop_stats/"
#define POP_STATE_BASE              "gen_%06"PRId32
#define POP_STATE_ROBUSTNESS_FILE   POP_STAT_DIR"/robustness_%06"PRId32".out"
#define POP_STATE_REPLICATION_FILE  POP_STAT_DIR"/replication_%06"PRId32".out"




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

/*!
  \brief Default constructor of population_statitistics class

  Initialisation of default values _nb_children (number of replications to compute proportion of neutral offsrpings), _wanted_rank and _wanted_index (indications of 
  individual of whom we want more replication information), and anaylsis type
  By default, the number of replications is 1,000 and the wanted individual is the best one (rank = 1) and the analysis is only one generation).
*/
population_statistics::population_statistics( void )
{
    _output_file = NULL;
    _replication_file = NULL;
    _robustness_file = NULL;
    
    _nb_children = 1000;
    _wanted_rank = 1;
    _wanted_index = -1;

    _f_nu_pop = NULL;
    _reprod_proba = NULL;
    _fitness = NULL;
    
    _type = ONE_GENERATION;	
}

/*!
  \brief Constructor of population_statitistics class with non-default values

  Initialisation of _nb_children, _wanted_rank and _wanted_index and _type with the given values
  
  \param type analysis type: ONE_GENERATION or MULTIPLE_GENERATIONS
  \param nb_children  number of replications needed to compute the proportion of neutral offsprings of each individual
  \param wanted_rank rank of the individual of whom we want more replication information
  \param wanted_index index of the individual of whom we want more replication information
*/
population_statistics::population_statistics( analysis_type type, int32_t nb_children, int32_t wanted_rank, int32_t wanted_index)
{   
    _nb_children = nb_children;
    _wanted_rank = wanted_rank;
    _wanted_index = wanted_index;
    _type = type;
    
    _f_nu_pop = NULL;
    _reprod_proba = NULL;
    _fitness = NULL;
    
    _output_file = NULL;
    _replication_file = NULL;
    _robustness_file = NULL;
     
    
    if (_type != ONE_GENERATION) 
    {
      _output_file = fopen(POP_STAT_FILE,"w");
      int status = mkdir( POP_STAT_DIR, 0755 );
      if ( (status == -1) && (errno != EEXIST) )
      {
        err( EXIT_FAILURE, POP_STAT_DIR );
      }
      assert(_output_file!=NULL);
    }
    
    
    fflush( stderr );    
}


// =================================================================
//                             Destructors
// =================================================================

/*!
  \brief Destructor of population_statistics class
*/
population_statistics::~population_statistics( void )
{
  if (_output_file != NULL) { fclose(_output_file); }
  
  delete [] _f_nu_pop;
  delete [] _reprod_proba;
  delete [] _fitness;
}


// =================================================================
//                            Public Methods
// =================================================================

/*!
  \brief Compute statistics of reproduction at population level at generation num_gener

  Compute at generation num_gener for each individual:
  * Fitness
  * Metabolic error
  * Genome size
  * Functional gene number
  * Reproduction probability
  * Proportion of neutral offsprings
  * Proportion of beneficial offsprings
  * Proportion of deleterious offsprings
  * Theoretical proportion of deleterious offsprings
  * Fitness mean of offsprings
  * Fitness variance of offsprings
  * Genome size mean of offsprings
  * Genome size variance of offsprings
  * Functional gene number mean of offsprings;
  * 17. Functional gene number variance of offsprings
  Write this information in a file robustness_numgener.out 
  
  Compute for each offspring of a choosen individual (by index or rank):
  * Fitness
  * Metabolic error
  * Genome size
  * Functional gene number
  * Number of coding bases
  * Number of transcribed but not translated bases
  * Number of non transcribed bases
  Write this information in a file replication_numgener.out 
  
  Complete _f_nu_pop, _reprod_proba, _fitness with proportion of neutral offspring, reproduction probability and fitness of each individual at num_gener
  
  \param exp_manager  current exp_manager
  \param num_gener    current generation number
*/
void population_statistics::compute_reproduction_stats(ae_exp_manager* exp_manager, int32_t num_gener)
{
  // ----------------------------------------
  //              Open output files
  // ----------------------------------------

  char* robustness_file_name  = new char[255];
  char* replication_file_name = new char[255];
  if (_type == ONE_GENERATION) 
  {
    sprintf( robustness_file_name,  ROBUSTNESS_FILE,  num_gener );
    sprintf( replication_file_name, REPLICATION_FILE, num_gener );
  }
  else
  {
    sprintf( robustness_file_name,  POP_STATE_ROBUSTNESS_FILE,  num_gener );
    sprintf( replication_file_name, POP_STATE_REPLICATION_FILE, num_gener );
  }
  _robustness_file = fopen(robustness_file_name,"w");
  _replication_file = fopen(replication_file_name,"w");  
  
  assert(_replication_file!=NULL);
  assert(_robustness_file!=NULL);
  delete [] robustness_file_name;
  delete [] replication_file_name;
  
  
  // -------------------------------------
  //            Write header
  // -------------------------------------
  fprintf(_robustness_file, "# ######################################################################\n" );
  fprintf(_robustness_file, "# Robustness data of individuals at generation %"PRId32"\n",num_gener );
  fprintf(_robustness_file, "# ######################################################################\n" );
  fprintf(_robustness_file,"#  1.  Rank\n");
  fprintf(_robustness_file,"#  2.  Index\n");
  fprintf(_robustness_file,"#  3.  Fitness\n");
  fprintf(_robustness_file,"#  4.  Metabolic error\n");
  fprintf(_robustness_file,"#  5.  Genome size\n");
  fprintf(_robustness_file,"#  6.  Functional gene number\n");
  fprintf(_robustness_file,"#  7.  Reproduction probability\n");
  fprintf(_robustness_file,"#  8.  Proportion of neutral offsprings\n");
  fprintf(_robustness_file,"#  9. Proportion of beneficial offsprings\n");
  fprintf(_robustness_file,"#  10. Proportion of deleterious offsprings\n");
  fprintf(_robustness_file,"#  11. Theoretical proportion of netural offsprings\n");
  fprintf(_robustness_file,"#  12. Fitness mean of offsprings\n");
  fprintf(_robustness_file,"#  13. Fitness variance of offsprings\n");
  fprintf(_robustness_file,"#  14. Genome size mean of offsprings\n");
  fprintf(_robustness_file,"#  15. Genome size variance of offsprings\n");
  fprintf(_robustness_file,"#  16. Functional gene number mean of offsprings\n");
  fprintf(_robustness_file,"#  17. Functional gene number variance of offsprings\n");
  fprintf(_robustness_file, "# ######################################################################\n" );
  
  fprintf(_replication_file, "# #######################################################################################################\n" );
  fprintf(_replication_file,"#  Offspring details of individual with rank %"PRId32" and index %"PRId32" at generation %"PRId32" \n",_wanted_rank, _wanted_index, num_gener );
  fprintf(_replication_file, "# #######################################################################################################\n" );
  fprintf(_replication_file,"#  1.  Fitness\n");
  fprintf(_replication_file,"#  2.  Metabolic error\n");
  fprintf(_replication_file,"#  3.  Genome size\n");
  fprintf(_replication_file,"#  4.  Functional gene number\n");
  fprintf(_replication_file,"#  5.  Number of coding bases\n");
  fprintf(_replication_file,"#  6.  Number of transcribed but not translated bases\n");
  fprintf(_replication_file,"#  7.  Number of non transcribed bases\n");
  fprintf(_replication_file, "# #######################################################################################################\n" );
  
  // --------------------------------------------------------------
  //       Get genome(s) of interest and compute Fv
  // --------------------------------------------------------------
  
  _pop_size = exp_manager->get_nb_indivs();
  
  ae_individual* initial_indiv = NULL;
  ae_list_node<ae_individual*>* node = exp_manager->get_pop()->get_indivs()->get_last();
  int32_t current_rank = 1;
  int32_t current_index = -1;
  
  double* reproduction_statistics =  new double[3];
  double* offsprings_statistics =  new double[6];
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
    
    // ------------------------------------
    //         Get initial individual
    // ------------------------------------
    #ifdef __NO_X
      #ifndef __REGUL
        initial_indiv = new ae_individual( tmpind, current_index, exp_manager->get_exp_s()->get_sel()->get_prng(), exp_manager->get_exp_s()->get_sel()->get_prng());
      #else
        initial_indiv = new ae_individual_R( (dynamic_cast<ae_individual_R*>(tmpind)), current_index, exp_manager->get_exp_s()->get_sel()->get_prng(), exp_manager->get_exp_s()->get_sel()->get_prng());
          #endif
    #elif defined __X11
      #ifndef __REGUL
        initial_indiv = new ae_individual_X11( (dynamic_cast<ae_individual_X11*>(tmpind)), current_index, exp_manager->get_exp_s()->get_sel()->get_prng(), exp_manager->get_exp_s()->get_sel()->get_prng() );
      #else
        initial_indiv = new ae_individual_R_X11( (dynamic_cast<ae_individual_R_X11*>(tmpind)), current_index, exp_manager->get_exp_s()->get_sel()->get_prng(), exp_manager->get_exp_s()->get_sel()->get_prng() );
      #endif
    #endif
    
    current_index = tmpind->get_id();
    current_rank = _pop_size - tmpind->get_rank();
    
    initial_indiv->evaluate(exp_manager->get_env());
    _fitness[current_rank] = initial_indiv->get_fitness();
    _reprod_proba[current_rank] = tmp_reprod[_pop_size-current_rank];
    
    // ------------------------------------
    //              Compute Fv
    // ------------------------------------	
    th_fv = initial_indiv->compute_theoritical_f_nu();
    
    if ( (_wanted_rank+1 == current_rank) || (_wanted_index+1 == current_index))
    {
      initial_indiv->compute_experimental_f_nu( _nb_children, reproduction_statistics, offsprings_statistics, _replication_file);
    }
    else
    {
      initial_indiv->compute_experimental_f_nu( _nb_children, reproduction_statistics, offsprings_statistics);
    }
    _f_nu_pop[current_rank] = reproduction_statistics[0];
    
    // ------------------------------------
    //            Write to file
    // ------------------------------------ 
    fprintf( _robustness_file, "%"PRId32" %"PRId32" %le %le %"PRId32" %"PRId32" %le %le %le %le %le %le %le %le %le %le %le\n",
                              current_rank, current_index, initial_indiv->get_fitness(),initial_indiv->get_dist_to_target_by_feature(METABOLISM ),initial_indiv->get_total_genome_size(),
                              initial_indiv->get_nb_functional_genes(), _reprod_proba[current_rank-1], reproduction_statistics[0], reproduction_statistics[1], reproduction_statistics[2], th_fv,
                              offsprings_statistics[0], offsprings_statistics[1], offsprings_statistics[2], offsprings_statistics[3],offsprings_statistics[4],offsprings_statistics[5]);
    delete initial_indiv;
    
  }
  //delete exp_manager;
  delete [] reproduction_statistics;
  delete [] offsprings_statistics;
  
  fclose(_replication_file); 
  fclose(_robustness_file);
}

/*!
  \brief Compute statistics at population level at generation num_gener

  Compute at generation num_gener for each individual:
  * Variation
  * Population_variability
  * Proportion of neutral offspring of the best individual
  * Variability of the best individual
  * Size of the first quasi species
  * Number of different fitness
  Write this information in _output_file
  
  Before this function, compute_reproduction_stats have to be called
  \see compute_reproduction_stats(ae_exp_manager* exp_manager, int32_t num_gener)
  
  \param num_gener    current generation number  
*/
void population_statistics::compute_population_stats(int32_t num_gener)
{
  printf("\n\nCompute population stats\n");
  double variation = 0;
  double population_variability = 0; 
  
  double best_variability = 1. - _f_nu_pop[0];
  
  int fitness_number = 1;
  int quasi_species_size = 1;
  
  population_variability = _f_nu_pop[0] * _reprod_proba[0];
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
    population_variability += _f_nu_pop[i] * _reprod_proba[i];
  }
  
  population_variability = 1. - population_variability;
  variation = (static_cast<double>(fitness_number)) / (static_cast<double>(_pop_size));
  
  printf("\n Global population stats : \n \tgeneration : %"PRId32"\n \tpopulation size : %"PRId32"\n \tvariation : %le\n \tpopulation variability : %le\n \tbest f nu : %le\n \tbest variability : %le\n \tbest quasi species size : %d\n \tfitness number : %d\n", num_gener, 
	   _pop_size, variation, population_variability, _f_nu_pop[0], best_variability, 
	   quasi_species_size, fitness_number);
  
  fprintf( _output_file, "%"PRId32" %"PRId32" %le %le %le %le %d %d\n", num_gener, _pop_size, variation, population_variability, _f_nu_pop[0], best_variability, quasi_species_size, fitness_number);
  fflush( _output_file );
}

// =================================================================
//                            Protected Methods
// =================================================================

// count how many proteins were modified after replication
/*int population_statistics::count_affected_genes( ae_individual* parent, ae_individual* child )
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
}*/


// print stats about replications of the individual
/*void population_statistics::print_replication_stats( ae_individual* initial_indiv, double* fitnesses, int* nb_aff_genes )
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
}*/


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
