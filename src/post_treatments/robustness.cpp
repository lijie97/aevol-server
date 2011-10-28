//
// This program takes individuals and computes f_nu
//
// Option are :
//     * -f / --file   backup_file
//            backup file where individuals should be extracted
//     * -o / --output   output_directory
//            output directory where statistics should be printed (default ".")
//     * -n / --nb-children  nb
//            perform nb replications (default 1000)
//     * -r / --rank  i
//            get individual with ith rank (default -1 -> "all individuals")
//            where 1st is the best individual !!!
//     * -i / --index  i
//            get individual with index i (default -1 -> "all individuals")
//     * -d / --with-details
//            print additional information in output_file/replications.out (see below)
//
// Data concerning individuals wanted is printed in one or two files.
//
// The first file (output_dir/fv.out) contains information about Fv:
//
//       "r_i_fit_pr_fvexp_fvexpnob_fvth_N\n"
//
// where:
//     * _ is a blank space
//     * r:        Rank of individual in population
//     * i:        Index of individual in population
//     * fit:      FITness of the individual
//     * pr:       Probability of Reproduction
//     * fvexp:    EXPerimental estimation of Fv
//     * fvexpnob: EXPerimental estimation of Fv counting children Neutral Or Better
//     * fvth:     THeoritical estimation of Fv
//     * N:        Number of individuals in population
//      
// The second file (output_dir/replications.out), if activated (-d option),
// contains information about children obtained when experimentally estimating Fv:
//
//       "#i_gl_nfg_cl_f_afc_sdfc_psp_pbp_pnga_anga\n"
//
// where:
//     * # is # and _ is a blank space
//     * i:     Index of individual in population
//     * gl:    total Genome Length
//     * nfg:   Number of Functional Genes
//     * cl:    total Coding Length
//     * f:     Fitness of individual
//     * afc:   Average Fitness of the Children
//     * sdfc:  Standard Deviation of Fitness values of the Children
//     * psp:   Proportion of children with Same fitness as Parent
//     * pbp:   Proportion of children Better than Parent
//     * pnga:  Proportion of children with No Gene Affected
//     * anga:  Average Number of Gene Affected by replication
//
// One additional line per child is printed with its fitness and the number of affected 
// genes ("fc_nga\n").
//
//
// Examples :
//
// For generation 20000, compute statistics for all the
// individuals and print them in directory out_020000 :
//
//    robustness -f backup/gen_020000.ae -o out_020000
//
// For generation 20000, write the best individual's statistics in
// out_020000_best with details about replications :
//
//    robustness -d -r 1 -f backup/gen_020000.ae -o out_020000_best
//



/** \class 
 *  \brief 
 *  \author  Carole Knibbe
 *  Input files :
 *               
 *  Output file :
 *               
 */


// =======================================================================
//                        Standard Libraries
// =======================================================================

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include <inttypes.h>
#include <getopt.h>
#include <math.h>


// =======================================================================
//                        Project Libraries
// =======================================================================

#include <ae_list.h>
#include <ae_common.h>
#include <ae_simulation.h>
#include <ae_population.h>
#include <ae_protein.h>

#ifdef __NO_X
#ifdef __REGUL
#include <ae_individual_R.h>
#else
#include <ae_individual.h>
#endif
#elif defined __X11
#ifdef __REGUL
#include <ae_individual_R_X11.h>
#else
#include <ae_individual_X11.h>
#endif
#endif



// =======================================================================
//                       Secondary Functions
// =======================================================================


// function copied from ae_individual's computation of experimental fv.
// In addition, it provides and prints information about replications
double compute_experimental_fv( ae_individual* indiv, int nb_children, double* neutral_or_better, FILE* replication_file );

// count how many proteins were modified after replication
int count_affected_genes( ae_individual* parent, ae_individual* child );

// print stats about replications of the individual
void print_replication_stats( FILE* output_file, ae_individual* initial_indiv, int nb_children, double* fitnesses, int* nb_aff_genes );

// TODO: update this function...
// reconstruct final individual from backup and lineage
// ae_individual * get_final_individual_using_dstory();     



// =====================================================================
//                         Main Function
// =====================================================================


#define FV_FILE "fv.out"
#define REP_FILE "replications.out"

int main( int argc, char* argv[] ) 
{
  // ----------------------------------------
  //     command-line option parsing
  // ----------------------------------------
  char* backup_file_name = NULL;
  char* output_dir = NULL;
  int nb_children = 1000;
  int wanted_rank = -1;
  int wanted_index = -1;
  bool details = false;

  const char * options_list = "f:o:n:r:i:d"; 
  static struct option long_options_list[] = {
    {"file", 1, NULL, 'f'},
    {"output", 1, NULL, 'o'},
    {"nb-children", 1, NULL, 'n'},
    {"rank", 1, NULL, 'r'},
    {"index", 1, NULL, 'i'},
    {"with-details", 0, NULL, 'd'},
    {0, 0, 0, 0}
  };

  int option = -1;
  while((option=getopt_long(argc,argv,options_list,long_options_list,NULL))!=-1) 
  {
    switch(option) 
    {
    case 'f' : 
      backup_file_name = new char[strlen(optarg) + 1];
      sprintf( backup_file_name, "%s", optarg);
      break;
    case 'o' : 
      output_dir = new char[strlen(optarg) + 1];
      sprintf( output_dir, "%s", optarg );
      break;
    case 'n' :
      nb_children = atoi(optarg);
      break;  
    case 'r' :
      wanted_rank = atoi(optarg);
      wanted_index = -1;
      break;  
    case 'i' :
      wanted_index = atoi(optarg);
      wanted_rank = -1;
      break;  
    case 'd' :
      details=true;
      break;  
    }
  }

  // ----------------------------------------
  //              Initializing....
  // ----------------------------------------
  FILE * fv_output = NULL;
  FILE * replication_output = NULL;
  if (output_dir==NULL) // default = .
  {
    fv_output = fopen(FV_FILE,"w");
    if (fv_output==NULL) { fprintf(stderr, "Warning: could not open file %s.\n", FV_FILE); }
    if ( details == true )
    { 
      replication_output=fopen(REP_FILE,"w");
      if (replication_output==NULL) { fprintf(stderr, "Warning: could not open file %s.\n", REP_FILE); }
    }
  }
  else
  {
    char* fv_file_name = new char [ strlen(output_dir) + strlen(FV_FILE) + 2 ];
    char* replication_file_name = new char [ strlen(output_dir) + strlen(REP_FILE) + 2 ];
    sprintf( fv_file_name, "%s/%s", output_dir, FV_FILE );
    sprintf( replication_file_name, "%s/%s", output_dir, REP_FILE );

    fv_output = fopen(fv_file_name,"w");
    if ( fv_output == NULL ) { fprintf(stderr, "Warning: could not open file %s.\n", fv_file_name); }
    if ( details == true )
    { 
      replication_output=fopen(replication_file_name,"w");
      if (replication_output==NULL) { fprintf(stderr, "Warning: could not open file %s.\n", replication_file_name); }
    }

    delete [] output_dir;
    delete [] fv_file_name;
    delete [] replication_file_name;
  }
  fflush( stderr );

  // Load the simulation from backup
  if ( backup_file_name != NULL )
  {
    // the second argument (false) indicates that we just want to do
    // post-processing, no need to initialize structures like 'tree'
    // that are modified at each generation
    ae_common::sim = new ae_simulation( backup_file_name, false, NULL );
  }
  else
  {
    fprintf(stderr, "You have to specify an existing backup file...\n");
    exit(EXIT_FAILURE);
  }
  delete [] backup_file_name;
   
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
  
  ae_common::sim->get_pop()->evaluate_individuals(ae_common::sim->get_env());
 
  // --------------------------------------------------------------
  //       Get genome(s) of interest and compute Fv
  // --------------------------------------------------------------
  
  ae_individual* initial_indiv = NULL;
  ae_list_node* node = ae_common::sim->get_pop()->get_indivs()->get_last();
  int current_rank = 1;
  int current_index = -1;
  
  // simply parse through individuals and keep those that are wanted
  while (node != NULL)
  {
    ae_individual* tmpind = (ae_individual*) node->get_obj();
    node = node->get_prev();

    if ( (wanted_rank==-1) || (wanted_rank == current_rank) )
    {
      current_index = tmpind->get_index_in_population();
      if ( (wanted_index == -1) || (wanted_index == current_index) )
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
	initial_indiv->evaluate(ae_common::sim->get_env());
	printf("found indiv %"PRId32" with fitness %lf (rank %d)\n", tmpind->get_index_in_population(), tmpind->get_fitness(), current_rank);

	// ------------------------------------
	//              Compute Fv
	// ------------------------------------
	double neutral_or_better = -1;
	double exp_fv = -1;
	if ( details == true )
	{ 
	  // 'overloaded' function
	  exp_fv = compute_experimental_fv( initial_indiv, nb_children, &neutral_or_better, replication_output );
	}
	else
	{
	  // aevol function
	  exp_fv = initial_indiv->compute_experimental_f_nu( nb_children, &neutral_or_better );
	}
	double th_fv = initial_indiv->compute_theoritical_f_nu();
	      
	// ------------------------------------
	//            Write to file
	// ------------------------------------     
	double* reprod_proba = ae_common::sim->get_pop()->get_prob_reprod();
	if ( fv_output != NULL )
	{
	  fprintf( fv_output, "%d %d %lf %lf %lf %lf %lf %"PRId32"\n",
		   current_rank, current_index, initial_indiv->get_fitness(),
		   reprod_proba[ae_common::init_pop_size-current_rank],
		   exp_fv, neutral_or_better, th_fv, ae_common::init_pop_size );
	  // fflush(fv_output);
	}
	      
	delete initial_indiv;
      }
    }
    current_rank++;
  }
  
  if (fv_output != NULL) { fclose(fv_output); }
  if (replication_output != NULL) { fclose(replication_output); }
  
  // TODO: update
  // rename dstory.bak.gz
  //rename("kept_dstory.bak.gz","dstory.bak.gz");
  //if(generation%param::back_step == 0) delete po
  
  ae_common::clean();

  return EXIT_SUCCESS;
}



//======================================================================
//              Implementation of secondary functions
//======================================================================


double compute_experimental_fv( ae_individual* indiv, int nb_children, double* neutral_or_better, FILE* replication_output )
{
  double initial_fitness = indiv->get_fitness();
  int32_t index = indiv->get_index_in_population();
  double Fv = 0;
  if ( neutral_or_better != NULL) *neutral_or_better = 0;
  
  // ------------------------------------------
  //      Simulate fitness degradation
  // ------------------------------------------
  double childfit = 0;
  double * fitnesses = new double [nb_children];
  int * nb_aff_genes = new int [nb_children];
	      
  // reproduce this individual to create 'nb_children' children 
  ae_individual * child = NULL;
  for (int i = 0; i < nb_children; i++)
  {
    child = ae_common::sim->get_pop()->do_replication(indiv, index); 
    childfit = child->get_fitness(); // child is automatically evaluated
    
    if (replication_output != NULL) // only interesting if we put that into a file...
    {
      fitnesses[i]    = childfit;
      nb_aff_genes[i] = count_affected_genes(indiv, child);
    }
    
    //count neutral offspring
    if ( fabs(initial_fitness - childfit) < 1e-15 )
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
  Fv /= nb_children;
  if ( neutral_or_better != NULL ) *neutral_or_better /= nb_children;
  	      
  // ------------------------------------------
  //      print statistics to file
  // ------------------------------------------
  if ( replication_output != NULL )
  {
    //general statistics of the individual
    fprintf( replication_output, "#%"PRId32" %"PRId32" %"PRId32" %d %.15lf",
	    indiv->get_index_in_population(),
	    indiv->get_total_genome_size(),
	    indiv->get_nb_functional_genes(),
	    (int) indiv->get_overall_size_coding_RNAs(),
	    indiv->get_fitness());
    //statistics of replication
    print_replication_stats( replication_output, indiv, nb_children, fitnesses, nb_aff_genes );    
  }  

  if ( fitnesses != NULL ) delete [] fitnesses;
  if ( nb_aff_genes != NULL ) delete [] nb_aff_genes;

  return Fv;
}


// ==============================================================================


int count_affected_genes( ae_individual * parent, ae_individual* child )
{  
  // ------------------------------------------------------------
  //       list all functional proteins of the child
  // ------------------------------------------------------------
  
  //copy protein list
  ae_list* child_protein_list = new ae_list(*(child->get_protein_list()));
  //delete all nodes containing degenerated proteins
  ae_list_node* current_prot_node = child_protein_list->get_first();
  ae_list_node* next_prot_node = NULL;
  ae_protein* current_prot = NULL;
  while ( current_prot_node != NULL )
  {
    current_prot = (ae_protein*) current_prot_node->get_obj();
    next_prot_node = current_prot_node->get_next();
    if ( current_prot->get_is_functional() == false ) 
    {
      //delete the node but not the protein
      child_protein_list->remove(current_prot_node, DELETE_OBJ, NO_DELETE);
    }
    current_prot_node = next_prot_node;
  }
  
  // ----------------------------------------------------------------
  //    for each functional protein of the parent, search for
  //    the same one in the child's proteins and pop it out of the  
  //    list if found
  // ----------------------------------------------------------------
  ae_list_node* parent_prot_node = parent->get_protein_list()->get_first();
  ae_protein* parent_prot = NULL;
  ae_list_node* child_prot_node = NULL;
  ae_protein* child_prot = NULL;
  bool found = false;

  int nb_proteins_left = 0;

  // parse parent protein list
  while ( parent_prot_node != NULL )
  {
    parent_prot = (ae_protein*) parent_prot_node->get_obj();
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
	  child_protein_list->remove(child_prot_node, DELETE_OBJ, NO_DELETE);
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


//==============================================================================

void print_replication_stats( FILE* output_file, ae_individual* initial_indiv, int nb_children, double* fitnesses, int* nb_aff_genes )
{
  if ( output_file == NULL ) return;
  
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

  for (int i = 0; i < nb_children; i++)
  {
    average_fit += fitnesses[i];
    avg_nb_aff_genes += nb_aff_genes[i];
    sd += fitnesses[i]*fitnesses[i];
    if ( fabs(fitnesses[i] - fitness_initial_indiv) < 1e-15 ) { prop_same_as_parent += 1; }
    else if ( fitnesses[i] > fitness_initial_indiv ) { prop_better_than_parent += 1; }
    if ( nb_aff_genes[i] == 0 ) { prop_no_genes_affected += 1; }
  }
  average_fit /= nb_children;
  avg_nb_aff_genes /= nb_children;
  prop_no_genes_affected /= nb_children;
  sd /= nb_children;
  sd -= (average_fit)*(average_fit);
  sd = sqrt(sd);
  prop_better_than_parent /= nb_children;
  prop_same_as_parent /= nb_children;

  // --------------------------------------------------------------
  //         print statistics on current line
  // --------------------------------------------------------------
  fprintf(output_file, " %.15lf %.15lf %.15lf %.15lf %.15lf %.15lf ", average_fit, sd, prop_same_as_parent, prop_better_than_parent, prop_no_genes_affected, avg_nb_aff_genes);
  fprintf(output_file, "\n");
  
  // --------------------------------------------------------------
  //         print fitness statistics for each child
  // --------------------------------------------------------------
  for (int i = 0; i < nb_children; i++)
  {
    fprintf(output_file, "%.15lf %d \n", fitnesses[i], nb_aff_genes[i]);
  }
}

//==============================================================================

/*
  ae_individual * get_final_individual_using_dstory()
  {
  // warning: delete the genome when you don't need it anymore

  gzFile * dstory_file = (gzFile *) gzopen("dstory.bak.gz","r");

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

  mut_list.erase(DELETE_OBJ);
  } 

  gzclose(dstory_file);
  
  indiv->compute_phenotype(true);
  return indiv;
  }
*/
