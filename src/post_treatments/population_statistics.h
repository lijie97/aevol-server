


/** \class
 *  \brief
 */


#ifndef __POPULATION_STATISTICS_H__
#define __POPULATION_STATISTICS_H__


// =================================================================
//                              Libraries
// =================================================================
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include <inttypes.h>
#include <getopt.h>
#include <math.h>

// =================================================================
//                            Project Files
// =================================================================
#include <ae_list.h>
//#include <ae_common.h>
#include <ae_exp_manager.h>
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



// =================================================================
//                          Class declarations
// =================================================================


enum analysis_type
{
  ONE_GENERATION = 0,
  MULTIPLE_GENERATIONS = 1
};



class population_statistics 
{
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    population_statistics( void );
    population_statistics( analysis_type type, int32_t nb_children = 1000, int32_t wanted_rank = 1, int32_t wanted_index = -1);

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~population_statistics( void );

    // =================================================================
    //                              Accessors
    // =================================================================


    // =================================================================
    //                            Public Methods
    // =================================================================
    
    void compute_reproduction_stats(ae_exp_manager* exp_manager, int32_t num_gener);
    void compute_population_stats(int32_t num_gener);


    // =================================================================
    //                           Public Attributes
    // =================================================================



  protected :
    


    // =================================================================
    //                         Forbidden Constructors
    // =================================================================



    // =================================================================
    //                           Protected Methods
    // =================================================================
    
    //int count_affected_genes( ae_individual* parent, ae_individual* child );
    //void print_replication_stats( ae_individual* initial_indiv, double* fitnesses, int* nb_aff_genes );
    
    
    // =================================================================
    //                          Protected Attributes
    // =================================================================
    
    FILE* _output_file;
    FILE*	_robustness_file;
    FILE*	_replication_file;
    
    int 	_nb_children;
    int 	_wanted_rank;
    int 	_wanted_index;
    int32_t	_pop_size;

    double*	_f_nu_pop;
    double*	_reprod_proba;
    double*	_fitness;
    analysis_type	_type;	
};

#endif // __POPULATION_STATISTICS_H__
