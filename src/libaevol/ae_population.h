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


#ifndef __AE_POPULATION_H__
#define __AE_POPULATION_H__


// =================================================================
//                              Libraries
// =================================================================
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>

#include <list>

// =================================================================
//                            Project Files
// =================================================================
#include "ae_list.h"
#include "ae_individual.h"
#include "ae_jumping_mt.h"
#include "ae_tree.h"

namespace aevol {


// =================================================================
//                          Class declarations
// =================================================================
class ae_exp_manager;

class ae_population
{
  friend class ae_selection;

  public :
    // =================================================================
    //                             Constructors
    // =================================================================
    ae_population( ae_exp_manager* exp_m );
    //~ ae_population( char* organism_file_name );
    //~ ae_population( gzFile backup_file );

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_population( void );

    // =================================================================
    //                        Accessors: Getters
    // =================================================================
    ae_list<ae_individual*>* get_indivs( void ) const;
    std::list<ae_individual*> get_indivs_std() const;
    int32_t                  get_nb_indivs( void ) const;
    ae_individual*           get_best( void ) const;
    ae_individual*                  get_indiv_by_id( int32_t id ) const;
    ae_individual*           get_indiv_by_rank( int32_t rank ) const;

    // PRNGs
    ae_jumping_mt* get_mut_prng( void ) const;
    ae_jumping_mt* get_stoch_prng( void ) const;

    // Spatial structure
    //~ double** get_secretion_present( void ) const;
    //~ double** get_secreted_amount( void ) const;
    //~ double** get_fitness_metabolic( void ) const;
    //~ double** get_fitness_total( void ) const;

    // =================================================================
    //                        Accessors: Setters
    // =================================================================
    void set_nb_indivs( size_t nb_indivs );

    void add_indiv( ae_individual* indiv );

    // PRNGs
    void set_mut_prng( ae_jumping_mt* prng );
    void set_stoch_prng( ae_jumping_mt* prng );

    // Mutation rates etc...
    void set_overall_point_mutation_rate( double point_mutation_rate);
    void set_overall_small_insertion_rate( double small_insertion_rate);
    void set_overall_small_deletion_rate( double small_deletion_rate);
    void set_overall_max_indel_size( int16_t max_indel_size);

    void set_overall_duplication_rate( double duplication_rate);
    void set_overall_deletion_rate( double deletion_rate);
    void set_overall_translocation_rate( double translocation_rate);
    void set_overall_inversion_rate( double inversion_rate);
    void set_overall_neighbourhood_rate( double neighbourhood_rate);
    void set_overall_duplication_proportion( double duplication_proportion);
    void set_overall_deletion_proportion( double deletion_proportion);
    void set_overall_translocation_proportion( double translocation_proportion);
    void set_overall_inversion_proportion( double inversion_proportion);
    void set_overall_transfer_ins_rate (double transfer_ins_rate);
    void set_overall_transfer_repl_rate (double transfer_repl_rate);

    void set_replication_reports( ae_tree* tree, int32_t num_gener);

    // =================================================================
    //                            Public Methods
    // =================================================================
    void replace_population(std::list<ae_individual*>&& new_indivs);
    //~ void            step_to_next_generation( void );
    //~ void            step_to_next_generation_grid( void );
    //~ ae_individual*  do_replication( ae_individual* parent, int32_t id, int16_t x = -1, int16_t y = -1 );
    //~ void            secretion_grid_update ( void );
    //~ ae_individual*  calculate_local_competition ( int16_t x, int16_t y );
    ae_individual*  calculate_GU_transfer ( int16_t x, int16_t y );
    void            do_random_migrations ( void );
    void     evaluate_individuals( Environment* envir );
    void            sort_individuals( void );
    void            update_best( void );

    void save( gzFile backup_file ) const;
    void load( gzFile backup_file, bool verbose );
    void load(const char* backup_file_name, bool verbose);

    #ifndef DISTRIBUTED_PRNG
      void backup_stoch_prng( void );
    #endif

    // =================================================================
    //                           Public Attributes
    // =================================================================





  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    ae_population( void )
    {
      printf( "%s:%d: error: call to forbidden constructor.\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
    ae_population( const ae_population &model )
    {
      printf( "%s:%d: error: call to forbidden constructor.\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };


    // =================================================================
    //                           Protected Methods
    // =================================================================
    void update_population(std::list<ae_individual*>&& new_indivs);
    ae_individual* create_random_individual( int32_t id );
    ae_individual* create_random_individual_with_good_gene( int32_t id );
    ae_individual* create_clone( ae_individual* dolly, int32_t id );
    ae_individual* create_individual_from_file( char* organism_file_name, int32_t id );

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    ae_exp_manager* _exp_m;

    #ifndef DISTRIBUTED_PRNG
      ae_jumping_mt* _mut_prng;
      ae_jumping_mt* _stoch_prng;
      ae_jumping_mt* _stoch_prng_bak;
    #endif

    // Individuals
    std::list<ae_individual*>  _indivs;
};

} // namespace aevol
#endif // __AE_POPULATION_H__
