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



// =================================================================
//                            Project Files
// =================================================================
#include <ae_list.h>
#include <ae_individual.h>
#include <ae_jumping_mt.h>
#include <ae_tree.h>

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
    inline ae_list<ae_individual*>* get_indivs( void ) const;
    inline int32_t                  get_nb_indivs( void ) const;
    inline ae_individual*           get_best( void ) const;
    ae_individual*                  get_indiv_by_id( int32_t id ) const;
    inline ae_individual*           get_indiv_by_rank( int32_t rank ) const;

    // PRNGs
    inline ae_jumping_mt* get_mut_prng( void ) const;
    inline ae_jumping_mt* get_stoch_prng( void ) const;

    // Spatial structure
    //~ inline double** get_secretion_present( void ) const;
    //~ inline double** get_secreted_amount( void ) const;
    //~ inline double** get_fitness_metabolic( void ) const;
    //~ inline double** get_fitness_total( void ) const;

    // =================================================================
    //                        Accessors: Setters
    // =================================================================
    void set_nb_indivs( int32_t nb_indivs );

    inline void add_indiv( ae_individual* indiv );

    // PRNGs
    inline void set_mut_prng( ae_jumping_mt* prng );
    inline void set_stoch_prng( ae_jumping_mt* prng );

    // Mutation rates etc...
    inline void set_overall_point_mutation_rate( double point_mutation_rate);
    inline void set_overall_small_insertion_rate( double small_insertion_rate);
    inline void set_overall_small_deletion_rate( double small_deletion_rate);
    inline void set_overall_max_indel_size( int16_t max_indel_size);

    inline void set_overall_duplication_rate( double duplication_rate);
    inline void set_overall_deletion_rate( double deletion_rate);
    inline void set_overall_translocation_rate( double translocation_rate);
    inline void set_overall_inversion_rate( double inversion_rate);
    inline void set_overall_neighbourhood_rate( double neighbourhood_rate);
    inline void set_overall_duplication_proportion( double duplication_proportion);
    inline void set_overall_deletion_proportion( double deletion_proportion);
    inline void set_overall_translocation_proportion( double translocation_proportion);
    inline void set_overall_inversion_proportion( double inversion_proportion);
    inline void set_overall_transfer_ins_rate (double transfer_ins_rate);
    inline void set_overall_transfer_repl_rate (double transfer_repl_rate);

    inline void set_replication_reports( ae_tree* tree, int32_t num_gener);

    // =================================================================
    //                            Public Methods
    // =================================================================
    void replace_population( ae_list<ae_individual*>* new_indivs );
    //~ void            step_to_next_generation( void );
    //~ void            step_to_next_generation_grid( void );
    //~ ae_individual*  do_replication( ae_individual* parent, int32_t id, int16_t x = -1, int16_t y = -1 );
    //~ void            secretion_grid_update ( void );
    //~ ae_individual*  calculate_local_competition ( int16_t x, int16_t y );
    ae_individual*  calculate_GU_transfer ( int16_t x, int16_t y );
    void            do_random_migrations ( void );
    inline void     evaluate_individuals( Environment* envir );
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
    void update_population(ae_list<ae_individual*>* new_indivs);
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
    int32_t                   _nb_indivs;
    ae_list<ae_individual*>*  _indivs;
    // new individual list targeted at replacing _nb_indivs and _indivs
    std::list<ae_individual> std_indivs;
};


// =====================================================================
//                           Getters' definitions
// =====================================================================
inline ae_list<ae_individual*>* ae_population::get_indivs( void ) const
{
  return _indivs;
}

inline int32_t ae_population::get_nb_indivs( void ) const
{
  return _nb_indivs;
}

inline ae_individual* ae_population::get_best( void ) const
{
  return _indivs->get_last()->get_obj();
}

/*!
  Get the indiv corresponding to the given rank (1 for the worst indiv, POP_SIZE for the best)

  Warning, be sure you call sort_individuals() before using get_indiv_by_rank
*/
inline ae_individual* ae_population::get_indiv_by_rank( int32_t rank ) const
{
  ae_list_node<ae_individual*>* indiv_node = _indivs->get_first();

  for ( int32_t i = 1 ; i < rank ; i++ )
  {
    indiv_node = indiv_node->get_next();
  }

  assert( indiv_node->get_obj()->get_rank() == rank );

  return indiv_node->get_obj();
}

inline ae_jumping_mt* ae_population::get_mut_prng( void ) const
{
  return _mut_prng;
}

inline ae_jumping_mt* ae_population::get_stoch_prng( void ) const
{
  return _stoch_prng;
}

// =====================================================================
//                           Setters' definitions
// =====================================================================
/*inline void ae_population::set_nb_indivs( int32_t nb_indivs )
{
  _nb_indivs = nb_indivs;
}*/

inline void ae_population::add_indiv( ae_individual* indiv )
{
  _indivs->add( indiv );
  _nb_indivs++;
}


inline void ae_population::set_mut_prng( ae_jumping_mt* prng )
{
  if (_mut_prng != NULL) delete _mut_prng;
  _mut_prng = prng;
  ae_list_node<ae_individual*>* indiv_node = _indivs->get_first();
  ae_individual*  indiv;
  for ( int32_t i = 0 ; i < _nb_indivs ; i++ )
  {
    indiv = indiv_node->get_obj();
    indiv->set_mut_prng( _mut_prng );
    indiv_node = indiv_node->get_next();
  }
}

inline void ae_population::set_stoch_prng( ae_jumping_mt* prng )
{
  if (_stoch_prng != NULL) delete _stoch_prng;
  _stoch_prng = prng;
  ae_list_node<ae_individual*>* indiv_node = _indivs->get_first();
  ae_individual*  indiv;
  for ( int32_t i = 0 ; i < _nb_indivs ; i++ )
  {
    indiv = indiv_node->get_obj();
    indiv->set_stoch_prng( _stoch_prng );
    indiv_node = indiv_node->get_next();
  }
}


// Mutation rates etc...
inline void ae_population::set_overall_point_mutation_rate( double point_mutation_rate )
{
  ae_list_node<ae_individual*>* indiv_node = _indivs->get_first();
  ae_individual*  indiv;
  for ( int32_t i = 0 ; i < _nb_indivs ; i++ )
  {
    indiv = indiv_node->get_obj();
    indiv->set_point_mutation_rate( point_mutation_rate );
    indiv_node = indiv_node->get_next();
  }
}

inline void ae_population::set_overall_small_insertion_rate( double small_insertion_rate )
{
  ae_list_node<ae_individual*>* indiv_node = _indivs->get_first();
  ae_individual*  indiv;
  for ( int32_t i = 0 ; i < _nb_indivs ; i++ )
  {
    indiv = indiv_node->get_obj();
    indiv->set_small_insertion_rate( small_insertion_rate );
    indiv_node = indiv_node->get_next();
  }
}

inline void ae_population::set_overall_small_deletion_rate( double small_deletion_rate )
{
  ae_list_node<ae_individual*>* indiv_node = _indivs->get_first();
  ae_individual*  indiv;
  for ( int32_t i = 0 ; i < _nb_indivs ; i++ )
  {
    indiv = indiv_node->get_obj();
    indiv->set_small_deletion_rate( small_deletion_rate );
    indiv_node = indiv_node->get_next();
  }
}

inline void ae_population::set_overall_max_indel_size( int16_t max_indel_size )
{
  ae_list_node<ae_individual*>* indiv_node = _indivs->get_first();
  ae_individual*  indiv;
  for ( int32_t i = 0 ; i < _nb_indivs ; i++ )
  {
    indiv = indiv_node->get_obj();
    indiv->set_max_indel_size( max_indel_size );
    indiv_node = indiv_node->get_next();
  }
}

inline void ae_population::set_overall_duplication_rate( double duplication_rate )
{
  ae_list_node<ae_individual*>* indiv_node = _indivs->get_first();
  ae_individual*  indiv;
  for ( int32_t i = 0 ; i < _nb_indivs ; i++ )
  {
    indiv = indiv_node->get_obj();
    indiv->set_duplication_rate( duplication_rate );
    indiv_node = indiv_node->get_next();
  }
}

inline void ae_population::set_overall_deletion_rate( double deletion_rate)
{
  ae_list_node<ae_individual*>* indiv_node = _indivs->get_first();
  ae_individual*  indiv;
  for ( int32_t i = 0 ; i < _nb_indivs ; i++ )
  {
    indiv = indiv_node->get_obj();
    indiv->set_deletion_rate( deletion_rate );
    indiv_node = indiv_node->get_next();
  }
}

inline void ae_population::set_overall_translocation_rate( double translocation_rate)
{
  ae_list_node<ae_individual*>* indiv_node = _indivs->get_first();
  ae_individual*  indiv;
  for ( int32_t i = 0 ; i < _nb_indivs ; i++ )
  {
    indiv = indiv_node->get_obj();
    indiv->set_translocation_rate( translocation_rate );
    indiv_node = indiv_node->get_next();
  }
}

inline void ae_population::set_overall_inversion_rate( double inversion_rate)
{
  ae_list_node<ae_individual*>* indiv_node = _indivs->get_first();
  ae_individual*  indiv;
  for ( int32_t i = 0 ; i < _nb_indivs ; i++ )
  {
    indiv = indiv_node->get_obj();
    indiv->set_inversion_rate( inversion_rate );
    indiv_node = indiv_node->get_next();
  }
}

inline void ae_population::set_overall_transfer_ins_rate (double transfer_ins_rate)
{
	ae_list_node<ae_individual*>* indiv_node = _indivs->get_first();
  ae_individual*  indiv;
  for ( int32_t i = 0 ; i < _nb_indivs ; i++ )
  {
    indiv = indiv_node->get_obj();
    indiv->set_HT_ins_rate( transfer_ins_rate );
    indiv_node = indiv_node->get_next();
  }
}

inline void ae_population::set_overall_transfer_repl_rate (double transfer_repl_rate)
{
	ae_list_node<ae_individual*>* indiv_node = _indivs->get_first();
  ae_individual*  indiv;
  for ( int32_t i = 0 ; i < _nb_indivs ; i++ )
  {
    indiv = indiv_node->get_obj();
    indiv->set_HT_repl_rate( transfer_repl_rate );
    indiv_node = indiv_node->get_next();
  }
}

inline void ae_population::set_overall_neighbourhood_rate( double neighbourhood_rate)
{
  ae_list_node<ae_individual*>* indiv_node = _indivs->get_first();
  ae_individual*  indiv;
  for ( int32_t i = 0 ; i < _nb_indivs ; i++ )
  {
    indiv = indiv_node->get_obj();
    indiv->set_neighbourhood_rate( neighbourhood_rate );
    indiv_node = indiv_node->get_next();
  }
}

inline void ae_population::set_overall_duplication_proportion( double duplication_proportion)
{
  ae_list_node<ae_individual*>* indiv_node = _indivs->get_first();
  ae_individual*  indiv;
  for ( int32_t i = 0 ; i < _nb_indivs ; i++ )
  {
    indiv = indiv_node->get_obj();
    indiv->set_duplication_proportion( duplication_proportion );
    indiv_node = indiv_node->get_next();
  }
}

inline void ae_population::set_overall_deletion_proportion( double deletion_proportion)
{
  ae_list_node<ae_individual*>* indiv_node = _indivs->get_first();
  ae_individual*  indiv;
  for ( int32_t i = 0 ; i < _nb_indivs ; i++ )
  {
    indiv = indiv_node->get_obj();
    indiv->set_deletion_proportion( deletion_proportion );
    indiv_node = indiv_node->get_next();
  }
}

inline void ae_population::set_overall_translocation_proportion( double translocation_proportion)
{
  ae_list_node<ae_individual*>* indiv_node = _indivs->get_first();
  ae_individual*  indiv;
  for ( int32_t i = 0 ; i < _nb_indivs ; i++ )
  {
    indiv = indiv_node->get_obj();
    indiv->set_translocation_proportion( translocation_proportion );
    indiv_node = indiv_node->get_next();
  }
}

inline void ae_population::set_overall_inversion_proportion( double inversion_proportion)
{
  ae_list_node<ae_individual*>* indiv_node = _indivs->get_first();
  ae_individual*  indiv;
  for ( int32_t i = 0 ; i < _nb_indivs ; i++ )
  {
    indiv = indiv_node->get_obj();
    indiv->set_inversion_proportion( inversion_proportion );
    indiv_node = indiv_node->get_next();
  }
}

inline void ae_population::set_replication_reports( ae_tree* tree, int32_t num_gener)
{
  ae_list_node<ae_individual*>* indiv_node = _indivs->get_first();
  ae_individual*  indiv;
  for ( int32_t i = 0 ; i < _nb_indivs ; i++ )
  {
    indiv = indiv_node->get_obj();
    indiv->set_replication_report( tree->get_report_by_index( num_gener, indiv->get_id()));
    indiv_node = indiv_node->get_next();
  }
}

// =====================================================================
//                       Inline functions' definition
// =====================================================================
inline void ae_population::evaluate_individuals( Environment* envir )
{
  ae_list_node<ae_individual*>*  indiv_node  = _indivs->get_first();
  ae_individual * indiv       = NULL;

  while ( indiv_node != NULL )
  {
    indiv = indiv_node->get_obj();
    indiv->evaluate( envir );
    indiv->compute_statistical_data();

    indiv_node = indiv_node->get_next();
  }
}
} // namespace aevol

#endif // __AE_POPULATION_H__
