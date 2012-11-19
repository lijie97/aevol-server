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


#ifndef __AE_POPULATION_H__
#define  __AE_POPULATION_H__


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



// =================================================================
//                          Class declarations
// =================================================================
class ae_exp_manager;






class ae_population : public ae_object
{
  public :
    // =================================================================
    //                             Constructors
    // =================================================================
    ae_population( ae_exp_manager* exp_m );
    //~ ae_population( char* organism_file_name );
    //~ ae_population( gzFile* backup_file );

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_population( void );

    // =================================================================
    //                        Accessors: Getters
    // =================================================================
    inline ae_list*         get_indivs( void ) const;
    inline int32_t          get_nb_indivs( void ) const;
    inline ae_individual*   get_best( void ) const;
    ae_individual *         get_indiv_by_id( int32_t id ) const;
    inline ae_individual *  get_indiv_by_rank( int32_t rank ) const;
  
    //~ inline double*          get_prob_reprod( void ) const;

    // Spatial structure
    //~ inline double** get_secretion_present( void ) const;
    //~ inline double** get_secreted_amount( void ) const;
    //~ inline double** get_fitness_metabolic( void ) const;
    //~ inline double** get_fitness_total( void ) const;

    // =================================================================
    //                        Accessors: Setters
    // =================================================================
    //~ inline void set_nb_indivs( int32_t nb_indivs );
    
    inline void add_indiv( ae_individual* indiv );
    
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

    // =================================================================
    //                            Public Methods
    // =================================================================
    void replace_population( ae_list* new_pop );
    //~ void            step_to_next_generation( void );
    //~ void            step_to_next_generation_grid( void );
    //~ ae_individual*  do_replication( ae_individual* parent, int32_t id, int16_t x = -1, int16_t y = -1 );
    //~ void            secretion_grid_update ( void ); 
    //~ ae_individual*  calculate_local_competition ( int16_t x, int16_t y );
    ae_individual*  calculate_GU_transfer ( int16_t x, int16_t y );
    void            do_random_migrations ( void );
    inline void     evaluate_individuals( ae_environment* envir );
    void            sort_individuals( void );
    void            update_best( void );

    void save( gzFile* backup_file ) const;
    void load( gzFile* backup_file, bool verbose );

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
  
    ae_individual* create_random_individual( int32_t id );
    ae_individual* create_random_individual_with_good_gene( int32_t id );
    ae_individual* create_clone( ae_individual* dolly, int32_t id );
    ae_individual* create_individual_from_file( char* organism_file_name, int32_t id );
    
    //~ void compute_prob_reprod( void );
    //~ void compute_local_prob_reprod( void );

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    ae_exp_manager* _exp_m;
    
    // Individuals
    int32_t           _nb_indivs;
    ae_list*          _indivs;
};


// =====================================================================
//                           Getters' definitions
// =====================================================================
inline ae_list* ae_population::get_indivs( void ) const
{
  return _indivs;
}

inline int32_t ae_population::get_nb_indivs( void ) const
{
  return _nb_indivs;
}

inline ae_individual* ae_population::get_best( void ) const
{
  return (ae_individual*)_indivs->get_last()->get_obj();
}

/*!
  Get the indiv corresponding to the given rank (1 for the worst indiv, POP_SIZE for the best)

  Warning, be sure you call sort_individuals() before using get_indiv_by_rank
*/
inline ae_individual * ae_population::get_indiv_by_rank( int32_t rank ) const
{
  ae_list_node* indiv_node = _indivs->get_first();
  
  for ( int32_t i = 1 ; i < rank ; i++ )
  {
    indiv_node = indiv_node->get_next();
  }
  
  // <DEBUG>
    assert( ((ae_individual*) indiv_node->get_obj())->get_rank() == rank );
  // </DEBUG>
  
  return (ae_individual*) indiv_node->get_obj();
}

/*inline double* ae_population::get_prob_reprod( void ) const
{
  if ( _prob_reprod == NULL )
  {
    printf( "ERROR, _prob_reprod has not been computed %s:%d\n", __FILE__, __LINE__ );
    exit( EXIT_FAILURE );
  }
  return _prob_reprod;
}*/


// =====================================================================
//                           Setters' definitions
// =====================================================================
//~ inline void ae_population::set_nb_indivs( int32_t nb_indivs )
//~ {
  //~ _nb_indivs = nb_indivs;
//~ }

inline void ae_population::add_indiv( ae_individual* indiv )
{
  _indivs->add( indiv );
  _nb_indivs++;
}


// Mutation rates etc...
inline void ae_population::set_overall_point_mutation_rate( double point_mutation_rate)
{
  ae_list_node*   indiv_node = _indivs->get_first();
  ae_individual*  indiv;
  for ( int32_t i = 0 ; i < _nb_indivs ; i++ )
  {
    indiv = ( ae_individual* ) indiv_node->get_obj();
    indiv->set_point_mutation_rate( point_mutation_rate );
    indiv_node = indiv_node->get_next();
  }
}

inline void ae_population::set_overall_small_insertion_rate( double small_insertion_rate)
{
  ae_list_node*   indiv_node = _indivs->get_first();
  ae_individual*  indiv;
  for ( int32_t i = 0 ; i < _nb_indivs ; i++ )
  {
    indiv = ( ae_individual* ) indiv_node->get_obj();
    indiv->set_small_insertion_rate( small_insertion_rate );
    indiv_node = indiv_node->get_next();
  }
}

inline void ae_population::set_overall_small_deletion_rate( double small_deletion_rate)
{
  ae_list_node*   indiv_node = _indivs->get_first();
  ae_individual*  indiv;
  for ( int32_t i = 0 ; i < _nb_indivs ; i++ )
  {
    indiv = ( ae_individual* ) indiv_node->get_obj();
    indiv->set_small_deletion_rate( small_deletion_rate );
    indiv_node = indiv_node->get_next();
  }
}

inline void ae_population::set_overall_max_indel_size( int16_t max_indel_size)
{
  ae_list_node*   indiv_node = _indivs->get_first();
  ae_individual*  indiv;
  for ( int32_t i = 0 ; i < _nb_indivs ; i++ )
  {
    indiv = ( ae_individual* ) indiv_node->get_obj();
    indiv->set_max_indel_size( max_indel_size );
    indiv_node = indiv_node->get_next();
  }
}

inline void ae_population::set_overall_duplication_rate( double duplication_rate)
{
  ae_list_node*   indiv_node = _indivs->get_first();
  ae_individual*  indiv;
  for ( int32_t i = 0 ; i < _nb_indivs ; i++ )
  {
    indiv = ( ae_individual* ) indiv_node->get_obj();
    indiv->set_duplication_rate( duplication_rate );
    indiv_node = indiv_node->get_next();
  }
}

inline void ae_population::set_overall_deletion_rate( double deletion_rate)
{
  ae_list_node*   indiv_node = _indivs->get_first();
  ae_individual*  indiv;
  for ( int32_t i = 0 ; i < _nb_indivs ; i++ )
  {
    indiv = ( ae_individual* ) indiv_node->get_obj();
    indiv->set_deletion_rate( deletion_rate );
    indiv_node = indiv_node->get_next();
  }
}

inline void ae_population::set_overall_translocation_rate( double translocation_rate)
{
  ae_list_node*   indiv_node = _indivs->get_first();
  ae_individual*  indiv;
  for ( int32_t i = 0 ; i < _nb_indivs ; i++ )
  {
    indiv = ( ae_individual* ) indiv_node->get_obj();
    indiv->set_translocation_rate( translocation_rate );
    indiv_node = indiv_node->get_next();
  }
}

inline void ae_population::set_overall_inversion_rate( double inversion_rate)
{
  ae_list_node*   indiv_node = _indivs->get_first();
  ae_individual*  indiv;
  for ( int32_t i = 0 ; i < _nb_indivs ; i++ )
  {
    indiv = ( ae_individual* ) indiv_node->get_obj();
    indiv->set_inversion_rate( inversion_rate );
    indiv_node = indiv_node->get_next();
  }
}

inline void ae_population::set_overall_neighbourhood_rate( double neighbourhood_rate)
{
  ae_list_node*   indiv_node = _indivs->get_first();
  ae_individual*  indiv;
  for ( int32_t i = 0 ; i < _nb_indivs ; i++ )
  {
    indiv = ( ae_individual* ) indiv_node->get_obj();
    indiv->set_neighbourhood_rate( neighbourhood_rate );
    indiv_node = indiv_node->get_next();
  }
}

inline void ae_population::set_overall_duplication_proportion( double duplication_proportion)
{
  ae_list_node*   indiv_node = _indivs->get_first();
  ae_individual*  indiv;
  for ( int32_t i = 0 ; i < _nb_indivs ; i++ )
  {
    indiv = ( ae_individual* ) indiv_node->get_obj();
    indiv->set_duplication_proportion( duplication_proportion );
    indiv_node = indiv_node->get_next();
  }
}

inline void ae_population::set_overall_deletion_proportion( double deletion_proportion)
{
  ae_list_node*   indiv_node = _indivs->get_first();
  ae_individual*  indiv;
  for ( int32_t i = 0 ; i < _nb_indivs ; i++ )
  {
    indiv = ( ae_individual* ) indiv_node->get_obj();
    indiv->set_deletion_proportion( deletion_proportion );
    indiv_node = indiv_node->get_next();
  }
}

inline void ae_population::set_overall_translocation_proportion( double translocation_proportion)
{
  ae_list_node*   indiv_node = _indivs->get_first();
  ae_individual*  indiv;
  for ( int32_t i = 0 ; i < _nb_indivs ; i++ )
  {
    indiv = ( ae_individual* ) indiv_node->get_obj();
    indiv->set_translocation_proportion( translocation_proportion );
    indiv_node = indiv_node->get_next();
  }
}

inline void ae_population::set_overall_inversion_proportion( double inversion_proportion)
{
  ae_list_node*   indiv_node = _indivs->get_first();
  ae_individual*  indiv;
  for ( int32_t i = 0 ; i < _nb_indivs ; i++ )
  {
    indiv = ( ae_individual* ) indiv_node->get_obj();
    indiv->set_inversion_proportion( inversion_proportion );
    indiv_node = indiv_node->get_next();
  }
}

// =====================================================================
//                       Inline functions' definition
// =====================================================================
inline void ae_population::evaluate_individuals( ae_environment* envir )
{
  ae_list_node *  indiv_node  = _indivs->get_first();
  ae_individual * indiv       = NULL;

  while ( indiv_node != NULL )
  {
    indiv = (ae_individual *) indiv_node->get_obj();
    indiv->evaluate( envir );
    indiv->compute_statistical_data();
    
    indiv_node = indiv_node->get_next();
  }
}

#endif // __AE_POPULATION_H__
