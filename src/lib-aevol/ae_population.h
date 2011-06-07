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






class ae_population : public ae_object
{
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    ae_population( void );
    ae_population( char* organism_file_name );
    ae_population( gzFile* backup_file );

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_population( void );

    // =================================================================
    //                              Accessors
    // =================================================================
    inline ae_list*         get_indivs( void )                  const;
    inline int32_t          get_nb_indivs( void )               const;
    inline ae_individual*   get_best( void )                    const;
    ae_individual *         get_indiv_by_index( int32_t index ) const;
    inline ae_individual *  get_indiv_by_rank( int32_t rank )   const;
  
    inline ae_grid_cell***  get_pop_grid( void );
    inline double**         get_secretion_present( void );
    inline double**         get_secreted_amount( void );
    inline double**         get_fitness_metabolic( void );
    inline double**         get_fitness_total( void );

    // =================================================================
    //                            Public Methods
    // =================================================================
    void            step_to_next_generation( void );
    void            step_to_next_generation_grid( void );
    ae_individual*  do_replication( ae_individual* parent, int32_t index, int16_t x = 0, int16_t y = 0 );
    void            secretion_grid_update ( void ); 
    ae_individual*  calculate_local_competition ( int16_t x, int16_t y );
    ae_individual*  calculate_GU_transfer ( int16_t x, int16_t y );
    void            do_random_migrations ( void );
    inline void     evaluate_individuals( ae_environment* envir );

    void    write_to_backup( gzFile* backup_file );

    // =================================================================
    //                           Public Attributes
    // =================================================================





  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    //~ ae_population( void )
    //~ {
      //~ printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      //~ exit( EXIT_FAILURE );
    //~ };
    ae_population( const ae_population &model )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };


    // =================================================================
    //                           Protected Methods
    // =================================================================
    void sort_individuals( void );
  
    ae_individual* create_random_individual( int32_t index );
    ae_individual* create_random_individual_with_good_gene( int32_t index );
    ae_individual* create_clone( ae_individual* dolly, int32_t index );
    ae_individual* create_individual_from_file( char* organism_file_name, int32_t index );
    
    void compute_prob_reprod( void );

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    ae_list*          _indivs;
    ae_grid_cell***   _pop_grid;
    double**          _output_grid;
    
    int32_t           _nb_indivs;
 
    double*           _prob_reprod;
    double            _prob_reprod_previous_best;
};


// =====================================================================
//                          Accessors definitions
// =====================================================================
ae_list* ae_population::get_indivs( void ) const
{
  return _indivs;
}

int32_t ae_population::get_nb_indivs( void ) const
{
  return _nb_indivs;
}

ae_individual* ae_population::get_best( void ) const
{
  return (ae_individual*)_indivs->get_last()->get_obj();
}

ae_individual * ae_population::get_indiv_by_rank( int32_t rank ) const
{
  ae_list_node* indiv_node = _indivs->get_first();
  
  for ( int32_t i = 1 ; i < rank ; i++ )
  {
    indiv_node = indiv_node->get_next();
  }
  
  // <DEBUG>
    assert( ((ae_individual*) indiv_node->get_obj())->get_rank_in_population() == rank );
  // </DEBUG>
  
  return (ae_individual*) indiv_node->get_obj();
}



ae_grid_cell*** ae_population::get_pop_grid( void )
{
  return _pop_grid;
}

double** ae_population::get_secretion_present( void )
{
  for (int16_t x = 0; x < ae_common::grid_x; x++) 
  {
    for (int16_t y = 0; y < ae_common::grid_y; y++) 
    { 
      _output_grid[x][y] = _pop_grid[x][y]->get_compound_amount(); 
    }
  }
  
  return _output_grid;
}

double** ae_population::get_secreted_amount( void )
{
  for (int16_t x = 0; x < ae_common::grid_x; x++) 
  {
    for (int16_t y = 0; y < ae_common::grid_y; y++) 
    { 
      _output_grid[x][y] = _pop_grid[x][y]->get_individual()->get_fitness_by_feature(SECRETION); 
    } 
  }
  
  return _output_grid;
}

double** ae_population::get_fitness_metabolic( void )
{
  for (int16_t x = 0; x < ae_common::grid_x; x++) 
  {
    for (int16_t y = 0; y < ae_common::grid_y; y++) 
    { 
       _output_grid[x][y] = _pop_grid[x][y]->get_individual()->get_fitness_by_feature(METABOLISM); 
    } 
  }
  
  return _output_grid;
}

double** ae_population::get_fitness_total( void )
{
  for (int16_t x = 0; x < ae_common::grid_x; x++) 
  {
    for (int16_t y = 0; y < ae_common::grid_y; y++) 
    { 
      _output_grid[x][y] = _pop_grid[x][y]->get_individual()->get_fitness(); 
    } 
  }
  
  return _output_grid;
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
