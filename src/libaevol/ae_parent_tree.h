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


/*! \class ae_parent_tree
  \brief
*/
 
 
#ifndef  __AE_PARENT_TREE_H__
#define  __AE_PARENT_TREE_H__
 
 
// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>
#include <stdio.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_object.h>
#include <ae_individual.h>
#include <ae_replication_report.h>



// =================================================================
//                          Class declarations
// =================================================================


class ae_parent_tree : public ae_object
{  
  public :
  
    // =================================================================
    //                             Constructors
    // =================================================================
    inline ae_parent_tree( void );
    // ae_dump( int32_t num_gener );
    
    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_parent_tree( void )
      {
        for( int32_t n = 0 ; n < ae_common::run_for_gens ; n++ )
        delete parent[n];
        delete parent;
      }
    
    // =================================================================
    //                              Accessors
    // =================================================================
    
    // =================================================================
    //                            Public Methods
    // =================================================================

    void insert_current_generation( void ); // insert
    inline int32_t get_LCA( int32_t n, int32_t a, int32_t b, int32_t* j );
    // return the generation of the last common ancestor of a and b (of generation n). stores its index in j if provided.
    
    // =================================================================
    //                           Public Attributes
    // =================================================================





  protected :
    
    
    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    ae_parent_tree( const ae_parent_tree &model )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };


    // =================================================================
    //                           Protected Methods
    // =================================================================
    inline void insert_indiv( int32_t n, int32_t i, int32_t j );
    // insert individual i from generation n with parent j from generation n-1
    /*
    void _dump(); // FOR DEBUG PURPOSE ONLY
    */

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    int32_t** parent; // parent[n][i] = j === (n-1)(j) is the parent of (n)(i)
    
    static int32_t NO_PARENT;
};


// =====================================================================
//                          Accessors' definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================



inline ae_parent_tree::ae_parent_tree( void )
{
  //  printf("init parent tree (run_for_gens = %ld, init_pop_size = %ld)\n", ae_common :: run_for_gens, ae_common :: init_pop_size);
  parent = new int32_t*[ae_common::run_for_gens];
  //  printf("init 2\n");
  for( int32_t n = 0 ; n < ae_common::run_for_gens ; n++ )
  {
    parent[n] = new int32_t[ae_common::init_pop_size];
  }
  //  printf("parent tree allocated\n");
  for( int32_t i = 0 ; i < ae_common::init_pop_size ; i++ )
  {
    parent[0][i] = NO_PARENT;
  }
  //  printf("parent tree initialized\n");
  return;
}
  
inline void ae_parent_tree::insert_indiv( int32_t n, int32_t i, int32_t j )
{
  // insert individual i from generation n with parent j from generation n-1
  parent[n][i] = j;
}

inline int32_t ae_parent_tree::get_LCA( int32_t n, int32_t a, int32_t b, int32_t* j = NULL )
{
  // return the generation of the last common ancestor of a and b. stores its index in j if provided.
  while ( a != b )
  {
    a = parent[n][a];
    b = parent[n][b];
    n = n-1;
  }
  if ( j != NULL )
  {
    *j = a;
  }
  
  return n;
}

#endif // __AE_PARENT_TREE_
