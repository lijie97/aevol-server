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
 
 
#ifndef __AE_TREE_H__
#define  __AE_TREE_H__


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_object.h>
#include <ae_common.h>
#include <ae_replication_report.h>




// =================================================================
//                          Class declarations
// =================================================================


class ae_tree : public ae_object
{  
 public :
  
  // =================================================================
  //                             Constructors
  // =================================================================
  ae_tree( void );
  ae_tree( char* backup_file_name, char* tree_file_name ); 
  // To be called when we want to INSPECT a tree, 
  // not when we want to run a simulation.
  
  // =================================================================
  //                             Destructors
  // =================================================================
  virtual ~ae_tree( void );
  
  // =================================================================
  //                              Accessors
  // =================================================================
  
  
  
  // =================================================================
  //                            Public Methods
  // =================================================================
  
  // Precondition for both methods: 
  // the tree was emptied every TREE_STEP generations ==> it contains
  // only the last generations since the last emptying ==> do not ask
  // something about an older generation 
  int32_t get_nb_indivs( int32_t generation ) const; 
  ae_replication_report * get_report_by_index( int32_t generation, int32_t index ) const;
  ae_replication_report * get_report_by_rank( int32_t generation, int32_t rank ) const;
  
  void fill_tree_with_cur_gener( void );
  void write_to_tree_file( gzFile* tree_file );
  
  inline int32_t get_LCA( int32_t num_gener, int32_t a, int32_t b, int32_t *j );
  // Returns the generation of the last common ancestor of individuals number a and b
  // of generation num_gener and stores its index in j if provided (not NULL)
  

  // =================================================================
  //                           Public Attributes
  // =================================================================
  static const int32_t NO_PARENT;
  
  
  
  
  
 protected :
  
  // =================================================================
  //                         Forbidden Constructors
  // =================================================================
  /*     ae_tree( void ) */
  /*     { */
  /*       printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ ); */
  /*       exit( EXIT_FAILURE ); */
  /*     }; */
  
  
  ae_tree( const ae_tree &model )
  {
    printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
    exit( EXIT_FAILURE );
  };
  
  // =================================================================
  //                           Protected Methods
  // =================================================================
  
  // =================================================================
  //                          Protected Attributes
  // =================================================================

  ae_tree_mode  _tree_mode;
  int32_t*      _nb_indivs;
  
  ae_replication_report*** _replics;
  // Two-dimensional table of ae_replication_report*
  //    dimension 1 (lines)   : generation
  //    dimension 2 (columns) : individual
  //
  // !!!!! WARNING !!!!!
  // The report at line 0, column 5 is for the
  // replication that created the indiv with index 5 of generation 1 (not generation 0)
  
  // light tree representation
  int32_t** _parent;
};


// =====================================================================
//                          Accessors' definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================

inline int32_t ae_tree::get_LCA( int32_t num_gener, int32_t a, int32_t b, int32_t* j = NULL )
// Returns the generation of the last common ancestor of individuals number a and b
// of generation num_gener and stores its index in j if provided (not NULL)
{
  switch( _tree_mode )
  {
    case NORMAL :
    {
      // TODO
      break;
    }
    case LIGHT :
    {
      while ( a != b )
      {
        a = _parent[num_gener][a];
        b = _parent[num_gener][b];
        num_gener = num_gener - 1;
      }
      
      if ( j != NULL )
      {
        *j = a;
      }
      
      break;
    }
  }
  
  return num_gener;
}


#endif // __AE_TREE_H__
