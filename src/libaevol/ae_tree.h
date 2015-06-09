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
 
 
#ifndef __AE_TREE_H__
#define __AE_TREE_H__


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>



// =================================================================
//                            Project Files
// =================================================================
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "ae_replication_report.h"
#include "ae_enums.h"


namespace aevol {


// =================================================================
//                          Class declarations
// =================================================================
class ae_exp_manager;


class ae_tree
{  
  public :
    
    // =================================================================
    //                             Constructors
    // =================================================================
    ae_tree(ae_exp_manager* exp_m, ae_tree_mode tree_mode, int64_t tree_step);
    // To be used when we want to run a simulation.
    ae_tree(ae_exp_manager* exp_m, char* tree_file_name);
    // To be used when we want to INSPECT a tree, 
    // not when we want to run a simulation.
    
    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_tree(void);

    // =================================================================
    //                        Accessors: getters
    // =================================================================
    inline int64_t       get_tree_step(void) const;
    inline ae_tree_mode  get_tree_mode(void) const;
    
    // Precondition for the following 3 methods: 
    // the tree was emptied every TREE_STEP generations ==> it contains
    // only the last generations since the last emptying ==> do not ask
    // something about an older generation 
    int32_t get_nb_indivs(int64_t t) const;
    ae_replication_report* get_report_by_index(int64_t t, int32_t index) const;
    ae_replication_report* get_report_by_rank(int64_t t, int32_t rank) const;
  

    // =================================================================
    //                        Accessors: setters
    // =================================================================
    void set_replic_report(int32_t id, ae_replication_report* replic_report);
    void set_replic_report(int64_t t, int32_t id, ae_replication_report* replic_report);
    void set_nb_indivs(int32_t nb_indivs, int64_t t);
    
    
    // =================================================================
    //                            Public Methods
    // =================================================================
    void fill_tree_with_cur_gener(void);
    void write_to_tree_file(gzFile tree_file);
    
    /** Returns the date of birth of the last common ancestor of individuals
     * of index a and b at time t and stores its index in j if provided
     * (not NULL)
     */
    inline int64_t get_LCA(int64_t t, int32_t a, int32_t b, int32_t *j);
    

    // =================================================================
    //                           Public Attributes
    // =================================================================
    static const int32_t NO_PARENT;
    
  
  
  
  
  protected :
    
    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    ae_tree( void )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
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
    ae_exp_manager* _exp_m;
    
    int64_t       _tree_step;
    ae_tree_mode  _tree_mode;
    int32_t*      _nb_indivs;
    
    ae_replication_report*** _replics;
    // Two-dimensional table of ae_replication_report*
    //    dimension 1 (lines)   : generation
    //    dimension 2 (columns) : individual
    //
    // !!!!! WARNING !!!!!
    // The report at line l, column c is for the
    // replication that created the indiv with index c of generation l+1
    
    // light tree representation
    int32_t** _parent;
};


// =====================================================================
//                           Getters' definitions
// =====================================================================
inline int64_t ae_tree::get_tree_step(void) const
{
  return _tree_step;
}

inline ae_tree_mode ae_tree::get_tree_mode(void) const
{
  return _tree_mode;
}

// =====================================================================
//                           Setters' definitions
// =====================================================================


// =====================================================================
//                       Inline functions' definition
// =====================================================================
inline int64_t ae_tree::get_LCA(int64_t t, int32_t a, int32_t b, int32_t* j = NULL)
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
        a = _parent[t][a];
        b = _parent[t][b];
        t--;
      }
      
      if (j != NULL)
      {
        *j = a;
      }
      
      break;
    }
  }
  
  return t;
}

} // namespace aevol

#endif // __AE_TREE_H__
