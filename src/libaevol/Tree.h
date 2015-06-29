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
 
 
#ifndef AEVOL_TREE_H__
#define AEVOL_TREE_H__


// =================================================================
//                              Includes
// =================================================================
#include <inttypes.h>

#include <cstdio>
#include <cstdlib>
#include <cassert>

#include "ReplicationReport.h"
#include "Observer.h"
#include "ae_enums.h"


namespace aevol {


// =================================================================
//                          Class declarations
// =================================================================
class ExpManager;


class Tree : public Observer
{  
  public :
    
    // =================================================================
    //                             Constructors
    // =================================================================
    Tree() = delete;
    Tree(const Tree &model) = delete;
    // To be used when we want to run a simulation.
    Tree(ExpManager* exp_m, TreeMode tree_mode, int64_t tree_step);
    // To be used when we want to INSPECT a tree,
    // not when we want to run a simulation.
    Tree(ExpManager* exp_m, char* tree_file_name);
    
    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~Tree(void);

    // =================================================================
    //                        Accessors: getters
    // =================================================================
    inline int64_t  get_tree_step(void) const;
    inline TreeMode get_tree_mode(void) const;
    
    // Precondition for the following methods:
    // the tree was emptied every TREE_STEP generations ==> it contains
    // only the last generations since the last emptying ==> do not ask
    // something about an older generation
    ReplicationReport** get_reports(int64_t t) const;
    ReplicationReport* get_report_by_index(int64_t t, int32_t index) const;
    ReplicationReport* get_report_by_rank(int64_t t, int32_t rank) const;
  

    // =================================================================
    //                        Accessors: setters
    // =================================================================
    
    // =================================================================
    //                            Public Methods
    // =================================================================
    void signal_end_of_generation();
    void write_to_tree_file(gzFile tree_file);
    
    /** Returns the date of birth of the last common ancestor of individuals
     * of index a and b at time t and stores its index in j if provided
     * (not NULL)
     */
    inline int64_t get_LCA(int64_t t, int32_t a, int32_t b, int32_t *j);

  void update(Observable& o, void* arg) override {

  };
    

    // =================================================================
    //                           Public Attributes
    // =================================================================
    static const int32_t NO_PARENT;
    
  
  
  
  
  protected :
    // =================================================================
    //                           Protected Methods
    // =================================================================
    
    // =================================================================
    //                          Protected Attributes
    // =================================================================
    ExpManager* _exp_m;
    
    int64_t _tree_step;
    TreeMode _tree_mode;
    
    ReplicationReport*** _replics;
    // Two-dimensional table of ReplicationReport*
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
inline int64_t Tree::get_tree_step(void) const
{
  return _tree_step;
}

inline TreeMode Tree::get_tree_mode(void) const
{
  return _tree_mode;
}

// =====================================================================
//                           Setters' definitions
// =====================================================================


// =====================================================================
//                       Inline functions' definition
// =====================================================================
inline int64_t Tree::get_LCA(int64_t t, int32_t a, int32_t b, int32_t* j = NULL)
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

#endif // AEVOL_TREE_H__