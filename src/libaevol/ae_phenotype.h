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
 
 
#ifndef __AE_PHENOTYPE_H__
#define __AE_PHENOTYPE_H__
 
 
// =================================================================
//                              Libraries
// =================================================================
#include <cinttypes>



// =================================================================
//                            Project Files
// =================================================================
#include "fuzzy.h"

namespace aevol {

// =================================================================
//                          Class declarations
// =================================================================
class ae_individual;




 
class ae_phenotype : public Fuzzy
{
  public :
  
    // =================================================================
    //                             Constructors
    // =================================================================
    ae_phenotype(void);
    ae_phenotype(const ae_phenotype&) = delete;
    ae_phenotype(ae_individual* indiv);
    ae_phenotype(ae_individual* indiv, const ae_phenotype &model);
  
    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_phenotype(void) = default;
  
    // =================================================================
    //                              Accessors
    // =================================================================
  
    // =================================================================
    //                            Public Methods
    // =================================================================
  
    // =================================================================
    //                           Public Attributes
    // =================================================================
  
  
  
  
  
  protected :
    // =================================================================
    //                           Protected Methods
    // =================================================================
  
    // =================================================================
    //                          Protected Attributes
    // =================================================================
    ae_individual* _indiv;
};


// =====================================================================
//                          Accessors definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================

} // namespace aevol

#endif // __AE_PHENOTYPE_H__
