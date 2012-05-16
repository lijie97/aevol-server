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
 
 
 #ifndef __AE_FUZZY_SET_X11_H__
#define  __AE_FUZZY_SET_X11_H__
 
 
// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_fuzzy_set.h>
#include <ae_X11_window.h>




// =================================================================
//                          Class declarations
// =================================================================





 
class ae_fuzzy_set_X11 : public ae_fuzzy_set
{  
  public :
  
    // =================================================================
    //                             Constructors
    // =================================================================
    ae_fuzzy_set_X11( void );
  
    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_fuzzy_set_X11( void );
  
    // =================================================================
    //                              Accessors
    // =================================================================
  
    // =================================================================
    //                            Public Methods
    // =================================================================
    virtual void display( ae_X11_window* win, color_map color, bool fill = false, bool bold = false );
  
    // =================================================================
    //                           Public Attributes
    // =================================================================
  
  
  
  
  
  protected :
  
    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    //~ ae_fuzzy_set_X11( void )
    //~ {
      //~ printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      //~ exit( EXIT_FAILURE );
    //~ };
    ae_fuzzy_set_X11( const ae_fuzzy_set &model )
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
};


// =====================================================================
//                          Accessors definitions
// =====================================================================


#endif // __AE_FUZZY_SET_X11_H__
