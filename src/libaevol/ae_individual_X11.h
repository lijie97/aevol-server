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
 
 
 #ifndef __AE_INDIVIDUAL_X11_H__
#define  __AE_INDIVIDUAL_X11_H__
 
 
// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>
#include <zlib.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_individual.h>
#include <ae_X11_window.h>




// =================================================================
//                          Class declarations
// =================================================================





 
class ae_individual_X11 : public virtual ae_individual
{  
  public :
  
    // =================================================================
    //                             Constructors
    // =================================================================
    ae_individual_X11(  ae_exp_manager* exp_m,
                        ae_jumping_mt* alea, 
                        ae_params_mut* param_mut,
                        double w_max,
                        int32_t min_genome_length,
                        int32_t max_genome_length,
                        bool allow_plasmids,
                        int32_t plasmid_minimal_length,
                        int32_t id,
                        int32_t age );
  
    ae_individual_X11( ae_exp_manager* exp_manager, gzFile* backup_file );
    ae_individual_X11( const ae_individual_X11 &model );
    ae_individual_X11( ae_individual_X11* const parent, int32_t id );
    //~ ae_individual_X11( char* genome, int32_t genome_size );
  
    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_individual_X11( void );
  
    // =================================================================
    //                              Accessors
    // =================================================================
  
    // =================================================================
    //                            Public Methods
    // =================================================================
    virtual void display( void );
    virtual void display_cdss( ae_X11_window* win );
    virtual void display_rnas( ae_X11_window* win );
  
    // =================================================================
    //                           Public Attributes
    // =================================================================
  
  
  
  
  
  protected :
  
    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    ae_individual_X11( void )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
  
    // =================================================================
    //                           Protected Methods
    // =================================================================
    void reset_sectors( void );
    void add_layer( void );
    void init_occupied_sectors( void );
  
    // =================================================================
    //                          Protected Attributes
    // =================================================================
    // These are used to manage overlapping CDS and RNA display
    int16_t _outmost_layer;
    bool*   _occupied_sectors[2][100];  // TODO : find a way to manage this table's size properly?
};


// =====================================================================
//                          Accessors definitions
// =====================================================================


#endif // __AE_INDIVIDUAL_X11_H__
