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


/*! \class ae_exp_setup
    \brief This class contains the experiment setup among which the environment
*/


// =================================================================
//                              Libraries
// =================================================================
#include <stdio.h>
#include <zlib.h>
#include <sys/types.h>
#include <sys/stat.h>




// =================================================================
//                            Project Files
// =================================================================
#include <ae_exp_setup.h>
#include <ae_population.h>

#ifdef __X11
  #include <ae_population_X11.h>
#endif




//##############################################################################
//                                                                             #
//                              Class ae_exp_setup                             #
//                                                                             #
//##############################################################################

// ===========================================================================
//                         Definition of static attributes
// ===========================================================================

// ===========================================================================
//                                  Constructors
// ===========================================================================
ae_exp_setup::ae_exp_setup( ae_exp_manager* exp_m )
{
  _exp_m  = exp_m;
  
  // -------------------------------------------------------------- Selection
  _sel = new ae_selection( exp_m );
  
  // --------------------------------------------------------------- Transfer
  _with_HT          = false;
  _HT_ins_rate      = 0.0;
  _HT_repl_rate     = 0.0;
  
  // --------------------------------------------------------------- Plasmids
  _with_plasmids    = false;
  _prob_plasmid_HT  = 0.0;
  _tune_donor_ability     = 0.0;
  _tune_recipient_ability = 0.0;
  _donor_cost       = 0.0;
  _recipient_cost   = 0.0;
  _swap_GUs         = false;
  
  // -------------------------------------------------------------- Secretion
  _with_secretion = false;
  _secretion_contrib_to_fitness = 0.0;
  _secretion_cost               = 0.0;
}
  

// ===========================================================================
//                                 Destructor
// ===========================================================================
ae_exp_setup::~ae_exp_setup( void )
{
  delete _sel;
}

// ===========================================================================
//                                 Public Methods
// ===========================================================================
/*!
*/
void ae_exp_setup::write_setup_file( gzFile exp_setup_file ) const
{
  // --------------------------------------------------------------- Transfer
  int8_t tmp_with_HT = _with_HT;
  gzwrite( exp_setup_file, &tmp_with_HT, sizeof(tmp_with_HT) );
  if ( _with_HT )
  {
    gzwrite( exp_setup_file, &_HT_ins_rate,  sizeof(_HT_ins_rate) );
    gzwrite( exp_setup_file, &_HT_repl_rate, sizeof(_HT_repl_rate) );
  }
  
  // --------------------------------------------------------------- Plasmids
  int8_t tmp_with_plasmids = get_with_plasmids();
  gzwrite( exp_setup_file, &tmp_with_plasmids, sizeof(tmp_with_plasmids) );
  if ( tmp_with_plasmids )
  {
    gzwrite( exp_setup_file, &_prob_plasmid_HT,  sizeof(_prob_plasmid_HT) );
    gzwrite( exp_setup_file, &_tune_donor_ability,  sizeof(_tune_donor_ability) );
    gzwrite( exp_setup_file, &_tune_recipient_ability,  sizeof(_tune_recipient_ability) );
    gzwrite( exp_setup_file, &_donor_cost,  sizeof(_donor_cost) );
    gzwrite( exp_setup_file, &_recipient_cost,  sizeof(_recipient_cost) );
    int8_t tmp_swap_GUs = _swap_GUs;
    gzwrite( exp_setup_file, &tmp_swap_GUs, sizeof(tmp_swap_GUs) );
  }
  
  // -------------------------------------------------------------- Secretion
  int8_t tmp_with_secretion = _with_secretion;
  gzwrite( exp_setup_file, &tmp_with_secretion, sizeof(tmp_with_secretion) );
  if ( _with_secretion )
  {
    gzwrite( exp_setup_file, &_secretion_contrib_to_fitness, sizeof(_secretion_contrib_to_fitness) );
    gzwrite( exp_setup_file, &_secretion_cost, sizeof(_secretion_cost) );
  }
  
  get_sel()->write_setup_file( exp_setup_file );
}

void ae_exp_setup::write_setup_file( FILE* exp_setup_file ) const
{
  // --------------------------------------------------------------- Transfer
  //...
  
  // -------------------------------------------------------------- Secretion
  //...
  
  get_sel()->write_setup_file( exp_setup_file );
}

void ae_exp_setup::load( gzFile setup_file, gzFile backup_file, bool verbose )
{
  // -------------------------------------------- Retrieve transfer parameters
  int8_t tmp_with_HT;
  gzread( setup_file, &tmp_with_HT, sizeof(tmp_with_HT) );
  _with_HT = tmp_with_HT ? 1 : 0;
  if ( _with_HT )
  {
    gzread( setup_file, &_HT_ins_rate,  sizeof(_HT_ins_rate) );
    gzread( setup_file, &_HT_repl_rate, sizeof(_HT_repl_rate) );
  }
  
  // -------------------------------------------- Retrieve plasmid parameters
  int8_t tmp_with_plasmids;
  gzread( setup_file, &tmp_with_plasmids, sizeof(tmp_with_plasmids) );
  _with_plasmids = tmp_with_plasmids ? 1 : 0;
  if ( _with_plasmids )
  {
    gzread( setup_file, &_prob_plasmid_HT,  sizeof(_prob_plasmid_HT) );
    gzread( setup_file, &_tune_donor_ability,  sizeof(_tune_donor_ability) );
    gzread( setup_file, &_tune_recipient_ability,  sizeof(_tune_recipient_ability) );
    gzread( setup_file, &_donor_cost,  sizeof(_donor_cost) );
    gzread( setup_file, &_recipient_cost,  sizeof(_recipient_cost) );
    int8_t tmp_swap_GUs;
    gzread( setup_file, &tmp_swap_GUs, sizeof(tmp_swap_GUs) );
    _swap_GUs = tmp_swap_GUs ? 1 : 0;
  }
  
  // ------------------------------------------ Retrieve secretion parameters
  int8_t tmp_with_secretion;
  gzread( setup_file, &tmp_with_secretion, sizeof(tmp_with_secretion) );
  _with_secretion = tmp_with_secretion ? true : false;
  if ( _with_secretion )
  {
    gzread( setup_file, &_secretion_contrib_to_fitness, sizeof(_secretion_contrib_to_fitness) );
    gzread( setup_file, &_secretion_cost, sizeof(_secretion_cost) );
  }
  
  // ---------------------------------------------- Retrieve selection context
  get_sel()->load( setup_file, backup_file, verbose );
}

void ae_exp_setup::load( FILE* setup_file, gzFile backup_file, bool verbose )
{
  printf( "Plain text setup file support not implemented yet (sorry)\n" );
  exit( EXIT_FAILURE );
  
  // -------------------------------------------- Retrieve transfer parameters
  int8_t tmp_with_HT;
  gzread( setup_file, &tmp_with_HT, sizeof(tmp_with_HT) );
  _with_HT = tmp_with_HT ? 1 : 0;
  if ( _with_HT )
  {
    gzread( setup_file, &_HT_ins_rate,  sizeof(_HT_ins_rate) );
    gzread( setup_file, &_HT_repl_rate, sizeof(_HT_repl_rate) );
  }
  int8_t tmp_with_plasmids;
  gzread( setup_file, &tmp_with_plasmids, sizeof(tmp_with_plasmids) );
  _with_plasmids = tmp_with_plasmids ? 1 : 0;
  if ( tmp_with_plasmids )
  {
    gzread( setup_file, &_prob_plasmid_HT,  sizeof(_prob_plasmid_HT) );
    gzread( setup_file, &_tune_donor_ability,  sizeof(_tune_donor_ability) );
    gzread( setup_file, &_tune_recipient_ability,  sizeof(_tune_recipient_ability) );
    gzread( setup_file, &_donor_cost,  sizeof(_donor_cost) );
    gzread( setup_file, &_recipient_cost,  sizeof(_recipient_cost) );
    int8_t tmp_swap_GUs;
    gzread( setup_file, &tmp_swap_GUs, sizeof(tmp_swap_GUs) );
    _swap_GUs = tmp_swap_GUs ? 1 : 0;
  }
  
  // ------------------------------------------ Retrieve secretion parameters
  int8_t tmp_with_secretion;
  gzread( setup_file, &tmp_with_secretion, sizeof(tmp_with_secretion) );
  _with_secretion = tmp_with_secretion ? true : false;
  if ( _with_secretion )
  {
    gzread( setup_file, &_secretion_contrib_to_fitness, sizeof(_secretion_contrib_to_fitness) );
    gzread( setup_file, &_secretion_cost, sizeof(_secretion_cost) );
  }
  
  // ---------------------------------------------- Retrieve selection context
  printf( "  Loading selection context..." );
  fflush( stdout );
  get_sel()->load( setup_file, backup_file, verbose );
  printf( "OK\n" );
}




// ===========================================================================
//                                Protected Methods
// ===========================================================================
