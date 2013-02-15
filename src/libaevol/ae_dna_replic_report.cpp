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
 
 
// =================================================================
//                              Libraries
// =================================================================



// =================================================================
//                            Project Files
// =================================================================
#include <ae_dna_replic_report.h>
#include <ae_mutation.h>




//##############################################################################
//                                                                             #
//                          Class ae_dna_replic_report                         #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
ae_dna_replic_report::ae_dna_replic_report( void )
{
  _mutations        = new ae_list<ae_mutation*>();
  _rearrangements   = new ae_list<ae_mutation*>();

  _nb_mut[SWITCH] = 0;
  _nb_mut[S_INS]  = 0;
  _nb_mut[S_DEL]  = 0;
  _nb_mut[DUPL]   = 0;
  _nb_mut[DEL]    = 0;
  _nb_mut[TRANS]  = 0;
  _nb_mut[INV]    = 0;
}

ae_dna_replic_report::ae_dna_replic_report( const ae_dna_replic_report &model )
{
  // Initialize trivial attributes
  _nb_mut[SWITCH] = model._nb_mut[SWITCH];
  _nb_mut[S_INS]  = model._nb_mut[S_INS];
  _nb_mut[S_DEL]  = model._nb_mut[S_DEL];
  _nb_mut[DUPL]   = model._nb_mut[DUPL];
  _nb_mut[DEL]    = model._nb_mut[DEL];
  _nb_mut[TRANS]  = model._nb_mut[TRANS];
  _nb_mut[INV]    = model._nb_mut[INV];
  
  
  // Copy _mutations and _rearrangements
  _mutations = new ae_list<ae_mutation*>();
  ae_list_node<ae_mutation*>* mut_node = (model._mutations)->get_first();
  ae_mutation* mut = NULL;
  while ( mut_node != NULL )
  {
    mut = mut_node->get_obj();
    _mutations->add( new ae_mutation( *mut ) );
    mut_node = mut_node->get_next();
  }
  
  _rearrangements = new ae_list<ae_mutation*>();
  ae_list_node<ae_mutation*>* rear_node = (model._rearrangements)->get_first();
  ae_mutation* rear = NULL;
  while ( rear_node != NULL )
  {
    rear = rear_node->get_obj();
    _rearrangements->add( new ae_mutation( *rear ) );
    rear_node = rear_node->get_next();
  }
}

// =================================================================
//                             Destructors
// =================================================================
ae_dna_replic_report::~ae_dna_replic_report( void )
{
  _mutations->erase( true );
  _rearrangements->erase( true );
  
  delete _mutations;
  delete _rearrangements;
}



// =================================================================
//                            Public Methods
// =================================================================

void ae_dna_replic_report::compute_stats( void )
{
  // useful when we inspect a tree file
  // because stats are not saved in the file 
  _nb_mut[SWITCH] = 0;
  _nb_mut[S_INS]  = 0;
  _nb_mut[S_DEL]  = 0;
  _nb_mut[DUPL]   = 0;
  _nb_mut[DEL]    = 0;
  _nb_mut[TRANS]  = 0;
  _nb_mut[INV]    = 0;
  
  ae_list_node<ae_mutation*>* event_node = _rearrangements->get_first();
  ae_mutation* event = NULL;
  
  while ( event_node != NULL )
  {
    event = event_node->get_obj();
    assert( event->get_mut_type() == DUPL || event->get_mut_type() == DEL || 
            event->get_mut_type() == TRANS || event->get_mut_type() == INV );
    
    _nb_mut[ event->get_mut_type() ]++;
    
    event_node = event_node->get_next();
  }

  
  event_node  = _mutations->get_first();
  event       = NULL;
  
  while ( event_node != NULL )
  {
    event = event_node->get_obj();
    assert( event->get_mut_type() == SWITCH || 
            event->get_mut_type() == S_INS || event->get_mut_type() == S_DEL );
    
    _nb_mut[ event->get_mut_type() ]++;
    
    event_node = event_node->get_next();
  }
}



// =================================================================
//                           Protected Methods
// =================================================================
