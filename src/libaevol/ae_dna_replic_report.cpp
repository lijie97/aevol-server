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




// =================================================================
//                              Libraries
// =================================================================



// =================================================================
//                            Project Files
// =================================================================
#include <ae_dna_replic_report.h>
#include <ae_mutation.h>

namespace aevol {

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
  _HT               = new ae_list<ae_mutation*>();

  _nb_mut[SWITCH] = 0;
  _nb_mut[S_INS]  = 0;
  _nb_mut[S_DEL]  = 0;
  _nb_mut[DUPL]   = 0;
  _nb_mut[DEL]    = 0;
  _nb_mut[TRANS]  = 0;
  _nb_mut[INV]    = 0;
  _nb_mut[INS_HT] = 0;
  _nb_mut[REPL_HT]= 0;
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
  _nb_mut[INS_HT] = model._nb_mut[INS_HT];
  _nb_mut[REPL_HT]= model._nb_mut[REPL_HT];
  
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
  
  _HT = new ae_list<ae_mutation*>();
  ae_list_node<ae_mutation*>* HT_node = (model._HT)->get_first();
  ae_mutation* HT = NULL;
  while ( HT_node != NULL )
  {
    HT = HT_node->get_obj();
    _HT->add( new ae_mutation( *HT ) );
    HT_node = HT_node->get_next();
  }
}

// =================================================================
//                             Destructors
// =================================================================
ae_dna_replic_report::~ae_dna_replic_report( void )
{
  _mutations->erase( true );
  _rearrangements->erase( true );
  _HT->erase( true );
  
  delete _mutations;
  delete _rearrangements;
  delete _HT;
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
  _nb_mut[INS_HT] = 0;
  _nb_mut[REPL_HT]= 0;
  
  ae_list_node<ae_mutation*>* event_node  = _HT->get_first();
  ae_mutation* event       = NULL;
  while ( event_node != NULL )
  {
    event = event_node->get_obj();
    assert( event->get_mut_type() == INS_HT || event->get_mut_type() == REPL_HT );
    
    _nb_mut[ event->get_mut_type() ]++;
    
    event_node = event_node->get_next();
  }
  
  event_node = _rearrangements->get_first();
  event = NULL;
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
} // namespace aevol
