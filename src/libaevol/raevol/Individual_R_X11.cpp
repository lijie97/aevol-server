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
#include "Individual_R_X11.h"

namespace aevol {

//##############################################################################
//                                                                             #
//                           Class ae_individual_R_X11                         #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
Individual_R_X11::Individual_R_X11( const Individual_R_X11 &model, bool replication_report_copy  ) :
  Individual( model, replication_report_copy  ), Individual_R( model, replication_report_copy  ), Individual_X11( model, replication_report_copy  )
{
  //printf("ae_individual_R_X11( model )");
}

Individual_R_X11::Individual_R_X11( void )  :
Individual(), Individual_R(), Individual_X11()
{
  //printf("ae_individual_R_X11( void )");
}

Individual_R_X11::Individual_R_X11( Individual_R_X11* parent, int32_t id,
                                          ae_jumping_mt* mut_prng, ae_jumping_mt* stoch_prng ) :
        Individual( parent, id, mut_prng, stoch_prng ),
        Individual_R( parent, id, mut_prng, stoch_prng  ),
        Individual_X11( parent, id, mut_prng, stoch_prng  )
{
  //printf("ae_individual_R_X11( parent )");
}

Individual_R_X11::Individual_R_X11( gzFile backup_file ) :
Individual( backup_file ), Individual_R( backup_file ), Individual_X11( backup_file )
{
}

// =================================================================
//                             Destructors
// =================================================================
Individual_R_X11::~Individual_R_X11( void )
{
}

// =================================================================
//                            Public Methods
// =================================================================

void ae_individual_R_X11::display_concentrations( ae_X11_window* win, const std::vector<ae_environment*>& env_list )
{
  char* color = new char[8];
  char* color2 = NULL;
  strcpy( color, "#FFFFFF" );

  //variable qui sert à stocker du texte à afficher
  char display_string[40];
  int16_t nb_prot = 0;
  int16_t life_time = ae_common::individual_evaluation_dates->get_value( ae_common::individual_evaluation_nbr - 1 );

  //deux pointeurs utilisés pour défiler dans la liste de protéines
  //  ae_list_node* prot_node  = NULL;
  //  ae_protein_R* prot       = NULL;

  //Draw the number of proteins
  if (ae_common::with_heredity)
  {
    nb_prot = _protein_list.size();
    sprintf( display_string, "Nb proteins: %"PRId32" (inherited: %"PRId32")",nb_prot,_inherited_protein_list.size());
  }
  else
  {
    nb_prot = _protein_list.size();
    sprintf( display_string, "Nb proteins: %"PRId32" (without heredity)",nb_prot);
  }
  win->draw_string( 15, 15, display_string );

  //Draw the life time
  sprintf( display_string, "Life duration: %"PRId32"",life_time);
  win->draw_string( 15, 30, display_string );

  //Draw the box
  win->draw_line( win->get_width() / 10 , 2* win->get_height() / 10 , 9 * win->get_width() / 10 , 2 * win->get_height() / 10 ,color);
  win->draw_line( win->get_width() / 10 , 9 * win->get_height() / 10 , 9 * win->get_width() / 10 , 9 * win->get_height() / 10 ,color);
  win->draw_line( win->get_width() / 10 , 2* win->get_height() / 10 ,  win->get_width() / 10 , 9 * win->get_height() / 10 ,color);
  win->draw_line( 9 * win->get_width() / 10 , 2 * win->get_height() / 10 , 9 * win->get_width() / 10 , 9 * win->get_height() / 10 ,color);

  // save the initial list of proteins
  std::vector<ae_protein*> init_prot_list = _protein_list;
  //  init_prot_list->add_list(_protein_list);
  // Add the signals to the list of proteins
  ae_environment* envir = NULL;
  //  std::vector<ae_protein_R*> _signals;
  std::map<int,std::vector<ae_protein_R*>*> _cloned_signals;
  ae_protein_R* cloned_signal = NULL;
  std::vector<int> env_switch(env_list.size());

  for(int i = 0; i < env_list.size(); i++)
  {
    envir = env_list[i];
    env_switch[i] = envir->get_id();
    //     printf("Env at age %d : %d (%d)\n",i,envir->get_id(),env_list.size());

    if (_cloned_signals.find(envir->get_id()) == _cloned_signals.end()) {
      std::vector<ae_protein_R*>* loc_signals = new std::vector<ae_protein_R*>();
      for ( int8_t j = 0; j < envir->get_signals().size(); j++)
      {
        cloned_signal = new ae_protein_R(NULL,*(envir->_signals[j]));
        cloned_signal->set_concentration(0.0);
        _protein_list.push_back(cloned_signal);
        loc_signals->push_back(cloned_signal);
      }
      _cloned_signals.insert(std::make_pair(envir->get_id(),loc_signals));
    }
  }

  //set the concentrations of proteins to their initial value
  double* concentrations = new double[_protein_list.size()]; // initialise le tableau de concentrations.
  //  int16_t prot_index = 0;
  for (int i = 0; i < _protein_list.size(); i++) {
    if (! ((ae_protein_R*)_protein_list[i])->is_signal() ) ((ae_protein_R*)_protein_list[i])->reset_concentration();
    concentrations[i] = ((ae_protein_R*)_protein_list[i])->get_concentration();
    //    prot_node = prot_node->get_next();
  }

  // compute steps
  double x_step = 0.8 * win->get_width() / (double)(life_time / ae_common::degradation_step);
  double y_step = 0.7 * win->get_height();

  // Go from an evaluation date to the next
  int8_t compteur_env   = 0;
  envir = env_list[compteur_env];
  int16_t nb_signals = 0;

  //Add the signals protein to the individual
  for ( int8_t i = 0; i < _cloned_signals[env_switch[compteur_env]]->size(); i++)
  {
    ((ae_protein_R*) _cloned_signals[env_switch[compteur_env]]->at(i))->set_concentration(0.9);
  }

//	printf("Nb proteins : %d\n",nb_prot);

//	  ae_list_node* rna_node  = _rna_list->get_first();
//	  ae_rna*       rna       = NULL;
//
//	  while ( rna_node != NULL )
//	  {
//	    rna = (ae_rna*) rna_node->get_obj();
//
//	    if ( rna->is_coding() == true )
//	    {
//	    	printf("RNA : %f\n",rna->get_basal_level());
//	    }
//
//	    rna_node = rna_node->get_next();
//	  }

//	for (int i = 0; i < _rna_list_coding.size(); i++) {
//
//	}

  for(int16_t indiv_age = 0 ; indiv_age < life_time ; indiv_age+=0 )
  {
    //Updating the concentrations in order to respect the degradation step.
    for( int16_t i = 0; i < 1/ae_common::degradation_step; i++ )
    {
      update_concentrations();

      //affichage des points n+1 dans la concentration
      //prot_index = 0;
      for (int proti = 0; proti < _protein_list.size(); proti++) {

        // morceau ajouté pour colorer les protéines en fonctions de leur paramètres
        if ( ((ae_protein_R*)_protein_list[proti])->get_is_functional() )
        {
          color2 = ae_X11_window::get_color( ((ae_protein_R*)_protein_list[proti])->get_mean() );
        }
        else
        {
          color2 = new char[8];
          strcpy( color2, "#FFFFFF" );
        }

//				if ( ((ae_protein_R*)_protein_list[proti])->is_signal() ) {
//					printf("Time %d :(%d): %f %f\n",indiv_age*10+i,concentrations[proti],
//							proti,
//							((ae_protein_R*)_protein_list[proti])->get_concentration());
//				}
        win->draw_line( (int16_t)((win->get_width() / 10) + (((indiv_age/ae_common::degradation_step)+i)*x_step))  ,
                        (int16_t)(( 9 * win->get_height() / 10)-(concentrations[proti]*y_step)) ,
                        (int16_t)((win->get_width() / 10) + ((((indiv_age / ae_common::degradation_step)+i) + 1) * x_step)) ,
                        (int16_t)((9 * win->get_height() / 10)-(((ae_protein_R*)_protein_list[proti])->get_concentration()*y_step)) ,color2);
        concentrations[proti]=((ae_protein_R*)_protein_list[proti])->get_concentration();
        //	    prot_node = prot_node->get_next();
        delete[] color2;
      }

    }
    indiv_age+=1;

    if( ae_common::individual_environment_dates->search(indiv_age) != -1)
    {
      //Remove the signals of this environment
      for ( int8_t i = 0; i < _cloned_signals[env_switch[compteur_env]]->size(); i++)
      {
        _cloned_signals[env_switch[compteur_env]]->at(i)->set_concentration(0.0);
      }
      //      for ( int8_t i = 0; i < _signals->get_nb_elts(); i++)
      //      {
      //	((ae_protein_R*) _signals->get_object(i))->set_concentration(0.);
      //      }

      // Change the environment at is next value
      compteur_env+=1;
      envir = env_list[compteur_env];

      // Add the signals of this new environment
      for ( int8_t i = 0; i < _cloned_signals[env_switch[compteur_env]]->size(); i++)
      {
        _cloned_signals[env_switch[compteur_env]]->at(i)->set_concentration(0.9);
      }
      //      for ( int8_t i = 0; i < _signals->get_nb_elts(); i++)
      //      {
      //	((ae_protein_R*) _signals->get_object(i))->set_concentration(0.9);
      //      }
    }
  }

  //Remove the signals of the last environment
  for ( int8_t i = 0; i < _cloned_signals[env_switch[compteur_env]]->size(); i++)
  {
    _cloned_signals[env_switch[compteur_env]]->at(i)->set_concentration(0.0);
  }
  //  for ( int8_t i = 0; i < _signals->get_nb_elts(); i++)
  //  {
  //    ((ae_protein_R*) _signals->get_object(i))->set_concentration(0.0);
  //  }

  _protein_list.clear();
  _protein_list = init_prot_list;
  init_prot_list.clear();

  delete[] concentrations;
  //  delete init_prot_list;
  delete[] color;
}
// =================================================================
//                           Protected Methods
// =================================================================
} // namespace aevol
