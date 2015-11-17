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
#include "../ExpManager.h"
#include "../ExpManager_X11.h"
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
Individual_R_X11::Individual_R_X11( const Individual_R_X11 &model  ) :
  Individual( model  ), Individual_R( model  ), Individual_X11( model  )
{
  //printf("ae_individual_R_X11( model )");
}

Individual_R_X11::Individual_R_X11(ExpManager* exp_m,
                           std::shared_ptr<JumpingMT> mut_prng,
                           std::shared_ptr<JumpingMT> stoch_prng,
                           std::shared_ptr<MutationParams> param_mut,
                           double w_max,
                           int32_t min_genome_length,
                           int32_t max_genome_length,
                           bool allow_plasmids,
                           int32_t id,
                           const char* strain_name,
                           int32_t age) :
                           Individual(exp_m,mut_prng,stoch_prng,param_mut,w_max,min_genome_length,
                            max_genome_length,allow_plasmids,id,strain_name,age),
                           Individual_R(exp_m,mut_prng,stoch_prng,param_mut,w_max,min_genome_length,
                            max_genome_length,allow_plasmids,id,strain_name,age),
                           Individual_X11(exp_m,mut_prng,stoch_prng,param_mut,w_max,min_genome_length,
                            max_genome_length,allow_plasmids,id,strain_name,age)
{

}


Individual_R_X11::Individual_R_X11( Individual_R_X11* parent, int32_t id,
                                    std::shared_ptr<JumpingMT> mut_prng, std::shared_ptr<JumpingMT> stoch_prng ) :
        Individual( parent, id, mut_prng, stoch_prng ),
        Individual_R( parent, id, mut_prng, stoch_prng  ),
        Individual_X11( parent, id, mut_prng, stoch_prng  )
{
  //printf("ae_individual_R_X11( parent )");
}

Individual_R_X11::Individual_R_X11( ExpManager* exp_m, gzFile backup_file ) :
Individual( exp_m, backup_file ), Individual_R( exp_m, backup_file ), Individual_X11( exp_m, backup_file )
{
}

// =================================================================
//                             Destructors
// =================================================================
Individual_R_X11::~Individual_R_X11( void ) noexcept
{
}

// =================================================================
//                            Public Methods
// =================================================================

void Individual_R_X11::display_concentrations( X11Window* win )
{
  char* color = new char[8];
  char* color2 = NULL;
  strcpy( color, "#FFFFFF" );

  //variable qui sert à stocker du texte à afficher
  char display_string[40];
  int16_t nb_prot = 0;
  int16_t life_time =  _exp_m->get_exp_s()->get_nb_indiv_age();

  //deux pointeurs utilisés pour défiler dans la liste de protéines
  //  ae_list_node* prot_node  = NULL;
  //  ae_protein_R* prot       = NULL;

  //Draw the number of proteins
  if ( _exp_m->get_exp_s()->get_with_heredity())
  {
    nb_prot = _protein_list.size();
    sprintf( display_string, "Nb proteins: %" PRId32 " (inherited: %ld)",nb_prot,_inherited_protein_list.size());
  }
  else
  {
    nb_prot = _protein_list.size();
    sprintf( display_string, "Nb proteins: %" PRId32 " (without heredity)",nb_prot);
  }
  win->draw_string( 15, 15, display_string );

  //Draw the life time
  sprintf( display_string, "Life duration: %" PRId32 "",life_time);
  win->draw_string( 15, 30, display_string );

  //Draw the box
  win->draw_line( win->get_width() / 10 , 2* win->get_height() / 10 , 9 * win->get_width() / 10 , 2 * win->get_height() / 10 ,color);
  win->draw_line( win->get_width() / 10 , 9 * win->get_height() / 10 , 9 * win->get_width() / 10 , 9 * win->get_height() / 10 ,color);
  win->draw_line( win->get_width() / 10 , 2* win->get_height() / 10 ,  win->get_width() / 10 , 9 * win->get_height() / 10 ,color);
  win->draw_line( 9 * win->get_width() / 10 , 2 * win->get_height() / 10 , 9 * win->get_width() / 10 , 9 * win->get_height() / 10 ,color);

  // save the initial list of proteins
  //std::list<Protein*> init_prot_list = _protein_list;

  //set the concentrations of proteins to their initial value
  double* concentrations = new double[_protein_list.size()]; // initialise le tableau de concentrations.
  //  int16_t prot_index = 0;
  int i = 0;
  for (const auto& prot : _protein_list) {
    ((Protein_R*)prot)->reset_concentration();
    concentrations[i] = ((Protein_R*)prot)->get_concentration();
    i++;
  }

  // compute steps
  double x_step = 0.8 * win->get_width() / (double)(get_exp_m()->get_exp_s()->get_nb_indiv_age());
  double y_step = 0.7 * win->get_height();

  for (int i = 1; i <= get_exp_m()->get_exp_s()->get_nb_indiv_age(); i++) {
    for (int j = 0; j < get_exp_m()->get_exp_s()->get_nb_degradation_step(); j++)
      update_concentrations();

      //affichage des points n+1 dans la concentration
      //prot_index = 0;
      int proti = 0;
      for (const auto& prot : _protein_list) {

        // morceau ajouté pour colorer les protéines en fonctions de leur paramètres
        if ( ((Protein_R*)prot)->get_is_functional() )
        {
          color2 = X11Window::get_color( ((Protein_R*)prot)->get_mean() );
        }
        else
        {
          color2 = new char[8];
          strcpy( color2, "#FFFFFF" );
        }

        win->draw_line( (int16_t)((win->get_width() / 10) + (i*x_step)),
                        (int16_t)(( 9 * win->get_height() / 10)-(concentrations[proti]*y_step)),
                        (int16_t)((win->get_width() / 10) + (( i + 1)  * x_step)),
                        (int16_t)((9 * win->get_height() / 10)-(((Protein_R*)prot)->get_concentration()*y_step)) ,color2);
        concentrations[proti]=((Protein_R*)prot)->get_concentration();
        printf("%d -- %f (%f) | ",proti,((Protein_R*)prot)->get_concentration(),((Protein_R*)prot)->get_mean());

        delete[] color2;
        proti++;
      }
    printf("\n");
  }

  //_protein_list.clear();
  //_protein_list = init_prot_list;
  //init_prot_list.clear();

  delete[] concentrations;
  delete[] color;
}

void Individual_R_X11::display_regulation( X11Window* win )
{
  int16_t nb_activators = 0;
  int16_t nb_operators = 0;

  double min_merged_activator_activity = 10;
  double min_merged_operator_activity = 10;
  double max_merged_activator_activity = 0;
  double max_merged_operator_activity = 0;
  double min_activator_activity = 10;
  double min_operator_activity = 10;
  double max_activator_activity = 0;
  double max_operator_activity = 0;
  double mean_activator_activity = 0;
  double mean_operator_activity = 0;
  // Retreive the genetic unit corresponding to the main chromosome
  GeneticUnit* gen_unit = &_genetic_unit_list.front();
  int32_t genome_length = gen_unit->get_dna()->get_length();

  // draw color scale
  char *color_bar = new char[8];
  for ( int16_t i = 0 ; i < (win->get_width()/2) ; i++ )
  {
    sprintf( color_bar, "#%02x%02x%02x", 255 - (i * 255 * 2) / win->get_width(),0,0 );
    win->draw_line( i, win->get_height() * 19 / 20, i, win->get_height(), color_bar );
  }

  for ( int16_t i = (win->get_width()/2) ; i < (win->get_width()) ; i++ )
  {
    sprintf( color_bar, "#%02x%02x%02x",0, ((i - (win->get_width() / 2)) * 255 * 2 ) / win->get_width(), 0 );
    win->draw_line( i, win->get_height() * 19 / 20, i, win->get_height(), color_bar );
  }
  delete [] color_bar;


  char display_string[80];
  sprintf( display_string, "-");
  win->draw_string( 15, (win->get_height() * 19 / 20) - 5, display_string );
  sprintf( display_string, "0");
  win->draw_string( win->get_width()/2,  (win->get_height() * 19 / 20) - 5, display_string );
  sprintf( display_string, "+");
  win->draw_string( win->get_width()-15,  (win->get_height() * 19 / 20) - 5, display_string );


  // Compute display diameter according to genome length and window size
  int16_t win_size      = std::min( win->get_width(), win->get_height() );
  int16_t diam          = round( win_size * log( (double)genome_length ) / 16 );

  // Prevent diameter from getting greater than 2/3 of the window size
  if ( diam > 2 * win_size / 3 )
  {
    diam = 2 * win_size / 3;
  }

  // Compute coordinates of the upper-left corner of the containing square
  int16_t pos_x = (win->get_width() - diam) / 2;
  int16_t pos_y = (win->get_height() - diam) / 2;
  double len_link;

  char* color = new char[8];
  strcpy( color, "#FFFFFF" );

  // Draw main circle
  win->draw_circle( pos_x, pos_y, diam );


  // ---------------
  //  Draw each RNA
  // ---------------

  // NB : As we want OriC to be at the "top" of the circle and the orientation
  //      to be clockwise, the drawing angle (theta) will be given as
  //      (90 - alpha), alpha being the "classical" trigonometric angle
  int16_t alpha_rna_first,alpha_prot_first; // Angles of first and last transcribed bases from OriC (degrees)
  int16_t theta_rna_first, theta_prot_first; // Transposed angles on the trigonometric circle (degrees)
  // Same as above with precision = 1/64 degree
  int16_t alpha_rna_first_64,alpha_prot_first_64;
  int16_t theta_rna_first_64,theta_prot_first_64;
  int16_t pos_rna_x,pos_rna_y,pos_prot_x,pos_prot_y;
  int8_t nb_signals = 0;

  // search for max and min regulation values (in order to scale colors)

  for (const auto& rna: _rna_list_coding) {
    // ---------------
    //  Draw each regulation link
    // ---------------

    for (unsigned int i = 0; i < rna->_protein_list.size(); i++) {
      //compute the activity
      if (rna->_enhancing_coef_list[i] > 0)
      {
        nb_activators++;
        mean_activator_activity += rna->_enhancing_coef_list[i];
        if (rna->_enhancing_coef_list[i]>max_activator_activity) max_activator_activity = rna->_enhancing_coef_list[i];
        if (rna->_enhancing_coef_list[i]<min_activator_activity) min_activator_activity = rna->_enhancing_coef_list[i];
      }

      if (rna->_operating_coef_list[i] > 0)
      {
        nb_operators++;
        mean_operator_activity += rna->_operating_coef_list[i];
        if (rna->_operating_coef_list[i]>max_operator_activity) max_operator_activity = rna->_operating_coef_list[i];
        if (rna->_operating_coef_list[i]<min_operator_activity) min_operator_activity = rna->_operating_coef_list[i];
      }
      double merged_activity = rna->_enhancing_coef_list[i] - rna->_operating_coef_list[i];
      if (merged_activity > 0)
      {
        if (merged_activity>max_merged_activator_activity) max_merged_activator_activity = merged_activity;
        if (merged_activity<min_activator_activity) min_merged_activator_activity = merged_activity;
      }
      if (merged_activity < 0)
      {
        if (merged_activity<max_merged_operator_activity) max_merged_operator_activity = merged_activity;
        if (merged_activity>min_merged_operator_activity) min_merged_operator_activity = merged_activity;
      }
    }
  }

  sprintf( display_string, "Activation links: %" PRId32 ", mean: %lf (max: %lf min:%lf)", nb_activators,mean_activator_activity/(double)nb_activators,max_activator_activity,min_activator_activity);
  win->draw_string( 15, 15, display_string );
  sprintf( display_string, "Inhibition links: %" PRId32 ", mean: %lf (max: %lf min: %lf)", nb_operators,mean_operator_activity/(double)nb_operators,max_operator_activity,min_operator_activity);
  win->draw_string( 15, 30, display_string );

  // end search for max

  // begin links drawing procedure
  for (const auto& rna: _rna_list_coding) {

    // Alpha : angles from OriC (in degrees)
    // Theta : angles on the trigonometric circle (in degrees)
    // nb_sect : "length" in degrees of the arc to be drawn
    alpha_rna_first   = (int16_t) round(  360 * ((double)rna->get_first_transcribed_pos() / (double)genome_length ));
    theta_rna_first   = std::fmod( 90 - alpha_rna_first, 360 );

    // These are the same as above but with a higher precision (1/64 degrees)
    alpha_rna_first_64   = (int16_t) round(64 * 360 * ((double)rna->get_first_transcribed_pos() / (double)genome_length ));
    theta_rna_first_64   = std::fmod( 64 * 90 - alpha_rna_first_64, 64 * 360 );

    pos_rna_x = (win->get_width() / 2.0) + (cos((theta_rna_first_64/(64*180.0)*M_PI)) * diam / 2.0);
    pos_rna_y = (win->get_height() / 2.0) - (sin((theta_rna_first_64/(64*180.0)*M_PI)) * diam / 2.0);

    // ---------------
    //  Draw each regulation link
    // ---------------
    nb_signals = 0;
    int i = 0;
    for (const auto& prot: rna->_protein_list) {
      if ( !(prot->is_signal()) )
      {
        alpha_prot_first   = (int16_t) round(  360 * ((double)prot->get_first_translated_pos() / (double)genome_length ));
        theta_prot_first   = std::fmod( 90 - alpha_prot_first, 360 );

        alpha_prot_first_64   = (int16_t) round(64 * 360 * ((double)prot->get_first_translated_pos() / (double)genome_length ));
        theta_prot_first_64   = std::fmod( 64 * 90 - alpha_prot_first_64, 64 * 360 );

        pos_prot_x = (win->get_width() / 2.0) + (cos((theta_prot_first_64/(64*180.0)*M_PI)) * diam / 2.0);
        pos_prot_y = (win->get_height() / 2.0) - (sin((theta_prot_first_64/(64*180.0)*M_PI)) * diam / 2.0);
      }
      else
      {
        nb_signals+=1;
        pos_prot_x = (win->get_width() / 10.0)*nb_signals;
        pos_prot_y = (win->get_height() * 0.9);
      }

      // compute the color of the link
      double merged_influence = rna->_enhancing_coef_list[i]-rna->_operating_coef_list[i];





      if ( merged_influence > 0 )
      {
        if (max_merged_activator_activity > 0) {
          //printf("ONE %lf %lf %d\n", merged_influence, max_merged_activator_activity,(int)((255 * merged_influence) / max_merged_activator_activity));
          //printf("COLOR : #%02x%02x%02x\n", 0,(int)((255 * merged_influence) / max_merged_activator_activity),0);
          //sprintf( color, "#%02x%02x%02x", 0,(int)((255 * merged_influence) / max_merged_activator_activity),0);
        }
      }
      else
      {
        if (max_merged_activator_activity > 0) {
          //printf("TWO %lf %lf %d\n", merged_influence, max_merged_activator_activity,(int)((255 * merged_influence) / max_merged_activator_activity));
          //printf("COLOR : #%02x%02x%02x\n", (int)((255 * merged_influence) / max_merged_operator_activity),0,0);
          //sprintf( color, "#%02x%02x%02x", (int)((255 * merged_influence) / max_merged_operator_activity),0,0);

        }

      }


      if (merged_influence != 0.0) {
        //compute the lenght of the line
        len_link = sqrt (((pos_rna_x - pos_prot_x) * (pos_rna_x - pos_prot_x)) + ((pos_rna_y - pos_prot_y) * (pos_rna_y - pos_prot_y)));

        //draw the link
        win->draw_line(pos_rna_x,pos_rna_y,pos_prot_x,pos_prot_y,color);

        //draw the arrow going to the protein to the rna regulated the arrow is well centered
        win->draw_line(
            (int16_t)(((pos_rna_x + pos_prot_x)/2.0)+(5.0*(pos_rna_x - pos_prot_x)/len_link)),
            (int16_t)(((pos_rna_y + pos_prot_y)/2.0)+(5.0*(pos_rna_y - pos_prot_y)/len_link)),
            (int16_t)(((pos_rna_x + pos_prot_x)/2.0)+(5.0*(pos_rna_y - pos_prot_y)/len_link)-(5.0*(pos_rna_x - pos_prot_x)/len_link)),
            (int16_t)(((pos_rna_y + pos_prot_y)/2.0)-(5.0*(pos_rna_x - pos_prot_x)/len_link)-(5.0*(pos_rna_y - pos_prot_y)/len_link)),
            color);

        win->draw_line(
            (int16_t)(((pos_rna_x + pos_prot_x)/2.0)+(5.0*(pos_rna_x - pos_prot_x)/len_link)),
            (int16_t)(((pos_rna_y + pos_prot_y)/2.0)+(5.0*(pos_rna_y - pos_prot_y)/len_link)),
            (int16_t)(((pos_rna_x + pos_prot_x)/2.0)-(5.0*(pos_rna_y - pos_prot_y)/len_link)-(5.0*(pos_rna_x - pos_prot_x)/len_link)),
            (int16_t)(((pos_rna_y + pos_prot_y)/2.0)+(5.0*(pos_rna_x - pos_prot_x)/len_link)-(5.0*(pos_rna_y - pos_prot_y)/len_link)),
            color);
      }

      i++;
    }
  }

  delete[] color;
}

void Individual_R_X11::display_phenotype( X11Window* win )
{
  double dist_temp = 0;
  char* color = new char[8];
  char* color2 = NULL;
  int nb_eval = 0;
  strcpy( color, "#FFFFFF" );

  //variable qui sert à stocker du texte à afficher
  char display_string[40];
  int16_t nb_prot = 0;
  int16_t life_time = get_exp_m()->get_exp_s()->get_nb_indiv_age();

  //deux pointeurs utilisés pour défiler dans la liste de protéines
  //  ae_list_node* prot_node  = NULL;
  //  ae_protein_R* prot       = NULL;

  // save the initial list of proteins
  //  ae_list* init_prot_list = new ae_list;
  //  init_prot_list->add_list(_protein_list);
  // Add the signals to the list of proteins
  //ae_environment* envir = NULL;
  //std::vector<ae_protein*> init_protein = _protein_list;
  //std::unordered_map<int,std::vector<ae_protein_R*>> _cloned_signals;


  //for(int i = 0; i < env_list.size(); i++)
  //{
    //envir = env_list[i];
    //    _signals = envir->get_signals();
    //initialise all concentrations of signals to 0
    //std::vector<ae_protein_R*> loc_signals;
    //for ( int8_t j = 0; j < envir->get_signals().size(); j++)
    //{
    //  ae_protein_R* cloned_signal = new ae_protein_R(NULL,*envir->get_signals()[j]);
     // cloned_signal->set_concentration(0.0);
     // _protein_list.push_back(cloned_signal);
     // loc_signals.push_back(cloned_signal);
      //      ((ae_protein_R*)  envir->get_signals()[j])->set_concentration(0.);
      //      _protein_list->add(new ae_list_node(_signals[j]));
    //}
    //_cloned_signals.insert(std::make_pair(envir->get_id(),loc_signals));
    //    _protein_list->add_list(_signals);
  //}

  //set the concentrations of proteins to their initial value
  double* concentrations = new double[_protein_list.size()]; // initialise le tableau de concentrations.
  int16_t prot_index = 0;
  for (auto& prot : _protein_list) {
    ((Protein_R*)prot)->reset_concentration();
    //    printf("prot %d : concentration initiale : %e\n",prot_index,prot->get_concentration());
    concentrations[prot_index++] = ((Protein_R*)prot)->get_concentration();
    //    prot_node = prot_node->get_next();
  }

  // compute steps
  double x_step = 0.8 * win->get_width() / (double)(life_time * get_exp_m()->get_exp_s()->get_nb_degradation_step());
  double y_step = 0.7 * win->get_height();

  // Go from an evaluation date to the next
  //int8_t compteur_env   = 0;
  //envir = env_list[compteur_env];
  //int16_t nb_signals = 0;

  //Add the signals protein to the individual
  //  _signals = envir->get_signals();
  //for ( int8_t i = 0; i < _cloned_signals[compteur_env].size(); i++)
  //{
  //  ((ae_protein_R*) _cloned_signals[compteur_env][i])->set_concentration(0.9);
  //}

  std::set<int>* eval = get_exp_m()->get_exp_s()->get_list_eval_step();

  for (int i = 1; i <= get_exp_m()->get_exp_s()->get_nb_indiv_age(); i++) {
    for (int j = 0; j < get_exp_m()->get_exp_s()->get_nb_degradation_step(); j++)
      update_concentrations();

    update_phenotype();

    //  if (!(indiv_age % 2))
    dynamic_cast<ExpManager_X11*>(get_exp_m())->display_3D(win, *_phenotype, WHITE, (life_time * 5) - (5 * i), (life_time * -3) + (3 * i),true);
    dynamic_cast<ExpManager_X11*>(get_exp_m())->display_3D(win, *(habitat().phenotypic_target().fuzzy()), RED, (life_time * 5) - (5 * i), (life_time * -3) + (3 * i), false);
    // if we have to change of environment
    /*if( ae_common::individual_environment_dates->search(indiv_age) != -1)
    {
      //Remove the signals of this environment
      for ( int8_t i = 0; i < _cloned_signals[compteur_env].size(); i++)
      {
        _cloned_signals[compteur_env][i]->set_concentration(0.);
      }
      // Change the environment at is next value
      compteur_env+=1;
      envir = env_list[compteur_env];

      // Add the signals of this new environment
      for ( int8_t i = 0; i < _cloned_signals[compteur_env].size(); i++)
      {
        _cloned_signals[compteur_env][i]->set_concentration(0.9);
      }
    }*/

    // if its an evaluation date
    if (eval->find(i) != eval->end()) {
      //  printf("indiv_age : %d\n",indiv_age);
      dynamic_cast<ExpManager_X11*>(get_exp_m())->display_3D(win, *_phenotype, RED,  (life_time * 5) - (5 * i), (life_time * -3) + (3 * i),false);
      //    compute_distance_to_target( envir );

      // Compute the difference between the (whole) phenotype and the environment
      AbstractFuzzy* delta = FuzzyFactory::fuzzyFactory->create_fuzzy( *_phenotype );
      delta->sub( *(habitat().phenotypic_target().fuzzy()) );
      double dist_to_target = 0;

      if ( habitat().phenotypic_target().nb_segments() == 1 )
      {
        dist_to_target = delta->get_geometric_area();
      } else {
        PhenotypicSegment ** segments = habitat().phenotypic_target().segments();
        for (size_t i = 0 ; i < static_cast<size_t>(habitat().phenotypic_target().nb_segments()) ; i++) {
          if (segments[i]->feature == METABOLISM) {
            dist_to_target += delta->get_geometric_area(
                segments[i]->start, segments[i]->stop);
          }
        }
      }

      delete delta;

      dist_temp += dist_to_target;
      nb_eval++;

      //Draw the intermediate evaluation result
      sprintf( display_string, " t = %" PRId32 ", dist_to_target =  %lf",i,dist_to_target);
      win->draw_string( 15, 15*nb_eval, display_string );
    }
  }

  //Draw the evaluation result
  sprintf( display_string, "Mean dist_to_target =  %lf",dist_temp/(double)nb_eval);
  win->draw_string( 15, 15*(nb_eval + 1), display_string );

  //Remove the signals of the last environment
  /*for ( int8_t i = 0; i < _cloned_signals[compteur_env].size(); i++)
  {
    _cloned_signals[compteur_env][i]->set_concentration(0.0);
  }*/

  //  _protein_list->erase(NO_DELETE);
  //  _protein_list->add_list(init_prot_list);
  //_protein_list.clear();
  //_protein_list = init_protein;
  //init_protein.clear();
  delete[] concentrations;
  //  delete init_prot_list;
  delete[] color;
}


// =================================================================
//                           Protected Methods
// =================================================================
} // namespace aevol