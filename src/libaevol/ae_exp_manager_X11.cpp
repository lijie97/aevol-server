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
#include <stdio.h>
#include <unistd.h>
#include <X11/Xlib.h>
#include <X11/keysym.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_exp_manager_X11.h>
#include <ae_spatial_structure.h>
#include <point.h>
#include <ae_population_X11.h>
#include <ae_individual_X11.h>
#include <ae_X11_window.h>
#ifdef __X11
#include "fuzzy_X11.h"
#endif

namespace aevol {


// XCheckMaskEvent() doesn't get ClientMessage Events so use XCheckIfEvent()
// with this Predicate function as a work-around ( ClientMessage events
// are needed in order to catch "WM_DELETE_WINDOW")
static Bool AlwaysTruePredicate ( Display*, XEvent*, char* ) { return True; }



//##############################################################################
//                                                                             #
//                          Class ae_exp_manager_X11                           #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================
#define NB_WIN INT32_C(7) // Number of windows that can be showed => CDS, RNA, phenotype, ...

// =================================================================
//                             Constructors
// =================================================================
ae_exp_manager_X11::ae_exp_manager_X11( void ) : ae_exp_manager()
{
  // Basic initializations
  _win      = NULL;
  _win_size = NULL;
  _win_pos  = NULL;
  _win_name = NULL;
  
  _display_on             = false;
  _handle_display_on_off  = false;
  
  // Initialize XLib stuff
  _display  = XOpenDisplay( NULL );
  _screen   = XDefaultScreen( _display );
  _atoms    = new Atom[2];
  _atoms[0] = XInternAtom( _display, "WM_DELETE_WINDOW", False );
  _atoms[1] = XInternAtom( _display, "WM_PROTOCOLS", False );
  set_codes();
}

//ae_exp_manager_X11::ae_exp_manager_X11( ae_param_overloader* param_overloader /* = NULL */ ) : ae_exp_manager( param_overloader )
/*{
  initialize();
}*/

//ae_exp_manager_X11::ae_exp_manager_X11( char* backup_file_name, bool to_be_run /* = true */, ae_param_overloader* param_overloader /* = NULL */ ) :
/*                   ae_exp_manager( backup_file_name, to_be_run, param_overloader  )
{
  initialize();
}*/

// =================================================================
//                             Destructors
// =================================================================
ae_exp_manager_X11::~ae_exp_manager_X11( void )
{
  delete [] _key_codes;
  delete [] _atoms;
  
  for ( int8_t i = 0 ; i < NB_WIN ; i++ )
  {
    if ( _win != NULL )
    {
      if ( _win[i] != NULL ) delete _win[i];
    }
    if ( _win_size != NULL )
    {
      if ( _win_size[i] != NULL ) delete [] _win_size[i];
    }
    if ( _win_pos != NULL )
    {
      if ( _win_pos[i] != NULL ) delete [] _win_pos[i];
    }
  }
  if ( _win != NULL ) delete [] _win;
  
  XCloseDisplay( _display );
  
  if ( _win_name != NULL ) delete [] _win_name;
  if ( _win_size != NULL ) delete [] _win_size;
  if ( _win_pos != NULL ) delete [] _win_pos;
  
}

// =================================================================
//                            Public Methods
// =================================================================


bool ae_exp_manager_X11::quit_signal_received( void )
{
  return _quit_signal_received;
}

void ae_exp_manager_X11::display( void )
{
  // ---------------------
  // 1) Handle user events
  // ---------------------
  if ( _display_on ) handle_events();


  // --------------------------------------------------
  // 2) Handle signal that toggle the display on or off
  // --------------------------------------------------
  if ( _handle_display_on_off )
  {
    _handle_display_on_off = false;
    
    if ( _display_on ) // Display was "on", close all windows (after saving their current size and position)
    {
      for ( int8_t num_win = 0 ; num_win < NB_WIN ; num_win++ )
      {
        if ( _win[num_win] != NULL )
        {
          // 1) Save current window position and size
          Window aWindow; // Unused
          int x_return, y_return;
          int dest_x_return, dest_y_return;
          unsigned int border_width_return, depth_return; // Unused
          XGetGeometry( _display, _win[num_win]->get_window(), &aWindow, &x_return, &y_return,
                        &_win_size[num_win][0], &_win_size[num_win][1], &border_width_return, &depth_return );
          XTranslateCoordinates( _display, _win[num_win]->get_window(), DefaultRootWindow(_display), 0, 0, &dest_x_return, &dest_y_return, &aWindow );
          
          _win_pos[num_win][0] = dest_x_return - x_return;
          _win_pos[num_win][1] = dest_y_return - y_return;
          
          // 2) Delete window
          delete _win[num_win];
          _win[num_win] = NULL;
        }
      }

      XFlush( _display );
      delete _win;
      _win = NULL;
      
      _display_on = false;
    }
    else // Display was "off", open windows
    {
      _show_window |= 1;
      _new_show_window = _show_window;
      
      // If it's the first time the display is switched on, initialize it.
      if ( _win == NULL ) initialize();
      
      for ( int8_t i = 0 ; i < NB_WIN ; i++ )
      {
        if ( get_show_window(i) )
        {
          _win[i] = new ae_X11_window( _display, _screen, _atoms, _win_pos[i][0], _win_pos[i][1], _win_size[i][0], _win_size[i][1], _win_name[i] );
        }
      }
      
      _display_on = true;
    }
  }

  
  // ----------
  // 3) Display
  // ----------
  if ( _display_on )
  {
    for ( int8_t i = 0 ; i < NB_WIN ; i++ )
    {
      if ( get_show_window( i ) )
      {
        if ( get_new_show_window( i ) )
        {
          draw_window( i );
        }
        else
        {
          refresh_window( i );
        }
      }
    }
    
    // Refresh all windows
    XFlush( _display );
  }
}

void ae_exp_manager_X11::handle_events( void )
{
  XEvent event;
  int8_t win_number;
  

  // XCheckMaskEvent() doesn't get ClientMessage Events so use XCheckIfEvent()
  // with the custom AlwaysTruePredicate function as a work-around (ClientMessage 
  // events are needed in order to catch "WM_DELETE_WINDOW")
  int iCurrEvent    = 0;
  int iIgnoreNoise  = 0;
  while( XCheckIfEvent( _display, &event, AlwaysTruePredicate, 0) )
  {
    iCurrEvent ++;
    win_number = identify_window( event.xany.window );
    
    if( win_number == -1) continue;
    // We discard this event because it occurred on a destroyed window
    // (e.g. the user pressed F3 and then moved or resized the window,
    // before the window was destroyed)

    switch( event.type )
    {
      case ConfigureNotify :
      {
        _win[win_number]->resize( event.xconfigure.width, event.xconfigure.height );
        //~ _win[win_number]->repos( event.xconfigure.x, event.xconfigure.y );
        
        // Mark window as having to be entirely redrawn
        _new_show_window |= 1 << win_number;
        break;
      }
      case MapNotify :
      {
        draw_window( win_number );
        break;
      }
      case Expose:
      {
        if( iCurrEvent > iIgnoreNoise )
        {
          draw_window( win_number );
          iIgnoreNoise = iCurrEvent + XQLength( _display );
        }
        break;
      }
      case ClientMessage :
      {
        if ( ((Atom) event.xclient.data.l[0]) == _atoms[0] ) // The user closed the window by clicking on the cross
        {
          // 1) Save current window position and size
          Window aWindow; // Unused
          int x_return, y_return;
          int dest_x_return, dest_y_return;
          unsigned int border_width_return, depth_return; // Unused
          XGetGeometry( _display, _win[win_number]->get_window(), &aWindow, &x_return, &y_return,
                        &_win_size[win_number][0], &_win_size[win_number][1], &border_width_return, &depth_return );
          XTranslateCoordinates( _display, _win[win_number]->get_window(), DefaultRootWindow(_display), 0, 0, &dest_x_return, &dest_y_return, &aWindow );
          
          _win_pos[win_number][0] = dest_x_return - x_return;
          _win_pos[win_number][1] = dest_y_return - y_return;
          
          // 2) Delete window and mark as "not to be shown"
          delete _win[win_number];
          _win[win_number] = NULL;
          _show_window &= ~(1 << win_number);
          
          // 3) If it was the main that was closed, turn display off.
          if ( win_number == 0 )
          {
            _handle_display_on_off = true;
          }
        }
        break;
      }
      case KeyPress :
      {
        // Not sure a switch would work on any platform => use ifs instead
        if ( event.xkey.keycode == _key_codes[KEY_F1] ||
             event.xkey.keycode == _key_codes[KEY_F2] ||
             event.xkey.keycode == _key_codes[KEY_F3] ||
             event.xkey.keycode == _key_codes[KEY_F4] ||
             event.xkey.keycode == _key_codes[KEY_F5] ||
             event.xkey.keycode == _key_codes[KEY_F6] )
        {
          int8_t num_win;
          
          // Not sure a switch would work on any platform => use ifs instead
          if ( event.xkey.keycode == _key_codes[KEY_F1] ) num_win = 1;
          else if ( event.xkey.keycode == _key_codes[KEY_F2] ) num_win = 2;
          else if ( event.xkey.keycode == _key_codes[KEY_F3] ) num_win = 3;
          else if ( event.xkey.keycode == _key_codes[KEY_F4] ) num_win = 4;
          else if ( event.xkey.keycode == _key_codes[KEY_F5] ) num_win = 5;
          else  num_win = 6; // case where ( event.xkey.keycode == _key_codes[KEY_F6] ) 
          
          if ( get_show_window( num_win ) )
          {
            // 1) Save current window position and size
            Window aWindow; // Unused
            int x_return, y_return;
            int dest_x_return, dest_y_return;
            unsigned int border_width_return, depth_return; // Unused
            XGetGeometry( _display, _win[num_win]->get_window(), &aWindow, &x_return, &y_return,
                          &_win_size[num_win][0], &_win_size[num_win][1], &border_width_return, &depth_return );
            XTranslateCoordinates( _display, _win[num_win]->get_window(), DefaultRootWindow(_display), 0, 0, &dest_x_return, &dest_y_return, &aWindow );
            
            _win_pos[num_win][0] = dest_x_return - x_return;
            _win_pos[num_win][1] = dest_y_return - y_return;
            
            // 2) Delete window and mark as "not to be shown"
            delete _win[num_win];
            _win[num_win] = NULL;
            _show_window &= ~(1 << num_win);
          }
          else
          {
            _win[num_win] = new ae_X11_window(  _display, _screen, _atoms, _win_pos[num_win][0], _win_pos[num_win][1],
                                                _win_size[num_win][0], _win_size[num_win][1], _win_name[num_win] );
            _new_show_window |= 1 << num_win;
            _show_window |= _new_show_window;
            draw_window( num_win );
          }
        }     
        else if ( event.xkey.keycode == _key_codes[KEY_P] )
        {
          printf(" P A U S E D \n");
          bool pause_key  = false;
          while ( ! pause_key )
          {
            if ( XCheckIfEvent( _display, &event, AlwaysTruePredicate, 0) )
            {
              if ( event.xkey.keycode == _key_codes[KEY_P] )
              {
                pause_key = true;               
              }
            }
          }
        }
        else if ( (event.xkey.state & ControlMask) && (event.xkey.keycode == _key_codes[KEY_Q]))
        {
          printf(" Q U I T   R E Q U E S T E D\n");
          _quit_signal_received = true;
          // We do not exit here, because it is cleaner to let the main program
          // interpret the signal and call the destructor of the simulation.
          // This ensures that the connection to the X server will be closed 
          // in a clean way, which is much better if we want to run other
          // X applications afterwards!
        }
        
        break;
      }
      default :
      {
        // printf("event not handled \n");
        break;
      }
    }
  }
}

void ae_exp_manager_X11::toggle_display_on_off( void )
{
  // Mark action to be done
  _handle_display_on_off = true;
}



// =================================================================
//                           Protected Methods
// =================================================================
void ae_exp_manager_X11::initialize( bool with_grid /*= false*/, bool with_plasmids /*= false*/ )
{
  // Initialize window structures
  _win      = new ae_X11_window* [NB_WIN];
  _win_size = new unsigned int* [NB_WIN];
  _win_pos  = new int* [NB_WIN];
  
  for ( int8_t i = 0 ; i < NB_WIN ; i++ )
  {
    _win[i] = NULL;
    
    // Default values
    _win_size[i] = new unsigned int[2];
    _win_size[i][0] = 300;
    _win_size[i][1] = 300;
    _win_pos[i] = new int[2];
    _win_pos[i][0]  = 0;
    _win_pos[i][1]  = 0;
  }
  
  // Set phenotype window width
  _win_size[1][0] = 600;
  
  // Set CDS and RNA window width
  if ( with_plasmids )
  {
    _win_size[2][0] = 600;
    _win_size[3][0] = 600;
  }
  
  // Set initial positions if screen is large enough
  if ( with_plasmids && with_grid )
  {
    //if ( XDisplayWidth( _display, _screen ) >= 900 && XDisplayHeight( _display, _screen ) >= 650 )
    {
      _win_pos[0][0]  = 0;
      _win_pos[0][1]  = 0;
      _win_pos[1][0]  = 300;
      _win_pos[1][1]  = 0;
      _win_pos[2][0]  = 0;
      _win_pos[2][1]  = 350;
      _win_pos[3][0]  = 0;
      _win_pos[3][1]  = 700;
    }
  }
  else if ( with_plasmids )
  {
    //if ( XDisplayWidth( _display, _screen ) >= 900 && XDisplayHeight( _display, _screen ) >= 650 )
    {
      _win_pos[0][0]  = 0;
      _win_pos[0][1]  = 0;
      _win_pos[1][0]  = 300;
      _win_pos[1][1]  = 0;
      _win_pos[2][0]  = 0;
      _win_pos[2][1]  = 350;
      _win_pos[3][0]  = 0;
      _win_pos[3][1]  = 700;
    }
  }
  else if ( with_grid )
  {
    //if ( XDisplayWidth( _display, _screen ) >= 900 && XDisplayHeight( _display, _screen ) >= 650 )
    {
      _win_pos[0][0]  = 0;
      _win_pos[0][1]  = 0;
      _win_pos[1][0]  = 300;
      _win_pos[1][1]  = 0;
      _win_pos[2][0]  = 0;
      _win_pos[2][1]  = 350;
      _win_pos[3][0]  = 300;
      _win_pos[3][1]  = 350;
    }
  }
  else // ( ! with_plasmids && ! with_grid )
  {
    //if ( XDisplayWidth( _display, _screen ) >= 900 && XDisplayHeight( _display, _screen ) >= 650 )
    {
      _win_pos[0][0]  = 0;
      _win_pos[0][1]  = 0;
      _win_pos[1][0]  = 300;
      _win_pos[1][1]  = 0;
      _win_pos[2][0]  = 0;
      _win_pos[2][1]  = 350;
      _win_pos[3][0]  = 300;
      _win_pos[3][1]  = 350;
    }
  }
  
  
  // Visible windows at the beginning of the run
  if ( with_grid )
  {
    _show_window  = 0x007F; // hex for bin 1111111  => show first 7 windows
  }
  else
  {
    _show_window  = 0x000F;   // hex for bin 1111   => show first 4 windows
  }
  _new_show_window = _show_window;
  
  
  _win_name = new char*[NB_WIN];
  if ( with_grid )  _win_name[0] = (char*) "Population grid";
  else              _win_name[0] = (char*) "Population";
  
  _win_name[1] = (char*) "Phenotypic profile";
  _win_name[2] = (char*) "Genes";
  _win_name[3] = (char*) "RNAs";
  _win_name[4] = (char*) "Secreted compound present";
  _win_name[5] = (char*) "Metabolic fitness";
  _win_name[6] = (char*) "Current secretion";
}

int8_t ae_exp_manager_X11::identify_window( Window winID )
{
  for ( int8_t i = 0 ; i < NB_WIN ; i++ )
  {
    if ( _win[i] != NULL )
    {
      if ( _win[i]->get_window() == winID ) return i;
    }
  }
  
  return -1;
}

void ae_exp_manager_X11::draw_window( int8_t win_number )
{
  if ( _win[win_number] == NULL)
  {
    fprintf(stderr, "Error: cannot draw this window, it doesn't exist.\n");
    return;
  }
  
  ae_X11_window* cur_win = _win[win_number];

  switch ( win_number )
  {
    case 0:
    {
      break;
    }

    case 1:
    {
      cur_win->blacken();
      
      // Display colour bar
      char* color;
      for ( int16_t i = 0 ; i < cur_win->get_width() ; i++ )
      {
        color = ae_X11_window::get_color( ((double)i / cur_win->get_width()) * (X_MAX - X_MIN) );
        //~ cur_win->draw_line( i, 0, i, cur_win->get_height() / 20, color );
        cur_win->draw_line( i, cur_win->get_height() * 19 / 20, i, cur_win->get_height(), color );
        delete [] color;
      }
      
      break;
    }

    case 2:
    {
      break;
    }

    case 3:
    {
      break;
    }

    // Amount of secreted compound present at each grid location
    case 4:
    {
      break;
    }
   
    // Metabolic fitness  grid
    case 5:
    {
      break;
    }

    // Current secretion (how much is secreted by each organism)
    case 6:
    {
      break;
    }
  }
  
  refresh_window( win_number );
  _new_show_window &= ~(1 << win_number);

  XFlush(_display);
}

void ae_exp_manager_X11::refresh_window( int8_t win_number )
{
  if ( _win[win_number] == NULL)
  {
    fprintf(stderr, "Error: cannot draw this window, it doesn't exist.\n");
    return;
  }
  
  ae_X11_window* cur_win = _win[win_number];

  switch ( win_number )
  {
    // Main window (population)
    case 0 :
    {
      cur_win->blacken();
      
      if ( is_spatially_structured() )
      {
        double** grid = get_spatial_structure()->get_total_fitness_grid();
        ((ae_population_X11*)_pop)->display_grid( cur_win, grid );
        
        // Has been allocated in ae_spatial_structure::get_total_fitness_grid()
        for ( int16_t x = 0 ; x < get_grid_width() ; x++ )
        {
          delete [] grid[x];
        }
        delete [] grid;
      }
      else
      {
        ((ae_population_X11*)_pop)->display( cur_win );
      }
      break;
    }
    
    // Display phenotypes and environment
    case 1 :
    {
      // Blacken all the window except the colour bar
      cur_win->fill_rectangle( 0, 0, cur_win->get_width(), cur_win->get_height() * 19 / 20, BLACK );
      
      // Mark all the non-metabolic segments (paint them in grey)
      if ( get_env()->get_nb_segments() > 1 )
      {
        ae_env_segment** segments = get_env()->get_segments();
     
        for ( size_t i = 0 ; i < get_env()->get_nb_segments() ; i++ )
        {
          if ( segments[i]->feature != METABOLISM )
          {
            if ( segments[i]->feature == NEUTRAL )
            {
              cur_win->fill_rectangle(  cur_win->get_width() * segments[i]->start / (X_MAX-X_MIN), 0.0,
                                        cur_win->get_width() * (segments[i]->stop - segments[i]->start) / (X_MAX-X_MIN),
                                        cur_win->get_height() * 19 / 20, DARKER_GREY );
            }
            else
            {
              cur_win->fill_rectangle(  cur_win->get_width() * segments[i]->start / (X_MAX-X_MIN), 0.0,
                                        cur_win->get_width() * (segments[i]->stop - segments[i]->start) / (X_MAX-X_MIN),
                                        cur_win->get_height() * 19 / 20, GREY );
            }
          }
        }
      }    

      // Display all the phenotypes (blue)
      ae_list_node<ae_individual*>* indiv_node = _pop->get_indivs()->get_first();
      ae_individual*  indiv;
      
      while ( indiv_node != NULL )
      {
        indiv = indiv_node->get_obj();
        
        indiv->get_phenotype()->display( cur_win, BLUE );
        
        if ( indiv->get_allow_plasmids() )
        {
          ((Fuzzy_X11*)indiv->get_genetic_unit( 0 )->get_phenotypic_contribution())->display( cur_win, YELLOW );
          ((Fuzzy_X11*)indiv->get_genetic_unit( 1 )->get_phenotypic_contribution())->display( cur_win, GREEN );
        }
        
        indiv_node = indiv_node->get_next();
      }
        
      // Display best indiv's phenotype (white)
      get_best_indiv()->get_phenotype()->display( cur_win, WHITE, true );
      
      // Display environment (red)
      get_env()->display( cur_win, RED, false, true );
    }
    break;

    // Display genes
    case 2 :
    {
      cur_win->blacken();
      
      ae_individual_X11* indiv = dynamic_cast<ae_individual_X11*>( get_best_indiv() );
      indiv->display_cdss( cur_win );
    }
    break;

    // Display RNAs
    case 3 :
    {
      cur_win->blacken();
      
      ae_individual_X11* indiv = dynamic_cast<ae_individual_X11*>( get_best_indiv() );
      indiv->display_rnas( cur_win );
    }
    break;

    // Display the amount of secreted compound present at each location
    case 4 :
    {
      cur_win->blacken();
      
      if ( is_spatially_structured() )
      {
        ((ae_population_X11*)_pop)->display_grid( cur_win, get_spatial_structure()->get_secretion_present_grid());
      }
    }
    break;
   
    // Display the metabolic fitness grid
    case 5 :
    {
      cur_win->blacken();
      
      ((ae_population_X11*)_pop)->display_grid( cur_win, get_spatial_structure()->get_metabolic_fitness_grid());
    }
    break;

    // display current secretion (how much is secreted by each organism)
    case 6:
    {
      cur_win->blacken();
      
      ((ae_population_X11*)_pop)->display_grid( cur_win, get_spatial_structure()->get_secreted_amount_grid());
    }
    break;
  }

  XFlush( _display );
}





void ae_exp_manager_X11::set_codes( void )
{
  _key_codes = new KeyCode[50];
  assert( _key_codes );
  
  _key_codes[KEY_ESCAPE]  = XKeysymToKeycode( _display, XK_Escape );
  _key_codes[KEY_F1]      = XKeysymToKeycode( _display, XK_F1 );
  _key_codes[KEY_F2]      = XKeysymToKeycode( _display, XK_F2 );
  _key_codes[KEY_F3]      = XKeysymToKeycode( _display, XK_F3 );
  _key_codes[KEY_F4]      = XKeysymToKeycode( _display, XK_F4 );
  _key_codes[KEY_F5]      = XKeysymToKeycode( _display, XK_F5 );
  _key_codes[KEY_F6]      = XKeysymToKeycode( _display, XK_F6 );
  _key_codes[KEY_F7]      = XKeysymToKeycode( _display, XK_F7 );
  _key_codes[KEY_F8]      = XKeysymToKeycode( _display, XK_F8 );
  _key_codes[KEY_F9]      = XKeysymToKeycode( _display, XK_F9 );
  _key_codes[KEY_F10]     = XKeysymToKeycode( _display, XK_F10 );
  _key_codes[KEY_F11]     = XKeysymToKeycode( _display, XK_F11 );
  _key_codes[KEY_F12]     = XKeysymToKeycode( _display, XK_F12 );
  _key_codes[KEY_A]       = XKeysymToKeycode( _display, XK_A );
  _key_codes[KEY_Q]       = XKeysymToKeycode( _display, XK_Q );
  _key_codes[KEY_W]       = XKeysymToKeycode( _display, XK_W );
  _key_codes[KEY_Z]       = XKeysymToKeycode( _display, XK_Z );
  _key_codes[KEY_S]       = XKeysymToKeycode( _display, XK_S );
  _key_codes[KEY_X]       = XKeysymToKeycode( _display, XK_X );
  _key_codes[KEY_E]       = XKeysymToKeycode( _display, XK_E );
  _key_codes[KEY_D]       = XKeysymToKeycode( _display, XK_D );
  _key_codes[KEY_C]       = XKeysymToKeycode( _display, XK_C );
  _key_codes[KEY_R]       = XKeysymToKeycode( _display, XK_R );
  _key_codes[KEY_F]       = XKeysymToKeycode( _display, XK_F );
  _key_codes[KEY_V]       = XKeysymToKeycode( _display, XK_V );
  _key_codes[KEY_T]       = XKeysymToKeycode( _display, XK_T );
  _key_codes[KEY_G]       = XKeysymToKeycode( _display, XK_G );
  _key_codes[KEY_B]       = XKeysymToKeycode( _display, XK_B );
  _key_codes[KEY_Y]       = XKeysymToKeycode( _display, XK_Y );
  _key_codes[KEY_H]       = XKeysymToKeycode( _display, XK_H );
  _key_codes[KEY_N]       = XKeysymToKeycode( _display, XK_N );
  _key_codes[KEY_U]       = XKeysymToKeycode( _display, XK_U );
  _key_codes[KEY_J]       = XKeysymToKeycode( _display, XK_J );
  _key_codes[KEY_I]       = XKeysymToKeycode( _display, XK_I );
  _key_codes[KEY_K]       = XKeysymToKeycode( _display, XK_K );
  _key_codes[KEY_O]       = XKeysymToKeycode( _display, XK_O );
  _key_codes[KEY_L]       = XKeysymToKeycode( _display, XK_L );
  _key_codes[KEY_P]       = XKeysymToKeycode( _display, XK_P );
  _key_codes[KEY_M]       = XKeysymToKeycode( _display, XK_M );
  _key_codes[KEY_1]       = XKeysymToKeycode( _display, XK_1 );
  _key_codes[KEY_2]       = XKeysymToKeycode( _display, XK_2 );
  _key_codes[KEY_3]       = XKeysymToKeycode( _display, XK_3 );
  _key_codes[KEY_4]       = XKeysymToKeycode( _display, XK_4 );
  _key_codes[KEY_5]       = XKeysymToKeycode( _display, XK_5 );
  _key_codes[KEY_6]       = XKeysymToKeycode( _display, XK_6 );
  _key_codes[KEY_7]       = XKeysymToKeycode( _display, XK_7 );
  _key_codes[KEY_8]       = XKeysymToKeycode( _display, XK_8 );
  _key_codes[KEY_9]       = XKeysymToKeycode( _display, XK_9 );
}
} // namespace aevol
