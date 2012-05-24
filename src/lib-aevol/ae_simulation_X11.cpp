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
#include <stdio.h>
#include <unistd.h>
#include <X11/Xlib.h>
#include <X11/keysym.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_simulation_X11.h>
#include <ae_population_X11.h>
#include <ae_individual_X11.h>
#include <ae_X11_window.h>
#include <ae_point_2d.h>



// XCheckMaskEvent() doesn't get ClientMessage Events so use XCheckIfEvent()
// with this Predicate function as a work-around ( ClientMessage events
// are needed in order to catch "WM_DELETE_WINDOW")
static Bool AlwaysTruePredicate (Display*, XEvent*, char*) { return True; }



//##############################################################################
//                                                                             #
//                           Class ae_simulation_X11                           #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================
#define NB_WIN INT32_C(7) // Number of windows that can be showed => CDS, RNA, phenotype, ...

// =================================================================
//                             Constructors
// =================================================================
ae_simulation_X11::ae_simulation_X11( ae_param_overloader* param_overloader /* = NULL */ ) : ae_simulation( param_overloader )
{
  _display_on             = false;
  _handle_display_on_off  = false;
  
  // Initialize XLib stuff
  _display  = XOpenDisplay( NULL );
  _screen   = XDefaultScreen( _display );
  _atoms    = new Atom[2];
  _atoms[0] = XInternAtom( _display, "WM_DELETE_WINDOW", False );
  _atoms[1] = XInternAtom( _display, "WM_PROTOCOLS", False );
  set_codes();
  
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
  
  // If screen is large enough, set initial positions
  if ( XDisplayWidth( _display, _screen ) >= 900 && XDisplayHeight( _display, _screen ) >= 650 )
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
  
  
  // Visible windows at the beginning of the run
  if ( ae_common::pop_structure == true )
  {
    _show_window  = 0x007F; // hex for bin 1111111  => show first 7 windows
  }
  else
  {
    _show_window  = 0x000F;   // hex for bin 1111   => show first 4 windows
  }
  _new_show_window = _show_window;
  
  
  _quit_signal = false;
  
  _window_name = new char*[NB_WIN];
  if ( ae_common::pop_structure == true )
  {
    _window_name[0] = (char*) "Population grid";
  }
  else
  {
    _window_name[0] = (char*) "Population";
  }
  _window_name[1] = (char*) "Phenotypic profile";
  _window_name[2] = (char*) "Genes";
  _window_name[3] = (char*) "RNAs";
  _window_name[4] = (char*) "Secreted compound present";
  _window_name[5] = (char*) "Metabolic fitness";
  _window_name[6] = (char*) "Current secretion";
}

ae_simulation_X11::ae_simulation_X11( char* backup_file_name, bool to_be_run /* = true */, ae_param_overloader* param_overloader /* = NULL */ ) : \
                   ae_simulation( backup_file_name, to_be_run, param_overloader  )
{
  _display_on             = false;
  _handle_display_on_off  = false;
  
  _display  = XOpenDisplay( NULL );
  _screen   = DefaultScreen( _display );
  _atoms    = new Atom[2];
  _atoms[0] = XInternAtom( _display, "WM_DELETE_WINDOW", False );
  _atoms[1] = XInternAtom( _display, "WM_PROTOCOLS", False );

  set_codes();
  
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
  
  // If screen is large enough, set initial positions
  if ( XDisplayWidth( _display, _screen ) >= 900 && XDisplayHeight( _display, _screen ) >= 650 )
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
  
  // Visible windows at the beginning of the run
  if ( ae_common::pop_structure == true )
  {
    _show_window  = 0x007F; // hex for bin 1111111  => show first 7 windows
  }
  else
  {
    _show_window  = 0x000F;   // hex for bin 1111   => show first 4 windows
  }
  _new_show_window = _show_window;
  
  
  _quit_signal = false;
  
  _window_name = new char*[NB_WIN];
  if ( ae_common::pop_structure == true )
  {
    _window_name[0] = (char*) "Population grid";
  }
  else
  {
    _window_name[0] = (char*) "Population";
  }
  _window_name[1] = (char*) "Phenotypic profile";
  _window_name[2] = (char*) "Genes";
  _window_name[3] = (char*) "RNAs";
  _window_name[4] = (char*) "Secreted compound present";
  _window_name[5] = (char*) "Metabolic fitness";
  _window_name[6] = (char*) "Current secretion";
}

// =================================================================
//                             Destructors
// =================================================================
ae_simulation_X11::~ae_simulation_X11( void )
{
  delete [] _key_codes;
  delete [] _atoms;
  
  for ( int8_t i = 0 ; i < NB_WIN ; i++ )
  {
    if ( _win[i] != NULL ) delete _win[i];
  }
  delete [] _win;
  
  XCloseDisplay( _display );
  
  delete [] _window_name;
}

// =================================================================
//                            Public Methods
// =================================================================


bool ae_simulation_X11::quit_signal_received( void )
{
  return _quit_signal;
}

void ae_simulation_X11::display( void )
{
  // --------------------------------------------------
  // 1) Handle signal that toggle the display on or off
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
      
      _display_on = false;
    }
    else // Display was "off", open windows
    {
      _show_window |= 1;
      _new_show_window = _show_window;
      
      for ( int8_t i = 0 ; i < NB_WIN ; i++ )
      {
        if ( get_show_window(i) )
        {
          _win[i] = new ae_X11_window( _display, _screen, _atoms, _win_pos[i][0], _win_pos[i][1], _win_size[i][0], _win_size[i][1], _window_name[i] );
        }
      }
      
      _display_on = true;
    }
  }


  // ---------------------
  // 2) Handle user events
  // ---------------------
  handle_events();

  
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

void ae_simulation_X11::handle_events( void )
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
          else if ( event.xkey.keycode == _key_codes[KEY_F6] ) num_win = 6;
          
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
                                                _win_size[num_win][0], _win_size[num_win][1], _window_name[num_win] );
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
          _quit_signal = true;
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

void ae_simulation_X11::toggle_display_on_off( void )
{
  // Mark action to be done
  _handle_display_on_off = true;
}



// =================================================================
//                           Protected Methods
// =================================================================
int8_t ae_simulation_X11::identify_window( Window winID )
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

void ae_simulation_X11::draw_window( int8_t win_number )
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
        color = ae_X11_window::get_color( ((double)i / cur_win->get_width()) * (MAX_X - MIN_X) );
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

void ae_simulation_X11::refresh_window( int8_t win_number )
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
      
      if ( ae_common::pop_structure )
      {
        ((ae_population_X11*)_pop)->display_grid( cur_win, _pop->get_fitness_total() );
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
      if ( ae_common::env_axis_is_segmented )
      {
        ae_env_segment** segments = _env->get_segments();
     
        for ( int16_t i = 0 ; i < ae_common::env_axis_nb_segments ; i++ )
        {
          if ( segments[i]->feature != METABOLISM )
          {
            if ( segments[i]->feature == NEUTRAL )
            {
              cur_win->fill_rectangle(  cur_win->get_width() * segments[i]->start / (MAX_X-MIN_X), 0.0,
                                        cur_win->get_width() * (segments[i]->stop - segments[i]->start) / (MAX_X-MIN_X),
                                        cur_win->get_height() * 19 / 20, DARKER_GREY );
            }
            else
            {
              cur_win->fill_rectangle(  cur_win->get_width() * segments[i]->start / (MAX_X-MIN_X), 0.0,
                                        cur_win->get_width() * (segments[i]->stop - segments[i]->start) / (MAX_X-MIN_X),
                                        cur_win->get_height() * 19 / 20, GREY );
            }
          }
        }
      }    

      // Display all the phenotypes (blue)
      ae_list_node*   indiv_node = _pop->get_indivs()->get_first();
      ae_individual*  indiv;
      
      while ( indiv_node != NULL )
      {
        indiv = (ae_individual*) indiv_node->get_obj();
        
        ((ae_fuzzy_set_X11*)indiv->get_phenotype())->display( cur_win, BLUE );
        
        if ( ae_common::allow_plasmids )
        {
          ((ae_fuzzy_set_X11*)indiv->get_genetic_unit( 0 )->get_phenotypic_contribution())->display( cur_win, YELLOW );
          ((ae_fuzzy_set_X11*)indiv->get_genetic_unit( 1 )->get_phenotypic_contribution())->display( cur_win, GREEN );
        }
        
        indiv_node = indiv_node->get_next();
      }
        
      // Display best indiv's phenotype (white)
      ((ae_fuzzy_set_X11*)_pop->get_best()->get_phenotype())->display( cur_win, WHITE, true );
      
      // Display environment (red)
      // ((ae_fuzzy_set_X11*)_env)->display( cur_win, RED ); // TODO : line replaced by next line
      _env->display( cur_win, RED, false, true );
    }
    break;

    // Display genes
    case 2 :
    {
      cur_win->blacken();
      
      dynamic_cast<ae_individual_X11*>(_pop->get_best())->display_cdss( cur_win );
    }
    break;

    // Display RNAs
    case 3 :
    {
      cur_win->blacken();
      
      ae_individual* indiv1     = _pop->get_best();
      ae_individual_X11* indiv2 = dynamic_cast<ae_individual_X11*>(indiv1);
      indiv2->display_rnas( cur_win );
    }
    break;

    // Display the amount of secreted compound present at each location
    case 4 :
    {
      cur_win->blacken();
      
      if ( ae_common::pop_structure )
      {
        ((ae_population_X11*)_pop)->display_grid( cur_win, _pop->get_secretion_present());
      }
    }
    break;
   
    // Display the metabolic fitness grid
    case 5 :
    {
      cur_win->blacken();
      
      ((ae_population_X11*)_pop)->display_grid( cur_win, _pop->get_fitness_metabolic());
    }
    break;

    // display current secretion (how much is secreted by each organism)
    case 6:
    {
      cur_win->blacken();
      
      ((ae_population_X11*)_pop)->display_grid( cur_win, _pop->get_secreted_amount());
    }
    break;
  }

  XFlush( _display );
}





void ae_simulation_X11::set_codes( void )
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
