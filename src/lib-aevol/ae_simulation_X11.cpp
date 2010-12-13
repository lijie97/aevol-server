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
  _display                = false;
  _handle_display_on_off  = false;
  
  _graphical_device = XOpenDisplay( NULL );
  _screen           = DefaultScreen( _graphical_device );
  _atoms            = new Atom[2];
  _atoms[0]         = XInternAtom( _graphical_device, "WM_DELETE_WINDOW", False );
  _atoms[1]         = XInternAtom( _graphical_device, "WM_PROTOCOLS", False );
  
  set_codes();
  
  _win = new ae_X11_window* [NB_WIN];
  for ( int8_t i = 0 ; i < NB_WIN ; i++ )
  {
    _win[i] = NULL;
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
  _display                = false;
  _handle_display_on_off  = false;
  
  _graphical_device = XOpenDisplay( NULL );
  _screen           = DefaultScreen( _graphical_device );
  _atoms            = new Atom[2];
  _atoms[0]         = XInternAtom( _graphical_device, "WM_DELETE_WINDOW", False );
  _atoms[1]         = XInternAtom( _graphical_device, "WM_PROTOCOLS", False );

  set_codes();

  _win = new ae_X11_window* [NB_WIN];
  for ( int8_t i = 0 ; i < NB_WIN ; i++ )
  {
    _win[i] = NULL;
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
  
  XCloseDisplay( _graphical_device );
  
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
    
    if ( _display ) // Display was "on", close all windows
    {
      for ( int8_t i = 0 ; i < NB_WIN ; i++ )
      {
        if ( _win[i] != NULL )
        {
          delete _win[i];
          _win[i] = NULL;
        }
      }
      
      _display = false;
    }
    else // Display was "off", open windows
    {
      _show_window |= 1;
      _new_show_window = _show_window;
      
      for ( int8_t i = 0 ; i < NB_WIN ; i++ )
      {
        if ( get_show_window(i) )
        {
          if ( i == 1 ) // Phenotype window
          {
            _win[i] = new ae_X11_window( _graphical_device, _screen, _atoms, 600, 300, _window_name[i] );
          }
          else
          {
            _win[i] = new ae_X11_window( _graphical_device, _screen, _atoms, 300, 300, _window_name[i] );
          }
        }
      }
      
      _display = true;
    }
  }


  // ---------------------
  // 2) Handle user events
  // ---------------------
  handle_events();

  
  // ----------
  // 3) Display
  // ----------
  if ( _display )
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
    XFlush( _graphical_device );
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
  while( XCheckIfEvent( _graphical_device, &event, AlwaysTruePredicate, 0) )
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
        
        // Mark window as having to be entirely redrawn
        _new_show_window |= 1 << win_number;
        break;
      }
      case MapNotify :
      {
        draw_window(win_number);
        break;
      }
      case Expose:
      {
        if( iCurrEvent > iIgnoreNoise )
        {
          draw_window(win_number);
          iIgnoreNoise = iCurrEvent + XQLength( _graphical_device );
        }
        break;
      }
      case ClientMessage :
      {
        if ( (event.xclient.data.l[0] == _atoms[0]) ) // TODO : Warning on this line
        {
          // the user closed the window by clicking on the cross
          delete _win[win_number];
          _win[win_number] = NULL;
          _show_window &= ~(1 << win_number);
          
          if ( _show_window == 0 )
          {
            _display = false;
          }
        }
        break;
      }
      case KeyPress :
      {
        if ( event.xkey.keycode == _key_codes[KEY_F1] )
        {
          if ( get_show_window(1) )
          {
            // Window 1 was shown, we must close it
            delete _win[1];
            _win[1] = NULL;
            _show_window &= ~(1 << 1);
          }
          else
          {
            _win[1] = new ae_X11_window( _graphical_device, _screen, _atoms, 300, 300, _window_name[1] );
            _new_show_window |= 1 << 1;
            _show_window |= _new_show_window;
            draw_window( 1 );
          }
        }
        else if ( event.xkey.keycode == _key_codes[KEY_F2] )
        {
          if ( get_show_window(2) )
          {
            // Window 2 was shown, we must close it
            delete _win[2];
            _win[2] = NULL;
            _show_window &= ~(1 << 2);
          }
          else
          {
            _win[2] = new ae_X11_window( _graphical_device, _screen, _atoms, 300, 300, _window_name[2] );
            _new_show_window |= 1 << 2;
            _show_window |= _new_show_window;
            draw_window( 2 );
          }
        }
        else if ( event.xkey.keycode == _key_codes[KEY_F3] )
        {
          if ( get_show_window(3) )
          {
            // Window 3 was shown, we must close it
            delete _win[3];
            _win[3] = NULL;
            _show_window &= ~(1 << 3);
          }
          else
          {
            _win[3] = new ae_X11_window( _graphical_device, _screen, _atoms, 300, 300, _window_name[3] );
            _new_show_window |= 1 << 3;
            _show_window |= _new_show_window;
            draw_window( 3 );
          }
        }          
        else if ( event.xkey.keycode == _key_codes[KEY_F4] )
        {
          if ( get_show_window(4) )
          {
            // Window 4 was shown, we must close it
            delete _win[4];
            _win[4] = NULL;
            _show_window &= ~(1 << 4);
          }
          else
          {
            _win[4] = new ae_X11_window( _graphical_device, _screen, _atoms, 300, 300, _window_name[4] );
            _new_show_window |= 1 << 4;
            _show_window |= _new_show_window;
            draw_window( 4 );
          }
        }
        else if ( event.xkey.keycode == _key_codes[KEY_F5] )
        {
          if ( get_show_window(5) )
          {
            // Window 5 was shown, we must close it
            delete _win[5];
            _win[5] = NULL;
            _show_window &= ~(1 << 5);
          }
          else
          {
            _win[5] = new ae_X11_window( _graphical_device, _screen, _atoms, 300, 300, _window_name[5] );
            _new_show_window |= 1 << 5;
            _show_window |= _new_show_window;
            draw_window( 5 );
          }
        }
        else if ( event.xkey.keycode == _key_codes[KEY_F6] )
        {
          if ( get_show_window(6) )
          {
            // Window 6 was shown, we must close it
            delete _win[6];
            _win[6] = NULL;
            _show_window &= ~(1 << 6);
          }
          else
          {
            _win[6] = new ae_X11_window( _graphical_device, _screen, _atoms, 300, 300, _window_name[6] );
            _new_show_window |= 1 << 6;
            _show_window |= _new_show_window;
            draw_window( 6 );
          }
        }          
        else if ( event.xkey.keycode == _key_codes[KEY_P] )
        {
          printf(" P A U S E D \n");
          bool pause_key  = false;
          while ( ! pause_key )
          {
            if ( XCheckIfEvent( _graphical_device, &event, AlwaysTruePredicate, 0) )
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

  XFlush(_graphical_device);
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
    case 0:
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
      
      // Display environment (red)
      // ((ae_fuzzy_set_X11*)_env)->display( cur_win, RED ); // TODO : line replaced by next line
      _env->display( cur_win, RED );      

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
      ((ae_fuzzy_set_X11*)_pop->get_best()->get_phenotype())->display( cur_win, WHITE );
      break;
    }

    case 2:
    {
      cur_win->blacken();
      
      dynamic_cast<ae_individual_X11*>(_pop->get_best())->display_cdss( cur_win );
      break;
    }

    case 3:
    {
      cur_win->blacken();
      
      ae_individual* indiv1     = _pop->get_best();
      ae_individual_X11* indiv2 = dynamic_cast<ae_individual_X11*>(indiv1);
      indiv2->display_rnas( cur_win );
      break;
    }

    // display the amount of secreted compound present at each location
    case 4:
    {
      cur_win->blacken();
      
      if ( ae_common::pop_structure )
      {
        ((ae_population_X11*)_pop)->display_grid( cur_win, _pop->get_secretion_present());
      }
      
      break;
    }
   
    // display the metabolic fitness  grid
    case 5:
    {
      cur_win->blacken();
      
      ((ae_population_X11*)_pop)->display_grid( cur_win, _pop->get_fitness_metabolic());
      break;
    }

    // display current secretion (how much is secreted by each organism)
    case 6:
    {
      cur_win->blacken();
      
      ((ae_population_X11*)_pop)->display_grid( cur_win, _pop->get_secreted_amount());
      break;
    }
  }

  XFlush(_graphical_device);
}





void ae_simulation_X11::set_codes( void )
{
  _key_codes = new KeyCode[50];
  assert( _key_codes );
  
  _key_codes[KEY_ESCAPE]  = XKeysymToKeycode( _graphical_device, XK_Escape );
  _key_codes[KEY_F1]      = XKeysymToKeycode( _graphical_device, XK_F1 );
  _key_codes[KEY_F2]      = XKeysymToKeycode( _graphical_device, XK_F2 );
  _key_codes[KEY_F3]      = XKeysymToKeycode( _graphical_device, XK_F3 );
  _key_codes[KEY_F4]      = XKeysymToKeycode( _graphical_device, XK_F4 );
  _key_codes[KEY_F5]      = XKeysymToKeycode( _graphical_device, XK_F5 );
  _key_codes[KEY_F6]      = XKeysymToKeycode( _graphical_device, XK_F6 );
  _key_codes[KEY_F7]      = XKeysymToKeycode( _graphical_device, XK_F7 );
  _key_codes[KEY_F8]      = XKeysymToKeycode( _graphical_device, XK_F8 );
  _key_codes[KEY_F9]      = XKeysymToKeycode( _graphical_device, XK_F9 );
  _key_codes[KEY_F10]     = XKeysymToKeycode( _graphical_device, XK_F10 );
  _key_codes[KEY_F11]     = XKeysymToKeycode( _graphical_device, XK_F11 );
  _key_codes[KEY_F12]     = XKeysymToKeycode( _graphical_device, XK_F12 );
  _key_codes[KEY_A]       = XKeysymToKeycode( _graphical_device, XK_A );
  _key_codes[KEY_Q]       = XKeysymToKeycode( _graphical_device, XK_Q );
  _key_codes[KEY_W]       = XKeysymToKeycode( _graphical_device, XK_W );
  _key_codes[KEY_Z]       = XKeysymToKeycode( _graphical_device, XK_Z );
  _key_codes[KEY_S]       = XKeysymToKeycode( _graphical_device, XK_S );
  _key_codes[KEY_X]       = XKeysymToKeycode( _graphical_device, XK_X );
  _key_codes[KEY_E]       = XKeysymToKeycode( _graphical_device, XK_E );
  _key_codes[KEY_D]       = XKeysymToKeycode( _graphical_device, XK_D );
  _key_codes[KEY_C]       = XKeysymToKeycode( _graphical_device, XK_C );
  _key_codes[KEY_R]       = XKeysymToKeycode( _graphical_device, XK_R );
  _key_codes[KEY_F]       = XKeysymToKeycode( _graphical_device, XK_F );
  _key_codes[KEY_V]       = XKeysymToKeycode( _graphical_device, XK_V );
  _key_codes[KEY_T]       = XKeysymToKeycode( _graphical_device, XK_T );
  _key_codes[KEY_G]       = XKeysymToKeycode( _graphical_device, XK_G );
  _key_codes[KEY_B]       = XKeysymToKeycode( _graphical_device, XK_B );
  _key_codes[KEY_Y]       = XKeysymToKeycode( _graphical_device, XK_Y );
  _key_codes[KEY_H]       = XKeysymToKeycode( _graphical_device, XK_H );
  _key_codes[KEY_N]       = XKeysymToKeycode( _graphical_device, XK_N );
  _key_codes[KEY_U]       = XKeysymToKeycode( _graphical_device, XK_U );
  _key_codes[KEY_J]       = XKeysymToKeycode( _graphical_device, XK_J );
  _key_codes[KEY_I]       = XKeysymToKeycode( _graphical_device, XK_I );
  _key_codes[KEY_K]       = XKeysymToKeycode( _graphical_device, XK_K );
  _key_codes[KEY_O]       = XKeysymToKeycode( _graphical_device, XK_O );
  _key_codes[KEY_L]       = XKeysymToKeycode( _graphical_device, XK_L );
  _key_codes[KEY_P]       = XKeysymToKeycode( _graphical_device, XK_P );
  _key_codes[KEY_M]       = XKeysymToKeycode( _graphical_device, XK_M );
  _key_codes[KEY_1]       = XKeysymToKeycode( _graphical_device, XK_1 );
  _key_codes[KEY_2]       = XKeysymToKeycode( _graphical_device, XK_2 );
  _key_codes[KEY_3]       = XKeysymToKeycode( _graphical_device, XK_3 );
  _key_codes[KEY_4]       = XKeysymToKeycode( _graphical_device, XK_4 );
  _key_codes[KEY_5]       = XKeysymToKeycode( _graphical_device, XK_5 );
  _key_codes[KEY_6]       = XKeysymToKeycode( _graphical_device, XK_6 );
  _key_codes[KEY_7]       = XKeysymToKeycode( _graphical_device, XK_7 );
  _key_codes[KEY_8]       = XKeysymToKeycode( _graphical_device, XK_8 );
  _key_codes[KEY_9]       = XKeysymToKeycode( _graphical_device, XK_9 );
}
