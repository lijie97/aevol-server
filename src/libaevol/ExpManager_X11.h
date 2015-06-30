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
// ****************************************************************************
 
 
#ifndef AEVOL_EXP_SETUP_X11_H__
#define AEVOL_EXP_SETUP_X11_H__
 
 
// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>
#include <X11/Xlib.h>
#include <X11/keysym.h>


// =================================================================
//                            Project Files
// =================================================================
#include "ExpManager.h"
#include "X11Window.h"

namespace aevol {

// =================================================================
//                          Class declarations
// =================================================================
class ExpSetup;
class X11Window;

enum key_map
{
  KEY_ESCAPE = 0,
  KEY_F1 = 1,
  KEY_F2 = 2,
  KEY_F3 = 3,
  KEY_F4 = 4,
  KEY_F5 = 5,
  KEY_F6 = 6,
  KEY_F7 = 7,
  KEY_F8 = 8,
  KEY_F9 = 9,
  KEY_F10 = 10,
  KEY_F11 = 11,
  KEY_F12 = 12,
  KEY_A = 13,
  KEY_Q = 14,
  KEY_W = 15,
  KEY_Z = 16,
  KEY_S = 17,
  KEY_X = 18,
  KEY_E = 19,
  KEY_D = 20,
  KEY_C = 21,
  KEY_R = 22,
  KEY_F = 23,
  KEY_V = 24,
  KEY_T = 25,
  KEY_G = 26,
  KEY_B = 27,
  KEY_Y = 28,
  KEY_H = 29,
  KEY_N = 30,
  KEY_U = 31,
  KEY_J = 32,
  KEY_I = 33,
  KEY_K = 34,
  KEY_O = 35,
  KEY_L = 36,
  KEY_P = 37,
  KEY_M = 38,
  KEY_1 = 41,
  KEY_2 = 42,
  KEY_3 = 43,
  KEY_4 = 44,
  KEY_5 = 45,
  KEY_6 = 46,
  KEY_7 = 47,
  KEY_8 = 48,
  KEY_9 = 49
};




 
class ExpManager_X11 : public ExpManager
{
  friend class ExpSetup;
  
  public :
    
    // =================================================================
    //                             Constructors
    // =================================================================
    ExpManager_X11(void);
  
    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ExpManager_X11(void);
  
    // =================================================================
    //                              Accessors
    // =================================================================
    inline bool             get_display_on(void);
    inline Display *        get_display(void);
    inline int8_t           get_screen(void);
    inline Atom *           get_atoms (void);
    inline bool             get_show_window(int8_t win);
    inline bool             get_new_show_window(int8_t win);
    inline X11Window *  get_window(int8_t win);
  
    // =================================================================
    //                            Public Methods
    // =================================================================
    KeyCode* get_key_codes(void)  { return _key_codes;  };
    virtual void display(void);
    void toggle_display_on_off(void);
    void handle_events(void);
    bool quit_signal_received(void);
    void display(X11Window * win, const Fuzzy& fuzzy, color_map color,
        bool fill = false, bool bold = false);
    void display_grid(X11Window * win, double** cell_grid);

    // =================================================================
    //                           Public Attributes
    // =================================================================
  
  
  
  
  
    protected :
  
    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    //~ ExpManager_X11(void)
    //~ {
      //~ printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      //~ exit( EXIT_FAILURE );
    //~ };
    ExpManager_X11( const ExpManager_X11 &model )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
  
  // =================================================================
  //                           Protected Methods
  // =================================================================
  void initialize(bool with_grid = false, bool with_plasmids = false);
  void compute_colormap(void);  
  void set_codes(void);
  int8_t identify_window(Window winID);
  void draw_window(int8_t win_number);
  void refresh_window(int8_t win_number);

  
  // =================================================================
  //                          Protected Attributes
  // =================================================================
  bool      _display_on;
  bool      _handle_display_on_off;
  uint32_t  _show_window;     // (bitmap) windows that have to be displayed (user switches value pressing F1, F2, ...)
  uint32_t  _new_show_window; // (bitmap) windows that have to be displayed but were not displayed at the last refresh
  Display*  _display;
  int8_t    _screen;
  Atom*     _atoms;
  KeyCode*  _key_codes;
  
  X11Window ** _win;       // Table containing the <nb_windows> windows
  char **         _win_name;  // window names
  unsigned int**  _win_size;  // window sizes
  int**           _win_pos;   // window positions

  std::vector<char*> _col_map;
};


// =====================================================================
//                          Accessors' definitions
// =====================================================================
bool ExpManager_X11::get_display_on(void)
{
  return _display_on;
}

Display *ExpManager_X11::get_display(void)
{
  return _display;
}


bool ExpManager_X11::get_show_window( int8_t win )
{
  return ( (_show_window >> win) & 1 );
}

bool ExpManager_X11::get_new_show_window( int8_t win )
{
  return ( (_new_show_window >> win) & 1 );
}


int8_t ExpManager_X11::get_screen(void)
{
  return _screen;
}


Atom*ExpManager_X11::get_atoms(void)
{
  return _atoms;
}


X11Window *ExpManager_X11::get_window( int8_t win )
{
  return _win[win];
}

// =====================================================================
//                       Inline functions' definition
// =====================================================================

} // namespace aevol

#endif // AEVOL_EXP_SETUP_X11_H__
