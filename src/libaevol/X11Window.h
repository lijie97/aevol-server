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
 
 
#ifndef  AEVOL_X11_WINDOW_H__
#define  AEVOL_X11_WINDOW_H__
 
 
// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>
#include <X11/Xlib.h>



// =================================================================
//                            Project Files
// =================================================================
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>


namespace aevol {


// =================================================================
//                          Class declarations
// =================================================================

enum color_map
{
  WHITE,
  BLACK,
  RED,
  GREEN,
  BLUE,
  ORANGE,
  YELLOW,
  GREY,
  LIGHT_GREY,
  DARK_GREY,
  DARKER_GREY
};





 
class X11Window
{  
  public :
  
    // =================================================================
    //                             Constructors
    // =================================================================
    X11Window( void );
    X11Window(  Display* display, int8_t screen, Atom * atoms,
                    uint16_t pos_x, uint16_t pos_y, 
                    uint16_t width, uint16_t height,
                    const char* caption );
  
    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~X11Window( void );
  
    // =================================================================
    //                              Accessors
    // =================================================================
    GC get_gc_white( void )     { return _gcWhite;      };
    GC get_gc_black( void )     { return _gcBlack;      };
    GC getGCBlue( void )        { return _gcBlue;       };
    GC getGCRed( void )         { return _gcRed;        };
    GC getGCLightGrey( void )   { return _gcLightGrey;  };
    GC getGCDarkGrey( void )    { return _gcDarkGrey;   };
    GC getGCDarkerGrey( void )  { return _gcDarkerGrey; };
    GC getGCGrey( void )        { return _gcGrey;       };
    GC getGCGreen( void )       { return _gcGreen;      };
    GC getGCOrange( void )      { return _gcOrange;     };
    GC getGCYellow( void )      { return _gcYellow;     };

    uint16_t  get_width( void )         { return _width;   };
    uint16_t  get_height( void )        { return _height;  };
    Window        get_window( void )    { return _window;  };
    Display *     get_display( void )   { return _display; };
  
    // =================================================================
    //                            Public Methods
    // =================================================================

    void resize( unsigned int width, unsigned int height );
    inline void blacken( void );
    void draw_string( int16_t x, int16_t y, char* str );
    void draw_line( int16_t x1, int16_t y1, int16_t x2, int16_t y2, color_map color, bool bold = false );
    void draw_line( int16_t x1, int16_t y1, int16_t x2, int16_t y2, char* color, bool bold = false );
    void draw_circle( int16_t x, int16_t y, int16_t diam );
    void draw_arc( int16_t x, int16_t y, int16_t diam, int16_t angle1, int16_t angle2 );
    void draw_arc( int16_t x, int16_t y, int16_t diam, int16_t angle1, int16_t angle2, char* color );
    void draw_arc_64( int16_t x, int16_t y, int16_t diam, int16_t angle1, int16_t angle2 );
    void draw_arc_64( int16_t x, int16_t y, int16_t diam, int16_t angle1, int16_t angle2, char* color );
    void fill_arc( int16_t x, int16_t y, int16_t diam, int16_t angle1, int16_t angle2 );
    void fill_arc( int16_t x, int16_t y, int16_t diam, int16_t angle1, int16_t angle2, char* color );
    void fill_arc_64( int16_t x, int16_t y, int16_t diam, int16_t angle1, int16_t angle2 );
    void fill_arc_64( int16_t x, int16_t y, int16_t diam, int16_t angle1, int16_t angle2, char* color );
    void fill_rectangle( int16_t x, int16_t y, int16_t width, int16_t height, color_map color );
    void fill_rectangle( int16_t x, int16_t y, int16_t width, int16_t height, char* color );
    
    static char* get_color( double mean );
  
    // =================================================================
    //                           Public Attributes
    // =================================================================
  
  
  
  
  
  protected :
  
    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    //~ X11Window( void )
    //~ {
      //~ printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      //~ exit( EXIT_FAILURE );
    //~ };
    X11Window( const X11Window &model )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };

    // =================================================================
    //                           Protected Methods
    // =================================================================
    uint32_t get_pixel( Display *display, int8_t screen, char *color_name, uint32_t default_color );
    void set_codes( void );
  
    // =================================================================
    //                          Protected Attributes
    // =================================================================
    Display * _display;
    int8_t    _screen;
    Window    _window;
    Cursor    _cursor;
    
    uint16_t _width;
    uint16_t _height;

  
    // Graphic Contexts
    GC _gcWhite;
    GC _gcBlack;
    GC _gcRed;
    GC _gcGreen;
    GC _gcBlue;
    GC _gcOrange;
    GC _gcYellow;
    GC _gcGrey;
    GC _gcLightGrey;
    GC _gcDarkGrey;
    GC _gcDarkerGrey;
};


// =====================================================================
//                          Accessors' definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================
void X11Window::blacken( void )
{
  fill_rectangle( 0, 0, _width, _height, BLACK );
}

} // namespace aevol
#endif // AEVOL_X11_WINDOW_H__
