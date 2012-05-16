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
#include <assert.h>
#include <string>
#include <X11/Xutil.h>
#include <X11/cursorfont.h>




// =================================================================
//                            Project Files
// =================================================================
#include <ae_X11_window.h>
#include <ae_common.h>
#include <ae_simulation.h>




// =================================================================
//                       Basic X11/Xlib notions
// =================================================================
//
// THE DISPLAY
//
// The major notion of using Xlib is the X display. This is a structure 
// representing the connection we have open with a given X server. It 
// hides a queue of messages coming from the server, and a queue of 
// pending requests that our client intends to send to the server. 
// In Xlib, this structure is named 'Display'. When we open a connection 
// to an X server, the library returns a pointer to such a structure. 
// Later, we supply this pointer to any Xlib function that should send 
// messages to the X server or receive messages from this server.
// 
//
// THE WINDOWS
//
// X11 relies on a hierarchical model of rectangular areas called "Windows".
//
// 1. Each Window can be included in another Window (its parent) and may include
//    other Windows (its children). Windows sharing the same owner are called 
//    siblings.
// 2. The screen itself is a Window (the Root Window) that contains all Windows.
// 3. A window can be above or behind a sibling Window. The Window which is above
//    hides partly or completely the other one.
// 4. Any drawing made in a Window is automatically "cut", meaning that only the 
//    part of the drawing which is inside the Window is drawn.
// 5. A Window can be hidden or displayed ("mapped"). The drawing instructions
//    made on an unmapped Window are ignored. By default, newly created windows 
//    are not mapped on the screen - they are invisible. In order to make a 
//    window visible, we must use the XMapWindow() function.
// 6. Each event (keyboard, mouse) is aimed at a specific Window.
// 7. A Window does not memorize its content. Each time it must be re-displayed, 
//    it gets an Expose event, and the content must be redrawn as a response to
//    this event. 
//
//
// THE GC (GRAPHICS CONTEXT)
//
// When we perform various drawing operations (graphics, text, etc), we may 
// specify various options for controlling how the data will be drawn - what 
// foreground and background colors to use, how line edges will be connected, 
// what font to use when drawing some text, etc). In order to avoid the need 
// to supply zillions of parameters to each drawing function, a graphical context 
// structure, of type 'GC' is used. We set the various drawing options in this 
// structure, and then pass a pointer to this structure to any drawing routines. 
// This is rather handy, as we often needs to perform several drawing requests 
// with the same options. Thus, we would initialize a graphical context, set the 
// desired options, and pass this GC structure to all drawing functions. 
// Allocating a new GC is done using the XCreateGC() function. 
//    GC XCreateGC( Display *display, Drawable d, uint32_t valuemask, 
//                  XGCValues *values )
// Since a graphics context has zillions of attributes, and since often we want 
// to define only few of them, we need to be able to tell the XCreateGC() which 
// attributes we want to set. This is what the "valuemask" variable is for. 
// We then use the "values" variable to specify actual values for the attributes 
// we defined in the "valuesmask". The rest of the attributes of this GC will 
// be set to their default values. Once we created a graphics context, we can 
// use it in drawing functions. We can also modify its parameters using various 
// functions (e.g. XSetForeground to change the foreground color of the GC).
//
//
// THE EVENTS
//
// A structure of type 'XEvent' is used to pass events received from the X server. 
// Xlib supports a large amount of event types. The XEvent structure contains the 
// type of event received, as well as the data associated with the event (e.g. 
// position on the screen where the event was generated, mouse button associated 
// with the event, region of screen associated with a 'redraw' event, etc). The way 
// to read the event's data depends on the event type. Thus, an XEvent structure 
// contains a C language union of all possible event types (if you're not sure what 
// C unions are, it is time to check your favourite C language manual...). Thus, 
// we could have an XExpose event, an XButton event, an XMotion event, etc.
// After a program creates a window (or several windows), it should tell the X 
// server what types of events it wishes to receive for this window. By default, 
// no events are sent to the program. This is done for optimizing the server-to-client
// connection (i.e. why send a program (that might even be running at the other 
// side of the globe) an event it is not interested in?). It may register for 
// various mouse (also called "pointer") events, keyboard events, expose events, etc.
// In Xlib, we use the XSelectInput() function to register for events. This function 
// accepts 3 parameters - the display structure, an ID of a window, and a mask of 
// the event types it wishes to get. 




//##############################################################################
//                                                                             #
//                             Class ae_X11_window                             #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
ae_X11_window::ae_X11_window( void )
{
}


ae_X11_window::ae_X11_window( Display* display, int8_t screen, Atom* atoms,
                              uint16_t x, uint16_t y, uint16_t width, uint16_t height,
                              const char* caption )
{
  _width    = width;
  _height   = height;
  _display  = display;
  _screen   = screen;


  XSetWindowAttributes win_attributes;
  win_attributes.event_mask = StructureNotifyMask | ExposureMask | KeyPressMask;
  win_attributes.background_pixel = XBlackPixel( _display, _screen );
  
  _window = XCreateWindow(  _display, DefaultRootWindow(_display), x, y, _width, _height, 0,
                            CopyFromParent, CopyFromParent, CopyFromParent,
                            CWBackPixel|CWEventMask, &win_attributes );
  // NB: the 7th parameter is the width of the window's border, it has nothing to do with 
  // the border appended by the window manager, so this is most often set to zero. 

  // Define the title & iconname of the window 
  XSetStandardProperties( _display, _window, caption, caption, None, NULL, 0, NULL );


  // We want to get MapNotify events, KeyPress events...
  XSelectInput( _display, _window, StructureNotifyMask | ExposureMask | KeyPressMask );


  // Create graphical contexts
  uint32_t whiteColor = WhitePixel( _display, _screen );
  XGCValues values;
  values.line_width = 1;

  values.foreground = get_pixel( _display, _screen, (char*)"white", whiteColor );
  values.background = get_pixel( _display, _screen, (char*)"black", whiteColor );
  _gcWhite = XCreateGC( _display, _window, GCForeground|GCBackground|GCLineWidth, &values );


  values.foreground = get_pixel( _display, _screen, (char*)"black",whiteColor );
  values.background = get_pixel( _display, _screen, (char*)"white",whiteColor );
  _gcBlack = XCreateGC( _display, _window, GCForeground|GCBackground|GCLineWidth, &values );

  values.foreground = get_pixel( _display, _screen, (char*)"red",  whiteColor );
  values.background = get_pixel( _display, _screen, (char*)"black",whiteColor );
  _gcRed = XCreateGC( _display, _window, GCForeground|GCBackground|GCLineWidth, &values );

  values.foreground = get_pixel( _display, _screen, (char*)"green",whiteColor );
  values.background = get_pixel( _display, _screen, (char*)"green",whiteColor );
  _gcGreen = XCreateGC( _display, _window, GCForeground|GCBackground|GCLineWidth, &values );

  values.foreground = get_pixel( _display, _screen, (char*)"blue", whiteColor );
  values.background = get_pixel( _display, _screen, (char*)"black",whiteColor );
  _gcBlue = XCreateGC( _display, _window, GCForeground|GCBackground|GCLineWidth, &values );

  values.foreground = get_pixel( _display, _screen, (char*)"orange",whiteColor );
  values.background = get_pixel( _display, _screen, (char*)"orange",whiteColor );
  _gcOrange = XCreateGC( _display, _window, GCForeground|GCBackground|GCLineWidth, &values );

  values.foreground = get_pixel( _display, _screen, (char*)"yellow",whiteColor );
  values.background = get_pixel( _display, _screen, (char*)"yellow",whiteColor );
  _gcYellow = XCreateGC( _display, _window, GCForeground|GCBackground|GCLineWidth, &values );

  values.foreground = get_pixel( _display, _screen, (char*)"lightgrey",whiteColor );
  values.background = get_pixel( _display, _screen, (char*)"lightgrey",whiteColor );
  _gcLightGrey = XCreateGC( _display, _window, GCForeground|GCBackground|GCLineWidth, &values );

  values.foreground = get_pixel( _display, _screen, (char*)"darkgrey", whiteColor );
  values.background = get_pixel( _display, _screen, (char*)"darkgrey", whiteColor );
  _gcDarkGrey = XCreateGC( _display, _window, GCForeground|GCBackground|GCLineWidth, &values );

  values.foreground = get_pixel( _display, _screen, (char*)"grey15", whiteColor );
  values.background = get_pixel( _display, _screen, (char*)"grey15", whiteColor );
  _gcDarkerGrey = XCreateGC( _display, _window, GCForeground|GCBackground|GCLineWidth, &values );

  values.foreground = get_pixel( _display, _screen, (char*)"grey",whiteColor );
  values.background = get_pixel( _display, _screen, (char*)"grey",whiteColor );
  _gcGrey = XCreateGC( _display, _window, GCForeground|GCBackground|GCLineWidth, &values );


  XMapWindow( _display, _window );
  XMoveWindow( _display, _window, x, y );
  XFlush( _display );


  // Necessary to handle window closing
  XSetWMProtocols( _display, _window, atoms, 2 );
}


// =================================================================
//                             Destructors
// =================================================================

ae_X11_window::~ae_X11_window( void )
{
  XFreeGC( _display, _gcWhite );
  XFreeGC( _display, _gcBlack );
  XFreeGC( _display, _gcRed );
  XFreeGC( _display, _gcGreen );
  XFreeGC( _display, _gcBlue );
  XFreeGC( _display, _gcOrange );
  XFreeGC( _display, _gcYellow );
  XFreeGC( _display, _gcGrey );
  XFreeGC( _display, _gcLightGrey );
  XFreeGC( _display, _gcDarkGrey );
  XFreeGC( _display, _gcDarkerGrey );

  XDestroyWindow( _display, _window );
}

// =================================================================
//                            Public Methods
// =================================================================



void ae_X11_window::resize( unsigned int width, unsigned int height )
{
  _width  = width;
  _height = height;
}


void ae_X11_window::draw_string( int16_t x, int16_t y, char * str )
{
  XDrawImageString( _display, _window, _gcWhite, x, y, str, strlen(str) );
}

void ae_X11_window::draw_line( int16_t x1, int16_t y1, int16_t x2, int16_t y2, color_map color, bool bold /*= false*/ )
{
  GC* gc = NULL;
  
  // Determine which GC to use
  switch ( color )
  {
    case WHITE :
      gc = &_gcWhite;
      break;
    case BLACK :
      gc = & _gcBlack;
      break;
    case RED :
      gc = & _gcRed;
      break;
    case GREEN :
      gc = & _gcGreen;
      break;
    case BLUE :
      gc = & _gcBlue;
      break;
    case ORANGE :
      gc = & _gcOrange;
      break;
    case YELLOW :
      gc = & _gcYellow;
      break;
    case GREY :
      gc = & _gcGrey;
      break;
    case LIGHT_GREY :
      gc = & _gcLightGrey;
      break;
    case DARK_GREY :
      gc = & _gcDarkGrey;
      break;
    case DARKER_GREY :
      gc = & _gcDarkerGrey;
      break;
  }
  
  
  // Draw line (lines if bold)
  XDrawLine( _display, _window, *gc, x1, y1, x2, y2 );
  if ( bold )
  {
    XDrawLine( _display, _window, *gc, x1-1, y1, x2-1, y2 );
    XDrawLine( _display, _window, *gc, x1+1, y1, x2+1, y2 );
  }
}

void ae_X11_window::draw_line( int16_t x1, int16_t y1, int16_t x2, int16_t y2, char* color, bool bold /*= false*/ )
{
  // Create custom GC
  XGCValues values;
  values.foreground = get_pixel( _display, _screen, color, WhitePixel(_display,_screen) );
  values.background = get_pixel( _display, _screen, color, WhitePixel(_display,_screen) );
  GC tmp_gc = XCreateGC( _display, _window, GCForeground|GCBackground, &values );

  // Draw line (lines if bold)
  XDrawLine( _display, _window, tmp_gc, x1, y1, x2, y2 );
  if ( bold )
  {
    XDrawLine( _display, _window, tmp_gc, x1-1, y1, x2-1, y2 );
    XDrawLine( _display, _window, tmp_gc, x1+1, y1, x2+1, y2 );
  }

  XFreeGC( _display, tmp_gc );
}

void ae_X11_window::draw_circle( int16_t x, int16_t y, int16_t diam )
{
  XDrawArc( _display, _window, _gcWhite, x, y, diam, diam, 0, 64*360 );
}

void ae_X11_window::draw_arc( int16_t x, int16_t y, int16_t diam, int16_t angle1, int16_t angle2 )
{
  XDrawArc( _display, _window, _gcWhite, x, y, diam, diam, 64*angle1, 64*angle2 );
}

void ae_X11_window::draw_arc( int16_t x, int16_t y, int16_t diam, int16_t angle1, int16_t angle2, char* color )
{
  XGCValues values;
  values.line_width = 2;
  values.foreground = get_pixel( _display, _screen, color, WhitePixel(_display,_screen) );
  values.background = get_pixel( _display, _screen, color, WhitePixel(_display,_screen) );
  GC tmp_gc = XCreateGC( _display, _window, GCForeground|GCBackground|GCLineWidth, &values );

  XDrawArc( _display, _window, tmp_gc, x, y, diam, diam, 64*angle1, 64*angle2 );

  XFreeGC( _display, tmp_gc );
}

void ae_X11_window::draw_arc_64( int16_t x, int16_t y, int16_t diam, int16_t angle1, int16_t angle2 )
{
  XDrawArc( _display, _window, _gcWhite, x, y, diam, diam, angle1, angle2 );
}

void ae_X11_window::draw_arc_64( int16_t x, int16_t y, int16_t diam, int16_t angle1, int16_t angle2, char* color )
{
  XGCValues values;
  values.line_width = 2;
  values.foreground = get_pixel( _display, _screen, color, WhitePixel(_display,_screen) );
  values.background = get_pixel( _display, _screen, color, WhitePixel(_display,_screen) );
  GC tmp_gc = XCreateGC( _display, _window, GCForeground|GCBackground|GCLineWidth, &values );

  XDrawArc( _display, _window, tmp_gc, x, y, diam, diam, angle1, angle2 );

  XFreeGC( _display, tmp_gc );
}

void ae_X11_window::fill_arc( int16_t x, int16_t y, int16_t diam, int16_t angle1, int16_t angle2 )
{
  XFillArc( _display, _window, _gcWhite, x, y, diam, diam, 64*angle1, 64*angle2 );
}

void ae_X11_window::fill_arc( int16_t x, int16_t y, int16_t diam, int16_t angle1, int16_t angle2, char* color )
{
  XGCValues values;
  values.line_width = 2;
  values.foreground = get_pixel( _display, _screen, color, WhitePixel(_display,_screen) );
  values.background = get_pixel( _display, _screen, color, WhitePixel(_display,_screen) );
  GC tmp_gc = XCreateGC( _display, _window, GCForeground|GCBackground|GCLineWidth, &values );

  XFillArc( _display, _window, tmp_gc, x, y, diam, diam, 64*angle1, 64*angle2 );

  XFreeGC( _display, tmp_gc );
}

void ae_X11_window::fill_arc_64( int16_t x, int16_t y, int16_t diam, int16_t angle1, int16_t angle2 )
{
  XFillArc( _display, _window, _gcWhite, x, y, diam, diam, angle1, angle2 );
}

void ae_X11_window::fill_arc_64( int16_t x, int16_t y, int16_t diam, int16_t angle1, int16_t angle2, char* color )
{
  XGCValues values;
  values.line_width = 2;
  values.foreground = get_pixel( _display, _screen, color, WhitePixel(_display,_screen) );
  values.background = get_pixel( _display, _screen, color, WhitePixel(_display,_screen) );
  GC tmp_gc = XCreateGC( _display, _window, GCForeground|GCBackground|GCLineWidth, &values );

  XFillArc( _display, _window, tmp_gc, x, y, diam, diam, angle1, angle2 );

  XFreeGC( _display, tmp_gc );
}

void ae_X11_window::fill_rectangle( int16_t x, int16_t y, int16_t width, int16_t height, color_map color )
{
  switch ( color )
  {
    case WHITE :
      XFillRectangle( _display, _window, _gcWhite, x, y, width, height );
      break;
    case BLACK :
      XFillRectangle( _display, _window, _gcBlack, x, y, width, height );
      break;
    case RED :
      XFillRectangle( _display, _window, _gcRed, x, y, width, height );
      break;
    case GREEN :
      XFillRectangle( _display, _window, _gcGreen, x, y, width, height );
      break;
    case BLUE :
      XFillRectangle( _display, _window, _gcBlue, x, y, width, height );
      break;
    case ORANGE :
      XFillRectangle( _display, _window, _gcOrange, x, y, width, height );
      break;
    case YELLOW :
      XFillRectangle( _display, _window, _gcYellow, x, y, width, height );
      break;
    case GREY :
      XFillRectangle( _display, _window, _gcGrey, x, y, width, height );
      break;
    case LIGHT_GREY :
      XFillRectangle( _display, _window, _gcLightGrey, x, y, width, height );
      break;
    case DARK_GREY :
      XFillRectangle( _display, _window, _gcDarkGrey, x, y, width, height );
      break;
    case DARKER_GREY :
      XFillRectangle( _display, _window, _gcDarkerGrey, x, y, width, height );
      break;
  }
}

void ae_X11_window::fill_rectangle( int16_t x, int16_t y, int16_t width, int16_t height, char* color )
{
  XGCValues values;
  values.foreground = get_pixel( _display, _screen, color, WhitePixel(_display,_screen) );
  values.background = get_pixel( _display, _screen, color, WhitePixel(_display,_screen) );
  GC tmp_gc = XCreateGC( _display, _window, GCForeground|GCBackground, &values );

  XFillRectangle( _display, _window, tmp_gc, x, y, width, height );

  XFreeGC( _display, tmp_gc );
}

char* ae_X11_window::get_color( double mean )
{
  int16_t red, green, blue;

  double  mean_range     = MAX_X - MIN_X;
  double  mean_range_5   = MIN_X + mean_range / 5;
  double  mean_range_2_5 = MIN_X + 2 * mean_range / 5;
  double  mean_range_3_5 = MIN_X + 3 * mean_range / 5;
  double  mean_range_4_5 = MIN_X + 4 * mean_range / 5;

  if ( mean < mean_range_5 )
  {
    red   = 0;
    green = 255 * ( 1.0 - ( (mean_range_5 - mean) / mean_range_5 ) );
    blue  = 255;
  }
  else if ( mean < mean_range_2_5 )
  {
    red   = 0;
    green = 255;
    blue  = 255 * ( (mean_range_2_5 - mean) / mean_range_5 );
  }
  else if ( mean < mean_range_3_5 )
  {
    red   = 255 * ( 1.0 - ( (mean_range_3_5 - mean) / mean_range_5 ) );
    green = 255;
    blue  = 0;
  }
  else if ( mean < mean_range_4_5 )
  {
    red   = 255;
    green = 255 * ( (mean_range_4_5 - mean) / mean_range_5 );
    blue  = 0;
  }
  else
  {
    red   = 255;
    green = 0;
    blue  = 255 * ( 1.0 - ( (mean_range - mean) / mean_range_5 ) );
  }

  char* color = new char[8];
  sprintf( color, "#%02x%02x%02x", red, green, blue );
  return color;
}


// =================================================================
//                           Protected Methods
// =================================================================

uint32_t ae_X11_window::get_pixel( Display *display, int8_t screen, char *color_name, uint32_t default_color )
{
  XColor color;

  if ( XParseColor( display, DefaultColormap(display,screen), color_name, &color ) == 0 )
  {
    fprintf( stderr, "Invalid colour : %s\n", color_name );
    return default_color;
  }

  if ( XAllocColor( display, DefaultColormap(display,screen), &color ) == 0 )
  {
    fprintf( stderr, "Could not allocate colour %s\n", color_name );
    return default_color;
  }


  return color.pixel;
}




