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



// =================================================================
//                            Project Files
// =================================================================
#include <ae_population_X11.h>
#include <ae_individual_X11.h>
#include <ae_simulation_X11.h>




//##############################################################################
//                                                                             #
//                           Class ae_population_X11                           #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
ae_population_X11::ae_population_X11( void ) : ae_population()
{
  compute_colormap();
}

ae_population_X11::ae_population_X11( gzFile backup_file ) : ae_population( backup_file )
{
  compute_colormap();
}

// =================================================================
//                             Destructors
// =================================================================
ae_population_X11::~ae_population_X11( void )
{
  delete [] _col_map;
}

// =================================================================
//                            Public Methods
// =================================================================

// Display unstructured population
void ae_population_X11::display( ae_X11_window* win )
{
  char generation[40];
  sprintf( generation, "Generation = %"PRId32, ae_common::sim->get_num_gener() );
  win->draw_string( 15, 15, generation );
}

// Display a grid of values
void ae_population_X11::display_grid( ae_X11_window* win, double** cell_grid )
{
  assert( ae_common::pop_structure );
  
  // printf("display grid\n");
  char gener[40];
  int num_colors = 50; 
  
  sprintf( gener, "Generation = %"PRId32, ae_common::sim->get_num_gener() );
  win->draw_string( 15, 15, gener );
  
  
  const int grid_x = ae_common::grid_x;
  const int grid_y = ae_common::grid_y;

  int nb_slots_in_a_row = (int) grid_y;
  int slot_width = 200/nb_slots_in_a_row;
  int x1 = 50 + 50 + slot_width/2;
  int y1 = 75 + 50 + slot_width/2;

  // create the colormap colors to be used for grid plotting
  int cell_size = 5;

  // draw the color scale for fitness
  int y_step_size = grid_y*cell_size/num_colors;
  for ( int i = 0; i  < num_colors; i++ )
  {
    win->fill_rectangle( x1 - 30, y1 - 80 + y_step_size * i,
                         cell_size * 5, y_step_size,
                         _col_map[num_colors-1-i] );
  }

  // find min/max of the matrix
  double grid_max = 0;
  double grid_min = 1000000;
  for (int x = 0; x < grid_x; x++) {
  for (int y = 0; y < grid_y; y++) {
     if (cell_grid[x][y] > grid_max) {grid_max = cell_grid[x][y];}
     if (cell_grid[x][y] < grid_min) {grid_min = cell_grid[x][y];}
   }
  }
  double col_sec_interval = (grid_max - grid_min)/49;

  char scale_txt[40];
  sprintf(scale_txt,"%f", grid_max);
  win->draw_string(x1-80, y1-80,scale_txt);
  sprintf(scale_txt,"%f", grid_min);
  win->draw_string(x1-80, y1-80+grid_y*cell_size,scale_txt);

  for (int x = 0; x < grid_x; x++)
  {
    for (int y = 0; y < grid_y; y++)
    {
      char * col_string;
      // calculate the color
      int new_col;
      if (col_sec_interval==0)
      {
        new_col = 0;
      }
      else
      {
        new_col = (int) floor((cell_grid[x][y] - grid_min) / col_sec_interval);
      }
      col_string = _col_map[new_col];

      // draw a colored rectangle for each cell
      win->fill_rectangle( x1 + 50 + x*cell_size, y1 - 80 + y*cell_size, cell_size, cell_size, col_string );
    }
  }
}


void ae_population_X11::compute_colormap( void )
{
    _col_map = new char* [50];
    
    _col_map[0] = (char*)"RGBi:1.0/0.0/0.0";
    _col_map[1] = (char*)"RGBi:1.0/0.1/0.0";   
    _col_map[2] = (char*)"RGBi:1.0/0.2/0.0";
    _col_map[3] = (char*)"RGBi:1.0/0.3/0.0";
    _col_map[4] = (char*)"RGBi:1.0/0.4/0.0";
    _col_map[5] = (char*)"RGBi:1.0/0.5/0.0";
    _col_map[6] = (char*)"RGBi:1.0/0.6/0.0";
    _col_map[7] = (char*)"RGBi:1.0/0.7/0.0";
    _col_map[8] = (char*)"RGBi:1.0/0.8/0.0";
    _col_map[9] = (char*)"RGBi:1.0/0.9/0.0";

    _col_map[10] = (char*)"RGBi:0.9/1.0/0.0";
    _col_map[11] = (char*)"RGBi:0.8/1.0/0.0";
    _col_map[12] = (char*)"RGBi:0.7/1.0/0.0";
    _col_map[13] = (char*)"RGBi:0.6/1.0/0.0";
    _col_map[14] = (char*)"RGBi:0.5/1.0/0.0";
    _col_map[15] = (char*)"RGBi:0.4/1.0/0.0";
    _col_map[16] = (char*)"RGBi:0.3/1.0/0.0";
    _col_map[17] = (char*)"RGBi:0.2/1.0/0.0";
    _col_map[18] = (char*)"RGBi:0.1/1.0/0.0";
    _col_map[19] = (char*)"RGBi:0.0/1.0/0.0";

    _col_map[20] = (char*)"RGBi:0.0/1.0/0.1";
    _col_map[21] = (char*)"RGBi:0.0/1.0/0.2";
    _col_map[22] = (char*)"RGBi:0.0/1.0/0.3";
    _col_map[23] = (char*)"RGBi:0.0/1.0/0.4";
    _col_map[24] = (char*)"RGBi:0.0/1.0/0.5";
    _col_map[25] = (char*)"RGBi:0.0/1.0/0.6";
    _col_map[26] = (char*)"RGBi:0.0/1.0/0.7";
    _col_map[27] = (char*)"RGBi:0.0/1.0/0.8";
    _col_map[28] = (char*)"RGBi:0.0/1.0/0.9";
    _col_map[29] = (char*)"RGBi:0.0/1.0/1.0";

    _col_map[30] = (char*)"RGBi:0.0/0.9/1.0";
    _col_map[31] = (char*)"RGBi:0.0/0.8/1.0";
    _col_map[32] = (char*)"RGBi:0.0/0.7/1.0";
    _col_map[33] = (char*)"RGBi:0.0/0.6/1.0";
    _col_map[34] = (char*)"RGBi:0.0/0.5/1.0";
    _col_map[35] = (char*)"RGBi:0.0/0.4/1.0";
    _col_map[36] = (char*)"RGBi:0.0/0.3/1.0";
    _col_map[37] = (char*)"RGBi:0.0/0.2/1.0";
    _col_map[38] = (char*)"RGBi:0.0/0.1/1.0";
    _col_map[39] = (char*)"RGBi:0.0/0.0/1.0";

    _col_map[40] = (char*)"RGBi:0.1/0.0/1.0";
    _col_map[41] = (char*)"RGBi:0.2/0.0/1.0";
    _col_map[42] = (char*)"RGBi:0.3/0.0/1.0";
    _col_map[43] = (char*)"RGBi:0.4/0.0/1.0";
    _col_map[44] = (char*)"RGBi:0.5/0.0/1.0";
    _col_map[45] = (char*)"RGBi:0.6/0.0/1.0";
    _col_map[46] = (char*)"RGBi:0.7/0.0/1.0";
    _col_map[47] = (char*)"RGBi:0.8/0.0/1.0";
    _col_map[48] = (char*)"RGBi:0.9/0.0/1.0";
    _col_map[49] = (char*)"RGBi:1.0/0.0/1.0";
}


// =================================================================
//                           Protected Methods
// =================================================================
