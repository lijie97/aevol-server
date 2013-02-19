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
#include <math.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_individual_X11.h>

#include <ae_exp_manager.h>
#include <ae_exp_setup.h>
#include <ae_utils.h>




//##############################################################################
//                                                                             #
//                           Class ae_individual_X11                           #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
ae_individual_X11::ae_individual_X11( ae_exp_manager* exp_m,
                                      ae_jumping_mt* alea,
                                      ae_params_mut* param_mut,
                                      double w_max,
                                      int32_t min_genome_length,
                                      int32_t max_genome_length,
                                      bool allow_plasmids,
                                      int32_t plasmid_minimal_length,
                                      int32_t id,
                                      int32_t age )
        : ae_individual( exp_m, alea, param_mut, w_max, min_genome_length,
                         max_genome_length, allow_plasmids,
                         plasmid_minimal_length, id, age )
{
  init_occupied_sectors();
}

ae_individual_X11::ae_individual_X11( ae_exp_manager* exp_manager, gzFile backup_file )
        : ae_individual( exp_manager, backup_file )
{
  init_occupied_sectors();
}

ae_individual_X11::ae_individual_X11( const ae_individual_X11 &model )
        : ae_individual( model )
{
  init_occupied_sectors();
}

ae_individual_X11::ae_individual_X11( ae_individual_X11* const parent, int32_t id , ae_jumping_mt* prng ) : ae_individual( parent, id, prng)
{
  init_occupied_sectors();
}

/*ae_individual_X11::ae_individual_X11( char* genome, int32_t genome_size ) : ae_individual( genome, genome_size )
{
  init_occupied_sectors();
}*/

// =================================================================
//                             Destructors
// =================================================================
ae_individual_X11::~ae_individual_X11( void )
{
  for ( int16_t layer = 0 ; layer < _outmost_layer ; layer++ )
  {
    delete [] _occupied_sectors[LEADING][layer];
    delete [] _occupied_sectors[LAGGING][layer];
  }
}

// =================================================================
//                            Public Methods
// =================================================================
void ae_individual_X11::display( void )
{
}

void ae_individual_X11::display_cdss( ae_X11_window* win )
{
  // Retreive the genetic unit corresponding to the main chromosome
  ae_genetic_unit* gen_unit = get_genetic_unit(0);
  int32_t genome_length = gen_unit->get_dna()->get_length();
  
  // Display the number of CDSs
  char display_string[40];
  sprintf( display_string, "Main chromosome size : %"PRId32"bp", genome_length );
  win->draw_string( 15, 25, display_string );
  sprintf( display_string, "Leading : %"PRId32" CDSs", gen_unit->get_protein_list()[LEADING]->get_nb_elts() );
  win->draw_string( 15, 40, display_string );
  sprintf( display_string, "Lagging : %"PRId32" CDSs", gen_unit->get_protein_list()[LAGGING]->get_nb_elts() );
  win->draw_string( 15, 55, display_string );

  // Compute display diameter according to genome length and window size
  int16_t canvas_width;
  if ( _allow_plasmids ) canvas_width = win->get_width() / 2;
  else canvas_width = win->get_width();
  int16_t canvas_height = win->get_height();
  
  int16_t canvas_size = ae_utils::min( canvas_width, canvas_width );
  int16_t diam        = round( canvas_size * log( (double)genome_length ) / 16 );

  // Prevent diameter from getting greater than 2/3 of the window size
  if ( diam > 2 * canvas_size / 3 )
  {
    diam = 2 * canvas_size / 3;
  }

  // Compute coordinates of the upper-left corner of the containing square
  int16_t pos_x = (canvas_width - diam) / 2;
  int16_t pos_y = (canvas_height - diam) / 2;

  // Draw main circle
  win->draw_circle( pos_x, pos_y, diam );

  // Sector occupation management
  reset_sectors();


  // ---------------
  //  Draw each CDS
  // ---------------
  ae_list_node<ae_protein*>* cds_node  = NULL;
  ae_protein*   cds       = NULL;

  // NB : As we want OriC to be at the "top" of the circle and the orientation
  //      to be clockwise, the drawing angle (theta) will be given as
  //      (90 - alpha), alpha being the "classical" trigonometric angle
  int16_t alpha_first, alpha_last; // Angles of first and last transcribed bases from OriC (degrees)
  int16_t theta_first, theta_last; // Transposed angles on the trigonometric circle (degrees)
  int16_t nb_sect;
  // Same as above with precision = 1/64 degree
  int16_t alpha_first_64, alpha_last_64;
  int16_t theta_first_64, theta_last_64;
  int16_t nb_sect_64;

  // ----------------
  //  LEADING strand
  // ----------------
  cds_node = gen_unit->get_protein_list()[LEADING]->get_first();

  while ( cds_node != NULL )
  {
    cds = cds_node->get_obj();

    // Alpha : angles from OriC (in degrees)
    // Theta : angles on the trigonometric circle (in degrees)
    // nb_sect : "length" in degrees of the arc to be drawn
    alpha_first   = (int16_t) round(  360 * ((double)cds->get_first_translated_pos() / (double)genome_length ));
    alpha_last    = (int16_t) round(  360 * ((double)cds->get_last_translated_pos()  / (double)genome_length ));
    theta_first   = ae_utils::mod( 90 - alpha_first, 360 );
    theta_last    = ae_utils::mod( 90 - alpha_last, 360 );
    nb_sect       = ae_utils::mod( alpha_last - alpha_first + 1,  360 );

    // These are the same as above but with a higher precision (1/64 degrees)
    alpha_first_64   = (int16_t) round(64 * 360 * ((double)cds->get_first_translated_pos() / (double)genome_length ));
    alpha_last_64    = (int16_t) round(64 * 360 * ((double)cds->get_last_translated_pos() / (double)genome_length ));    
    theta_first_64   = ae_utils::mod( 64 * 90 - alpha_first_64, 64 * 360 );
    theta_last_64    = ae_utils::mod( 64 * 90 - alpha_last_64, 64 * 360 );
    nb_sect_64       = ae_utils::mod( alpha_last_64 - alpha_first_64 + 1,  64 * 360 );


    // Look for the inmost layer that has all the sectors between
    // theta_first and theta_last free
    int16_t layer = 0;
    bool sectors_free = false;
    while ( ! sectors_free )
    {
      sectors_free = true;

      for ( int16_t rho = 0 ; rho < nb_sect ; rho++ )
      {
        if ( _occupied_sectors[LEADING][layer][ae_utils::mod(theta_first-rho, 360)] )
        {
          sectors_free = false;
          break;
        }
      }

      if ( sectors_free )
      {
        break; // All the needed sectors are free on the current layer
      }
      else
      {
        layer++;

        if ( layer >= _outmost_layer )
        {
          add_layer();
          break; // An added layer is necessarily free, no need to look further
        }
      }
    }

    // Mark sectors to be drawn as occupied
    for ( int16_t rho = 0 ; rho < nb_sect ; rho++ )
    {
      _occupied_sectors[LEADING][layer][ae_utils::mod(theta_first-rho, 360)] = true;
    }
    // Mark flanking sectors as occupied
    _occupied_sectors[LEADING][layer][ae_utils::mod(theta_first+1, 360)] = true;
    _occupied_sectors[LEADING][layer][ae_utils::mod(theta_first-nb_sect, 360)] = true;


    // Draw
    const int8_t cds_layer_spacing = 5; // TODO : param?
    layer++; // index starting at 0 but needed to start at 1

    int16_t diam2 = diam + (layer * 2 * cds_layer_spacing);
    pos_x         = (int16_t) round( (double)(canvas_width  - diam2 ) / 2.0 );
    pos_y         = (int16_t) round( (double)(canvas_height - diam2 ) / 2.0 );

    char* color = ae_X11_window::get_color( cds->get_mean() );
    win->draw_arc_64( pos_x, pos_y, diam2, theta_first_64 - nb_sect_64, nb_sect_64, color );
    delete [] color;

    cds_node = cds_node->get_next();
  }

  // ----------------
  //  LAGGING strand
  // ----------------
  cds_node = gen_unit->get_protein_list()[LAGGING]->get_first();

  while ( cds_node != NULL )
  {
    cds = cds_node->get_obj();

    // Alpha : angles from OriC (in degrees)
    // Theta : angles on the trigonometric circle (in degrees)
    // nb_sect : "length" in degrees of the arc to be drawn
    alpha_first   = (int16_t) round(  360 * ((double)cds->get_first_translated_pos() / (double)genome_length ));
    alpha_last    = (int16_t) round(  360 * ((double)cds->get_last_translated_pos()  / (double)genome_length ));
    theta_first   = ae_utils::mod( 90 - alpha_first, 360 );
    theta_last    = ae_utils::mod( 90 - alpha_last, 360 );
    nb_sect = ae_utils::mod( alpha_first - alpha_last + 1,  360 );

    // These are the same as above but with a higher precision (1/64 degrees)
    alpha_first_64   = (int16_t) round(64 * 360 * ((double)cds->get_first_translated_pos() / (double)genome_length ));
    alpha_last_64    = (int16_t) round(64 * 360 * ((double)cds->get_last_translated_pos() / (double)genome_length ));
    theta_first_64   = ae_utils::mod( 64 * 90 - alpha_first_64, 64 * 360 );
    theta_last_64    = ae_utils::mod( 64 * 90 - alpha_last_64, 64 * 360 );
    nb_sect_64 = ae_utils::mod( alpha_first_64 - alpha_last_64 + 1,  64 * 360 );


    // Look for the inmost layer that has all the sectors between
    // theta_first and theta_last free
    int16_t layer = 0;
    bool sectors_free = false;
    while ( ! sectors_free )
    {
      sectors_free = true;

      for ( int16_t rho = 0 ; rho < nb_sect ; rho++ )
      {
        if ( _occupied_sectors[LAGGING][layer][ae_utils::mod(theta_first+rho, 360)] )
        {
          sectors_free = false;
          break;
        }
      }

      if ( sectors_free )
      {
        break; // All the needed sectors are free on the current layer
      }
      else
      {
        layer++;

        if ( layer >= _outmost_layer )
        {
          add_layer();
          break; // An added layer is necessarily free, no need to look further
        }
      }
    }

    // Mark sectors to be drawn as occupied
    for ( int16_t rho = 0 ; rho < nb_sect ; rho++ )
    {
      _occupied_sectors[LAGGING][layer][ae_utils::mod(theta_first+rho, 360)] = true;
    }
    // Mark flanking sectors as occupied
    _occupied_sectors[LAGGING][layer][ae_utils::mod(theta_first-1, 360)] = true;
    _occupied_sectors[LAGGING][layer][ae_utils::mod(theta_first+nb_sect, 360)] = true;


    // Draw
    const int8_t cds_layer_spacing = 5; // TODO : param?
    layer++; // index starting at 0 but needed to start at 1

    int16_t diam2 = diam - (layer * 2 * cds_layer_spacing);
    pos_x         = (int16_t) round( (double)(canvas_width  - diam2 ) / 2.0 );
    pos_y         = (int16_t) round( (double)(canvas_height - diam2 ) / 2.0 );

    char* color = ae_X11_window::get_color( cds->get_mean() );
    win->draw_arc_64( pos_x, pos_y, diam2, theta_first_64, nb_sect_64, color );
    delete [] color;

    cds_node = cds_node->get_next();
  }
  
  
  
  
  
  
  
  
  
  // --------------------------------------------------------------------------This is temporary, it is a big copy-paste of what's above.
  if ( _allow_plasmids )
  {
    // Retreive the genetic unit corresponding to the plasmid
    ae_genetic_unit* gen_unit = get_genetic_unit(1);
    if ( gen_unit == NULL ) return;
    
    int32_t genome_length = gen_unit->get_dna()->get_length();
    

    // Compute display diameter according to genome length and window size
    int16_t canvas_width;
    int16_t canvas_height;
    int16_t canvas_size ;
    if ( _allow_plasmids )
    {
      canvas_width  = win->get_width() / 2;
      canvas_size   = canvas_width;
    }
    else
    {
      canvas_width  = win->get_width();
      canvas_size   = ae_utils::min( canvas_width, canvas_width );
    }
    canvas_height = win->get_height();
    
    int16_t diam  = round( canvas_size * log( (double)genome_length ) / 16 );

    // Prevent diameter from getting greater than 2/3 of the window size
    if ( diam > 2 * canvas_size / 3 )
    {
      diam = 2 * canvas_size / 3;
    }

    // Compute coordinates of the upper-left corner of the containing square
    int16_t pos_x = canvas_width + (canvas_width - diam) / 2;
    int16_t pos_y = (canvas_height - diam) / 2;
    
    
    // Draw main circle
    win->draw_circle( pos_x, pos_y, diam );

    // Sector occupation management
    reset_sectors();


    // ---------------
    //  Draw each CDS
    // ---------------
    ae_list_node<ae_protein*>* cds_node  = NULL;
    ae_protein*   cds       = NULL;

    // NB : As we want OriC to be at the "top" of the circle and the orientation
    //      to be clockwise, the drawing angle (theta) will be given as
    //      (90 - alpha), alpha being the "classical" trigonometric angle
    int16_t alpha_first, alpha_last; // Angles of first and last transcribed bases from OriC (degrees)
    int16_t theta_first, theta_last; // Transposed angles on the trigonometric circle (degrees)
    int16_t nb_sect;
    // Same as above with precision = 1/64 degree
    int16_t alpha_first_64, alpha_last_64;
    int16_t theta_first_64, theta_last_64;
    int16_t nb_sect_64;

    // ----------------
    //  LEADING strand
    // ----------------
    cds_node = gen_unit->get_protein_list()[LEADING]->get_first();

    while ( cds_node != NULL )
    {
      cds = cds_node->get_obj();

      // Alpha : angles from OriC (in degrees)
      // Theta : angles on the trigonometric circle (in degrees)
      // nb_sect : "length" in degrees of the arc to be drawn
      alpha_first   = (int16_t) round(  360 * ((double)cds->get_first_translated_pos() / (double)genome_length ));
      alpha_last    = (int16_t) round(  360 * ((double)cds->get_last_translated_pos()  / (double)genome_length ));
      theta_first   = ae_utils::mod( 90 - alpha_first, 360 );
      theta_last    = ae_utils::mod( 90 - alpha_last, 360 );
      nb_sect       = ae_utils::mod( alpha_last - alpha_first + 1,  360 );

      // These are the same as above but with a higher precision (1/64 degrees)
      alpha_first_64   = (int16_t) round(64 * 360 * ((double)cds->get_first_translated_pos() / (double)genome_length ));
      alpha_last_64    = (int16_t) round(64 * 360 * ((double)cds->get_last_translated_pos() / (double)genome_length ));    
      theta_first_64   = ae_utils::mod( 64 * 90 - alpha_first_64, 64 * 360 );
      theta_last_64    = ae_utils::mod( 64 * 90 - alpha_last_64, 64 * 360 );
      nb_sect_64       = ae_utils::mod( alpha_last_64 - alpha_first_64 + 1,  64 * 360 );


      // Look for the inmost layer that has all the sectors between
      // theta_first and theta_last free
      int16_t layer = 0;
      bool sectors_free = false;
      while ( ! sectors_free )
      {
        sectors_free = true;

        for ( int16_t rho = 0 ; rho < nb_sect ; rho++ )
        {
          if ( _occupied_sectors[LEADING][layer][ae_utils::mod(theta_first-rho, 360)] )
          {
            sectors_free = false;
            break;
          }
        }

        if ( sectors_free )
        {
          break; // All the needed sectors are free on the current layer
        }
        else
        {
          layer++;

          if ( layer >= _outmost_layer )
          {
            add_layer();
            break; // An added layer is necessarily free, no need to look further
          }
        }
      }

      // Mark sectors to be drawn as occupied
      for ( int16_t rho = 0 ; rho < nb_sect ; rho++ )
      {
        _occupied_sectors[LEADING][layer][ae_utils::mod(theta_first-rho, 360)] = true;
      }
      // Mark flanking sectors as occupied
      _occupied_sectors[LEADING][layer][ae_utils::mod(theta_first+1, 360)] = true;
      _occupied_sectors[LEADING][layer][ae_utils::mod(theta_first-nb_sect, 360)] = true;


      // Draw
      const int8_t cds_layer_spacing = 5; // TODO : param?
      layer++; // index starting at 0 but needed to start at 1

      int16_t diam2 = diam + (layer * 2 * cds_layer_spacing);
      pos_x         = canvas_width + (int16_t) round( (double)(canvas_width  - diam2 ) / 2.0 );
      pos_y         = (int16_t) round( (double)(canvas_height - diam2 ) / 2.0 );

      char* color = ae_X11_window::get_color( cds->get_mean() );
      win->draw_arc_64( pos_x, pos_y, diam2, theta_first_64 - nb_sect_64, nb_sect_64, color );
      delete [] color;

      cds_node = cds_node->get_next();
    }

    // ----------------
    //  LAGGING strand
    // ----------------
    cds_node = gen_unit->get_protein_list()[LAGGING]->get_first();

    while ( cds_node != NULL )
    {
      cds = cds_node->get_obj();

      // Alpha : angles from OriC (in degrees)
      // Theta : angles on the trigonometric circle (in degrees)
      // nb_sect : "length" in degrees of the arc to be drawn
      alpha_first   = (int16_t) round(  360 * ((double)cds->get_first_translated_pos() / (double)genome_length ));
      alpha_last    = (int16_t) round(  360 * ((double)cds->get_last_translated_pos()  / (double)genome_length ));
      theta_first   = ae_utils::mod( 90 - alpha_first, 360 );
      theta_last    = ae_utils::mod( 90 - alpha_last, 360 );
      nb_sect = ae_utils::mod( alpha_first - alpha_last + 1,  360 );

      // These are the same as above but with a higher precision (1/64 degrees)
      alpha_first_64   = (int16_t) round(64 * 360 * ((double)cds->get_first_translated_pos() / (double)genome_length ));
      alpha_last_64    = (int16_t) round(64 * 360 * ((double)cds->get_last_translated_pos() / (double)genome_length ));
      theta_first_64   = ae_utils::mod( 64 * 90 - alpha_first_64, 64 * 360 );
      theta_last_64    = ae_utils::mod( 64 * 90 - alpha_last_64, 64 * 360 );
      nb_sect_64 = ae_utils::mod( alpha_first_64 - alpha_last_64 + 1,  64 * 360 );


      // Look for the inmost layer that has all the sectors between
      // theta_first and theta_last free
      int16_t layer = 0;
      bool sectors_free = false;
      while ( ! sectors_free )
      {
        sectors_free = true;

        for ( int16_t rho = 0 ; rho < nb_sect ; rho++ )
        {
          if ( _occupied_sectors[LAGGING][layer][ae_utils::mod(theta_first+rho, 360)] )
          {
            sectors_free = false;
            break;
          }
        }

        if ( sectors_free )
        {
          break; // All the needed sectors are free on the current layer
        }
        else
        {
          layer++;

          if ( layer >= _outmost_layer )
          {
            add_layer();
            break; // An added layer is necessarily free, no need to look further
          }
        }
      }

      // Mark sectors to be drawn as occupied
      for ( int16_t rho = 0 ; rho < nb_sect ; rho++ )
      {
        _occupied_sectors[LAGGING][layer][ae_utils::mod(theta_first+rho, 360)] = true;
      }
      // Mark flanking sectors as occupied
      _occupied_sectors[LAGGING][layer][ae_utils::mod(theta_first-1, 360)] = true;
      _occupied_sectors[LAGGING][layer][ae_utils::mod(theta_first+nb_sect, 360)] = true;


      // Draw
      const int8_t cds_layer_spacing = 5; // TODO : param?
      layer++; // index starting at 0 but needed to start at 1

      int16_t diam2 = diam - (layer * 2 * cds_layer_spacing);
      pos_x         = canvas_width + (int16_t) round( (double)(canvas_width  - diam2 ) / 2.0 );
      pos_y         = (int16_t) round( (double)(canvas_height - diam2 ) / 2.0 );

      char* color = ae_X11_window::get_color( cds->get_mean() );
      win->draw_arc_64( pos_x, pos_y, diam2, theta_first_64, nb_sect_64, color );
      delete [] color;

      cds_node = cds_node->get_next();
    }
  }
}

void ae_individual_X11::display_rnas( ae_X11_window* win )
{
  // Retreive the genetic unit corresponding to the main chromosome
  ae_genetic_unit* gen_unit = _genetic_unit_list->get_first()->get_obj();
  int32_t genome_length = gen_unit->get_dna()->get_length();
  
  // Display the number of RNAs
  char nb_rna[40];
  sprintf( nb_rna, "Leading : %"PRId32" RNAs", gen_unit->get_rna_list()[LEADING]->get_nb_elts() );
  win->draw_string( 15, 15, nb_rna );
  sprintf( nb_rna, "Lagging : %"PRId32" RNAs", gen_unit->get_rna_list()[LAGGING]->get_nb_elts() );
  win->draw_string( 15, 30, nb_rna );

  // Compute display diameter according to genome length and window size
  int16_t win_size      = ae_utils::min( win->get_width(), win->get_height() );
  int16_t diam          = round( win_size * log( (double)genome_length ) / 16 );

  // Prevent diameter from getting greater than 2/3 of the window size
  if ( diam > 2 * win_size / 3 )
  {
    diam = 2 * win_size / 3;
  }

  // Compute coordinates of the upper-left corner of the containing square
  int16_t pos_x = (win->get_width() - diam) / 2;
  int16_t pos_y = (win->get_height() - diam) / 2;

  // Draw main circle
  win->draw_circle( pos_x, pos_y, diam );

  // Sector occupation management
  reset_sectors();


  // ---------------
  //  Draw each RNA
  // ---------------
  ae_list_node<ae_rna*>* rna_node  = NULL;
  ae_rna*       rna       = NULL;

  // NB : As we want OriC to be at the "top" of the circle and the orientation
  //      to be clockwise, the drawing angle (theta) will be given as
  //      (90 - alpha), alpha being the "classical" trigonometric angle
  int16_t alpha_first, alpha_last; // Angles of first and last transcribed bases from OriC (degrees)
  int16_t theta_first, theta_last; // Transposed angles on the trigonometric circle (degrees)
  int16_t nb_sect;
  // Same as above with precision = 1/64 degree
  int16_t alpha_first_64, alpha_last_64;
  int16_t theta_first_64, theta_last_64;
  int16_t nb_sect_64;

  // ----------------
  //  LEADING strand
  // ----------------
  rna_node = gen_unit->get_rna_list()[LEADING]->get_first();

  while ( rna_node != NULL )
  {
    rna = rna_node->get_obj();

    // Alpha : angles from OriC (in degrees)
    // Theta : angles on the trigonometric circle (in degrees)
    // nb_sect : "length" in degrees of the arc to be drawn
    alpha_first   = (int16_t) round(  360 * ((double)rna->get_first_transcribed_pos() / (double)genome_length ));
    alpha_last    = (int16_t) round(  360 * ((double)rna->get_last_transcribed_pos() / (double)genome_length ));
    theta_first   = ae_utils::mod( 90 - alpha_first, 360 );
    theta_last    = ae_utils::mod( 90 - alpha_last, 360 );
    nb_sect       = ae_utils::mod( alpha_last - alpha_first + 1,  360 );

    // These are the same as above but with a higher precision (1/64 degrees)
    alpha_first_64   = (int16_t) round(64 * 360 * ((double)rna->get_first_transcribed_pos() / (double)genome_length ));
    alpha_last_64    = (int16_t) round(64 * 360 * ((double)rna->get_last_transcribed_pos() / (double)genome_length ));
    theta_first_64   = ae_utils::mod( 64 * 90 - alpha_first_64, 64 * 360 );
    theta_last_64    = ae_utils::mod( 64 * 90 - alpha_last_64, 64 * 360 );
    nb_sect_64       = ae_utils::mod( alpha_last_64 - alpha_first_64 + 1,  64 * 360 );

    //~ printf( "    LEADING RNA %"PRId32" => %"PRId32" :: %"PRId16" %"PRId16"\n", rna->get_first_transcribed_pos(), rna->get_last_transcribed_pos(), theta_first, theta_last );


    // Look for the inmost layer that has all the sectors between
    // theta_first and theta_last free
    int16_t layer = 0;
    bool sectors_free = false;
    while ( ! sectors_free )
    {
      sectors_free = true;

      for ( int16_t rho = 0 ; rho < nb_sect ; rho++ )
      {
        if ( _occupied_sectors[LEADING][layer][ae_utils::mod(theta_first-rho, 360)] )
        {
          sectors_free = false;
          break;
        }
      }

      if ( sectors_free )
      {
        break; // All the needed sectors are free on the current layer
      }
      else
      {
        layer++;

        if ( layer >= _outmost_layer )
        {
          add_layer();
          break; // An added layer is necessarily free, no need to look further
        }
      }
    }

    // Mark sectors to be drawn as occupied
    for ( int16_t rho = 0 ; rho < nb_sect ; rho++ )
    {
      _occupied_sectors[LEADING][layer][ae_utils::mod(theta_first-rho, 360)] = true;
    }
    // Mark flanking sectors as occupied
    _occupied_sectors[LEADING][layer][ae_utils::mod(theta_first+1, 360)] = true;
    _occupied_sectors[LEADING][layer][ae_utils::mod(theta_first-nb_sect, 360)] = true;

    
    // Determine drawing color
    char* color;
    if ( rna->is_coding() )
    {
      color = ae_X11_window::get_color( rna->get_basal_level() );
    }
    else
    {
      color = new char[8];
      strcpy( color, "#FFFFFF" );
    }

    // Draw arrow body
    const int8_t rna_layer_spacing = 5; // TODO : param?
    layer++; // index starting at 0 but needed to start at 1

    int16_t diam2 = diam + (layer * 2 * rna_layer_spacing);
    pos_x         = (int16_t) round( (double)(win->get_width()  - diam2 ) / 2.0 );
    pos_y         = (int16_t) round( (double)(win->get_height() - diam2 ) / 2.0 );

    win->draw_arc_64( pos_x, pos_y, diam2, theta_first_64 - nb_sect_64, nb_sect_64, color );

    // Draw arrow head
    int8_t arrow_thick = 6; // Must be an even value
    pos_x = (win->get_width() / 2.0) + (cos((theta_first_64 - nb_sect_64)/(64*180.0)*M_PI) * diam2 / 2.0) - (arrow_thick / 2.0);
    pos_y = (win->get_height() / 2.0) - (sin((theta_first_64 - nb_sect_64)/(64*180.0)*M_PI) * diam2 / 2.0) - (arrow_thick / 2.0);

    win->fill_arc( pos_x, pos_y, arrow_thick, ae_utils::mod(180+theta_last, 360), 180, color );
    
    delete [] color;

    rna_node = rna_node->get_next();
  }

  // ----------------
  //  LAGGING strand
  // ----------------
  rna_node = gen_unit->get_rna_list()[LAGGING]->get_first();

  while ( rna_node != NULL )
  {
    rna = rna_node->get_obj();

    // Alpha : angles from OriC (in degrees)
    // Theta : angles on the trigonometric circle (in degrees)
    // nb_sect : "length" in degrees of the arc to be drawn
    alpha_first   = (int16_t) round(  360 * ((double)rna->get_first_transcribed_pos() / (double)genome_length ));
    alpha_last    = (int16_t) round(  360 * ((double)rna->get_last_transcribed_pos()  / (double)genome_length ));
    theta_first   = ae_utils::mod( 90 - alpha_first, 360 );
    theta_last    = ae_utils::mod( 90 - alpha_last, 360 );
    nb_sect = ae_utils::mod( alpha_first - alpha_last + 1,  360 );

    // These are the same as above but with a higher precision (1/64 degrees)
    alpha_first_64   = (int16_t) round(64 * 360 * ((double)rna->get_first_transcribed_pos() / (double)genome_length ));
    alpha_last_64    = (int16_t) round(64 * 360 * ((double)rna->get_last_transcribed_pos()  / (double)genome_length ));
    theta_first_64   = ae_utils::mod( 64 * 90 - alpha_first_64, 64 * 360 );
    theta_last_64    = ae_utils::mod( 64 * 90 - alpha_last_64, 64 * 360 );
    nb_sect_64 = ae_utils::mod( alpha_first_64 - alpha_last_64 + 1,  64 * 360 );

    //~ printf( "    LAGGING RNA %"PRId32" => %"PRId32" :: %"PRId16" %"PRId16"\n", rna->get_first_transcribed_pos(), rna->get_last_transcribed_pos(), theta_first, theta_last );


    // Look for the inmost layer that has all the sectors between
    // theta_first and theta_last free
    int16_t layer = 0;
    bool sectors_free = false;
    while ( ! sectors_free )
    {
      sectors_free = true;

      for ( int16_t rho = 0 ; rho < nb_sect ; rho++ )
      {
        if ( _occupied_sectors[LAGGING][layer][ae_utils::mod(theta_first+rho, 360)] )
        {
          sectors_free = false;
          break;
        }
      }

      if ( sectors_free )
      {
        break; // All the needed sectors are free on the current layer
      }
      else
      {
        layer++;

        if ( layer >= _outmost_layer )
        {
          add_layer();
          break; // An added layer is necessarily free, no need to look further
        }
      }
    }

    // Mark sectors to be drawn as occupied
    for ( int16_t rho = 0 ; rho < nb_sect ; rho++ )
    {
      _occupied_sectors[LAGGING][layer][ae_utils::mod(theta_first+rho, 360)] = true;
    }
    // Mark flanking sectors as occupied
    _occupied_sectors[LAGGING][layer][ae_utils::mod(theta_first-1, 360)] = true;
    _occupied_sectors[LAGGING][layer][ae_utils::mod(theta_first+nb_sect, 360)] = true;

    
    // Determine drawing color
    char* color;
    if ( rna->is_coding() )
    {
      color = ae_X11_window::get_color( rna->get_basal_level() );
    }
    else
    {
      color = new char[8];
      strcpy( color, "#FFFFFF" );
    }

    // Draw arrow body
    const int8_t rna_layer_spacing = 5; // TODO : param?
    layer++; // index starting at 0 but needed to start at 1

    int16_t diam2 = diam - (layer * 2 * rna_layer_spacing);
    pos_x         = (int16_t) round( (double)(win->get_width()  - diam2 ) / 2.0 );
    pos_y         = (int16_t) round( (double)(win->get_height() - diam2 ) / 2.0 );

    win->draw_arc_64( pos_x, pos_y, diam2, theta_first_64, nb_sect_64, color );

    // Draw arrow head
    int8_t arrow_thick = 6; // Must be an even value
    pos_x = (win->get_width() / 2.0) + (cos((theta_last_64+1)/(64*180.0)*M_PI) * diam2 / 2.0) - (arrow_thick / 2.0);
    pos_y = (win->get_height() / 2.0) - (sin((theta_last_64+1)/(64*180.0)*M_PI) * diam2 / 2.0) - (arrow_thick / 2.0);

    win->fill_arc( pos_x, pos_y, arrow_thick, theta_last, 180, color );
    
    delete [] color;

    rna_node = rna_node->get_next();
  }
}

// =================================================================
//                           Protected Methods
// =================================================================
void ae_individual_X11::add_layer( void )
{
  _occupied_sectors[LEADING][_outmost_layer] = new bool[360];
  _occupied_sectors[LAGGING][_outmost_layer] = new bool[360];

  for ( int16_t angle = 0 ; angle < 360 ; angle++ )
  {
    _occupied_sectors[LEADING][_outmost_layer][angle] = false;
    _occupied_sectors[LAGGING][_outmost_layer][angle] = false;
  }

  _outmost_layer++;
}

void ae_individual_X11::init_occupied_sectors( void )
{
  _outmost_layer = 1;

  for ( int16_t layer = 0 ; layer < _outmost_layer ; layer++ )
  {
    _occupied_sectors[LEADING][layer] = new bool[360];
    _occupied_sectors[LAGGING][layer] = new bool[360];
  }
}

void ae_individual_X11::reset_sectors( void )
{
  for ( int16_t layer = 0 ; layer < _outmost_layer ; layer++ )
  {
    for ( int16_t angle = 0 ; angle < 360 ; angle++ )
    {
      _occupied_sectors[LEADING][layer][angle] = false;
      _occupied_sectors[LAGGING][layer][angle] = false;
    }
  }
}
