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
#include <math.h>



// =================================================================
//                            Project Files
// =================================================================
#include "Individual_X11.h"

#include "ExpManager.h"
#include "ExpSetup.h"
#include "Utils.h"

namespace aevol {

//##############################################################################
//                                                                             #
//                           Class Individual_X11                           #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
Individual_X11::Individual_X11(ExpManager * exp_manager, gzFile backup_file) :
    Individual(exp_manager, backup_file)
{
  init_occupied_sectors();
}

//Individual_X11::Individual_X11(const Individual_X11 &model,
//                                     bool replication_report_copy) :
//    Individual(model, replication_report_copy)
//{
//  init_occupied_sectors();
//}

Individual_X11::Individual_X11(Individual_X11 * const parent,
                                     int32_t id,
                                     std::shared_ptr<JumpingMT> mut_prng,
                                     std::shared_ptr<JumpingMT> stoch_prng) :
    Individual(parent, id, mut_prng, stoch_prng)
{
  init_occupied_sectors();
}

/*Individual_X11::Individual_X11(char* genome, int32_t genome_size) : Individual(genome, genome_size)
{
  init_occupied_sectors();
}*/

// =================================================================
//                             Destructors
// =================================================================
Individual_X11::~Individual_X11(void)
{
  for (int16_t layer = 0 ; layer < _outmost_layer ; layer++)
  {
    delete [] _occupied_sectors[LEADING][layer];
    delete [] _occupied_sectors[LAGGING][layer];
  }
}

// =================================================================
//                            Public Methods
// =================================================================
void Individual_X11::display(void)
{
}

void Individual_X11::display_cdss(X11Window * win)
{
  // Retreive the genetic unit corresponding to the main chromosome
  GeneticUnit* gen_unit = &_genetic_unit_list.front();
  int32_t genome_length = gen_unit->get_dna()->length();

  // Display the number of CDSs
  char display_string[40];
  sprintf(display_string, "Main chromosome size : %" PRId32 "bp", genome_length);
  win->draw_string(15, 25, display_string);
  sprintf(display_string, "Leading : %" PRId32 " CDSs", static_cast<int32_t>(gen_unit->get_protein_list(LEADING).size()));
  win->draw_string(15, 40, display_string);
  sprintf(display_string, "Lagging : %" PRId32 " CDSs", static_cast<int32_t>(gen_unit->get_protein_list(LAGGING).size()));
  win->draw_string(15, 55, display_string);

  // Compute display diameter according to genome length and window size
  int16_t canvas_width;
  if (_allow_plasmids) canvas_width = win->get_width() / 2;
  else canvas_width = win->get_width();
  int16_t canvas_height = win->get_height();

  int16_t canvas_size = Utils::min(canvas_width, canvas_width);
  int16_t diam        = round(canvas_size * log((double)genome_length) / 16);

  // Prevent diameter from getting greater than 2/3 of the window size
  if (diam > 2 * canvas_size / 3)
  {
    diam = 2 * canvas_size / 3;
  }

  // Compute coordinates of the upper-left corner of the containing square
  int16_t pos_x = (canvas_width - diam) / 2;
  int16_t pos_y = (canvas_height - diam) / 2;

  // Draw main circle
  win->draw_circle(pos_x, pos_y, diam);

  // Sector occupation management
  reset_sectors();


  // ---------------
  //  Draw each CDS
  // ---------------
  // NB : As we want OriC to be at the "top" of the circle and the orientation
  //      to be clockwise, the drawing angle (theta) will be given as
  //      (90 - alpha), alpha being the "classical" trigonometric angle
  int16_t alpha_first, alpha_last; // Angles of first and last transcribed bases from OriC (degrees)
  int16_t theta_first; //, theta_last; // Transposed angles on the trigonometric circle (degrees)
  int16_t nb_sect;
  // Same as above with precision = 1/64 degree
  int16_t alpha_first_64, alpha_last_64;
  int16_t theta_first_64; //, theta_last_64;
  int16_t nb_sect_64;

  // ----------------
  //  LEADING strand
  // ----------------
  for (const auto& cds: gen_unit->get_protein_list(LEADING))
  {
    // Alpha : angles from OriC (in degrees)
    // Theta : angles on the trigonometric circle (in degrees)
    // nb_sect : "length" in degrees of the arc to be drawn
    alpha_first   = (int16_t) round( 360 * ((double)cds.get_first_translated_pos() / (double)genome_length));
    alpha_last    = (int16_t) round( 360 * ((double)cds.get_last_translated_pos()  / (double)genome_length));
    theta_first   = Utils::mod(90 - alpha_first, 360);
    // theta_last    = Utils::mod(90 - alpha_last, 360);
    nb_sect       = Utils::mod(alpha_last - alpha_first + 1,  360);

    // These are the same as above but with a higher precision (1/64 degrees)
    alpha_first_64   = (int16_t) round(64 * 360 * ((double)cds.get_first_translated_pos() / (double)genome_length));
    alpha_last_64    = (int16_t) round(64 * 360 * ((double)cds.get_last_translated_pos() / (double)genome_length));
    theta_first_64   = Utils::mod(64 * 90 - alpha_first_64, 64 * 360);
    // theta_last_64    = Utils::mod(64 * 90 - alpha_last_64, 64 * 360);
    nb_sect_64       = Utils::mod(alpha_last_64 - alpha_first_64 + 1,  64 * 360);


    // Look for the inmost layer that has all the sectors between
    // theta_first and theta_last free
    int16_t layer = 0;
    bool sectors_free = false;
    while (! sectors_free)
    {
      sectors_free = true;

      for (int16_t rho = 0 ; rho < nb_sect ; rho++)
      {
        if (_occupied_sectors[LEADING][layer][Utils::mod(theta_first-rho, 360)])
        {
          sectors_free = false;
          break;
        }
      }

      if (sectors_free)
      {
        break; // All the needed sectors are free on the current layer
      }
      else
      {
        layer++;

        if (layer >= _outmost_layer)
        {
          add_layer();
          break; // An added layer is necessarily free, no need to look further
        }
      }
    }

    // Mark sectors to be drawn as occupied
    for (int16_t rho = 0 ; rho < nb_sect ; rho++)
    {
      _occupied_sectors[LEADING][layer][Utils::mod(theta_first-rho, 360)] = true;
    }
    // Mark flanking sectors as occupied
    _occupied_sectors[LEADING][layer][Utils::mod(theta_first+1, 360)] = true;
    _occupied_sectors[LEADING][layer][Utils::mod(theta_first-nb_sect, 360)] = true;


    // Draw
    const int8_t cds_layer_spacing = 5; // TODO : param?
    layer++; // index starting at 0 but needed to start at 1

    int16_t diam2 = diam + (layer * 2 * cds_layer_spacing);
    pos_x         = (int16_t) round((double)(canvas_width  - diam2) / 2.0);
    pos_y         = (int16_t) round((double)(canvas_height - diam2) / 2.0);

    char* color = X11Window::get_color(cds.get_mean());
    win->draw_arc_64(pos_x, pos_y, diam2, theta_first_64 - nb_sect_64, nb_sect_64, color);
    delete [] color;
  }

  // ----------------
  //  LAGGING strand
  // ----------------
  for (const auto& cds: gen_unit->get_protein_list(LAGGING))
  {
    // Alpha : angles from OriC (in degrees)
    // Theta : angles on the trigonometric circle (in degrees)
    // nb_sect : "length" in degrees of the arc to be drawn
    alpha_first   = (int16_t) round( 360 * ((double)cds.get_first_translated_pos() / (double)genome_length));
    alpha_last    = (int16_t) round( 360 * ((double)cds.get_last_translated_pos()  / (double)genome_length));
    theta_first   = Utils::mod(90 - alpha_first, 360);
    // theta_last    = Utils::mod(90 - alpha_last, 360);
    nb_sect = Utils::mod(alpha_first - alpha_last + 1,  360);

    // These are the same as above but with a higher precision (1/64 degrees)
    alpha_first_64   = (int16_t) round(64 * 360 * ((double)cds.get_first_translated_pos() / (double)genome_length));
    alpha_last_64    = (int16_t) round(64 * 360 * ((double)cds.get_last_translated_pos() / (double)genome_length));
    theta_first_64   = Utils::mod(64 * 90 - alpha_first_64, 64 * 360);
    // theta_last_64    = Utils::mod(64 * 90 - alpha_last_64, 64 * 360);
    nb_sect_64 = Utils::mod(alpha_first_64 - alpha_last_64 + 1,  64 * 360);


    // Look for the inmost layer that has all the sectors between
    // theta_first and theta_last free
    int16_t layer = 0;
    bool sectors_free = false;
    while (! sectors_free)
    {
      sectors_free = true;

      for (int16_t rho = 0 ; rho < nb_sect ; rho++)
      {
        if (_occupied_sectors[LAGGING][layer][Utils::mod(theta_first+rho, 360)])
        {
          sectors_free = false;
          break;
        }
      }

      if (sectors_free)
      {
        break; // All the needed sectors are free on the current layer
      }
      else
      {
        layer++;

        if (layer >= _outmost_layer)
        {
          add_layer();
          break; // An added layer is necessarily free, no need to look further
        }
      }
    }

    // Mark sectors to be drawn as occupied
    for (int16_t rho = 0 ; rho < nb_sect ; rho++)
    {
      _occupied_sectors[LAGGING][layer][Utils::mod(theta_first+rho, 360)] = true;
    }
    // Mark flanking sectors as occupied
    _occupied_sectors[LAGGING][layer][Utils::mod(theta_first-1, 360)] = true;
    _occupied_sectors[LAGGING][layer][Utils::mod(theta_first+nb_sect, 360)] = true;


    // Draw
    const int8_t cds_layer_spacing = 5; // TODO : param?
    layer++; // index starting at 0 but needed to start at 1

    int16_t diam2 = diam - (layer * 2 * cds_layer_spacing);
    pos_x         = (int16_t) round((double)(canvas_width  - diam2) / 2.0);
    pos_y         = (int16_t) round((double)(canvas_height - diam2) / 2.0);

    char* color = X11Window::get_color(cds.get_mean());
    win->draw_arc_64(pos_x, pos_y, diam2, theta_first_64, nb_sect_64, color);
    delete [] color;
  }









  // --------------------------------------------------------------------------This is temporary, it is a big copy-paste of what's above.
  if (_allow_plasmids)
  {
    // Retreive the genetic unit corresponding to the plasmid (i.e. index 1 in _genetic_unit_list)
    GeneticUnit* gen_unit = &*std::next(_genetic_unit_list.begin());
    if (gen_unit == NULL) return;

    int32_t genome_length = gen_unit->get_dna()->length();


    // Compute display diameter according to genome length and window size
    int16_t canvas_width;
    int16_t canvas_height;
    int16_t canvas_size ;
    if (_allow_plasmids)
    {
      canvas_width  = win->get_width() / 2;
      canvas_size   = canvas_width;
    }
    else
    {
      canvas_width  = win->get_width();
      canvas_size   = Utils::min(canvas_width, canvas_width);
    }
    canvas_height = win->get_height();

    int16_t diam  = round(canvas_size * log((double)genome_length) / 16);

    // Prevent diameter from getting greater than 2/3 of the window size
    if (diam > 2 * canvas_size / 3)
    {
      diam = 2 * canvas_size / 3;
    }

    // Compute coordinates of the upper-left corner of the containing square
    int16_t pos_x = canvas_width + (canvas_width - diam) / 2;
    int16_t pos_y = (canvas_height - diam) / 2;


    // Draw main circle
    win->draw_circle(pos_x, pos_y, diam);

    // Sector occupation management
    reset_sectors();


    // ---------------
    //  Draw each CDS
    // ---------------
    // NB : As we want OriC to be at the "top" of the circle and the orientation
    //      to be clockwise, the drawing angle (theta) will be given as
    //      (90 - alpha), alpha being the "classical" trigonometric angle
    int16_t alpha_first, alpha_last; // Angles of first and last transcribed bases from OriC (degrees)
    int16_t theta_first; //, theta_last; // Transposed angles on the trigonometric circle (degrees)
    int16_t nb_sect;
    // Same as above with precision = 1/64 degree
    int16_t alpha_first_64, alpha_last_64;
    int16_t theta_first_64; //, theta_last_64;
    int16_t nb_sect_64;

    // ----------------
    //  LEADING strand
    // ----------------
    for (const auto& cds: gen_unit->get_protein_list(LEADING))
    {
      // Alpha : angles from OriC (in degrees)
      // Theta : angles on the trigonometric circle (in degrees)
      // nb_sect : "length" in degrees of the arc to be drawn
      alpha_first   = (int16_t) round( 360 * ((double)cds.get_first_translated_pos() / (double)genome_length));
      alpha_last    = (int16_t) round( 360 * ((double)cds.get_last_translated_pos()  / (double)genome_length));
      theta_first   = Utils::mod(90 - alpha_first, 360);
      // theta_last    = Utils::mod(90 - alpha_last, 360);
      nb_sect       = Utils::mod(alpha_last - alpha_first + 1,  360);

      // These are the same as above but with a higher precision (1/64 degrees)
      alpha_first_64   = (int16_t) round(64 * 360 * ((double)cds.get_first_translated_pos() / (double)genome_length));
      alpha_last_64    = (int16_t) round(64 * 360 * ((double)cds.get_last_translated_pos() / (double)genome_length));
      theta_first_64   = Utils::mod(64 * 90 - alpha_first_64, 64 * 360);
      // theta_last_64    = Utils::mod(64 * 90 - alpha_last_64, 64 * 360);
      nb_sect_64       = Utils::mod(alpha_last_64 - alpha_first_64 + 1,  64 * 360);


      // Look for the inmost layer that has all the sectors between
      // theta_first and theta_last free
      int16_t layer = 0;
      bool sectors_free = false;
      while (! sectors_free)
      {
        sectors_free = true;

        for (int16_t rho = 0 ; rho < nb_sect ; rho++)
        {
          if (_occupied_sectors[LEADING][layer][Utils::mod(theta_first-rho, 360)])
          {
            sectors_free = false;
            break;
          }
        }

        if (sectors_free)
        {
          break; // All the needed sectors are free on the current layer
        }
        else
        {
          layer++;

          if (layer >= _outmost_layer)
          {
            add_layer();
            break; // An added layer is necessarily free, no need to look further
          }
        }
      }

      // Mark sectors to be drawn as occupied
      for (int16_t rho = 0 ; rho < nb_sect ; rho++)
      {
        _occupied_sectors[LEADING][layer][Utils::mod(theta_first-rho, 360)] = true;
      }
      // Mark flanking sectors as occupied
      _occupied_sectors[LEADING][layer][Utils::mod(theta_first+1, 360)] = true;
      _occupied_sectors[LEADING][layer][Utils::mod(theta_first-nb_sect, 360)] = true;


      // Draw
      const int8_t cds_layer_spacing = 5; // TODO : param?
      layer++; // index starting at 0 but needed to start at 1

      int16_t diam2 = diam + (layer * 2 * cds_layer_spacing);
      pos_x         = canvas_width + (int16_t) round((double)(canvas_width  - diam2) / 2.0);
      pos_y         = (int16_t) round((double)(canvas_height - diam2) / 2.0);

      char* color = X11Window::get_color(cds.get_mean());
      win->draw_arc_64(pos_x, pos_y, diam2, theta_first_64 - nb_sect_64, nb_sect_64, color);
      delete [] color;
    }

    // ----------------
    //  LAGGING strand
    // ----------------
    for (const auto& cds: gen_unit->get_protein_list(LAGGING))
    {
      // Alpha : angles from OriC (in degrees)
      // Theta : angles on the trigonometric circle (in degrees)
      // nb_sect : "length" in degrees of the arc to be drawn
      alpha_first   = (int16_t) round( 360 * ((double)cds.get_first_translated_pos() / (double)genome_length));
      alpha_last    = (int16_t) round( 360 * ((double)cds.get_last_translated_pos()  / (double)genome_length));
      theta_first   = Utils::mod(90 - alpha_first, 360);
      // theta_last    = Utils::mod(90 - alpha_last, 360);
      nb_sect = Utils::mod(alpha_first - alpha_last + 1,  360);

      // These are the same as above but with a higher precision (1/64 degrees)
      alpha_first_64   = (int16_t) round(64 * 360 * ((double)cds.get_first_translated_pos() / (double)genome_length));
      alpha_last_64    = (int16_t) round(64 * 360 * ((double)cds.get_last_translated_pos() / (double)genome_length));
      theta_first_64   = Utils::mod(64 * 90 - alpha_first_64, 64 * 360);
      // theta_last_64    = Utils::mod(64 * 90 - alpha_last_64, 64 * 360);
      nb_sect_64 = Utils::mod(alpha_first_64 - alpha_last_64 + 1,  64 * 360);


      // Look for the inmost layer that has all the sectors between
      // theta_first and theta_last free
      int16_t layer = 0;
      bool sectors_free = false;
      while (! sectors_free)
      {
        sectors_free = true;

        for (int16_t rho = 0 ; rho < nb_sect ; rho++)
        {
          if (_occupied_sectors[LAGGING][layer][Utils::mod(theta_first+rho, 360)])
          {
            sectors_free = false;
            break;
          }
        }

        if (sectors_free)
        {
          break; // All the needed sectors are free on the current layer
        }
        else
        {
          layer++;

          if (layer >= _outmost_layer)
          {
            add_layer();
            break; // An added layer is necessarily free, no need to look further
          }
        }
      }

      // Mark sectors to be drawn as occupied
      for (int16_t rho = 0 ; rho < nb_sect ; rho++)
      {
        _occupied_sectors[LAGGING][layer][Utils::mod(theta_first+rho, 360)] = true;
      }
      // Mark flanking sectors as occupied
      _occupied_sectors[LAGGING][layer][Utils::mod(theta_first-1, 360)] = true;
      _occupied_sectors[LAGGING][layer][Utils::mod(theta_first+nb_sect, 360)] = true;


      // Draw
      const int8_t cds_layer_spacing = 5; // TODO : param?
      layer++; // index starting at 0 but needed to start at 1

      int16_t diam2 = diam - (layer * 2 * cds_layer_spacing);
      pos_x         = canvas_width + (int16_t) round((double)(canvas_width  - diam2) / 2.0);
      pos_y         = (int16_t) round((double)(canvas_height - diam2) / 2.0);

      char* color = X11Window::get_color(cds.get_mean());
      win->draw_arc_64(pos_x, pos_y, diam2, theta_first_64, nb_sect_64, color);
      delete [] color;
    }
  }
}

void Individual_X11::display_rnas(X11Window * win)
{
  // Retreive the genetic unit corresponding to the main chromosome
  const GeneticUnit* gen_unit = &_genetic_unit_list.front();
  int32_t genome_length = gen_unit->get_dna()->length();

  // Display the number of RNAs
  char nb_rna[40];
  sprintf(nb_rna, "Leading : %" PRId32 " RNAs", static_cast<int32_t>(gen_unit->get_rna_list()[LEADING].size()));
  win->draw_string(15, 15, nb_rna);
  sprintf(nb_rna, "Lagging : %" PRId32 " RNAs", static_cast<int32_t>(gen_unit->get_rna_list()[LAGGING].size()));
  win->draw_string(15, 30, nb_rna);

  // Compute display diameter according to genome length and window size
  int16_t win_size      = Utils::min(win->get_width(), win->get_height());
  int16_t diam          = round(win_size * log((double)genome_length) / 16);

  // Prevent diameter from getting greater than 2/3 of the window size
  if (diam > 2 * win_size / 3)
  {
    diam = 2 * win_size / 3;
  }

  // Compute coordinates of the upper-left corner of the containing square
  int16_t pos_x = (win->get_width() - diam) / 2;
  int16_t pos_y = (win->get_height() - diam) / 2;

  // Draw main circle
  win->draw_circle(pos_x, pos_y, diam);

  // Sector occupation management
  reset_sectors();


  // ---------------
  //  Draw each RNA
  // ---------------
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
  const auto& rna_lists = gen_unit->get_rna_list();
  for (const auto& rna: rna_lists[LEADING])
  {
    // Alpha : angles from OriC (in degrees)
    // Theta : angles on the trigonometric circle (in degrees)
    // nb_sect : "length" in degrees of the arc to be drawn
    alpha_first   = (int16_t) round( 360 * ((double)rna.get_first_transcribed_pos() / (double)genome_length));
    alpha_last    = (int16_t) round( 360 * ((double)rna.get_last_transcribed_pos() / (double)genome_length));
    theta_first   = Utils::mod(90 - alpha_first, 360);
    theta_last    = Utils::mod(90 - alpha_last, 360);
    nb_sect       = Utils::mod(alpha_last - alpha_first + 1,  360);

    // These are the same as above but with a higher precision (1/64 degrees)
    alpha_first_64   = (int16_t) round(64 * 360 * ((double)rna.get_first_transcribed_pos() / (double)genome_length));
    alpha_last_64    = (int16_t) round(64 * 360 * ((double)rna.get_last_transcribed_pos() / (double)genome_length));
    theta_first_64   = Utils::mod(64 * 90 - alpha_first_64, 64 * 360);
    theta_last_64    = Utils::mod(64 * 90 - alpha_last_64, 64 * 360);
    nb_sect_64       = Utils::mod(alpha_last_64 - alpha_first_64 + 1,  64 * 360);

    //~ printf("    LEADING RNA %"PRId32" => %"PRId32" :: %"PRId16" %"PRId16"\n", rna->get_first_transcribed_pos(), rna->get_last_transcribed_pos(), theta_first, theta_last);


    // Look for the inmost layer that has all the sectors between
    // theta_first and theta_last free
    int16_t layer = 0;
    bool sectors_free = false;
    while (! sectors_free)
    {
      sectors_free = true;

      for (int16_t rho = 0 ; rho < nb_sect ; rho++)
      {
        if (_occupied_sectors[LEADING][layer][Utils::mod(theta_first-rho, 360)])
        {
          sectors_free = false;
          break;
        }
      }

      if (sectors_free)
      {
        break; // All the needed sectors are free on the current layer
      }
      else
      {
        layer++;

        if (layer >= _outmost_layer)
        {
          add_layer();
          break; // An added layer is necessarily free, no need to look further
        }
      }
    }

    // Mark sectors to be drawn as occupied
    for (int16_t rho = 0 ; rho < nb_sect ; rho++)
    {
      _occupied_sectors[LEADING][layer][Utils::mod(theta_first-rho, 360)] = true;
    }
    // Mark flanking sectors as occupied
    _occupied_sectors[LEADING][layer][Utils::mod(theta_first+1, 360)] = true;
    _occupied_sectors[LEADING][layer][Utils::mod(theta_first-nb_sect, 360)] = true;


    // Determine drawing color
    char* color;
    if (rna.is_coding())
    {
      color = X11Window::get_color(rna.get_basal_level());
    }
    else
    {
      color = new char[8];
      strcpy(color, "#FFFFFF");
    }

    // Draw arrow body
    const int8_t rna_layer_spacing = 5; // TODO : param?
    layer++; // index starting at 0 but needed to start at 1

    int16_t diam2 = diam + (layer * 2 * rna_layer_spacing);
    pos_x         = (int16_t) round((double)(win->get_width()  - diam2) / 2.0);
    pos_y         = (int16_t) round((double)(win->get_height() - diam2) / 2.0);

    win->draw_arc_64(pos_x, pos_y, diam2, theta_first_64 - nb_sect_64, nb_sect_64, color);

    // Draw arrow head
    int8_t arrow_thick = 6; // Must be an even value
    pos_x = (win->get_width() / 2.0) + (cos((theta_first_64 - nb_sect_64)/(64*180.0)*M_PI) * diam2 / 2.0) - (arrow_thick / 2.0);
    pos_y = (win->get_height() / 2.0) - (sin((theta_first_64 - nb_sect_64)/(64*180.0)*M_PI) * diam2 / 2.0) - (arrow_thick / 2.0);

    win->fill_arc(pos_x, pos_y, arrow_thick, Utils::mod(180+theta_last, 360), 180, color);

    delete [] color;
  }

  // ----------------
  //  LAGGING strand
  // ----------------
  for (const auto& rna: rna_lists[LAGGING]) {
    // Alpha : angles from OriC (in degrees)
    // Theta : angles on the trigonometric circle (in degrees)
    // nb_sect : "length" in degrees of the arc to be drawn
    alpha_first   = (int16_t) round( 360 * ((double)rna.get_first_transcribed_pos() / (double)genome_length));
    alpha_last    = (int16_t) round( 360 * ((double)rna.get_last_transcribed_pos()  / (double)genome_length));
    theta_first   = Utils::mod(90 - alpha_first, 360);
    theta_last    = Utils::mod(90 - alpha_last, 360);
    nb_sect = Utils::mod(alpha_first - alpha_last + 1,  360);

    // These are the same as above but with a higher precision (1/64 degrees)
    alpha_first_64   = (int16_t) round(64 * 360 * ((double)rna.get_first_transcribed_pos() / (double)genome_length));
    alpha_last_64    = (int16_t) round(64 * 360 * ((double)rna.get_last_transcribed_pos()  / (double)genome_length));
    theta_first_64   = Utils::mod(64 * 90 - alpha_first_64, 64 * 360);
    theta_last_64    = Utils::mod(64 * 90 - alpha_last_64, 64 * 360);
    nb_sect_64 = Utils::mod(alpha_first_64 - alpha_last_64 + 1,  64 * 360);

    //~ printf("    LAGGING RNA %"PRId32" => %"PRId32" :: %"PRId16" %"PRId16"\n", rna->get_first_transcribed_pos(), rna->get_last_transcribed_pos(), theta_first, theta_last);


    // Look for the inmost layer that has all the sectors between
    // theta_first and theta_last free
    int16_t layer = 0;
    bool sectors_free = false;
    while (! sectors_free)
    {
      sectors_free = true;

      for (int16_t rho = 0 ; rho < nb_sect ; rho++)
      {
        if (_occupied_sectors[LAGGING][layer][Utils::mod(theta_first+rho, 360)])
        {
          sectors_free = false;
          break;
        }
      }

      if (sectors_free)
      {
        break; // All the needed sectors are free on the current layer
      }
      else
      {
        layer++;

        if (layer >= _outmost_layer)
        {
          add_layer();
          break; // An added layer is necessarily free, no need to look further
        }
      }
    }

    // Mark sectors to be drawn as occupied
    for (int16_t rho = 0 ; rho < nb_sect ; rho++)
    {
      _occupied_sectors[LAGGING][layer][Utils::mod(theta_first+rho, 360)] = true;
    }
    // Mark flanking sectors as occupied
    _occupied_sectors[LAGGING][layer][Utils::mod(theta_first-1, 360)] = true;
    _occupied_sectors[LAGGING][layer][Utils::mod(theta_first+nb_sect, 360)] = true;


    // Determine drawing color
    char* color;
    if (rna.is_coding())
    {
      color = X11Window::get_color(rna.get_basal_level());
    }
    else
    {
      color = new char[8];
      strcpy(color, "#FFFFFF");
    }

    // Draw arrow body
    const int8_t rna_layer_spacing = 5; // TODO : param?
    layer++; // index starting at 0 but needed to start at 1

    int16_t diam2 = diam - (layer * 2 * rna_layer_spacing);
    pos_x         = (int16_t) round((double)(win->get_width()  - diam2) / 2.0);
    pos_y         = (int16_t) round((double)(win->get_height() - diam2) / 2.0);

    win->draw_arc_64(pos_x, pos_y, diam2, theta_first_64, nb_sect_64, color);

    // Draw arrow head
    int8_t arrow_thick = 6; // Must be an even value
    pos_x = (win->get_width() / 2.0) + (cos((theta_last_64+1)/(64*180.0)*M_PI) * diam2 / 2.0) - (arrow_thick / 2.0);
    pos_y = (win->get_height() / 2.0) - (sin((theta_last_64+1)/(64*180.0)*M_PI) * diam2 / 2.0) - (arrow_thick / 2.0);

    win->fill_arc(pos_x, pos_y, arrow_thick, theta_last, 180, color);

    delete [] color;
  }
}

// =================================================================
//                           Protected Methods
// =================================================================
void Individual_X11::add_layer(void)
{
  _occupied_sectors[LEADING][_outmost_layer] = new bool[360];
  _occupied_sectors[LAGGING][_outmost_layer] = new bool[360];

  for (int16_t angle = 0 ; angle < 360 ; angle++)
  {
    _occupied_sectors[LEADING][_outmost_layer][angle] = false;
    _occupied_sectors[LAGGING][_outmost_layer][angle] = false;
  }

  _outmost_layer++;
}

void Individual_X11::init_occupied_sectors(void)
{
  _outmost_layer = 1;

  for (int16_t layer = 0 ; layer < _outmost_layer ; layer++)
  {
    _occupied_sectors[LEADING][layer] = new bool[360];
    _occupied_sectors[LAGGING][layer] = new bool[360];
  }
}

void Individual_X11::reset_sectors(void)
{
  for (int16_t layer = 0 ; layer < _outmost_layer ; layer++)
  {
    for (int16_t angle = 0 ; angle < 360 ; angle++)
    {
      _occupied_sectors[LEADING][layer][angle] = false;
      _occupied_sectors[LAGGING][layer][angle] = false;
    }
  }
}
} // namespace aevol
