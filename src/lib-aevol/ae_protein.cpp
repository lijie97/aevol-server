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
#include <math.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_codon.h>
#include <ae_protein.h>
#include <ae_genetic_unit.h>
#include <ae_rna.h>
#include <ae_utils.h>




//##############################################################################
//                                                                             #
//                              Class ae_protein                               #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
/**
 * Copy constructor.
 *
 * Copies the protein but does nothing regarding the RNAs transcribing it (creates an empty list).
 */
ae_protein::ae_protein( ae_genetic_unit* gen_unit, const ae_protein &model )
{
  _gen_unit  = gen_unit;
  
  _strand                 = model._strand;
  _shine_dal_pos          = model._shine_dal_pos;
  _first_translated_pos   = model._first_translated_pos;
  _last_translated_pos    = model._last_translated_pos;
  _length                 = model._length;
  _concentration          = model._concentration;
  _is_functional          = model._is_functional;
  
  _rna_list = new ae_list();
  
  // Copy the list of amino-acids
  _AA_list  = new ae_list();
  
  ae_list_node* AA_node = model._AA_list->get_first();
  ae_codon*     AA      = NULL;

  
  while ( AA_node != NULL )
  {
    AA = (ae_codon*) AA_node->get_obj();
    
    _AA_list->add( new ae_codon( *AA ) );
    
    AA_node = AA_node->get_next();
  }
  
  // Copy triangle parameters
  _mean   = model._mean;
  _width  = model._width;
  _height = model._height;
}

ae_protein::ae_protein( ae_genetic_unit* gen_unit, ae_list* codon_list, ae_strand strand, int32_t shine_dal_pos, ae_rna* rna )
{
  _gen_unit       = gen_unit;
  _strand         = strand;
  _shine_dal_pos  = shine_dal_pos;
  _length         = codon_list->get_nb_elts();
  
  // In Aevol the concentration of a new protein is set at the basal level
  #ifndef __REGUL
  _concentration  = rna->get_basal_level();
  // In Raevol, there is two case, depending on the heredity
  #else
  if ( ae_common::with_heredity )
  {
    // With heredity the new protein has a concentration set at 0, because there are inherited proteins which allow the regulation
    _concentration = 0;
  }
  else
  {
    // Without heredity, we use the same concentration as in Aevol (No inherited proteins)
    _concentration = rna->get_basal_level();
  }
  #endif
  
  // TODO : make this cleaner...
  _AA_list = codon_list;
  
  _rna_list = new ae_list();
  _rna_list->add( rna );

  if ( _strand == LEADING )
  {
    _first_translated_pos = ae_utils::mod( _shine_dal_pos + (SHINE_DAL_SIZE + SHINE_START_SPACER + CODON_SIZE), _gen_unit->get_dna()->get_length() );
    _last_translated_pos  = ae_utils::mod( _first_translated_pos + (_length * CODON_SIZE - 1), _gen_unit->get_dna()->get_length() );
  }
  else
  {
    _first_translated_pos = ae_utils::mod( _shine_dal_pos - (SHINE_DAL_SIZE + SHINE_START_SPACER + CODON_SIZE), _gen_unit->get_dna()->get_length() );
    _last_translated_pos = ae_utils::mod( _first_translated_pos - (_length * CODON_SIZE - 1), _gen_unit->get_dna()->get_length() );
  }



  // ============================================================================
  // Folding process (compute phenotypic contribution parameters from codon list)
  // ============================================================================
  //  1) Compute values for M, W and H
  //  2) Normalize M, W and H values according to number of codons of each kind
  //  3) Normalize M, W and H values according to the allowed ranges (defined in macros.h)


  //  --------------------------------
  //  1) Compute values for M, W and H
  //  --------------------------------
  long double M = 0.0;
  long double W = 0.0;
  long double H = 0.0;

  int32_t nb_m = 0;
  int32_t nb_w = 0;
  int32_t nb_h = 0;

  bool bin_m = false; // Initializing to false will yield a conservation of the high weight bit
  bool bin_w = false; // when applying the XOR operator for the Gray to standard conversion
  bool bin_h = false;

  ae_list_node* node = codon_list->get_first();
  ae_codon* codon = NULL;


  while ( node != NULL )
  {
    codon = (ae_codon*) node->get_obj();

    switch ( codon->get_value() )
    {
      case CODON_M0 :
      {
        // M codon found
        nb_m++;

        // Convert Gray code to "standard" binary code
        bin_m ^= false; // as bin_m was initialized to false, the XOR will have no effect on the high weight bit

        // A lower-than-the-previous-lowest weight bit was found, make a left bitwise shift
        //~ M <<= 1;
        M *= 2;

        // Add this nucleotide's contribution to M
        if ( bin_m ) M += 1;

        break;
      }
      case CODON_M1 :
      {
        // M codon found
        nb_m++;

        // Convert Gray code to "standard" binary code
        bin_m ^= true; // as bin_m was initialized to false, the XOR will have no effect on the high weight bit

        // A lower-than-the-previous-lowest bit was found, make a left bitwise shift
        //~ M <<= 1;
        M *= 2;

        // Add this nucleotide's contribution to M
        if ( bin_m ) M += 1;

        break;
      }
      case CODON_W0 :
      {
        // W codon found
        nb_w++;

        // Convert Gray code to "standard" binary code
        bin_w ^= false; // as bin_m was initialized to false, the XOR will have no effect on the high weight bit

        // A lower-than-the-previous-lowest weight bit was found, make a left bitwise shift
        //~ W <<= 1;
        W *= 2;

        // Add this nucleotide's contribution to M
        if ( bin_w ) W += 1;

        break;
      }
      case CODON_W1 :
      {
        // W codon found
        nb_w++;

        // Convert Gray code to "standard" binary code
        bin_w ^= true; // as bin_m was initialized to false, the XOR will have no effect on the high weight bit

        // A lower-than-the-previous-lowest weight bit was found, make a left bitwise shift
        //~ W <<= 1;
        W *= 2;

        // Add this nucleotide's contribution to M
        if ( bin_w ) W += 1;

        break;
      }
      case CODON_H0 :
      case CODON_START : // Start codon codes for the same amino-acid as H0 codon
      {
        // H codon found
        nb_h++;

        // Convert Gray code to "standard" binary code
        bin_h ^= false; // as bin_m was initialized to false, the XOR will have no effect on the high weight bit

        // A lower-than-the-previous-lowest weight bit was found, make a left bitwise shift
        //~ H <<= 1;
        H *= 2;

        // Add this nucleotide's contribution to M
        if ( bin_h ) H += 1;

        break;
      }
      case CODON_H1 :
      {
        // H codon found
        nb_h++;

        // Convert Gray code to "standard" binary code
        bin_h ^= true; // as bin_m was initialized to false, the XOR will have no effect on the high weight bit

        // A lower-than-the-previous-lowest weight bit was found, make a left bitwise shift
        //~ H <<= 1;
        H *= 2;

        // Add this nucleotide's contribution to M
        if ( bin_h ) H += 1;

        break;
      }
    }

    node = node->get_next();
  }



  //  ----------------------------------------------------------------------------------
  //  2) Normalize M, W and H values in [0;1] according to number of codons of each kind
  //  ----------------------------------------------------------------------------------
  if ( nb_m != 0 )  _mean = M / (pow(2, nb_m) - 1);
  else              _mean = 0.5;
  if ( nb_w != 0 )  _width = W / (pow(2, nb_w) - 1);
  else              _width = 0.0;
  if ( nb_h != 0 )  _height = H / (pow(2, nb_h) - 1);
  else              _height = 0.5;

  assert( _mean >= 0.0 && _mean <= 1.0 );
  assert( _width >= 0.0 && _width <= 1.0 );
  assert( _height >= 0.0 && _height <= 1.0 );



  //  ------------------------------------------------------------------------------------
  //  3) Normalize M, W and H values according to the allowed ranges (defined in macros.h)
  //  ------------------------------------------------------------------------------------
  // min_x <= M <= max_x
  // min_w <= W <= max_w
  // min_h <= H <= max_h
  _mean   = (MAX_X - MIN_X) * _mean + MIN_X;
  _width  = (MAX_W - MIN_W) * _width + MIN_W;
  _height = (MAX_H - MIN_H) * _height + MIN_H;
  
  if ( nb_m == 0 || nb_w == 0 || nb_h == 0 || _width == 0.0 || _height == 0.0 )
  {
    _is_functional = false;
  }
  else
  {
    _is_functional = true;
  }

  assert( _mean >= MIN_X && _mean <= MAX_X );
  assert( _width >= MIN_W && _width <= MAX_W );
  assert( _height >= MIN_H && _height <= MAX_H );
}

/*
ae_protein::ae_protein( ae_protein* parent )
{
  _gen_unit             = parent->_gen_unit;
  _strand               = parent->_strand;
  _shine_dal_pos        = parent->_shine_dal_pos;
  _first_translated_pos = parent->_first_translated_pos;
  _last_translated_pos  = parent->_last_translated_pos;
  _length               = parent->_length;
  _concentration        = parent->_concentration;
  _mean                 = parent->_mean;
  _width                = parent->_width;
  _height               = parent->_height;
}
*/

ae_protein::ae_protein( gzFile* backup_file )
{
  // the Rna_list is empty, and is not uselfull there, thus it is not save.
  gzread( backup_file, &_gen_unit,   			      sizeof(_gen_unit)             );
  gzread( backup_file, &_strand,    	 		      sizeof(_strand)               );
  gzread( backup_file, &_shine_dal_pos,			    sizeof(_shine_dal_pos)        );
  gzread( backup_file, &_first_translated_pos, 	sizeof(_first_translated_pos) );
  gzread( backup_file, &_last_translated_pos,  	sizeof(_last_translated_pos)  );
  gzread( backup_file, &_length,     			      sizeof(_length)               );
  gzread( backup_file, &_concentration,     		sizeof(_concentration)        );
  gzread( backup_file, &_is_functional,         sizeof(_is_functional)         );
  gzread( backup_file, &_mean,  			          sizeof(_mean)                 );
  gzread( backup_file, &_width,    			        sizeof(_width)                );
  gzread( backup_file, &_height,                sizeof(_height)               );
  
  _rna_list = new ae_list();

  // Retreive the AA
  _AA_list = new ae_list();
  int16_t nb_AA = 0;
  gzread( backup_file, &nb_AA,  sizeof(nb_AA) );
  
  for ( int16_t i = 0 ; i < nb_AA ; i++ )
  {
    _AA_list->add( new ae_codon( backup_file ) );
  }
  
}

// =================================================================
//                             Destructors
// =================================================================
ae_protein::~ae_protein( void )
{
  _rna_list->erase( NO_DELETE );
  delete _rna_list;
  
  _AA_list->erase( DELETE_OBJ );
  delete _AA_list;
}

// =================================================================
//                            Public Methods
// =================================================================
int32_t ae_protein::get_last_STOP_base_pos( void ) const
{
  if ( _strand == LEADING )
  {
    return ae_utils::mod( _last_translated_pos + 3, _gen_unit->get_dna()->get_length() );
  }
  else
  {
    return ae_utils::mod( _last_translated_pos - 3, _gen_unit->get_dna()->get_length() );
  }
}

void ae_protein::add_RNA( ae_rna* rna )
{
  _rna_list->add( rna );
  _concentration += rna->get_basal_level();
}

char* ae_protein::get_AA_sequence( void ) const
{
  char* seq = new char[3*_length]; // + 1 (for the '\0')  - 1 (_length - 1 spaces)
  
  ae_list_node* codon_node = _AA_list->get_first();
  
  int32_t i = 0;
  while ( codon_node != NULL )
  {
    ae_codon* codon = (ae_codon*) codon_node->get_obj();
    
    if ( i != 0 ) seq[i++] = ' ';
    
    switch ( codon->get_value() )
    {
      case CODON_START :
      {
        seq[i++] = 'S';
        seq[i++] = 'T';
        break;
      }
      case CODON_M0 :
      {
        seq[i++] = 'M';
        seq[i++] = '0';
        break;
      }
      case CODON_M1 :
      {
        seq[i++] = 'M';
        seq[i++] = '1';
        break;
      }
      case CODON_W0 :
      {
        seq[i++] = 'W';
        seq[i++] = '0';
        break;
      }
      case CODON_W1 :
      {
        seq[i++] = 'W';
        seq[i++] = '1';
        break;
      }
      case CODON_H0 :
      {
        seq[i++] = 'H';
        seq[i++] = '0';
        break;
      }
      case CODON_H1 :
      {
        seq[i++] = 'H';
        seq[i++] = '1';
        break;
      }
    }
    
    codon_node = codon_node->get_next();
  }
  
  seq[3*_length-1] = '\0';
  return seq;
}

void ae_protein::write_to_backup( gzFile* backup_file )
{
  // The rna_list is not write because there is no need to, it is an empty list.
  gzwrite( backup_file, &_gen_unit,             sizeof(_gen_unit)             );
  gzwrite( backup_file, &_strand,               sizeof(_strand)               );
  gzwrite( backup_file, &_shine_dal_pos,        sizeof(_shine_dal_pos)        );
  gzwrite( backup_file, &_first_translated_pos, sizeof(_first_translated_pos) );
  gzwrite( backup_file, &_last_translated_pos,  sizeof(_last_translated_pos)  );
  gzwrite( backup_file, &_length,     			    sizeof(_length)               );
  gzwrite( backup_file, &_concentration,     		sizeof(_concentration)        );
  gzwrite( backup_file, &_is_functional,        sizeof(_is_functional)         );
  gzwrite( backup_file, &_mean,  			          sizeof(_mean)                 );
  gzwrite( backup_file, &_width,    			      sizeof(_width)                );
  gzwrite( backup_file, &_height,		     	      sizeof(_height)               );

  // Write the Acide Amino in the backup file
  int16_t nb_AA = _AA_list->get_nb_elts();
  gzwrite( backup_file, &nb_AA,  sizeof(nb_AA) );

  ae_list_node*     AA_node = _AA_list->get_first();
  ae_codon*  AA;

  for ( int16_t i = 0 ; i < nb_AA ; i++ )
  {
  AA = (ae_codon*)AA_node->get_obj();
   
  AA->write_to_backup( backup_file );
    
  AA_node = AA_node->get_next();
  }
}

// =================================================================
//                           Protected Methods
// =================================================================
