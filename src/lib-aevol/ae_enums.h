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



#ifndef __AE_ENUMS_H__
#define  __AE_ENUMS_H__

enum ae_align_fun_shape
{
  LINEAR  = 0,
  SIGMOID = 1
};

enum ae_env_axis_feature
{
  NEUTRAL     = 0,
  METABOLISM  = 1,
  SECRETION   = 2,
  TRANSFER    = 3,
  NB_FEATURES = 4 // This is used to know how many possible features exist to make them easy to parse.
};

enum ae_env_var
{
  NONE                    = 0,
  AUTOREGRESSIVE_MEAN_VAR = 1,
  LOCAL_GAUSSIANS_VAR     = 2
};

enum ae_init_method
{
  ONE_GOOD_GENE   = 0x01,
  CLONE           = 0x02,
  WITH_INS_SEQ    = 0x04
};

enum ae_sense
{
  DIRECT      = 0,
  INDIRECT    = 1,
  BOTH_SENSES = 2
};
  
enum ae_selection_scheme
{
  RANK_LINEAR           = 0,
  RANK_EXPONENTIAL      = 1,
  FITNESS_PROPORTIONATE = 2,
  FITTEST               = 3
};

enum ae_strand
{
  LEADING = 0,
  LAGGING = 1
};
  
enum ae_tree_mode
{
  LIGHT   = 0,
  NORMAL  = 1
};

enum ae_log_type
{
  LOG_TRANSFER  = 0x01,
  LOG_REAR      = 0x02,
  LOG_BARRIER   = 0x04,
  LOG_LOADS     = 0x08
};


#endif // __AE_ENUMS_H__
