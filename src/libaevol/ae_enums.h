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


#ifndef __AE_ENUMS_H__
#define __AE_ENUMS_H__

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
  DONOR       = 3,
  RECIPIENT   = 4
};

#define NB_FEATURES 5 // This is used to know how many possible features exist to make them easy to parse.

enum ae_env_var
{
  NO_VAR                    = 0,
  AUTOREGRESSIVE_MEAN_VAR   = 1,
  AUTOREGRESSIVE_HEIGHT_VAR = 2,
  LOCAL_GAUSSIANS_VAR       = 3
};

enum ae_env_noise
{
  NO_NOISE  = 0,
  FRACTAL   = 1
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
  //LOG_LOADS     = 0x08
};


#endif // __AE_ENUMS_H__
