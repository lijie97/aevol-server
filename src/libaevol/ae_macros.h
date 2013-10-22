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


#include <inttypes.h>

#define NB_BASE INT8_C(2) // WARNING :  A lot of stuff has been optimized for binary genomes
                          //            Changing the value of NB_BASE implies verifying the existing code
                          //            and make changes where necessary

// Backup directories and file name formats
// Experimental Setup
#define EXP_S_DIR                 "exp_setup"
#define EXP_S_FNAME_BASE          "exp_setup_%06"PRId32
#define EXP_S_FNAME_FORMAT        EXP_S_DIR"/"EXP_S_FNAME_BASE".ae"
//~ #define EXP_S_TXT_FNAME_FORMAT    EXP_S_DIR"/"EXP_S_FNAME_BASE".in"
#define EXP_S_CONST_FNAME_BASE        "exp_setup_const"
#define EXP_S_CONST_FNAME_FORMAT      EXP_S_DIR"/"EXP_S_CONST_FNAME_BASE".ae"
#define EXP_S_TXT_CONST_FNAME_FORMAT  EXP_S_DIR"/"EXP_S_CONST_FNAME_BASE".in"
// Output Profile
#define OUT_P_DIR                 "output_profile"
#define OUT_P_FNAME_BASE          "output_profile"
#define OUT_P_FNAME_FORMAT        OUT_P_DIR"/"OUT_P_FNAME_BASE".ae"
#define OUT_P_TXT_FNAME_FORMAT    OUT_P_DIR"/"OUT_P_FNAME_BASE".in"
#define OUT_P_CUR_FNAME           "output_profile.ae"
#define OUT_P_TXT_CUR_FNAME       "output_profile.in"
// Environment
#define ENV_DIR                   "environment"
#define ENV_FNAME_BASE            "env_%06"PRId32
#define ENV_FNAME_FORMAT          ENV_DIR"/"ENV_FNAME_BASE".ae"
// Population
#define POP_DIR                   "populations"
#define POP_FNAME_BASE            "pop_%06"PRId32
#define POP_FNAME_FORMAT          POP_DIR"/"POP_FNAME_BASE".ae"
// Spatial Structure
#define SP_STRUCT_DIR             "spatial_structure"
#define SP_STRUCT_FNAME_BASE      "sp_struct_%06"PRId32
#define SP_STRUCT_FNAME_FORMAT    SP_STRUCT_DIR"/"SP_STRUCT_FNAME_BASE".ae"
// Stats
#define STATS_DIR   "stats"
// Tree
#define TREE_DIR    "tree"
// Last gener file
#define LAST_GENER_FNAME "last_gener.txt"
// Best last organism file
#define BEST_LAST_ORG_FNAME "best_last_org.txt"


#define FIXED_POPULATION_SIZE // Some calculation can be spared if we know that the size of the population is fixed

#define PROM_SIZE       INT8_C(22)
#define PROM_SEQ        "0101011001110010010110"

#define PROM_MAX_DIFF   INT8_C(4)
#define TERM_STEM_SIZE  INT8_C(4)
#define TERM_LOOP_SIZE  INT8_C(3)
#define TERM_SIZE       INT8_C(2 * TERM_STEM_SIZE + TERM_LOOP_SIZE)

#define SHINE_DAL_SIZE      INT8_C(6)
#define SHINE_DAL_SEQ       "011011"
#define SHINE_START_SPACER  INT8_C(4)

#define CODON_SIZE    INT8_C(3)
#define CODON_START   INT8_C(0)   // "000"
#define CODON_STOP    INT8_C(1)   // "001"
#define CODON_M0      INT8_C(4)   // "100"
#define CODON_M1      INT8_C(5)   // "101"
#define CODON_W0      INT8_C(2)   // "010"
#define CODON_W1      INT8_C(3)   // "011"
#define CODON_H0      INT8_C(6)   // "110"
#define CODON_H1      INT8_C(7)   // "111"

#ifdef __REGUL
  #define MAX_CODON   INT8_C(1 << CODON_SIZE)
  #define QUADON_SIZE 4
  #define MAX_QUADON  INT8_C(1 << QUADON_SIZE)
#endif

#define X_MIN   0.0
#define X_MAX   1.0
#define Y_MIN   0.0
#define Y_MAX   1.0
#define H_MIN   -1.0
#define H_MAX   1.0
#define W_MIN   0.0
// W_MAX is defined through a parameter

#define SC_MATCH_BONUS    INT8_C(1)
#define SC_MISMATCH_COST  INT8_C(2)

