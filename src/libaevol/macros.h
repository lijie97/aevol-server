// Aevol - An in silico experimental evolution platform
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

#include <cinttypes>

#ifndef AEVOL_MACROS_H__
#define AEVOL_MACROS_H__

constexpr int8_t NB_BASE = 2; // WARNING :  A lot of stuff has been optimized for binary genomes
                              //            Changing the value of NB_BASE implies verifying the existing code
                              //            and make changes where necessary

// NB The following strings are not easily replaced with `constexpr
// const char*` because they are meant to be concatenated by the
// preprocessor.

// Backup directories and file name formats
// Experimental Setup
#define EXP_S_DIR                 "exp_setup"
#define EXP_S_FNAME_BASE          "exp_setup_%06" PRId64
#define EXP_S_FNAME_FORMAT        EXP_S_DIR "/" EXP_S_FNAME_BASE ".ae"
#define EXP_S_CONST_FNAME_BASE        "exp_setup_const"
#define EXP_S_CONST_FNAME_FORMAT      EXP_S_DIR "/" EXP_S_CONST_FNAME_BASE ".ae"
// Output Profile
#define OUT_P_DIR                 "output_profile"
#define OUT_P_FNAME_BASE          "output_profile"
#define OUT_P_FNAME_FORMAT        OUT_P_DIR "/" OUT_P_FNAME_BASE ".ae"
#define OUT_P_CUR_FNAME           "output_profile.ae"
// Spatial Structure
#define WORLD_DIR             "world"
#define WORLD_FNAME_BASE      "world_%06" PRId64
#define WORLD_FNAME_FORMAT    WORLD_DIR "/" WORLD_FNAME_BASE".ae"
// Stats
#define STATS_DIR   "stats"
// Tree
#define TREE_DIR    "tree"
// Last gener file
#define LAST_GENER_FNAME "last_gener.txt"
// Best last organism file
#define BEST_LAST_ORG_FNAME "best_last_org.txt"

#define FIXED_POPULATION_SIZE // Some calculation can be spared if we know that the size of the population is fixed

constexpr int8_t PROM_SIZE = 22;
constexpr const char* PROM_SEQ = "0101011001110010010110";

constexpr int8_t PROM_MAX_DIFF  = 4;
constexpr int8_t TERM_STEM_SIZE = 4;
constexpr int8_t TERM_LOOP_SIZE = 3;
constexpr int8_t TERM_SIZE      = 2 * TERM_STEM_SIZE + TERM_LOOP_SIZE;

constexpr int8_t SHINE_DAL_SIZE = 6;
constexpr const char* SHINE_DAL_SEQ = "011011";
constexpr int8_t SHINE_START_SPACER = 4;

constexpr int8_t CODON_SIZE  = 3;
constexpr int8_t CODON_START = 0b000;
constexpr int8_t CODON_STOP  = 0b001;
constexpr int8_t CODON_M0    = 0b100;
constexpr int8_t CODON_M1    = 0b101;
constexpr int8_t CODON_W0    = 0b010;
constexpr int8_t CODON_W1    = 0b011;
constexpr int8_t CODON_H0    = 0b110;
constexpr int8_t CODON_H1    = 0b111;

#ifdef __REGUL
constexpr int8_t MAX_CODON   = 1 << CODON_SIZE;
constexpr int8_t QUADON_SIZE = 4;
constexpr int8_t MAX_QUADON  = 1 << QUADON_SIZE;
#endif

constexpr double X_MIN = 0.0;
constexpr double X_MAX = 1.0;
constexpr double Y_MIN = 0.0;
constexpr double Y_MAX = 1.0;
constexpr double H_MIN = -1.0;
constexpr double H_MAX = 1.0;
constexpr double W_MIN = 0.0;
// W_MAX is defined through a parameter

constexpr int8_t SC_MATCH_BONUS   = 1;
constexpr int8_t SC_MISMATCH_COST = 2;

#endif // AEVOL_MACROS_H__
