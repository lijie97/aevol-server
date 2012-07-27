#include <inttypes.h>

#define NB_BASE INT8_C(2) // WARNING :  A lot of stuff has been optimized for binary genomes
                          //            Changing the value of NB_BASE implies verifying the existing code
                          //            and make changes where necessary

// Backup directories and file name formats
#define ENV_DIR               "environment"
#define POP_DIR               "populations"
#define SP_STRUCT_DIR         "spatial_structure"
#define ENV_FNAME_BASE        "env_%06"PRId32
#define POP_FNAME_BASE        "pop_%06"PRId32
#define SP_STRUCT_FNAME_BASE  "sp_struct_%06"PRId32

#ifdef __REGUL
  #define ENV_FNAME_FORMAT        ENV_DIR"/"ENV_FNAME_BASE".rae"
  #define POP_FNAME_FORMAT        POP_DIR"/"POP_FNAME_BASE".rae"
  #define SP_STRUCT_FNAME_FORMAT  SP_STRUCT_DIR"/"SP_STRUCT_FNAME_BASE".rae"
#else
  #define ENV_FNAME_FORMAT        ENV_DIR"/"ENV_FNAME_BASE".ae"
  #define POP_FNAME_FORMAT        POP_DIR"/"POP_FNAME_BASE".ae"
  #define SP_STRUCT_FNAME_FORMAT  SP_STRUCT_DIR"/"SP_STRUCT_FNAME_BASE".ae"
#endif

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
