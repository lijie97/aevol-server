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
#include <errno.h>
#include <inttypes.h>
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <signal.h>
#include <sys/stat.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// =================================================================
//                            Project Files
// =================================================================
#ifdef __X11
  #include "ExpManager_X11.h"
#else
  #include "ExpManager.h"
#endif

#include "macros.h"

using namespace aevol;

#ifndef __NO_X
void catch_usr1(int sig_num);
#endif

void print_help(char* prog_path);

static ExpManager* exp_manager = NULL;




int main(int argc, char* argv[])
{
  // Set handlers for signals to be caught
  #ifndef __NO_X
    signal(SIGUSR1, catch_usr1);
  #endif
  
  
  // Print warning for debug mode
  #ifdef DEBUG
    printf("aevol is being run in DEBUG mode\n");
  #endif
  
 
  
  // =========================================================================
  //                           Get command-line options
  // =========================================================================
  // 1) Initialize command-line option variables with default values
  // 2) Define allowed options
  // 3) Get actual values of the command-line options
  // 4) Check the consistancy of the command-line options
  // 5) Set file names according to options
  // =========================================================================
  
  // -------------------------------------------------------------------------
  // 1) Initialize command-line option variables with default values
  // -------------------------------------------------------------------------
  // bool pause_on_startup = false;
  bool verbose          = false;
  
  int64_t t0 = -1;
  int64_t t_end = -1;
  int64_t nb_steps = -1;
  
  #ifndef __NO_X
    bool show_display_on_startup = true;
  #endif

  bool run_in_parallel = false;
  
  
  // -------------------------------------------------------------------------
  // 2) Define allowed options
  // -------------------------------------------------------------------------
  const char * options_list = "he:n:r:vVwxp:";
  static struct option long_options_list[] = {
    // Print help
    { "help",     no_argument,        NULL, 'h' },
    // time up to which to simulate (use either -e or -n)
    { "end",  required_argument,  NULL, 'e' },
    // Number of generations to be run (use either -e or -n)
    { "nsteps",  required_argument,  NULL, 'n' },
    // Resume from generation X (default: 0)
    { "resume",   required_argument,  NULL, 'r' },
    // Be verbose
    { "verbose",  no_argument,        NULL, 'v' },
    // Print version
    { "version",  no_argument,        NULL, 'V' },
    // Pause after loading
    { "wait",     no_argument,        NULL, 'w' },
    // Don't display X outputs on start
    { "noX",      no_argument,        NULL, 'x' },
    // Run in parallel on x threads (0 or negative value yields system default)
    { "parallel", required_argument,  NULL, 'p' },
    { 0, 0, 0, 0 }
  };
  
  
  // -------------------------------------------------------------------------
  // 3) Get actual values of the command-line options
  // -------------------------------------------------------------------------
  int option;
  while ((option =
              getopt_long(argc, argv, options_list, long_options_list, NULL))
         != -1) {
    switch (option) {
      case 'h' : {
        print_help(argv[0]);
        exit(EXIT_SUCCESS);
      }
      case 'V' : {
        Utils::PrintAevolVersion();
        exit(EXIT_SUCCESS);
      }
      case 'e' : {
        if (nb_steps != -1)
          Utils::ExitWithUsrMsg("use either option -n or -e, not both");

        t_end = atol(optarg);
        break;
      }
      case 'n' : {
        if (t_end != -1)
          Utils::ExitWithUsrMsg("use either option -n or -e, not both");

        nb_steps = atol(optarg);
        break;
      }
      case 'r' : {
        t0 = atol(optarg);
        break;      
      }
      case 'v' : {
        verbose = true;
        break;
      }
      case 'w' : {
        // pause_on_startup = true;
        break;
      }
      case 'x' : {
        #ifdef __NO_X
          printf("%s: error: Program was compiled with __NO_X option, "
                 "no visualisation available.\n", argv[0]);
          exit(EXIT_FAILURE);
        #else
          show_display_on_startup = false;
        #endif
        
        break;
      }
      case 'p' : {
        #ifdef _OPENMP
        run_in_parallel = true;
        if (atoi(optarg) > 0)
          omp_set_num_threads(atoi(optarg));
        #endif
        break;
      }
      default : {
        // An error message is printed in getopt_long, we just need to exit
        exit(EXIT_FAILURE);
      }
    }
  }

  // If t0 wasn't provided, use default
  if (t0 < 0)
    t0 = OutputManager::get_last_gener();

  // If t_end wasn't provided, set it according to nb_steps or use default (run
  // for 1000 timesteps)
  if (t_end < 0)
    if (nb_steps >= 0)
      t_end = t0 + nb_steps;
    else
      t_end = t0 + 1000;

  // It the user didn't ask for a parallel run, set number of threads to 1
  #ifdef _OPENMP
  if (not run_in_parallel)
    omp_set_num_threads(1);
  #endif
  
  // =================================================================
  //                          Load the simulation
  // =================================================================
  #ifndef __NO_X
    exp_manager = new ExpManager_X11();
  #else
    exp_manager = new ExpManager();
  #endif
  
  exp_manager->load(t0, verbose, true);
  exp_manager->set_t_end(t_end);
  
 

  // Make a numbered copy of each static input file (dynamic files are saved elsewhere)
  // TODO (?)

  #ifndef __NO_X
    if (show_display_on_startup)
    {
      ((ExpManager_X11 *) exp_manager)->toggle_display_on_off();
    }
  #endif
  
  
  // =================================================================
  //                         Run the simulation
  // =================================================================
  exp_manager->run_evolution();
  
  
  delete exp_manager;
  //~ return EXIT_SUCCESS;
}


#ifndef __NO_X
void catch_usr1(int sig_num)
{
  signal(SIGUSR1, catch_usr1);
  
  printf("display on/off\n");
  if (exp_manager != NULL)
  {
    ((ExpManager_X11 *) exp_manager)->toggle_display_on_off();
  }
}
#endif







/*!
  \brief 
  
*/
void print_help(char* prog_path) 
{
  // Get the program file-name in prog_name (strip prog_path of the path)
  char* prog_name; // No new, it will point to somewhere inside prog_path
  if ((prog_name = strrchr(prog_path, '/'))) prog_name++;
  else prog_name = prog_path;
  
	printf("******************************************************************************\n");
	printf("*                                                                            *\n");
	printf("*                        aevol - Artificial Evolution                        *\n");
	printf("*                                                                            *\n");
	printf("* Aevol is a simulation platform that allows one to let populations of       *\n");
  printf("* digital organisms evolve in different conditions and study experimentally  *\n");
  printf("* the mechanisms responsible for the structuration of the genome and the     *\n");
  printf("* transcriptome.                                                             *\n");
	printf("*                                                                            *\n");
	printf("******************************************************************************\n");
  printf("\n");
	printf("%s: run an aevol simulation.\n", prog_name);
  printf("\n");
	printf("Usage : %s -h or --help\n", prog_name);
	printf("   or : %s -V or --version\n", prog_name);
	printf("   or : %s [-r TIME] [-e TIME|-n NB_TIMESTEPS] [-tvwx]\n", prog_name);
	printf("\nOptions\n");
	printf("  -h, --help\n\tprint this help, then exit\n\n");
	printf("  -V, --version\n\tprint version number, then exit\n\n");
  printf("  -r, --resume TIME\n\tspecify time t0 to resume simulation at (default value read in last_gener.txt)\n\n");
  printf("  -e, --end TIME\n\tspecify time of the end of the simulation\n\n");
  printf("  -n, --nsteps NB_TIMESTEPS\n\tspecify number of timesteps to be simulated (default 1000)\n\n");
  printf("  -t, --text\n\tuse text files instead of binary files when possible\n\n");
	printf("  -v, --verbose\n\tbe verbose\n\n");
  printf("  -w, --wait\n\tpause after loading\n\n");
  printf("  -x, --noX\n\tdon't display X outputs upon start\n\tsend SIGUSR1 to switch X output on/off\n");
}
