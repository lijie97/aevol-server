//*****************************************************************************
//
//                         aevol - Artificial Evolution
//
// Copyright (C) 2004 LIRIS.
// Web: https://liris.cnrs.fr/
// E-mail: carole.knibbe@liris.cnrs.fr
// Original Authors : Guillaume Beslon, Carole Knibbe, David Parsons
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

const char* DEFAULT_PARAM_FILE_NAME = "param.in";

// =================================================================
//                              Includes
// =================================================================
#include <cstdlib>
#include <cstdio>
#include <cstring>

#include <getopt.h>

#ifdef __X11
  #include "ExpManager_X11.h"
#else
  #include "ExpManager.h"
#endif

#include "ExpManager.h"
#include "ParamLoader.h"


using namespace aevol;

// =================================================================
//                         Function declarations
// =================================================================
void print_help(char* prog_path);





int main(int argc, char* argv[])
{
  // 1) Initialize command-line option variables with default values
  char* param_file_name = NULL;
  char* output_dir = NULL;
  char* chromosome_file_name = NULL;
  char* plasmid_file_name = NULL;


  // 2) Define allowed options
  const char * options_list = "hVf:o:c:p:";
  static struct option long_options_list[] = {
    { "help",     no_argument,        NULL, 'h' },
    { "version",  no_argument,        NULL, 'V' },
    { "file",     required_argument,  NULL, 'f' },
    { "out",      required_argument,  NULL, 'o' },
    { "chromosome",   required_argument,  NULL, 'c' },
    { "plasmid",   required_argument,  NULL, 'p' },
    { 0, 0, 0, 0 }
  };


  // 3) Get actual values of the command-line options
  int option;
  while ((option = getopt_long(argc, argv, options_list, long_options_list, NULL)) != -1)
  {
    switch (option)
    {
      case 'h' :
      {
        print_help(argv[0]);
        exit(EXIT_SUCCESS);
      }
      case 'V' :
      {
        Utils::PrintAevolVersion();
        exit(EXIT_SUCCESS);
      }
      case 'f' :
      {
        param_file_name = new char[strlen(optarg)+1];
        strcpy(param_file_name, optarg);
        break;
      }
      case 'o' :
      {
        output_dir = new char[strlen(optarg)+1];
        strcpy(output_dir, optarg);
        break;
      }
      case 'c':
      {
        chromosome_file_name = new char [strlen(optarg)+1];
        strcpy(chromosome_file_name, optarg);
        break;
      }
      case 'p':
      {
        plasmid_file_name = new char [strlen(optarg)+1];
        strcpy(plasmid_file_name, optarg);
        break;
      }
      default :
      {
        // An error message is printed in getopt_long, we just need to exit
        exit(EXIT_FAILURE);
      }
    }
  }


  // 4) Set undefined command line parameters to default values
  if (param_file_name == NULL)
  {
    param_file_name = new char[strlen(DEFAULT_PARAM_FILE_NAME)+1];
    sprintf(param_file_name, "%s", DEFAULT_PARAM_FILE_NAME);
  }


  // 5) Create a param loader for the parameter file
  ParamLoader * my_param_loader = new ParamLoader(param_file_name);
  delete param_file_name;


  // 6) Initialize the experiment manager
  ExpManager * exp_manager = new ExpManager();


  // 7) Initialize the simulation from the parameter file
  int32_t lchromosome = -1;
  char* chromosome;

  if (chromosome_file_name != NULL)
  {
    const int max_input_chrom_size = 1000000;
    char raw_chromosome[max_input_chrom_size];
    FILE* chromosome_file = fopen(chromosome_file_name,"r");
    if (chromosome_file==NULL)
    {
      printf("ERROR: failed to open source chromosome file %s\n",chromosome_file_name);
      exit(EXIT_FAILURE);
    }
    if (fgets(raw_chromosome, max_input_chrom_size, chromosome_file) == NULL)
    {
      printf("ERROR: failed to read from chromosome file %s\n",chromosome_file_name);
      exit(EXIT_FAILURE);
    }
    lchromosome = strlen(raw_chromosome)-1;
    chromosome = new char[lchromosome]; // Warning: will become the DNA of the first individual created -> no not delete, will be freed in Dna   strncpy(chromosome, raw_chromosome, lchromosome);
    printf("Loading chromosome from text file %s (%" PRId32 " base pairs) \n",chromosome_file_name,lchromosome);
    fclose(chromosome_file);
  }

  int32_t lplasmid = -1;
  char* plasmid;

  if (plasmid_file_name != NULL)
  {
    const int max_input_plasmid_size = 1000000;
    char raw_plasmid[max_input_plasmid_size];
    FILE* plasmid_file = fopen(plasmid_file_name,"r");
    if (plasmid_file==NULL)
    {
      printf("ERROR: failed to open source chromosome file %s\n",plasmid_file_name);
      exit(EXIT_FAILURE);
    }
    if (fgets(raw_plasmid, max_input_plasmid_size, plasmid_file) == NULL)
    {
      printf("ERROR: failed to read from chromosome file %s\n",chromosome_file_name);
      exit(EXIT_FAILURE);
    }
    lplasmid = strlen(raw_plasmid)-1;
    plasmid = new char[lplasmid]; // Warning: will become the DNA of the first individual created -> no not delete, will be freed in Dna   strncpy(plasmid, raw_plasmid, lplasmid);
    printf("Loading plasmid from text file %s (%" PRId32 " base pairs) \n",plasmid_file_name,lplasmid);
    fclose(plasmid_file);
  }

  if (lchromosome > -1)
  {
    if (lplasmid > -1)
    {
      my_param_loader->load(exp_manager, true, chromosome, lchromosome, plasmid, lplasmid);
    }
    else
    {
      my_param_loader->load(exp_manager, true, chromosome, lchromosome);
    }
  }
  else
  {
    my_param_loader->load(exp_manager, true);
  }
  delete my_param_loader;


  //~ ((ExpManager_X11*)exp_manager)->toggle_display_on_off();
  //~ exp_manager->display();
  //~ getchar();

  // 8) Save the experiment
  if (output_dir == NULL)
  {
    exp_manager->Save();
  }
  else
  {
    // Save everything in the provided directory
    exp_manager->save_copy(output_dir);

    delete output_dir;
  }

  delete exp_manager;
}







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
  printf("%s: create an experiment with setup as specified in param_file.\n", prog_name);
  printf("\n");
  printf("Usage : %s -h or --help\n", prog_name);
  printf("   or : %s -V or --version\n", prog_name);
  printf("   or : %s [-f PARAM_FILE] [-o OUTDIR] [-c CFILE] [-p PFILE]\n", prog_name);
  printf("\nOptions\n");
  printf("  -h, --help\n\tprint this help, then exit\n\n");
  printf("  -V, --version\n\tprint version number, then exit\n\n");
  printf("  -f, --file PARAM_FILE\n\tspecify parameter file (default: param.in)\n");
  printf("  -o, --out OUTDIR\n\tspecify output directory (default \"./\")\n\n");
  printf("  -c, --chromosome CFILE\n\tload chromosome from given text file instead of generating it\n");
  printf("  -p, --plasmid PFILE\n\tload plasmid from given text file instead of generating it\n");
}
