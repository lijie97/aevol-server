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

 
const char* DEFAULT_PARAM_FILE_NAME = "param.in";


// =================================================================
//                              Libraries
// =================================================================
#include <err.h>
#include <errno.h>
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>




// =================================================================
//                            Project Files
// =================================================================
#ifdef __X11
  #include <ae_exp_manager_X11.h>
#else
  #include <ae_exp_manager.h>
#endif
#include <param_loader.h>



// =================================================================
//                         Function declarations
// =================================================================
void print_help( char* prog_path );
void print_version( void );





int main( int argc, char* argv[] )
{
  // 1) Initialize command-line option variables with default values
  int32_t num_gener = -1;
  char* input_dir   = NULL;
  char* output_dir  = NULL;
  bool  verbose     = false;
  
  // 2) Define allowed options
  const char * options_list = "g:hi:o:vV";
  static struct option long_options_list[] = {
    { "gener",    required_argument,  NULL, 'g' },
    { "help",     no_argument,        NULL, 'h' },
    { "in",       required_argument,  NULL, 'i' },
    { "out",      required_argument,  NULL, 'o' },
    { "verbose",  no_argument,        NULL, 'v' },
    { "version",  no_argument,        NULL, 'V' },
    { 0, 0, 0, 0 }
  };
      
  // 3) Get actual values of the command-line options
  int option;
  while ( ( option = getopt_long(argc, argv, options_list, long_options_list, NULL) ) != -1 ) 
  {
    switch ( option ) 
    {
      case 'h' :
      {
        print_help( argv[0] );
        exit( EXIT_SUCCESS );
      }
      case 'V' :
      {
        print_version();
        exit( EXIT_SUCCESS );
      }
      case 'g' :
      {
        num_gener = atoi( optarg );
        break;
      }
      case 'i' :
      {
        input_dir = new char[strlen(optarg)+1];
        strcpy( input_dir, optarg );
        break;
      }
      case 'o' :
      {
        output_dir = new char[strlen(optarg)+1];
        strcpy( output_dir, optarg );
        break;
      }
      case 'v' :
      {
        verbose = true;
        break;
      }
      default :
      {
        // An error message is printed in getopt_long, we just need to exit
        exit( EXIT_FAILURE );
      }
    }
  }
  
  
  // 4) Set undefined command line parameters to default values
  if ( input_dir == NULL )
  {
    input_dir = new char[255];
    sprintf( input_dir, "%s", "." );
  }
  if ( output_dir == NULL )
  {
    output_dir = new char[255];
    sprintf( output_dir, "%s", "output" );
  }
  if ( num_gener == -1 )
  {
    // Set num_gener to the content of the LAST_GENER file if it exists.
    // If it doesn't, print help and exit
    FILE* lg_file = fopen( LAST_GENER_FNAME, "r" );
    if ( lg_file != NULL )
    {
      fscanf( lg_file, "%"PRId32"\n", &num_gener );
      fclose( lg_file );
    }
    else
    {
      printf( "aevol_propagate: no generation number provided.\n" );
      print_help( argv[0] );
      exit( EXIT_FAILURE );
    }
  }
  
  
  // 5) Check whether the output directory is missing
  struct stat stat_buf;
  if ( (stat( output_dir, &stat_buf ) == -1) && (errno == ENOENT) )
  {
    printf( "Directory \"%s\" does not exist. Create it ? [Y/n]\n", output_dir );
    char answer = getchar();
    while ( answer != 'y' and answer != 'n' and answer != '\n' )
    {
      printf( "Please answer by 'y' or 'n'. Create output directory ? [Y/n]\n" );
      while( answer != '\n' && answer != EOF) answer = getchar(); // "flush" stdin
      answer = getchar();
    }
    char flush = answer;
    while( flush != '\n' && flush != EOF) flush = getchar(); // "flush" stdin
    if ( answer == '\n' ) answer = 'y';
    
    if ( answer == 'n' ) exit( EXIT_SUCCESS );
    
    if ( mkdir( output_dir, 0755 ) )
    {
      err( EXIT_FAILURE, output_dir, errno );
    }
  }
  
  
  
  // =================================================================
  //                    Load the model experiment
  // =================================================================
  #ifndef __NO_X
    ae_exp_manager* exp_manager = new ae_exp_manager_X11();
  #else
    ae_exp_manager* exp_manager = new ae_exp_manager();
  #endif
  
  exp_manager->load( input_dir, num_gener, false, verbose, false );
  
  exp_manager->save_copy( output_dir );
}







/*!
  \brief 
  
*/
void print_help( char* prog_path ) 
{
  // Get the program file-name in prog_name (strip prog_path of the path)
  char* prog_name; // No new, it will point to somewhere inside prog_path
  if ( prog_name = strrchr( prog_path, '/' ) ) prog_name++;
  else prog_name = prog_path;
  
	printf( "******************************************************************************\n" );
	printf( "*                                                                            *\n" );
	printf( "*                        aevol - Artificial Evolution                        *\n" );
	printf( "*                                                                            *\n" );
	printf( "* Aevol is a simulation platform that allows one to let populations of       *\n" );
  printf( "* digital organisms evolve in different conditions and study experimentally  *\n" );
  printf( "* the mechanisms responsible for the structuration of the genome and the     *\n" );
  printf( "* transcriptome.                                                             *\n" );
	printf( "*                                                                            *\n" );
	printf( "******************************************************************************\n" );
  printf( "\n" );
	printf( "%s:\n", prog_name );
  printf( "\tCreate a fresh copy of the experiment as it was at the given generation.\n" );
  printf( "\tThe generation number of the copy will be reset to 0.\n" );
  printf( "\n" );
	printf( "Usage : %s -h or --help\n", prog_name );
	printf( "   or : %s -V or --version\n", prog_name );
	printf( "   or : %s [-g GENER] [-i in_dir] [-o out_dir] [-v]\n", prog_name );
	printf( "\nOptions\n" );
	printf( "  -h, --help\n\tprint this help, then exit\n\n" );
	printf( "  -V, --version\n\tprint version number, then exit\n\n" );
  printf( "  -g, --gener GENER\n\tspecify generation number\n" );
  printf( "\t(default: that contained in file last_gener.txt, if any)\n\n" );
	printf( "  -i, --in in_dir\n\tspecify input directory (default \".\")\n\n" );
	printf( "  -o, --out out_dir\n\tspecify output directory (default \"./output\")\n\n" );
	printf( "  -v, --verbose\n\tbe verbose\n" );
}

/*!
  \brief 
  
*/
void print_version( void ) 
{
	printf( "aevol %s\n", VERSION );
}
