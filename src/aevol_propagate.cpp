//*****************************************************************************
//
//                         aevol - Artificial Evolution
//
// Copyright (C) 2004 LIRIS.
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





int main( int argc, char* argv[] )
{
  // 1) Initialize command-line option variables with default values
  int32_t num_gener = -1;
  char* input_dir   = NULL;
  char* output_dir  = NULL;
  bool  verbose     = false;
  
  // 2) Define allowed options
  const char * options_list = "g:hi:o:v";
  static struct option long_options_list[] = {
    { "gener",    required_argument,  NULL, 'g' },
    { "help",     no_argument,        NULL, 'h' },
    { "in",       required_argument,  NULL, 'i' },
    { "out",      required_argument,  NULL, 'o' },
    { "verbose",  no_argument,        NULL, 'v' },
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
  
  exp_manager->load( num_gener, false, verbose );
  
  exp_manager->save_copy( output_dir );
}







void print_help( char* prog_path ) 
{
  // Get the program file-name in prog_name (strip prog_path of the path)
  char* prog_name; // No new, it will point to somewhere inside prog_path
  if ( prog_name = strrchr( prog_path, '/' ) ) prog_name++;
  else prog_name = prog_path;
  
	printf( "******************************************************************************\n" );
	printf( "*                        aevol - Artificial Evolution                        *\n" );
	printf( "******************************************************************************\n" );
	printf( "Usage : %s -h\n", prog_name );
	printf( "   or : %s [-f param_file]\n", prog_name );
	printf( "  -h, --help       Display this screen\n" );
	printf( "  -g, --gener      Specify generation number (input)\n" );
	printf( "  -i, --in         Specify input directory (default \".\")\n" );
	printf( "  -o, --out        Specify output directory (default \"./output\")\n" );
	printf( "  -v, --verbose    Be verbose\n" );
	printf( "Create a fresh copy of the provided experiment.\n" );
}
