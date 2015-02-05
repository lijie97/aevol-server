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



//
// This program extracts some data about the individuals and write
// them into text files easy to parse with e.g. matlab.
//
// Two kinds of data can be extracted :
//
//  * data about the phenotype (option -t) : write information about
//    the proteins in a text file. A space delimits two proteins, a
//    new line delimits two individuals. For each protein, the output
//    is "m_h_w_c_r_s_f_l_z_g" where :
//        * m, h, w and c are the mean, height, width and concentration of the
//            protein
//        * r is an identifier of the rna it belongs (usefull to
//            know if several proteins are on the same rna)
//        * s indicates the strand (LEADING/LAGGING)
//        * f and l are the first and last translated base
//        * z indicates the feature (at the center of the protein)
//        * g indicates the genetic unit to which the protein belongs
//            (0=chromosome, 1=plasmid)
//
//  * sequences of the individuals (option -s) : write the sequences
//    in a text file. A new line delimits two individuals. In case
//    there are several GU, they are separated by spaces.
//
// The option -b only treats the best individual
//
// The input can be either a generation number, in which case we
// will attempt to load a full backup tree, or a population file,
// in which case features of the proteins won't be outputed as we
// need to know the environment to infer them.
//
// Examples :
//
// For generation 20000, write infos about the phenotypes of all the
// individuals in phe_020000 and the sequences of all the
// individuals in seq_020000 :
//
//    extract -r 20000 -t phe_020000 -s seq_020000
//
// For generation 20000, write the best individual's sequence in
// seq_020000_best :
//
//    extract -b -r 20000 -s seq_020000_best
// or extract -b -p populations/pop_020000.ae -s seq_020000_best
//




// =================================================================
//                              Libraries
// =================================================================
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <signal.h>


// =================================================================
//                            Project Files
// =================================================================
#include <ae_population.h>
#include <ae_individual.h>
#include <environment.h>
#include <ae_list.h>
#include <ae_exp_manager.h>

using namespace aevol;

// =================================================================
//                         Function declarations
// =================================================================
void print_help(char* prog_path);
void print_version( void );

void analyse_indiv( ae_individual* indiv, FILE* triangles_file, FILE* sequence_file, int16_t gu, Environment* env );
void analyse_gu( ae_genetic_unit* gen_unit, int32_t gen_unit_number, FILE* triangles_file, Environment* env );




int main( int argc, char* argv[] )
{
  // Initialize command-line option variables with default values
  char* pop_file_name  = NULL;
  char* triangles_file_name  = NULL;
  char* sequence_file_name  = NULL;
  bool best_only = false;
  int16_t gu = -1;
  int32_t num_gener = -1;

  // Define allowed options
  const char * options_list = "hVp:r:t:s:bg:";
  static struct option long_options_list[] = {
    {"help",      no_argument,        NULL, 'h'},
    {"version",   no_argument,        NULL, 'V'},
    {"popfile",   required_argument,  NULL, 'p'},
    {"resume",    required_argument,  NULL, 'r'},
    {"triangles", required_argument,  NULL, 't'},
    {"sequence",  required_argument,  NULL, 's'},
    {"best",      no_argument,        NULL, 'b'},
    {"gu",        required_argument,  NULL, 'g'},
    {0, 0, 0, 0}
  };

  // Get actual values of the command-line options
  int option;
  while ( ( option = getopt_long(argc, argv, options_list, long_options_list, NULL) ) != -1 )
  {
    switch ( option )
    {
      case 'h' :
      {
        print_help(argv[0]);
        exit( EXIT_SUCCESS );
      }
      case 'V' :
      {
        print_version();
        exit( EXIT_SUCCESS );
      }
      case 'p' :
        pop_file_name = new char[strlen(optarg) + 1];
        sprintf( pop_file_name, "%s", optarg );
        break;
      case 'r':
        num_gener = atol( optarg );
        break;
      case 't' :
        triangles_file_name = new char[strlen(optarg) + 1];
        sprintf( triangles_file_name, "%s", optarg );
        break;
      case 's' :
        sequence_file_name = new char[strlen(optarg) + 1];
        sprintf( sequence_file_name, "%s", optarg );
        break;
      case 'g' :
        gu = atoi( optarg );
        break;
      case 'b' :
        best_only = true;
        break;
    }
  }

  // Open the files
  FILE* triangles_file = NULL;
  FILE* sequence_file = NULL;

  if ( triangles_file_name != NULL )
  {
    triangles_file = fopen(triangles_file_name,"w");
  }
  if ( sequence_file_name != NULL )
  {
    sequence_file = fopen(sequence_file_name,"w");
  }

  ae_population* pop = NULL;
  Environment* env = NULL;
  ae_exp_manager* exp_manager = new ae_exp_manager();

  // Two possible sources: either the user provided a "full" simulation via a generation number (option '-r'), either he just provided a population file (option '-p').
  if ( num_gener != -1 )
  {
    exp_manager->load( num_gener, false, false, false );
    pop = exp_manager->get_pop();
    env = exp_manager->get_env();
  }
  else
  {
    if ( pop_file_name == NULL )
    {
      printf("You must specify either a generation number or a source population file");
      exit(EXIT_FAILURE);
    }

    // Load the simulation from population file
    pop = new ae_population(exp_manager);

    gzFile pop_file = gzopen( pop_file_name, "r" );
    if ( pop_file == Z_NULL )
    {
      printf( "%s:%d: error: could not open backup file %s\n", __FILE__, __LINE__, pop_file_name );
      exit( EXIT_FAILURE );
    }
    pop->load( pop_file, false );
    gzclose( pop_file );
  }

  // The best individual is already known because it is the last in the list
  // Thus we do not need to know anything about the environment and to evaluate the individuals

  // Parse the individuals
  if (best_only)
  {
    ae_individual* best = pop->get_best();
    best->do_transcription_translation_folding(); // We need to recompute proteins if not already done (ie if using a population file and not a full backup)
    analyse_indiv(best, triangles_file, sequence_file, gu, env);
  }
  else
  {
    ae_list_node<ae_individual*>* indiv_node = pop->get_indivs()->get_first();
    ae_individual* indiv      = NULL;
    while( indiv_node != NULL )
    {
      indiv = (ae_individual *) indiv_node->get_obj();
      indiv->do_transcription_translation_folding(); // We need to recompute proteins if not already done (ie if using a population file and not a full backup)
      analyse_indiv(indiv, triangles_file, sequence_file, gu, env);
      indiv_node = indiv_node->get_next();
    }
  }

  if (sequence_file_name != NULL)
  {
    fclose(sequence_file);
  }
  if (triangles_file_name != NULL)
  {
    fclose(triangles_file);
  }

  if (pop_file_name != NULL) {delete [] pop_file_name;}
  if (triangles_file_name != NULL) {delete [] triangles_file_name;}
  if (sequence_file_name != NULL) {delete [] sequence_file_name;}

  delete exp_manager;
  if ((num_gener == -1)&&(pop!=NULL)) {delete pop;}

  return EXIT_SUCCESS;
}

// Parsing an individual
inline void analyse_indiv( ae_individual* indiv, FILE* triangles_file, FILE* sequence_file, int16_t gu, Environment* env )
{
  if ( gu == -1 ) // We want to treat all genetic units
  {
    int32_t gen_unit_number = 0;
    ae_list_node<ae_genetic_unit*>* gen_unit_node = indiv->get_genetic_unit_list()->get_first();
    while ( gen_unit_node != NULL )
    {
      ae_genetic_unit* gen_unit = (ae_genetic_unit*) gen_unit_node->get_obj();

      if ( triangles_file != NULL )
      {
        analyse_gu(gen_unit, gen_unit_number, triangles_file, env); // We call the triangle parser for each GU successively
      }
      if ( sequence_file != NULL )
      {
        const char* dna = gen_unit->get_dna()->get_data();
        int32_t length = gen_unit->get_dna()->get_length();
        fprintf(sequence_file,"%.*s ",length,dna); // We output the sequences of each GU separated by a space
      }

      gen_unit_node = gen_unit_node->get_next();
      gen_unit_number++;
    }
  }
  else // User specified a genetic unit
  {
    ae_genetic_unit* gen_unit = indiv->get_genetic_unit(gu);
    if ( triangles_file != NULL )
    {
      analyse_gu(gen_unit, gu, triangles_file, env); // We call the triangle parser
    }
    if ( sequence_file != NULL )
    {
      const char* dna = gen_unit->get_dna()->get_data();
      int32_t length = gen_unit->get_dna()->get_length();
      fprintf(sequence_file,"%.*s",length,dna); // We output the sequence
    }
  }

  // We go to next line in each file
  if ( triangles_file != NULL )
  {
    fprintf(triangles_file,"\n");
  }
  if ( sequence_file != NULL )
  {
    fprintf(sequence_file,"\n");
  }
}

// Parsing a GU
inline void analyse_gu( ae_genetic_unit* gen_unit, int32_t gen_unit_number, FILE* triangles_file, Environment* env )
{
  // Construct the list of all rnas
  ae_list<ae_rna*>** llrnas = gen_unit->get_rna_list();
  ae_list<ae_rna*>* lrnas = new ae_list<ae_rna*>();
  lrnas->add_list(llrnas[LEADING]);
  lrnas->add_list(llrnas[LAGGING]);

  // Parse this list
  ae_list_node<ae_rna*>* rna_node = lrnas->get_first();;

  ae_rna* rna = NULL;
  int rna_nb = 0;

  while( rna_node != NULL )
  {
    rna = (ae_rna *) rna_node->get_obj();

    ae_list<ae_protein*>* lprot = rna->get_transcribed_proteins();
    ae_list_node<ae_protein*>* prot_node = lprot->get_first();
    ae_protein* prot = NULL;

    while( prot_node != NULL )
    {
      prot = (ae_protein*) prot_node->get_obj();

      double height = prot->get_height();
      double width = prot->get_width();
      double mean = prot->get_mean();
      double concentration=rna->get_basal_level();
      int32_t fpos = prot->get_first_translated_pos();
      int32_t lpos = prot->get_last_translated_pos();

      int nfeat = -1;
      // Retrieving the feature of the protein also necessitates the an environment file.
      if ( env != NULL )
      {
        for (size_t i = 0; i <= env->get_nb_segments() - 1; ++i)
        {
          if ( (mean > env->get_segment_boundaries(i) ) && (mean < env->get_segment_boundaries(i+1)) )
          {
            nfeat = env->get_axis_feature(i);
            break;
          }
        }
      }

      fprintf(triangles_file,"%f_%f_%f_%f_%d_%d_%i_%i_%d_%d ",mean,height,width,concentration,rna_nb,rna->get_strand(),fpos,lpos,nfeat,gen_unit_number);

      prot_node = prot_node->get_next();
    }

    rna_node = rna_node->get_next();
    rna_nb++;
  }

  delete lrnas;
}


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
  printf("%s: extracts the genotype and/or data about the phenotype of individuals in the provided population and write them into text files easy to parse with e.g. matlab.\n", prog_name);
  printf("\n");
  printf("Usage : %s -h\n", prog_name);
  printf("   or : %s -V or --version\n", prog_name);
  printf("   or :    %s [-r GENER | -p POP_FILE] [-t PHEN_FILE] [-s SEQ_FILE] [-g NUM_GU] [-b]\n", prog_name);
  printf("\nOptions\n");
  printf("  -h, --help\n\tprint this help, then exit\n\n");
  printf("  -V, --version\n\tprint version number, then exit\n\n");
  printf("  -r GENER  :\n");
  printf("\tread generation GENER from a full aevol backup\n");
  printf("  -p POP_FILE:\n");
  printf("\tread the population saved in population file POP_FILE\n");
  printf("  -t PHEN_FILE:\n");
  printf("\textract and save some infos about the phenotypes of the individuals to file PHEN_FILE\n");
  printf("  -s SEQ_FILE:\n");
  printf("\textract and save the sequences of the individuals to file SEQ_FILE\n");
  printf("  -g NUM_GU:\n");
  printf("\tonly treat this genetic unit (by default: treat all genetic units)\n");
  printf("  -b:\n");
  printf("\tonly treat the best individual\n");
  printf("\n\n");

  printf("\n\
This program extracts some data about the individuals and write\n\
them into text files easy to parse with e.g. matlab.\n\
\n\
Two kinds of data can be extracted :\n\
\n\
 * data about the phenotype (option -t) : write information about\n\
   the proteins in a text file. A space delimits two proteins, a\n\
   new line delimits two individuals. For each protein, the output\n\
   is \"m_h_w_c_r_s_f_l_z_g\" where :\n\
       * m, h, w and c are the mean, height, width and concentration of the protein\n\
       * r is an identifier of the rna it belongs (useful to\n\
           know if several proteins are on the same rna)\n\
       * s indicates the strand (LEADING/LAGGING)\n\
       * f and l are the first and last translated base\n\
       * z indicates the feature (at the center of the protein)\n\
       * g indicates the genetic unit to which the protein belongs (0=chromosome, 1=plasmid)\n\
\n\
 * sequences of the individuals (option -s) : write the sequences\n\
   in a text file. A new line delimits two individuals. In case\n\
   there are several GU, they are separated by whitespaces.\n\
\n\
With option -b, only the best individual is treated.\n\
\n\
The input can be either a generation number, in which case we\n\
will attempt to load a full backup tree, or a population file,\n\
in which case features of the proteins won't be outputed as we\n\
need to know the environment to infer them.\n\
\n\
Examples :\n\
\n\
For generation 20000, write infos about the phenotypes of all the\n\
individuals in phe_020000 and the sequences of all the\n\
individuals in seq_020000 :\n\
\n\
   extract -r 20000 -t phe_020000 -s seq_020000\n\
\n\
For generation 20000, write the best individual's sequence in\n\
seq_020000_best :\n\
\n\
   extract -b -r 20000 -s seq_020000_best\n\
or extract -b -p populations/pop_020000.ae -s seq_020000_best\n");

}


/*!
  \brief Print aevol version number

*/
void print_version( void )
{
  printf( "aevol %s\n", VERSION );
}
