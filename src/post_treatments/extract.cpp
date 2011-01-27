// This program extracts some data about the individuals and write
// them into text files easy to parse with matlab.
//
// Two kinds of data can be extracted :
//
//  * data about the phenotype (option -p) : write information about
//    the proteins in a text file. A space delimits two proteins, a
//    new line delimits two individuals. For each protein, the output
//    is "m_h_w_r_f_l_z" where : 
//        * m, h and w are the mean, height and width of the protein
//        * r is an identifier of the rna it belongs (usefull to
//            know if several proteins are on the same rna)
//        * f and l are the first and last translated base
//        * z indicates the feature (at the center of the protein)
//    At the end of each line, "@n" is printed, where n is the number
//    of bases of the individual.
//
//  * sequence of the individuals (option -s) : write the sequences 
//    in a text file. A new line delimits two individuals. It can be
//    used to do sequence analysis in matlab.
//
// The option -b only treats the best individual
//
// Examples :
//
// For generation 20000, write infos about the phenotypes of all the
// individuals in phe_020000 and the sequences of all the
// individuals in seq_020000 :
//
//    extract -f backup/gen_020000.ae -p phe_020000 -s seq_020000
//
// For generation 20000, write the best individual's sequence in
// seq_020000_best :
//
//    extract -b -f backup/gen_020000.ae -s seq_020000_best
//


/** \class
 *  \brief
 */
 
 
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
#include <ae_common.h>
#include <ae_population.h>
#include <ae_individual.h>
#include <ae_list.h>
#include <ae_simulation.h>





void print_help( char* prog_name );
void analyse_indiv( ae_individual* indiv, FILE* phenotype_file, FILE* sequence_file );




int main( int argc, char* argv[] )
{
  // Initialize command-line option variables with default values
  char* backup_file_name  = NULL;
  char* phenotype_file_name  = NULL;
  char* sequence_file_name  = NULL;
  bool best_only = false;
  
  // Define allowed options
  const char * options_list = "hf:p:s:b";
  static struct option long_options_list[] = {
    { "file", 1, NULL, 'f' },
    { "phenotype", 1, NULL, 'p' },
    { "sequence", 1, NULL, 's' },
    { "best", 0, NULL, 'b' },
    { 0, 0, 0, 0 }
  };

  // Get actual values of the command-line options
  int option;
  while ( ( option = getopt_long(argc, argv, options_list, long_options_list, NULL) ) != -1 ) 
  {
          switch ( option ) 
            {
            case 'h' :
              print_help( argv[0] );
              exit( EXIT_SUCCESS );
              break;
            case 'f' :
              backup_file_name = new char[strlen(optarg) + 1];
              sprintf( backup_file_name, "%s", optarg );
              break;
            case 'p' :
              phenotype_file_name = new char[strlen(optarg) + 1];
              sprintf( phenotype_file_name, "%s", optarg );
              break;
            case 's' :
              sequence_file_name = new char[strlen(optarg) + 1];
              sprintf( sequence_file_name, "%s", optarg );
              break;
            case 'b' :
              best_only = true;
              break;
            }
  }
  
  // Open the files
  FILE* phenotype_file = NULL;
  FILE* sequence_file = NULL;
  
  if ( backup_file_name == NULL )
  {
    printf("You must specify a backup file. Please use the option -f or --file.\n");
    exit(EXIT_FAILURE);
  }
  if ( phenotype_file_name != NULL )
  {
    phenotype_file = fopen(phenotype_file_name,"w");
  }
  if ( sequence_file_name != NULL )
  {
    sequence_file = fopen(sequence_file_name,"w");
  }

  
  printf( "Reading backup file <%s>... \n", backup_file_name );
  fflush(stdout);

  // Load the simulation from backup
  ae_common::sim = new ae_simulation( backup_file_name, false );
  printf("done\n");
  
  printf( "Computing phenotypes... \n" );
  fflush( stdout );

  // Evaluate the individuals
  (ae_common::sim->get_pop())->evaluate_individuals(ae_common::sim->get_env());
  
  int i = 0;
  int nb_indiv = (ae_common::sim->get_pop())->get_nb_indivs();

  
  // Parse the individuals
  if (best_only)
  {
    ae_individual* best = ae_common::sim->get_pop()->get_best();
    analyse_indiv(best, phenotype_file, sequence_file);
  }
  else
  {
    if (ae_common::pop_structure)
    {
      ae_grid_cell*** _pop_grid = (ae_common::sim->get_pop())->get_pop_grid();
      for ( int16_t x = 0 ; x < ae_common::grid_x ; x++ )
      {
        for ( int16_t y = 0 ; y < ae_common::grid_y ; y++ )
        {
          ae_individual* indiv = (_pop_grid[x][y]->get_individual());
          analyse_indiv(indiv, phenotype_file, sequence_file);
          i++;
        }  
      }
    }
    else
    {
      ae_list_node  * indiv_node = ((ae_common::sim->get_pop())->get_indivs())->get_first();
      ae_individual * indiv      = NULL;

      while( indiv_node != NULL )
      {
        indiv = (ae_individual *) indiv_node->get_obj();
        analyse_indiv(indiv, phenotype_file, sequence_file);
        indiv_node = indiv_node->get_next();
        i++;
      }
    }
  }

  if (sequence_file_name != NULL)
  {
    fclose(sequence_file);
  }
  if (phenotype_file_name != NULL)
  {
    fclose(phenotype_file);
  }
  if (backup_file_name != NULL) {delete [] backup_file_name;}
  if (phenotype_file_name != NULL) {delete [] phenotype_file_name;}
  if (sequence_file_name != NULL) {delete [] sequence_file_name;}
  
  delete ae_common::sim;

  return EXIT_SUCCESS;
}

// The export fonction
inline void analyse_indiv( ae_individual* indiv, FILE* phenotype_file , FILE* sequence_file)
{
  ae_list* lrnas = indiv->get_rna_list();
  ae_list_node* rna_node = lrnas->get_first();
  ae_rna* rna = NULL;
  int32_t nb_bases =  indiv->get_genetic_unit(0)->get_dna()->get_length();
  int rna_nb=0;
  while( rna_node != NULL )
  {
    rna = (ae_rna *) rna_node->get_obj();
    
    ae_list* lprot = rna->get_transcribed_proteins();
    ae_list_node* prot_node = lprot->get_first();
    ae_protein* prot = NULL;
    while( prot_node != NULL )
    {
      prot = (ae_protein*) prot_node->get_obj();
      
      double height = prot->get_height();
      double width = prot->get_width();
      double mean = prot->get_mean();
      int32_t fpos = prot->get_first_translated_pos();
      int32_t lpos = prot->get_last_translated_pos();
      int32_t zone = 0;
      
      if ( mean != 1 )
      {
        zone = (int32_t) floor( mean * ae_common::env_axis_nb_segments );
      }
      else
      {
        zone = (int32_t) floor( mean * ae_common::env_axis_nb_segments ) - 1;
      }
      
      ae_env_axis_feature feat = ae_common::env_axis_features[zone];
      int nfeat = -1;
      switch ( feat )
      {
        case NEUTRAL :
          nfeat=1;
          break;
        case METABOLISM :
          nfeat=2;
          break;
        case SECRETION :
          nfeat=3;
          break;
        case TRANSFER :
          nfeat=4;
          break;
        case NB_FEATURES : 
          nfeat=-1;
          break;
      }

      if (phenotype_file!=NULL)
      {
        fprintf(phenotype_file,"%f_%f_%f_%d_%i_%i_%d ",mean,height,width,rna_nb,fpos,lpos,nfeat);
      }
      prot_node = prot_node->get_next();
    }
    
    rna_node = rna_node->get_next();
    rna_nb++;
  }
  if (phenotype_file!=NULL)
  {
    fprintf(phenotype_file,"@%d\n",nb_bases);
  }
  if (sequence_file!=NULL)
  {
    const char* dna = indiv->get_genetic_unit(0)->get_dna()->get_data();
    int32_t length = indiv->get_genetic_unit(0)->get_dna()->get_length();
    fprintf(sequence_file,"%.*s\n",length,dna);
  }
}

void print_help( char* prog_name ) 
{
  printf( "\n\
Usage : extract -h\n\
   or : extract -f source [-p phenotype_file] [-s sequence_file] [-b] \n\n\
\t-h : display this screen\n\
\t--file source : read from the backup file source\n\
\t--phenotype pf : extract and save some infos about the phenotypes of the individuals to file pf\
\t--sequence sf : extract and save the sequences of the individuals to file sf\
\t--best : only treat the best individual\n");
}
