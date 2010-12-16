//*****************************************************************************
//
//                         aevol - Artificial Evolution
//
// Copyright (C) 2009  LIRIS.
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
 
 
 
 
// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>
#include <sys/stat.h>  // for the permission symbols used with mkdir



// =================================================================
//                            Project Files
// =================================================================
#include <ae_common.h>
#include <ae_population.h>
#include <ae_individual.h>
#include <ae_environment.h>
#include <ae_protein.h>
#include <ae_rna.h>
#include <ae_list.h>
#include <ae_simulation.h>



// ========
//  TO DO 
// ========
// 
//  * option --color ?
//  * Raevol-specific output (EPS file with the network) ?





void print_help( void );


// The height of each triangle is proportional to the product c*m, where c is the 
// concentration of the protein and m its intrinsic efficiency (depending on its 
// aminoacid sequence). In the case of Raevol, the concentration used here is the 
// final one, i.e. the one reached after all the time steps of the lifetime. 
// If a coding sequence has several promoters, only one triangle is drawn.
void draw_triangles( ae_individual* indiv, ae_environment* env, char * directoryName );



// In the case of Raevol, the profile is drawn using the final concentrations
// of the proteins, i.e. the ones reached after all the time steps of the lifetime.
void draw_pos_neg_profiles( ae_individual * indiv, ae_environment* env, char * directoryName );



// In the case of Raevol, the phenotype is drawn using the final concentrations
// of the proteins, i.e. the ones reached after all the time steps of the lifetime.
void draw_phenotype( ae_individual * indiv, ae_environment * envir, char * directoryName );




// The chromosome is drawn as a circle. Coding sequences are represented by arcs starting 
// at the start codon and ending at the stop codon. Coding sequences that are on the 
// leading (resp. lagging) strand are drawn outside (resp. inside) the circle. If a coding 
// sequence has several promoters, only one arc is drawn.  
void draw_genetic_unit_with_CDS( ae_genetic_unit* gen_unit, char * directoryName );


// The chromosome is drawn as a circle. Transcribed sequences are represented by arcs starting 
// at the first trasncribed position codon and ending at the last transcribed position. mRNAs 
// that are on the leading (resp. lagging) strand are drawn outside (resp. inside) the circle. 
// mRNAs which include at least one coding sequence are black, the others are gray.
 void draw_genetic_unit_with_mRNAs( ae_genetic_unit* gen_unit, char * directoryName );


int main( int argc, char* argv[] )
{
  // =================================================================
  //                      Get command-line options
  // =================================================================
  //
  // 1) Initialize command-line option variables with default values

  bool  verbose               = false;
  bool  use_single_indiv_file = false;
  char* backup_file_name      = NULL;
  
  // 2) Define allowed options
  const char * options_list = "hvf:";
  static struct option long_options_list[] = {
    { "file", 1, NULL, 'f' },
    { 0, 0, 0, 0 }
  };

  // 3) Get actual values of the command-line options
  int option;
  while ( ( option = getopt_long(argc, argv, options_list, long_options_list, NULL) ) != -1 ) 
  {
    switch ( option ) 
    {
      case 'h' :
        print_help();
        exit( EXIT_SUCCESS );
        break;
      case 'v' :
        verbose = true;
        break;
      case 'f' :
        if ( strstr( optarg, "best" ) != NULL )
        {
          use_single_indiv_file = true;
        }
        else
        {
          printf( "%s|\n", optarg );
        }
        backup_file_name = new char[strlen(optarg) + 1];
        sprintf( backup_file_name, "%s", optarg );
        break;
    }
  }
  
  
  // =================================================================
  //                       Read the backup file
  // =================================================================
  
  ae_individual*  best_indiv;
  ae_environment* env;
  int32_t         num_gener;
  
  if ( backup_file_name == NULL )
  {
    printf("You must specify a backup file. Please use the option -f or --file.\n");
    exit(EXIT_FAILURE);
  }
  else
  {
    if ( use_single_indiv_file )
    {
      // TODO : best* backups don't look right...
      printf( "Reading single individual backup file <%s>... ", backup_file_name );
      gzFile* backup_file = (gzFile*) gzopen( backup_file_name, "r" );
      ae_common::read_from_backup( backup_file, verbose );
      env = new ae_environment(); // Uses the ae_common data
      best_indiv = new ae_individual( backup_file );
      
      num_gener = -1; // TODO!!!
      printf("done\n");
    }
    else
    {
      printf( "Reading backup file <%s>... ", backup_file_name );
      fflush(stdout);

      // Load simulation from backup
      ae_common::sim  = new ae_simulation ( backup_file_name, false );
      best_indiv      = ae_common::sim->get_pop()->get_best();
      env             = ae_common::sim->get_env();
      num_gener       = ae_common::sim->get_num_gener();
      printf("done\n");
    }
  }
  
  
  // The constructor of the ae_simulation has read the genomes of the individuals
  // and located their promoters, but has not performed the translation nor the
  // phenotype computation. We must do it now.
  // However, as the individuals in the backups are sorted, we don't need to evaluate
  // all the individuals, only those we are interested in (here only the best one)
    
  
  best_indiv->evaluate( env );


    
  // =================================================================
  //                      Create the EPS files 
  // =================================================================

  char directory_name[30];
  sprintf( directory_name, "files-generation%06"PRId32, num_gener );

  // Check whether the directory already exists and is writable
  if ( access( directory_name, F_OK ) == 0 )
  {
//       struct stat status;
//       stat( directory_name, &status );
//       if ( status.st_mode & S_IFDIR ) cout << "The directory exists." << endl;
//       else cout << "This path is a file." << endl;

    if ( access( directory_name, X_OK | W_OK) != 0 )
    {
      fprintf(stderr, "Error: cannot enter or write in directory %s.\n", directory_name);
      exit( EXIT_FAILURE );
    }
  }
  else
  {
    // Create the directory with permissions : rwx r-x r-x
    if ( mkdir( directory_name, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH) != 0 )
    {
      fprintf(stderr, "Error: cannot create directory %s.\n", directory_name);
      exit( EXIT_FAILURE );
    }
  }

  
  
  // =================================================================
  //                  Write the data in the EPS files 
  // =================================================================

  ae_genetic_unit*  best_indiv_main_genome = best_indiv->get_genetic_unit( 0 );
    
  printf( "Creating the EPS file with the triangles of the best individual... " );
  fflush(stdout);
  draw_triangles( best_indiv, env, directory_name );
  printf( "OK\n" );

  printf( "Creating the EPS file with the positive and negatives profiles of the best individual... " );
  fflush(stdout);
  draw_pos_neg_profiles( best_indiv, env, directory_name );
  printf( "OK\n" );

    
  printf( "Creating the EPS file with the phenotype of the best individual... " );
  fflush(stdout);  
  draw_phenotype( best_indiv, env, directory_name );
  printf( "OK\n" );

  printf( "Creating the EPS file with the CDS of the best individual... " );
  fflush(stdout);  
  draw_genetic_unit_with_CDS( best_indiv_main_genome, directory_name );
  printf( "OK\n" );

  printf( "Creating the EPS file with the mRNAs of the best individual... " );
  fflush(stdout);  
  draw_genetic_unit_with_mRNAs( best_indiv_main_genome, directory_name );
  printf( "OK\n" );


  delete [] backup_file_name;
  
  if ( use_single_indiv_file )
  {
    delete best_indiv;
    delete env;
  }
  else
  {
    delete ae_common::sim;
  }

  return EXIT_SUCCESS;
}




void print_help( void ) 
{
}




void draw_triangles( ae_individual* indiv, ae_environment* env, char * directoryName )
{
  const uint8_t bbsize = 200;  // a4 paper: 595*842 
  double margin = 0.1;
  double scale = 0.8*(1 - 2*margin);
  
  char filename[50];
  strncpy( filename, directoryName, 29 );
  strcat(  filename, "/best_triangles.eps" );
  FILE * drawingfile = fopen( filename, "w" );


  fprintf(drawingfile, "%%!PS-Adobe-3.0 EPSF-3.0\n");
  fprintf(drawingfile, "%%%%BoundingBox: 0 0 %d %d\n", bbsize, bbsize);
  fprintf(drawingfile, "%d %d scale\n", bbsize, bbsize);
  fprintf(drawingfile, "%d %d 8 [100 0 0 -100 0 100]\n",bbsize, bbsize);
  fprintf(drawingfile, "{currentfile 3 100 mul string readhexstring pop} bind\n");


  // -----------------------------
  //  paint neutral zones in grey
  // -----------------------------
  if ( ae_common::env_axis_is_segmented )
  {
    int16_t nb_segments = env->get_nb_segments();
    ae_env_segment** segments = env->get_segments();
    
    for ( int16_t i = 0 ; i < nb_segments ; i++ )
    {
      if ( segments[i]->feature == NEUTRAL )
      {    
        fprintf( drawingfile, "%lf 0 moveto\n", margin + scale * segments[i]->start );
        fprintf( drawingfile, "%lf 1 lineto\n", margin + scale * segments[i]->start );
        fprintf( drawingfile, "%lf 1 lineto\n", margin + scale * segments[i]->stop );
        fprintf( drawingfile, "%lf 0 lineto\n", margin + scale * segments[i]->stop );
        fprintf( drawingfile, "closepath\n" );
        fprintf( drawingfile, "0.8 setgray\n" );
        fprintf( drawingfile, "fill\n" );
        fprintf( drawingfile, "0 setgray\n" );
      }
    }
  }
 
  // -----------
  //  draw axis
  // -----------

  double arrowsize = 0.03;
  double arrowangle = 3.14/6;

  fprintf(drawingfile, "0.001 setlinewidth\n");
  fprintf(drawingfile, "0 0 0 setrgbcolor\n");

  // axis X + arrow
  fprintf(drawingfile, "%lf %lf moveto\n", margin/2, 0.5);
  fprintf(drawingfile, "%lf %lf lineto\n", 1-margin, 0.5);
  fprintf(drawingfile, "stroke\n");
  fprintf(drawingfile, "%lf %lf moveto\n", 1-margin, 0.5);
  fprintf(drawingfile, "%lf %lf lineto\n", 1-margin-arrowsize*cos(arrowangle), 0.5 + arrowsize*sin(arrowangle) );
  fprintf(drawingfile, "stroke\n");
  fprintf(drawingfile, "%lf %lf moveto\n", 1-margin, 0.5);
  fprintf(drawingfile, "%lf %lf lineto\n", 1-margin-arrowsize*cos(arrowangle), 0.5 - arrowsize*sin(arrowangle) );
  fprintf(drawingfile, "stroke\n");

  // axis Y + arrow
  fprintf( drawingfile, "%lf %lf moveto\n", margin, margin/2);
  fprintf( drawingfile, "%lf %lf lineto\n", margin, 1-margin);
  fprintf( drawingfile, "%lf %lf moveto\n", margin, 1-margin);
  fprintf( drawingfile, "%lf %lf lineto\n", margin-arrowsize*sin(arrowangle), 1 - margin - arrowsize*cos(arrowangle) );
  fprintf( drawingfile, "stroke\n");
  fprintf( drawingfile, "%lf %lf moveto\n", margin, 1-margin);
  fprintf( drawingfile, "%lf %lf lineto\n", margin+arrowsize*sin(arrowangle), 1 - margin - arrowsize*cos(arrowangle) );
  fprintf( drawingfile, "stroke\n");



  // ----------------
  //  draw triangles
  // ----------------

  fprintf(drawingfile,"[ ] 0 setdash\n");

  double h;
  ae_list_node*     gen_unit_node = indiv->get_genetic_unit_list()->get_first();
  ae_genetic_unit*  gen_unit      = NULL;
  
  while ( gen_unit_node != NULL )
  {
    gen_unit = (ae_genetic_unit*) gen_unit_node->get_obj();
    
    ae_list_node* prot_node  = NULL;
    ae_protein*   prot       = NULL;
    
    prot_node = (gen_unit->get_protein_list())[LEADING]->get_first();
    while ( prot_node != NULL )
    {
      prot = (ae_protein*) prot_node->get_obj();
      h = prot->get_height() * prot->get_concentration();
      fprintf(drawingfile, "%lf %lf moveto\n", margin, 0.5);
      fprintf(drawingfile, "%lf %lf lineto\n", margin + scale*(prot->get_mean() - prot->get_width()), 0.5);
      fprintf(drawingfile, "%lf %lf lineto\n", margin + scale*(prot->get_mean()), 0.5 + scale*(h));
      fprintf(drawingfile, "%lf %lf lineto\n", margin + scale*(prot->get_mean() + prot->get_width()), 0.5);
      fprintf(drawingfile, "%lf %lf moveto\n", margin + scale*(1), 0.5);
      fprintf(drawingfile, "stroke\n");
      prot_node = prot_node->get_next();
    }


    prot_node = (gen_unit->get_protein_list())[LAGGING]->get_first();
    while ( prot_node != NULL )
    {
      prot = (ae_protein*) prot_node->get_obj();
      h = prot->get_height() * prot->get_concentration();
      fprintf(drawingfile, "%lf %lf moveto\n", margin, 0.5);
      fprintf(drawingfile, "%lf %lf lineto\n", margin + scale*(prot->get_mean() - prot->get_width()), 0.5);
      fprintf(drawingfile, "%lf %lf lineto\n", margin + scale*(prot->get_mean()), 0.5 + scale*(h));
      fprintf(drawingfile, "%lf %lf lineto\n", margin + scale*(prot->get_mean() + prot->get_width()), 0.5);
      fprintf(drawingfile, "%lf %lf moveto\n", margin + scale*(1), 0.5);
      fprintf(drawingfile, "stroke\n");
      
      prot_node = prot_node->get_next();
    }
    
    gen_unit_node = gen_unit_node->get_next();
  }



  fprintf(drawingfile,"%%%%EOF\n");
  fclose(drawingfile);

}



void draw_pos_neg_profiles( ae_individual * indiv, ae_environment* env, char * directoryName )
{
  const uint8_t bbsize = 200;  // a4 paper: 595*842 
  double margin = 0.1;
  double scale = 0.8*(1 - 2*margin);
  
  char filename[50];
  strncpy( filename, directoryName, 29 );
  strcat(  filename, "/best_pos_neg_profiles.eps" );
  FILE * drawingfile = fopen( filename, "w" );


  fprintf(drawingfile, "%%!PS-Adobe-3.0 EPSF-3.0\n");
  fprintf(drawingfile, "%%%%BoundingBox: 0 0 %d %d\n", bbsize, bbsize);
  fprintf(drawingfile, "%d %d scale\n", bbsize, bbsize);
  fprintf(drawingfile, "%d %d 8 [100 0 0 -100 0 100]\n",bbsize, bbsize);
  fprintf(drawingfile, "{currentfile 3 100 mul string readhexstring pop} bind\n");


  // -----------------------------
  //  paint neutral zones in grey
  // -----------------------------
  if ( ae_common::env_axis_is_segmented )
  {
    int16_t nb_segments = env->get_nb_segments();
    ae_env_segment** segments = env->get_segments();
    
    for ( int16_t i = 0 ; i < nb_segments ; i++ )
    {
      if ( segments[i]->feature == NEUTRAL )
      {    
        fprintf( drawingfile, "%lf 0 moveto\n", margin + scale * segments[i]->start );
        fprintf( drawingfile, "%lf 1 lineto\n", margin + scale * segments[i]->start );
        fprintf( drawingfile, "%lf 1 lineto\n", margin + scale * segments[i]->stop );
        fprintf( drawingfile, "%lf 0 lineto\n", margin + scale * segments[i]->stop );
        fprintf( drawingfile, "closepath\n" );
        fprintf( drawingfile, "0.8 setgray\n" );
        fprintf( drawingfile, "fill\n" );
        fprintf( drawingfile, "0 setgray\n" );
      }
    }
  }
  
 
  // -----------
  //  draw axis
  // -----------

  double arrowsize = 0.03;
  double arrowangle = 3.14/6;

  fprintf(drawingfile, "0.001 setlinewidth\n");
  fprintf(drawingfile, "0 0 0 setrgbcolor\n");

  // axis X + arrow
  fprintf(drawingfile, "%lf %lf moveto\n", margin/2, 0.5);
  fprintf(drawingfile, "%lf %lf lineto\n", 1-margin, 0.5);
  fprintf(drawingfile, "stroke\n");
  fprintf(drawingfile, "%lf %lf moveto\n", 1-margin, 0.5);
  fprintf(drawingfile, "%lf %lf lineto\n", 1-margin-arrowsize*cos(arrowangle), 0.5 + arrowsize*sin(arrowangle) );
  fprintf(drawingfile, "stroke\n");
  fprintf(drawingfile, "%lf %lf moveto\n", 1-margin, 0.5);
  fprintf(drawingfile, "%lf %lf lineto\n", 1-margin-arrowsize*cos(arrowangle), 0.5 - arrowsize*sin(arrowangle) );
  fprintf(drawingfile, "stroke\n");

  // axis Y + arrow
  fprintf(drawingfile, "%lf %lf moveto\n", margin, margin/2);
  fprintf(drawingfile, "%lf %lf lineto\n", margin, 1-margin);
  fprintf(drawingfile, "%lf %lf moveto\n", margin, 1-margin);
  fprintf(drawingfile, "%lf %lf lineto\n", margin-arrowsize*sin(arrowangle), 1 - margin - arrowsize*cos(arrowangle) );
  fprintf(drawingfile, "stroke\n");
  fprintf(drawingfile, "%lf %lf moveto\n", margin, 1-margin);
  fprintf(drawingfile, "%lf %lf lineto\n", margin+arrowsize*sin(arrowangle), 1 - margin - arrowsize*cos(arrowangle) );
  fprintf(drawingfile, "stroke\n");


  // -----------------------
  //  draw positive profile
  // -----------------------

  fprintf(drawingfile,"[ ] 0 setdash\n");
  fprintf(drawingfile, "0.002 setlinewidth\n");
  fprintf(drawingfile, "%lf %lf moveto\n", margin, 0.5);

  ae_list_node * node = ((indiv->get_phenotype_activ())->get_points())->get_first();
  ae_point_2d * pt = NULL;
  while (node != NULL)
  {
    pt = (ae_point_2d *) node->get_obj();
    fprintf(drawingfile, "%lf %lf lineto\n", margin + scale*pt->x, 0.5 + scale*pt->y);
    node = node->get_next();
  }
  fprintf(drawingfile, "stroke\n");
    

  // -----------------------
  //  draw negative profile
  // -----------------------

  fprintf(drawingfile,"[ ] 0 setdash\n");
  fprintf(drawingfile, "0.002 setlinewidth\n");
  fprintf(drawingfile, "%lf %lf moveto\n", margin, 0.5);
    
  node = ((indiv->get_phenotype_inhib())->get_points())->get_first();
  pt = NULL;
  while (node != NULL)
  {
    pt = (ae_point_2d *) node->get_obj();
    fprintf(drawingfile, "%lf %lf lineto\n", margin + scale*pt->x, 0.5 + scale*pt->y);
    node = node->get_next();
  }
  fprintf(drawingfile, "stroke\n");

  fprintf(drawingfile,"%%%%EOF\n");
  fclose(drawingfile);

}




void draw_phenotype( ae_individual * indiv, ae_environment * envir, char * directoryName )
{
  const uint8_t bbsize = 200;  // a4 paper: 595*842 
  double margin = 0.1;
  double scale = 0.8*(1 - 2*margin);
  

  char filename[50];
  strncpy( filename, directoryName, 29 );
  strcat(  filename, "/best_phenotype.eps" );
  FILE * drawingfile = fopen( filename, "w" );
 
  if (drawingfile == NULL)
    {
      fprintf(stderr, "Error: could not create the file %s\n", filename);
      return;
    }
  

  fprintf(drawingfile, "%%!PS-Adobe-3.0 EPSF-3.0\n");
  fprintf(drawingfile, "%%%%BoundingBox: 0 0 %d %d\n", bbsize, bbsize);
  fprintf(drawingfile, "%d %d scale\n", bbsize, bbsize);
  fprintf(drawingfile, "%d %d 8 [100 0 0 -100 0 100]\n",bbsize, bbsize);
  fprintf(drawingfile, "{currentfile 3 100 mul string readhexstring pop} bind\n");


  // -----------------------------
  //  paint neutral zones in grey
  // -----------------------------
  if ( ae_common::env_axis_is_segmented )
  {
    int16_t nb_segments = envir->get_nb_segments();
    ae_env_segment** segments = envir->get_segments();
    
    for ( int16_t i = 0 ; i < nb_segments ; i++ )
    {
      if ( segments[i]->feature == NEUTRAL )
      {    
        fprintf( drawingfile, "%lf 0 moveto\n", margin + scale * segments[i]->start );
        fprintf( drawingfile, "%lf 1 lineto\n", margin + scale * segments[i]->start );
        fprintf( drawingfile, "%lf 1 lineto\n", margin + scale * segments[i]->stop );
        fprintf( drawingfile, "%lf 0 lineto\n", margin + scale * segments[i]->stop );
        fprintf( drawingfile, "closepath\n" );
        fprintf( drawingfile, "0.8 setgray\n" );
        fprintf( drawingfile, "fill\n" );
        fprintf( drawingfile, "0 setgray\n" );
      }
    }
  }
 
    
  // -----------
  //  draw axis
  // -----------

  double arrowsize = 0.03;
  double arrowangle = 3.14/6;

  fprintf(drawingfile, "0.001 setlinewidth\n");
  fprintf(drawingfile, "0 0 0 setrgbcolor\n");

  // axis X + arrow
  fprintf(drawingfile, "%lf %lf moveto\n", margin/2, margin);
  fprintf(drawingfile, "%lf %lf lineto\n", 1-margin, margin);
  fprintf(drawingfile, "stroke\n");
  fprintf(drawingfile, "%lf %lf moveto\n", 1-margin, margin);
  fprintf(drawingfile, "%lf %lf lineto\n", 1-margin-arrowsize*cos(arrowangle), margin + arrowsize*sin(arrowangle) );
  fprintf(drawingfile, "stroke\n");
  fprintf(drawingfile, "%lf %lf moveto\n", 1-margin, margin);
  fprintf(drawingfile, "%lf %lf lineto\n", 1-margin-arrowsize*cos(arrowangle), margin - arrowsize*sin(arrowangle) );
  fprintf(drawingfile, "stroke\n");

  // axis Y + arrow
  fprintf(drawingfile, "%lf %lf moveto\n", margin, margin/2);
  fprintf(drawingfile, "%lf %lf lineto\n", margin, 1-margin);
  fprintf(drawingfile, "%lf %lf moveto\n", margin, 1-margin);
  fprintf(drawingfile, "%lf %lf lineto\n", margin-arrowsize*sin(arrowangle), 1 - margin - arrowsize*cos(arrowangle) );
  fprintf(drawingfile, "stroke\n");
  fprintf(drawingfile, "%lf %lf moveto\n", margin, 1-margin);
  fprintf(drawingfile, "%lf %lf lineto\n", margin+arrowsize*sin(arrowangle), 1 - margin - arrowsize*cos(arrowangle) );
  fprintf(drawingfile, "stroke\n");

  // max degree = 1
  fprintf(drawingfile, "[0.02 0.02] 0 setdash\n");
  fprintf(drawingfile, "%lf %lf moveto\n", margin, margin + 1*scale);
  fprintf(drawingfile, "%lf %lf lineto\n", 1-margin, margin + 1*scale);
  fprintf(drawingfile, "stroke\n");


  // ----------------
  //  draw phenotype
  // ----------------
  fprintf(drawingfile,"[ ] 0 setdash\n");
  fprintf(drawingfile, "0.002 setlinewidth\n");
  fprintf(drawingfile, "%lf %lf moveto\n", margin, margin);
  ae_list_node * node = ((indiv->get_phenotype())->get_points())->get_first();
  ae_point_2d * pt = NULL;
  while (node != NULL)
    {
      pt = (ae_point_2d *) node->get_obj();
      fprintf(drawingfile, "%lf %lf lineto\n", margin + scale*pt->x, margin + scale*pt->y);
      node = node->get_next();
    }
  fprintf(drawingfile, "stroke\n");


  // ------------------
  //  draw environment
  // ------------------
  fprintf(drawingfile,"[ ] 0 setdash\n");
  fprintf(drawingfile, "0.001 setlinewidth\n");
  fprintf(drawingfile, "%lf %lf moveto\n", margin, margin);
  node = (envir->get_points())->get_first();
  pt = NULL;
  while (node != NULL)
    {
      pt = (ae_point_2d *) node->get_obj();
      fprintf(drawingfile, "%lf %lf lineto\n", margin + scale*pt->x, margin + scale*pt->y);
      node = node->get_next();
    }
  fprintf(drawingfile, "stroke\n");



  fprintf(drawingfile,"%%%%EOF\n");
  fclose(drawingfile);
 
}







void draw_genetic_unit_with_CDS( ae_genetic_unit* gen_unit, char * directoryName )
{
  const uint8_t bbsize = 200;  // a4 paper: 595*842 
  int32_t gen_length = (gen_unit->get_dna())->get_length();
  double r = 0.35;
  double scale = 2*M_PI*r/gen_length;

  char filename[50];
  strncpy( filename, directoryName, 29 );
  strcat(  filename, "/best_genome_with_CDS.eps" );
  FILE * drawingfile = fopen( filename, "w" );

  fprintf(drawingfile, "%%!PS-Adobe-3.0 EPSF-3.0\n");
  fprintf(drawingfile, "%%%%BoundingBox: 0 0 %d %d\n", bbsize, bbsize);
  fprintf(drawingfile, "%d %d scale\n", bbsize, bbsize);
  fprintf(drawingfile, "%d %d 8 [100 0 0 -100 0 100]\n",bbsize, bbsize);
  fprintf(drawingfile, "{currentfile 3 100 mul string readhexstring pop} bind\n");
 
  // -----------
  //  chromosome
  // -----------
  
  fprintf(drawingfile, "0.001 setlinewidth\n");
  fprintf(drawingfile, "0 0 0 setrgbcolor\n");
  fprintf(drawingfile, "%lf %lf %lf 0 360 arc\n", 0.5, 0.5, r); // arcn = clockwise arc
  fprintf(drawingfile, "stroke\n");

  // -----------
  //  scale
  // -----------

  double scalesize = 0.15;
  fprintf(drawingfile, "%lf %lf moveto\n", 0.5-scalesize/2, 0.5);
  fprintf(drawingfile, "%lf %lf lineto\n", 0.5+scalesize/2, 0.5);
  fprintf(drawingfile, "stroke\n");
  fprintf(drawingfile, "/Helvetica findfont\n");   
  fprintf(drawingfile, "0.035 scalefont\n");           
  fprintf(drawingfile, "setfont\n"); 
  fprintf(drawingfile, "newpath\n"); 
  fprintf(drawingfile, "%lf %lf moveto\n", 0.5-scalesize/3, 0.52);
  fprintf(drawingfile, "(scale : %.0lf bp) show\n", scalesize/scale); 

  // -----------
  //  genes
  // -----------

  ae_list_node* node  = NULL;
  ae_protein*   prot  = NULL;

  int32_t first;
  int32_t last;
  int8_t  layer = 0;
  int8_t  outmost_layer = 1;
  int16_t nb_sect;
  int16_t rho;
  int16_t angle;
  bool    sectors_free;
  int16_t alpha_first;
  int16_t alpha_last;
  int16_t theta_first;
  int16_t theta_last;
  

  // To handle gene overlaps, we remember where we have aldready drawn
  // something, with a precision of 1 degree. We handle up to 100 layers:
  //  - 50 layers inside the circle (lagging strand), 
  //  - 50 layers outside the cricle (leading strand).
  // At first, only one layer is created, we create new layers when we 
  // need them. 
  bool* occupied_sectors[2][50];
  occupied_sectors[LEADING][0] = new bool[360];
  occupied_sectors[LAGGING][0] = new bool[360];
  for ( int16_t angle = 0 ; angle < 360 ; angle++ )
  {
    occupied_sectors[LEADING][0][angle] = false;
    occupied_sectors[LAGGING][0][angle] = false;
  }


  printf("LEADING\n");
  node = (gen_unit->get_protein_list())[LEADING]->get_first();
  while (node != NULL)
  {
    prot = (ae_protein *) node->get_obj();
    first = prot->get_first_translated_pos();
    last = prot->get_last_translated_pos();
    // h = prot->get_height() * prot->get_concentration();

    alpha_first   = (int16_t) round(  (double)(360 * first) / (double)gen_length );  //  == sect1 == alphaB
    alpha_last    = (int16_t) round(  (double)(360 * last)  / (double)gen_length );  //  == sect2 == alphaA
    theta_first   = utils::mod( 90 - alpha_first, 360 );  //  == tetaB
    theta_last    = utils::mod( 90 - alpha_last, 360 );  //   == tetaA
    if ( theta_first == theta_last ) theta_first = utils::mod( theta_first + 1, 360 ); 

    nb_sect = utils::mod( theta_first - theta_last + 1, 360 );


    // Outside the circle, look for the inmost layer that has all the sectors between
    // theta_first and theta_last free
    layer = 0;
    sectors_free = false;
    while ( ! sectors_free )
    {
      sectors_free = true;
      for ( rho = 0 ; rho < nb_sect ; rho++ )
      {
        if ( occupied_sectors[LEADING][layer][utils::mod(theta_first - rho, 360)] )
        {
          sectors_free = false;
          break;
        }
      }
        
      if ( sectors_free )
      {
        break; // All the needed sectors are free on the current layer
      }
      else
      {
        layer++;

        if ( (layer >= outmost_layer) && (layer < 49) )
        {
          // Add a new layer (actually, 2 layers, to maintain the symmetry)
          occupied_sectors[LEADING][outmost_layer] = new bool[360];
          occupied_sectors[LAGGING][outmost_layer] = new bool[360];
          for ( int16_t angle = 0 ; angle < 360 ; angle++ )
          {
            occupied_sectors[LEADING][outmost_layer][angle] = false;
            occupied_sectors[LAGGING][outmost_layer][angle] = false;
          }
                
          outmost_layer++;
          break; // A new layer is necessarily free, no need to look further
        }
        if ( layer == 49 ) 
        {
          // We shall not create a 51th layer, the CDS will be drawn on the
          // layer, probably over another CDS
          break;
        }
      }
    }
    
    printf("f %d, l %d, af %d, al %d, tf %d, tl %d, nbsect %d, layer %d\n", first, last, alpha_first, alpha_last, theta_first, theta_last, nb_sect, layer);

    // Mark sectors to be drawn as occupied
    for ( rho = 0 ; rho < nb_sect ; rho++ )
    {
      occupied_sectors[LEADING][layer][utils::mod(theta_first - rho, 360)] = true;
    }
    // Mark flanking sectors as occupied
    occupied_sectors[LEADING][layer][utils::mod(theta_first + 1, 360)] = true;
    occupied_sectors[LEADING][layer][utils::mod(theta_first - nb_sect, 360)] = true;
    
    
    // draw !
    fprintf(drawingfile, "0.018 setlinewidth\n");
    // fprintf(drawingfile, "%lf %lf %lf setrgbcolor\n",  1-(0.8*h/max_height + 0.2), 1-(0.8*h/max_height + 0.2),1-(0.8*h/max_height + 0.2));
    layer++; // index starting at 0 but needed to start at 1
    
    if (theta_last > theta_first) 
    {
      fprintf(drawingfile, "newpath\n");
      fprintf(drawingfile, "%lf %lf %lf %d %d arc\n", 0.5, 0.5, r + layer*0.02, theta_last, 360);
      fprintf(drawingfile, "stroke\n");
      fprintf(drawingfile, "%lf %lf %lf %d %d arc\n", 0.5, 0.5, r + layer*0.02, 0, theta_first);
      fprintf(drawingfile, "stroke\n");
    }
    else
    {
      fprintf(drawingfile, "newpath\n");
      fprintf(drawingfile, "%lf %lf %lf %d %d arc\n", 0.5, 0.5, r + layer*0.02, theta_last, theta_first);
      fprintf(drawingfile, "stroke\n");
    }
    
    node = node->get_next();
  }
  

  printf("LAGGING\n");
  node = (gen_unit->get_protein_list())[LAGGING]->get_first();
  while (node != NULL)
  {
    prot = (ae_protein *) node->get_obj();
    first = prot->get_first_translated_pos();
    last = prot->get_last_translated_pos();
    // h = prot->get_height() * prot->get_concentration();

    alpha_first   = (int16_t) round(  (double)(360 * first) / (double)gen_length );  
    alpha_last    = (int16_t) round(  (double)(360 * last)  / (double)gen_length );
    theta_first   = utils::mod( 90 - alpha_first, 360 );  
    theta_last    = utils::mod( 90 - alpha_last, 360 );     
    if ( theta_first == theta_last ) theta_last = utils::mod( theta_last + 1, 360 ); 

    nb_sect = utils::mod( theta_last - theta_first + 1, 360 );


    // Inside the circle, look for the inmost layer that has all the sectors between
    // theta_first and theta_last free
    layer = 0;
    sectors_free = false;
    while ( ! sectors_free )
    {
      sectors_free = true;
      for ( rho = 0 ; rho < nb_sect ; rho++ )
      {
        if ( occupied_sectors[LAGGING][layer][utils::mod(theta_first + rho, 360)] )
        {
          sectors_free = false;
          break;
        }
      }
      
      if ( sectors_free )
      {
        break; // All the needed sectors are free on the current layer
      }
      else
      {
        layer++;

        if ( (layer >= outmost_layer) && (layer < 49) )
        {
          // Add a new layer (actually, 2 layers, to maintain the symmetry)
          occupied_sectors[LEADING][outmost_layer] = new bool[360];
          occupied_sectors[LAGGING][outmost_layer] = new bool[360];
          for ( angle = 0 ; angle < 360 ; angle++ )
          {
            occupied_sectors[LEADING][outmost_layer][angle] = false;
            occupied_sectors[LAGGING][outmost_layer][angle] = false;
          }
          
          outmost_layer++;
          break; // A new layer is necessarily free, no need to look further
        }
        if ( layer == 49 ) 
        {
          // We shall not create a 51th layer, the CDS will be drawn on the
          // layer, probably over another CDS
          break;
        }
      }
    }

    printf("f %d, l %d, af %d, al %d, tf %d, tl %d, nbsect %d, layer %d\n", first, last, alpha_first, alpha_last, theta_first, theta_last, nb_sect, layer);

    // Mark sectors to be drawn as occupied
    for ( rho = 0 ; rho < nb_sect ; rho++ )
    {
      occupied_sectors[LAGGING][layer][utils::mod(theta_first + rho, 360)] = true;
    }
    // Mark flanking sectors as occupied
    occupied_sectors[LAGGING][layer][utils::mod(theta_first - 1, 360)] = true;
    occupied_sectors[LAGGING][layer][utils::mod(theta_first + nb_sect, 360)] = true;

  
    // draw !
    fprintf(drawingfile, "0.018 setlinewidth\n");
    // fprintf(drawingfile, "%lf %lf %lf setrgbcolor\n",  1-(0.8*h/max_height + 0.2), 1-(0.8*h/max_height + 0.2),1-(0.8*h/max_height + 0.2));
    layer++; // index starting at 0 but needed to start at 1

    if (theta_first > theta_last) 
    {
      fprintf(drawingfile, "newpath\n");
      fprintf(drawingfile, "%lf %lf %lf %d %d arc\n", 0.5, 0.5, r - layer*0.02, theta_first, 360);
      fprintf(drawingfile, "stroke\n");
      fprintf(drawingfile, "%lf %lf %lf %d %d arc\n", 0.5, 0.5, r - layer*0.02, 0, theta_last);
      fprintf(drawingfile, "stroke\n");
    }
    else
    {
      fprintf(drawingfile, "newpath\n");
      fprintf(drawingfile, "%lf %lf %lf %d %d arc\n", 0.5, 0.5, r - layer*0.02, theta_first, theta_last);
      fprintf(drawingfile, "stroke\n");
    }
      
    node = node->get_next();
  }


  fprintf(drawingfile,"showpage\n");
  fprintf(drawingfile,"%%%%EOF\n");
  fclose(drawingfile);

  for ( layer = 0 ; layer < outmost_layer ; layer++ )
  {
    delete occupied_sectors[LEADING][layer];
    delete occupied_sectors[LAGGING][layer];
  }

}




void draw_genetic_unit_with_mRNAs( ae_genetic_unit* gen_unit, char * directoryName )
{
  const uint8_t bbsize = 200;  // a4 paper: 595*842 
  int32_t gen_length = (gen_unit->get_dna())->get_length();
  double r = 0.35;
  double scale = 2*M_PI*r/gen_length;

  char filename[50];
  strncpy( filename, directoryName, 29 );
  strcat(  filename, "/best_genome_with_mRNAs.eps" );
  FILE * drawingfile = fopen( filename, "w" );

  fprintf(drawingfile, "%%!PS-Adobe-3.0 EPSF-3.0\n");
  fprintf(drawingfile, "%%%%BoundingBox: 0 0 %d %d\n", bbsize, bbsize);
  fprintf(drawingfile, "%d %d scale\n", bbsize, bbsize);
  fprintf(drawingfile, "%d %d 8 [100 0 0 -100 0 100]\n",bbsize, bbsize);
  fprintf(drawingfile, "{currentfile 3 100 mul string readhexstring pop} bind\n");
 
  // -----------
  //  chromosome
  // -----------
  
  fprintf(drawingfile, "0.001 setlinewidth\n");
  fprintf(drawingfile, "0 0 0 setrgbcolor\n");
  fprintf(drawingfile, "%lf %lf %lf 0 360 arc\n", 0.5, 0.5, r); // arcn = clockwise arc
  fprintf(drawingfile, "stroke\n");

  // -----------
  //  scale
  // -----------

  double scalesize = 0.15;
  fprintf(drawingfile, "%lf %lf moveto\n", 0.5-scalesize/2, 0.5);
  fprintf(drawingfile, "%lf %lf lineto\n", 0.5+scalesize/2, 0.5);
  fprintf(drawingfile, "stroke\n");
  fprintf(drawingfile, "/Helvetica findfont\n");   
  fprintf(drawingfile, "0.035 scalefont\n");           
  fprintf(drawingfile, "setfont\n"); 
  fprintf(drawingfile, "newpath\n"); 
  fprintf(drawingfile, "%lf %lf moveto\n", 0.5-scalesize/3, 0.52);
  fprintf(drawingfile, "(scale : %.0lf bp) show\n", scalesize/scale); 

  // -----------
  //  mRNAs
  // -----------

  ae_list_node * node  = NULL;
  ae_rna *  rna  = NULL;

  int32_t first;
  int32_t last;
  int8_t  layer = 0;
  int8_t  outmost_layer = 1;
  int16_t nb_sect;
  int16_t rho;
  int16_t angle;
  bool    sectors_free;
  int16_t alpha_first;
  int16_t alpha_last;
  int16_t theta_first;
  int16_t theta_last;
  

  // To handle gene overlaps, we remember where we have aldready drawn
  // something, with a precision of 1 degree. We handle up to 100 layers:
  //  - 50 layers inside the circle (lagging strand), 
  //  - 50 layers outside the cricle (leading strand).
  // At first, only one layer is created, we create new layers when we 
  // need them. 
  bool* occupied_sectors[2][50];
  occupied_sectors[LEADING][0] = new bool[360];
  occupied_sectors[LAGGING][0] = new bool[360];
  for ( int16_t angle = 0 ; angle < 360 ; angle++ )
  {
    occupied_sectors[LEADING][0][angle] = false;
    occupied_sectors[LAGGING][0][angle] = false;
  }



  node = (gen_unit->get_rna_list())[LEADING]->get_first();
  while (node != NULL)
  {
    rna = (ae_rna *) node->get_obj();
    first = rna->get_first_transcribed_pos();
    last = rna->get_last_transcribed_pos();
    
    
    alpha_first   = (int16_t) round(  (double)(360 * first) / (double)gen_length );  //  == sect1 == alphaB
    alpha_last    = (int16_t) round(  (double)(360 * last)  / (double)gen_length );  //  == sect2 == alphaA
    theta_first   = utils::mod( 90 - alpha_first, 360 );  //  == tetaB
    theta_last    = utils::mod( 90 - alpha_last, 360 );  //   == tetaA
    if ( theta_first == theta_last ) theta_first = utils::mod( theta_first + 1, 360 ); 

    nb_sect = utils::mod( theta_first - theta_last + 1, 360 );


    // Outside the circle, look for the inmost layer that has all the sectors between
    // theta_first and theta_last free
    layer = 0;
    sectors_free = false;
    while ( ! sectors_free )
    {
      sectors_free = true;
      for ( rho = 0 ; rho < nb_sect ; rho++ )
      {
        if ( occupied_sectors[LEADING][layer][utils::mod(theta_first - rho, 360)] )
        {
          sectors_free = false;
          break;
        }
      }
          
      if ( sectors_free )
      {
        break; // All the needed sectors are free on the current layer
      }
      else
      {
        layer++;

        if ( (layer >= outmost_layer) && (layer < 49) )
        {
          // Add a new layer (actually, 2 layers, to maintain the symmetry)
          occupied_sectors[LEADING][outmost_layer] = new bool[360];
          occupied_sectors[LAGGING][outmost_layer] = new bool[360];
          for ( int16_t angle = 0 ; angle < 360 ; angle++ )
          {
            occupied_sectors[LEADING][outmost_layer][angle] = false;
            occupied_sectors[LAGGING][outmost_layer][angle] = false;
          }
                  
          outmost_layer++;
          break; // A new layer is necessarily free, no need to look further
        }
        if ( layer == 49 ) 
        {
          // We shall not create a 51th layer, the CDS will be drawn on the
          // layer, probably over another CDS
          break;
        }
      }
    }

    // Mark sectors to be drawn as occupied
    for ( rho = 0 ; rho < nb_sect ; rho++ )
    {
      occupied_sectors[LEADING][layer][utils::mod(theta_first - rho, 360)] = true;
    }

    // Mark flanking sectors as occupied
    occupied_sectors[LEADING][layer][utils::mod(theta_first + 1, 360)] = true;
    occupied_sectors[LEADING][layer][utils::mod(theta_first - nb_sect, 360)] = true;

    
    // draw !
    fprintf(drawingfile, "0.018 setlinewidth\n");
    if ( rna->is_coding() ) fprintf(drawingfile, "0 0 0 setrgbcolor\n");
    else fprintf(drawingfile, "0.7 0.7 0.7 setrgbcolor\n");
    layer++; // index starting at 0 but needed to start at 1

    if (theta_last > theta_first) 
    {
      fprintf(drawingfile, "newpath\n");
      fprintf(drawingfile, "%lf %lf %lf %d %d arc\n", 0.5, 0.5, r + layer*0.02, theta_last, 360);
      fprintf(drawingfile, "stroke\n");
      fprintf(drawingfile, "%lf %lf %lf %d %d arc\n", 0.5, 0.5, r + layer*0.02, 0, theta_first);
      fprintf(drawingfile, "stroke\n");
    }
    else
    {
      fprintf(drawingfile, "newpath\n");
      fprintf(drawingfile, "%lf %lf %lf %d %d arc\n", 0.5, 0.5, r + layer*0.02, theta_last, theta_first);
      fprintf(drawingfile, "stroke\n");
    }
        
    node = node->get_next();
  }
  


  node = (gen_unit->get_rna_list())[LAGGING]->get_first();
  while (node != NULL)
  {
    rna = (ae_rna *) node->get_obj();
    first = rna->get_first_transcribed_pos();
    last = rna->get_last_transcribed_pos();


    alpha_first   = (int16_t) round(  (double)(360 * first) / (double)gen_length );  
    alpha_last    = (int16_t) round(  (double)(360 * last)  / (double)gen_length );
    theta_first   = utils::mod( 90 - alpha_first, 360 );  
    theta_last    = utils::mod( 90 - alpha_last, 360 );     
    nb_sect = utils::mod( alpha_first - alpha_last + 1,  360 );


    // Inside the circle, look for the inmost layer that has all the sectors between
    // theta_first and theta_last free
    layer = 0;
    sectors_free = false;
    while ( ! sectors_free )
    {
      sectors_free = true;
      for ( rho = 0 ; rho < nb_sect ; rho++ )
      {
        if ( occupied_sectors[LAGGING][layer][utils::mod(theta_first + rho, 360)] )
        {
          sectors_free = false;
          break;
        }
      }
          
      if ( sectors_free )
      {
        break; // All the needed sectors are free on the current layer
      }
      else
      {
        layer++;

        if ( (layer >= outmost_layer) && (layer < 49) )
        {
          // Add a new layer (actually, 2 layers, to maintain the symmetry)
          occupied_sectors[LEADING][outmost_layer] = new bool[360];
          occupied_sectors[LAGGING][outmost_layer] = new bool[360];
          for ( angle = 0 ; angle < 360 ; angle++ )
          {
            occupied_sectors[LEADING][outmost_layer][angle] = false;
            occupied_sectors[LAGGING][outmost_layer][angle] = false;
          }
                  
          outmost_layer++;
          break; // A new layer is necessarily free, no need to look further
        }
        if ( layer == 49 ) 
        {
          // We shall not create a 51th layer, the CDS will be drawn on the
          // layer, probably over another CDS
          break;
        }
      }
    }

    // Mark sectors to be drawn as occupied
    for ( rho = 0 ; rho < nb_sect ; rho++ )
    {
      occupied_sectors[LAGGING][layer][utils::mod(theta_first + rho, 360)] = true;
    }

    // Mark flanking sectors as occupied
    occupied_sectors[LAGGING][layer][utils::mod(theta_first - 1, 360)] = true;
    occupied_sectors[LAGGING][layer][utils::mod(theta_first + nb_sect, 360)] = true;

    
    // draw !
    fprintf(drawingfile, "0.018 setlinewidth\n");
    if ( rna->is_coding() ) fprintf(drawingfile, "0 0 0 setrgbcolor\n");
    else fprintf(drawingfile, "0.7 0.7 0.7 setrgbcolor\n");
    layer++; // index starting at 0 but needed to start at 1

    if (theta_first > theta_last) 
    {
      fprintf(drawingfile, "newpath\n");
      fprintf(drawingfile, "%lf %lf %lf %d %d arc\n", 0.5, 0.5, r - layer*0.02, theta_first, 360);
      fprintf(drawingfile, "stroke\n");
      fprintf(drawingfile, "%lf %lf %lf %d %d arc\n", 0.5, 0.5, r - layer*0.02, 0, theta_last);
      fprintf(drawingfile, "stroke\n");
    }
    else
    {
      fprintf(drawingfile, "newpath\n");
      fprintf(drawingfile, "%lf %lf %lf %d %d arc\n", 0.5, 0.5, r - layer*0.02, theta_first, theta_last);
      fprintf(drawingfile, "stroke\n");
    }
        
    node = node->get_next();
  }


  fprintf(drawingfile,"showpage\n");
  fprintf(drawingfile,"%%%%EOF\n");
  fclose(drawingfile);

  for ( layer = 0 ; layer < outmost_layer ; layer++ )
  {
    delete occupied_sectors[LEADING][layer];
    delete occupied_sectors[LAGGING][layer];
  }

}
  