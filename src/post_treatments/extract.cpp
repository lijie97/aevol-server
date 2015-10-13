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
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <signal.h>


// =================================================================
//                            Project Files
// =================================================================
#include "aevol.h"

using namespace aevol;

// =================================================================
//                         Function declarations
// =================================================================
void print_help(char* prog_path);

void analyse_indiv(Individual* indiv, FILE* triangles_file, FILE* sequence_file,
                   int16_t gu, const PhenotypicTarget& phenotypicTarget);
void analyse_gu(GeneticUnit* gen_unit, int32_t gen_unit_number,
                FILE* triangles_file, const PhenotypicTarget& phenotypicTarget);




int main(int argc, char* argv[])
{
  // Initialize command-line option variables with default values
  char* pop_file_name  = NULL;
  char* triangles_file_name  = NULL;
  char* sequence_file_name  = NULL;
  bool best_only = false;
  int16_t gu = -1;
  int32_t num_gener = -1;

  // Define allowed options
  const char * options_list = "hVr:t:s:bg:";
  static struct option long_options_list[] = {
      {"help",      no_argument,        NULL, 'h'},
      {"version",   no_argument,        NULL, 'V'},
      {"resume",    required_argument,  NULL, 'r'},
      {"triangles", required_argument,  NULL, 't'},
      {"sequence",  required_argument,  NULL, 's'},
      {"best",      no_argument,        NULL, 'b'},
      {"gu",        required_argument,  NULL, 'g'},
      {0, 0, 0, 0}
  };

  // Get actual values of the command-line options
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
      case 'r':
        num_gener = atol(optarg);
        break;
      case 't' :
        triangles_file_name = new char[strlen(optarg) + 1];
        sprintf(triangles_file_name, "%s", optarg);
        break;
      case 's' :
        sequence_file_name = new char[strlen(optarg) + 1];
        sprintf(sequence_file_name, "%s", optarg);
        break;
      case 'g' :
        gu = atoi(optarg);
        break;
      case 'b' :
        best_only = true;
        break;
    }
  }

  // If num_gener is not provided, assume last gener
  if (num_gener == -1) {
    num_gener = OutputManager::get_last_gener();
  }

  if (triangles_file_name == NULL && sequence_file_name == NULL) {
    Utils::ExitWithDevMsg("Use option -s or -t (-h for more info)",
                          __FILE__, __LINE__);
  }

  // Open the files
  FILE* triangles_file = NULL;
  FILE* sequence_file = NULL;

  if (triangles_file_name != NULL)
  {
    triangles_file = fopen(triangles_file_name,"w");

    // Write file headers
    int key = 1;
    fprintf(triangles_file, "# %2.d individual's identifier (id)\n", key++);
    fprintf(triangles_file, "# %2.d chromosome or plasmid (c_or_p)\n", key++);
    fprintf(triangles_file, "# %2.d strand\n", key++);
    fprintf(triangles_file, "# %2.d protein position (pos)\n", key++);
    fprintf(triangles_file, "# %2.d length (len)\n", key++);
    fprintf(triangles_file, "# %2.d position of last translated nucleotide (lpos)\n", key++);
    fprintf(triangles_file, "# %2.d primary sequence (sequence)\n", key++);
    fprintf(triangles_file, "# %2.d mean (m)\n", key++);
    fprintf(triangles_file, "# %2.d width (w)\n", key++);
    fprintf(triangles_file, "# %2.d height (h)\n", key++);
    fprintf(triangles_file, "# %2.d concentration (c)\n", key++);
    fprintf(triangles_file, "# %2.d feature (f)\n", key++);
    fprintf(triangles_file, "# %2.d promoter position (prom_pos)\n", key++);
    fprintf(triangles_file, "# %2.d RNA length (rna_len)\n", key++);
    fprintf(triangles_file, "# %2.d basal level (basal_level)\n", key++);
    fprintf(triangles_file, "\n");
    fprintf(triangles_file, "id c_or_p strand pos len lpos sequence m w h c f prom_pos rna_len basal_level\n");
  }
  if (sequence_file_name != NULL)
  {
    sequence_file = fopen(sequence_file_name,"w");
  }

  auto exp_manager = new ExpManager();
  exp_manager->load(num_gener, false, false);


  // The best individual is already known because it is the last in the list
  // Thus we do not need to know anything about the environment and to evaluate the individuals

  // Parse the individuals
  if (best_only)
  {
    Individual* best = exp_manager->get_best_indiv();
    best->do_transcription_translation_folding(); // We need to recompute proteins if not already done (ie if using a population file and not a full backup)
    analyse_indiv(best, triangles_file, sequence_file, gu, best->habitat().phenotypic_target());
  }
  else
  {
    for (const auto& indiv: exp_manager->get_indivs()) {
      indiv->do_transcription_translation_folding(); // We need to recompute proteins if not already done (ie if using a population file and not a full backup)
      analyse_indiv(indiv, triangles_file, sequence_file, gu, indiv->habitat().phenotypic_target());
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

  return EXIT_SUCCESS;
}

// Parsing an individual
inline void analyse_indiv(Individual* indiv, FILE* triangles_file,
                          FILE* sequence_file, int16_t gu,
                          const PhenotypicTarget & phenotypicTarget)
{
  if (gu == -1) // We want to treat all genetic units
  {
    int32_t gen_unit_number = 0;
    for (auto& gen_unit: indiv->get_genetic_unit_list_nonconst()) {
      if (triangles_file != NULL)
      {
        analyse_gu(&gen_unit, gen_unit_number, triangles_file, phenotypicTarget); // We call the triangle parser for each GU successively
      }
      if (sequence_file != NULL)
      {
        const char* dna = gen_unit.get_dna()->data();
        int32_t length = gen_unit.get_dna()->length();
        fprintf(sequence_file,"%.*s ",length,dna); // We output the sequences of each GU separated by a space
      }

      gen_unit_number++;
    }
  }
  else // User specified a genetic unit
  {
    GeneticUnit* gen_unit = &indiv->get_genetic_unit_nonconst(gu);
    if (triangles_file != NULL)
    {
      analyse_gu(gen_unit, gu, triangles_file, phenotypicTarget); // We call the triangle parser
    }
    if (sequence_file != NULL)
    {
      const char* dna = gen_unit->get_dna()->data();
      int32_t length = gen_unit->get_dna()->length();
      fprintf(sequence_file,"%.*s",length,dna); // We output the sequence
    }
  }

  // We go to next line in each file
  if (triangles_file != NULL)
  {
    fprintf(triangles_file,"\n");
  }
  if (sequence_file != NULL)
  {
    fprintf(sequence_file,"\n");
  }
}

// Parsing a GU
inline void analyse_gu(GeneticUnit* gen_unit, int32_t gen_unit_number,
                       FILE* triangles_file,
                       const PhenotypicTarget& phenotypicTarget)
{
  // Construct the list of all rnas
  auto llrnas = gen_unit->get_rna_list();
  auto lrnas = llrnas[LEADING];
  lrnas.splice(lrnas.end(), llrnas[LAGGING]);

  // Parse this list
  int rna_nb = 0;
  for (const auto& rna: lrnas) {
    for (const auto& protein: rna.get_transcribed_proteins()) {
      double mean = protein->get_mean();

      int nfeat = -1;

      // Retrieving the feature of the protein also necessitates the an environment file.
      for (size_t i = 0; i <= static_cast<size_t>(phenotypicTarget.nb_segments()) - 1; ++i)
        if ((mean > phenotypicTarget.segments()[i]->start) and (mean < phenotypicTarget.segments()[i]->stop)) {
          nfeat = phenotypicTarget.segments()[i]->feature;
          break;
        }

      char *dummy;
      fprintf(triangles_file,
              "%" PRId32 " %s %s %" PRId32 " %" PRId32 " %" PRId32
                  " %s %f %f %f %f %d %" PRId32 " %" PRId32 " %f\n",
              gen_unit->get_indiv()->get_id(),
              gen_unit_number != 0 ? "PLASMID" :
              "CHROM  ",
              protein->get_strand() == LEADING ? "LEADING" :
              "LAGGING",
              protein->get_first_translated_pos(),
              protein->get_length(),
              protein->get_last_translated_pos(),
              dummy = protein->get_AA_sequence('_'),
              mean,
              protein->get_width(),
              protein->get_height(),
              protein->get_concentration(),
              nfeat,
              rna.get_promoter_pos(),
              rna.get_transcript_length(),
              rna.get_basal_level());
      delete dummy;
    }
    rna_nb++;
  }
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
  printf("   or :    %s [-r GENER] [-t PHEN_FILE] [-s SEQ_FILE] [-g NUM_GU] [-b]\n", prog_name);
  printf("\nOptions\n");
  printf("  -h, --help\n\tprint this help, then exit\n\n");
  printf("  -V, --version\n\tprint version number, then exit\n\n");
  printf("  -r GENER  :\n");
  printf("\tread generation GENER from a full aevol backup\n\t(default: reads from last_gener.txt)\n");
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
