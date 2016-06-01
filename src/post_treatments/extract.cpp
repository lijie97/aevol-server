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
// ****************************************************************************

// ============================================================================
//                                   Includes
// ============================================================================
#include <getopt.h>
#include <cstdlib>
#include <cstdio>

#include "aevol.h"

using namespace aevol;

// Command-line option variables
static char* triangles_file_name  = nullptr;
static char* sequence_file_name  = nullptr;
static bool best_only = true;
static int16_t gu = -1;
static int32_t timestep = -1;

// Helper functions
void print_help(char* prog_path);
void interpret_cmd_line_options(int argc, char* argv[]);

void analyse_indiv(Individual* indiv, FILE* triangles_file, FILE* sequence_file,
                   int16_t gu, const PhenotypicTarget& phenotypicTarget);

void analyse_gu(GeneticUnit* gen_unit, int32_t gen_unit_number,
                FILE* triangles_file, const PhenotypicTarget& phenotypicTarget);




int main(int argc, char* argv[]) {
  interpret_cmd_line_options(argc, argv);

  // Open the files
  FILE* triangles_file = nullptr;
  FILE* sequence_file = nullptr;

  if (triangles_file_name != nullptr) {
    triangles_file = fopen(triangles_file_name, "w");

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
    fprintf(triangles_file,
            "id c_or_p strand pos len lpos sequence m w h c f "
                "prom_pos rna_len basal_level\n");
  }
  if (sequence_file_name != nullptr) {
    sequence_file = fopen(sequence_file_name,"w");
  }

  auto exp_manager = new ExpManager();
  exp_manager->load(timestep, false, false);


  // The best individual is already known because it is the last in the list
  // Thus we do not need to know anything about the environment and to evaluate
  // the individuals

  // Parse the individuals
  if (best_only) {
    Individual* best = exp_manager->best_indiv();
    analyse_indiv(best, triangles_file, sequence_file, gu,
                  best->habitat().phenotypic_target());
  }
  else {
    for (const auto& indiv: exp_manager->indivs()) {
      analyse_indiv(indiv, triangles_file, sequence_file, gu,
                    indiv->habitat().phenotypic_target());
    }
  }

  if (sequence_file_name != nullptr) {
    fclose(sequence_file);
  }
  if (triangles_file_name != nullptr) {
    fclose(triangles_file);
  }

  delete [] triangles_file_name;
  delete [] sequence_file_name;

  delete exp_manager;

  return EXIT_SUCCESS;
}

// Parsing an individual
inline void analyse_indiv(Individual* indiv, FILE* triangles_file,
                          FILE* sequence_file, int16_t gu,
                          const PhenotypicTarget & phenotypicTarget) {
  if (gu == -1) { // We want to treat all genetic units
    int32_t gen_unit_number = 0;
    for (auto& gen_unit: indiv->genetic_unit_list_nonconst()) {
      if (triangles_file != nullptr) {
        analyse_gu(&gen_unit, gen_unit_number, triangles_file,
                   phenotypicTarget);
      }
      if (sequence_file != nullptr) {
        // The sequences of different GUs are separated by a space
        if (gen_unit_number > 0) fprintf(sequence_file, " ");

        const char* dna = gen_unit.dna()->data();
        int32_t length = gen_unit.dna()->length();
        fprintf(sequence_file, "%.*s", length, dna);
      }

      gen_unit_number++;
    }
  }
  else { // User has specified a genetic unit
    GeneticUnit* gen_unit = &indiv->genetic_unit_nonconst(gu);
    if (triangles_file != nullptr) {
      analyse_gu(gen_unit, gu, triangles_file, phenotypicTarget);
    }
    if (sequence_file != nullptr) {
      const char* dna = gen_unit->dna()->data();
      int32_t length = gen_unit->dna()->length();
      fprintf(sequence_file, "%.*s", length, dna);
    }
  }

  // We go to next line in each file
  if (triangles_file != nullptr) {
    fprintf(triangles_file, "\n");
  }
  if (sequence_file != nullptr) {
    fprintf(sequence_file, "\n");
  }
}

// Parsing a GU
inline void analyse_gu(GeneticUnit* gen_unit, int32_t gen_unit_number,
                       FILE* triangles_file,
                       const PhenotypicTarget& phenotypicTarget) {
  // Construct the list of all rnas
  auto llrnas = gen_unit->rna_list();
  auto lrnas = llrnas[LEADING];
  lrnas.splice(lrnas.end(), llrnas[LAGGING]);

  // Parse this list
  int rna_nb = 0;
  for (const auto& rna: lrnas) {
    for (const auto& protein: rna.transcribed_proteins()) {
      double mean = protein->mean();

      int nfeat = -1;
      for (size_t i = 0 ;
           i <= static_cast<size_t>(phenotypicTarget.nb_segments()) - 1 ;
           ++i) {
        if ((mean > phenotypicTarget.segments()[i]->start) and
            (mean < phenotypicTarget.segments()[i]->stop)) {
          nfeat = phenotypicTarget.segments()[i]->feature;
          break;
        }
      }

      char *dummy;
      fprintf(triangles_file,
              "%" PRId32 " %s %s %" PRId32 " %" PRId32 " %" PRId32
                  " %s %f %f %f %f %d %" PRId32 " %" PRId32 " %f\n",
              gen_unit->indiv()->id(),
              gen_unit_number != 0 ? "PLASMID" :
              "CHROM  ",
              protein->strand() == LEADING ? "LEADING" :
              "LAGGING",
              protein->first_translated_pos(),
              protein->length(),
              protein->last_translated_pos(),
              dummy = protein->AA_sequence('_'),
              mean,
              protein->width(),
              protein->height(),
              protein->concentration(),
              nfeat,
              rna.promoter_pos(),
              rna.transcript_length(),
              rna.basal_level());
      delete dummy;
    }
    rna_nb++;
  }
}


void print_help(char* prog_path) {
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
  printf("%s:\n", prog_name);
  printf("\tExtracts the genotype and/or data about the phenotype of individuals\n");
  printf("\tin the provided population and write them into text files easy to parse\n");
  printf("\twith e.g. matlab.\n");
  printf("\n");
  printf("Usage : %s -h\n", prog_name);
  printf("   or : %s -V or --version\n", prog_name);
  printf("   or : %s [-t TIMESTEP] [-S SEQ_FILE] [-T TRIANGLE_FILE] [-U NUM_GU] [-a]\n",
         prog_name);
  printf("\nOptions\n");
  printf("  -h, --help\n\tprint this help, then exit\n");
  printf("  -V, --version\n\tprint version number, then exit\n");
  printf("  -t TIMESTEP\n");
  printf("\tspecify timestep of the individual(s) of interest\n");
  printf("  -S SEQ_FILE\n");
  printf("\textract sequences into file SEQ_FILE\n");
  printf("  -T TRIANGLE_FILE\n");
  printf("\textract phenotypic data into file TRIANGLE_FILE\n");
  printf("  -U NUM_GU\n");
  printf("\tonly treat genetic unit #NUM_GU (default: treat all genetic units)\n");
  printf("  -a\n");
  printf("\ttreat all the individuals (default: treat only the best)\n");
}


void interpret_cmd_line_options(int argc, char* argv[]) {
  // Define allowed options
  const char * options_list = "hVt:aU:S:T:";
  static struct option long_options_list[] = {
      {"help",      no_argument,        nullptr, 'h'},
      {"version",   no_argument,        nullptr, 'V'},
      {"timestep",  required_argument,  nullptr, 't'},
      {"all",       no_argument,        nullptr, 'a'},
      {"gu",        required_argument,  nullptr, 'U'},
      {"sequence",  required_argument,  nullptr, 'S'},
      {"triangles", required_argument,  nullptr, 'T'},
      {0, 0, 0, 0}
  };

  // Get actual values of the command-line options
  int option;
  while ((option = getopt_long(argc, argv, options_list,
                               long_options_list, nullptr)) != -1) {
    switch (option) {
      case 'h' : {
        print_help(argv[0]);
        exit(EXIT_SUCCESS);
      }
      case 'V' : {
        Utils::PrintAevolVersion();
        exit(EXIT_SUCCESS);
      }
      case 't' : {
        timestep = atol(optarg);
        break;
      }
      case 'a' : {
        best_only = false;
        break;
      }
      case 'U' : {
        gu = atoi(optarg);
        break;
      }
      case 'S' : {
        sequence_file_name = new char[strlen(optarg) + 1];
        sprintf(sequence_file_name, "%s", optarg);
        break;
      }
      case 'T' : {
        triangles_file_name = new char[strlen(optarg) + 1];
        sprintf(triangles_file_name, "%s", optarg);
        break;
      }
    }
  }

  // If timestep wasn't provided, use default
  if (timestep < 0) {
    timestep = OutputManager::last_gener();
  }

  // If neither the sequence_file_name nor the triangles_file_name was provided,
  // we will output only the sequence in a default-named file
  if (sequence_file_name == nullptr && triangles_file_name == nullptr) {
    sequence_file_name = new char[255];
    strcpy(sequence_file_name, "sequence");
  }
}
