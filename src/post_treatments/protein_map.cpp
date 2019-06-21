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

// The input file is produced by the lineage post-treatment, please refer to it
// for e.g. the file format/content

// ============================================================================
//                                   Includes
// ============================================================================
#include <cinttypes>
#include <getopt.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <zlib.h>
#include <err.h>
#include <cerrno>
#include <sys/stat.h>
#include <unistd.h>
#include <list>
#include <iostream>
#include <fstream>

#include "aevol.h"

using namespace aevol;

class ProteinMap {
public:
    ProteinMap(int nb_prot, int dna_length) {
        nb_prot_ = nb_prot;
        start_pos_.resize(nb_prot_);
        length_.resize(nb_prot_);
        basal_level_.resize(nb_prot_);
        hamming_dist_.resize(nb_prot_);
        dist_next_prot_.resize(nb_prot_);
        dna_length_ = dna_length;
    }

    void add_protein(int32_t start_pos, int32_t length, int32_t basal_level, int32_t hamming_dist, int32_t dist_next_prot) {
        start_pos_[cpt_] = start_pos;
        length_[cpt_] = length;
        basal_level_[cpt_] = basal_level;
        hamming_dist_[cpt_] = hamming_dist;
        dist_next_prot_[cpt_] = dist_next_prot;
        cpt_++;
    }

    std::vector<int32_t> start_pos_;
    std::vector<int32_t> length_;
    std::vector<double> basal_level_;
    std::vector<int32_t> hamming_dist_;
    std::vector<int32_t> dist_next_prot_;
    int nb_prot_;
    int cpt_ = 0;
    int dna_length_;
};

// Helper functions
void interpret_cmd_line_options(int argc, char* argv[]);
void print_help(char* prog_path);
ProteinMap* compute_protein_map(Individual* indiv);

// Command-line option variables
static char* lineage_file_name = nullptr;
static bool verbose = false;

static long pt_begin = -1;
static long pt_end = -1;

int main(int argc, char* argv[]) {
  interpret_cmd_line_options(argc, argv);

  printf("\n"
         "WARNING : Parameter change during simulation is not managed in general.\n"
         "          Only changes in environmental target done with aevol_modify are handled.\n"
         "\n");

  // =======================
  //  Open the lineage file
  // =======================
  gzFile lineage_file = gzopen(lineage_file_name, "r");
  if (lineage_file == Z_NULL) {
    fprintf(stderr, "ERROR : Could not read the lineage file %s\n", lineage_file_name);
    exit(EXIT_FAILURE);
  }

  int64_t t0 = 0;
  int64_t t_end = 0;
  int32_t final_indiv_index = 0;
  int32_t final_indiv_rank  = 0;

  gzread(lineage_file, &t0, sizeof(t0));
  gzread(lineage_file, &t_end, sizeof(t_end));
  gzread(lineage_file, &final_indiv_index, sizeof(final_indiv_index));
  gzread(lineage_file, &final_indiv_rank,  sizeof(final_indiv_rank));

  if (pt_begin == -1) pt_begin = t0;
  if (pt_end == -1) pt_end = t_end;

  if (verbose) {
    printf("\n\n""===============================================================================\n");
    printf(" Statistics of the ancestors of indiv. %" PRId32
           " (rank %" PRId32 ") from time %" PRId64 " to %" PRId64 "\n",
           final_indiv_index, final_indiv_rank, t0, t_end);
    printf("================================================================================\n");
  }



  // =============================
  //  Open the experiment manager
  // =============================
  ExpManager* exp_manager = new ExpManager();
  exp_manager->load(t0, true, false);

  // The current version doesn't allow for phenotypic variation nor for
  // different phenotypic targets among the grid
  if (not exp_manager->world()->phenotypic_target_shared())
    Utils::ExitWithUsrMsg("sorry, ancestor stats has not yet been implemented "
                              "for per grid-cell phenotypic target");
  auto phenotypicTargetHandler =
      exp_manager->world()->phenotypic_target_handler();
  if (not (phenotypicTargetHandler->var_method() == NO_VAR))
    Utils::ExitWithUsrMsg("sorry, ancestor stats has not yet been implemented "
                              "for variable phenotypic targets");

  int64_t backup_step = exp_manager->backup_step();


  // =========================
  //  Open the output file(s)
  // =========================
  // Create missing directories
  int status;
  status = mkdir("stats/ancestor_stats/", 0755);
  if ((status == -1) && (errno != EEXIST)) {
    err(EXIT_FAILURE, "stats/ancestor_stats/");
  }

    // =========================
    //  Create data structure
    // =========================
    ProteinMap** list_prot_map = new ProteinMap*[t_end-t0+1];

    // ==================================================
  //  Prepare the initial ancestor and write its stats
  // ==================================================
  GridCell* grid_cell = new GridCell(lineage_file, exp_manager, nullptr);
  auto* indiv = grid_cell->individual();
  indiv->Evaluate();
  indiv->compute_statistical_data();
  indiv->compute_non_coding();

  list_prot_map[time()] = compute_protein_map(indiv);

  // ==========================================================================
  //  Replay the mutations to get the successive ancestors and analyze them
  // ==========================================================================
  ReplicationReport* rep = nullptr;
  int32_t index;
  ExpManager* exp_manager_backup = nullptr;
  int32_t unitlen_before;
  double metabolic_error_before;
  double impact_on_metabolic_error;
  char mut_descr_string[255];


  aevol::AeTime::plusplus();
  while (time() <= t_end)
  {
    rep = new ReplicationReport(lineage_file, indiv);
    index = rep->id(); // who we are building...

    if (verbose)
        printf("Rebuilding ancestor at generation %" PRId64
            " (index %" PRId32 ")...", time(), index);

    indiv->Reevaluate();

    // 2) Replay replication (create current individual's child)
    GeneticUnit& gen_unit = indiv->genetic_unit_nonconst(0);
    GeneticUnit* stored_gen_unit = nullptr;
    Individual* stored_indiv = nullptr;

    // For each genetic unit, replay the replication (undergo all mutations)
    // TODO <david.parsons@inria.fr> disabled for multiple GUs
    const auto& dnarep = rep->dna_replic_report();

    // TODO(dpa) The following 3 for loops should be factorized.
    // However, this is not as easy as it sounds :-D
    // see std::list::splice
    for (const auto& mut: dnarep.HT())
      gen_unit.dna()->undergo_this_mutation(*mut);

    for (const auto& mut: dnarep.rearrangements()) {
      // Apply mutation
      gen_unit.dna()->undergo_this_mutation(*mut);
    }

    for (const auto& mut: dnarep.mutations()) {
      // Apply mutation
      gen_unit.dna()->undergo_this_mutation(*mut);

    }

    // 3) All the mutations have been replayed, we can now evaluate the new individual
    indiv->Reevaluate();
    indiv->compute_statistical_data();
    indiv->compute_non_coding();

    list_prot_map[time()] = compute_protein_map(indiv);

    if (verbose) printf(" OK\n");

    delete rep;

    aevol::AeTime::plusplus();
  }

  gzclose(lineage_file);

  std::ofstream proteins_map_file;
    proteins_map_file.open("proteins_map_gen.csv",std::ofstream::trunc);
    proteins_map_file<<"generation,protein_id,shine_dal,length,concentration,hamming_dist,dist_next_protein,nb_proteins,dna_length"<<std::endl;

  for (int i = t0; i <= t_end; i++) {
      for (int i_prot = 0; i_prot < list_prot_map[i]->nb_prot_; i_prot++) {
          proteins_map_file<<i<<","<<i_prot<<","<<list_prot_map[i]->start_pos_[i_prot]
                  <<","<<list_prot_map[i]->length_[i_prot]
                  <<","<<list_prot_map[i]->basal_level_[i_prot]
                  <<","<<list_prot_map[i]->hamming_dist_[i_prot]
                  <<","<<list_prot_map[i]->dist_next_prot_[i_prot]
                  <<","<<list_prot_map[i]->nb_prot_<<","<<list_prot_map[i]->dna_length_<<std::endl;
      }
  }

    proteins_map_file.flush();
    proteins_map_file.close();

  // Additional outputs

  delete exp_manager;
  delete indiv;

  return EXIT_SUCCESS;
}

void interpret_cmd_line_options(int argc, char* argv[]) {
  // =====================
  //  Parse command line
  // =====================
  const char * short_options = "hVb:e:v";
  static struct option long_options[] = {
    {"help",                no_argument, NULL, 'h'},
    {"version",             no_argument, NULL, 'V'},
    {"verbose",             no_argument, NULL, 'v'},
    {"begin",     required_argument, nullptr, 'b'},
    {"end",       required_argument, nullptr, 'e'},
    {0, 0, 0, 0}
  };

  int option;
  while((option = getopt_long(argc, argv, short_options,
                              long_options, nullptr)) != -1) {
    switch(option) {
      case 'h':
        print_help(argv[0]);
        exit(EXIT_SUCCESS);
      case 'V':
        Utils::PrintAevolVersion();
        exit(EXIT_SUCCESS);
      case 'v':
        verbose = true;
        break;
      case 'b':
            pt_begin = atol(optarg);
            break;
      case 'e':
            pt_end = atol(optarg);
            break;
      default:
        // An error message is printed in getopt_long, we just need to exit
        exit(EXIT_FAILURE);
    }
  }

  // There should be only one remaining arg: the lineage file
  if (optind != argc - 1) {
    Utils::ExitWithUsrMsg("please specify a lineage file");
  }

  lineage_file_name = new char[strlen(argv[optind]) + 1];
  sprintf(lineage_file_name, "%s", argv[optind]);
}

void print_help(char* prog_path) {
  // Get the program file-name in prog_name (strip prog_path of the path)
  char* prog_name; // No new, it will point to somewhere inside prog_path
  if ((prog_name = strrchr(prog_path, '/'))) {
    prog_name++;
  }
  else {
    prog_name = prog_path;
  }

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
  printf("%s: create an experiment with setup as specified in PARAM_FILE.\n",
  prog_name);
  printf("\n");
  printf("Usage : %s -h or --help\n", prog_name);
  printf("   or : %s -V or --version\n", prog_name);
  printf("   or : %s LINEAGE_FILE [-FMv]\n",
  prog_name);
  printf("\nOptions\n");
  printf("  -h, --help\n\tprint this help, then exit\n");
  printf("  -V, --version\n\tprint version number, then exit\n");
  printf("  -v, --verbose\n\tbe verbose\n");
}

ProteinMap* compute_protein_map(Individual* indiv) {
    ProteinMap* pmap = new ProteinMap(indiv->protein_list().size(),indiv->amount_of_dna());

    // Make a copy of each genetic unit's protein list
    for (auto& gen_unit: indiv->genetic_unit_list_nonconst()) {
        // append all proteins from `gen_unit` to `protein_list_`
        for (auto& strand_id: {LEADING, LAGGING}) {
            auto& strand = gen_unit.protein_list(strand_id);
            int pos_next = std::prev(gen_unit.protein_list(strand_id).end())->shine_dal_pos();
            bool first = true;
            int pos_first = gen_unit.protein_list(strand_id).begin()->shine_dal_pos();
            int pos_prev = -1;
            for (auto& p: strand) {
                int dist = -1;
                if (first) {
                    if (strand_id == LEADING)
                        dist = p.shine_dal_pos() + (indiv->amount_of_dna() - pos_next);
                    else
                        dist = (indiv->amount_of_dna() - p.shine_dal_pos()) + pos_next;

                    first = false;
                } else if (p.shine_dal_pos() == std::prev(gen_unit.protein_list(strand_id).end())->shine_dal_pos()) {
                    if (strand_id == LEADING)
                        dist = (indiv->amount_of_dna() - p.shine_dal_pos()) + pos_first;
                    else
                        dist = p.shine_dal_pos() + (indiv->amount_of_dna() - pos_first);
                } else {
                    if (strand_id == LEADING)
                        dist = p.shine_dal_pos() - pos_prev;
                    else
                        dist = pos_prev - p.shine_dal_pos();
                }
                pos_prev = p.shine_dal_pos();
                int8_t prom_dist;
                gen_unit.is_promoter(LEADING, (*p.rna_list().begin())->promoter_pos(),
                                                      prom_dist);

                pmap->add_protein(p.shine_dal_pos(),p.length(),p.concentration(),prom_dist,dist);
            }
        }
    }

    return pmap;
}