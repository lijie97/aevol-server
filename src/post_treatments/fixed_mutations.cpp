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

/// The input file (lineage.ae or lineage.rae) must contain the following information:
///
/// - common data                                                (ae_common::write_to_backup)
/// - begin gener                                                (int32_t)
/// - end gener                                                  (int32_t)
/// - final individual index                                     (int32_t)
/// - initial genome size                                        (int32_t)
/// - initial ancestor (nb genetic units + sequences)            (ae_individual::write_to_backup)
/// - replication report of ancestor at generation begin_gener+1 (ae_replic_report::write_to_backup)
/// - replication report of ancestor at generation begin_gener+2 (ae_replic_report::write_to_backup)
/// - replication report of ancestor at generation begin_gener+3 (ae_replic_report::write_to_backup)
/// - ...
/// - replication report of ancestor at generation end_gener     (ae_replic_report::write_to_backup)

#include <list>
#include <cstdlib>
#include <string>
#include <cinttypes>
#include <cstring>
#include <iostream>
#include <exception>
#include <getopt.h>

#include "aevol.h"
#include "GzHelpers.h"

using std::string;
using std::list;
using std::cout;
using std::cerr;
using std::endl;
using std::exit;
using std::ofstream;
using std::to_string;

using namespace aevol;

/// At each generation, the mutations are recorded on the one hand
/// and the resulting generation is recorded on the other hand.
/// This enum determines when both recordings match.
enum class Check {
  FULL,  ///< at every generation: fail early, for debugging purposes
  LIGHT, ///< only at the last generation: sanity check, for when things should go well
  ENV,   ///< when environment changed
  NONE   ///< never: for the believers
};

void print_help(char* prog_path);

int main(int argc, char** argv) {
  // =====================
  //  Parse command line
  // =====================

  // Default values
  Check check(Check::LIGHT);
  string lineage_file_name;
  bool verbose(false);
  double tolerance(0.0);

  const char* short_options = "hVvncf:t:";
  static struct option long_options[] = {
    {"help",      no_argument,       nullptr, 'h'},
    {"version",   no_argument,       nullptr, 'V'},
    {"verbose",   no_argument,       nullptr, 'v'},
    {"nocheck",   no_argument,       nullptr, 'n'},
    {"fullcheck", no_argument,       nullptr, 'c'},
    {"file",      required_argument, nullptr, 'f'},
    {"tolerance", required_argument, nullptr, 't'},
    {0, 0, 0, 0}
  };

  while((int option = getopt_long(argc, argv, short_options, long_options, nullptr)) != -1) {
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
      case 'n':
        check = Check::NONE;
        break;
      case 'c':
        check = Check::FULL;
        break;
      case 'f':
        if (strcmp(optarg, "") == 0)
          throw std::invalid_argument("Option -f or --file missing argument.");
        lineage_file_name = optarg;
        break;
      case 't':
        if (strcmp(optarg, "") == 0)
          throw std::invalid_argument("Option -t or --tolerance missing argument.");
        check = Check::ENV;
        tolerance = atof(optarg);
        break;
      default:
        cerr << "ERROR: Unknown option provided." << endl;
        print_help(argv[0]);
        exit(EXIT_FAILURE);
    }
  }

  if (lineage_file_name.empty()) {
    cerr << "ERROR : Option -f or --file missing." << endl;
    exit(EXIT_FAILURE);
  }

  cout << "\n"
          "WARNING : Parameter change during simulation is not managed in general.\n"
          "          Only changes in environmental target done with aevol_modify are handled.\n"
          "\n";

  // =======================
  //  Open the lineage file
  // =======================

  gzFile lineage_file = gzopen(lineage_file_name.c_str(), "r");
  if (lineage_file == Z_NULL) {
    cerr << "ERROR: Could not read the lineage file " << lineage_file_name << endl;
    exit(EXIT_FAILURE);
  }

  int64_t t0;
  int64_t t_end;
  int32_t final_index;
  int32_t final_indiv_rank;
  gzread(lineage_file, t0, t_end, final_index, final_indiv_rank);

  if (verbose)
    cout << "\n\n"
            "================================================================================\n"
            " Statistics of the ancestors of indiv. #" << final_index << " (t=" << t0 << " to " << t_end << ")\n"
            "================================================================================\n";

  // =========================
  //  Open the experience manager
  // =========================

  // Open the experiment manager
  auto exp_manager = new ExpManager();
  exp_manager->load(t0, true, false);
  Environment* env = new Environment(*(exp_manager->env())); // independent copy

  int64_t backup_step = exp_manager->backup_step();


  // =========================
  //  Open the output file(s)
  // =========================
  char output_file_name[60];
  snprintf(output_file_name,
           60,
           "stats/fixedmut-b%06" PRId64 "-e%06" PRId64 "-i%" PRId32 "-r%" PRId32 ".out",
           t0, t_end, final_index, final_indiv_rank);

  ofstream output(output_file_name);
  if (not output.is_open()) {
    cerr << "ERROR: Could not create the output file " << output_file_name << endl;
    exit(EXIT_FAILURE);
  }


  // Write the header
  fprintf(output, "# #################################################################\n");
  fprintf(output, "#  Mutations in the lineage of the best indiv at generation %" PRId64 "\n", t_end);
  fprintf(output, "# #################################################################\n");
  fprintf(output, "#  1.  Generation       (mut. occurred when producing the indiv. of this generation)\n");
  fprintf(output, "#  2.  Genetic unit     (which underwent the mutation, 0 = chromosome) \n");
  fprintf(output, "#  3.  Mutation type    (0: switch, 1: smallins, 2: smalldel, 3:dupl, 4: del, 5:trans, 6:inv, 7:insert, 8:ins_HT, 9:repl_HT) \n");
  fprintf(output, "#  4.  pos_0            (position for the small events, begin_segment for the rearrangements, begin_segment of the inserted segment for ins_HT, begin_segment of replaced segment for repl_HT) \n");
  fprintf(output, "#  5.  pos_1            (-1 for the small events, end_segment for the rearrangements, end_segment of the inserted segment for ins_HT, begin_segment of donor segment for repl_HT) \n");
  fprintf(output, "#  6.  pos_2            (reinsertion point for duplic., cutting point in segment for transloc., insertion point in the receiver for ins_HT, end_segment of the replaced segment for repl_HT, -1 for other events)\n");
  fprintf(output, "#  7.  pos_3            (reinsertion point for transloc., breakpoint in the donor for ins_HT, end_segment of the donor segment for repl_HT, -1 for other events)\n");
  fprintf(output, "#  8.  invert           (transloc, was the segment inverted (0/1)?, sense of insertion for ins_HT (0=DIRECT, 1=INDIRECT), sense of the donor segment for repl_HT (0=DIRECT, 1=INDIRECT),-1 for other events)\n");
  fprintf(output, "#  9.  align_score      (score that was needed for the rearrangement to occur, score of the first alignment for ins_HT and repl_HT)\n");
  fprintf(output, "#  10. align_score2     (score for the reinsertion for transloc, score of the second alignment for ins_HT and repl_HT)\n");
  fprintf(output, "#  11. seg_len          (segment length for rearrangement, donor segment length for ins_HT and repl_HT)\n");
  fprintf(output, "#  12. repl_seg_len     (replaced segment length for repl_HT, -1 for the others)\n");
  fprintf(output, "#  13. GU_length        (before the event)\n");
  fprintf(output, "#  14. Impact of the mutation on the metabolic error (negative value = smaller gap after = beneficial mutation) \n");
  fprintf(output, "#  15. Number of coding RNAs possibly disrupted by the breakpoints \n");
  fprintf(output, "#  16. Number of coding RNAs completely included in the segment (donor segment in the case of a transfer) \n");
  fprintf(output, "#  17. Number of coding RNAs that were completely included in the replaced segment (meaningful only for repl_HT) \n");
  fprintf(output, "####################################################################################################################\n");
  fprintf(output, "#\n");
  fprintf(output, "# Header for R\n");
  fprintf(output, "gener gen_unit mut_type pos_0 pos_1 pos_2 pos_3 invert align_score align_score_2 seg_len repl_seg_len GU_len impact nbgenesatbreak nbgenesinseg nbgenesinreplseg\n");

  // ==============================
  //  Prepare the initial ancestor
  // ==============================

  auto indiv = new Individual(exp_manager, lineage_file);
  indiv->Evaluate(habitat);
  indiv->compute_statistical_data();

  if (verbose)
  {
    printf("Initial fitness     = %f\n", indiv->fitness());
    printf("Initial genome size = %" PRId32 "\n", indiv->total_genome_size());
  }



  // ===============================================================================
  //  Replay the mutation to get the successive ancestors and analyze them
  //  (and, optionally, check that the rebuilt envir and genome are correct each
  //   time a backup is available)
  // ===============================================================================

  Individual* stored_indiv = NULL;
  list<GeneticUnit>::const_iterator stored_unit;

  int32_t index, genetic_unit_number, unitlen_before;
  int32_t nb_genes_at_breakpoints, nb_genes_in_segment, nb_genes_in_replaced_segment;
  double metabolic_error_before, metabolic_error_after, impact_on_metabolic_error;
  char mut_descr_string[80];

  ExpManager* exp_manager_backup = nullptr;
  Environment* backup_env = nullptr;

  for (AeTime::plusplus(); time() <= t_end; AeTime::plusplus()) {
    ReplicationReport report(lineage_file, indiv);
    index = report.id(); // who are we building...
    indiv->set_replication_report(report);

    // Check now?
    bool check_now = ((check == Check::FULL or check == Check::ENV) and Utils::mod(time(), backup_step) == 0) or
                      (check == Check::LIGHT and time() == t_end);


    if ( verbose )
      printf("Rebuilding ancestor at generation %" PRId64 " (index %" PRId32 ")...", get_time(), index);

    env->build();
    env->apply_variation();
    indiv->Reevaluate();

    // Check, and possibly update, the environment according to the backup files (update necessary if the env. was modified by aevol_modify at some point)
    if (Utils::mod(time(), backup_step) == 0) {
        char env_file_name[255];
        sprintf(env_file_name, "./" ENV_FNAME_FORMAT, time());
        gzFile env_file = gzopen(env_file_name, "r");
        backup_env = new Environment();
        backup_env->load(env_file);

        if (! env->is_identical_to(*backup_env, tolerance)) {
            printf("Warning: At t=%" PRId64 ", the replayed environment is not the same\n", time());
            printf("         as the one saved at generation %" PRId64 "... \n", time());
            printf("         with tolerance of %lg\n", tolerance);
            printf("Replacing the replayed environment by the one stored in the backup.\n");
            delete env;
            env = new Environment(*backup_env);
          }
        delete backup_env;
      }


    // Warning: this portion of code won't work if the number of units changes
    // during the evolution, or if some translocations occurred between different genetic units

    genetic_unit_number = 0;
    auto unit = indiv->genetic_unit_list_nonconst().begin();

    if (check_now and Utils::mod(time(), backup_step) == 0) { // TODO vld: redundant with check_now definition?
      exp_manager_backup = new ExpManager();
      exp_manager_backup->load(time(), true, false);
      // TODO: disabled tmp
      // stored_indiv = new ae_individual(* (ae_individual *)exp_manager_backup->indiv_by_id(index), false);
      stored_unit = stored_indiv->genetic_unit_list().begin();
    }

    for (const auto& dnareport: report.dna_replic_report()) {
      assert(unit != indiv->genetic_unit_list().end());

      unit->dna()->set_replic_report(dnareport);

      // because the reports are stored in 3 different member-lists,
      // iterate among the DnaReplicationReport getters:
      // - horizontal transfer events (HT);
      // - rearrangements events (rearrangements);
      // - local  events (point mutations & small indels) (mutations).
      for (auto mutation_getter: dnareport->HT, dnareport->rearrangements, dnareport->mutations)
        for (const auto& mutation: mutation_getter()) {
          metabolic_error_before = indiv->dist_to_target_by_feature(METABOLISM);
          unitlen_before = unit->dna()->length();
          unit->compute_nb_of_affected_genes(&mutation, nb_genes_at_breakpoints, nb_genes_in_segment, nb_genes_in_replaced_segment);


          unit->dna()->undergo_this_mutation(&mutation);
          indiv->reevaluate(env);


          metabolic_error_after = indiv->dist_to_target_by_feature(METABOLISM);
          impact_on_metabolic_error = metabolic_error_after - metabolic_error_before;


          mutation.generic_description_string(mut_descr_string);
          fprintf(output,
                  "%" PRId64 " %" PRId32 " %s %" PRId32 " %.15f  %" PRId32
                    " %" PRId32 " %" PRId32 " \n",
                  time(), genetic_unit_number, mut_descr_string, unitlen_before,
                  impact_on_metabolic_error, nb_genes_at_breakpoints, nb_genes_in_segment, nb_genes_in_replaced_segment);
        }

      if (check_now && Utils::mod(time(), backup_step) == 0) {
        if (verbose) {
          printf("Checking the sequence of the unit...");
          fflush(NULL);
        }

        assert(stored_unit != stored_indiv->genetic_unit_list().end());

        string dna_computed(unit->dna()->data(), static_cast<size_t>(unit->dna()->length()));
        string dna_backuped((stored_unit->dna())->data(), static_cast<size_t>((stored_unit->dna())->length()));

        if (dna_computed == dna_backuped) { // TODO(vld): why does static analyzer state that condition is always true?
          if (verbose)
            printf(" OK\n");
        } else {
          if (verbose)
            printf(" ERROR!\n");
          fprintf(stderr, "Error: the rebuilt unit is not the same as \n");
          fprintf(stderr, "the one saved at generation %" PRId64 "... ", time());
          fprintf(stderr, "Rebuilt unit : %zu bp\n %s\n", strlen(dna_computed), dna_computed.c_str());
          fprintf(stderr, "Stored unit  : %zu bp\n %s\n", strlen(dna_backuped), dna_backuped.c_str());
          gzclose(lineage_file);
          delete indiv;
          delete stored_indiv;
          delete exp_manager_backup;
          delete exp_manager;
          exit(EXIT_FAILURE);
        }

        ++stored_unit;
      }

      ++unit;
      ++genetic_unit_number;
    }

    assert(unit == indiv->genetic_unit_list().end());


    if (verbose)
      printf(" OK\n");

    if (check_now and Utils::mod(get_time(), backup_step) == 0) {
      assert(stored_unit == stored_indiv->get_genetic_unit_list().end());
      delete stored_indiv;
      delete exp_manager_backup;
    }
  }

  gzclose(lineage_file);
  output.close();
  delete exp_manager;
  delete indiv;
  delete env;

  return EXIT_SUCCESS;
}

void print_help(char* prog_path)
{
  printf("\n");
  printf("*********************** aevol - Artificial Evolution ******************* \n");
  printf("*                                                                      * \n");
  printf("*               Fixed mutations post-treatment program                 * \n");
  printf("*                                                                      * \n");
  printf("************************************************************************ \n");
  printf("\n\n");
  printf("This program is Free Software. No Warranty.\n");
  printf("Copyright (C) 2009  LIRIS.\n");
  printf("\n");
#ifdef __REGUL
  printf("Usage : rfixed_mutations -h\n");
  printf("or :    rfixed_mutations [-vn] -f lineage_file \n");
#else
  printf("Usage : fixed_mutations -h\n");
  printf("or :    fixed_mutations [-vn] -f lineage_file \n");
#endif
  printf("\n");
  printf("This program computes the fixed mutations of the individuals within the lineage\n");
  printf("of lineage_file\n");
  printf("\n");
  printf("WARNING: This program should not be used for simulations run with several genetic units\n");
  printf("(eg chromosome + plasmids). It works only for simulations with only a single chromosome.\n");
  printf("\n");
  printf("\t-h or --help       : Display this help.\n");
  printf("\n");
  printf("\t-v or --verbose    : Be verbose, listing generations as they are \n");
  printf("\t                       treated.\n");
  printf("\n");
  printf("\t-n or --nocheck    : Disable genome sequence checking. Makes the \n");
  printf("\t                       program faster, but it is not recommended. \n");
  printf("\t                       It is better to let the program check that \n");
  printf("\t                       when we rebuild the genomes of the ancestors\n");
  printf("\t                       from the lineage file, we get the same sequences\n");
  printf("\t                       as those stored in the backup files.\n");
  printf("\n");
  printf("\t-c or --fullcheck  : Will perform the genome and environment checks every\n");
  printf("\t                       <BACKUP_STEP> generations. Default behaviour is\n");
  printf("\t                       lighter as it only perform sthese checks at the\n");
  printf("\t                       ending generation.\n");
  printf("\n");
  printf("\t-f lineage_file or --file lineage_file : \n");
  printf("\t                       Compute the fixed mutations of the individuals within lineage_file.\n");
  printf("\n");
  printf("\t-t tolerance or --tolerance tolerance : \n");
  printf("\t                       Tolerance used to compare the replayed environment to environment in backup\n");
  printf("\n");
}
