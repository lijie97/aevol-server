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


#include "ExpManager_7.h"

#include "7/Stats_7.h"
#include "Abstract_Metadata.h"
#include "DnaMutator.h"
#include "ExpManager.h"
#include "Promoter.h"
#include "Protein_7.h"
#include "Rna_7.h"
#include "List_Metadata.h"


#include <algorithm>
#include <err.h>
#include <sys/stat.h>

namespace aevol {

#ifndef WITH_STANDALONE_SIMD
bool ExpManager_7::standalone_simd = false;
#else
bool ExpManager_7::standalone_simd = true;
#endif

ExpManager_7::ExpManager_7(ExpManager* exp_m) {
  printf("  Loading SIMD Controller...");

  exp_m_ = exp_m;

  nb_indivs_ = exp_m_->nb_indivs();

  current_individuals       = new Individual_7*[exp_m_->nb_indivs()];
  previous_individuals      = new Individual_7*[exp_m_->nb_indivs()];

  // Allocate Dna_7
  int max_size_dna = -1;
  for (int indiv_id = 0; indiv_id < exp_m_->nb_indivs(); indiv_id++) {
    int x = indiv_id / exp_m_->world()->height();
    int y = indiv_id % exp_m_->world()->height();

    max_size_dna = max_size_dna < exp_m->world()->grid(x, y)->individual()->genetic_unit_seq_length(0) ?
                   exp_m->world()->grid(x, y)->individual()->genetic_unit_seq_length(0) : max_size_dna;
  }

  dna_factory_ = new DnaFactory(DnaFactory_Policy::FIRSTFIT,exp_m_->nb_indivs()*3,max_size_dna);


  for (int indiv_id = 0; indiv_id < exp_m_->nb_indivs(); indiv_id++) {
    int x = indiv_id / exp_m_->world()->height();
    int y = indiv_id % exp_m_->world()->height();

    current_individuals[indiv_id] = new Individual_7(exp_m_, exp_m_->best_indiv()->w_max(),dna_factory_);
    current_individuals[indiv_id]->dna_ = dna_factory_->get_dna(exp_m->world()->grid(x, y)->individual()->genetic_unit_seq_length(0));
    current_individuals[indiv_id]->dna_->set_indiv(exp_m->world()->grid(x, y)->individual()->genetic_unit(0).dna(),dna_factory_);
    current_individuals[indiv_id]->indiv_id = indiv_id;
    current_individuals[indiv_id]->parent_id = indiv_id;
    previous_individuals[indiv_id] = current_individuals[indiv_id];
    current_individuals[indiv_id]->global_id = AeTime::time() * 1024 + indiv_id;
  }

  dna_size = new int[exp_m_->nb_indivs()];
  for (int indiv_id = 0; indiv_id < exp_m_->nb_indivs(); indiv_id++) {
    dna_size[indiv_id] = current_individuals[indiv_id]->dna_->length();
  }

#ifdef __REGUL
  phenotypic_target_handler_ = new SIMD_PhenotypicTargetHandler_R(
        dynamic_cast<PhenotypicTargetHandler_R*>(exp_m->world()->phenotypic_target_handler()),
        exp_m->exp_s());
#else
#ifdef PHENOTYPE_VECTOR
  target = new double[PHENOTYPE_VECTOR_SIZE];
    for (int i = 0; i < PHENOTYPE_VECTOR_SIZE; i++) {
        double tmp = ((Fuzzy *) exp_m->world()->phenotypic_target_handler()->phenotypic_target().fuzzy())->y(((double)i)/D_PHENOTYPE_VECTOR_SIZE);
        target[i] = tmp;
    }
#else
  target = new Vector_Fuzzy(*(Fuzzy*)(exp_m->world()->phenotypic_target_handler()->phenotypic_target().fuzzy()));
#endif
#endif

#ifdef WITH_PERF_TRACES
  std::ofstream perf_traces_file_;
        perf_traces_file_.open("simd_perf_traces.csv", std::ofstream::trunc);
        perf_traces_file_ << "Generation,Indiv_ID,Runtime" << std::endl;
        perf_traces_file_.close();
#endif

  printf(" OK\n");

}


void ExpManager_7::selection(int indiv_id) {
  int32_t selection_scope_x = exp_m_->sel()->selection_scope_x();
  int32_t selection_scope_y = exp_m_->sel()->selection_scope_y();

  int32_t grid_width = exp_m_->grid_width();
  int32_t grid_height = exp_m_->grid_height();

  int16_t neighborhood_size = selection_scope_x * selection_scope_y;

  FitnessFunction fitness_function = exp_m_->sel()->fitness_func();
  int32_t fitness_function_scope_x = exp_m_->sel()->fitness_function_scope_x();
  int32_t fitness_function_scope_y = exp_m_->sel()->fitness_function_scope_y();

  double *  local_fit_array   = new double[neighborhood_size];
  double *  local_meta_array   = new double[neighborhood_size];
  double *  probs             = new double[neighborhood_size];
  int16_t   count             = 0;
  double    sum_local_fit     = 0.0;

  int32_t x = indiv_id / grid_height;
  int32_t y = indiv_id % grid_height;

  int cur_x,cur_y;

#ifdef __REGUL
  double ** fitness_sum_local_tab_;
#endif

  if (fitness_function == FITNESS_LOCAL_SUM) {
#ifdef __REGUL
    fitness_sum_local_tab_ = new double*[fitness_function_scope_x*fitness_function_scope_y];
        for (int tab_id = 0; tab_id < fitness_function_scope_x*fitness_function_scope_y; tab_id++)
          fitness_sum_local_tab_[tab_id] = new double[phenotypic_target_handler_->nb_env_];

        for (int env_id = 0; env_id < phenotypic_target_handler_->nb_env_; env_id++) {

          int tab_id = 0;
          fitness_sum_local_tab_[tab_id][env_id] = 0;

          for (int8_t i = -1; i < selection_scope_x - 1; i++) {
            for (int8_t j = -1; j < selection_scope_y - 1; j++) {
              cur_x = (x + i + grid_width) % grid_width;
              cur_y = (y + j + grid_height) % grid_height;

              int16_t new_x,new_y;
              for (int8_t ii = -1; ii < fitness_function_scope_x - 1; ii++) {
                for (int8_t jj = -1; jj < fitness_function_scope_y - 1; jj++) {
                  //TODO: Check values HERE !

                  new_x = (cur_x + ii + grid_width) % grid_width;
                  new_y = (cur_y + jj + grid_height) % grid_height;

                  fitness_sum_local_tab_[tab_id][env_id] +=
                      previous_individuals[new_x * grid_height + new_y]->fitness_by_env_id_[env_id];
                }
              }
              tab_id++;
            }
          }
        }
#else
    printf("Fitness local sum is not supported for Aevol (only R-Aevol)\n");
    exit(-1);
#endif
  }

  int tab_id = 0;

  for (int8_t i = -1 ; i < selection_scope_x-1 ; i++) {
    for (int8_t j = -1; j < selection_scope_y - 1; j++) {
      cur_x = (x + i + grid_width) % grid_width;
      cur_y = (y + j + grid_height) % grid_height;

      if (fitness_function == FITNESS_EXP) {
        local_fit_array[count] =
            previous_individuals[cur_x * grid_height + cur_y]->fitness;
        local_meta_array[count] =
            previous_individuals[cur_x * grid_height + cur_y]
                ->metaerror;
      } else if (fitness_function == FITNESS_GLOBAL_SUM) {
#ifdef __REGUL
        double composed_fitness = 0;
            for (int env_id = 0; env_id < phenotypic_target_handler_->nb_env_;
                 env_id++) {
              composed_fitness +=
                  previous_individuals[cur_x * grid_height + cur_y]->fitness_by_env_id_[env_id] /
                  fitness_sum_tab_[env_id];
            }
            composed_fitness /= phenotypic_target_handler_->nb_env_;
            local_fit_array[count] = composed_fitness;
#else
        printf("Fitness local sum is not supported for Aevol (only R-Aevol)\n");
        exit(-1);
#endif
      } else if (fitness_function == FITNESS_LOCAL_SUM) {
#ifdef __REGUL
        double composed_fitness = 0;
            for (int env_id = 0; env_id < phenotypic_target_handler_->nb_env_;
                 env_id++) {
              composed_fitness +=
                  previous_individuals[cur_x * grid_height + cur_y]->fitness_by_env_id_[env_id] /
                  fitness_sum_local_tab_[tab_id][env_id];
            }
            composed_fitness /= phenotypic_target_handler_->nb_env_;
            local_fit_array[count] = composed_fitness;
#else
        printf("Fitness local sum is not supported for Aevol (only R-Aevol)\n");
        exit(-1);
#endif
      }

      sum_local_fit += local_fit_array[count];

      count++;
      tab_id++;
    }
  }


  if (fitness_function == FITNESS_LOCAL_SUM) {
#ifdef __REGUL
    for (int tab_id = 0; tab_id < fitness_function_scope_x * fitness_function_scope_y; tab_id++)
          delete[] fitness_sum_local_tab_[tab_id];
        delete[] fitness_sum_local_tab_;
#else
    printf("Fitness local sum is not supported for Aevol (only R-Aevol)\n");
    exit(-1);
#endif
  }

  for(int16_t i = 0 ; i < neighborhood_size ; i++) {
    probs[i] = local_fit_array[i]/sum_local_fit;
  }

  int16_t found_org = exp_m_->world()->grid(x,y)->reprod_prng_simd_->roulette_random(probs, neighborhood_size);

  int16_t x_offset = (found_org / selection_scope_x) - 1;
  int16_t y_offset = (found_org % selection_scope_y) - 1;

  delete [] local_fit_array;
  delete [] local_meta_array;
  delete [] probs;


  exp_m_->next_generation_reproducer_[indiv_id] = ((x+x_offset+grid_width)  % grid_width)*grid_height+
                                          ((y+y_offset+grid_height) % grid_height);
  dna_size[indiv_id] =
      previous_individuals[cur_x*grid_height+cur_y]->dna_->length();
}

void ExpManager_7::check_selection(int indiv_id) {
  int32_t selection_scope_x = exp_m_->sel()->selection_scope_x();
  int32_t selection_scope_y = exp_m_->sel()->selection_scope_y();

  int32_t grid_width = exp_m_->grid_width();
  int32_t grid_height = exp_m_->grid_height();

  int16_t neighborhood_size = selection_scope_x * selection_scope_y;

  double *  local_fit_array   = new double[neighborhood_size];
  double *  local_meta_array   = new double[neighborhood_size];
  double *  probs             = new double[neighborhood_size];
  int16_t   count             = 0;
  double    sum_local_fit     = 0.0;

  int * indiv_index = new int[neighborhood_size];

  int32_t x = indiv_id / grid_height;
  int32_t y = indiv_id % grid_height;

  int cur_x,cur_y;

  for (int8_t i = -1 ; i < selection_scope_x-1 ; i++) {
    for (int8_t j = -1; j < selection_scope_y - 1; j++) {
      cur_x = (x + i + grid_width) % grid_width;
      cur_y = (y + j + grid_height) % grid_height;


      local_fit_array[count] =
          previous_individuals[cur_x * grid_height + cur_y]->fitness;

      local_meta_array[count] =
          previous_individuals[cur_x * grid_height + cur_y]->metaerror;

      if (local_meta_array[count] !=
          previous_individuals[cur_x * grid_height + cur_y]->metaerror) {
        printf("NONONNONONONN\n");
        exit(-1);
      }
      sum_local_fit += local_fit_array[count];


      indiv_index[count] = cur_x * grid_height + cur_y;
      count++;
    }
  }

  for(int16_t i = 0 ; i < neighborhood_size ; i++) {
    probs[i] = local_fit_array[i]/sum_local_fit;
  }

  int16_t found_org = exp_m_->world()->grid(x,y)->reprod_prng_simd_->roulette_random(probs, neighborhood_size);

  int16_t x_offset = (found_org / selection_scope_x) - 1;
  int16_t y_offset = (found_org % selection_scope_y) - 1;


  int found_id = ((x+x_offset+grid_width)  % grid_width)*grid_height+
                 ((y+y_offset+grid_height) % grid_height);


  if (found_id != exp_m_->next_generation_reproducer_[indiv_id]) {


    printf("For individual %d: Selection is diff SIMD %d CPU %d (Meta error %f -- %f || Fitness %e -- %e) \n",
           indiv_id, found_id, exp_m_->next_generation_reproducer_[indiv_id],
        previous_individuals[exp_m_->next_generation_reproducer_[indiv_id]]->metaerror,
        previous_individuals[found_id]->metaerror,
        previous_individuals[exp_m_->next_generation_reproducer_[indiv_id]]->fitness,
        previous_individuals[found_id]->fitness);



    for (int i = 0; i < neighborhood_size; i++) {
      if (i==8) {
        int v_x = indiv_index[i] / grid_height;
        int v_y = indiv_index[i] % grid_height;

        printf(
            "A-A-ERROR -- Individual %d (%d,%d): Metaerror (CPU/GPU) : %e/%e || Fitness (CPU/GPU) : %e/%e \n",
            indiv_index[i],v_x,v_y,
            exp_m_->world()->grid(v_x, v_y)->individual()->dist_to_target_by_feature(
                METABOLISM),
               previous_individuals[indiv_index[i]]->metaerror,
            exp_m_->world()->grid(v_x, v_y)->individual()->fitness(),
               previous_individuals[indiv_index[i]]->fitness);

        printf("ID CPU %lld SIMD %d -- PARENT ID CPU %d SIMD %d\n",
               exp_m_->world()->grid(v_x, v_y)->individual()->id(),
               previous_individuals[indiv_index[i]]->indiv_id,
               exp_m_->world()->grid(v_x, v_y)->individual()->parent_id_,
               previous_individuals[indiv_index[i]]->parent_id);

        printf(
            "Nb RNA SIMD/CPU %ud/%ld Protein %ud/%ld\n",
            previous_individuals[indiv_index[i]]->metadata_->rna_count(),
            exp_m_->world()->grid(v_x, v_y)->individual()->rna_list().size(),
            previous_individuals[indiv_index[i]]->metadata_->proteins_count(),
            exp_m_->world()->grid(v_x, v_y)->individual()->protein_list().size());

        int idx = 0;

        for (auto rna : exp_m_->world()->grid(v_x, v_y)->individual()->rna_list()) {
          printf("RNA CPU %d Start %d Stop %d Leading/Lagging %d Length %d\n", idx,
                 rna->promoter_pos(), rna->last_transcribed_pos(), rna->strand(), rna->transcript_length());
          idx++;
        }

        idx = 0;
        for (idx = 0; idx < (int) (previous_individuals[indiv_index[i]]->metadata_->rna_count()); idx++) {
          printf("RNA SIMD %d Start %d Stop %d Leading/Lagging %d Length %d\n", idx,
              previous_individuals[indiv_index[i]]->metadata_->rnas(idx)->begin,
              previous_individuals[indiv_index[i]]->metadata_->rnas(idx)->end,
              previous_individuals[indiv_index[i]]->metadata_->rnas(idx)->leading_lagging,
              previous_individuals[indiv_index[i]]->metadata_->rnas(idx)->length);
        }


      }

      printf("%d -- (Probs %e %e -- Fit Array %e %e -- Sum Fit %e %e\n",i,
             exp_m_->world()->grid(x,y)->probs[i],probs[i],
             exp_m_->world()->grid(x,y)->local_fit_array[i],local_fit_array[i],
             exp_m_->world()->grid(x,y)->sum_local_fit,sum_local_fit);
    }

    exit(-44);
  }



  delete [] local_fit_array;
  delete [] local_meta_array;
  delete [] probs;
}


void ExpManager_7::do_mutation(int indiv_id) {
  if (ExpManager_7::standalone() && !exp_m_->check_simd()) {

    int x = indiv_id / exp_m_->world()->height();
    int y = indiv_id % exp_m_->world()->height();
    delete exp_m_->dna_mutator_array_[indiv_id];

    exp_m_->dna_mutator_array_[indiv_id] = new DnaMutator(
        exp_m_->world()->grid(x, y)->mut_prng(),
        previous_individuals[exp_m_->next_generation_reproducer_[indiv_id]]->dna_->length(),
        exp_m_->exp_s()->mut_params()->duplication_rate(),
        exp_m_->exp_s()->mut_params()->deletion_rate(),
        exp_m_->exp_s()->mut_params()->translocation_rate(),
        exp_m_->exp_s()->mut_params()->inversion_rate(),
        exp_m_->exp_s()->mut_params()->point_mutation_rate(),
        exp_m_->exp_s()->mut_params()->small_insertion_rate(),
        exp_m_->exp_s()->mut_params()->small_deletion_rate(),
        exp_m_->exp_s()->mut_params()->max_indel_size(),
        exp_m_->exp_s()->min_genome_length(),
        exp_m_->exp_s()->max_genome_length(), indiv_id,x,y);
    exp_m_->dna_mutator_array_[indiv_id]->generate_mutations();
  }

  if (exp_m_->dna_mutator_array_[indiv_id]->hasMutate()) {
    current_individuals[indiv_id] =
        new Individual_7(exp_m_, previous_individuals[exp_m_->next_generation_reproducer_[indiv_id]],dna_factory_);

    current_individuals[indiv_id]->global_id = AeTime::time()*1024+indiv_id;
    current_individuals[indiv_id]->indiv_id = indiv_id;
    current_individuals[indiv_id]->parent_id =
        exp_m_->next_generation_reproducer_[indiv_id];
    if (ExpManager_7::standalone() && exp_m_->record_tree()) {
      int x = indiv_id / exp_m_->world()->height();
      int y = indiv_id % exp_m_->world()->height();
      NewIndivEvent *eindiv = new NewIndivEvent(
          current_individuals[indiv_id],
          previous_individuals[exp_m_->next_generation_reproducer_[indiv_id]],
                                                x, y,indiv_id,exp_m_->next_generation_reproducer_[indiv_id]);

      exp_m_->tree()->update_new_indiv(eindiv);
      delete eindiv;
    }

    auto size_before = current_individuals[indiv_id]->dna_->length_;
#ifdef WITH_PERF_TRACES
    auto t_start = std::chrono::steady_clock::now();
#endif
    if (ExpManager_7::standalone() && !exp_m_->check_simd())
      current_individuals[indiv_id]->dna_->apply_mutations_standalone();
    else
      current_individuals[indiv_id]->dna_->apply_mutations();
#ifdef WITH_PERF_TRACES
    auto t_end = std::chrono::steady_clock::now();
                apply_mutation[indiv_id] = t_end.time_since_epoch().count() - t_start.time_since_epoch().count();
#endif
    auto size_after = current_individuals[indiv_id]->dna_->length_;

#pragma omp atomic
    cumulate_size += size_after;

#pragma omp atomic
    cumulate_diff += std::abs(size_after-size_before);
  } else {



#pragma omp atomic
    nb_clones_++;

    int32_t parent_id = exp_m_->next_generation_reproducer_[indiv_id];

    current_individuals[indiv_id] = previous_individuals[parent_id];

    if (ExpManager_7::standalone() && exp_m_->record_tree()) {
      int x = indiv_id / exp_m_->world()->height();
      int y = indiv_id % exp_m_->world()->height();
      NewIndivEvent *eindiv = new NewIndivEvent(
          current_individuals[indiv_id],
          previous_individuals[exp_m_->next_generation_reproducer_[indiv_id]],
                                                x, y,indiv_id,exp_m_->next_generation_reproducer_[indiv_id]);
      exp_m_->tree()->update_new_indiv(eindiv);
      delete eindiv;
    }

#pragma omp atomic
    current_individuals[indiv_id]->usage_count_++;
#ifdef WITH_PERF_TRACES
    apply_mutation[indiv_id] = -1;
#endif
  }

  dna_size[indiv_id] = current_individuals[indiv_id]->dna_->length();
}

ExpManager_7::~ExpManager_7() {

  printf("Destroy SIMD Controller\n");

  delete stats_best;
  delete stats_mean;

  /* No need to delete current_individuals, at the end of a generation all the
   * element of the table are nullptr
   */

  for (int indiv_id = 0; indiv_id < (int) exp_m_->nb_indivs(); indiv_id++) {
      if (previous_individuals[indiv_id]->usage_count_ > 0)
        previous_individuals[indiv_id]->usage_count_--;
      else {
          for (int rn = 0; rn < previous_individuals[indiv_id]->metadata_->rna_count(); rn++) {
            delete previous_individuals[indiv_id]->metadata_->rnas(rn);
          }

          previous_individuals[indiv_id]->metadata_->rnas_clear();
          for (int rn = 0; rn < previous_individuals[indiv_id]->metadata_->proteins_count(); rn++) {
            delete previous_individuals[indiv_id]->metadata_->proteins(rn);
          }
          previous_individuals[indiv_id]->metadata_->proteins_clear();

          delete previous_individuals[indiv_id];
          previous_individuals[indiv_id] = nullptr;
      }

    delete exp_m_->dna_mutator_array_[indiv_id];
  }

  delete[] previous_individuals;
  delete[] current_individuals;

  delete[] dna_size;

  delete dna_factory_;

#ifdef __REGUL
#else
  delete target;
#endif
}


void ExpManager_7::start_stop_RNA(int indiv_id) {
  for (int dna_pos = 0; dna_pos < dna_size[indiv_id]; dna_pos++) {


    int len = dna_size[indiv_id];
    if (len >= PROM_SIZE) {
      int prom_dist_leading[26];
      int prom_dist_lagging[26];

      int term_dist_leading[4];
      int term_dist_lagging[4];

      for (int motif_id = 0; motif_id < 52; motif_id++) {
        if (motif_id >= 26 && motif_id < 48) {
          // LAGGING
          int t_motif_id = motif_id - 26;
          prom_dist_lagging[t_motif_id] =
              PROM_SEQ_LAG[t_motif_id] ==
                      current_individuals[indiv_id]->dna_->data_[
                  dna_pos - t_motif_id < 0 ? len +
                                             dna_pos -
                                             t_motif_id :
                  dna_pos - t_motif_id]
              ? 0 : 1;
        } else if (motif_id < 22) {
          // LEADING
          prom_dist_leading[motif_id] =
              PROM_SEQ_LEAD[motif_id] ==
                      current_individuals[indiv_id]->dna_->data_[
                  dna_pos + motif_id >= len ? dna_pos +
                                              motif_id -
                                              len
                                            : dna_pos +
                                              motif_id]
              ? 0
              : 1;
        } else if (motif_id >= 22 && motif_id < 26) {
          int t_motif_id = motif_id - 22;
          // LEADING
          term_dist_leading[t_motif_id] =
              current_individuals[indiv_id]->dna_->data_[
                  dna_pos + t_motif_id >= len ? dna_pos +
                                                t_motif_id -
                                                len :
                  dna_pos + t_motif_id] !=
                      current_individuals[indiv_id]->dna_->data_[
                  dna_pos - t_motif_id + 10 >= len ?
                  dna_pos - t_motif_id + 10 - len :
                  dna_pos -
                  t_motif_id +
                  10] ? 1
                      : 0;
        } else {
          int t_motif_id = motif_id - 48;
          term_dist_lagging[t_motif_id] =
              current_individuals[indiv_id]->dna_->data_[
                  dna_pos - t_motif_id < 0 ? dna_pos -
                                             t_motif_id +
                                             len
                                           : dna_pos -
                                             t_motif_id] !=
                      current_individuals[indiv_id]->dna_->data_[
                  dna_pos + t_motif_id - 10 < 0 ? dna_pos +
                                                  t_motif_id -
                                                  10 + len
                                                :
                  dna_pos + t_motif_id - 10] ? 1 : 0;
        }
      }


      int dist_lead = prom_dist_leading[0] +
                      prom_dist_leading[1] +
                      prom_dist_leading[2] +
                      prom_dist_leading[3] +
                      prom_dist_leading[4] +
                      prom_dist_leading[5] +
                      prom_dist_leading[6] +
                      prom_dist_leading[7] +
                      prom_dist_leading[8] +
                      prom_dist_leading[9] +
                      prom_dist_leading[10] +
                      prom_dist_leading[11] +
                      prom_dist_leading[12] +
                      prom_dist_leading[13] +
                      prom_dist_leading[14] +
                      prom_dist_leading[15] +
                      prom_dist_leading[16] +
                      prom_dist_leading[17] +
                      prom_dist_leading[18] +
                      prom_dist_leading[19] +
                      prom_dist_leading[20] +
                      prom_dist_leading[21];

      if (dist_lead <= 4) {
        PromoterStruct* nprom = new PromoterStruct(dna_pos, dist_lead,
                                                   true);
        int prom_idx;

        prom_idx = current_individuals[indiv_id]->metadata_->promoter_count();
        current_individuals[indiv_id]->metadata_->set_promoters_count(
            current_individuals[indiv_id]->metadata_->promoter_count()+ 1);

        current_individuals[indiv_id]->metadata_->promoter_add(prom_idx, nprom);

        delete nprom;
      }



      int dist_term_lead = term_dist_leading[0] +
                           term_dist_leading[1] +
                           term_dist_leading[2] +
                           term_dist_leading[3];

      if (dist_term_lead == 4) {
        current_individuals[indiv_id]->metadata_->terminator_add(LEADING, dna_pos);
      }

      int dist_lag = prom_dist_lagging[0] +
                     prom_dist_lagging[1] +
                     prom_dist_lagging[2] +
                     prom_dist_lagging[3] +
                     prom_dist_lagging[4] +
                     prom_dist_lagging[5] +
                     prom_dist_lagging[6] +
                     prom_dist_lagging[7] +
                     prom_dist_lagging[8] +
                     prom_dist_lagging[9] +
                     prom_dist_lagging[10] +
                     prom_dist_lagging[11] +
                     prom_dist_lagging[12] +
                     prom_dist_lagging[13] +
                     prom_dist_lagging[14] +
                     prom_dist_lagging[15] +
                     prom_dist_lagging[16] +
                     prom_dist_lagging[17] +
                     prom_dist_lagging[18] +
                     prom_dist_lagging[19] +
                     prom_dist_lagging[20] +
                     prom_dist_lagging[21];

      if (dist_lag <= 4) {
        PromoterStruct* nprom = new PromoterStruct(dna_pos, dist_lag,
                                                   false);
        int prom_idx;
        prom_idx = current_individuals[indiv_id]->metadata_->promoter_count();
        current_individuals[indiv_id]->metadata_->set_promoters_count(
            current_individuals[indiv_id]->metadata_->promoter_count() + 1);

        current_individuals[indiv_id]->metadata_->promoter_add(prom_idx, nprom);
        delete nprom;
      }

      int dist_term_lag = term_dist_lagging[0] +
                          term_dist_lagging[1] +
                          term_dist_lagging[2] +
                          term_dist_lagging[3];


      if (dist_term_lag == 4) {
        current_individuals[indiv_id]->metadata_->terminator_add(LAGGING, dna_pos);
      }
    }
  }
}


void ExpManager_7::opt_prom_compute_RNA(int indiv_id) {
  if (exp_m_->dna_mutator_array_[indiv_id]->hasMutate()) {
    current_individuals[indiv_id]->metadata_->proteins_clear();
    current_individuals[indiv_id]->metadata_->rnas_clear();
    current_individuals[indiv_id]->metadata_->terminators_clear();
  }
  if (exp_m_->dna_mutator_array_[indiv_id]->hasMutate()) {

    current_individuals[indiv_id]->metadata_->rnas_resize(
        current_individuals[indiv_id]->metadata_->promoter_count());

    current_individuals[indiv_id]->metadata_->promoter_begin();

    for (int prom_idx = 0;
         prom_idx < (int)current_individuals[indiv_id]->metadata_->promoter_count(); prom_idx++) {
      PromoterStruct*prom =
          current_individuals[indiv_id]->metadata_->promoter_next();

      if (prom != nullptr) {
        Dna_7*dna = current_individuals[indiv_id]->dna_;
        int32_t dna_length = dna->length();
        int prom_pos;
        bool lead_lag;
        double prom_error;

        prom_pos = prom->pos;
        lead_lag = prom->leading_or_lagging;
        prom_error = fabs(
            ((float) prom->error));

        if (lead_lag) {
          /* Search for terminators */
          int cur_pos =
              Utils::mod(prom_pos + 22,dna_length);
          int start_pos = cur_pos;

          bool terminator_found = false;
          bool no_terminator = false;
          int term_dist_leading = 0;

          int loop_size = 0;

          while (!terminator_found) {
            for (int t_motif_id = 0; t_motif_id < 4; t_motif_id++) {

              term_dist_leading +=
                  dna->get_lead(cur_pos + t_motif_id) !=
                  dna->get_lead(cur_pos - t_motif_id + 10)
                  ? 1 : 0;
            }

            if (term_dist_leading == 4) {
              terminator_found = true;
            }
            else {
              cur_pos = Utils::mod(cur_pos + 1,dna_length);

              term_dist_leading = 0;
              if (cur_pos == start_pos) {
                no_terminator = true;
                terminator_found = true;
              }
            }

            loop_size++;

          }

          if (!no_terminator) {

            int32_t rna_end =
                cur_pos;
            int32_t rna_length = 0;

            rna_length = loop_size + TERM_SIZE -1;

            if (rna_length > 0) {


              int glob_rna_idx = -1;
              glob_rna_idx =
                  current_individuals[indiv_id]->metadata_->rna_count_++;
              current_individuals[indiv_id]->metadata_->rna_add(glob_rna_idx,
                                                                 prom_pos,
                                                                 rna_end,
                                                                 !lead_lag,
                                                                 1.0 -
                                                                 prom_error /
                                                                 5.0, rna_length);
            }
          }
        } else {
          /* Search for terminator */
          int cur_pos =
              prom_pos - 22;
          cur_pos =
              cur_pos < 0 ? dna_length + (cur_pos) : cur_pos;
          int start_pos = cur_pos;
          bool terminator_found = false;
          bool no_terminator = false;
          int term_dist_lagging = 0;

          int loop_size = 0;

          while (!terminator_found) {
            for (int32_t t_motif_id = 0; t_motif_id < 4; t_motif_id++) {
              term_dist_lagging +=
                  dna->get_lag(cur_pos - t_motif_id) !=
                  dna->get_lag(cur_pos + t_motif_id - 10)
                  ? 1 : 0;
            }

            if (term_dist_lagging == 4) {
              terminator_found = true;
            }
            else {
              cur_pos =
                  cur_pos - 1 < 0 ? dna_length + (cur_pos - 1)
                                  : cur_pos - 1;
              term_dist_lagging = 0;
              if (cur_pos == start_pos) {
                no_terminator = true;
                terminator_found = true;
              }
            }
            loop_size++;
          }

          if (!no_terminator) {


            int32_t rna_end = Utils::mod(cur_pos - 10,dna_length);

            int32_t rna_length = 0;
            rna_length = loop_size + TERM_SIZE -1;

            if (rna_length >= 0) {
              int glob_rna_idx = -1;
              glob_rna_idx =
                  current_individuals[indiv_id]->metadata_->rna_count_++;

              current_individuals[indiv_id]->metadata_->rna_add(glob_rna_idx,
                                                                 prom_pos,
                                                                 rna_end,
                                                                 !lead_lag,
                                                                 1.0 -
                                                                 prom_error /
                                                                 5.0, rna_length);
            }

          }
        }

      }
    }
  }
}

void ExpManager_7::compute_RNA(int indiv_id) {

  {
    current_individuals[indiv_id]->metadata_->rnas_resize(
        current_individuals[indiv_id]->metadata_->promoter_count());
#ifdef WITH_FINETASKLOOP
#pragma omp taskloop grainsize(rna_grain_size)
#endif
    for (int rna_idx = 0; rna_idx <
                          (int)current_individuals[indiv_id]->metadata_->promoter_count(); rna_idx++) {
      {
        if (current_individuals[indiv_id]->metadata_->promoters(rna_idx) != nullptr) {
          if (current_individuals[indiv_id]->metadata_->promoters(rna_idx)->leading_or_lagging) {
            if (current_individuals[indiv_id]->metadata_->terminator_count(LEADING) != 0) {


              int k = current_individuals[indiv_id]->metadata_->promoters(rna_idx)->pos + 22;
              k = k >= dna_size[indiv_id] ? k - dna_size[indiv_id] : k;

              int32_t next_rna_end =
                  current_individuals[indiv_id]->metadata_->next_terminator(LEADING,k);

              int32_t rna_end =
                  next_rna_end + 10 >= dna_size[indiv_id] ?
                  next_rna_end + 10 - dna_size[indiv_id] :
                  next_rna_end + 10;

              int32_t rna_length = 0;

              if (current_individuals[indiv_id]->metadata_->promoters(rna_idx)->pos
                  > rna_end)
                rna_length = dna_size[indiv_id] -
                             current_individuals[indiv_id]->metadata_->promoters(rna_idx)->pos
                             + rna_end;
              else
                rna_length = rna_end - current_individuals[indiv_id]->
                    metadata_->promoters(rna_idx)->pos;

              rna_length -= 21;

              if (rna_length >= 0) {


                int glob_rna_idx = -1;
#pragma omp critical
                {
                  glob_rna_idx =
                      current_individuals[indiv_id]->metadata_->rna_count();
                  current_individuals[indiv_id]->metadata_->set_rna_count(
                      current_individuals[indiv_id]->metadata_->rna_count() + 1);
                }

                current_individuals[indiv_id]->metadata_->rna_add(glob_rna_idx, new Rna_7(current_individuals[indiv_id]->metadata_->promoters(rna_idx)->pos,
                    rna_end,
                    !current_individuals[indiv_id]->metadata_->promoters(rna_idx)->leading_or_lagging,
                    1.0 -
                    fabs(
                        ((float)current_individuals[indiv_id]->metadata_->promoters(rna_idx)->error)) /
                    5.0, rna_length));
              }
            }
          } else {
            // LAGGING
            if (current_individuals[indiv_id]->metadata_->terminator_count(LAGGING) != 0) {




              // Search for terminator
              int k = current_individuals[indiv_id]->metadata_->promoters(rna_idx)->pos - 22;
              k = k < 0 ? dna_size[indiv_id] + k : k;

              int32_t next_rna_end =
                  current_individuals[indiv_id]->metadata_->next_terminator(LAGGING,k);


              int32_t rna_end =
                  next_rna_end - 10 < 0 ? dna_size[indiv_id] + (next_rna_end - 10)
                                        :
                  next_rna_end - 10;

              int32_t rna_length = 0;

              if (current_individuals[indiv_id]->metadata_->promoters(rna_idx)->pos <
                  rna_end)
                rna_length = current_individuals[indiv_id]->metadata_->promoters(rna_idx)->pos +
                    dna_size[indiv_id] - rna_end;
              else
                rna_length = current_individuals[indiv_id]->metadata_->promoters(rna_idx)->pos -
                    rna_end;

              rna_length -= 21;

              if (rna_length >= 0) {


                int glob_rna_idx;
#pragma omp critical
                {
                  glob_rna_idx =
                      current_individuals[indiv_id]->metadata_->rna_count();
                  current_individuals[indiv_id]->metadata_->set_rna_count(
                      current_individuals[indiv_id]->metadata_->rna_count() + 1);
                }

                current_individuals[indiv_id]->metadata_->rna_add(rna_idx, new Rna_7(current_individuals[indiv_id]->metadata_->promoters(rna_idx)->pos,
                    rna_end,
                    !current_individuals[indiv_id]->metadata_->promoters(rna_idx)->leading_or_lagging,
                    1.0 -
                    fabs(
                        ((float)current_individuals[indiv_id]->metadata_->promoters(rna_idx)->error)) /
                    5.0, rna_length));

              }
            }
          }
        }
      }
    }
  }

}

void ExpManager_7::start_protein(int indiv_id) {
  current_individuals[indiv_id]->metadata_->rna_begin();

  for (int rna_idx = 0; rna_idx <
                        (int)current_individuals[indiv_id]->metadata_->rna_count(); rna_idx++) {
    {

      Rna_7* rna = current_individuals[indiv_id]->metadata_->rna_next();
      int32_t dna_length = dna_size[indiv_id];

      if (rna->is_init_) {
        int c_pos = rna->begin;

        if (rna->length >= 21) {
          if (rna->leading_lagging ==
              0) {
            c_pos += 22;
            c_pos =
                c_pos >= dna_length ? c_pos - dna_length
                                    : c_pos;
          } else {
            c_pos -= 22;
            c_pos =
                c_pos < 0 ? ((int) dna_length) + c_pos : c_pos;
          }


          int loop_size = 0;
          while (loop_size+DO_TRANSLATION_LOOP < rna->length) {
            bool start = false;
            int k_t;

            if (rna->leading_lagging ==
                0) {
              // Search for Shine Dalgarro + START codon on LEADING
              for (int k = 0; k < 9; k++) {
                k_t = k >= 6 ? k + 4 : k;

                if (current_individuals[indiv_id]->dna_->data_[Utils::mod((c_pos + k_t),dna_length)] ==
                    SHINE_DAL_SEQ_LEAD[k]) {
                  start = true;
                } else {
                  start = false;
                  break;
                }

              }
            } else {
              // Search for Shine Dalgarro + START codon on LAGGING
              for (int k = 0; k < 9; k++) {
                k_t = k >= 6 ? k + 4 : k;

                if (current_individuals[indiv_id]->dna_->data_[Utils::mod((c_pos - k_t),dna_length)] ==
                    SHINE_DAL_SEQ_LAG[k]) {
                  start = true;
                } else {
                  start = false;
                  break;
                }
              }

            }

            if (start) {
              rna->start_prot_count_++;

              rna->start_prot.push_back(c_pos);
            }

            if (rna->leading_lagging == 0) {
              c_pos++;
              c_pos =
                  c_pos >= dna_length ? c_pos - dna_length
                                      : c_pos;
            } else {
              c_pos--;
              c_pos = c_pos < 0 ? dna_length + c_pos : c_pos;
            }
            loop_size++;
          }
        }
      }
    }
  }
}


void ExpManager_7::compute_protein(int indiv_id) {
  int resize_to = 0;

  current_individuals[indiv_id]->metadata_->rna_begin();
  for (int rna_idx = 0; rna_idx <
                        (int)current_individuals[indiv_id]->metadata_->rna_count(); rna_idx++) {
    Rna_7* rna = current_individuals[indiv_id]->metadata_->rna_next();

    if (rna->is_init_)
      resize_to += rna->start_prot_count_;
  }

  current_individuals[indiv_id]->
      metadata_->proteins_resize(resize_to);

  Dna_7* dna = current_individuals[indiv_id]->dna_;
  int32_t dna_length = dna->length();

  current_individuals[indiv_id]->metadata_->rna_begin();
  for (int rna_idx = 0; rna_idx <
                        (int)current_individuals[indiv_id]->metadata_->rna_count(); rna_idx++) {

    Rna_7* rna = current_individuals[indiv_id]->metadata_->rna_next();

    if (rna->is_init_) {
      for (auto it_start_pos = rna->start_prot.begin(); it_start_pos != rna->start_prot.end(); it_start_pos++) {

        int start_prot = *it_start_pos;
        int start_protein_pos = rna->leading_lagging == 0 ?
                                start_prot +
                                13 :
                                start_prot -
                                13;
        int length;

        if (rna->leading_lagging == 0) {
          start_protein_pos = start_protein_pos >= dna_length ?
                              start_protein_pos - dna_length
                                                              : start_protein_pos;

          if (start_prot < rna->end) {
            length = rna->end - start_prot;
          } else {
            length = dna_length -
                     start_prot +
                     rna->end;

          }

          length -= 13;
        } else {


          start_protein_pos = start_protein_pos < 0 ?
                              dna_length + start_protein_pos
                                                    : start_protein_pos;

          if (start_prot > rna->end) {
            length = (*it_start_pos) - rna->end;
          } else {
            length = *it_start_pos +
                     dna_length -
                     rna->end;
          }

          length -= 13;
        }

        bool is_protein = false;

        length += 1;
        length = length - (length % 3);

        int j = 0;
        int32_t transcribed_start = 0;

        if (rna->leading_lagging == 0) {
          transcribed_start = rna->begin + 22;
          transcribed_start = transcribed_start >= dna_length ?
                              transcribed_start - dna_length
                                                              : transcribed_start;

          if (transcribed_start <= start_prot) {
            j = start_prot - transcribed_start;
          } else {
            j = dna_length -
                transcribed_start +
                start_prot;

          }
        } else {
          transcribed_start = rna->begin - 22;
          transcribed_start = transcribed_start < 0 ?
                              dna_length + transcribed_start
                                                    : transcribed_start;

          if (transcribed_start >=
              start_prot) {
            j = transcribed_start -
                start_prot;
          } else {
            j = transcribed_start +
                dna_length - start_prot;
          }
        }

        j += 13;

        while (rna->length - j >= 3) {
          int t_k;

          if (rna->leading_lagging == 0) {
            start_protein_pos =
                start_protein_pos >= dna_length ?
                start_protein_pos - dna_length
                                                : start_protein_pos;

            is_protein = false;

            for (int k = 0; k < 3; k++) {
              t_k = start_protein_pos + k;

              if (dna->get_lead(t_k) ==
                  PROTEIN_END_LEAD[k]) {
                is_protein = true;
              } else {
                is_protein = false;
                break;
              }
            }

            if (is_protein) {
              int prot_length = -1;

              if (start_prot + 13 < t_k) {
                prot_length = t_k -
                              (start_prot +
                               13);
              } else {
                prot_length = dna_length -
                              (start_prot +
                               13) + t_k;
              }

              if (prot_length >= 3) {
                int32_t glob_prot_idx = -1;
                glob_prot_idx =
                    current_individuals[indiv_id]->metadata_->proteins_count();
                current_individuals[indiv_id]->metadata_->set_proteins_count(
                    current_individuals[indiv_id]->metadata_->proteins_count() +
                    1);

                current_individuals[indiv_id]->
                    metadata_->protein_add(glob_prot_idx, new Protein_7(
                    Utils::mod(start_prot+13,dna_length), Utils::mod(t_k,dna_length),
                    prot_length/3,
                    rna->leading_lagging,
                    rna->e,rna
                ));

                rna->is_coding_ = true;
              }

              break;
            }

            start_protein_pos += 3;
            start_protein_pos =
                start_protein_pos >= dna_length ?
                start_protein_pos - dna_length
                                                : start_protein_pos;
          } else {


            start_protein_pos = start_protein_pos < 0 ?
                                dna_length + start_protein_pos
                                                      : start_protein_pos;

            is_protein = false;
            for (int k = 0; k < 3; k++) {
              t_k = start_protein_pos - k;

              if (dna->get_lag(t_k) ==
                  PROTEIN_END_LAG[k]) {
                is_protein = true;
              } else {
                is_protein = false;
                break;
              }
            }

            if (is_protein) {
              int prot_length = -1;
              if (start_prot - 13 > t_k) {
                prot_length =
                    (start_prot - 13) -
                    t_k;
              } else {
                prot_length =
                    (start_prot - 13) +
                    dna_length - t_k;
              }
              if (prot_length >= 3) {
                int32_t glob_prot_idx = -1;

                glob_prot_idx =
                    current_individuals[indiv_id]->metadata_->proteins_count();
                current_individuals[indiv_id]->metadata_->set_proteins_count(
                    current_individuals[indiv_id]->metadata_->proteins_count() +
                    1);
                ((List_Metadata*)current_individuals[indiv_id]->metadata_)->protein_add(
                    glob_prot_idx,
                        new Protein_7(Utils::mod(start_prot-13,dna_length), Utils::mod(t_k,dna_length),
                        prot_length/3,
                        rna->leading_lagging,
                        rna->e,rna
                    ));
                rna->is_coding_ = true;
              }
              break;
            }
            start_protein_pos = start_protein_pos - 3;
            start_protein_pos = start_protein_pos < 0 ?
                                dna_length + start_protein_pos
                                                      : start_protein_pos;
          }
          j += 3;
        }
      }
    }
  }

#ifdef __REGUL
  ((List_Metadata*)current_individuals[indiv_id]->metadata_)->add_inherited_proteins();
#endif
}


void ExpManager_7::translate_protein(int indiv_id, double w_max) {
  current_individuals[indiv_id]->metadata_->protein_begin();
  for (int protein_idx = 0; protein_idx <
                            (int)current_individuals[indiv_id]->metadata_->proteins_count(); protein_idx++) {
    Protein_7* prot = current_individuals[indiv_id]->metadata_->protein_next();

    if (prot->is_init_) {
      int c_pos = prot->protein_start, t_pos;
      int end_pos = prot->protein_end;
      if (prot->leading_lagging ==
          0) {
        end_pos -= 3;

        c_pos =
            c_pos >= dna_size[indiv_id] ? c_pos - dna_size[indiv_id]
                                        : c_pos;
        end_pos = end_pos < 0 ? dna_size[indiv_id] + end_pos : end_pos;
      } else {
        end_pos += 3;

        end_pos =
            end_pos >= dna_size[indiv_id] ? end_pos - dna_size[indiv_id]
                                          : end_pos;
        c_pos = c_pos < 0 ? dna_size[indiv_id] + c_pos : c_pos;
      }

      int8_t value = 0;
      int8_t codon_idx = 0;
      int32_t count_loop = 0;

      if (prot->leading_lagging ==
          0) {
        // LEADING

        while (count_loop <
               prot->protein_length &&
               codon_idx < 64*3) {
          value = 0;
          for (int8_t i = 0; i < 3; i++) {
            t_pos = c_pos + i;
            if (current_individuals[indiv_id]->dna_->get_lead(t_pos) ==
                '1')
              value += 1 << (CODON_SIZE - i - 1);
          }

          // if (indiv_id == 660) printf("Protein %d :: Add codon %d : %d\n",prot->protein_start,codon_idx,value);
          prot->codon_list[codon_idx] = value;

          codon_idx++;

          count_loop++;
          c_pos += 3;
          c_pos =
              c_pos >= dna_size[indiv_id] ? c_pos - dna_size[indiv_id]
                                          : c_pos;
        }
      } else {
        // LAGGING
        while (count_loop <
               prot->protein_length &&
               codon_idx < 64*3) {
          value = 0;
          for (int8_t i = 0; i < 3; i++) {
            t_pos = c_pos - i;
            if (current_individuals[indiv_id]->dna_->get_lag(t_pos) !=
                '1')
              value += 1 << (CODON_SIZE - i - 1);
          }

          // if (indiv_id == 660) printf("Protein %d :: Add codon %d : %d\n",prot->protein_start,codon_idx,value);
          prot->codon_list[codon_idx] = value;
          codon_idx++;

          count_loop++;

          c_pos -= 3;
          c_pos = c_pos < 0 ? c_pos + dna_size[indiv_id] : c_pos;
        }
      }

      if (count_loop < prot->protein_length && codon_idx==64*3) {
        std::ofstream last_gener_file("aevol_run.log",
                                    std::ofstream::out);
        last_gener_file << "Stop translating protein before end (length is greater than 196" << std::endl;
        last_gener_file.close();
      }

      double M = 0.0;
      double W = 0.0;
      double H = 0.0;

      int32_t nb_m = 0;
      int32_t nb_w = 0;
      int32_t nb_h = 0;

      bool bin_m = false; // Initializing to false will yield a conservation of the high weight bit
      bool bin_w = false; // when applying the XOR operator for the Gray to standard conversion
      bool bin_h = false;

      prot->nb_codons_ = codon_idx-1;

      for (int i = 0; i < codon_idx; i++) {
        switch (prot->codon_list[i]) {
        case CODON_M0 : {
          // M codon found
          nb_m++;

          // Convert Gray code to "standard" binary code
          bin_m ^= false; // as bin_m was initialized to false, the XOR will have no effect on the high weight bit

          // A lower-than-the-previous-lowest weight bit was found, make a left bitwise shift
          //~ M <<= 1;
          M *= 2;

          // Add this nucleotide's contribution to M
          if (bin_m) M += 1;

          break;
        }
        case CODON_M1 : {
          // M codon found
          nb_m++;

          // Convert Gray code to "standard" binary code
          bin_m ^= true; // as bin_m was initialized to false, the XOR will have no effect on the high weight bit

          // A lower-than-the-previous-lowest bit was found, make a left bitwise shift
          //~ M <<= 1;
          M *= 2;

          // Add this nucleotide's contribution to M
          if (bin_m) M += 1;

          break;
        }
        case CODON_W0 : {
          // W codon found
          nb_w++;

          // Convert Gray code to "standard" binary code
          bin_w ^= false; // as bin_m was initialized to false, the XOR will have no effect on the high weight bit

          // A lower-than-the-previous-lowest weight bit was found, make a left bitwise shift
          //~ W <<= 1;
          W *= 2;

          // Add this nucleotide's contribution to W
          if (bin_w) W += 1;

          break;
        }
        case CODON_W1 : {
          // W codon found
          nb_w++;

          // Convert Gray code to "standard" binary code
          bin_w ^= true; // as bin_m was initialized to false, the XOR will have no effect on the high weight bit

          // A lower-than-the-previous-lowest weight bit was found, make a left bitwise shift
          //~ W <<= 1;
          W *= 2;

          // Add this nucleotide's contribution to W
          if (bin_w) W += 1;

          break;
        }
        case CODON_H0 :
        case CODON_START : // Start codon codes for the same amino-acid as H0 codon
        {
          // H codon found
          nb_h++;

          // Convert Gray code to "standard" binary code
          bin_h ^= false; // as bin_m was initialized to false, the XOR will have no effect on the high weight bit

          // A lower-than-the-previous-lowest weight bit was found, make a left bitwise shift
          //~ H <<= 1;
          H *= 2;

          // Add this nucleotide's contribution to H
          if (bin_h) H += 1;

          break;
        }
        case CODON_H1 : {
          // H codon found
          nb_h++;

          // Convert Gray code to "standard" binary code
          bin_h ^= true; // as bin_m was initialized to false, the XOR will have no effect on the high weight bit

          // A lower-than-the-previous-lowest weight bit was found, make a left bitwise shift
          //~ H <<= 1;
          H *= 2;

          // Add this nucleotide's contribution to H
          if (bin_h) H += 1;

          break;
        }
        }
      }



      //  ----------------------------------------------------------------------------------
      //  2) Normalize M, W and H values in [0;1] according to number of codons of each kind
      //  ----------------------------------------------------------------------------------
      prot->m =
          nb_m != 0 ? M / (pow(2, nb_m) - 1) : 0.5;
      prot->w =
          nb_w != 0 ? W / (pow(2, nb_w) - 1) : 0.0;
      prot->h =
          nb_h != 0 ? H / (pow(2, nb_h) - 1) : 0.5;

      //  ------------------------------------------------------------------------------------
      //  3) Normalize M, W and H values according to the allowed ranges (defined in macros.h)
      //  ------------------------------------------------------------------------------------
      // x_min <= M <= x_max
      // w_min <= W <= w_max
      // h_min <= H <= h_max
      prot->m =
          (X_MAX - X_MIN) *
          prot->m +
          X_MIN;
      prot->w =
          (w_max - W_MIN) *
          prot->w +
          W_MIN;
      prot->h =
          (H_MAX - H_MIN) *
          prot->h +
          H_MIN;

      if (nb_m == 0 || nb_w == 0 || nb_h == 0 ||
          prot->w ==
          0.0 ||
          prot->h ==
          0.0) {
        prot->is_functional = false;
      } else {
        prot->is_functional = true;
      }
    }
  }


  std::map<int32_t, Protein_7*> lookup;

  current_individuals[indiv_id]->metadata_->protein_begin();

  for (int protein_idx = 0; protein_idx <
                            (int)current_individuals[indiv_id]->metadata_->proteins_count(); protein_idx++) {
    {

      Protein_7* prot  =
          current_individuals[indiv_id]->metadata_->protein_next();
      if (prot->is_init_ && prot->leading_lagging==0) {
        if (lookup.find(prot->protein_start) ==
            lookup.end()) {
          lookup[prot->protein_start] = prot;
        } else {
          lookup[prot->protein_start]->e += prot->e;
          lookup[prot->protein_start]->initial_e_ += prot->initial_e_;
          lookup[prot->protein_start]->rna_list_.insert(
              lookup[prot->protein_start]->rna_list_.end(),
              prot->rna_list_.begin(),prot->rna_list_.end());
          prot->is_init_ = false;
        }
      }
    }
  }

  lookup.clear();

  current_individuals[indiv_id]->metadata_->protein_begin();

  for (int protein_idx = 0; protein_idx <
                            (int)current_individuals[indiv_id]->metadata_->proteins_count(); protein_idx++) {
    {

      Protein_7* prot  =
          current_individuals[indiv_id]->metadata_->protein_next();
      if (prot->is_init_ && prot->leading_lagging==1) {
        if (lookup.find(prot->protein_start) ==
            lookup.end()) {
          lookup[prot->protein_start] = prot;
        } else {
          lookup[prot->protein_start]->e += prot->e;
          lookup[prot->protein_start]->initial_e_ += prot->initial_e_;
          lookup[prot->protein_start]->rna_list_.insert(
              lookup[prot->protein_start]->rna_list_.end(),
              prot->rna_list_.begin(),prot->rna_list_.end());
          prot->is_init_ = false;
        }
      }
    }
  }
}

void ExpManager_7::compute_phenotype(int indiv_id) {
#ifdef PHENOTYPE_VECTOR
  double activ_phenotype[PHENOTYPE_VECTOR_SIZE];
        double inhib_phenotype[PHENOTYPE_VECTOR_SIZE];
        {
            for (int fuzzy_idx = 0; fuzzy_idx < PHENOTYPE_VECTOR_SIZE; fuzzy_idx++) {
                activ_phenotype[fuzzy_idx] = 0;
                inhib_phenotype[fuzzy_idx] = 0;
            }
        }
    Fuzzy* f_activ_phenotype = new Fuzzy();
    Fuzzy* f_inhib_phenotype = new Fuzzy();
#else
  Vector_Fuzzy* activ_phenotype = new Vector_Fuzzy();
  Vector_Fuzzy* inhib_phenotype = new Vector_Fuzzy();
#endif

  std::vector<Protein_7*> protein_vector;
  current_individuals[indiv_id]->metadata_->protein_begin();
  for (int protein_idx = 0; protein_idx < current_individuals[indiv_id]->metadata_->proteins_count(); protein_idx++) {
    Protein_7* prot = current_individuals[indiv_id]->metadata_->protein_next();
    if (prot->is_init_) {
      if (prot->leading_lagging==0)
        protein_vector.push_back(prot);
    }
  }

  current_individuals[indiv_id]->metadata_->protein_begin();
  for (int protein_idx = 0; protein_idx < current_individuals[indiv_id]->metadata_->proteins_count(); protein_idx++) {
    Protein_7* prot = current_individuals[indiv_id]->metadata_->protein_next();
    if (prot->is_init_) {
      if (prot->leading_lagging==1)
        protein_vector.push_back(prot);
    }
  }

  std::sort(protein_vector.begin(), protein_vector.end(),
            [](Protein_7*a, Protein_7*b) { return *a < *b;});
  for (auto prot : protein_vector) {
    if (fabs(
        prot->w) >=
        1e-15 &&
        fabs(
            prot->h) >=
        1e-15) {

      if (prot->is_functional) {
        // Compute triangle points' coordinates
#ifdef PHENOTYPE_VECTOR
        double x0 =
                            prot->m -
                            prot->w;
                    double x1 = prot->m;
                    double x2 =
                            prot->m +
                            prot->w;

                    double height = (prot->h *
                                     prot->e);

                    int loop_A_start = (int) std::ceil(x0 * D_PHENOTYPE_VECTOR_SIZE);
                    loop_A_start = loop_A_start < 0 ? 0 : loop_A_start;
                    loop_A_start = loop_A_start > PHENOTYPE_VECTOR_SIZE ? PHENOTYPE_VECTOR_SIZE : loop_A_start;

                    int loop_A_end = (int) std::ceil(x1 * D_PHENOTYPE_VECTOR_SIZE);
                    loop_A_end = loop_A_end < 0 ? 0 : loop_A_end;
                    loop_A_end = loop_A_end > PHENOTYPE_VECTOR_SIZE ? PHENOTYPE_VECTOR_SIZE : loop_A_end;

                    for (int i = loop_A_start; i < loop_A_end; i++) {
                        if (prot->h > 0)
                            activ_phenotype[i] += (((i / D_PHENOTYPE_VECTOR_SIZE) - x0) / (x1 - x0)) * height;
                        else
                            inhib_phenotype[i] += (((i / D_PHENOTYPE_VECTOR_SIZE) - x0) / (x1 - x0)) * height;
                    }

                    // Compute the second equation of the triangle
                    // Updating value between x1 and x2
                    int loop_B_start = (int) std::ceil(x1 * D_PHENOTYPE_VECTOR_SIZE);
                    loop_B_start = loop_B_start < 0 ? 0 : loop_B_start;
                    loop_B_start = loop_B_start > PHENOTYPE_VECTOR_SIZE ? PHENOTYPE_VECTOR_SIZE : loop_B_start;

                    int loop_B_end = (int) std::ceil(x2 * D_PHENOTYPE_VECTOR_SIZE);
                    if (loop_B_end > PHENOTYPE_VECTOR_SIZE) {
                        if (prot->h > 0)
                            activ_phenotype[PHENOTYPE_VECTOR_SIZE-1] += height * ((x2 - 1.0) / (x2 - x1));
                        else
                            inhib_phenotype[PHENOTYPE_VECTOR_SIZE-1] += height * ((x2 - 1.0) / (x2 - x1));
                    }

                    loop_B_end = loop_B_end < 0 ? 0 : loop_B_end;
                    loop_B_end = loop_B_end > PHENOTYPE_VECTOR_SIZE ? PHENOTYPE_VECTOR_SIZE : loop_B_end;

                    for (int i = loop_B_start; i < loop_B_end; i++) {
                        if (prot->h > 0)
                            activ_phenotype[i] += height * ((x2 - (i / D_PHENOTYPE_VECTOR_SIZE)) / (x2 - x1));
                        else
                            inhib_phenotype[i] += height * ((x2 - (i / D_PHENOTYPE_VECTOR_SIZE)) / (x2 - x1));
                    }

                    if (prot->h > 0)
                        f_activ_phenotype->add_triangle(prot->m, prot->w, prot->h *
                                                                        prot->e);
                    else
                        f_inhib_phenotype->add_triangle(prot->m, prot->w, prot->h *
                                                                        prot->e);

//                    printf("SIMD_FUZZY ACTIV\n");
//                    for (int i = 0; i < PHENOTYPE_VECTOR_SIZE; i++)
//                        if (activ_phenotype[i]!=0) printf("[%lf : %lf]\n",((double)i)/D_PHENOTYPE_VECTOR_SIZE,activ_phenotype[i]);
//                    printf("\n");
//
//                    printf("CLASSIC FUZZY ACTIV\n");
//                    f_activ_phenotype->print();
//
//                    printf("SIMD FUZZY INHIB\n");
//                    for (int i = 0; i < PHENOTYPE_VECTOR_SIZE; i++)
//                        if (inhib_phenotype[i]!=0) printf("[%d : %f]\n",((double)i)/D_PHENOTYPE_VECTOR_SIZE,inhib_phenotype[i]);
//                    printf("\n");
//
//                    printf("CLASSIC FUZZY INHIB\n");
//                    f_inhib_phenotype->print();
         double geom_active = 0, geom_inhib = 0;
//
        for (int fuzzy_idx = 0; fuzzy_idx < PHENOTYPE_VECTOR_SIZE-1; fuzzy_idx++) {
            geom_active +=
                    ((std::fabs(activ_phenotype[fuzzy_idx]) +
                      std::fabs(activ_phenotype[fuzzy_idx + 1])) /
                     (D_PHENOTYPE_VECTOR_SIZE*2));
            geom_inhib += ((std::fabs(inhib_phenotype[fuzzy_idx]) +
                      std::fabs(inhib_phenotype[fuzzy_idx + 1])) /
                     (D_PHENOTYPE_VECTOR_SIZE*2));
        }

        double fgeom_active = f_activ_phenotype->get_geometric_area();
        double fgeom_active_round = roundf(fgeom_active * 10000);
        double geom_active_round = roundf(geom_active * 10000);

        double fgeom_inhib = f_inhib_phenotype->get_geometric_area();
        double fgeom_inhib_round = roundf(fgeom_inhib * 10000);
        double geom_inhib_round = roundf(geom_inhib * 10000);

        if ((geom_active_round != fgeom_active_round) || (geom_inhib_round != fgeom_inhib_round)) {
        printf("After adding triangle (%lf %lf %lf) : Active %.10lf/%.10lf Inhib %.10lf/%.10lf\n",prot->m,prot->w,prot->h*prot->e,
        geom_active,fgeom_active,
        geom_inhib,fgeom_inhib);
        exit(5);
        }
#else
        if (prot->h > 0)
          activ_phenotype->add_triangle(prot->m, prot->w, prot->h *
                                                          prot->e,false);
        else
          inhib_phenotype->add_triangle(prot->m, prot->w, prot->h *
                                                          prot->e,false);
#endif
      }
    }
  }


#ifdef PHENOTYPE_VECTOR
  for (int fuzzy_idx = 0; fuzzy_idx < PHENOTYPE_VECTOR_SIZE; fuzzy_idx++) {
            if (activ_phenotype[fuzzy_idx] > 1)
                activ_phenotype[fuzzy_idx] = 1;
            if (inhib_phenotype[fuzzy_idx] < -1)
                inhib_phenotype[fuzzy_idx] = -1;
        }

        for (int fuzzy_idx = 0; fuzzy_idx < PHENOTYPE_VECTOR_SIZE; fuzzy_idx++) {
            current_individuals[indiv_id]->phenotype[fuzzy_idx] = activ_phenotype[fuzzy_idx] + inhib_phenotype[fuzzy_idx];
            if (current_individuals[indiv_id]->phenotype[fuzzy_idx] < 0)
                current_individuals[indiv_id]->phenotype[fuzzy_idx] = 0;
        }
#else

  activ_phenotype->clip(AbstractFuzzy::max,   Y_MAX);
  inhib_phenotype->clip(AbstractFuzzy::min, - Y_MAX);

  // activ_phenotype->simplify();
  // inhib_phenotype->simplify();
  current_individuals[indiv_id]->phenotype = new Vector_Fuzzy();
  current_individuals[indiv_id]->phenotype->add(activ_phenotype);
  current_individuals[indiv_id]->phenotype->add(inhib_phenotype);
  current_individuals[indiv_id]->phenotype->clip(AbstractFuzzy::min, Y_MIN);
  current_individuals[indiv_id]->phenotype->simplify();

  delete activ_phenotype;
  delete inhib_phenotype;
#endif
}


void ExpManager_7::build_phenotypic_target(PhenotypicTargetHandler* phenotypic_target_handler) {
#ifdef PHENOTYPE_VECTOR
  for (int16_t i = 0; i <= 300; i++) {
    Point new_point = Point(
        X_MIN + (double) i * (X_MAX - X_MIN) / (double) 300, 0.0);
    //printf("Gaussians %d\n",phenotypic_target_handler->current_gaussians().size());

    for (const Gaussian& g: phenotypic_target_handler->current_gaussians())
      new_point.y += g.compute_y(new_point.x);

    //printf("Sampling %d = %lf\n",(int)(new_point.x * 300),new_point.y);

    if (i < 300)
      target[(int)(new_point.x * 300)] = (float) new_point.y;
  }

  for (int i = 1; i < 300; i++) {
      if (target[i] == 0.0) {
        int minL = i - 1;
        int maxL = i + 1;
        int dist = 1;

        while (target[maxL] == 0.0) {
          maxL++;
          dist++;
        }
        float inc = 0.0;
        if (target[maxL] > target[minL]) {
          inc = (target[maxL] - target[minL]) / dist;
        } else {
          inc = (target[minL] - target[maxL]) / dist;
          minL = maxL;
        }

        for (int j = i; j < maxL; j++) {
          target[j] = target[minL] + inc;
          inc += inc;
        }

      }
    }

  for (int fuzzy_idx = 0; fuzzy_idx < 300; fuzzy_idx++) {

    if (target[fuzzy_idx] > 1)
      target[fuzzy_idx] = 1;
    if (target[fuzzy_idx] < 0)
      target[fuzzy_idx] = 0;
  }
#endif
}

void ExpManager_7::compute_fitness(int indiv_id, double selection_pressure, int env_id) {
#ifdef __REGUL
  Vector_Fuzzy* delta = new Vector_Fuzzy(*current_individuals[indiv_id]->phenotype);
      if (phenotypic_target_handler_->var_method_ == SWITCH_IN_A_LIST)
        delta->sub(phenotypic_target_handler_->targets_fuzzy_[env_id]);
      else if (phenotypic_target_handler_->var_method_ == ONE_AFTER_ANOTHER)
        delta->sub(phenotypic_target_handler_->targets_fuzzy_by_id_[env_id]);

      if (phenotypic_target_handler_->var_method_ == SWITCH_IN_A_LIST)
        current_individuals[indiv_id]->metaerror_by_env_id_[env_id] = delta->get_geometric_area();
      else if (phenotypic_target_handler_->var_method_ == ONE_AFTER_ANOTHER)
        current_individuals[indiv_id]->metaerror_by_env_id_[env_id] += delta->get_geometric_area();

      delete delta;

      if (phenotypic_target_handler_->var_method_ == SWITCH_IN_A_LIST)
        current_individuals[indiv_id]->fitness_by_env_id_[env_id] = exp(
            -selection_pressure *
            ((double) current_individuals[indiv_id]->metaerror_by_env_id_[env_id]));
#else
#ifdef PHENOTYPE_VECTOR
  for (int fuzzy_idx = 0; fuzzy_idx < PHENOTYPE_VECTOR_SIZE; fuzzy_idx++) {

            if (current_individuals[indiv_id]->phenotype[fuzzy_idx] > 1)
                current_individuals[indiv_id]->phenotype[fuzzy_idx] = 1;
            if (current_individuals[indiv_id]->phenotype[fuzzy_idx] < 0)
                current_individuals[indiv_id]->phenotype[fuzzy_idx] = 0;

            current_individuals[indiv_id]->delta[fuzzy_idx] =
                    current_individuals[indiv_id]->phenotype[fuzzy_idx] -
                    target[fuzzy_idx];
        }

        current_individuals[indiv_id]->metaerror = 0;

        for (int fuzzy_idx = 0; fuzzy_idx < PHENOTYPE_VECTOR_SIZE; fuzzy_idx++) {
            current_individuals[indiv_id]->metaerror +=
                    ((std::fabs(current_individuals[indiv_id]->delta[fuzzy_idx]) +
                      std::fabs(current_individuals[indiv_id]->delta[fuzzy_idx + 1])) /
                     (D_PHENOTYPE_VECTOR_SIZE*2));
        }
#else
  Vector_Fuzzy* delta = new Vector_Fuzzy(*current_individuals[indiv_id]->phenotype);
  delta->sub(target);

  current_individuals[indiv_id]->metaerror = delta->get_geometric_area();
  delete delta;
#endif

        current_individuals[indiv_id]->fitness = exp(
      -selection_pressure *
      ((double)current_individuals[indiv_id]->metaerror));
#endif
}


#ifdef __REGUL
void ExpManager_7::compute_network(int indiv_id, double selection_pressure) {
  // Allocate fitness and metarror array
  if (phenotypic_target_handler_->var_method_ == SWITCH_IN_A_LIST) {
    current_individuals[indiv_id]->fitness_by_env_id_ = new double[phenotypic_target_handler_->nb_eval_];
    current_individuals[indiv_id]->metaerror_by_env_id_ = new double[phenotypic_target_handler_->nb_eval_];
  } else if (phenotypic_target_handler_->var_method_ == ONE_AFTER_ANOTHER) {
    current_individuals[indiv_id]->fitness_by_env_id_ = new double[phenotypic_target_handler_->nb_env_];
    current_individuals[indiv_id]->metaerror_by_env_id_ = new double[phenotypic_target_handler_->nb_env_];

    for (int env_id = 0; env_id < phenotypic_target_handler_->nb_env_; env_id++) {
      current_individuals[indiv_id]->fitness_by_env_id_[env_id] = 0;
      current_individuals[indiv_id]->metaerror_by_env_id_[env_id] = 0;
    }
  }

  // Add signals
  ((List_Metadata*)current_individuals[indiv_id]->metadata_)->signal_proteins_.resize(phenotypic_target_handler_->signals_models_.size());
  int i = 0;
  for (auto signal : phenotypic_target_handler_->signals_models_) {
    int glob_prot_idx = current_individuals[indiv_id]->metadata_->proteins_count();
    current_individuals[indiv_id]->metadata_->set_proteins_count(
        current_individuals[indiv_id]->metadata_->proteins_count() +
        1);

    Protein_7* prot = new Protein_7(signal);
    current_individuals[indiv_id]->metadata_->protein_add(glob_prot_idx, prot);
    ((List_Metadata*)current_individuals[indiv_id]->metadata_)->signal_proteins_[i] = (prot);
    i++;
  }

  // Set influences
  current_individuals[indiv_id]->metadata_->rna_begin();
  for (int i = 0; i < current_individuals[indiv_id]->metadata_->rna_count(); i++) {
    Rna_7* rna = current_individuals[indiv_id]->metadata_->rna_next();
    if (rna != nullptr) {
      if (rna->is_coding_) {
        rna->nb_influences_ = 0;

        int32_t enhancer_position = rna->enhancer_position(current_individuals[indiv_id]->dna_->length());
        int32_t operator_position = rna->operator_position(current_individuals[indiv_id]->dna_->length());

        double enhance,operate;
        current_individuals[indiv_id]->metadata_->protein_begin();
        for (int j = 0; j < current_individuals[indiv_id]->metadata_->proteins_count(); j++) {
          Protein_7* prot =
              current_individuals[indiv_id]->metadata_->protein_next();
          if (prot != nullptr) {
            if (prot->is_init_) {
              enhance = rna->affinity_with_protein(
                  enhancer_position, prot, current_individuals[indiv_id],
                  exp_m_);
              operate = rna->affinity_with_protein(
                  operator_position, prot, current_individuals[indiv_id],
                  exp_m_);

              if (enhance != 0.0 || operate != 0.0) {

                rna->affinity_list.push_back(
                    AffinityFactor(prot, enhance, operate));
                prot->is_TF_ = true;

                rna->nb_influences_++;

                 if (indiv_id==543 && AeTime::time() == 5895)  
                  printf("SIMD -- Affinity between RNA %d and Protein %d : %lf %lf\n",
                         rna->begin, prot->protein_start, enhance, operate);
              }
            }
          }
        }
      }
    }
  }
}

void ExpManager_7::update_network(int indiv_id, double selection_pressure) {
//
//  current_individuals[indiv_id]->metadata_->rna_begin();
//  for (int x = 0; x < current_individuals[indiv_id]->metadata_->rna_count(); x++) {
//    Rna_7* rna = current_individuals[indiv_id]->metadata_->rna_next();
//    if (rna != nullptr) {
//      if (rna->is_coding_) {
//
//        double enhancer_activity = 0;
//        double operator_activity = 0;
//
//        for (auto affinity: rna->affinity_list) {
//          enhancer_activity +=
//              affinity.enhancer_factor * affinity.protein_concentration;
//          operator_activity +=
//              affinity.operator_factor * affinity.protein_concentration;
//        }
//
//        ProteinConcentration enhancer_activity_pow_n =
//            enhancer_activity == 0
//                ? 0
//                : pow(enhancer_activity, exp_m_->exp_s()->get_hill_shape_n());
//        ProteinConcentration operator_activity_pow_n =
//            operator_activity == 0
//                ? 0
//                : pow(operator_activity, exp_m_->exp_s()->get_hill_shape_n());
//        rna->synthesis_rate =
//            rna->e *
//            (exp_m_->exp_s()->get_hill_shape() /
//             (operator_activity_pow_n + exp_m_->exp_s()->get_hill_shape())) *
//            (1 + ((1 / rna->e) - 1) * (enhancer_activity_pow_n /
//                                       (enhancer_activity_pow_n +
//                                        exp_m_->exp_s()->get_hill_shape())));
//      }
//    }
//  }

  current_individuals[indiv_id]->metadata_->protein_begin();
  for (int j = 0; j < current_individuals[indiv_id]->metadata_->proteins_count(); j++) {
    Protein_7* prot =
        current_individuals[indiv_id]->metadata_->protein_next();
    if (!prot->signal_) {
      if (prot->is_init_) {
        prot->delta_concentration_ = 0;

        for (auto rna: prot->rna_list_) {
          double synthesis_rate = rna->compute_synthesis_rate(current_individuals[indiv_id]);

          if (indiv_id==543 && AeTime::time() == 5895)  printf("SIMD -- Protein %d synthesis by RNA %d at rate %lf : DELTA BEFORE %f :: ",prot->protein_start,
                            rna->begin,
                            synthesis_rate,prot->delta_concentration_);
          prot->delta_concentration_ += synthesis_rate;

          if (indiv_id==543 && AeTime::time() == 5895)  {
            printf("DELTA AFTER %lf : %lf\n",prot->delta_concentration_,synthesis_rate);
          }
        }

        prot->delta_concentration_ -=
            exp_m_->exp_s()->get_degradation_rate() * prot->e;
        prot->delta_concentration_ *=
            1.0 / exp_m_->exp_s()->get_nb_degradation_step();
      }
    }
  }

  // Apply the changes in concentrations we have just computed


//  if (indiv_id==137)
//    ((List_Metadata*)current_individuals[indiv_id]->metadata_)->proteins_print();

  current_individuals[indiv_id]->metadata_->protein_begin();
  for (int j = 0; j < current_individuals[indiv_id]->metadata_->proteins_count(); j++) {
    Protein_7* prot =
        current_individuals[indiv_id]->metadata_->protein_next();
//    printf("SIMD -- Protein %d : %d %d\n",prot->protein_start,prot->signal_,prot->is_init_);
    if (!prot->signal_) {
      if (prot->is_init_) {
        // if (indiv_id == 70 && AeTime::time() == 1595)
        //   printf("SIMD -- Protein %d : %lf + %lf\n", prot->protein_start,
        //          prot->e, prot->delta_concentration_);
        prot->e += prot->delta_concentration_;
      }
    }
  }
}

void ExpManager_7::evaluate_network(int indiv_id, double selection_pressure, int env_id) {
  update_phenotype(indiv_id);

  //   current_individuals[indiv_id]->metadata_->protein_begin();
  // for (int j = 0; j < current_individuals[indiv_id]->metadata_->proteins_count(); j++) {
  //   Protein_7* prot =
  //       current_individuals[indiv_id]->metadata_->protein_next();

  //   if (!prot->signal_)
  //     if (prot->is_init_) {
  //       double old_e = prot->e;
  //       prot->e += prot->delta_concentration_;
  //       printf("SIMD -- Evaluate/Update %d : %lf + %lf\n", prot->protein_start,
  //                prot->e, prot->delta_concentration_);
  //     }
  // }

  compute_fitness(indiv_id,selection_pressure,env_id);

  if (phenotypic_target_handler_->var_method_ == SWITCH_IN_A_LIST)
    current_individuals[indiv_id]->metaerror += current_individuals[indiv_id]->metaerror_by_env_id_[env_id];
}

void ExpManager_7::finalize_network(int indiv_id, double selection_pressure) {
   double sum_meta = current_individuals[indiv_id]->metaerror;

  if (phenotypic_target_handler_->var_method_ == SWITCH_IN_A_LIST) {
    current_individuals[indiv_id]->metaerror =
        current_individuals[indiv_id]->metaerror /
        (double)(phenotypic_target_handler_->nb_eval_);
  }
  else if (phenotypic_target_handler_->var_method_ == ONE_AFTER_ANOTHER) {
    for (int env_id = 0; env_id < phenotypic_target_handler_->nb_eval_; env_id++) {
      current_individuals[indiv_id]->metaerror_by_env_id_[env_id] =
          current_individuals[indiv_id]->metaerror_by_env_id_[env_id] / 10.0;
      current_individuals[indiv_id]->metaerror +=
          current_individuals[indiv_id]->metaerror_by_env_id_[env_id];

      current_individuals[indiv_id]->fitness_by_env_id_[env_id] = exp(
          -selection_pressure *
          ((double) current_individuals[indiv_id]->metaerror_by_env_id_[env_id]));
    }
    current_individuals[indiv_id]->metaerror =
        current_individuals[indiv_id]->metaerror / (double) (phenotypic_target_handler_->nb_env_);
  }

  current_individuals[indiv_id]->fitness = exp(
      -selection_pressure *
      ((double) current_individuals[indiv_id]->metaerror));
}

void ExpManager_7::solve_network(int indiv_id, double selection_pressure) {
  current_individuals[indiv_id]->metadata_->protein_begin();
  for (int j = 0;
       j < current_individuals[indiv_id]->metadata_->proteins_count(); j++) {
    Protein_7* prot =
        current_individuals[indiv_id]->metadata_->protein_next();
    if (!prot->signal_) {
      if (prot->is_init_) {
        prot->e = prot->initial_e_;
      }
    }
  }

  compute_network(indiv_id, selection_pressure);

  current_individuals[indiv_id]->metaerror = 0;

  if (indiv_id==190 && AeTime::time() == 1936) 
   ((List_Metadata*)current_individuals[indiv_id]->metadata_)->proteins_print(0);

  if (phenotypic_target_handler_->var_method_ == ONE_AFTER_ANOTHER) {
    for (int16_t env_i = 0; env_i < phenotypic_target_handler_->nb_env_; env_i++) {

      //Set the concentration of signals for this age
      for (auto prot1: ((List_Metadata*)current_individuals[indiv_id]->metadata_)->signal_proteins_) {
        prot1->e = 0;
      }

      for (auto prot_id : phenotypic_target_handler_->env_signals_list_[env_i]) {
        ((List_Metadata*)current_individuals[indiv_id]->metadata_)->signal_proteins_[prot_id]->e = 0.9;
      }

      for (int16_t i = 0; i < 10; i++) {

        for (int j = 0; j < 10; j++) {
          update_network(indiv_id,selection_pressure);
        }

        // If we have to evaluate the individual at this age
        evaluate_network(indiv_id,selection_pressure,env_i);
                if (indiv_id==190 && AeTime::time() == 1936)  printf("%d -- Evaluate Network at %d :: %lf %lf -- %lf\n",indiv_id,i+1,
                         current_individuals[indiv_id]->metaerror,
               current_individuals[indiv_id]->metaerror_by_env_id_[0],
                         phenotypic_target_handler_->targets_fuzzy_by_id_[0]->get_geometric_area());
      }
    }

    finalize_network(indiv_id,selection_pressure);
  } else {
    std::set<int>* eval = exp_m_->exp_s()->get_list_eval_step();
    // i is thus the age of the individual
    for (int16_t i = 0; i < phenotypic_target_handler_->nb_indiv_age_; i++) {
      //Set the concentration of signals for this age
      for (auto prot1: ((List_Metadata*)current_individuals[indiv_id]->metadata_)->signal_proteins_) {
        prot1->e = 0;
      }

      for (auto prot_id : phenotypic_target_handler_->env_signals_list_[phenotypic_target_handler_->list_env_id_[i]]) {
        ((List_Metadata*)current_individuals[indiv_id]->metadata_)->signal_proteins_[prot_id]->e = 0.9;
      }

      for (int j = 0; j < exp_m_->exp_s()->get_nb_degradation_step(); j++) {
        update_network(indiv_id,selection_pressure);
      }


  if (indiv_id==543 && AeTime::time() == 5895) 
   ((List_Metadata*)current_individuals[indiv_id]->metadata_)->proteins_print(i+1);

      // If we have to evaluate the individual at this age
      if (eval->find(i+1) != eval->end() || (indiv_id==543 && AeTime::time() == 5895)) {// ||( (indiv_id == 70) && (AeTime::time()>=1570))) {
        evaluate_network(indiv_id,selection_pressure, phenotypic_target_handler_->list_env_id_[i]);
          if (indiv_id==543 && AeTime::time() == 5895) 
            printf("%d -- Evaluate Network at %d :: %lf %lf -- %lf\n",indiv_id,i+1,
                         current_individuals[indiv_id]->metaerror,
               current_individuals[indiv_id]->metaerror_by_env_id_[0],
                         phenotypic_target_handler_->targets_fuzzy_by_id_[0]->get_geometric_area());
      }
    }

    finalize_network(indiv_id,selection_pressure);
  }
}

void ExpManager_7::update_phenotype( int indiv_id ) {
  delete current_individuals[indiv_id]->phenotype;
  compute_phenotype(indiv_id);
}
#endif


void ExpManager_7::run_a_step(double w_max, double selection_pressure,bool optim_prom) {
#pragma omp single
  {
    nb_clones_ = 0;
    cumulate_size = 0;
    cumulate_diff = 0;


    if (exp_m_->sel()->fitness_func() == FITNESS_GLOBAL_SUM) {
#ifdef __REGUL
      fitness_sum_tab_ = new double[phenotypic_target_handler_->nb_env_];
        for (int env_id = 0; env_id < phenotypic_target_handler_->nb_env_; env_id++) {
          fitness_sum_tab_[env_id] = 0;
          for (int indiv_id = 0; indiv_id < exp_m_->nb_indivs(); indiv_id++)
            fitness_sum_tab_[env_id] += previous_individuals[indiv_id]->fitness_by_env_id_[env_id];
        }
#else
      printf("Fitness local sum is not supported for Aevol (only R-Aevol)\n");
      exit(-1);
#endif
    }
  }

#pragma omp for schedule(dynamic)
  for (int g_indiv_id = 0; g_indiv_id < exp_m_->nb_indivs(); g_indiv_id += 1) {
    {
      for (int indiv_id = g_indiv_id; indiv_id < g_indiv_id + 1; indiv_id++) {
        if (ExpManager_7::standalone() && optim_prom && !exp_m_->check_simd()) {
          selection(indiv_id);
        } else if (!ExpManager_7::standalone() && optim_prom) {
        } else if (optim_prom && standalone() && exp_m_->check_simd()) {
        }



        if (optim_prom) {
          do_mutation(indiv_id);
          opt_prom_compute_RNA(indiv_id);
        } else {
          int x = indiv_id / exp_m_->world()->height();
          int y = indiv_id % exp_m_->world()->height();
          delete exp_m_->dna_mutator_array_[indiv_id];

          exp_m_->dna_mutator_array_[indiv_id] = new DnaMutator(
              exp_m_->world()->grid(x, y)->mut_prng(),
              previous_individuals[exp_m_->next_generation_reproducer_[indiv_id]]->dna_->length(),
              exp_m_->exp_s()->mut_params()->duplication_rate(),
              exp_m_->exp_s()->mut_params()->deletion_rate(),
              exp_m_->exp_s()->mut_params()->translocation_rate(),
              exp_m_->exp_s()->mut_params()->inversion_rate(),
              exp_m_->exp_s()->mut_params()->point_mutation_rate(),
              exp_m_->exp_s()->mut_params()->small_insertion_rate(),
              exp_m_->exp_s()->mut_params()->small_deletion_rate(),
              exp_m_->exp_s()->mut_params()->max_indel_size(),
              exp_m_->exp_s()->min_genome_length(),
              exp_m_->exp_s()->max_genome_length(), indiv_id, x, y);
          exp_m_->dna_mutator_array_[indiv_id]->setMutate(true);

          start_stop_RNA(indiv_id);
          compute_RNA(indiv_id);
        }

        if (exp_m_->dna_mutator_array_[indiv_id]->hasMutate()) {
          start_protein(indiv_id);
          compute_protein(indiv_id);
          translate_protein(indiv_id, w_max);
#ifdef __REGUL
          current_individuals[indiv_id]->metadata_->protein_begin();
                      for (int j = 0; j < current_individuals[indiv_id]
                                              ->metadata_->proteins_count();
                           j++) {
                        Protein_7* prot = current_individuals[indiv_id]
                                             ->metadata_->protein_next();

//                        if (indiv_id == 22)
//                          printf("INACT -- Protein %d :: %lf\n", j,
//                                 prot->e);
                      }


                        solve_network(indiv_id,selection_pressure);
#else
          compute_phenotype(indiv_id);
          compute_fitness(indiv_id, selection_pressure);
#endif
        }

        if (ExpManager_7::standalone() && optim_prom && exp_m_->record_tree()) {
          int x = indiv_id / exp_m_->world()->height();
          int y = indiv_id % exp_m_->world()->height();

          EndReplicationEvent *eindiv = new EndReplicationEvent(current_individuals[indiv_id], x, y);
          // Tell observers the replication is finished
          exp_m_->tree()->update_end_replication(eindiv);
          delete eindiv;
        }
      }
    }
  }


#pragma omp single
  {
    if (optim_prom && exp_m_->record_tree())
      exp_m_->tree()->update_end_generation();
  }

  if (optim_prom) {
#pragma omp for schedule(static)
    for (int indiv_id = 0; indiv_id < (int) exp_m_->nb_indivs(); indiv_id++) {
      bool toDelete = false;

#pragma omp critical(indiv_list)
      {
        if (previous_individuals[indiv_id]->usage_count_ == 1)
          toDelete = true;
        else
          previous_individuals[indiv_id]->usage_count_--;
      }
      if (toDelete) {
        delete previous_individuals[indiv_id];
      }

      previous_individuals[indiv_id] = current_individuals[indiv_id];
      current_individuals[indiv_id] = nullptr;
    }
  }

#pragma omp single
  {
    // Search for the best
    double best_fitness = previous_individuals[0]->fitness;
    int idx_best = 0;
    for (int indiv_id = 1; indiv_id < (int) exp_m_->nb_indivs(); indiv_id++) {
      if (previous_individuals[indiv_id]->fitness > best_fitness) {
        idx_best = indiv_id;
        best_fitness = previous_individuals[indiv_id]->fitness;
      }
    }

    best_indiv = previous_individuals[idx_best];


    // Traces
#ifdef WITH_PERF_TRACES
    std::ofstream perf_traces_file_;
            perf_traces_file_.open("simd_perf_traces.csv", std::ofstream::app);
            for (int indiv_id = 0; indiv_id < (int) exp_m_->nb_indivs(); indiv_id++) {
                perf_traces_file_ << AeTime::time() << "," << indiv_id << "," << apply_mutation[indiv_id] << std::endl;
            }
            perf_traces_file_.close();
#endif
  }

  bool without_stats = true;
  if (!without_stats) {

#pragma omp for schedule(static)
    for (int indiv_id = 0; indiv_id < (int) exp_m_->nb_indivs(); indiv_id++) {
      previous_individuals[indiv_id]->reset_stats();
      previous_individuals[indiv_id]->metadata_->rna_begin();
      for (int i = 0; i < previous_individuals[indiv_id]->metadata_->rna_count(); i++) {
        Rna_7*rna = previous_individuals[indiv_id]->metadata_->rna_next();
        if (rna != nullptr) {
          if (rna->is_coding_)
            previous_individuals[indiv_id]->nb_coding_RNAs++;
          else
            previous_individuals[indiv_id]->nb_non_coding_RNAs++;
        }
      }


      for (int i = 0; i < previous_individuals[indiv_id]->metadata_->proteins_count(); i++) {
        Protein_7*prot =
            previous_individuals[indiv_id]->metadata_->proteins(i);
        if (prot != nullptr) {
          if (prot->is_functional) {
            previous_individuals[indiv_id]->nb_func_genes++;
          } else {
            previous_individuals[indiv_id]->nb_non_func_genes++;
          }
          if (prot->h > 0) {
            previous_individuals[indiv_id]->nb_genes_activ++;
          } else {
            previous_individuals[indiv_id]->nb_genes_inhib++;
          }
        }
      }
    }


#pragma omp single
    {
      // Stats
      if (!optim_prom) {
        stats_best = new Stats_7(this, AeTime::time(), true);
        stats_mean = new Stats_7(this, AeTime::time(), false);
      } else {
        stats_best->reinit(AeTime::time());
        stats_mean->reinit(AeTime::time());
      }

      stats_best->write_best();
      stats_mean->write_average();
    }
  } else {
    // ONLY BEST
    // Stats
#pragma omp single
    {
      if (!optim_prom) {
        stats_best = new Stats_7(this, AeTime::time(), true);
      } else {
        stats_best->reinit(AeTime::time());
      }

      best_indiv->reset_stats();
      best_indiv->metadata_->rna_begin();
      for (int i = 0; i < best_indiv->metadata_->rna_count(); i++) {
        Rna_7*rna = best_indiv->metadata_->rna_next();
        if (rna != nullptr) {
          if (rna->is_coding_)
            best_indiv->nb_coding_RNAs++;
          else
            best_indiv->nb_non_coding_RNAs++;
        }
      }


      for (int i = 0; i < best_indiv->metadata_->proteins_count(); i++) {
        Protein_7*prot = best_indiv->metadata_->proteins(i);
        if (prot != nullptr) {
          if (prot->is_functional) {
            best_indiv->nb_func_genes++;
          } else {
            best_indiv->nb_non_func_genes++;
          }
          if (prot->h > 0) {
            best_indiv->nb_genes_activ++;
          } else {
            best_indiv->nb_genes_inhib++;
          }
        }
      }


      stats_best->write_best();
    }
  }



  if (!first_gener_) {
#pragma omp single
    {

      if (ExpManager_7::standalone() && exp_m_->record_light_tree()) {
        if (ExpManager_7::standalone() && exp_m_->record_light_tree() && AeTime::time() % exp_m_->backup_step() == 0 &&
            AeTime::time() > 0) {
        }

        if (ExpManager_7::standalone() && exp_m_->record_light_tree() && AeTime::time() > 0) {
          exp_m_->output_m()->light_tree()->update_tree(AeTime::time(),
                                                        previous_individuals);

          if (AeTime::time() % exp_m_->backup_step() == 0) {
            std::cout << "writing light tree for gen : " << AeTime::time() << '\n';
            exp_m_->output_m()->write_light_tree(AeTime::time());
          }
        }

        if (ExpManager_7::standalone() && exp_m_->record_light_tree() && AeTime::time() % exp_m_->backup_step() == 0 &&
            AeTime::time() > 0) {
        }

      }


      if (ExpManager_7::standalone() && exp_m_->record_tree() && AeTime::time() %  exp_m_->tree_step() == 0 &&
          AeTime::time() > 0) {
        int status;
        status = mkdir(TREE_DIR, 0755);
        if ((status == -1) && (errno != EEXIST))
        {
          err(EXIT_FAILURE, "Impossible to create the directory %s", TREE_DIR);
        }

        char tree_file_name[50];

        sprintf(tree_file_name, "tree/tree_" TIMESTEP_FORMAT ".ae", AeTime::time());
        exp_m_->output_m()->tree()->write_to_tree_file(tree_file_name);
      }

    }
  }


  if (!first_gener_ && ExpManager_7::standalone() && !exp_m_->check_simd() && AeTime::time() % exp_m_->backup_step() == 0) {
#pragma omp single
    {
      for (int indiv_id = 0; indiv_id < (int) exp_m_->nb_indivs(); indiv_id++) {
        int x = indiv_id / exp_m_->world()->height();
        int y = indiv_id % exp_m_->world()->height();

        exp_m_->world()->grid(x, y)->individual()->clear_everything_except_dna_and_promoters();
        exp_m_->world()->grid(x, y)->individual()->genetic_unit_list_nonconst().clear();
        delete exp_m_->world()->grid(x, y)->individual();

#ifdef __REGUL
        Individual_R *indiv = new Individual_R(exp_m_,
                                                         exp_m_->world()->grid(x, y)->mut_prng(),
                                                         exp_m_->world()->grid(x, y)->stoch_prng(),
                                                         exp_m_->exp_s()->mut_params(),
                                                         w_max,
                                                         exp_m_->exp_s()->min_genome_length(),
                                                         exp_m_->exp_s()->max_genome_length(),
                                                         false,
                                                         indiv_id,
                                                         "",
                                                         0);
#else
        Individual *indiv = new Individual(exp_m_,
                                           exp_m_->world()->grid(x, y)->mut_prng(),
                                           exp_m_->world()->grid(x, y)->stoch_prng(),
                                           exp_m_->exp_s()->mut_params(),
                                           w_max,
                                           exp_m_->exp_s()->min_genome_length(),
                                           exp_m_->exp_s()->max_genome_length(),
                                           false,
                                           indiv_id,
                                           "",
                                           0);
#endif

        int32_t nb_blocks_ =
            previous_individuals[indiv_id]->dna_->nb_block();
        char *dna_string = new char[nb_blocks_ * BLOCK_SIZE];
        memset(dna_string, 0,
               (previous_individuals[indiv_id]->dna_->length() + 1) * sizeof(char));


        char *to_copy = previous_individuals[indiv_id]->dna_->to_char();


        memcpy(dna_string, to_copy,
               (previous_individuals[indiv_id]->dna_->length() + 1) * sizeof(char));


        indiv->add_GU(dna_string,
                      previous_individuals[indiv_id]->dna_->length());
        indiv->genetic_unit_nonconst(0).set_min_gu_length(exp_m_->exp_s()->min_genome_length());
        indiv->genetic_unit_nonconst(0).set_max_gu_length(exp_m_->exp_s()->max_genome_length());
        indiv->compute_statistical_data();
        indiv->EvaluateInContext(exp_m_->world()->grid(x, y)->habitat());


        exp_m_->world()->grid(x, y)->set_individual(indiv);
      }

      // Create missing directories
      exp_m_->WriteDynamicFiles();

      std::ofstream last_gener_file(LAST_GENER_FNAME,
                                    std::ofstream::out);

      last_gener_file << AeTime::time() << std::endl;
      last_gener_file.close();
    }
  } else {
#pragma omp single
    {
      first_gener_ = false;
    }
  }

  if (AeTime::time() == exp_m_->end_step() && ExpManager_7::standalone() && !exp_m_->check_simd()) {
#pragma omp for schedule(static)
    for (int indiv_id = 0; indiv_id < (int) exp_m_->nb_indivs(); indiv_id++) {
      int x = indiv_id / exp_m_->world()->height();
      int y = indiv_id % exp_m_->world()->height();

      exp_m_->world()->grid(x, y)->individual()->clear_everything_except_dna_and_promoters();
      exp_m_->world()->grid(x, y)->individual()->genetic_unit_list_nonconst().clear();
      delete exp_m_->world()->grid(x, y)->individual();
    }
  }



}


void ExpManager_7::check_dna() {

  int x, y;
  for (int i = 0; i < (int) exp_m_->nb_indivs(); i++) {
    x = i / exp_m_->world()->height();
    y = i % exp_m_->world()->height();

    for (int dna_pos = 0; dna_pos < dna_size[i]; dna_pos++) {
      if (exp_m_->world()->grid(x, y)->individual()->genetic_unit(
          0).dna()->data()[dna_pos] !=
          current_individuals[i]->dna_->data_[dna_pos]) {

        printf("Check DNA indiv %d %d %d --- NB Mutation %ld\n",i,dna_size[i],exp_m_->world()->grid(x, y)->individual()->genetic_unit(
            0).dna()->length(),exp_m_->dna_mutator_array_[i]->mutation_list_.size());

        printf("Divergence between classic DNA and SIMD DNA %d %d at pos %d\n",
               exp_m_->world()->grid(x, y)->individual()->genetic_unit(
                   0).dna()->data()[dna_pos],
               current_individuals[i]->dna_->data_[dna_pos],
               dna_pos);
        for (auto mute : exp_m_->dna_mutator_array_[i]->mutation_list_) {
          printf("Mutation type %d\n",mute->type());
        }

        break;
      }
    }
  }
}

void ExpManager_7::check_individual(int i, int x, int y) {
  exp_m_->world()->grid(x, y)->set_individual(exp_m_->world()->grid(x, y)->old_one);
  exp_m_->world()->grid(x, y)->old_one->Reevaluate();

  printf("%d %d %d -- ",i,x,y);

  printf(
      "Nb RNA SIMD/CPU %ud/%ld Protein %ud/%ld Metaerror %f/%f Fitness %e/%e DNA Size %d/%d\n",
      previous_individuals[i]->metadata_->rna_count(),
      exp_m_->world()->grid(x, y)->individual()->rna_list().size(),
      previous_individuals[i]->metadata_->proteins_count(),
      exp_m_->world()->grid(x, y)->individual()->protein_list().size(),
      previous_individuals[i]->metaerror,
      exp_m_->world()->grid(x, y)->individual()->dist_to_target_by_feature(
          METABOLISM),
      previous_individuals[i]->fitness,
      exp_m_->world()->grid(x, y)->individual()->fitness(),
      previous_individuals[i]->dna_->length(),
      exp_m_->world()->grid(x, y)->individual()->genetic_unit(
          0).seq_length());

  int idx = 0;

  for (auto rna : exp_m_->world()->grid(x, y)->old_one->rna_list()) {
    printf("RNA CPU %d Start %d Stop %d Leading/Lagging %d Length %d\n", idx,
           rna->promoter_pos(), rna->last_transcribed_pos(), rna->strand(), rna->transcript_length());
    idx++;
  }
  idx = 0;
  for (idx = 0; idx < (previous_individuals[i]->metadata_->promoter_count()); idx++) {
    if (previous_individuals[i]->metadata_->promoters(idx) != nullptr)
      printf("Promoters found at %d\n",
             previous_individuals[i]->metadata_->promoters(idx)->pos);
  }

  idx = 0;
  for (idx = 0; idx < (previous_individuals[i]->metadata_->rna_count()); idx++) {
    if (previous_individuals[i]->metadata_->rnas(idx) != nullptr) {
      printf("RNA SIMD %d Start %d Stop %d Leading/Lagging %d Length %d\n", idx,
             previous_individuals[i]->metadata_->rnas(idx)->begin,
             previous_individuals[i]->metadata_->rnas(idx)->end,
             previous_individuals[i]->metadata_->rnas(idx)->leading_lagging,
             previous_individuals[i]->metadata_->rnas(idx)->length);
    }
  }

  int prot_cpt_b=0;
  idx = 0;
  for (auto prot : exp_m_->world()->grid(x, y)->old_one->protein_list()) {
    bool found = false;

    for (int pidx = 0; pidx <
                       (int)previous_individuals[i]->metadata_->proteins_count(); pidx++) {
      if (previous_individuals[i]->metadata_->proteins(pidx)->is_init_) {
        if ((previous_individuals[i]->metadata_->proteins(pidx)->e == prot->concentration()) &&
            (previous_individuals[i]->metadata_->proteins(pidx)->protein_end == prot->last_STOP_base_pos())) {
          found = true;
          break;
        }
      }
    }

    if (!found) {
      printf("Proteins CPU %d Start %d (end %d stop %d) Length %d Leading/Lagging %d M/W/H %f/%f/%f Func %d -- Concentration %f\n",
             idx,
             prot->first_translated_pos(), prot->last_translated_pos(), prot->last_STOP_base_pos(),
             prot->length(), prot->strand(),
             prot->mean(), prot->width(), prot->height(), prot->is_functional(), prot->concentration());
    }
    idx++;
  }

  for (int idx = 0; idx <
                    (int)previous_individuals[i]->metadata_->proteins_count(); idx++) {
    if (previous_individuals[i]->metadata_->proteins(idx)->is_init_) {


      bool found = false;

      for (auto prot : exp_m_->world()->grid(x, y)->old_one->protein_list()) {
        if ((previous_individuals[i]->metadata_->proteins(idx)->e ==  prot->concentration()) &&
            (previous_individuals[i]->metadata_->proteins(idx)->protein_end ==  prot->last_STOP_base_pos())) {
          found = true;
          break;
        }
      }

      if (!found)
        printf("Proteins SIMD %d Start %d (end %d) Length %d Leading/Lagging %d M/W/H %f/%f/%f Func %d -- Concentration %f\n", idx,
            previous_individuals[i]->metadata_->proteins(idx)->protein_start,
            previous_individuals[i]->metadata_->proteins(idx)->protein_end,
            previous_individuals[i]->metadata_->proteins(idx)->protein_length,
            previous_individuals[i]->metadata_->proteins(idx)->leading_lagging,
            previous_individuals[i]->metadata_->proteins(idx)->m,
            previous_individuals[i]->metadata_->proteins(idx)->w,
            previous_individuals[i]->metadata_->proteins(idx)->h,
            previous_individuals[i]->metadata_->proteins(idx)->is_functional,
            previous_individuals[i]->metadata_->proteins(idx)->e
        );
      prot_cpt_b++;
    }
  }

}

void ExpManager_7::check_struct() {
  for (int i = 0; i < (int) exp_m_->nb_indivs(); i++) {

  }
}

void ExpManager_7::check_result() {


#pragma omp single
  {
    printf("Check results !!!!!\n");
    int x, y;

    bool validated_generation = true;

    for (int i = 0; i < (int) exp_m_->nb_indivs(); i++) {
      //if (i != 905) continue;

      x = i / exp_m_->world()->height();
      y = i % exp_m_->world()->height();

      double fit_1 = exp_m_->world()->grid(x,
                                           y)->individual()->dist_to_target_by_feature(
          METABOLISM);
      double fit_2 = previous_individuals[i]->metaerror;
      float i_fit_1 = roundf(fit_1 * 100);
      float i_fit_2 = roundf(fit_2 * 100);

      int count_prot = 0;

      for (int pidx = 0; pidx < previous_individuals[i]->metadata_->proteins_count(); pidx++) {
        if (previous_individuals[i]->metadata_->proteins(pidx)->is_init_) {
          count_prot++;
        }
      }

      int count_rna_cpu = 0;

      for (auto rna : exp_m_->world()->grid(x, y)->individual()->rna_list()) {
        if (rna->transcript_length() >= 0) {
          count_rna_cpu++;
        }
      }


//      int idx = 0, fidx = 0;
//      for (auto prot : exp_m_->world()->grid(x, y)->individual()->protein_list()) {
//        bool found = false;
//        fidx = 0;
//
//        for (int pidx = 0; pidx < previous_individuals[i]->metadata_->proteins_count(); pidx++) {
//          if (previous_individuals[i]->metadata_->proteins(pidx)->is_init_) {
//            if ((previous_individuals[i]->metadata_->proteins(pidx)->e ==
//                 prot->concentration()) &&
//                (previous_individuals[i]->metadata_->proteins(pidx)->protein_end ==
//                 prot->last_STOP_base_pos()))
//              if ((previous_individuals[i]->metadata_->proteins(pidx)->protein_length ==
//                   prot->length()) &&
//                  (previous_individuals[i]->metadata_->proteins(pidx)->protein_start ==
//                   prot->first_translated_pos())) {
//                found = true;
//                fidx = pidx;
//                break;
//              } else {
//                fidx = pidx;
//              }
//
//          }
//        }
//
//        if (!found) {
//          printf("==================-------------------------======================\n");
//          printf("Proteins CPU %d Start %d (end %d stop %d) Length %d Leading/Lagging %d M/W/H %f/%f/%f Func %d -- Concentration %f\n",
//                 idx,
//                 prot->first_translated_pos(), prot->last_translated_pos(),
//                 prot->last_STOP_base_pos(),
//                 prot->length(), prot->strand(),
//                 prot->mean(), prot->width(), prot->height(), prot->is_functional(),
//                 prot->concentration());
//
//          if (fidx < previous_individuals[i]->metadata_->proteins_count())
//            printf("Proteins SIMD %d Start %d (end %d) Length %d Leading/Lagging %d M/W/H %f/%f/%f Func %d -- Concentration %f\n",
//                   fidx,
//                previous_individuals[i]->metadata_->proteins(fidx)->protein_start,
//                previous_individuals[i]->metadata_->proteins(fidx)->protein_end,
//                previous_individuals[i]->metadata_->proteins(fidx)->protein_length,
//                previous_individuals[i]->metadata_->proteins(fidx)->leading_lagging,
//                previous_individuals[i]->metadata_->proteins(fidx)->m,
//                previous_individuals[i]->metadata_->proteins(fidx)->w,
//                previous_individuals[i]->metadata_->proteins(fidx)->h,
//                previous_individuals[i]->metadata_->proteins(fidx)->is_functional,
//                previous_individuals[i]->metadata_->proteins(fidx)->e
//            );
//          printf("==================-------------------------======================\n");
//        }
//        idx++;
//      }
//
//
//      for (int j = 0; j < previous_individuals[i]->dna_->length(); j++) {
//        if (previous_individuals[i]->dna_->data_[j] !=
//            exp_m_->world()->grid(x, y)->individual()->genetic_unit(0).dna()->data()[j]) {
//          printf("%ld -- %d -- DNA is different at %d !!!\n", AeTime::time(), i, j);
//
//          exit(-1);
//        }
//      }

      int prot_size = (int) exp_m_->world()->grid(x, y)->individual()->protein_list().size();
      if ((((previous_individuals[i]->metadata_->rna_count() !=
            count_rna_cpu) ||
           (count_prot != prot_size))
          || ((i_fit_1 != i_fit_2))) && dna_size[i] > 300) {
        validated_generation = false;


        printf(
            "X-X-ERROR -- %ld -- Individual %d  -- %d %d --(P %d / %d): Metaerror (CPU/GPU) : %e/%e || Fitness (CPU/GPU) : %e/%e \n",
            AeTime::time(), i, x,y,exp_m_->world()->grid(x, y)->individual()->parent_id_,
               previous_individuals[i]->parent_id,
            exp_m_->world()->grid(x, y)->individual()->dist_to_target_by_feature(
                METABOLISM),
               previous_individuals[i]->metaerror,
            exp_m_->world()->grid(x, y)->individual()->fitness(),
               previous_individuals[i]->fitness);

        printf(
            "Nb RNA SIMD/CPU %ud/%ld Protein %ud/%ld Metaerror %f/%f Fitness %e/%e DNA Size %d/%d\n",
               previous_individuals[i]->metadata_->rna_count(),
            exp_m_->world()->grid(x, y)->individual()->rna_list().size(),
            count_prot,
            exp_m_->world()->grid(x, y)->individual()->protein_list().size(),
               previous_individuals[i]->metaerror,
            exp_m_->world()->grid(x, y)->individual()->dist_to_target_by_feature(
                METABOLISM),
               previous_individuals[i]->fitness,
            exp_m_->world()->grid(x, y)->individual()->fitness(), dna_size[i],
            exp_m_->world()->grid(x, y)->individual()->genetic_unit(
                0).seq_length());

        printf("Promoters LEADING : ");
        for (auto& prom : ((List_Metadata*)previous_individuals[i]->
            metadata_)->promoters_list_[LEADING]) {
          printf("%d ",prom.pos);
        }
        printf("\n");
        printf("Promoters LAGGING : ");
        for (auto& prom : ((List_Metadata*)previous_individuals[i]->
            metadata_)->promoters_list_[LAGGING]) {
          printf("%d ",prom.pos);
        }
        printf("\n");

        int idx = 0;

        for (auto rna : exp_m_->world()->grid(x, y)->individual()->rna_list()) {
          bool found = false;
          // for (int pidx = 0; pidx < (int) (previous_individuals[i]->metadata_->rna_count());
          //      pidx++) {
          //   if ((rna->promoter_pos() ==
          //        previous_individuals[i]->metadata_->rnas(pidx)->begin) &&
          //       (rna->transcript_length() ==
          //        previous_individuals[i]->metadata_->rnas(pidx)->length)){
          //     found = true;
          //     break;
          //   }
          // }

          // if (i == 392) found = false;

          if (!found)
            printf("RNA CPU %d Start %d Stop %d Leading/Lagging %d Length %d Basal %lf\n", idx,
                   rna->promoter_pos(), rna->last_transcribed_pos(), rna->strand(),
                   rna->transcript_length(),rna->basal_level());
          idx++;
        }

        idx = 0;
        for (idx = 0; idx < (int) (previous_individuals[i]->metadata_->rna_count()); idx++) {
          bool found = false;
          // for (auto rna : exp_m_->world()->grid(x, y)->individual()->rna_list()) {
          //   if ((rna->promoter_pos() ==
          //        previous_individuals[i]->metadata_->rnas(idx)->begin) &&
          //       (rna->transcript_length() ==
          //        previous_individuals[i]->metadata_->rnas(idx)->length)) {
          //     found = true;
          //     break;
          //   }
          // }

          // if (i == 392) found = false;

          if (!found)
            printf("RNA SIMD %d Start %d Stop %d Leading/Lagging %d Length %d  Basal %lf\n", idx, previous_individuals[i]->metadata_->rnas(idx)->begin,
                previous_individuals[i]->metadata_->rnas(idx)->end,
                previous_individuals[i]->metadata_->rnas(idx)->leading_lagging,
                previous_individuals[i]->metadata_->rnas(idx)->length,
                previous_individuals[i]->metadata_->rnas(idx)->e);
        }


        idx = 0;
        int prot_cpt_b = 0;

        for (auto prot : exp_m_->world()->grid(x, y)->individual()->protein_list()) {
          bool found = false;

          // for (int pidx = 0; pidx < previous_individuals[i]->metadata_->proteins_count(); pidx++) {
          //   if (previous_individuals[i]->metadata_->proteins(pidx)->is_init_) {
          //     if ((previous_individuals[i]->metadata_->proteins(pidx)->e ==
          //          prot->concentration()) &&
          //         (previous_individuals[i]->metadata_->proteins(pidx)->protein_end ==
          //          prot->last_STOP_base_pos())) {
          //       found = true;
          //       break;
          //     }
          //   }
          // }

          // if (i == 392) found = false;
#ifdef __REGUL
          if (!found && !((Protein_R*)prot)->is_signal()) {
#else
          if (!found) {
#endif
            printf("Proteins CPU %d Start %d (end %d stop %d) Length %d Leading/Lagging %d M/W/H %f/%f/%f Func %d -- Concentration %f RNA : \n",
                   idx,
                   prot->first_translated_pos(), prot->last_translated_pos(),
                   prot->last_STOP_base_pos(),
                   prot->length(), prot->strand(),
                   prot->mean(), prot->width(), prot->height(), prot->is_functional(),
                   prot->concentration());

            for (auto rna : prot->rna_list()) {
              printf("[%d => %d]\n",rna->promoter_pos(),rna->last_transcribed_pos());
            }
          }
          idx++;
        }

        for (int idx = 0; idx < previous_individuals[i]->metadata_->proteins_count(); idx++) {
          if (previous_individuals[i]->metadata_->proteins(idx)->is_init_) {


            bool found = false;

            // for (auto prot : exp_m_->world()->grid(x, y)->individual()->protein_list()) {
            //   if ((previous_individuals[i]->metadata_->proteins(idx)->e ==
            //        prot->concentration()) &&
            //       (previous_individuals[i]->metadata_->proteins(idx)->protein_end ==
            //        prot->last_STOP_base_pos())) {
            //     found = true;
            //     break;
            //   }
            // }

            // if (i == 392) found = false;

            //for (idx = 0; idx < (int) (current_individuals[i]->proteins.size()); idx++) {
            if (!found) {
              printf("Proteins SIMD %d Start %d (end %d) Length %d Leading/Lagging %d M/W/H %f/%f/%f Func %d -- Concentration %f\n",
                     idx,
                     previous_individuals[i]->metadata_->proteins(idx)->protein_start,
                     previous_individuals[i]->metadata_->proteins(idx)->protein_end,
                     previous_individuals[i]->metadata_->proteins(idx)->protein_length,
                     previous_individuals[i]->metadata_->proteins(idx)->leading_lagging,
                     previous_individuals[i]->metadata_->proteins(idx)->m,
                     previous_individuals[i]->metadata_->proteins(idx)->w,
                     previous_individuals[i]->metadata_->proteins(idx)->h,
                     previous_individuals[i]->metadata_->proteins(idx)->is_functional,
                     previous_individuals[i]->metadata_->proteins(idx)->e
              );

                    for (auto rna : previous_individuals[i]->metadata_->proteins(idx)->rna_list_) {
                      printf("[%d => %d]\n",rna->begin,rna->end);
                    }
            
              validated_generation = false;
            }
            prot_cpt_b++;
          }
        }

        // printf("Start prot LEAD : ");
        // for (int pidx = 0; pidx < (int) (previous_individuals[i]->metadata_->rna_count());
        //      pidx++) {
        //   if (previous_individuals[i]->metadata_->rnas(pidx)->leading_lagging == 0) {
        //     for (int pos :
        //          previous_individuals[i]->metadata_->rnas(pidx)->start_prot) {
        //       printf("%d ",pos);
        //     }
        //   }
        // }
        // printf("\n");


        // printf("Start prot LAG : ");
        // for (int pidx = 0; pidx < (int) (previous_individuals[i]->metadata_->rna_count());
        //      pidx++) {
        //   if (previous_individuals[i]->metadata_->rnas(pidx)->leading_lagging == 1) {
        //     for (int pos :
        //          previous_individuals[i]->metadata_->rnas(pidx)->start_prot) {
        //       printf("%d ",pos);
        //     }
        //   }
        // }
        // printf("\n");

        // exp_m_->world()->grid(x, y)->individual()->phenotype()->print();

        // previous_individuals[i]->phenotype->print();

        exit(-1);
      }

    }


    if (validated_generation)
      printf("Generation %ld is replicated with SIMD without diff\n", AeTime::time());


  }
}
}