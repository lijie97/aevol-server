//
// Created by arrouan on 10/02/16.
//

#include "cuda_struct.h"

namespace aevol {

void cuda_struct::init_struct(int max_protein, int max_rna, int max_influence,
                              int nb_signals, int life_time, int nb_eval,
                              float selection_pressure) {
  max_protein_ = max_protein;
  max_rna_ = max_rna;
  max_influence_ = max_influence;
  nb_signals_ = nb_signals;
  life_time_ = life_time;
  nb_eval_ = nb_eval;
  selection_pressure_ = selection_pressure;

  phenotype_inhib = new float[300 * 1024];
  phenotype_activ = new float[300 * 1024];
  phenotype = new float[300 * 1024];
  dist_sum = new float[1024];

  environment = new float[300 * life_time];
  delta = new float[300 * 1024];

  eval_step = new int[life_time];

  protein_concentration = new float[max_protein * 1024];

  signals_concentration = new float[nb_signals * life_time];

  enhance_coef = new float[max_influence * max_rna * 1024];
  operate_coef = new float[max_influence * max_rna * 1024];
  rna_synthesis = new float[max_rna * 1024];
  basal_level = new float[max_rna * 1024];

  protein_influence = new int[max_influence * max_rna * 1024];

  protein_influenced = new float[max_protein * max_rna * 1024];

  protein_triangle_ix0 = new int[max_protein * 1024];
  protein_triangle_ix1 = new int[max_protein * 1024];
  protein_triangle_ix2 = new int[max_protein * 1024];

  protein_triangle_height = new float[max_protein * 1024];
}

void cuda_struct::transfert_to_gpu(ExpManager* exp_m) {
  std::set<int>* eval = exp_m->exp_s()->get_list_eval_step();
  const Habitat_R& hab = dynamic_cast<const Habitat_R&>(exp_m->world()->grid(0,0)->habitat());

  for (int8_t i = 1; i <= exp_m->exp_s()->get_nb_indiv_age(); i++) {
    //Set the concentration of signals for this age
    for (Protein_R* prot1 : hab.signals()) {
      signals_concentration[i * life_time_ + prot1->get_local_id()] = 0.0;
    }

    for (Protein_R* prot2 : hab.phenotypic_target(i).signals()) {
      signals_concentration[i * life_time_ + prot2->get_local_id()] = 0.9;
    }

    if (eval->find(i) != eval->end())
      eval_step[i] = 1;
    else
      eval_step[i] = 0;

    for (int j = 0; j < 300; j++)
      environment[(i-1) * 300 + j] = ((HybridFuzzy*) hab.phenotypic_target(
          i).fuzzy())->points()[j];
  }


  int i = 0;
  for (int ki = 0; ki < 32; ki++)
    for (int kj = 0; kj < 32; kj++) {
      Individual_R* indiv = dynamic_cast<Individual_R*>(exp_m->world()->indiv_at(ki, kj));

      i = ki * 32 + kj;

      for (int j = 0; j < 300; j++) {
        phenotype[i * 300 + j] = 0;
        phenotype_activ[i * 300 + j] = 0;
        phenotype_inhib[i * 300 + j] = 0;
        delta[i * 300 + j] = 0;
      }

      for (int je = 0; je < max_protein_; je++) {
        protein_concentration[i * max_protein_ + je] + 0;

        protein_triangle_ix0[i * max_protein_ + je] = 0;
        protein_triangle_ix1[i * max_protein_ + je] = 0;
        protein_triangle_ix2[i * max_protein_ + je] = 0;

        protein_triangle_height[i * max_protein_ + je] = 0;


        for (int k = 0; k < max_protein_ * max_rna_; k++)
          protein_influenced[i * max_protein_ * max_rna_ + je + k] = 0;
      }

      for (auto prot_a : indiv->_initial_protein_list) {
        Protein_R* prot = dynamic_cast<Protein_R*>(prot_a);

        protein_concentration[i * max_protein_ +
                              prot->get_local_id()] = prot->concentration();

        double x0 = prot->mean() - prot->width();
        double x1 = prot->mean();
        double x2 = prot->mean() + prot->width();

        int ix0 = (int) (x0 * 300);
        int ix1 = (int) (x1 * 300);
        int ix2 = (int) (x2 * 300);

        if (ix0 < 0) ix0 = 0; else if (ix0 > (300 - 1)) ix0 = 300 - 1;
        if (ix1 < 0) ix1 = 0; else if (ix1 > (300 - 1)) ix1 = 300 - 1;
        if (ix2 < 0) ix2 = 0; else if (ix2 > (300 - 1)) ix2 = 300 - 1;

        protein_triangle_ix0[i * max_protein_ + prot->get_local_id()] = ix0;
        protein_triangle_ix1[i * max_protein_ + prot->get_local_id()] = ix1;
        protein_triangle_ix2[i * max_protein_ + prot->get_local_id()] = ix2;

        protein_triangle_height[i * max_protein_ +
                                prot->get_local_id()] = prot->height();


        for (auto rna : prot->_rna_R_list)
          protein_influenced[i * max_protein_ * max_rna_ +
                             prot->get_local_id() * max_protein_
                             + rna->get_local_id()] = 1;
      }

      for (int je = 0; je < max_rna_; je++) {
        for (int ke = 0; ke < max_influence_; ke++) {
          enhance_coef[i * max_rna_ * max_influence_ +
                       je * max_influence_ + ke] = 0;
          operate_coef[i * max_rna_ * max_influence_ +
                       je * max_influence_ + ke] = 0;
          protein_influence[i * max_rna_ * max_influence_ +
                            je * max_influence_ + ke] = 0;
        }
      }


      for (auto rna : indiv->_rna_list_coding) {

        rna_synthesis[i * max_rna_ +
                      rna->get_local_id()] = rna->get_synthesis_rate();

        basal_level[i * max_rna_ +
                      rna->get_local_id()] = rna->basal_level();

        for (int k = 0; k < rna->_nb_influences; k++) {
          enhance_coef[i * max_rna_ * max_influence_ +
                       rna->get_local_id() * max_influence_ +
                       k] = rna->_enhancing_coef_list[k];
          operate_coef[i * max_rna_ * max_influence_ +
                       rna->get_local_id() * max_influence_ +
                       k] = rna->_operating_coef_list[k];
          if (rna->_protein_list[k]->is_signal())
            protein_influence[i * max_rna_ * max_influence_ +
                              rna->get_local_id() * max_influence_ +
                              k] =
                rna->_protein_list[k]->get_local_id() + max_protein_;
          else
            protein_influence[i * max_rna_ * max_influence_ +
                              rna->get_local_id() * max_influence_ +
                              k] = rna->_protein_list[k]->get_local_id();
        }
      }

    }

  /** TRANSFER Vector **/

}

void cuda_struct::compute_a_generation(ExpManager* exp_m) {

  int degradation_step = exp_m->exp_s()->get_nb_degradation_step();
  float degradation_rate = exp_m->exp_s()->get_degradation_rate();
  float hill_shape_n = exp_m->exp_s()->get_hill_shape_n();
  float hill_shape = exp_m->exp_s()->get_hill_shape();

  for (int ki = 0; ki < 32; ki++)
    for (int kj = 0; kj < 32; kj++) {
      int indiv_id = ki * 32 + kj;

      for (int8_t i = 1; i <= life_time_; i++) {
        for (int j = 0; j < degradation_step; j++) {
          // Compute synthesis rate of RNA
          for (int m = 0; m < max_rna_; m++) {
            float enhancer_activity = 0, operator_activity = 0;
            for (int n = 0; n < max_influence_; n++) {
              int prot_id = protein_influence
              [ki * 32 + kj * max_rna_ * max_influence_ + m * max_influence_ +
               n];

              if (prot_id > max_protein_) {
                enhancer_activity +=
                    enhance_coef[indiv_id * max_rna_ * max_influence_ +
                                 m * max_influence_ + n] *
                    signals_concentration[i * life_time_ + prot_id -
                                          max_protein_];

                operator_activity +=
                    operate_coef[indiv_id * max_rna_ * max_influence_ +
                                 m * max_influence_ + n] *
                    signals_concentration[i * life_time_ + prot_id -
                                          max_protein_];
              } else {
                enhancer_activity +=
                    enhance_coef[indiv_id * max_rna_ * max_influence_ +
                                 m * max_influence_ + n] *
                    protein_concentration[indiv_id * max_protein_ + prot_id];

                operator_activity +=
                    operate_coef[indiv_id * max_rna_ * max_influence_ +
                                 m * max_influence_ + n] *
                    protein_concentration[indiv_id * max_protein_ + prot_id];
              }
            }

            float enhancer_activity_pow_n = pow(enhancer_activity,
                                                hill_shape_n);
            float operator_activity_pow_n = pow(operator_activity,
                                                hill_shape_n);

            rna_synthesis[indiv_id * max_rna_ + m] = basal_level[indiv_id * max_rna_ + m]
                                                         * (hill_shape
                                                            /
                                                            (operator_activity_pow_n +
                                                             hill_shape))
                                                         * (1 + ((1 /
                basal_level[indiv_id * max_rna_ + m]) -
                                                                 1)
                                                                *
                                                                (enhancer_activity_pow_n /
                                                                 (enhancer_activity_pow_n +
                                                                  hill_shape)));
          }

          // Compute concentration
          for (int l = 0; l < max_protein_; l++) {
            float _delta_concentration;
            for (int m = 0; m < max_rna_; m++) {
              if (protein_influenced[indiv_id * max_protein_ * max_rna_ +
                                     l * max_protein_
                                     + m])
                _delta_concentration += rna_synthesis[indiv_id * max_rna_ +
                                                      m];
            }

            _delta_concentration -= degradation_rate * protein_concentration[indiv_id * max_protein_ +
                                                                             l];
            _delta_concentration *= 1 / ((float) degradation_step);

            protein_concentration[indiv_id * max_protein_ +
                                  l] += _delta_concentration;
          }
        }

        // If we have to evaluate the individual at this age
        if (eval_step[i]) {
          /** update phenotype **/
          for (int j = 0; j < 300; j++) {
            phenotype[indiv_id * 300 + j] = 0;
            phenotype_activ[indiv_id * 300 + j] = 0;
            phenotype_inhib[indiv_id * 300 + j] = 0;
          }

          for (int l = 0; l < max_protein_; l++) {
            float height = (protein_triangle_height[indiv_id * max_protein_ +
                                                    l] *
                            protein_concentration[indiv_id * max_protein_ + l]);
            float incY =
                protein_triangle_height[indiv_id * max_protein_ + l] *
                protein_concentration[indiv_id * max_protein_ + l] /
                (protein_triangle_ix1[indiv_id * max_protein_ + l] -
                 protein_triangle_ix0[indiv_id * max_protein_ + l]);

            if (protein_triangle_height[indiv_id * max_protein_ + l] > 0) {
              for (int j = 0; j < 300; j++) {
                if (j > protein_triangle_ix0[indiv_id * max_protein_ + l]
                    and
                    j < protein_triangle_ix1[indiv_id * max_protein_ + l])
                  phenotype_activ[indiv_id * 300 + j] +=
                      incY * (j -
                              protein_triangle_ix0[indiv_id * max_protein_ +
                                                   l]);
                else if (
                    j > protein_triangle_ix1[indiv_id * max_protein_ + l]
                    and
                    j < protein_triangle_ix2[indiv_id * max_protein_ + l])
                  phenotype_activ[indiv_id * 300 + j] += height
                                                         -
                                                         incY * (j -
                                                                 protein_triangle_ix1[
                                                                     indiv_id *
                                                                     max_protein_ +
                                                                     l]);
                else if (j ==
                         protein_triangle_ix1[indiv_id * max_protein_ + l])
                  phenotype_activ[indiv_id * 300 + j] +=
                      height;
              }
            } else {
              for (int j = 0; j < 300; j++) {
                if (j > protein_triangle_ix0[indiv_id * max_protein_ + l]
                    and
                    j < protein_triangle_ix1[indiv_id * max_protein_ + l])
                  phenotype_inhib[indiv_id * 300 + j] += incY * (j -
                                                                 protein_triangle_ix0[
                                                                     indiv_id *
                                                                     max_protein_ +
                                                                     l]);
                else if (
                    j > protein_triangle_ix1[indiv_id * max_protein_ + l]
                    and
                    j < protein_triangle_ix2[indiv_id * max_protein_ + l])
                  phenotype_inhib[indiv_id * 300 + j] +=
                      height -
                      incY * (j -
                              protein_triangle_ix1[indiv_id * max_protein_ +
                                                   l]);
                else if (j ==
                         protein_triangle_ix1[indiv_id * max_protein_ + l])
                  phenotype_inhib[indiv_id * 300 + j] +=
                      height;
              }
            }
          }

          for (int j = 0; j < 300; j++) {
            phenotype_activ[indiv_id * 300 + j] =
                phenotype_activ[indiv_id * 300 + j] > Y_MAX ? Y_MAX :
                phenotype_activ[indiv_id * 300 + j];
            phenotype_inhib[indiv_id * 300 + j] =
                phenotype_inhib[indiv_id * 300 + j] < -Y_MAX ? -Y_MAX :
                phenotype_inhib[indiv_id * 300 + j];
            phenotype[indiv_id * 300 + j] =
                phenotype_activ[indiv_id * 300 + j] +
                phenotype_inhib[indiv_id * 300 + j];
            phenotype[indiv_id * 300 + j] =
                phenotype[indiv_id * 300 + j] > Y_MIN ?
                Y_MIN : phenotype[indiv_id * 300 + j];
            delta[indiv_id * 300 + j] =
                phenotype[indiv_id * 300 + j] -
                environment[i * life_time_ + j];
          }

          /** compute distance to target **/

          float area = 0;

          for (int j = 0; j < 299; j++) {
            dist_sum[indiv_id] += ((delta[indiv_id * 300 + j]
                                    - delta[indiv_id * 300 + j + 1]) / 600.0);
          }
        }
      }

      dist_sum[indiv_id] = exp(
          selection_pressure_ * (dist_sum[indiv_id] / nb_eval_));
    }
}

void cuda_struct::print_dist(ExpManager* exp_m) {

}

}
