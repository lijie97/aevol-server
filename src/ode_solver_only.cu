#include <cuda.h>

#include "ode_solver_only_gpu.h"

__global__
void process_delta(int nb_signal, int degradstep, int degradrate, int ***rna_produce_protein_array,
                   int **nb_rna_produce_protein, int *nb_rna_produce,   double **protein_concentration_array,
                   double **rna_basal_concentration_array, int *nb_protein_array, int *nb_rna_array,
                   double ***rna_influence_enhancing_coef_array, double ***rna_influence_operating_coef_array,
                   int **nb_rna_influence_enhancing_coef, int  **nb_rna_influence_operating_coef,
                   double **env_concentration_array, double hill_shape, double hill_shape_n) {



  double delta = 0;

  int indiv_id = blockIdx.x;
  int prot_id = threadIdx.x;

  if (prot_id < nb_protein_array[indiv_id] - nb_signal) {
    for (int j = 0; j < degradstep; j++) {
      for (int j = 0;
           j < nb_rna_produce_protein[indiv_id][prot_id]; j++) {
        double enhancer_activity = 0;
        double operator_activity = 0;

        int rna_id = rna_produce_protein_array[indiv_id][prot_id][j];

        for (int i = 0; i <
                        nb_rna_influence_enhancing_coef[indiv_id][rna_id]; i++) {

          enhancer_activity +=
              rna_influence_enhancing_coef_array[indiv_id][rna_id][i]
              * protein_concentration_array[indiv_id][i];
          operator_activity +=
              rna_influence_operating_coef_array[indiv_id][rna_id][i]
              * protein_concentration_array[indiv_id][i];
        }

        double enhancer_activity_pow_n = enhancer_activity == 0 ? 0 :
                                         powf(enhancer_activity, hill_shape_n);
        double operator_activity_pow_n = operator_activity == 0 ? 0 :
                                         powf(operator_activity, hill_shape_n);
        delta += rna_basal_concentration_array[indiv_id][rna_id]
                 * (hill_shape
                    / (operator_activity_pow_n + hill_shape))
                 * (1 +
                    ((1 / rna_basal_concentration_array[indiv_id][rna_id]
                     ) -
                     1)
                    * (enhancer_activity_pow_n /
                       (enhancer_activity_pow_n + hill_shape)));
      }

      delta -=
          degradrate *
          protein_concentration_array[indiv_id][prot_id];
      delta *= 1 / (double) degradstep;

      __syncthreads();

      protein_concentration_array[indiv_id][prot_id] += delta[prot_id];
    }
  }
}
