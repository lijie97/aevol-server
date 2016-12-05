
#include "ode_solver_only_gpu.h"

void call_kernel_ode_cuda(int multiply_population, int max_protein,
                          int nb_signal, int degradationstep,
                          int degradation_rate) {
  process_delta<<<1024*multiply_population,max_protein>>>(nb_signal,degradationstep,degradation_rate,
      dev_rna_produce_protein_array, dev_nb_rna_produce_protein, dev_nb_rna_produce, dev_protein_concentration_array,
      dev_rna_basal_concentration_array, dev_nb_protein_array, dev_nb_rna_array,
      dev_rna_influence_enhancing_coef_array, dev_rna_influence_operating_coef_array,
      dev_nb_rna_influence_enhancing_coef, dev_nb_rna_influence_operating_coef,
      dev_env_concentration_array,hill_shape,hill_shape_n);
}


int transfert_data_to_gpu(int pop_size, int lifestep) {
  protein_concentration_array = (double**)malloc(pop_size * sizeof(double*));
  rna_basal_concentration_array = (double**)malloc(pop_size * sizeof(double*));
  rna_produce_protein_array = (int***)malloc(pop_size * sizeof(int**));
  rna_influence_enhancing_coef_array = (double***)malloc(pop_size * sizeof(double**));
  rna_influence_operating_coef_array = (double***)malloc(pop_size * sizeof(double**));

  nb_protein_array = (int*)malloc(pop_size * sizeof(int));
  nb_rna_array = (int*)malloc(pop_size * sizeof(int));
  nb_rna_produce_protein = (int**)malloc(pop_size * sizeof(int*));
  nb_rna_produce = (int*)malloc(pop_size * sizeof(int));
  nb_rna_influence_enhancing_coef = (int**) malloc(pop_size * sizeof(int*));
  nb_rna_influence_operating_coef = (int**) malloc(pop_size * sizeof(int*));
  nb_rna_influence_enhancing_coef_l1 = (int*) malloc(pop_size * sizeof(int));
  nb_rna_influence_operating_coef_l1 = (int*) malloc(pop_size * sizeof(int));

  int max_prot = 0;

  for (int i = 0; i < pop_size; i++){
    protein_concentration_array[i] = (double*)
        malloc(protein_concentration_list[i]->size() * sizeof(double));
    nb_protein_array[i] = (int) protein_concentration_list[i]->size();

    max_prot = protein_concentration_list[i]->size() > max_prot ?
               protein_concentration_list[i]->size() : max_prot;

    for (int prot_id = 0; prot_id < protein_concentration_list[i]->size(); prot_id++) {
      protein_concentration_array[i][prot_id] =
          protein_concentration_list[i]->at(prot_id);
    }

    rna_basal_concentration_array[i] = (double*)
        malloc(rna_basal_concentration_list[i]->size() * sizeof(double));
    nb_rna_array[i] = (int)rna_basal_concentration_list[i]->size();

    for (int rna_id = 0; rna_id < rna_basal_concentration_list[i]->size(); rna_id++) {
      rna_basal_concentration_array[i][rna_id] =
          rna_basal_concentration_list[i]->at(rna_id);
    }

    rna_produce_protein_array[i] = (int**)malloc(
        rna_produce_protein_list.size() * sizeof(int*));
    nb_rna_produce_protein[i] = (int*)malloc(
        rna_produce_protein_list.size() * sizeof(int));
    nb_rna_produce[i] = rna_produce_protein_list.size();

    for (int prot_id = 0; prot_id < rna_produce_protein_list.size(); prot_id++) {
      rna_produce_protein_array[i][prot_id] = (int*)malloc(
          rna_produce_protein_list[i]->at(prot_id)->size()*sizeof(int));
      nb_rna_produce_protein[i][prot_id] =
          rna_produce_protein_list[i]->at(prot_id)->size();
      for (int ix = 0; ix < rna_produce_protein_list[i]->at(prot_id)->size(); ix++) {
        rna_produce_protein_array[i][prot_id][ix] =
            rna_produce_protein_list[i]->at(prot_id)->at(ix);

      }
    }

    rna_influence_enhancing_coef_array[i] =
        (double**)malloc(rna_influence_enhancing_coef_list.size() * sizeof(double*));
    rna_influence_operating_coef_array[i] = (
        double**)malloc(rna_influence_operating_coef_list.size() * sizeof(double*));

    nb_rna_influence_enhancing_coef_l1[i] = rna_influence_enhancing_coef_list.size();
    nb_rna_influence_operating_coef_l1[i] = rna_influence_operating_coef_list.size();

    nb_rna_influence_enhancing_coef[i] =
        (int*)malloc(rna_influence_enhancing_coef_list.size() * sizeof(int));
    nb_rna_influence_operating_coef[i] =
        (int*)malloc(rna_influence_enhancing_coef_list.size() * sizeof(int));

    for (int rna_id = 0; rna_id < rna_influence_enhancing_coef_list.size(); rna_id++) {
      rna_influence_enhancing_coef_array[i][rna_id] =
          (double*)malloc(rna_influence_enhancing_coef_list[i]->at(rna_id)->size());
      nb_rna_influence_enhancing_coef[i][rna_id] =
          rna_influence_enhancing_coef_list[i]->at(rna_id)->size();

      for (int prot_id = 0; prot_id <
                            rna_influence_enhancing_coef_list[i]->at(rna_id)->size(); prot_id++) {
        rna_influence_enhancing_coef_array[i][rna_id][prot_id] =
            rna_influence_enhancing_coef_list[i]->at(rna_id)->at(prot_id);
      }
    }

    for (int rna_id = 0; rna_id < rna_influence_operating_coef_list.size(); rna_id++) {
      rna_influence_operating_coef_array[i][rna_id] =
          (double*)malloc(rna_influence_operating_coef_list[i]->at(rna_id)->size());
      nb_rna_influence_operating_coef[i][rna_id] =
          rna_influence_operating_coef_list[i]->at(rna_id)->size();

      for (int prot_id = 0; prot_id <
                            rna_influence_operating_coef_list[i]->at(rna_id)->size(); prot_id++) {
        rna_influence_operating_coef_array[i][rna_id][prot_id] =
            rna_influence_operating_coef_list[i]->at(rna_id)->at(prot_id);
      }
    }
  }

  env_concentration_array = (double**)malloc(lifestep * sizeof(double*));
  for (int i = 0; i < lifestep; i++) {
    env_concentration_array[i] = (double*)malloc(nb_signal * sizeof(double));
    for (int j=0; j < nb_signal; j++) {
      env_concentration_array[i][j] = env_concentration_list[i]->at(j);
    }
  }

  cudaMalloc((void**)&dev_protein_concentration_array, pop_size * sizeof(double *));
  cudaMalloc((void**)&dev_rna_basal_concentration_array, pop_size * sizeof(double *));

  cudaMalloc((void**)&dev_nb_protein_array, pop_size * sizeof(int));
  cudaMalloc((void**)&dev_nb_rna_array, pop_size * sizeof(int));

  cudaMemcpy(dev_nb_protein_array,
             nb_protein_array, pop_size * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_nb_rna_array,
             nb_rna_array, pop_size * sizeof(int), cudaMemcpyHostToDevice);

  cudaMalloc((void**)&dev_rna_produce_protein_array, pop_size * sizeof(int **));
  cudaMalloc((void**)&dev_nb_rna_produce_protein, pop_size * sizeof(int *));
  cudaMalloc((void**)&dev_nb_rna_produce, pop_size * sizeof(int));

  cudaMemcpy(dev_nb_rna_produce,
             nb_rna_produce, pop_size * sizeof(int), cudaMemcpyHostToDevice);

  cudaMalloc((void**)&dev_rna_influence_enhancing_coef_array,
             pop_size * sizeof(double **));
  cudaMalloc((void**)&dev_rna_influence_operating_coef_array,
             pop_size * sizeof(double **));

  cudaMalloc((void**)&dev_nb_rna_influence_enhancing_coef,
             pop_size * sizeof(int *));
  cudaMalloc((void**)&dev_nb_rna_influence_operating_coef,
             pop_size * sizeof(int *));
  cudaMalloc((void**)&dev_nb_rna_influence_enhancing_coef_l1,
             pop_size * sizeof(int));
  cudaMalloc((void**)&dev_nb_rna_influence_operating_coef_l1,
             pop_size * sizeof(int));

  cudaMemcpy(dev_nb_rna_influence_enhancing_coef_l1,
             nb_rna_influence_enhancing_coef_l1, pop_size * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_nb_rna_influence_operating_coef_l1,
             nb_rna_influence_operating_coef_l1, pop_size * sizeof(int), cudaMemcpyHostToDevice);


  for (int i = 0; i < pop_size; i++){
    cudaMalloc((void **)&dev_protein_concentration_array[i],
               nb_protein_array[i] * sizeof(double));
    cudaMemcpy(dev_protein_concentration_array[i],
               protein_concentration_array[i],
               nb_protein_array[i] * sizeof(double), cudaMemcpyHostToDevice);

    cudaMalloc((void **)&dev_rna_basal_concentration_array[i],
               rna_basal_concentration_list[i]->size() * sizeof(double));
    cudaMemcpy(dev_rna_basal_concentration_array[i],
               rna_basal_concentration_array[i],
               nb_rna_array[i] * sizeof(double), cudaMemcpyHostToDevice);

    cudaMalloc((void **)&dev_rna_produce_protein_array[i],
               nb_rna_produce[i] * sizeof(int*));
    cudaMalloc((void **)&dev_nb_rna_produce_protein[i],
               nb_rna_produce[i] * sizeof(int));
    cudaMemcpy(dev_nb_rna_produce_protein[i],
               nb_rna_produce_protein[i],
               nb_rna_produce[i] * sizeof(int), cudaMemcpyHostToDevice);

    for (int prot_id = 0; prot_id < nb_rna_produce[i]; prot_id++) {
      cudaMalloc((void **)&dev_rna_produce_protein_array[i][prot_id],
                 nb_rna_produce_protein[i][prot_id] * sizeof(int));
      cudaMemcpy(dev_rna_produce_protein_array[i][prot_id],
                 rna_produce_protein_array[i][prot_id],
                 nb_rna_produce_protein[i][prot_id] * sizeof(int),
                 cudaMemcpyHostToDevice);
    }


    cudaMalloc((void**)&dev_rna_influence_enhancing_coef_array[i],
               nb_rna_influence_enhancing_coef_l1[i] * sizeof(double *));
    cudaMalloc((void**)&dev_rna_influence_operating_coef_array[i],
               nb_rna_influence_operating_coef_l1[i] * sizeof(double *));

    cudaMalloc((void **)&dev_nb_rna_influence_enhancing_coef[i],
               nb_rna_influence_enhancing_coef_l1[i] * sizeof(int));
    cudaMalloc((void **)&dev_nb_rna_influence_operating_coef[i],
               nb_rna_influence_operating_coef_l1[i] * sizeof(int));

    cudaMemcpy(dev_nb_rna_influence_enhancing_coef[i],
               nb_rna_influence_enhancing_coef[i],
               nb_rna_influence_enhancing_coef_l1[i] * sizeof(int),
               cudaMemcpyHostToDevice);
    cudaMemcpy(dev_nb_rna_influence_operating_coef[i],
               nb_rna_influence_operating_coef[i],
               nb_rna_influence_operating_coef_l1[i] * sizeof(int),
               cudaMemcpyHostToDevice);

    for (int rna_id = 0; rna_id < nb_rna_influence_enhancing_coef_l1[i]; rna_id++) {
      cudaMalloc((void**)&dev_rna_influence_enhancing_coef_array[i][rna_id],
                 nb_rna_influence_enhancing_coef[i][rna_id] * sizeof(double));
      cudaMemcpy(dev_rna_influence_enhancing_coef_array[i][rna_id],
                 rna_influence_enhancing_coef_array[i][rna_id],
                 nb_rna_influence_enhancing_coef[i][rna_id] * sizeof(double),
                 cudaMemcpyHostToDevice);
    }


    for (int rna_id = 0; rna_id < nb_rna_influence_enhancing_coef_l1[i]; rna_id++) {
      cudaMalloc((void**)&dev_rna_influence_enhancing_coef_array[i][rna_id],
                 nb_rna_influence_enhancing_coef[i][rna_id] * sizeof(double));
      cudaMemcpy(dev_rna_influence_enhancing_coef_array[i][rna_id],
                 rna_influence_enhancing_coef_array[i][rna_id],
                 nb_rna_influence_enhancing_coef[i][rna_id] * sizeof(double),
                 cudaMemcpyHostToDevice);
    }

    for (int rna_id = 0; rna_id < nb_rna_influence_operating_coef_l1[i]; rna_id++) {
      cudaMalloc((void**)&dev_rna_influence_operating_coef_array[i][rna_id],
                 nb_rna_influence_operating_coef[i][rna_id] * sizeof(double));
      cudaMemcpy(dev_rna_influence_operating_coef_array[i][rna_id],
                 rna_influence_operating_coef_array[i][rna_id],
                 nb_rna_influence_operating_coef[i][rna_id] * sizeof(double),
                 cudaMemcpyHostToDevice);
    }

  }


  cudaMalloc((void**)&dev_env_concentration_array, lifestep * sizeof(double *));
  for (int i = 0; i < lifestep; i++) {
    cudaMalloc((void**)&dev_env_concentration_array[i], nb_signal * sizeof(double));
    cudaMemcpy(dev_env_concentration_array[i],
               env_concentration_array[i],
               nb_signal * sizeof(double),
               cudaMemcpyHostToDevice);
  }
}

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

      protein_concentration_array[indiv_id][prot_id] = delta;
    }
  }
}
