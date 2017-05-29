#include <vector>
#include <stdio.h>
#include <unistd.h>
#include<cuda_profiler_api.h>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include <thrust/copy.h>

#include "ode_solver_only_gpu.h"

void call_kernel_ode_cuda(int nb_gen, int multiply_population, int max_protein,
                          int nb_signal, int degradationstep,
                          int degradation_rate, double hill_shape_n, double hill_shape) {
  for (int gen = 0; gen < nb_gen; gen++) {
    process_delta << < 1024 * multiply_population, max_protein >> >
                                                   (nb_signal, degradationstep, degradation_rate,
                                                       dev_rna_produce_protein_array, dev_nb_rna_produce_protein, dev_nb_rna_produce, dev_protein_concentration_array,
                                                       dev_rna_basal_concentration_array, dev_nb_protein_array, dev_nb_rna_array,
                                                       dev_rna_influence_enhancing_coef_array, dev_rna_influence_operating_coef_array,
                                                       dev_nb_rna_influence_enhancing_coef, dev_nb_rna_influence_operating_coef,
                                                       dev_env_concentration_array, hill_shape, hill_shape_n);
  }

  cudaDeviceSynchronize();
  cudaProfilerStop();
}


void call_kernel_ode_cuda_float(int nb_gen, int multiply_population, int max_protein,
                          int nb_signal, int degradationstep,
                          int degradation_rate, double hill_shape_n, double hill_shape) {
  for (int gen = 0; gen < nb_gen; gen++) {
    process_delta_float << < 1024 * multiply_population, max_protein >> >
                                                         (nb_signal, degradationstep, degradation_rate,
                                                             dev_rna_produce_protein_array, dev_nb_rna_produce_protein, dev_nb_rna_produce, f_dev_protein_concentration_array,
                                                             f_dev_rna_basal_concentration_array, dev_nb_protein_array, dev_nb_rna_array,
                                                             f_dev_rna_influence_enhancing_coef_array, f_dev_rna_influence_operating_coef_array,
                                                             dev_nb_rna_influence_enhancing_coef, dev_nb_rna_influence_operating_coef,
                                                             f_dev_env_concentration_array, (float) hill_shape, (float) hill_shape_n);
  }

  cudaDeviceSynchronize();
  cudaProfilerStop();
}

void call_kernel_ode_cuda_thrust(int multiply_population,
                                int nb_signal, int degradationstep,
                                int degradation_rate, double hill_shape_n, double hill_shape) {





  /*process_delta_thrust<<<1024*multiply_population,g_max_protein>>>(nb_signal,degradationstep,degradation_rate,
      g_max_protein,g_max_rna,
      r_gpu_thrust_nb_rna_produce_protein, r_gpu_thrust_rna_produce_protein,
      r_gpu_thrust_protein_concentration,
      r_gpu_thrust_rna_basal_concentration, r_gpu_thrust_nb_protein,
      r_gpu_thrust_rna_influence_enhancing_coef, r_gpu_thrust_rna_influence_operating_coef,
      r_gpu_thrust_nb_influence,
      r_gpu_thrust_environment_concentration, hill_shape, hill_shape_n);*/
}

int transfert_data_to_gpu(int pop_size, int lifestep,
         std::vector<std::vector<double>*> const &protein_concentration_list,
         std::vector<std::vector<double>*> const &rna_basal_concentration_list,
         std::vector<std::vector<std::vector<int>*>*> const &rna_produce_protein_list,
         std::vector<std::vector<std::vector<double>*>*> const &rna_influence_enhancing_coef_list,
         std::vector<std::vector<std::vector<double>*>*> const &rna_influence_operating_coef_list,
         int nb_signal,
         std::vector<std::vector<double>*> const &env_concentration_list) {
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

  int max_prot = 0, min_prot = 100000000;

  for (int i = 0; i < pop_size; i++){
    protein_concentration_array[i] = (double*)
        malloc(protein_concentration_list[i]->size() * sizeof(double));
    nb_protein_array[i] = (int) protein_concentration_list[i]->size();

    max_prot = protein_concentration_list[i]->size() > max_prot ?
               protein_concentration_list[i]->size() : max_prot;

    min_prot = protein_concentration_list[i]->size() < min_prot ?
               protein_concentration_list[i]->size() : min_prot;



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
        rna_produce_protein_list[i]->size() * sizeof(int*));
    nb_rna_produce_protein[i] = (int*)malloc(
        rna_produce_protein_list[i]->size() * sizeof(int));
    nb_rna_produce[i] = rna_produce_protein_list[i]->size();

    for (int prot_id = 0; prot_id < rna_produce_protein_list[i]->size(); prot_id++) {
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
        (double**)malloc(rna_influence_enhancing_coef_list[i]->size() * sizeof(double*));
    rna_influence_operating_coef_array[i] = (
        double**)malloc(rna_influence_operating_coef_list[i]->size() * sizeof(double*));

    nb_rna_influence_enhancing_coef_l1[i] = rna_influence_enhancing_coef_list[i]->size();
    nb_rna_influence_operating_coef_l1[i] = rna_influence_operating_coef_list[i]->size();

    nb_rna_influence_enhancing_coef[i] =
        (int*)malloc(rna_influence_enhancing_coef_list[i]->size() * sizeof(int));
    nb_rna_influence_operating_coef[i] =
        (int*)malloc(rna_influence_enhancing_coef_list[i]->size() * sizeof(int));

    for (int rna_id = 0; rna_id < rna_influence_enhancing_coef_list[i]->size(); rna_id++) {

      rna_influence_enhancing_coef_array[i][rna_id] =
          (double*)malloc(rna_influence_enhancing_coef_list[i]->at(rna_id)->size() * sizeof(double));

      nb_rna_influence_enhancing_coef[i][rna_id] =
          rna_influence_enhancing_coef_list[i]->at(rna_id)->size();

      for (int prot_id = 0; prot_id <
                            rna_influence_enhancing_coef_list[i]->at(rna_id)->size(); prot_id++) {
        rna_influence_enhancing_coef_array[i][rna_id][prot_id] =
            rna_influence_enhancing_coef_list[i]->at(rna_id)->at(prot_id);
      }
    }

    for (int rna_id = 0; rna_id < rna_influence_operating_coef_list[i]->size(); rna_id++) {
      rna_influence_operating_coef_array[i][rna_id] =
          (double*)malloc(rna_influence_operating_coef_list[i]->at(rna_id)->size()* sizeof(double));
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

  ///// COPY TO GPU
  // env

  double **e_d = (double**)malloc(lifestep*sizeof(double*));

  for(int i=0; i<lifestep; i++) {
    cudaMalloc((void**) &e_d[i],
               nb_signal * sizeof(double));
    cudaMemcpy(e_d[i], env_concentration_array[i],
               nb_signal*sizeof(double),
               cudaMemcpyHostToDevice);
  }

  cudaMalloc((void***)&dev_env_concentration_array,lifestep*sizeof(double*));
  cudaMemcpy(dev_env_concentration_array,e_d,lifestep*sizeof(double*),cudaMemcpyHostToDevice);

  // rna_influence_enhancing_coef_array
  double ***h_c = (double***)malloc(pop_size*sizeof(double**));

  for(int i=0; i<pop_size; i++) {
    h_c[i] = (double**) malloc(nb_rna_influence_enhancing_coef_l1[i]*sizeof(double*));
    for(int j=0; j<nb_rna_influence_enhancing_coef_l1[i]; j++) {
      cudaMalloc((void**) &h_c[i][j],
                 nb_rna_influence_enhancing_coef[i][j] * sizeof(double));
      cudaMemcpy(h_c[i][j], rna_influence_enhancing_coef_array[i][j],
                 nb_rna_influence_enhancing_coef[i][j]*sizeof(double),
                 cudaMemcpyHostToDevice);
    }

  }

  double ***h_c1 = (double ***) malloc(pop_size*sizeof(double **));
  for (int i=0; i<pop_size; i++){
    cudaMalloc((void***)&(h_c1[i]), pop_size*sizeof(double*));
    cudaMemcpy(h_c1[i], h_c[i], pop_size*sizeof(double*), cudaMemcpyHostToDevice);
  }

  cudaMalloc((void****)&dev_rna_influence_enhancing_coef_array,pop_size*sizeof(double**));
  cudaMemcpy(dev_rna_influence_enhancing_coef_array,h_c1,pop_size*sizeof(double**),cudaMemcpyHostToDevice);

  // rna_influence_operating_coef_array

  h_c = (double***)malloc(pop_size*sizeof(double**));

  for(int i=0; i<pop_size; i++) {
    h_c[i] = (double**) malloc(nb_rna_influence_operating_coef_l1[i]*sizeof(double*));
    for(int j=0; j<nb_rna_influence_operating_coef_l1[i]; j++) {
      cudaMalloc((void**) &h_c[i][j],
                 nb_rna_influence_operating_coef[i][j] * sizeof(double));
      cudaMemcpy(h_c[i][j], rna_influence_operating_coef_array[i][j],
                 nb_rna_influence_operating_coef[i][j]*sizeof(double),
                 cudaMemcpyHostToDevice);
    }

  }

  h_c1 = (double ***) malloc(pop_size*sizeof(double **));
  for (int i=0; i<pop_size; i++){
    cudaMalloc((void***)&(h_c1[i]), pop_size*sizeof(double*));
    cudaMemcpy(h_c1[i], h_c[i], pop_size*sizeof(double*), cudaMemcpyHostToDevice);
  }

  cudaMalloc((void****)&dev_rna_influence_operating_coef_array,pop_size*sizeof(double**));
  cudaMemcpy(dev_rna_influence_operating_coef_array,h_c1,pop_size*sizeof(double**),cudaMemcpyHostToDevice);
  
  // protein_concentration_array (double**)malloc(pop_size * sizeof(double*));
  double **h_d = (double**)malloc(pop_size*sizeof(double*));

  for(int i=0; i<pop_size; i++) {
    cudaMalloc((void**) &h_d[i],
                 nb_protein_array[i] * sizeof(double));
    cudaMemcpy(h_d[i], protein_concentration_array[i],
                 nb_protein_array[i]*sizeof(double),
                 cudaMemcpyHostToDevice);
  }

  cudaMalloc((void***)&dev_protein_concentration_array,pop_size*sizeof(double*));
  cudaMemcpy(dev_protein_concentration_array,h_d,pop_size*sizeof(double*),cudaMemcpyHostToDevice);

  // rna_basal_concentration_array
  h_d = (double**)malloc(pop_size*sizeof(double*));

  for(int i=0; i<pop_size; i++) {
    cudaMalloc((void**) &h_d[i],
               nb_rna_array[i] * sizeof(double));
    cudaMemcpy(h_d[i], rna_basal_concentration_array[i],
               nb_rna_array[i]*sizeof(double),
               cudaMemcpyHostToDevice);
  }

  cudaMalloc((void***)&dev_rna_basal_concentration_array,pop_size*sizeof(double*));
  cudaMemcpy(dev_rna_basal_concentration_array,h_d,pop_size*sizeof(double*),cudaMemcpyHostToDevice);

  // rna_produce_protein_array

  int ***i_c = (int***)malloc(pop_size*sizeof(int**));

  for(int i=0; i<pop_size; i++) {
    i_c[i] = (int**) malloc(nb_rna_produce[i]*sizeof(int*));
    for(int j=0; j<nb_rna_produce[i]; j++) {
      cudaMalloc((void**) &i_c[i][j],
                 nb_rna_produce_protein[i][j] * sizeof(int));
      cudaMemcpy(i_c[i][j], rna_produce_protein_array[i][j],
                 nb_rna_produce_protein[i][j]*sizeof(int),
                 cudaMemcpyHostToDevice);
    }

  }

  int ***i_c1 = (int ***) malloc(pop_size*sizeof(int **));
  for (int i=0; i<pop_size; i++){
    cudaMalloc((void***)&(i_c1[i]), pop_size*sizeof(int*));
    cudaMemcpy(i_c1[i], i_c[i], pop_size*sizeof(int*), cudaMemcpyHostToDevice);
  }

  cudaMalloc((void****)&dev_rna_produce_protein_array,pop_size*sizeof(int**));
  cudaMemcpy(dev_rna_produce_protein_array,i_c1,pop_size*sizeof(int**),cudaMemcpyHostToDevice);

  // nb_protein_array
  cudaMalloc((void**)&dev_nb_protein_array, pop_size * sizeof(int));
  cudaMemcpy(dev_nb_protein_array,
             nb_protein_array, pop_size * sizeof(int), cudaMemcpyHostToDevice);

  // nb_rna_array
  cudaMalloc((void**)&dev_nb_rna_array, pop_size * sizeof(int));
  cudaMemcpy(dev_nb_rna_array,
             nb_rna_array, pop_size * sizeof(int), cudaMemcpyHostToDevice);

  // nb_rna_produce_protein
  int **i_d = (int**)malloc(pop_size*sizeof(int*));

  for(int i=0; i<pop_size; i++) {
    cudaMalloc((void**) &i_d[i],
               nb_rna_produce[i] * sizeof(int));
    cudaMemcpy(i_d[i], nb_rna_produce_protein[i],
               nb_rna_produce[i]*sizeof(int),
               cudaMemcpyHostToDevice);
  }

  cudaMalloc((void***)&dev_nb_rna_produce_protein,pop_size*sizeof(int*));
  cudaMemcpy(dev_nb_rna_produce_protein,i_d,pop_size*sizeof(int*),cudaMemcpyHostToDevice);

  // nb_rna_produce
  cudaMalloc((void**)&dev_nb_rna_produce, pop_size * sizeof(int));

  cudaMemcpy(dev_nb_rna_produce,
             nb_rna_produce, pop_size * sizeof(int), cudaMemcpyHostToDevice);

  // nb_rna_influence_enhancing_coef

  i_d = (int**)malloc(pop_size*sizeof(int*));

  for(int i=0; i<pop_size; i++) {
    cudaMalloc((void**) &i_d[i],
               nb_rna_influence_enhancing_coef_l1[i] * sizeof(int));
    cudaMemcpy(i_d[i], nb_rna_influence_enhancing_coef[i],
               nb_rna_influence_enhancing_coef_l1[i]*sizeof(int),
               cudaMemcpyHostToDevice);
  }

  cudaMalloc((void***)&dev_nb_rna_influence_enhancing_coef,pop_size*sizeof(int*));
  cudaMemcpy(dev_nb_rna_influence_enhancing_coef,i_d,pop_size*sizeof(int*),cudaMemcpyHostToDevice);

  // nb_rna_influence_operating_coef

  i_d = (int**)malloc(pop_size*sizeof(int*));

  for(int i=0; i<pop_size; i++) {
    cudaMalloc((void**) &i_d[i],
               nb_rna_influence_operating_coef_l1[i] * sizeof(int));
    cudaMemcpy(i_d[i], nb_rna_influence_operating_coef[i],
               nb_rna_influence_operating_coef_l1[i]*sizeof(int),
               cudaMemcpyHostToDevice);
  }

  cudaMalloc((void***)&dev_nb_rna_influence_operating_coef,pop_size*sizeof(int*));
  cudaMemcpy(dev_nb_rna_influence_operating_coef,i_d,pop_size*sizeof(int*),cudaMemcpyHostToDevice);

  // nb_rna_influence_enhancing_coef_l1
  cudaMalloc((void**)&dev_nb_rna_influence_enhancing_coef_l1,
             pop_size * sizeof(int));
  cudaMemcpy(dev_nb_rna_influence_enhancing_coef_l1,
             nb_rna_influence_enhancing_coef_l1, pop_size * sizeof(int), cudaMemcpyHostToDevice);

  // nb_rna_influence_operating_coef_l1
  cudaMalloc((void**)&dev_nb_rna_influence_operating_coef_l1,
             pop_size * sizeof(int));
  cudaMemcpy(dev_nb_rna_influence_operating_coef_l1,
             nb_rna_influence_operating_coef_l1, pop_size * sizeof(int), cudaMemcpyHostToDevice);
  ////

  return max_prot;
}



int transfert_data_to_gpu_float(int pop_size, int lifestep,
                          std::vector<std::vector<double>*> const &protein_concentration_list,
                          std::vector<std::vector<double>*> const &rna_basal_concentration_list,
                          std::vector<std::vector<std::vector<int>*>*> const &rna_produce_protein_list,
                          std::vector<std::vector<std::vector<double>*>*> const &rna_influence_enhancing_coef_list,
                          std::vector<std::vector<std::vector<double>*>*> const &rna_influence_operating_coef_list,
                          int nb_signal,
                          std::vector<std::vector<double>*> const &env_concentration_list) {
  f_protein_concentration_array = (float**)malloc(pop_size * sizeof(float*));
  f_rna_basal_concentration_array = (float**)malloc(pop_size * sizeof(float*));
  rna_produce_protein_array = (int***)malloc(pop_size * sizeof(int**));
  f_rna_influence_enhancing_coef_array = (float***)malloc(pop_size * sizeof(float**));
  f_rna_influence_operating_coef_array = (float***)malloc(pop_size * sizeof(float**));

  nb_protein_array = (int*)malloc(pop_size * sizeof(int));
  nb_rna_array = (int*)malloc(pop_size * sizeof(int));
  nb_rna_produce_protein = (int**)malloc(pop_size * sizeof(int*));
  nb_rna_produce = (int*)malloc(pop_size * sizeof(int));
  nb_rna_influence_enhancing_coef = (int**) malloc(pop_size * sizeof(int*));
  nb_rna_influence_operating_coef = (int**) malloc(pop_size * sizeof(int*));
  nb_rna_influence_enhancing_coef_l1 = (int*) malloc(pop_size * sizeof(int));
  nb_rna_influence_operating_coef_l1 = (int*) malloc(pop_size * sizeof(int));

  int max_prot = 0, min_prot = 100000000;

  for (int i = 0; i < pop_size; i++){
    f_protein_concentration_array[i] = (float*)
        malloc(protein_concentration_list[i]->size() * sizeof(float));
    nb_protein_array[i] = (int) protein_concentration_list[i]->size();

    max_prot = protein_concentration_list[i]->size() > max_prot ?
               protein_concentration_list[i]->size() : max_prot;

    min_prot = protein_concentration_list[i]->size() < min_prot ?
               protein_concentration_list[i]->size() : min_prot;



    for (int prot_id = 0; prot_id < protein_concentration_list[i]->size(); prot_id++) {
      f_protein_concentration_array[i][prot_id] = (float)
          protein_concentration_list[i]->at(prot_id);
    }

    f_rna_basal_concentration_array[i] = (float*)
        malloc(rna_basal_concentration_list[i]->size() * sizeof(float));
    nb_rna_array[i] = (int)rna_basal_concentration_list[i]->size();

    for (int rna_id = 0; rna_id < rna_basal_concentration_list[i]->size(); rna_id++) {
      f_rna_basal_concentration_array[i][rna_id] = (float)
          rna_basal_concentration_list[i]->at(rna_id);
    }

    rna_produce_protein_array[i] = (int**)malloc(
        rna_produce_protein_list[i]->size() * sizeof(int*));
    nb_rna_produce_protein[i] = (int*)malloc(
        rna_produce_protein_list[i]->size() * sizeof(int));
    nb_rna_produce[i] = rna_produce_protein_list[i]->size();

    for (int prot_id = 0; prot_id < rna_produce_protein_list[i]->size(); prot_id++) {
      rna_produce_protein_array[i][prot_id] = (int*)malloc(
          rna_produce_protein_list[i]->at(prot_id)->size()*sizeof(int));
      nb_rna_produce_protein[i][prot_id] =
          rna_produce_protein_list[i]->at(prot_id)->size();
      for (int ix = 0; ix < rna_produce_protein_list[i]->at(prot_id)->size(); ix++) {
        rna_produce_protein_array[i][prot_id][ix] =
            rna_produce_protein_list[i]->at(prot_id)->at(ix);

      }
    }

    f_rna_influence_enhancing_coef_array[i] =
        (float**)malloc(rna_influence_enhancing_coef_list[i]->size() * sizeof(float*));
    f_rna_influence_operating_coef_array[i] = (
        float**)malloc(rna_influence_operating_coef_list[i]->size() * sizeof(float*));

    nb_rna_influence_enhancing_coef_l1[i] = rna_influence_enhancing_coef_list[i]->size();
    nb_rna_influence_operating_coef_l1[i] = rna_influence_operating_coef_list[i]->size();

    nb_rna_influence_enhancing_coef[i] =
        (int*)malloc(rna_influence_enhancing_coef_list[i]->size() * sizeof(int));
    nb_rna_influence_operating_coef[i] =
        (int*)malloc(rna_influence_enhancing_coef_list[i]->size() * sizeof(int));

    for (int rna_id = 0; rna_id < rna_influence_enhancing_coef_list[i]->size(); rna_id++) {

      f_rna_influence_enhancing_coef_array[i][rna_id] =
          (float*)malloc(rna_influence_enhancing_coef_list[i]->at(rna_id)->size() * sizeof(float));

      nb_rna_influence_enhancing_coef[i][rna_id] =
          rna_influence_enhancing_coef_list[i]->at(rna_id)->size();

      for (int prot_id = 0; prot_id <
                            rna_influence_enhancing_coef_list[i]->at(rna_id)->size(); prot_id++) {
        f_rna_influence_enhancing_coef_array[i][rna_id][prot_id] = (float)
            rna_influence_enhancing_coef_list[i]->at(rna_id)->at(prot_id);
      }
    }

    for (int rna_id = 0; rna_id < rna_influence_operating_coef_list[i]->size(); rna_id++) {
      f_rna_influence_operating_coef_array[i][rna_id] =
          (float*)malloc(rna_influence_operating_coef_list[i]->at(rna_id)->size()* sizeof(float));
      nb_rna_influence_operating_coef[i][rna_id] =
          rna_influence_operating_coef_list[i]->at(rna_id)->size();

      for (int prot_id = 0; prot_id <
                            rna_influence_operating_coef_list[i]->at(rna_id)->size(); prot_id++) {
        f_rna_influence_operating_coef_array[i][rna_id][prot_id] = (float)
            rna_influence_operating_coef_list[i]->at(rna_id)->at(prot_id);
      }
    }
  }

  f_env_concentration_array = (float**)malloc(lifestep * sizeof(float*));
  for (int i = 0; i < lifestep; i++) {
    f_env_concentration_array[i] = (float*)malloc(nb_signal * sizeof(float));
    for (int j=0; j < nb_signal; j++) {
      f_env_concentration_array[i][j] = (float) env_concentration_list[i]->at(j);
    }
  }

  ///// COPY TO GPU
  // env

  float **e_d = (float**)malloc(lifestep*sizeof(float*));

  for(int i=0; i<lifestep; i++) {
    cudaMalloc((void**) &e_d[i],
               nb_signal * sizeof(float));
    cudaMemcpy(e_d[i], f_env_concentration_array[i],
               nb_signal*sizeof(float),
               cudaMemcpyHostToDevice);
  }

  cudaMalloc((void***)&f_dev_env_concentration_array,lifestep*sizeof(float*));
  cudaMemcpy(f_dev_env_concentration_array,e_d,lifestep*sizeof(float*),cudaMemcpyHostToDevice);

  // rna_influence_enhancing_coef_array
  float ***h_c = (float***)malloc(pop_size*sizeof(float**));

  for(int i=0; i<pop_size; i++) {
    h_c[i] = (float**) malloc(nb_rna_influence_enhancing_coef_l1[i]*sizeof(float*));
    for(int j=0; j<nb_rna_influence_enhancing_coef_l1[i]; j++) {
      cudaMalloc((void**) &h_c[i][j],
                 nb_rna_influence_enhancing_coef[i][j] * sizeof(float));
      cudaMemcpy(h_c[i][j], f_rna_influence_enhancing_coef_array[i][j],
                 nb_rna_influence_enhancing_coef[i][j]*sizeof(float),
                 cudaMemcpyHostToDevice);
    }

  }

  float ***h_c1 = (float ***) malloc(pop_size*sizeof(float **));
  for (int i=0; i<pop_size; i++){
    cudaMalloc((void***)&(h_c1[i]), pop_size*sizeof(float*));
    cudaMemcpy(h_c1[i], h_c[i], pop_size*sizeof(float*), cudaMemcpyHostToDevice);
  }

  cudaMalloc((void****)&f_dev_rna_influence_enhancing_coef_array,pop_size*sizeof(float**));
  cudaMemcpy(f_dev_rna_influence_enhancing_coef_array,h_c1,pop_size*sizeof(float**),cudaMemcpyHostToDevice);

  // rna_influence_operating_coef_array

  h_c = (float***)malloc(pop_size*sizeof(float**));

  for(int i=0; i<pop_size; i++) {
    h_c[i] = (float**) malloc(nb_rna_influence_operating_coef_l1[i]*sizeof(float*));
    for(int j=0; j<nb_rna_influence_operating_coef_l1[i]; j++) {
      cudaMalloc((void**) &h_c[i][j],
                 nb_rna_influence_operating_coef[i][j] * sizeof(float));
      cudaMemcpy(h_c[i][j], f_rna_influence_operating_coef_array[i][j],
                 nb_rna_influence_operating_coef[i][j]*sizeof(float),
                 cudaMemcpyHostToDevice);
    }

  }

  h_c1 = (float ***) malloc(pop_size*sizeof(float **));
  for (int i=0; i<pop_size; i++){
    cudaMalloc((void***)&(h_c1[i]), pop_size*sizeof(float*));
    cudaMemcpy(h_c1[i], h_c[i], pop_size*sizeof(float*), cudaMemcpyHostToDevice);
  }

  cudaMalloc((void****)&f_dev_rna_influence_operating_coef_array,pop_size*sizeof(float**));
  cudaMemcpy(f_dev_rna_influence_operating_coef_array,h_c1,pop_size*sizeof(float**),cudaMemcpyHostToDevice);

  // protein_concentration_array (double**)malloc(pop_size * sizeof(double*));
  float **h_d = (float**)malloc(pop_size*sizeof(float*));

  for(int i=0; i<pop_size; i++) {
    cudaMalloc((void**) &h_d[i],
               nb_protein_array[i] * sizeof(float));
    cudaMemcpy(h_d[i], f_protein_concentration_array[i],
               nb_protein_array[i]*sizeof(float),
               cudaMemcpyHostToDevice);
  }

  cudaMalloc((void***)&f_dev_protein_concentration_array,pop_size*sizeof(float*));
  cudaMemcpy(f_dev_protein_concentration_array,h_d,pop_size*sizeof(float*),cudaMemcpyHostToDevice);

  // rna_basal_concentration_array
  h_d = (float**)malloc(pop_size*sizeof(float*));

  for(int i=0; i<pop_size; i++) {
    cudaMalloc((void**) &h_d[i],
               nb_rna_array[i] * sizeof(float));
    cudaMemcpy(h_d[i], f_rna_basal_concentration_array[i],
               nb_rna_array[i]*sizeof(float),
               cudaMemcpyHostToDevice);
  }

  cudaMalloc((void***)&f_dev_rna_basal_concentration_array,pop_size*sizeof(float*));
  cudaMemcpy(f_dev_rna_basal_concentration_array,h_d,pop_size*sizeof(float*),cudaMemcpyHostToDevice);

  // rna_produce_protein_array

  int ***i_c = (int***)malloc(pop_size*sizeof(int**));

  for(int i=0; i<pop_size; i++) {
    i_c[i] = (int**) malloc(nb_rna_produce[i]*sizeof(int*));
    for(int j=0; j<nb_rna_produce[i]; j++) {
      cudaMalloc((void**) &i_c[i][j],
                 nb_rna_produce_protein[i][j] * sizeof(int));
      cudaMemcpy(i_c[i][j], rna_produce_protein_array[i][j],
                 nb_rna_produce_protein[i][j]*sizeof(int),
                 cudaMemcpyHostToDevice);
    }

  }

  int ***i_c1 = (int ***) malloc(pop_size*sizeof(int **));
  for (int i=0; i<pop_size; i++){
    cudaMalloc((void***)&(i_c1[i]), pop_size*sizeof(int*));
    cudaMemcpy(i_c1[i], i_c[i], pop_size*sizeof(int*), cudaMemcpyHostToDevice);
  }

  cudaMalloc((void****)&dev_rna_produce_protein_array,pop_size*sizeof(int**));
  cudaMemcpy(dev_rna_produce_protein_array,i_c1,pop_size*sizeof(int**),cudaMemcpyHostToDevice);

  // nb_protein_array
  cudaMalloc((void**)&dev_nb_protein_array, pop_size * sizeof(int));
  cudaMemcpy(dev_nb_protein_array,
             nb_protein_array, pop_size * sizeof(int), cudaMemcpyHostToDevice);

  // nb_rna_array
  cudaMalloc((void**)&dev_nb_rna_array, pop_size * sizeof(int));
  cudaMemcpy(dev_nb_rna_array,
             nb_rna_array, pop_size * sizeof(int), cudaMemcpyHostToDevice);

  // nb_rna_produce_protein
  int **i_d = (int**)malloc(pop_size*sizeof(int*));

  for(int i=0; i<pop_size; i++) {
    cudaMalloc((void**) &i_d[i],
               nb_rna_produce[i] * sizeof(int));
    cudaMemcpy(i_d[i], nb_rna_produce_protein[i],
               nb_rna_produce[i]*sizeof(int),
               cudaMemcpyHostToDevice);
  }

  cudaMalloc((void***)&dev_nb_rna_produce_protein,pop_size*sizeof(int*));
  cudaMemcpy(dev_nb_rna_produce_protein,i_d,pop_size*sizeof(int*),cudaMemcpyHostToDevice);

  // nb_rna_produce
  cudaMalloc((void**)&dev_nb_rna_produce, pop_size * sizeof(int));

  cudaMemcpy(dev_nb_rna_produce,
             nb_rna_produce, pop_size * sizeof(int), cudaMemcpyHostToDevice);

  // nb_rna_influence_enhancing_coef

  i_d = (int**)malloc(pop_size*sizeof(int*));

  for(int i=0; i<pop_size; i++) {
    cudaMalloc((void**) &i_d[i],
               nb_rna_influence_enhancing_coef_l1[i] * sizeof(int));
    cudaMemcpy(i_d[i], nb_rna_influence_enhancing_coef[i],
               nb_rna_influence_enhancing_coef_l1[i]*sizeof(int),
               cudaMemcpyHostToDevice);
  }

  cudaMalloc((void***)&dev_nb_rna_influence_enhancing_coef,pop_size*sizeof(int*));
  cudaMemcpy(dev_nb_rna_influence_enhancing_coef,i_d,pop_size*sizeof(int*),cudaMemcpyHostToDevice);

  // nb_rna_influence_operating_coef

  i_d = (int**)malloc(pop_size*sizeof(int*));

  for(int i=0; i<pop_size; i++) {
    cudaMalloc((void**) &i_d[i],
               nb_rna_influence_operating_coef_l1[i] * sizeof(int));
    cudaMemcpy(i_d[i], nb_rna_influence_operating_coef[i],
               nb_rna_influence_operating_coef_l1[i]*sizeof(int),
               cudaMemcpyHostToDevice);
  }

  cudaMalloc((void***)&dev_nb_rna_influence_operating_coef,pop_size*sizeof(int*));
  cudaMemcpy(dev_nb_rna_influence_operating_coef,i_d,pop_size*sizeof(int*),cudaMemcpyHostToDevice);

  // nb_rna_influence_enhancing_coef_l1
  cudaMalloc((void**)&dev_nb_rna_influence_enhancing_coef_l1,
             pop_size * sizeof(int));
  cudaMemcpy(dev_nb_rna_influence_enhancing_coef_l1,
             nb_rna_influence_enhancing_coef_l1, pop_size * sizeof(int), cudaMemcpyHostToDevice);

  // nb_rna_influence_operating_coef_l1
  cudaMalloc((void**)&dev_nb_rna_influence_operating_coef_l1,
             pop_size * sizeof(int));
  cudaMemcpy(dev_nb_rna_influence_operating_coef_l1,
             nb_rna_influence_operating_coef_l1, pop_size * sizeof(int), cudaMemcpyHostToDevice);
  ////

  return max_prot;
}


int transfert_data_to_gpu_dense(int pop_size, int lifestep,
                          std::vector<std::vector<double>*> const &protein_concentration_list,
                          std::vector<std::vector<double>*> const &rna_basal_concentration_list,
                          std::vector<std::vector<std::vector<int>*>*> const &rna_produce_protein_list,
                          std::vector<std::vector<std::vector<double>*>*> const &rna_influence_enhancing_coef_list,
                          std::vector<std::vector<std::vector<double>*>*> const &rna_influence_operating_coef_list,
                          int nb_signal,
                          std::vector<std::vector<double>*> const &env_concentration_list) {

  int max_prot = 0, max_rna = 0;

  /* Computing max rna and protein */
  for (int i = 0; i < pop_size; i++) {
    max_prot = protein_concentration_list[i]->size() > max_prot ?
               protein_concentration_list[i]->size() : max_prot;
    max_rna  = rna_basal_concentration_list[i]->size() > max_rna ?
               rna_basal_concentration_list[i]->size() : max_rna;
  }

  /* Storing as dense tab */

  protein_tab = (double*)malloc(pop_size * max_prot * sizeof(double));
  memset (protein_tab, -1, pop_size * max_prot * sizeof(double));


  rna_basal_tab = (double*)malloc(pop_size * max_rna * sizeof(double));
  memset (rna_basal_tab, -1, pop_size * max_rna * sizeof(double));

  rna_produce_protein_tab = (int*)malloc(pop_size * max_prot * max_rna * sizeof(int));
  memset (rna_produce_protein_tab, -1, pop_size * max_prot * max_rna * sizeof(int));


  rna_influence_enhancing_tab = (double*)malloc(pop_size * max_prot * max_rna * sizeof(double));
  memset (rna_influence_enhancing_tab, 0, pop_size * max_prot * max_rna * sizeof(double));

  rna_influence_operating_tab = (double*)malloc(pop_size * max_prot * max_rna * sizeof(double));
  memset (rna_influence_operating_tab, 0, pop_size * max_prot * max_rna * sizeof(double));

  for (int i = 0; i < pop_size; i++) {
    for (int prot_id = 0; prot_id < protein_concentration_list[i]->size(); prot_id++) {
      protein_tab[i*max_prot+prot_id] =
          protein_concentration_list[i]->at(prot_id);
    }

    for (int rna_id = 0; rna_id < rna_basal_concentration_list[i]->size(); rna_id++) {
      rna_basal_tab[i*max_rna+rna_id] =
          rna_basal_concentration_list[i]->at(rna_id);
    }

    for (int prot_id = 0; prot_id < rna_produce_protein_list[i]->size(); prot_id++) {
      for (int ix = 0; ix < rna_produce_protein_list[i]->at(prot_id)->size(); ix++) {
        rna_produce_protein_tab[i*max_prot+prot_id*max_rna+ix] =
            rna_produce_protein_list[i]->at(prot_id)->at(ix);

      }
    }

    for (int rna_id = 0; rna_id < rna_influence_enhancing_coef_list[i]->size(); rna_id++) {
      for (int prot_id = 0; prot_id <
                            rna_influence_enhancing_coef_list[i]->at(rna_id)->size(); prot_id++) {
        rna_influence_enhancing_tab[i*max_rna+rna_id*max_prot+prot_id] =
            rna_influence_enhancing_coef_list[i]->at(rna_id)->at(prot_id);
      }
    }

    for (int rna_id = 0; rna_id < rna_influence_operating_coef_list[i]->size(); rna_id++) {
      for (int prot_id = 0; prot_id <
                            rna_influence_operating_coef_list[i]->at(rna_id)->size(); prot_id++) {
        rna_influence_operating_tab[i*max_rna+rna_id*max_prot+prot_id] =
            rna_influence_operating_coef_list[i]->at(rna_id)->at(prot_id);
      }
    }
  }

  env_concentration_tab = (double*)malloc(lifestep * nb_signal * sizeof(double));
  for (int i = 0; i < lifestep; i++) {
    for (int j=0; j < nb_signal; j++) {
      env_concentration_tab[i*nb_signal+j] = env_concentration_list[i]->at(j);
    }
  }

  ///// COPY TO GPU
  cudaMallocPitch(&dev_protein_tab, &pitch_protein_tab, sizeof(double)*max_prot, pop_size);
  cudaMemcpy2D(dev_protein_tab,pitch_protein_tab,protein_tab,
               sizeof(double)*max_prot,sizeof(double)*max_prot,pop_size,
               cudaMemcpyHostToDevice);

  cudaMallocPitch(&dev_rna_basal_tab, &pitch_rna_basal_tab, sizeof(double)*max_rna, pop_size);
  cudaMemcpy2D(dev_rna_basal_tab,pitch_rna_basal_tab,rna_basal_tab,
               sizeof(double)*max_rna,sizeof(double)*max_rna,pop_size,
               cudaMemcpyHostToDevice);

  // Allocate 3D memory on the device
  cudaExtent volumeSizeBytes_rna_produce =
      make_cudaExtent(sizeof(int) * max_rna, max_prot, pop_size);

  cudaPitchedPtr devicePitchedPointer_rna_produce;

  cudaMalloc3D(&devicePitchedPointer_rna_produce, volumeSizeBytes_rna_produce);

  cudaMemcpy3DParms p_rna_produce = { 0 };

  p_rna_produce.srcPtr.ptr = rna_produce_protein_tab;
  p_rna_produce.srcPtr.pitch = max_rna * sizeof(int);
  p_rna_produce.srcPtr.xsize = max_rna;
  p_rna_produce.srcPtr.ysize = max_prot;
  p_rna_produce.dstPtr.ptr = devicePitchedPointer_rna_produce.ptr;
  p_rna_produce.dstPtr.pitch = devicePitchedPointer_rna_produce.pitch;
  p_rna_produce.dstPtr.xsize = max_rna;
  p_rna_produce.dstPtr.ysize = max_prot;
  p_rna_produce.extent.width = max_rna * sizeof(int);
  p_rna_produce.extent.height = max_prot;
  p_rna_produce.extent.depth = pop_size;
  p_rna_produce.kind = cudaMemcpyHostToDevice;

  cudaMemcpy3D(&p_rna_produce);


  cudaExtent volumeSizeBytes_influence_enhancing =
      make_cudaExtent(sizeof(double) * max_prot, max_rna, pop_size);

  cudaPitchedPtr devicePitchedPointer_influence_enhancing;

  cudaMalloc3D(&devicePitchedPointer_influence_enhancing, volumeSizeBytes_influence_enhancing);

  cudaMemcpy3DParms p_influence_enhancing = { 0 };

  p_influence_enhancing.srcPtr.ptr = rna_influence_enhancing_tab;
  p_influence_enhancing.srcPtr.pitch = max_prot * sizeof(double);
  p_influence_enhancing.srcPtr.xsize = max_prot;
  p_influence_enhancing.srcPtr.ysize = max_rna;
  p_influence_enhancing.dstPtr.ptr = devicePitchedPointer_influence_enhancing.ptr;
  p_influence_enhancing.dstPtr.pitch = devicePitchedPointer_influence_enhancing.pitch;
  p_influence_enhancing.dstPtr.xsize = max_prot;
  p_influence_enhancing.dstPtr.ysize = max_rna;
  p_influence_enhancing.extent.width = max_prot * sizeof(double);
  p_influence_enhancing.extent.height = max_rna;
  p_influence_enhancing.extent.depth = pop_size;
  p_influence_enhancing.kind = cudaMemcpyHostToDevice;

  cudaMemcpy3D(&p_influence_enhancing);


  cudaExtent volumeSizeBytes_influence_operating =
      make_cudaExtent(sizeof(double) * max_prot, max_rna, pop_size);

  cudaPitchedPtr devicePitchedPointer_influence_operating;

  cudaMalloc3D(&devicePitchedPointer_influence_operating, volumeSizeBytes_influence_operating);

  cudaMemcpy3DParms p_influence_operating = { 0 };

  p_influence_operating.srcPtr.ptr = rna_influence_operating_tab;
  p_influence_operating.srcPtr.pitch = max_prot * sizeof(double);
  p_influence_operating.srcPtr.xsize = max_prot;
  p_influence_operating.srcPtr.ysize = max_rna;
  p_influence_operating.dstPtr.ptr = devicePitchedPointer_influence_enhancing.ptr;
  p_influence_operating.dstPtr.pitch = devicePitchedPointer_influence_enhancing.pitch;
  p_influence_operating.dstPtr.xsize = max_prot;
  p_influence_operating.dstPtr.ysize = max_rna;
  p_influence_operating.extent.width = max_prot * sizeof(double);
  p_influence_operating.extent.height = max_rna;
  p_influence_operating.extent.depth = pop_size;
  p_influence_operating.kind = cudaMemcpyHostToDevice;

  cudaMemcpy3D(&p_influence_operating);
  return max_prot;
}


int transfert_data_to_gpu_thrust(
    int nb_gen,
    int multiply_population,
     int nb_signal, int degradationstep,
     int degradation_rate, double hill_shape_n, double hill_shape,
    int pop_size, int lifestep,
                          std::vector<std::vector<double>*> const &protein_concentration_list,
                          std::vector<std::vector<double>*> const &rna_basal_concentration_list,
                          std::vector<std::vector<std::vector<int>*>*> const &rna_produce_protein_list,
                          std::vector<std::vector<std::vector<double>*>*> const &rna_influence_enhancing_coef_list,
                          std::vector<std::vector<std::vector<double>*>*> const &rna_influence_operating_coef_list,
                                                    std::vector<std::vector<double>*> const &env_concentration_list) {

  int max_prot = 0, max_rna = 0;

  thrust::device_vector<double> gpu_thrust_protein_concentration;
  thrust::device_vector<double> gpu_thrust_rna_basal_concentration;
  thrust::device_vector<int>    gpu_thrust_rna_produce_protein;
  thrust::device_vector<double> gpu_thrust_rna_influence_enhancing_coef;
  thrust::device_vector<double> gpu_thrust_rna_influence_operating_coef;

  thrust::device_vector<double> gpu_thrust_environment_concentration;


  thrust::device_vector<int>    gpu_thrust_nb_protein;
  thrust::device_vector<int>    gpu_thrust_nb_rna_produce_protein;
  thrust::device_vector<int>    gpu_thrust_nb_influence;

  int g_max_rna;
  int g_max_protein;

  double* r_gpu_thrust_protein_concentration;
  double* r_gpu_thrust_rna_basal_concentration;
  int*    r_gpu_thrust_rna_produce_protein;
  double* r_gpu_thrust_rna_influence_enhancing_coef;
  double* r_gpu_thrust_rna_influence_operating_coef;

  double* r_gpu_thrust_environment_concentration;


  int* r_gpu_thrust_nb_protein;
  int* r_gpu_thrust_nb_rna_produce_protein;
  int* r_gpu_thrust_nb_influence;

  for (int i = 0; i < pop_size; i++) {
    max_prot = protein_concentration_list[i]->size() > max_prot ?
               protein_concentration_list[i]->size() : max_prot;
    max_rna = rna_basal_concentration_list[i]->size() > max_rna ?
              rna_basal_concentration_list[i]->size() : max_rna;
  }

  g_max_protein = max_prot;
  g_max_rna = max_rna;

  thrust::host_vector<double> thrust_protein_concentration(max_prot*pop_size);
  thrust::host_vector<double> thrust_rna_basal_concentration(max_rna*pop_size);
  thrust::host_vector<int> thrust_rna_produce_protein(max_rna*max_prot*pop_size);
  thrust::host_vector<double> thrust_rna_influence_enhancing_coef(max_rna*max_prot*pop_size);
  thrust::host_vector<double> thrust_rna_influence_operating_coef(max_rna*max_prot*pop_size);

  thrust::host_vector<int> thrust_nb_protein(pop_size);
  thrust::host_vector<int> thrust_nb_rna_produce_protein(max_prot*pop_size);
  thrust::host_vector<int> thrust_nb_influence(max_rna*pop_size);

  for (int i = 0; i < pop_size; i++){

    thrust_nb_protein[i] = protein_concentration_list[i]->size();


    thrust::copy(
                                        protein_concentration_list[i]->begin(),
                                        protein_concentration_list[i]->end(),
                                        thrust_protein_concentration.begin()+max_prot*i
                                        );

    thrust::copy(
                                        rna_basal_concentration_list[i]->begin(),
                                        rna_basal_concentration_list[i]->end(),
                                          thrust_rna_basal_concentration.begin()+max_rna*i
    );

    for (int prot_id = 0; prot_id < rna_produce_protein_list[i]->size(); prot_id++) {
      thrust_nb_rna_produce_protein[i*max_prot+prot_id] = rna_produce_protein_list[i]->at(prot_id)->size();
      thrust::copy(
                                          rna_produce_protein_list[i]->at(prot_id)->begin(),
                                          rna_produce_protein_list[i]->at(prot_id)->end(),
                                          thrust_rna_produce_protein.begin()+max_prot*max_rna*i+max_rna*prot_id
      );
    }

    for (int rna_id = 0; rna_id < rna_influence_enhancing_coef_list[i]->size(); rna_id++) {
      thrust_nb_influence[i*max_rna+rna_id] = rna_influence_enhancing_coef_list[i]->at(rna_id)->size();
      thrust::copy(
                                          rna_influence_enhancing_coef_list[i]->at(rna_id)->begin(),
                                          rna_influence_enhancing_coef_list[i]->at(rna_id)->end(),
                                          thrust_rna_influence_enhancing_coef.begin()+
                                                     max_prot*max_rna*i+max_prot*rna_id
      );
    }

    for (int rna_id = 0; rna_id < rna_influence_operating_coef_list[i]->size(); rna_id++) {
      thrust::copy(
                                                 rna_influence_operating_coef_list[i]->at(rna_id)->begin(),
                                                 rna_influence_operating_coef_list[i]->at(rna_id)->end(),
                                                 thrust_rna_influence_operating_coef.begin()+
                                                     max_prot*max_rna*i+max_prot*rna_id
      );

    }
  }

  thrust::host_vector<double> thrust_environment_concentration(lifestep*nb_signal);

  for (int i = 0; i < lifestep; i++) {
    thrust::copy(
                                            env_concentration_list[i]->begin(),
                                            env_concentration_list[i]->end(),
                                            thrust_environment_concentration.begin()+
                                               nb_signal*i
    );
  }

  printf("Copying to device");

  gpu_thrust_protein_concentration = thrust_protein_concentration;
  gpu_thrust_rna_basal_concentration = thrust_rna_basal_concentration;
  gpu_thrust_rna_produce_protein = thrust_rna_produce_protein;
  gpu_thrust_rna_influence_enhancing_coef = thrust_rna_influence_enhancing_coef;
  gpu_thrust_rna_influence_operating_coef = thrust_rna_influence_operating_coef;

  gpu_thrust_environment_concentration = thrust_environment_concentration;


  gpu_thrust_nb_protein = thrust_nb_protein;
  gpu_thrust_nb_rna_produce_protein = thrust_nb_rna_produce_protein;
  gpu_thrust_nb_influence = thrust_nb_influence;

  r_gpu_thrust_protein_concentration = thrust::raw_pointer_cast(gpu_thrust_protein_concentration.data());
  r_gpu_thrust_rna_basal_concentration = thrust::raw_pointer_cast(gpu_thrust_rna_basal_concentration.data());
  r_gpu_thrust_rna_produce_protein = thrust::raw_pointer_cast(gpu_thrust_rna_produce_protein.data());
  r_gpu_thrust_rna_influence_enhancing_coef = thrust::raw_pointer_cast(gpu_thrust_rna_influence_enhancing_coef.data());
  r_gpu_thrust_rna_influence_operating_coef = thrust::raw_pointer_cast(gpu_thrust_rna_influence_operating_coef.data());

  r_gpu_thrust_environment_concentration = thrust::raw_pointer_cast(gpu_thrust_environment_concentration.data());


  r_gpu_thrust_nb_protein = thrust::raw_pointer_cast(gpu_thrust_nb_protein.data());
  r_gpu_thrust_nb_rna_produce_protein = thrust::raw_pointer_cast(gpu_thrust_nb_rna_produce_protein.data());
  r_gpu_thrust_nb_influence = thrust::raw_pointer_cast(gpu_thrust_nb_influence.data());

  printf("Launching kernel\n");


  for (int gen = 0; gen < nb_gen; gen++) {

    process_delta_thrust << < 1024 * multiply_population, g_max_protein >> >
                                                          (nb_signal, degradationstep, degradation_rate,
                                                              g_max_protein, g_max_rna,
                                                              r_gpu_thrust_nb_rna_produce_protein, r_gpu_thrust_rna_produce_protein,
                                                              r_gpu_thrust_protein_concentration,
                                                              r_gpu_thrust_rna_basal_concentration, r_gpu_thrust_nb_protein,
                                                              r_gpu_thrust_rna_influence_enhancing_coef, r_gpu_thrust_rna_influence_operating_coef,
                                                              r_gpu_thrust_nb_influence,
                                                              r_gpu_thrust_environment_concentration, hill_shape, hill_shape_n);
  }
  return max_prot;
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



__global__
void process_delta_float(int nb_signal, int degradstep, int degradrate, int ***rna_produce_protein_array,
                   int **nb_rna_produce_protein, int *nb_rna_produce,   float **protein_concentration_array,
                   float **rna_basal_concentration_array, int *nb_protein_array, int *nb_rna_array,
                   float ***rna_influence_enhancing_coef_array, float ***rna_influence_operating_coef_array,
                   int **nb_rna_influence_enhancing_coef, int  **nb_rna_influence_operating_coef,
                   float **env_concentration_array, float hill_shape, float hill_shape_n) {



  float delta = 0;

  int indiv_id = blockIdx.x;
  int prot_id = threadIdx.x;

  if (prot_id < nb_protein_array[indiv_id] - nb_signal) {
    for (int j = 0; j < degradstep; j++) {
      for (int j = 0;
           j < nb_rna_produce_protein[indiv_id][prot_id]; j++) {
        float enhancer_activity = 0;
        float operator_activity = 0;

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

        float enhancer_activity_pow_n = enhancer_activity == 0 ? 0 :
                                         powf(enhancer_activity, hill_shape_n);
        float operator_activity_pow_n = operator_activity == 0 ? 0 :
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
      delta *= 1 / (float) degradstep;

      __syncthreads();

      protein_concentration_array[indiv_id][prot_id] = delta;
    }
  }
}


__global__
void process_delta_dense(int nb_signal, int degradstep, int degradrate, int ***rna_produce_protein_array,
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


__global__
void process_delta_thrust(int nb_signal, int degradstep, int degradrate,
                          int max_prot, int max_rna,
                   int *nb_rna_produce_protein, int *rna_produce_protein,
                   double *protein_concentration_array,
                   double *rna_basal_concentration_array, int *nb_protein_array,
                   double *rna_influence_enhancing_coef_array, double *rna_influence_operating_coef_array,
                   int *nb_rna_influence_enhancing_coef,
                   double *env_concentration_array, double hill_shape, double hill_shape_n) {



  double delta = 0;

  int indiv_id = blockIdx.x;
  int prot_id = threadIdx.x;

  if (prot_id < nb_protein_array[indiv_id] - nb_signal) {
    for (int j = 0; j < degradstep; j++) {
      for (int j = 0;
           j < nb_rna_produce_protein[indiv_id*max_prot+prot_id]; j++) {
        double enhancer_activity = 0;
        double operator_activity = 0;

        int rna_id = rna_produce_protein[indiv_id*max_prot*max_rna+prot_id*max_rna+j];

        for (int i = 0; i <
                        nb_rna_influence_enhancing_coef[indiv_id*max_rna+rna_id]; i++) {

          enhancer_activity +=
              rna_influence_enhancing_coef_array[indiv_id*max_rna*max_prot+rna_id*max_prot+i]
              * protein_concentration_array[indiv_id*max_prot+i];
          operator_activity +=
              rna_influence_operating_coef_array[indiv_id*max_rna*max_prot+rna_id*max_prot+i]
              * protein_concentration_array[indiv_id*max_prot+i];
        }

        double enhancer_activity_pow_n = enhancer_activity == 0 ? 0 :
                                         powf(enhancer_activity, hill_shape_n);
        double operator_activity_pow_n = operator_activity == 0 ? 0 :
                                         powf(operator_activity, hill_shape_n);
        delta += rna_basal_concentration_array[indiv_id*max_rna+rna_id]
                 * (hill_shape
                    / (operator_activity_pow_n + hill_shape))
                 * (1 +
                    ((1 / rna_basal_concentration_array[indiv_id*max_rna+rna_id]
                     ) -
                     1)
                    * (enhancer_activity_pow_n /
                       (enhancer_activity_pow_n + hill_shape)));
      }

      delta -=
          degradrate *
          protein_concentration_array[indiv_id*max_prot+prot_id];
      delta *= 1 / (double) degradstep;

      __syncthreads();

      protein_concentration_array[indiv_id*max_prot+prot_id] = delta;
    }
  }
}
