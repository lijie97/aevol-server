//
// Created by arrouan on 13/12/16.
//

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <getopt.h>
#include <chrono>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <zlib.h>

#include "ode_solver_only_omp_gpu.h"
#include "ode_solver_only.h"
using namespace std::chrono;

//#include "aevol.h"


#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char* argv[])
{
  bool verbose = false;

  std::string line;

  const char * short_options = "Vve:p:s:d:m:o:n:";
  static struct option long_options[] = {
      {"version",   no_argument,       NULL,  'V'},
      {"verbose",   no_argument,       NULL,  'v'},
      {"executionmode",       required_argument,  NULL, 'e'},
      {"multipopulation",       required_argument,  NULL, 'p'},
      {"lifestep",       required_argument,  NULL, 's'},
      {"degradationstep",       required_argument,  NULL, 'd'},
      {"merge",       required_argument,  NULL, 'm'},
      {"odesolver",       required_argument,  NULL, 'o'},
      {"parallel",       required_argument,  NULL, 'n'},
      {0, 0, 0, 0}
  };

  int execution_mode = 0; // 0: sequential, 1: sequential_merge, 2: parallel, 3: parallel:merge
  int multiply_population = 0;
  int lifestep = 0;
  int merge = 0;
  int odesolver = 0;
  int parallelthread = 1;

  int option;
  while((option = getopt_long(argc, argv, short_options, long_options, NULL)) != -1)
  {
    switch(option)
    {
      case 'V' :
      {
        exit(EXIT_SUCCESS);
      }
      case 'v' : verbose = true;                    break;
      case 'e' :
      {
        if (strcmp(optarg, "") == 0)
        {
          printf("%s: error: Option -e or --executionmode : missing argument.\n", argv[0]);
          exit(EXIT_FAILURE);
        }

        execution_mode = atoi(optarg);

        break;
      }
      case 'p' :
      {
        if (strcmp(optarg, "") == 0)
        {
          printf("%s: error: Option -p or --multipopulation : missing argument.\n", argv[0]);
          exit(EXIT_FAILURE);
        }

        multiply_population = atoi(optarg);

        break;
      }
      case 's' :
      {
        if (strcmp(optarg, "") == 0)
        {
          printf("%s: error: Option -s or --lifestep : missing argument.\n", argv[0]);
          exit(EXIT_FAILURE);
        }

        lifestep = atoi(optarg);

        break;
      }
      case 'd' :
      {
        if (strcmp(optarg, "") == 0)
        {
          printf("%s: error: Option -d or --degradationstep : missing argument.\n", argv[0]);
          exit(EXIT_FAILURE);
        }

        degradationstep = atoi(optarg);

        break;
      }
      case 'm' :
      {
        if (strcmp(optarg, "") == 0)
        {
          printf("%s: error: Option -m or --merge : missing argument.\n", argv[0]);
          exit(EXIT_FAILURE);
        }

        merge = atoi(optarg);

        break;
      }
      case 'o' :
      {
        if (strcmp(optarg, "") == 0)
        {
          printf("%s: error: Option -o or --odesolver : missing argument.\n", argv[0]);
          exit(EXIT_FAILURE);
        }

        odesolver = atoi(optarg);

        break;
      }
      case 'n' :
      {
        if (strcmp(optarg, "") == 0)
        {
          printf("%s: error: Option -n or --parallel : missing argument.\n", argv[0]);
          exit(EXIT_FAILURE);
        }

        parallelthread = atoi(optarg);

        break;
      }
    }
  }

//  omp_set_num_threads(parallelthread);

  std::ifstream protein_concentration_file("protein_concentration.csv");
  std::ifstream rna_basal_concentration_file("rna_basal_concentration.csv");
  std::ifstream rna_produce_protein_file("rna_produce_protein.csv");
  std::ifstream rna_influence_enhancing_coef("rna_influence_enhancing_coef.csv");
  std::ifstream rna_influence_operating_coef("rna_influence_operating_coef.csv");
  std::ifstream nb_protein_signal("nb_protein_signal.csv");
  std::ifstream env_concentration("env_concentration_list.csv");

  getline (nb_protein_signal,line);
  nb_signal = atoi(line.c_str());
  nb_protein_signal.close();

  protein_concentration_list.reserve(1024);
  rna_basal_concentration_list.reserve(1024);
  rna_produce_protein_list.reserve(1024);
  rna_influence_enhancing_coef_list.reserve(1024);
  rna_influence_operating_coef_list.reserve(1024);

  std::vector<std::vector<double>*>::iterator itprot = protein_concentration_list.begin();
  std::vector<std::vector<std::vector<int>*>*>::iterator itproduce = rna_produce_protein_list.begin();
  std::vector<std::vector<double>*>::iterator itrna = rna_basal_concentration_list.begin();
  std::vector<std::vector<std::vector<double>*>*>::
  iterator itrnaenh = rna_influence_enhancing_coef_list.begin();
  std::vector<std::vector<std::vector<double>*>*>::
  iterator itrnaop = rna_influence_operating_coef_list.begin();

  for (int i = 0; i < 1024; i++) {
    protein_concentration_list.insert(itprot+i,new std::vector<double>());
    rna_basal_concentration_list.insert(itrna+i,new std::vector<double>());
    rna_produce_protein_list.insert(itproduce+i,new std::vector<std::vector<int>*>());
    rna_influence_enhancing_coef_list.insert(itrnaenh+i,new std::vector<std::vector<double>*>());
    rna_influence_operating_coef_list.insert(itrnaop+i,new std::vector<std::vector<double>*>());
  }


  while ( getline (protein_concentration_file,line) ) {
    std::vector<double> vect;

    std::stringstream ss(line);

    double i;

    while (ss >> i) {
      vect.push_back(i);

      if (ss.peek() == ',')
        ss.ignore();
    }

    std::vector<double>::iterator itprot2 = protein_concentration_list[(int)vect[0]]->begin();
    protein_concentration_list[(int)vect[0]]->insert(itprot2+(int)vect[1],vect[2]);

    std::vector<std::vector<int>*>::iterator itprod2 = rna_produce_protein_list[(int)vect[0]]->begin();
    rna_produce_protein_list[(int)vect[0]]->insert(itprod2+(int)vect[1],new std::vector<int>());
  }

  while ( getline (rna_produce_protein_file,line) ) {
    std::vector<double> vect;

    std::stringstream ss(line);

    int i;

    while (ss >> i) {
      vect.push_back(i);

      if (ss.peek() == ',')
        ss.ignore();
    }

    std::vector<int>::iterator itproduce2 = rna_produce_protein_list[(int)vect[0]]->at((int)vect[1])->begin();
    rna_produce_protein_list[(int)vect[0]]->at((int)vect[1])->push_back((int)vect[2]);
  }

  while ( getline (rna_basal_concentration_file,line) ) {
    std::vector<double> vect;

    std::stringstream ss(line);

    double i;

    while (ss >> i) {
      vect.push_back(i);

      if (ss.peek() == ',')
        ss.ignore();
    }

    std::vector<double>::iterator itrna2 = rna_basal_concentration_list[(int)vect[0]]->begin();
    rna_basal_concentration_list[(int)vect[0]]->insert(itrna2+(int)vect[1],vect[2]);

    std::vector<std::vector<double>*>::iterator itrnaenh2 =
        rna_influence_enhancing_coef_list[(int)vect[0]]->begin();
    std::vector<std::vector<double>*>::iterator itrnaop2 =
        rna_influence_operating_coef_list[(int)vect[0]]->begin();

    rna_influence_enhancing_coef_list[(int)vect[0]]->insert(itrnaenh2+(int)vect[1],new std::vector<double>());
    rna_influence_operating_coef_list[(int)vect[0]]->insert(itrnaop2+(int)vect[1],new std::vector<double>());
  }

  while ( getline (rna_influence_enhancing_coef,line) ) {
    std::vector<double> vect;

    std::stringstream ss(line);

    double i;

    while (ss >> i) {
      vect.push_back(i);

      if (ss.peek() == ',')
        ss.ignore();
    }

    std::vector<double>::iterator itrna2 =
        rna_influence_enhancing_coef_list[(int)vect[0]]->at((int)vect[1])->begin();

    rna_influence_enhancing_coef_list[(int)vect[0]]->at((int)vect[1])->insert(itrna2+(int)vect[2],vect[3]);
  }


  while ( getline (rna_influence_operating_coef,line) ) {
    std::vector<double> vect;

    std::stringstream ss(line);

    double i;

    while (ss >> i) {
      vect.push_back(i);

      if (ss.peek() == ',')
        ss.ignore();
    }


    std::vector<double>::iterator itrna2 =
        rna_influence_operating_coef_list[(int)vect[0]]->at((int)vect[1])->begin();

    rna_influence_operating_coef_list[(int)vect[0]]->at((int)vect[1])->insert(itrna2+(int)vect[2],vect[3]);
  }

  env_concentration_list.reserve(lifestep);
  std::vector<std::vector<double>*>::iterator itenv = env_concentration_list.begin();
  for (int i = 0; i < lifestep; i++) {
    env_concentration_list.insert(itenv+i,new std::vector<double>());
  }

  while ( getline (env_concentration,line) ) {
    std::vector<double> vect;

    std::stringstream ss(line);

    double i;

    while (ss >> i) {
      vect.push_back(i);

      if (ss.peek() == ',')
        ss.ignore();
    }

    std::vector<double>::iterator itenv2 = env_concentration_list[(int)vect[0]]->begin();
    env_concentration_list[(int)vect[0]]->insert(itenv2+(int)vect[1],vect[2]);
  }

  printf("Lifestep %d -- Population %d (%ld)\n",lifestep,
         multiply_population*1024,protein_concentration_list.size());
  high_resolution_clock::time_point t_t1 = high_resolution_clock::now();
  high_resolution_clock::time_point t_t2,t1,t2;
  printf("Starting lookup\n");

  char* lookup_table_file_name = new char[100];

      sprintf( lookup_table_file_name, "lookup_table.ae" );

      gzFile lookup_table_file = gzopen( lookup_table_file_name, "r" );

      if ( lookup_table_file == Z_NULL )
      {
        printf( "ERROR : Could not read lookup table file %s\n", lookup_table_file_name );
        exit( EXIT_FAILURE );
      }

  printf("Starting filling\n");

  double value;
      for (int i=0; i < LOOKUP_TABLE_SIZE; i++) {
        gzread( lookup_table_file, &value, sizeof(double));
        lookup_table_pow[i] = (double) value;
      }

  printf("Closingr\n");

  gzclose( lookup_table_file );

      delete[] lookup_table_file_name;

      printf("Starting transfer\n");
      transfer_to_tab(1024, lifestep);
      compute_openmp_gpu(lifestep,degradationstep,nb_signal,hill_shape_n,hill_shape,degradation_rate);

}


void compute_openmp_gpu(int lifestep, int degradationstep, int nb_signal,
                        double hill_shape,double hill_shape_n,double degradation_rate ) {
  int pop_size = 1024;

  double* delta = (double*) malloc(1024 * max_prot * sizeof(double));

  double* protein_tab;
  double* rna_basal_tab;
  int* rna_produce_protein_tab;
  double* rna_influence_enhancing_tab;
  double* rna_influence_operating_tab;
  double* env_concentration_tab;
  int* nb_protein;

  for (int i = 0; i < pop_size; i++) {
    max_prot = protein_concentration_list[i]->size() > max_prot ?
               protein_concentration_list[i]->size() : max_prot;
    max_rna = rna_basal_concentration_list[i]->size() > max_rna ?
              rna_basal_concentration_list[i]->size() : max_rna;
  }

  /* Storing as dense tab */

  protein_tab = (double*) malloc(pop_size * max_prot * sizeof(double));
  memset(protein_tab, -1, pop_size * max_prot * sizeof(double));


  nb_protein = (int*) malloc(pop_size * sizeof(int));
  memset(nb_protein, -1, pop_size * sizeof(int));


  rna_basal_tab = (double*) malloc(pop_size * max_rna * sizeof(double));
  memset(rna_basal_tab, -1, pop_size * max_rna * sizeof(double));

  rna_produce_protein_tab = (int*) malloc(
      pop_size * max_prot * max_rna * sizeof(int));
  memset(rna_produce_protein_tab, -1,
         pop_size * max_prot * max_rna * sizeof(int));


  rna_influence_enhancing_tab = (double*) malloc(
      pop_size * max_prot * max_rna * sizeof(double));
  memset(rna_influence_enhancing_tab, 0,
         pop_size * max_prot * max_rna * sizeof(double));

  rna_influence_operating_tab = (double*) malloc(
      pop_size * max_prot * max_rna * sizeof(double));
  memset(rna_influence_operating_tab, 0,
         pop_size * max_prot * max_rna * sizeof(double));

  env_concentration_tab = (double*) malloc(
      lifestep * nb_signal * sizeof(double));
  memset(env_concentration_tab, 0, lifestep * nb_signal * sizeof(double));

  for (int i = 0; i < pop_size; i++) {
    nb_protein[i] = protein_concentration_list[i]->size();

    for (int prot_id = 0;
         prot_id < protein_concentration_list[i]->size(); prot_id++) {
      protein_tab[i * max_prot + prot_id] =
          protein_concentration_list[i]->at(prot_id);
    }

    for (int rna_id = 0;
         rna_id < rna_basal_concentration_list[i]->size(); rna_id++) {
      rna_basal_tab[i * max_rna + rna_id] =
          rna_basal_concentration_list[i]->at(rna_id);
    }

    for (int prot_id = 0;
         prot_id < rna_produce_protein_list[i]->size(); prot_id++) {
      for (int ix = 0;
           ix < rna_produce_protein_list[i]->at(prot_id)->size(); ix++) {
        rna_produce_protein_tab[i * max_prot + prot_id * max_rna + ix] =
            rna_produce_protein_list[i]->at(prot_id)->at(ix);

      }
    }

    for (int rna_id = 0;
         rna_id < rna_influence_enhancing_coef_list[i]->size(); rna_id++) {
      for (int prot_id = 0; prot_id <
                            rna_influence_enhancing_coef_list[i]->at(
                                rna_id)->size(); prot_id++) {
        rna_influence_enhancing_tab[i * max_rna + rna_id * max_prot + prot_id] =
            rna_influence_enhancing_coef_list[i]->at(rna_id)->at(prot_id);
      }
    }

    for (int rna_id = 0;
         rna_id < rna_influence_operating_coef_list[i]->size(); rna_id++) {
      for (int prot_id = 0; prot_id <
                            rna_influence_operating_coef_list[i]->at(
                                rna_id)->size(); prot_id++) {
        rna_influence_operating_tab[i * max_rna + rna_id * max_prot + prot_id] =
            rna_influence_operating_coef_list[i]->at(rna_id)->at(prot_id);
      }
    }
  }

  for (int i = 0; i < lifestep; i++) {
    for (int j = 0; j < nb_signal; j++) {
      env_concentration_tab[i * nb_signal + j] = env_concentration_list[i]->at(
          j);
    }
  }

#pragma omp target data map(to:delta[0:1024*max_prot],protein_tab[0:1024*max_prot],\
      rna_basal_tab[0:1024*max_rna],\
      rna_produce_protein_tab[0:1024*max_prot*max_rna],\
      rna_influence_enhancing_tab[0:1024*max_prot*max_rna],\
      rna_influence_operating_tab[0:1024*max_prot*max_rna],\
      env_concentration_tab[0:lifestep*nb_signal],\
      nb_protein[0:1024],\
      lookup_table_pow[0:LOOKUP_TABLE_SIZE])
  {
    //for (int lstep = 0; lstep < lifestep; lstep++) {
//teams distribute parallel for collapse(2) schedule(static,1)
    for (int l = 0; l < lifestep; l++) {
#pragma omp target teams distribute parallel for schedule(static,1) num_teams(1024) num_threads(max_prot)
      for (int idx = 0; idx < 1024 * max_prot; idx++) {
        int indiv_id = idx / 1024;
        int prot_id = idx % 1024;
        //for (int indiv_id = 0; indiv_id < 1024; indiv_id++) {
        //for (int prot_id = nb_protein[indiv_id] - nb_signal;
        //     prot_id < nb_protein[indiv_id]; prot_id++) {
        protein_tab[indiv_id * max_prot + prot_id] =
            env_concentration_tab[lifestep * nb_signal + prot_id -
                                  nb_protein[indiv_id] + nb_signal];
        //}
      }

      for (int j = 0; j < degradationstep; j++) {
#pragma omp target teams distribute parallel for schedule(static,1) num_teams(1024) num_threads(max_prot)
        for (int idx = 0; idx < 1024 * max_prot; idx++) {
          int indiv_id = idx / 1024;
          int prot_id = idx % 1024;
          delta[indiv_id * max_prot + prot_id] = 0;
        }

//#pragma omp distribute parallel for collapse(3) schedule(static,1)
#pragma omp target teams distribute parallel for schedule(static,1) num_teams(1024) num_threads(max_prot)
        for (int idx = 0; idx < 1024 * max_prot * max_rna; idx++) {
          int indiv_id = idx / (1024 * max_prot);
          int prot_id = (idx % (1024 * max_prot)) / max_rna;
          int j = (idx % (1024 * max_prot)) % max_rna;


          int rna_id = rna_produce_protein_tab[indiv_id * max_prot +
                                               prot_id * max_rna + j];
          if (rna_id != -1) {

            double enhancer_activity = 0;
            double operator_activity = 0;
            for (int i = 0; i <
                            max_rna * max_prot; i++) {


              enhancer_activity +=
                  rna_influence_enhancing_tab[indiv_id * max_rna +
                                              rna_id * max_prot + i]
                  * protein_tab[indiv_id * max_prot + i];
              operator_activity +=
                  rna_influence_operating_tab[indiv_id * max_rna +
                                              rna_id * max_prot + i]
                  * protein_tab[indiv_id * max_prot + i];
            }

            double enhancer_activity_pow_n = enhancer_activity == 0 ? 0 :
                                             enhancer_activity <= 1 ?
                                             (lookup_table_pow[(int) (
                                                 enhancer_activity *
                                                 LOOKUP_TABLE_SIZE)] +
                                              lookup_table_pow[((int) (
                                                  enhancer_activity *
                                                  LOOKUP_TABLE_SIZE)) +
                                                               1]) / 2 : 0;
            double operator_activity_pow_n = operator_activity == 0 ? 0 :
                                             operator_activity <= 1 ?
                                             (lookup_table_pow[(int) (
                                                 operator_activity *
                                                 LOOKUP_TABLE_SIZE)] +
                                              lookup_table_pow[((int) (
                                                  operator_activity *
                                                  LOOKUP_TABLE_SIZE)) +
                                                               1]) / 2 : 0;
            delta[indiv_id * max_prot + prot_id] +=
                rna_basal_tab[indiv_id * max_rna + rna_id]
                * (hill_shape
                   / (operator_activity_pow_n + hill_shape))
                * (1 +
                   ((1 /
                     rna_basal_tab[indiv_id * max_rna + rna_id]) -
                    1)
                   * (enhancer_activity_pow_n /
                      (enhancer_activity_pow_n + hill_shape)));

          }
        }

#pragma omp target teams distribute parallel for schedule(static,1) num_teams(1024) num_threads(max_prot)
        for (int idx = 0; idx < 1024 * max_prot; idx++) {
          int indiv_id = idx / 1024;
          int prot_id = idx % 1024;

          delta[indiv_id * max_prot + prot_id] -=
              degradation_rate *
              protein_tab[indiv_id * max_prot + prot_id];
          delta[indiv_id * max_prot + prot_id] *=
              1 / (double) degradationstep;
        }

        // Apply the changes in concentrations we have just computed
#pragma omp target teams distribute parallel for schedule(static,1) num_teams(1024) num_threads(max_prot)
        for (int idx = 0; idx < 1024 * max_prot; idx++) {
          int indiv_id = idx / 1024;
          int prot_id = idx % 1024;

          protein_tab[indiv_id * max_prot +
                      prot_id] += delta[indiv_id * max_prot + prot_id];

        }
      }

    }
  }
}




void transfer_to_tab(int pop_size, int lifestep) {
/*  for (int i = 0; i < pop_size; i++) {
    max_prot = protein_concentration_list[i]->size() > max_prot ?
               protein_concentration_list[i]->size() : max_prot;
    max_rna  = rna_basal_concentration_list[i]->size() > max_rna ?
               rna_basal_concentration_list[i]->size() : max_rna;
  }

  omp_protein_tab = (double*)malloc(pop_size * max_prot * sizeof(double));
  memset (omp_protein_tab, -1, pop_size * max_prot * sizeof(double));


  omp_rna_basal_tab = (double*)malloc(pop_size * max_rna * sizeof(double));
  memset (omp_rna_basal_tab, -1, pop_size * max_rna * sizeof(double));

  omp_rna_produce_protein_tab = (int*)malloc(pop_size * max_prot * max_rna * sizeof(int));
  memset (omp_rna_produce_protein_tab, -1, pop_size * max_prot * max_rna * sizeof(int));


  omp_rna_influence_enhancing_tab = (double*)malloc(pop_size * max_prot * max_rna * sizeof(double));
  memset (omp_rna_influence_enhancing_tab, 0, pop_size * max_prot * max_rna * sizeof(double));

  omp_rna_influence_operating_tab = (double*)malloc(pop_size * max_prot * max_rna * sizeof(double));
  memset (omp_rna_influence_operating_tab, 0, pop_size * max_prot * max_rna * sizeof(double));

  for (int i = 0; i < pop_size; i++) {
    for (int prot_id = 0; prot_id < protein_concentration_list[i]->size(); prot_id++) {
      omp_protein_tab[i*max_prot+prot_id] =
          protein_concentration_list[i]->at(prot_id);
    }

    for (int rna_id = 0; rna_id < rna_basal_concentration_list[i]->size(); rna_id++) {
      omp_rna_basal_tab[i*max_rna+rna_id] =
          rna_basal_concentration_list[i]->at(rna_id);
    }

    for (int prot_id = 0; prot_id < rna_produce_protein_list[i]->size(); prot_id++) {
      for (int ix = 0; ix < rna_produce_protein_list[i]->at(prot_id)->size(); ix++) {
        omp_rna_produce_protein_tab[i*max_prot+prot_id*max_rna+ix] =
            rna_produce_protein_list[i]->at(prot_id)->at(ix);

      }
    }

    for (int rna_id = 0; rna_id < rna_influence_enhancing_coef_list[i]->size(); rna_id++) {
      for (int prot_id = 0; prot_id <
                            rna_influence_enhancing_coef_list[i]->at(rna_id)->size(); prot_id++) {
        omp_rna_influence_enhancing_tab[i*max_rna+rna_id*max_prot+prot_id] =
            rna_influence_enhancing_coef_list[i]->at(rna_id)->at(prot_id);
      }
    }

    for (int rna_id = 0; rna_id < rna_influence_operating_coef_list[i]->size(); rna_id++) {
      for (int prot_id = 0; prot_id <
                            rna_influence_operating_coef_list[i]->at(rna_id)->size(); prot_id++) {
        omp_rna_influence_operating_tab[i*max_rna+rna_id*max_prot+prot_id] =
            rna_influence_operating_coef_list[i]->at(rna_id)->at(prot_id);
      }
    }
  }

  omp_env_concentration_tab = (double*)malloc(lifestep * nb_signal * sizeof(double));
  for (int i = 0; i < lifestep; i++) {
    for (int j=0; j < nb_signal; j++) {
      omp_env_concentration_tab[i*nb_signal+j] = env_concentration_list[i]->at(j);
    }
  }*/

  /* Computing max rna and protein */
  /*for (int i = 0; i < pop_size; i++) {
    max_prot = protein_concentration_list[i]->size() > max_prot ?
               protein_concentration_list[i]->size() : max_prot;
    max_rna  = rna_basal_concentration_list[i]->size() > max_rna ?
               rna_basal_concentration_list[i]->size() : max_rna;
  }

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
*/
}
