//
// Created by arrouan on 13/12/16.
//

#include <vector>
#include <map>
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
#include <omp.h>
#include <zlib.h>
#include "tbb/tbb.h"

using namespace tbb;
constexpr int32_t LOOKUP_TABLE_SIZE = 10000000;

#include "ode_solver_only_omp_gpu.h"
#include "ode_solver_only.h"
using namespace std::chrono;

//#include "aevol.h"

#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void loop_one_user(int lifestep, int i);
void solve_one_indiv_one_step(int indiv_id);
void update_env_indiv(int lifestep, int indiv_id);
void one_prot(std::vector<double> &delta, int indiv_id, int prot_id);
void solve_one_indiv_one_step_parallel(int indiv_id);
void loop_one_user_parallel(int lifestep, int i);
void update_env_indiv_parallel(int lifestep, int indiv_id);
void update_one_env(int lifestep, int indiv_id, int prot_id);
void update_one_prot(std::vector<double> &delta, int indiv_id, int prot_id);

int main(int argc, char* argv[])
{
  bool verbose = false;

  std::string line;

  const char * short_options = "Vve:p:s:d:m:o:n:g:";
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
      {"generation",       required_argument,  NULL, 'g'},
      {0, 0, 0, 0}
  };

  int execution_mode = 0; // 0: sequential, 1: sequential_merge, 2: parallel, 3: parallel:merge
  int multiply_population = 0;
  int lifestep = 0;
  int merge = 0;
  int odesolver = 0;
  int parallelthread = 1;
  int nb_gen = 1;

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
      case 'g' :
      {
        if (strcmp(optarg, "") == 0)
        {
          printf("%s: error: Option -g or --generation : missing argument.\n", argv[0]);
          exit(EXIT_FAILURE);
        }

        nb_gen = atoi(optarg);

        break;
      }
    }
  }

  omp_set_num_threads(parallelthread);

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


  if (multiply_population > 1) {
    std::vector<std::vector<double>*> m_protein_concentration_list;
    std::vector<std::vector<double>*> m_rna_basal_concentration_list;
    std::vector<std::vector<std::vector<int>*>*> m_rna_produce_protein;
    std::vector<std::vector<std::vector<double>*>*> m_rna_influence_enhancing_coef_list;
    std::vector<std::vector<std::vector<double>*>*> m_rna_influence_operating_coef_list;


    m_protein_concentration_list.reserve(1024*multiply_population);
    m_rna_basal_concentration_list.reserve(1024*multiply_population);
    m_rna_produce_protein.reserve(1024*multiply_population);
    m_rna_influence_enhancing_coef_list.reserve(1024*multiply_population);
    m_rna_influence_operating_coef_list.reserve(1024*multiply_population);


    std::vector<std::vector<double>*>::iterator itprot = m_protein_concentration_list.begin();
    std::vector<std::vector<std::vector<int>*>*>::iterator itproduce = m_rna_produce_protein.begin();
    std::vector<std::vector<double>*>::iterator itrna = m_rna_basal_concentration_list.begin();
    std::vector<std::vector<std::vector<double>*>*>::
    iterator itrnaenh = m_rna_influence_enhancing_coef_list.begin();
    std::vector<std::vector<std::vector<double>*>*>::
    iterator itrnaop = m_rna_influence_operating_coef_list.begin();

    for (int i = 0; i < 1024*multiply_population; i++) {
      m_protein_concentration_list.insert(itprot+i,new std::vector<double>());
      m_rna_basal_concentration_list.insert(itrna+i,new std::vector<double>());
      m_rna_produce_protein.insert(itproduce+i,new std::vector<std::vector<int>*>());
      m_rna_influence_enhancing_coef_list.insert(itrnaenh+i,new std::vector<std::vector<double>*>());
      m_rna_influence_operating_coef_list.insert(itrnaop+i,new std::vector<std::vector<double>*>());
    }

    for (int m = 0; m < multiply_population; m++) {
      for (int i = 0+m*multiply_population; i < 1024+m*multiply_population; i++) {
        int old_i = i - m*multiply_population;

        int prot_id = 0;
        std::vector<double>::iterator itprot = m_protein_concentration_list[i]->begin();
        for (auto concentration : *(protein_concentration_list[old_i])) {
          m_protein_concentration_list[i]->push_back(concentration);
          prot_id++;
        }

        int rna_id = 0;
        std::vector<double>::iterator itprod = m_rna_basal_concentration_list[i]->begin();
        for (auto concentration : *(rna_basal_concentration_list[old_i])) {
          m_rna_basal_concentration_list[i]->push_back(concentration);
          rna_id++;
        }


        std::vector<std::vector<int>*>::iterator itproduce =
            m_rna_produce_protein[i]->begin();
        prot_id = 0;

        for (auto rna_produce : *(rna_produce_protein_list[old_i])) {
          m_rna_produce_protein[i]->push_back(new std::vector<int>());

          std::vector<int>::iterator itx = m_rna_produce_protein[i]->at(prot_id)->begin();
          rna_id = 0;

          for (int xi = 0; xi < rna_produce_protein_list[old_i]->at(prot_id)->size(); xi++) {
            m_rna_produce_protein[i]->at(prot_id)->push_back(
                rna_produce_protein_list[old_i]->at(prot_id)->at(xi)
            );
            rna_id++;
          }
          prot_id++;
        }

        rna_id=0;
        std::vector<std::vector<double>*>::iterator itrnae =
            m_rna_influence_enhancing_coef_list[i]->begin();

        for (auto rna_enh_list : *(rna_influence_enhancing_coef_list[old_i])) {
          m_rna_influence_enhancing_coef_list[i]->push_back(new std::vector<double>());
          int prot_id = 0;

          std::vector<double>::iterator itx = m_rna_influence_enhancing_coef_list[i]->at(rna_id)->begin();

          for (int xi = 0; xi < rna_influence_enhancing_coef_list[old_i]->at(rna_id)->size(); xi++) {
            m_rna_influence_enhancing_coef_list[i]->at(rna_id)->push_back(
                rna_influence_enhancing_coef_list[old_i]->at(rna_id)->at(xi)
            );
            prot_id++;
          }
          rna_id++;
        }

        rna_id=0;
        std::vector<std::vector<double>*>::iterator itrnao =
            m_rna_influence_operating_coef_list[i]->begin();
        for (auto rna_enh_list : *(rna_influence_operating_coef_list[old_i])) {
          m_rna_influence_operating_coef_list[i]->push_back(new std::vector<double>());
          int prot_id = 0;
          std::vector<double>::iterator itx = m_rna_influence_operating_coef_list[i]->at(rna_id)->begin();

          for (int xi = 0; xi < rna_influence_operating_coef_list[old_i]->at(rna_id)->size(); xi++) {
            m_rna_influence_operating_coef_list[i]->at(rna_id)->push_back(
                rna_influence_operating_coef_list[old_i]->at(rna_id)->at(xi)
            );
            prot_id++;
          }
          rna_id++;
        }
      }
    }

    for (int i = 0; i < 1024; i++) {
      protein_concentration_list[i]->clear();
      rna_basal_concentration_list[i]->clear();
      rna_produce_protein_list[i]->clear();

      for (auto rna_enh_list : *(rna_influence_enhancing_coef_list[i])) {
        rna_enh_list->clear();
      }

      for (auto rna_enh_list : *(rna_influence_operating_coef_list[i])) {
        rna_enh_list->clear();
      }
    }
    protein_concentration_list.clear();
    rna_basal_concentration_list.clear();
    rna_produce_protein_list.clear();
    rna_influence_enhancing_coef_list.clear();
    rna_influence_operating_coef_list.clear();

    protein_concentration_list.swap(m_protein_concentration_list);
    rna_basal_concentration_list.swap(m_rna_basal_concentration_list);
    rna_produce_protein_list.swap(m_rna_produce_protein);
    rna_influence_enhancing_coef_list.swap(m_rna_influence_enhancing_coef_list);
    rna_influence_operating_coef_list.swap(m_rna_influence_operating_coef_list);
  }


  env_concentration_list.reserve(lifestep);
  std::vector<std::vector<double>*>::iterator itenv = env_concentration_list.begin();
  for (int i = 0; i < lifestep; i++) {
    env_concentration_list.insert(itenv+i,new std::vector<double>(nb_signal));
  }

  int lstep = 0;
  while ( getline (env_concentration,line) && lstep < lifestep) {
    std::vector<double> vect;

    std::stringstream ss(line);

    double i;

    while (ss >> i) {
      vect.push_back(i);
      if (ss.peek() == ',')
        ss.ignore();
    }
    lstep++;

    std::vector<double>::iterator itenv2 = env_concentration_list[(int)vect[0]]->begin();
    env_concentration_list[(int)vect[0]]->insert(itenv2+(int)vect[1],vect[2]);
  }

  printf("Lifestep %d -- Population %d (%ld)\n",lifestep,
         multiply_population*1024,protein_concentration_list.size());

  char* lookup_table_file_name = new char[100];

      sprintf( lookup_table_file_name, "lookup_table.ae" );

      gzFile lookup_table_file = gzopen( lookup_table_file_name, "r" );

      if ( lookup_table_file == Z_NULL )
      {
        printf( "ERROR : Could not read lookup table file %s\n", lookup_table_file_name );
        exit( EXIT_FAILURE );
      }

  double value;
      for (int i=0; i < LOOKUP_TABLE_SIZE; i++) {
        gzread( lookup_table_file, &value, sizeof(double));
        lookup_table_pow[i] = (double) value;
      }

  gzclose( lookup_table_file );

      delete[] lookup_table_file_name;

  if (odesolver == 0) {
    // Explicit Euler custom Solver
    if (execution_mode == 6) {
      transfer_to_tab(1024, lifestep);

      compute_openmp_gpu(lifestep, degradationstep, nb_signal, hill_shape_n,
                         hill_shape, degradation_rate);
    } else if (execution_mode == 7) {
      transfer_to_tab(1024, lifestep);

      compute_openmp_gpu_less_tasks(lifestep, degradationstep, nb_signal, hill_shape_n,
                         hill_shape, degradation_rate);
    } else if (execution_mode == 8) {
      tbb::task_group tgroup;

      for (int gen = 0; gen < nb_gen; gen++) {
        for (int i = 0; i < multiply_population * 1024; i++) {
          //t1 = high_resolution_clock::now();
          tgroup.run([=] { loop_one_user(lifestep,i); });
        }
        tgroup.wait();
      }
    } else if (execution_mode == 9) {
      tbb::task_group tgroup;

      for (int gen = 0; gen < nb_gen; gen++) {
        for (int i = 0; i < multiply_population * 1024; i++) {
          //t1 = high_resolution_clock::now();
          tgroup.run([=] { loop_one_user_parallel(lifestep,i); });
        }
        tgroup.wait();
      }
    }
  }
}

void loop_one_user(int lifestep, int i) {
  for (int lstep = 0; lstep < lifestep; lstep++) {
    update_env_indiv(lstep, i);
    solve_one_indiv_one_step(i);
  }
}

void loop_one_user_parallel(int lifestep, int i) {
  for (int lstep = 0; lstep < lifestep; lstep++) {
    update_env_indiv_parallel(lstep, i);
    solve_one_indiv_one_step_parallel(i);
  }
}

void compute_openmp_gpu_less_tasks(int lifestep, int degradationstep, int nb_signal,
                        double hill_shape,double hill_shape_n,double degradation_rate ) {
  int pop_size = 1024;


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

  double* delta = (double*) malloc(1024 * max_prot * sizeof(double));

  high_resolution_clock::time_point t_t1 = high_resolution_clock::now();

  {
    //for (int lstep = 0; lstep < lifestep; lstep++) {
//teams distribute parallel for collapse(2) schedule(static,1)
    for (int l = 0; l < lifestep; l++) {
#pragma omp parallel for
      for (int idx = 0; idx < 1024 * max_prot; idx++) {
        int indiv_id = idx / max_prot;
        int prot_id = idx % max_prot;
        //for (int indiv_id = 0; indiv_id < 1024; indiv_id++) {
        //for (int prot_id = nb_protein[indiv_id] - nb_signal;
        //     prot_id < nb_protein[indiv_id]; prot_id++) {
        protein_tab[indiv_id * max_prot + prot_id] =
            env_concentration_tab[lifestep * nb_signal + prot_id -
                                  nb_protein[indiv_id] + nb_signal];
        //}
      }

      for (int j = 0; j < degradationstep; j++) {
#pragma omp parallel for
        for (int idx = 0; idx < 1024 * max_prot; idx++) {
          int indiv_id = idx / max_prot;
          int prot_id = idx % max_prot;
          delta[indiv_id * max_prot + prot_id] = 0;
        }

//#pragma omp distribute parallel for collapse(3) schedule(static,1)
#pragma omp parallel for
        for (int idx = 0; idx < 1024 * max_prot; idx++) {
          int indiv_id = idx / max_prot;
          int prot_id = idx % max_prot;

          for (int j = 0; j < max_rna; j++) {

            int rna_id = rna_produce_protein_tab[indiv_id * max_prot +
                                                 prot_id * max_rna + j];
            if (rna_id != -1) {

              double enhancer_activity = 0;
              double operator_activity = 0;

#pragma omp simd
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
        }

#pragma omp parallel for
        for (int idx = 0; idx < 1024 * max_prot; idx++) {
          int indiv_id = idx / max_prot;
          int prot_id = idx % max_prot;

          delta[indiv_id * max_prot + prot_id] -=
              degradation_rate *
              protein_tab[indiv_id * max_prot + prot_id];
          delta[indiv_id * max_prot + prot_id] *=
              1 / (double) degradationstep;
        }

        // Apply the changes in concentrations we have just computed
#pragma omp parallel for
        for (int idx = 0; idx < 1024 * max_prot; idx++) {
          int indiv_id = idx / max_prot;
          int prot_id = idx % max_prot;

          protein_tab[indiv_id * max_prot +
                      prot_id] += delta[indiv_id * max_prot + prot_id];

        }
      }

    }
  }

  high_resolution_clock::time_point t_t2 = high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t_t2 - t_t1 ).count();

  std::cout<<"DURATION: "<<duration<<std::endl;
}

void compute_openmp_gpu(int lifestep, int degradationstep, int nb_signal,
                        double hill_shape,double hill_shape_n,double degradation_rate ) {
  int pop_size = 1024;
/*

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
*/
  //double* delta = (double*) malloc(1024 * max_prot * sizeof(double));
  std::map<int,std::vector<double>*> delta;

  for (int i = 0; i < pop_size; i++)
    delta[i] = new std::vector<double>(protein_concentration_list[i]->size());

  high_resolution_clock::time_point t_t1 = high_resolution_clock::now();

#pragma omp parallel
#pragma omp single
  {
    //for (int lstep = 0; lstep < lifestep; lstep++) {
//teams distribute parallel for collapse(2) schedule(static,1)
    for (int l = 0; l < lifestep; l++) {
      #pragma omp taskloop
      for (int idx = 0; idx < 1024 * max_prot; idx++) {
        int indiv_id = idx / max_prot;
        int prot_id = idx % max_prot;

        if ((prot_id >=
             protein_concentration_list[indiv_id]->size() - nb_signal) &&
            (prot_id < protein_concentration_list[indiv_id]->size()))
          protein_concentration_list[indiv_id]->at(prot_id) =
              env_concentration_list[lifestep]->at(
                  prot_id - protein_concentration_list[indiv_id]->size() +
                  nb_signal);
      }

      for (int j = 0; j < degradationstep; j++) {
        #pragma omp taskloop
        for (int idx = 0; idx < 1024 * max_prot; idx++) {
          int indiv_id = idx / max_prot;
          int prot_id = idx % max_prot;

          if (prot_id <
              protein_concentration_list[indiv_id]->size() - nb_signal)
            delta[indiv_id]->at(prot_id) = 0;
        }

//#pragma omp distribute parallel for collapse(3) schedule(static,1)

#pragma omp taskloop untied grainsize(max_rna)
        //num_tasks(65536)
        for (int idx = 0; idx < 1024 * max_prot * max_rna; idx++) {
          int indiv_id = idx / (max_prot * max_rna);
          int tmp2 = (idx % (max_rna * max_prot));

          int prot_id = tmp2 / max_rna;
          int j = tmp2 % max_rna;

          if ((prot_id <
               protein_concentration_list[indiv_id]->size() - nb_signal) &&
              (j < rna_produce_protein_list[indiv_id]->size())) {

            int rna_id = rna_produce_protein_list[indiv_id]->at(prot_id)->at(
                j);
            if (rna_id != -1) {

              double enhancer_activity = 0;
              double operator_activity = 0;

#pragma omp simd
              for (int i = 0; i <
                              rna_influence_enhancing_coef_list[indiv_id]->at(
                                  rna_id)->size(); i++) {

                enhancer_activity +=
                    rna_influence_enhancing_coef_list[indiv_id]->at(
                        rna_id)->at(
                        i)
                    * protein_concentration_list[indiv_id]->at(i);
                operator_activity +=
                    rna_influence_operating_coef_list[indiv_id]->at(
                        rna_id)->at(
                        i)
                    * protein_concentration_list[indiv_id]->at(i);
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
              delta[indiv_id]->at(prot_id) +=
                  rna_basal_concentration_list[indiv_id]->at(rna_id)
                  * (hill_shape
                     / (operator_activity_pow_n + hill_shape))
                  * (1 +
                     ((1 /
                       rna_basal_concentration_list[indiv_id]->at(rna_id)) -
                      1)
                     * (enhancer_activity_pow_n /
                        (enhancer_activity_pow_n + hill_shape)));

            }
          }
        }


        #pragma omp taskloop
        for (int idx = 0; idx < 1024 * max_prot; idx++) {
          int indiv_id = idx / max_prot;
          int prot_id = idx % max_prot;

          if (prot_id <
              protein_concentration_list[indiv_id]->size() - nb_signal) {
            delta[indiv_id]->at(prot_id) -=
                degradation_rate *
                protein_concentration_list[indiv_id]->at(prot_id);
            delta[indiv_id]->at(prot_id) *= 1 / (double) degradationstep;
          }
        }

        // Apply the changes in concentrations we have just computed
        #pragma omp taskloop
        for (int idx = 0; idx < 1024 * max_prot; idx++) {
          int indiv_id = idx / max_prot;
          int prot_id = idx % max_prot;

          if (prot_id <
              protein_concentration_list[indiv_id]->size() - nb_signal) {
            protein_concentration_list[indiv_id]->at(
                prot_id) += delta[indiv_id]->at(prot_id);

          }
        }
      }
    }

  }
  //}

  high_resolution_clock::time_point t_t2 = high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t_t2 - t_t1 ).count();

  std::cout<<"DURATION: "<<duration<<std::endl;
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


void update_env_indiv(int lifestep, int indiv_id) {
  for (int prot_id = protein_concentration_list[indiv_id]->size() - nb_signal;
       prot_id < protein_concentration_list[indiv_id]->size(); prot_id++) {
    protein_concentration_list[indiv_id]->at(prot_id) =
        env_concentration_list[lifestep]->at(prot_id-protein_concentration_list[indiv_id]->size()+nb_signal);
  }
}

void solve_one_indiv_one_step(int indiv_id) {
  std::vector<double> delta;

  for (int j = 0; j < degradationstep; j++) {
    delta.clear();

    for (int prot_id = 0;
         prot_id < protein_concentration_list[indiv_id]->size(); prot_id++) {

      if (prot_id < protein_concentration_list[indiv_id]->size() - nb_signal) {
        delta.insert(delta.begin() + prot_id, 0);

        for (int j = 0;
             j < rna_produce_protein_list[indiv_id]->at(prot_id)->size(); j++) {
          double enhancer_activity = 0;
          double operator_activity = 0;

          int rna_id = rna_produce_protein_list[indiv_id]->at(prot_id)->at(j);

          for (int i = 0; i <
                          rna_influence_enhancing_coef_list[indiv_id]->at(rna_id)->size(); i++) {

            enhancer_activity +=
                rna_influence_enhancing_coef_list[indiv_id]->at(rna_id)->at(i)
                * protein_concentration_list[indiv_id]->at(i);
            operator_activity +=
                rna_influence_operating_coef_list[indiv_id]->at(rna_id)->at(i)
                * protein_concentration_list[indiv_id]->at(i);
          }

          double enhancer_activity_pow_n = enhancer_activity == 0 ? 0 :
                                           pow(enhancer_activity, hill_shape_n);
          double operator_activity_pow_n = operator_activity == 0 ? 0 :
                                           pow(operator_activity, hill_shape_n);
          delta[prot_id] += rna_basal_concentration_list[indiv_id]->at(rna_id)
                            * (hill_shape
                               / (operator_activity_pow_n + hill_shape))
                            * (1 +
                               ((1 / rna_basal_concentration_list[indiv_id]->at(rna_id)) -
                                1)
                               * (enhancer_activity_pow_n /
                                  (enhancer_activity_pow_n + hill_shape)));
        }

        delta[prot_id] -=
            degradation_rate * protein_concentration_list[indiv_id]->at(prot_id);
        delta[prot_id] *= 1 / (double) degradationstep;
      }
    }

    // Apply the changes in concentrations we have just computed
    for (int prot_id = 0;
         prot_id < protein_concentration_list[indiv_id]->size(); prot_id++) {
      if (prot_id < protein_concentration_list[indiv_id]->size() - nb_signal) {
        protein_concentration_list[indiv_id]->at(prot_id) += delta[prot_id];
      }
    }
  }
}


void solve_one_indiv_one_step_parallel(int indiv_id) {
  std::vector<double> delta;

  for (int j = 0; j < degradationstep; j++) {
    delta.clear();

    parallel_for(size_t(0), protein_concentration_list[indiv_id]->size(), size_t(1) , [&](size_t prot_id) {one_prot(delta,indiv_id,prot_id);});

    parallel_for(size_t(0), protein_concentration_list[indiv_id]->size(), size_t(1) , [&](size_t prot_id) { update_one_prot(delta,indiv_id,prot_id);});
  }
}

void update_env_indiv_parallel(int lifestep, int indiv_id) {
  parallel_for(protein_concentration_list[indiv_id]->size() - nb_signal,
               protein_concentration_list[indiv_id]->size(), size_t(1) , [&](size_t prot_id) {
          update_one_env(lifestep,indiv_id,prot_id);
      });
}

void update_one_env(int lifestep, int indiv_id, int prot_id) {
  protein_concentration_list[indiv_id]->at(prot_id) =
      env_concentration_list[lifestep]->at(prot_id-protein_concentration_list[indiv_id]->size()+nb_signal);

}
void update_one_prot(std::vector<double> &delta, int indiv_id, int prot_id) {
  if (prot_id < protein_concentration_list[indiv_id]->size() - nb_signal) {
    protein_concentration_list[indiv_id]->at(prot_id) += delta[prot_id];
  }
}

void one_prot(std::vector<double> &delta, int indiv_id, int prot_id) {
  if (prot_id < protein_concentration_list[indiv_id]->size() - nb_signal) {
    delta.insert(delta.begin() + prot_id, 0);

    for (int j = 0;
         j < rna_produce_protein_list[indiv_id]->at(prot_id)->size(); j++) {
      double enhancer_activity = 0;
      double operator_activity = 0;

      int rna_id = rna_produce_protein_list[indiv_id]->at(prot_id)->at(j);

      for (int i = 0; i <
                      rna_influence_enhancing_coef_list[indiv_id]->at(rna_id)->size(); i++) {

        enhancer_activity +=
            rna_influence_enhancing_coef_list[indiv_id]->at(rna_id)->at(i)
            * protein_concentration_list[indiv_id]->at(i);
        operator_activity +=
            rna_influence_operating_coef_list[indiv_id]->at(rna_id)->at(i)
            * protein_concentration_list[indiv_id]->at(i);
      }

      double enhancer_activity_pow_n = enhancer_activity == 0 ? 0 :
                                       pow(enhancer_activity, hill_shape_n);
      double operator_activity_pow_n = operator_activity == 0 ? 0 :
                                       pow(operator_activity, hill_shape_n);
      delta[prot_id] += rna_basal_concentration_list[indiv_id]->at(rna_id)
                        * (hill_shape
                           / (operator_activity_pow_n + hill_shape))
                        * (1 +
                           ((1 / rna_basal_concentration_list[indiv_id]->at(rna_id)) -
                            1)
                           * (enhancer_activity_pow_n /
                              (enhancer_activity_pow_n + hill_shape)));
    }

    delta[prot_id] -=
        degradation_rate * protein_concentration_list[indiv_id]->at(prot_id);
    delta[prot_id] *= 1 / (double) degradationstep;
  }
}
