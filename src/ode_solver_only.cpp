//
// Created by arrouan on 14/11/16.
//

//
// Created by arrouan on 14/11/16.
//
// =================================================================
//                              Libraries
// =================================================================
#include <errno.h>
#include <inttypes.h>
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include <sys/stat.h>
#include <omp.h>
#include <list>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include<chrono>
// =================================================================
//                            Project Files
// =================================================================
#include "ode_solver_only.h"
#include "ae_logger.h"
using namespace std::chrono;


unordered_map<int,unordered_multiset<string>> ae_logger::logMap;
string ae_logger::logFile = "logger_csv.log";
mutex ae_logger::loggerMtx;

#include "aevol.h"

using namespace aevol;




int main(int argc, char* argv[]) {
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
        Utils::PrintAevolVersion();
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

  AeTime::set_time(0);

  printf("Lifestep %d -- Population %d (%ld)\n",lifestep,
         multiply_population*1024,protein_concentration_list.size());
  high_resolution_clock::time_point t_t1 = high_resolution_clock::now();
  high_resolution_clock::time_point t_t2,t1,t2;
  if (odesolver == 0) {
    // Explicit Euler custom Solver
    if (execution_mode == 0) {
      // Sequential solver
      for (int i = 0; i < multiply_population*1024; i++) {
        for (int lstep = 0; lstep < lifestep; lstep++) {
          update_env_indiv(lstep,i);
          solve_one_indiv_one_step(i);

        }
      }
    } else if (execution_mode == 1) {
      // Sequential merge

      for (int i = 0; i < multiply_population*1024; i+=merge) {
        for (int lstep = 0; lstep < lifestep; lstep++) {
          update_env_list_indiv(lstep,i,i+merge-1);
          solve_list_indiv_one_step(i,i+merge-1);

        }
      }
    } else if (execution_mode == 2) {
      // Parallel

      #pragma omp parallel for num_threads(parallelthread) schedule(dynamic)
      for (int i = 0; i < multiply_population*1024; i++) {
        //t1 = high_resolution_clock::now();
        for (int lstep = 0; lstep < lifestep; lstep++) {
          update_env_indiv(lstep,i);
          solve_one_indiv_one_step(i);
        }
/*
        t2 = high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
        printf("%d %d %d\n",t1,t2,duration);

        ae_logger::addLog(SELECTION,duration);
        ae_logger::flush(AeTime::time());
        AeTime::plusplus();*/
      }
    } else if (execution_mode == 3) {
      // Parallel merge
      #pragma omp parallel for num_threads(parallelthread) schedule(dynamic)
      for (int i = 0; i < multiply_population*1024; i+=merge) {
        //t1 = high_resolution_clock::now();

        for (int lstep = 0; lstep < lifestep; lstep++) {
          update_env_list_indiv(lstep,i,i+merge-1);
          solve_list_indiv_one_step(i,i+merge-1);

        }
/*
        t2 = high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
        printf("%d %d %d\n",t1,t2,duration);

        ae_logger::addLog(SELECTION,duration);
        ae_logger::flush(AeTime::time());
        AeTime::plusplus();*/
      }
    } else if (execution_mode == 4) {
      int max_protein = transfert_data_to_gpu(1024*multiply_population,lifestep);

      process_delta<<<1024*multiply_population,max_protein>>>(nb_signal,degradationstep,degradation_rate,
          dev_rna_produce_protein_array, dev_nb_rna_produce_protein, dev_nb_rna_produce, dev_protein_concentration_array,
          dev_rna_basal_concentration_array, dev_nb_protein_array, dev_nb_rna_array,
          dev_rna_influence_enhancing_coef_array, dev_rna_influence_operating_coef_array,
          dev_nb_rna_influence_enhancing_coef, dev_nb_rna_influence_operating_coef,
          dev_env_concentration_array,hill_shape,hill_shape_n);
    }
  }
  t_t2 = high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t_t2 - t_t1 ).count();
  ae_logger::addLog(TOTAL,duration);
  ae_logger::flush(AeTime::time());

  std::ofstream finished("finished.lock");
  finished<<"ok"<<std::endl;
  finished.close();
}

void update_env_list_indiv(int lifestep, int start_indiv_id, int end_indiv_id) {
  for (int indiv_id = start_indiv_id; indiv_id <= end_indiv_id; indiv_id++) {
    for (int prot_id = protein_concentration_list[indiv_id]->size() - nb_signal;
         prot_id < protein_concentration_list[indiv_id]->size(); prot_id++) {

      protein_concentration_list[indiv_id]->at(prot_id) =
          env_concentration_list[lifestep]->at(prot_id-protein_concentration_list[indiv_id]->size()+nb_signal);
    }
  }
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


void solve_list_indiv_one_step(int start_indiv_id, int end_indiv_id) {
  for (int indiv_id = start_indiv_id; indiv_id <= end_indiv_id; indiv_id++) {
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
      rna_produce_protein_array[i][rna_id] = (int*)malloc(
          rna_produce_protein_list[i]->at(rna_id)->size()*sizeof(int));
      nb_rna_produce_protein[i][prot_id] =
          rna_produce_protein_list[i]->at(rna_id)->size();
      for (int ix = 0; ix < rna_produce_protein_list[i]->at(prot_id)->size(); ix++) {
        rna_produce_protein_array[i][rna_id][prot_id] =
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

  cudaMemcpy(dev_nb_rna_influence_enhacing_coef_l1,
             nb_rna_influence_enhacing_coef_l1, pop_size * sizeof(int), cudaMemcpyHostToDevice);
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
