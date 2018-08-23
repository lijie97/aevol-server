//
// Created by arrouan on 27/07/17.
//

#include "SIMD_Individual.h"
#include "Dna_SIMD.h"
#include "DnaMutator.h"
#include "HybridFuzzy.h"
#include "Stats_SIMD.h"
#include "ExpManager.h"


#include <omp.h>
namespace aevol {

#ifndef WITH_STANDALONE_SIMD
     bool SIMD_Individual::standalone_simd = false;
#else
     bool SIMD_Individual::standalone_simd = true;
#endif

SIMD_Individual::SIMD_Individual(ExpManager* exp_m) {

  printf("Create SIMD Controller\n");
    standalone_ = standalone_simd;
  exp_m_ = exp_m;

  nb_indivs_ = exp_m_->nb_indivs();

  internal_simd_struct = new Internal_SIMD_Struct* [exp_m_->nb_indivs()];
  prev_internal_simd_struct = new Internal_SIMD_Struct* [exp_m_->nb_indivs()];

  next_generation_reproducer_ = new int32_t[exp_m_->nb_indivs()];

  for (int indiv_id = 0; indiv_id < exp_m_->nb_indivs(); indiv_id++) {
    int x = indiv_id / exp_m_->world()->height();
    int y = indiv_id % exp_m_->world()->height();

    internal_simd_struct[indiv_id] = new Internal_SIMD_Struct(exp_m);
    //printf("DNA SIMD %d\n",indiv_id);
    internal_simd_struct[indiv_id]->dna_ = new Dna_SIMD(exp_m->world()->grid(x,y)->individual()->genetic_unit(0).dna());
    //printf("SetIndiv\n");
    internal_simd_struct[indiv_id]->indiv_id = indiv_id;
    internal_simd_struct[indiv_id]->parent_id = indiv_id;
    //printf("Set Prev Indiv\n");
    prev_internal_simd_struct[indiv_id] = internal_simd_struct[indiv_id];
      internal_simd_struct[indiv_id]->global_id = AeTime::time()*1024+indiv_id;
      //printf("Building indiv %d %d\n",indiv_id,internal_simd_struct[indiv_id]->dna_->length_);
  }

  dna_size = new int[exp_m_->nb_indivs()];
  int x, y;
  for (int indiv_id = 0; indiv_id < exp_m_->nb_indivs(); indiv_id++) {
    x = indiv_id / exp_m_->world()->height();
    y = indiv_id % exp_m_->world()->height();


    dna_size[indiv_id] = internal_simd_struct[indiv_id]->dna_->length();
  }

  target = new double[300];
  for (int i = 0; i < 300; i++) {
      double tmp =((HybridFuzzy*) exp_m->world()->phenotypic_target_handler()->phenotypic_target().fuzzy())->points()[i];


    target[i] = tmp;
      //printf("AT %d value is %e -- %e\n",i,tmp,target[i]);
  }

    /*for (int ip = 0; ip < 300; ip++) {
        printf("PH[%d] = %e -- %e\n",ip,

              // ((HybridFuzzy*)exp_m_->world()->grid(0, 0)->phenotypic_target().fuzzy())[ip],
               ((HybridFuzzy*) exp_m->world()->phenotypic_target_handler()->phenotypic_target().fuzzy())->points()[ip],
               target[ip]
        );
    }*/

  //build_phenotypic_target(exp_m->world()->phenotypic_target_handler());



  /*for (int i = 0; i < 300; i++) {
    printf("%i  -- %f / %lf\n", i, target[i],
           ((HybridFuzzy*) exp_m->world()->phenotypic_target_handler()->phenotypic_target().fuzzy())->points()[i]);
  }*/
}

void SIMD_Individual::clear_struct_before_next_step() {
/*  for (int indiv_id = 0; indiv_id < (int) exp_m_->nb_indivs(); indiv_id++) {
    delete prev_internal_simd_struct[indiv_id];
  }

  delete[] prev_internal_simd_struct; */

  // Empty all struct except promoters
  for (int indiv_id = 0; indiv_id < exp_m_->nb_indivs(); indiv_id++) {
    internal_simd_struct[indiv_id]->proteins.clear();
    internal_simd_struct[indiv_id]->rnas.clear();
    internal_simd_struct[indiv_id]->terminator_lead.clear();
    internal_simd_struct[indiv_id]->terminator_lag.clear();
  }
}

void SIMD_Individual::selection() {
  int32_t selection_scope_x = exp_m_->sel()->selection_scope_x();
  int32_t selection_scope_y = exp_m_->sel()->selection_scope_y();

  int32_t grid_width = exp_m_->grid_width();
  int32_t grid_height = exp_m_->grid_height();

  int16_t neighborhood_size = selection_scope_x * selection_scope_y;

  for (int indiv_id = 0; indiv_id < exp_m_->nb_indivs(); indiv_id++) {
    double *  local_fit_array   = new double[neighborhood_size];
    double *  local_meta_array   = new double[neighborhood_size];
    double *  probs             = new double[neighborhood_size];
    int16_t   count             = 0;
    double    sum_local_fit     = 0.0;

    int32_t x = indiv_id / grid_height;
    int32_t y = indiv_id % grid_height;

    int cur_x,cur_y;

    for (int8_t i = -1 ; i < selection_scope_x-1 ; i++) {
      for (int8_t j = -1; j < selection_scope_y - 1; j++) {
        cur_x = (x + i + grid_width)  % grid_width;
        cur_y = (y + j + grid_height) % grid_height;

        local_fit_array[count]  = prev_internal_simd_struct[cur_x*grid_height+cur_y]->fitness;
        local_meta_array[count]  = prev_internal_simd_struct[cur_x*grid_height+cur_y]->metaerror;
        sum_local_fit += local_fit_array[count];
        count++;
      }
    }

    for(int16_t i = 0 ; i < neighborhood_size ; i++) {
      probs[i] = local_fit_array[i]/sum_local_fit;
      //printf("%d -- prob[%d] : %e : fitness %e (%f) sum %e\n",indiv_id,i,probs[i],
      //       local_fit_array[i],local_meta_array[i],sum_local_fit);
    }

    int16_t found_org = exp_m_->world()->grid(x,y)->reprod_prng_->roulette_random(probs, neighborhood_size);

    int16_t x_offset = (found_org / selection_scope_x) - 1;
    int16_t y_offset = (found_org % selection_scope_y) - 1;

    delete [] local_fit_array;
    delete [] local_meta_array;
    delete [] probs;

    exp_m_->simd_individual->internal_simd_struct[indiv_id] =
        new Internal_SIMD_Struct(exp_m_,prev_internal_simd_struct
                [((x+x_offset+grid_width)  % grid_width)*grid_height+
                    ((y+y_offset+grid_height) % grid_height)],false);

    exp_m_->simd_individual->internal_simd_struct[indiv_id]->indiv_id = indiv_id;
    exp_m_->simd_individual->internal_simd_struct[indiv_id]->parent_id =
        ((x+x_offset+grid_width)  % grid_width)*grid_height+
        ((y+y_offset+grid_height) % grid_height);

    //printf("New indiv %d parent %d\n",exp_m_->simd_individual->internal_simd_struct[indiv_id]->indiv_id,
    //       exp_m_->simd_individual->internal_simd_struct[indiv_id]->parent_id);
  }
}



    void SIMD_Individual::selection(int indiv_id) {
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

            int32_t x = indiv_id / grid_height;
            int32_t y = indiv_id % grid_height;

            int cur_x,cur_y;

            for (int8_t i = -1 ; i < selection_scope_x-1 ; i++) {
                for (int8_t j = -1; j < selection_scope_y - 1; j++) {
                    cur_x = (x + i + grid_width)  % grid_width;
                    cur_y = (y + j + grid_height) % grid_height;

                    local_fit_array[count]  = prev_internal_simd_struct[cur_x*grid_height+cur_y]->fitness;
                    local_meta_array[count]  = prev_internal_simd_struct[cur_x*grid_height+cur_y]->metaerror;
                    sum_local_fit += local_fit_array[count];
                    count++;
                }
            }

            for(int16_t i = 0 ; i < neighborhood_size ; i++) {
                probs[i] = local_fit_array[i]/sum_local_fit;
                //printf("%d -- prob[%d] : %e : fitness %e (%f) sum %e\n",indiv_id,i,probs[i],
                //       local_fit_array[i],local_meta_array[i],sum_local_fit);
            }

            int16_t found_org = exp_m_->world()->grid(x,y)->reprod_prng_->roulette_random(probs, neighborhood_size);

            int16_t x_offset = (found_org / selection_scope_x) - 1;
            int16_t y_offset = (found_org % selection_scope_y) - 1;

            delete [] local_fit_array;
            delete [] local_meta_array;
            delete [] probs;

            next_generation_reproducer_[indiv_id] = ((x+x_offset+grid_width)  % grid_width)*grid_height+
                                          ((y+y_offset+grid_height) % grid_height);
            dna_size[indiv_id] = prev_internal_simd_struct[cur_x*grid_height+cur_y]->dna_->length();



             /*exp_m_->simd_individual->internal_simd_struct[indiv_id] =
                    new Internal_SIMD_Struct(exp_m_,prev_internal_simd_struct
                    [((x+x_offset+grid_width)  % grid_width)*grid_height+
                     ((y+y_offset+grid_height) % grid_height)],false);


             exp_m_->simd_individual->internal_simd_struct[indiv_id]->indiv_id = indiv_id;

             exp_m_->simd_individual->internal_simd_struct[indiv_id]->parent_id =
                    ((x+x_offset+grid_width)  % grid_width)*grid_height+
                    ((y+y_offset+grid_height) % grid_height);*/



            //printf("New indiv %d parent %d\n",exp_m_->simd_individual->internal_simd_struct[indiv_id]->indiv_id,
            //       exp_m_->simd_individual->internal_simd_struct[indiv_id]->parent_id);
    }



    void SIMD_Individual::check_selection(int indiv_id) {
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


                local_fit_array[count] = prev_internal_simd_struct[cur_x * grid_height + cur_y]->fitness;

                local_meta_array[count] = prev_internal_simd_struct[cur_x * grid_height + cur_y]->metaerror;

                //printf("Local fit %e %d %d %d -> %e\n",local_fit_array[count],cur_x,cur_y,cur_x*grid_height+cur_y,
                //       prev_internal_simd_struct[cur_x*grid_height+cur_y]->fitness);


                /*printf("Metaerror : %e -- %e\n", local_meta_array[count],
                       internal_simd_struct[cur_x * grid_height + cur_y]->metaerror);*/
                if (local_meta_array[count] != prev_internal_simd_struct[cur_x * grid_height + cur_y]->metaerror) {
                    printf("NONONNONONONN\n");
                    exit(-1);
                }
                sum_local_fit += local_fit_array[count];


                indiv_index[count] = cur_x * grid_height + cur_y;
                /*if (indiv_id == 0)
                    printf("SIMD Local SUM Fit %e -- Fitness %e (CPU %e)\n",sum_local_fit,local_fit_array[count],exp_m_->world()->grid(x,y)->local_fit_array[count]);*/
                count++;
            }
        }

        for(int16_t i = 0 ; i < neighborhood_size ; i++) {
            probs[i] = local_fit_array[i]/sum_local_fit;

                /*printf("Local fit X %e %d %d %d -> %e\n",local_fit_array[i],cur_x,cur_y,cur_x*grid_height+cur_y,
                       prev_internal_simd_struct[cur_x*grid_height+cur_y]->fitness);*/
            //printf("%d -- prob[%d] : %e : fitness %e (%f) sum %e\n",indiv_id,i,probs[i],
            //       local_fit_array[i],local_meta_array[i],sum_local_fit);
        }

        //printf("SIMD PRNG\n");
        /*        bool verbose = false;
        if (indiv_id == 2)  {
            printf("SIMD PRNG\n");
            verbose = true;
        }*/
        int16_t found_org = exp_m_->world()->grid(x,y)->reprod_prng_simd_->roulette_random(probs, neighborhood_size);

        int16_t x_offset = (found_org / selection_scope_x) - 1;
        int16_t y_offset = (found_org % selection_scope_y) - 1;


        int found_id = ((x+x_offset+grid_width)  % grid_width)*grid_height+
                       ((y+y_offset+grid_height) % grid_height);


            if (found_id != next_generation_reproducer_[indiv_id]) {


                printf("For individual %d: Selection is diff SIMD %d CPU %d (Meta error %f -- %f || Fitness %e -- %e) \n",
                       indiv_id, found_id, next_generation_reproducer_[indiv_id],
                       prev_internal_simd_struct[next_generation_reproducer_[indiv_id]]->metaerror,
                       prev_internal_simd_struct[found_id]->metaerror,
                       prev_internal_simd_struct[next_generation_reproducer_[indiv_id]]->fitness,
                       prev_internal_simd_struct[found_id]->fitness);



                for (int i = 0; i < neighborhood_size; i++) {
/*                    if (i==0) {
                        for (int ip = 0; ip < 300; ip++) {
                            printf("PH[%d] = %e -- %e (%e) || %e -- %e\n",ip,
                                   exp_m_->world()->grid(x, y)->loc_phenotype[ip],
                                   prev_internal_simd_struct[indiv_index[i]]->phenotype[ip],
                                   exp_m_->world()->grid(x, y)->loc_phenotype[ip] -
                                   prev_internal_simd_struct[indiv_index[i]]->phenotype[ip],
                                   ((HybridFuzzy*) exp_m_->world()->phenotypic_target_handler()->phenotypic_target().fuzzy())->points()[ip],
                                   target[ip]
                                );
                        }

                    }*/

                    if (i==5) {
                        int v_x = indiv_index[i] / grid_height;
                        int v_y = indiv_index[i] % grid_height;

                        printf(
                                "ERROR -- Individual %d (%d,%d): Metaerror (CPU/GPU) : %e/%e || Fitness (CPU/GPU) : %e/%e \n",
                                indiv_index[i],v_x,v_y,
                                exp_m_->world()->grid(v_x, v_y)->individual()->dist_to_target_by_feature(
                                        METABOLISM),
                                prev_internal_simd_struct[indiv_index[i]]->metaerror,
                                exp_m_->world()->grid(v_x, v_y)->individual()->fitness(),
                                prev_internal_simd_struct[indiv_index[i]]->fitness);

                        printf(
                                "Nb RNA SIMD/CPU %ld/%ld Protein %ld/%ld Metaerror %f/%f Fitness %e/%e DNA Size %d/%d\n",
                                prev_internal_simd_struct[indiv_index[i]]->rnas.size(),
                                exp_m_->world()->grid(v_x, v_y)->individual()->rna_list().size(),
                                prev_internal_simd_struct[indiv_index[i]]->proteins.size(),
                                exp_m_->world()->grid(v_x, v_y)->individual()->protein_list().size(),
                                prev_internal_simd_struct[indiv_index[i]]->metaerror,
                                exp_m_->world()->grid(v_x, v_y)->individual()->dist_to_target_by_feature(
                                        METABOLISM), internal_simd_struct[indiv_index[i]]->fitness,
                                exp_m_->world()->grid(v_x, v_y)->individual()->fitness(), dna_size[i],
                                exp_m_->world()->grid(v_x, v_y)->individual()->genetic_unit(
                                        0).seq_length());

                        check_individual(indiv_index[i],v_x,v_y);

                    }


                    printf("%d -- %d / %d-- (Probs %e %e -- Fit Array %e %e -- Sum Fit %e %e -- Metaerror %e %e (%e) -- LFIT %e %e)\n",i,
                           exp_m_->world()->grid(x, y)->indiv_index[i], indiv_index[i],
                           exp_m_->world()->grid(x,y)->probs[i],probs[i],
                           exp_m_->world()->grid(x,y)->local_fit_array[i],local_fit_array[i],
                           exp_m_->world()->grid(x,y)->sum_local_fit,sum_local_fit,
                           exp_m_->world()->grid(x,y)->local_meta_array[i],local_meta_array[i],
                           exp_m_->world()->grid(x,y)->local_meta_array[i]-local_meta_array[i],
                           exp(-exp_m_->selection_pressure()*exp_m_->world()->grid(x,y)->local_meta_array[i]),
                           exp(-exp_m_->selection_pressure()*local_meta_array[i]));
                }

                exit(-44);
            }



        delete [] local_fit_array;
        delete [] local_meta_array;
        delete [] probs;


/*        delete [] exp_m_->world()->grid(x,y)->probs;
        delete [] exp_m_->world()->grid(x,y)->local_fit_array;*/

    }


    void SIMD_Individual::do_mutation() {
//  for (int indiv_id = 0; indiv_id < exp_m_->nb_indivs(); indiv_id++) {
//    if (internal_simd_struct[indiv_id]->indiv_id == 212) {
//      printf("APPLY_MUTATION : Leading promoters lists : ");
//      for (auto it : internal_simd_struct[indiv_id]->leading_prom_pos) {
//        printf("%d (%d) || ", it.first, it.second);
//      }
//      printf("\n");
//
//      printf("APPLY_MUTATION : Leading promoters lists (promoters): ");
//      for (auto it : internal_simd_struct[indiv_id]->promoters) {
//        printf("%d (%d) -- ", it.second->pos, it.first);
//      }
//
//      printf("\n");
//    }
//  }
  if (standalone_) {

    for (int indiv_id = 0; indiv_id < exp_m_->nb_indivs(); indiv_id++) {
      int x = indiv_id / exp_m_->world()->height();
      int y = indiv_id % exp_m_->world()->height();
      delete exp_m_->dna_mutator_array_[indiv_id];

      exp_m_->dna_mutator_array_[indiv_id] = new DnaMutator(
          exp_m_->world()->grid(x, y)->mut_prng(),
          internal_simd_struct[indiv_id]->dna_->length(),
          exp_m_->exp_s()->mut_params()->duplication_rate(),
          exp_m_->exp_s()->mut_params()->deletion_rate(),
          exp_m_->exp_s()->mut_params()->translocation_rate(),
          exp_m_->exp_s()->mut_params()->inversion_rate(),
          exp_m_->exp_s()->mut_params()->point_mutation_rate(),
          exp_m_->exp_s()->mut_params()->small_insertion_rate(),
          exp_m_->exp_s()->mut_params()->small_deletion_rate(),
          exp_m_->exp_s()->mut_params()->max_indel_size(),
          exp_m_->exp_s()->min_genome_length(),
          exp_m_->exp_s()->max_genome_length());
      exp_m_->dna_mutator_array_[indiv_id]->generate_mutations();
    }
  }
  nb_clones_ = 0;
  for (int indiv_id = 0; indiv_id < exp_m_->nb_indivs(); indiv_id++) {

//    dna_size[indiv_id] = internal_simd_struct[indiv_id]->dna_->length();
//    if (indiv_id == 43)
//      printf("DNA BEFORE SIZE of %d is %d (%d)\n",indiv_id,dna_size[indiv_id],internal_simd_struct[indiv_id]->dna_->length());

      if (exp_m_->dna_mutator_array_[indiv_id]->hasMutate()) {
#ifdef WITH_BITSET
        internal_simd_struct[indiv_id]->dna_->bitset_ =
            new BitSet_SIMD(prev_internal_simd_struct
              [internal_simd_struct[indiv_id]->parent_id]->dna_->bitset_);
#endif

 //       int x = indiv_id / exp_m_->world()->height();
 //       int y = indiv_id % exp_m_->world()->height();
//        printf("Before mutation %d (%d %d)-- %d %d %d\n",indiv_id,x,y,
//               internal_simd_struct[indiv_id]->dna_->length(),dna_size[indiv_id],
//               exp_m_->world()->grid(x, y)->individual()->genetic_unit(0).seq_length());

        /*if (indiv_id == 49) {
          printf("Size before mutation %d\n",internal_simd_struct[indiv_id]->dna_->length());
        }*/
        if (standalone_)
          internal_simd_struct[indiv_id]->dna_->apply_mutations_standalone();
        else
          internal_simd_struct[indiv_id]->dna_->apply_mutations();
      } else {
        nb_clones_++;
        int32_t parent_id = internal_simd_struct[indiv_id]->parent_id;
        delete internal_simd_struct[indiv_id];
        internal_simd_struct[indiv_id] = prev_internal_simd_struct[parent_id];
        internal_simd_struct[indiv_id]->usage_count_++;
      }

/*
    if (indiv_id == 49) {
      printf("Size after mutation %d\n",internal_simd_struct[indiv_id]->dna_->length());
    }
*/
//    if (indiv_id == 43)
//      printf("DNA AFTER SIZE of %d is %d (%d)\n",indiv_id,dna_size[indiv_id],internal_simd_struct[indiv_id]->dna_->length());

/*
      if (internal_simd_struct[indiv_id]->dna_->mutation_list.size() > 0) {
        int x = indiv_id / exp_m_->world()->height();
        int y = indiv_id % exp_m_->world()->height();

        bool same = true;
        if (internal_simd_struct[indiv_id]->dna_->length() !=
              exp_m_->world()->grid(x, y)->individual()->genetic_unit(0).seq_length())
          same=false;

        if (same)
          for (int i = 0; i < internal_simd_struct[indiv_id]->dna_->length(); i++) {
            if (internal_simd_struct[indiv_id]->dna_->data_[i] != exp_m_->world()->grid(x, y)->individual()->genetic_unit(0).dna()->data()[i]) {
              same = false;
              break;
            }
          }

        // TODO Add promoters verification
        for (auto rna : exp_m_->world()->grid(x, y)->individual()->rna_list()) {
          bool found = false;
          if (rna->strand() == LEADING) {
            for (auto it : internal_simd_struct[indiv_id]->leading_prom_pos) {
              if (it.first == rna->promoter_pos()) {
                found = true;
                break;
              }
            }
          } else {
            for (auto it : internal_simd_struct[indiv_id]->lagging_prom_pos) {
              if (it.first == rna->promoter_pos()) {
                found = true;
                break;
              }
            }
          }

          if (!found) {
            printf("Indiv %d -- Promoter at position %d not found\n",indiv_id,rna->promoter_pos());
            same = false;
            break;
          }
        }

        if (!same) {
          printf("Incorrect mutations replay %d (%d %d)\n",indiv_id,x,y);
          //printf("Aevol (%d) %s\n",exp_m_->world()->grid(x, y)->individual()->genetic_unit(0).dna()->length(),
          //       exp_m_->world()->grid(x, y)->individual()->genetic_unit(0).dna()->data());
          //printf("SIMD (%d) %s\n",internal_simd_struct[indiv_id]->dna_->length(),
          //       internal_simd_struct[indiv_id]->dna_->data_);

          printf("Leading promoters lists : ");
          for (auto it : internal_simd_struct[indiv_id]->leading_prom_pos) {
            printf("%d (%d) || ", it.first, it.second);
          }
          printf("\n");

          printf("Lagging promoters lists : ");
          for (auto it : internal_simd_struct[indiv_id]->lagging_prom_pos) {
            printf("%d (%d) || ", it.first, it.second);
          }
          printf("\n");

          printf("Promoters lists (promoters): ");
          for (auto it : internal_simd_struct[indiv_id]->promoters) {
            printf("%d (%d) -- ", it.second->pos, it.first);
          }
          printf("\n");

          printf("RNA Promoters lists : ");
          for (auto rna : exp_m_->world()->grid(x, y)->individual()->rna_list()) {
            printf("%d ",rna->promoter_pos());
          }

          printf("\n");

          exit(40);
        }

      }*/

    dna_size[indiv_id] = internal_simd_struct[indiv_id]->dna_->length();
  }
  //exit(12);
}



    void SIMD_Individual::do_mutation(int indiv_id) {
//  for (int indiv_id = 0; indiv_id < exp_m_->nb_indivs(); indiv_id++) {
//    if (internal_simd_struct[indiv_id]->indiv_id == 212) {
//      printf("APPLY_MUTATION : Leading promoters lists : ");
//      for (auto it : internal_simd_struct[indiv_id]->leading_prom_pos) {
//        printf("%d (%d) || ", it.first, it.second);
//      }
//      printf("\n");
//
//      printf("APPLY_MUTATION : Leading promoters lists (promoters): ");
//      for (auto it : internal_simd_struct[indiv_id]->promoters) {
//        printf("%d (%d) -- ", it.second->pos, it.first);
//      }
//
//      printf("\n");
//    }
//  }
        if (standalone_) {

                int x = indiv_id / exp_m_->world()->height();
                int y = indiv_id % exp_m_->world()->height();
                delete exp_m_->dna_mutator_array_[indiv_id];

                exp_m_->dna_mutator_array_[indiv_id] = new DnaMutator(
                        exp_m_->world()->grid(x, y)->mut_prng(),
                        //prev_internal_simd_struct[next_generation_reproducer_[indiv_id]]->dna_->length(),
                        dna_size[indiv_id],
                        exp_m_->exp_s()->mut_params()->duplication_rate(),
                        exp_m_->exp_s()->mut_params()->deletion_rate(),
                        exp_m_->exp_s()->mut_params()->translocation_rate(),
                        exp_m_->exp_s()->mut_params()->inversion_rate(),
                        exp_m_->exp_s()->mut_params()->point_mutation_rate(),
                        exp_m_->exp_s()->mut_params()->small_insertion_rate(),
                        exp_m_->exp_s()->mut_params()->small_deletion_rate(),
                        exp_m_->exp_s()->mut_params()->max_indel_size(),
                        exp_m_->exp_s()->min_genome_length(),
                        exp_m_->exp_s()->max_genome_length());
                exp_m_->dna_mutator_array_[indiv_id]->generate_mutations();
        }

//    dna_size[indiv_id] = internal_simd_struct[indiv_id]->dna_->length();
//    if (indiv_id == 43)
//      printf("DNA BEFORE SIZE of %d is %d (%d)\n",indiv_id,dna_size[indiv_id],internal_simd_struct[indiv_id]->dna_->length());

            if (exp_m_->dna_mutator_array_[indiv_id]->hasMutate()) {
/*             if (internal_simd_struct[indiv_id] != nullptr) {
                printf("Warning will replace something that is not null !!!\n");

             }*/

             internal_simd_struct[indiv_id] =
                        new Internal_SIMD_Struct(exp_m_,prev_internal_simd_struct
                        [next_generation_reproducer_[indiv_id]],false);

                internal_simd_struct[indiv_id]->global_id = AeTime::time()*1024+indiv_id;
                internal_simd_struct[indiv_id]->indiv_id = indiv_id;
                internal_simd_struct[indiv_id]->parent_id =
                        next_generation_reproducer_[indiv_id];

                if (standalone_) {
#pragma omp critical
                    {
                        int x = indiv_id / exp_m_->world()->height();
                        int y = indiv_id % exp_m_->world()->height();
                        NewIndivEvent *eindiv = new NewIndivEvent(internal_simd_struct[indiv_id],
                                                                  prev_internal_simd_struct[next_generation_reproducer_[indiv_id]],
                                                                  x, y,indiv_id,next_generation_reproducer_[indiv_id]);
                        notifyObservers(NEW_INDIV, eindiv);
                        delete eindiv;
                    }
                }

#ifdef WITH_BITSET
                internal_simd_struct[indiv_id]->dna_->bitset_ =
            new BitSet_SIMD(prev_internal_simd_struct
              [internal_simd_struct[indiv_id]->parent_id]->dna_->bitset_);
#endif

                //       int x = indiv_id / exp_m_->world()->height();
                //       int y = indiv_id % exp_m_->world()->height();
//        printf("Before mutation %d (%d %d)-- %d %d %d\n",indiv_id,x,y,
//               internal_simd_struct[indiv_id]->dna_->length(),dna_size[indiv_id],
//               exp_m_->world()->grid(x, y)->individual()->genetic_unit(0).seq_length());

                /*if (indiv_id == 49) {
                  printf("Size before mutation %d\n",internal_simd_struct[indiv_id]->dna_->length());
                }*/
                if (standalone_)
                    internal_simd_struct[indiv_id]->dna_->apply_mutations_standalone();
                else
                    internal_simd_struct[indiv_id]->dna_->apply_mutations();
            } else {



                #pragma omp atomic
                nb_clones_++;

                int32_t parent_id;
                if (standalone_)
                    parent_id = next_generation_reproducer_[indiv_id];
                else
                    parent_id = next_generation_reproducer_[indiv_id];

                //delete internal_simd_struct[indiv_id];
                internal_simd_struct[indiv_id] = prev_internal_simd_struct[parent_id];

                if (standalone_) {
#pragma omp critical
                    {
                        int x = indiv_id / exp_m_->world()->height();
                        int y = indiv_id % exp_m_->world()->height();
                        NewIndivEvent *eindiv = new NewIndivEvent(internal_simd_struct[indiv_id],
                                                                  prev_internal_simd_struct[next_generation_reproducer_[indiv_id]],
                                                                  x, y,indiv_id,next_generation_reproducer_[indiv_id]);
                        notifyObservers(NEW_INDIV, eindiv);
                        delete eindiv;
                    }
                }

                #pragma omp atomic
                internal_simd_struct[indiv_id]->usage_count_++;

            }

/*
    if (indiv_id == 49) {
      printf("Size after mutation %d\n",internal_simd_struct[indiv_id]->dna_->length());
    }
*/
//    if (indiv_id == 43)
//      printf("DNA AFTER SIZE of %d is %d (%d)\n",indiv_id,dna_size[indiv_id],internal_simd_struct[indiv_id]->dna_->length());

/*
      if (internal_simd_struct[indiv_id]->dna_->mutation_list.size() > 0) {
        int x = indiv_id / exp_m_->world()->height();
        int y = indiv_id % exp_m_->world()->height();

        bool same = true;
        if (internal_simd_struct[indiv_id]->dna_->length() !=
              exp_m_->world()->grid(x, y)->individual()->genetic_unit(0).seq_length())
          same=false;

        if (same)
          for (int i = 0; i < internal_simd_struct[indiv_id]->dna_->length(); i++) {
            if (internal_simd_struct[indiv_id]->dna_->data_[i] != exp_m_->world()->grid(x, y)->individual()->genetic_unit(0).dna()->data()[i]) {
              same = false;
              break;
            }
          }

        // TODO Add promoters verification
        for (auto rna : exp_m_->world()->grid(x, y)->individual()->rna_list()) {
          bool found = false;
          if (rna->strand() == LEADING) {
            for (auto it : internal_simd_struct[indiv_id]->leading_prom_pos) {
              if (it.first == rna->promoter_pos()) {
                found = true;
                break;
              }
            }
          } else {
            for (auto it : internal_simd_struct[indiv_id]->lagging_prom_pos) {
              if (it.first == rna->promoter_pos()) {
                found = true;
                break;
              }
            }
          }

          if (!found) {
            printf("Indiv %d -- Promoter at position %d not found\n",indiv_id,rna->promoter_pos());
            same = false;
            break;
          }
        }

        if (!same) {
          printf("Incorrect mutations replay %d (%d %d)\n",indiv_id,x,y);
          //printf("Aevol (%d) %s\n",exp_m_->world()->grid(x, y)->individual()->genetic_unit(0).dna()->length(),
          //       exp_m_->world()->grid(x, y)->individual()->genetic_unit(0).dna()->data());
          //printf("SIMD (%d) %s\n",internal_simd_struct[indiv_id]->dna_->length(),
          //       internal_simd_struct[indiv_id]->dna_->data_);

          printf("Leading promoters lists : ");
          for (auto it : internal_simd_struct[indiv_id]->leading_prom_pos) {
            printf("%d (%d) || ", it.first, it.second);
          }
          printf("\n");

          printf("Lagging promoters lists : ");
          for (auto it : internal_simd_struct[indiv_id]->lagging_prom_pos) {
            printf("%d (%d) || ", it.first, it.second);
          }
          printf("\n");

          printf("Promoters lists (promoters): ");
          for (auto it : internal_simd_struct[indiv_id]->promoters) {
            printf("%d (%d) -- ", it.second->pos, it.first);
          }
          printf("\n");

          printf("RNA Promoters lists : ");
          for (auto rna : exp_m_->world()->grid(x, y)->individual()->rna_list()) {
            printf("%d ",rna->promoter_pos());
          }

          printf("\n");

          exit(40);
        }

      }*/


            dna_size[indiv_id] = internal_simd_struct[indiv_id]->dna_->length();
        //exit(12);
    }


SIMD_Individual::~SIMD_Individual() {

    printf("Destroy SIMD Controller\n");

    delete stats_best;
    delete stats_mean;

  for (int indiv_id = 0; indiv_id < (int) exp_m_->nb_indivs(); indiv_id++) {

/*
      if (internal_simd_struct[indiv_id] != nullptr) {
          printf("Internal %d : %d \n", indiv_id, internal_simd_struct[indiv_id]->usage_count_);
      }
*/

      if (internal_simd_struct[indiv_id] != nullptr) {


          if (internal_simd_struct[indiv_id]->usage_count_ > 0)
              internal_simd_struct[indiv_id]->usage_count_--;
          else {
              if (internal_simd_struct[indiv_id]->usage_count_ != -1) {
                  internal_simd_struct[indiv_id]->usage_count_ = -1;

                  for (auto rn : internal_simd_struct[indiv_id]->rnas) {
                      delete rn;
                  }
                  internal_simd_struct[indiv_id]->rnas.clear();
                  for (auto prot : internal_simd_struct[indiv_id]->proteins) {
                      delete prot;
                  }
                  internal_simd_struct[indiv_id]->proteins.clear();

                  delete internal_simd_struct[indiv_id];
                  internal_simd_struct[indiv_id] = nullptr;
              }
          }
      }
  }

      for (int indiv_id = 0; indiv_id < (int) exp_m_->nb_indivs(); indiv_id++) {
    if (prev_internal_simd_struct[indiv_id] != nullptr) {
      if (prev_internal_simd_struct[indiv_id]->usage_count_ > 0)
        prev_internal_simd_struct[indiv_id]->usage_count_--;
      else {
          if (prev_internal_simd_struct[indiv_id]->usage_count_ != -1) {
              prev_internal_simd_struct[indiv_id]->usage_count_ = -1;

              /*printf("DELETE --- PREV %d (%d) : %d \n", indiv_id, prev_internal_simd_struct[indiv_id]->usage_count_,
                     prev_internal_simd_struct[indiv_id]->global_id);*/
              for (auto rn : prev_internal_simd_struct[indiv_id]->rnas) {
                  delete rn;
              }
              prev_internal_simd_struct[indiv_id]->rnas.clear();
              for (auto prot : prev_internal_simd_struct[indiv_id]->proteins) {
                  delete prot;
              }
              prev_internal_simd_struct[indiv_id]->proteins.clear();

              delete prev_internal_simd_struct[indiv_id];
              prev_internal_simd_struct[indiv_id] = nullptr;
          }
      }
    }

    delete exp_m_->dna_mutator_array_[indiv_id];
  }

    for (int indiv_id = 0; indiv_id < (int) exp_m_->nb_indivs(); indiv_id++) {
     /*   if (internal_simd_struct[indiv_id] != nullptr) {
            printf("Still here %d : %d\n",internal_simd_struct[indiv_id]->global_id,
                        internal_simd_struct[indiv_id]->usage_count_);
        //}
        //    delete internal_simd_struct[indiv_id];
        }
*/
        if (prev_internal_simd_struct[indiv_id] != nullptr) {
            if (prev_internal_simd_struct[indiv_id]->usage_count_ != -1) {
                prev_internal_simd_struct[indiv_id]->usage_count_ = -1;
                /*printf("Still here %d -- %d : %d\n", indiv_id, prev_internal_simd_struct[indiv_id]->global_id,
                       prev_internal_simd_struct[indiv_id]->usage_count_);*/

                for (auto rn : prev_internal_simd_struct[indiv_id]->rnas) {
                    delete rn;
                }

                prev_internal_simd_struct[indiv_id]->rnas.clear();
                for (auto prot : prev_internal_simd_struct[indiv_id]->proteins) {
                    delete prot;
                }
                prev_internal_simd_struct[indiv_id]->proteins.clear();

                delete prev_internal_simd_struct[indiv_id];
                prev_internal_simd_struct[indiv_id] = nullptr;
            }
        }
    }

/*    for (int indiv_id = 0; indiv_id < (int) exp_m_->nb_indivs(); indiv_id++) {
        *//*   if (internal_simd_struct[indiv_id] != nullptr) {
               printf("Still here %d : %d\n",internal_simd_struct[indiv_id]->global_id,
                           internal_simd_struct[indiv_id]->usage_count_);
           //}
           //    delete internal_simd_struct[indiv_id];
           }
   *//*
        if (prev_internal_simd_struct[indiv_id] != nullptr) {
            printf("WHAT Still here %d -- %d : %d\n",indiv_id,prev_internal_simd_struct[indiv_id]->global_id,
                   prev_internal_simd_struct[indiv_id]->usage_count_);
        }
    }*/

  delete[] prev_internal_simd_struct;
  delete[] internal_simd_struct;

  delete[] dna_size;



  if (standalone_) {
      /*for (int indiv_id = 0; indiv_id < (int) exp_m_->nb_indivs(); indiv_id++) {
          int x = indiv_id / exp_m_->world()->height();
          int y = indiv_id % exp_m_->world()->height();

          if (exp_m_->world()->grid(x, y)->individual()!= nullptr) {
              if (exp_m_->world()->grid(x, y)->individual()->transcribed()) {
                  exp_m_->world()->grid(x, y)->individual()->clear_everything_except_dna_and_promoters();
                  exp_m_->world()->grid(x, y)->individual()->genetic_unit_list_nonconst().clear();
                  delete exp_m_->world()->grid(x, y)->individual();
              }
          }
      }*/
  }
}


void SIMD_Individual::start_stop_RNA() {
  int nb_indiv = exp_m_->nb_indivs();
  //int x, y;

  //
  ExpManager* exp_m = exp_m_;

  //#pragma omp parallel for collapse(2) default(shared)
//#pragma omp parallel for

//#pragma omp parallel
//#pragma omp single
//{
//  #pragma omp parallel for //collapse(2) default(shared)
  for (int indiv_id = 0; indiv_id < nb_indiv; indiv_id++) {
    //internal_simd_struct[indiv_id]->promoters.resize(internal_simd_struct[indiv_id]->dna_->length()/3);

//#pragma omp parallel for firstprivate(indiv_id)
    for (int dna_pos = 0; dna_pos < dna_size[indiv_id]; dna_pos++) {
//#pragma omp task firstprivate(indiv_id,dna_pos)
//      {
#ifdef WITH_BITSET
    //printf("Looking for DNA of %d\n",indiv_id);

      if (internal_simd_struct[indiv_id]->dna_->bitset_->length_ > 0) {

        int dist_lead = internal_simd_struct[indiv_id]->dna_->bitset_->is_promoter(
            true, dna_pos);
        int dist_lag = internal_simd_struct[indiv_id]->dna_->bitset_->is_promoter(
            false, dna_pos);
        bool is_terminator_lead = internal_simd_struct[indiv_id]->dna_->bitset_->is_terminator(
            true, dna_pos);
        bool is_terminator_lag = internal_simd_struct[indiv_id]->dna_->bitset_->is_terminator(
            false, dna_pos);
#else
      int x = indiv_id / exp_m->world()->height();
      int y = indiv_id % exp_m->world()->height();

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
                internal_simd_struct[indiv_id]->dna_->data_[
                    dna_pos - t_motif_id < 0 ? len +
                                               dna_pos -
                                               t_motif_id :
                    dna_pos - t_motif_id]
                ? 0 : 1;
          } else if (motif_id < 22) {
            // LEADING
            prom_dist_leading[motif_id] =
                PROM_SEQ_LEAD[motif_id] ==
                internal_simd_struct[indiv_id]->dna_->data_[
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
                internal_simd_struct[indiv_id]->dna_->data_[
                    dna_pos + t_motif_id >= len ? dna_pos +
                                                  t_motif_id -
                                                  len :
                    dna_pos + t_motif_id] !=
                internal_simd_struct[indiv_id]->dna_->data_[
                    dna_pos - t_motif_id + 10 >= len ?
                    dna_pos - t_motif_id + 10 - len :
                    dna_pos -
                    t_motif_id +
                    10] ? 1
                        : 0;
          } else {
            int t_motif_id = motif_id - 48;
            term_dist_lagging[t_motif_id] =
                internal_simd_struct[indiv_id]->dna_->data_[
                    dna_pos - t_motif_id < 0 ? dna_pos -
                                               t_motif_id +
                                               len
                                             : dna_pos -
                                               t_motif_id] !=
                internal_simd_struct[indiv_id]->dna_->data_[
                    dna_pos + t_motif_id - 10 < 0 ? dna_pos +
                                                    t_motif_id -
                                                    10 + len
                                                  :
                    dna_pos + t_motif_id - 10] ? 1 : 0;

            /*if (dna_pos == 22201 && indiv_id == 309 && AeTime::time() == 105) {
              printf("Term @ 7 Motif ID %d : %c %c (%d %d) Error %d\n",t_motif_id,
                     exp_m_->world()->grid(x, y)->individual()->genetic_unit(
                         0).dna()->data()[dna_pos - t_motif_id < 0 ? dna_pos - t_motif_id +len : dna_pos - t_motif_id],
                     exp_m_->world()->grid(x, y)->individual()->genetic_unit(
                         0).dna()->data()[dna_pos + t_motif_id - 10 < 0 ? dna_pos + t_motif_id - 10 +len : dna_pos + t_motif_id - 10],
                     dna_pos - t_motif_id < 0 ? dna_pos - t_motif_id +len : dna_pos - t_motif_id,
                     dna_pos + t_motif_id - 10 < 0 ? dna_pos + t_motif_id - 10 +len : dna_pos + t_motif_id - 10,
                     term_dist_lagging[t_motif_id]);
            }*/
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
#endif

        if (dist_lead <= 4) {
          promoterStruct* nprom = new promoterStruct(dna_pos, dist_lead,
                                                     true);

//#pragma omp critical(add_to_promoters)
          {
            int prom_idx;

            prom_idx = internal_simd_struct[indiv_id]->count_prom;
            internal_simd_struct[indiv_id]->count_prom =
                internal_simd_struct[indiv_id]->count_prom + 1;
/*
              if (indiv_id == 6)
                printf("Adding promoters %d at %d\n",dna_pos,prom_idx);*/

            internal_simd_struct[indiv_id]->promoters[prom_idx] = nprom;
            internal_simd_struct[indiv_id]->leading_prom_pos[dna_pos] = prom_idx;
          }
        }

#ifndef WITH_BITSET
        int dist_term_lead = term_dist_leading[0] +
                             term_dist_leading[1] +
                             term_dist_leading[2] +
                             term_dist_leading[3];
        /*if (dna_pos >= 2536 && indiv_id == 915 && AeTime::time() == 26) {
          printf("Distance for %d : %d\n",dna_pos,dist_term_lead);
        }*/


        if (dist_term_lead == 4) {
#else
          if (is_terminator_lead) {
#endif
          internal_simd_struct[indiv_id]->terminator_lead.insert(
              dna_pos);
        }

#ifndef WITH_BITSET
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
#endif

        if (dist_lag <= 4) {
          promoterStruct* nprom = new promoterStruct(dna_pos, dist_lag,
                                                     false);
//#pragma omp critical(add_to_promoters)
          {
            int prom_idx;
            prom_idx = internal_simd_struct[indiv_id]->count_prom;
            internal_simd_struct[indiv_id]->count_prom =
                internal_simd_struct[indiv_id]->count_prom + 1;
            internal_simd_struct[indiv_id]->promoters[prom_idx] = nprom;
            internal_simd_struct[indiv_id]->lagging_prom_pos[dna_pos] = prom_idx;
          }
        }

#ifndef WITH_BITSET
        int dist_term_lag = term_dist_lagging[0] +
                            term_dist_lagging[1] +
                            term_dist_lagging[2] +
                            term_dist_lagging[3];


        if (dist_term_lag == 4) {
#else
          if (is_terminator_lag) {
#endif
          internal_simd_struct[indiv_id]->terminator_lag.insert(
              dna_pos);
        }
      }
    }
  }
//}
//}
//#pragma omp taskwait
}



    void SIMD_Individual::start_stop_RNA(int indiv_id) {
        int nb_indiv = exp_m_->nb_indivs();
        //int x, y;

        //
        ExpManager* exp_m = exp_m_;

        //#pragma omp parallel for collapse(2) default(shared)
//#pragma omp parallel for

//#pragma omp parallel
//#pragma omp single
//{
//  #pragma omp parallel for //collapse(2) default(shared)
            //internal_simd_struct[indiv_id]->promoters.resize(internal_simd_struct[indiv_id]->dna_->length()/3);

//#pragma omp parallel for firstprivate(indiv_id)
            for (int dna_pos = 0; dna_pos < dna_size[indiv_id]; dna_pos++) {
//#pragma omp task firstprivate(indiv_id,dna_pos)
//      {
#ifdef WITH_BITSET
                //printf("Looking for DNA of %d\n",indiv_id);

      if (internal_simd_struct[indiv_id]->dna_->bitset_->length_ > 0) {

        int dist_lead = internal_simd_struct[indiv_id]->dna_->bitset_->is_promoter(
            true, dna_pos);
        int dist_lag = internal_simd_struct[indiv_id]->dna_->bitset_->is_promoter(
            false, dna_pos);
        bool is_terminator_lead = internal_simd_struct[indiv_id]->dna_->bitset_->is_terminator(
            true, dna_pos);
        bool is_terminator_lag = internal_simd_struct[indiv_id]->dna_->bitset_->is_terminator(
            false, dna_pos);
#else
                int x = indiv_id / exp_m->world()->height();
                int y = indiv_id % exp_m->world()->height();

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
                                    internal_simd_struct[indiv_id]->dna_->data_[
                                            dna_pos - t_motif_id < 0 ? len +
                                                                       dna_pos -
                                                                       t_motif_id :
                                            dna_pos - t_motif_id]
                                    ? 0 : 1;
                        } else if (motif_id < 22) {
                            // LEADING
                            prom_dist_leading[motif_id] =
                                    PROM_SEQ_LEAD[motif_id] ==
                                    internal_simd_struct[indiv_id]->dna_->data_[
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
                                    internal_simd_struct[indiv_id]->dna_->data_[
                                            dna_pos + t_motif_id >= len ? dna_pos +
                                                                          t_motif_id -
                                                                          len :
                                            dna_pos + t_motif_id] !=
                                    internal_simd_struct[indiv_id]->dna_->data_[
                                            dna_pos - t_motif_id + 10 >= len ?
                                            dna_pos - t_motif_id + 10 - len :
                                            dna_pos -
                                            t_motif_id +
                                            10] ? 1
                                                : 0;
                        } else {
                            int t_motif_id = motif_id - 48;
                            term_dist_lagging[t_motif_id] =
                                    internal_simd_struct[indiv_id]->dna_->data_[
                                            dna_pos - t_motif_id < 0 ? dna_pos -
                                                                       t_motif_id +
                                                                       len
                                                                     : dna_pos -
                                                                       t_motif_id] !=
                                    internal_simd_struct[indiv_id]->dna_->data_[
                                            dna_pos + t_motif_id - 10 < 0 ? dna_pos +
                                                                            t_motif_id -
                                                                            10 + len
                                                                          :
                                            dna_pos + t_motif_id - 10] ? 1 : 0;

                            /*if (dna_pos == 22201 && indiv_id == 309 && AeTime::time() == 105) {
                              printf("Term @ 7 Motif ID %d : %c %c (%d %d) Error %d\n",t_motif_id,
                                     exp_m_->world()->grid(x, y)->individual()->genetic_unit(
                                         0).dna()->data()[dna_pos - t_motif_id < 0 ? dna_pos - t_motif_id +len : dna_pos - t_motif_id],
                                     exp_m_->world()->grid(x, y)->individual()->genetic_unit(
                                         0).dna()->data()[dna_pos + t_motif_id - 10 < 0 ? dna_pos + t_motif_id - 10 +len : dna_pos + t_motif_id - 10],
                                     dna_pos - t_motif_id < 0 ? dna_pos - t_motif_id +len : dna_pos - t_motif_id,
                                     dna_pos + t_motif_id - 10 < 0 ? dna_pos + t_motif_id - 10 +len : dna_pos + t_motif_id - 10,
                                     term_dist_lagging[t_motif_id]);
                            }*/
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
#endif

                    if (dist_lead <= 4) {
                        promoterStruct* nprom = new promoterStruct(dna_pos, dist_lead,
                                                                   true);

//#pragma omp critical(add_to_promoters)
                        {
                            int prom_idx;

                            prom_idx = internal_simd_struct[indiv_id]->count_prom;
                            internal_simd_struct[indiv_id]->count_prom =
                                    internal_simd_struct[indiv_id]->count_prom + 1;
/*
              if (indiv_id == 6)
                printf("Adding promoters %d at %d\n",dna_pos,prom_idx);*/

                            internal_simd_struct[indiv_id]->promoters[prom_idx] = nprom;
                            internal_simd_struct[indiv_id]->leading_prom_pos[dna_pos] = prom_idx;
                        }
                    }

#ifndef WITH_BITSET
                    int dist_term_lead = term_dist_leading[0] +
                                         term_dist_leading[1] +
                                         term_dist_leading[2] +
                                         term_dist_leading[3];
                    /*if (dna_pos >= 2536 && indiv_id == 915 && AeTime::time() == 26) {
                      printf("Distance for %d : %d\n",dna_pos,dist_term_lead);
                    }*/


                    if (dist_term_lead == 4) {
#else
                        if (is_terminator_lead) {
#endif
                        internal_simd_struct[indiv_id]->terminator_lead.insert(
                                dna_pos);
                    }

#ifndef WITH_BITSET
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
#endif

                    if (dist_lag <= 4) {
                        promoterStruct* nprom = new promoterStruct(dna_pos, dist_lag,
                                                                   false);
//#pragma omp critical(add_to_promoters)
                        {
                            int prom_idx;
                            prom_idx = internal_simd_struct[indiv_id]->count_prom;
                            internal_simd_struct[indiv_id]->count_prom =
                                    internal_simd_struct[indiv_id]->count_prom + 1;
                            internal_simd_struct[indiv_id]->promoters[prom_idx] = nprom;
                            internal_simd_struct[indiv_id]->lagging_prom_pos[dna_pos] = prom_idx;
                        }
                    }

#ifndef WITH_BITSET
                    int dist_term_lag = term_dist_lagging[0] +
                                        term_dist_lagging[1] +
                                        term_dist_lagging[2] +
                                        term_dist_lagging[3];


                    if (dist_term_lag == 4) {
#else
                        if (is_terminator_lag) {
#endif
                        internal_simd_struct[indiv_id]->terminator_lag.insert(
                                dna_pos);
                    }
                }
            }

//}
//}
//#pragma omp taskwait
    }


void SIMD_Individual::opt_prom_compute_RNA() {

  int nb_indiv = exp_m_->nb_indivs();


  for (int indiv_id = 0; indiv_id < nb_indiv; indiv_id++) {
    if (exp_m_->dna_mutator_array_[indiv_id]->hasMutate()) {
      internal_simd_struct[indiv_id]->proteins.clear();
      internal_simd_struct[indiv_id]->rnas.clear();
      internal_simd_struct[indiv_id]->terminator_lead.clear();
      internal_simd_struct[indiv_id]->terminator_lag.clear();
    }
  }


//#pragma omp parallel
//#pragma omp single
  //{
//#pragma omp parallel for schedule(dynamic)
    for (int indiv_id = 0; indiv_id < nb_indiv; indiv_id++) {
      if (exp_m_->dna_mutator_array_[indiv_id]->hasMutate()) {
        internal_simd_struct[indiv_id]->rnas.resize(
            internal_simd_struct[indiv_id]->promoters.size());
        /*internal_simd_struct[indiv_id]->proteins.clear();
        internal_simd_struct[indiv_id]->rnas.clear();
        internal_simd_struct[indiv_id]->terminator_lead.clear();
        internal_simd_struct[indiv_id]->terminator_lag.clear();*/
//#pragma omp parallel for schedule(dynamic)
//#pragma omp task firstprivate(indiv_id)
//          {

        for (auto promx : internal_simd_struct[indiv_id]->promoters) {
          int rna_idx = promx.first;
            promoterStruct* prom;
            prom = internal_simd_struct[indiv_id]->promoters[rna_idx];
            //internal_simd_struct[indiv_id]->rnas[rna_idx] = nullptr;


            if (prom != nullptr) {
              int prom_pos;
              bool lead_lag;
              double prom_error;
              prom_pos = internal_simd_struct[indiv_id]->promoters[rna_idx]->pos;
              lead_lag = internal_simd_struct[indiv_id]->promoters[rna_idx]->leading_or_lagging;
              prom_error = fabs(
                  ((float) internal_simd_struct[indiv_id]->promoters[rna_idx]->error));


              if (lead_lag) {
//        if (indiv_id == 152) printf("Searching for RNA (OPT) for indiv %d RNA %d LEAD\n",indiv_id,rna_idx);
                /* Search for terminators */
                int cur_pos =
                    prom_pos + 22;
                cur_pos = cur_pos >= dna_size[indiv_id] ? cur_pos -
                                                          dna_size[indiv_id] :
                          cur_pos;
                int start_pos = cur_pos;

                bool terminator_found = false;
                bool no_terminator = false;
                int term_dist_leading = 0;

                int loop_size = 0;

                while (!terminator_found) {
                  loop_size++;
#ifdef WITH_BITSET
                    if (cur_pos >= dna_size[indiv_id] || cur_pos < 0) {
                      printf("-----LEAD-----> Position %d Length %d\n",cur_pos,dna_size[indiv_id]);
                  }

                  bool is_term = internal_simd_struct[indiv_id]->dna_->bitset_->is_terminator(
                      true, cur_pos);

                  if (is_term)
#else
                  //#pragma omp simd aligned(internal_simd_struct[indiv_id]->dna_->data_:64)
                  for (int t_motif_id = 0; t_motif_id < 4; t_motif_id++)
                    term_dist_leading +=
                        internal_simd_struct[indiv_id]->dna_->data_[
                            cur_pos + t_motif_id >= dna_size[indiv_id] ?
                            cur_pos +
                            t_motif_id -
                            dna_size[indiv_id] :
                            cur_pos + t_motif_id] !=
                        internal_simd_struct[indiv_id]->dna_->data_[
                            cur_pos - t_motif_id + 10 >= dna_size[indiv_id] ?
                            cur_pos - t_motif_id + 10 - dna_size[indiv_id] :
                            cur_pos -
                            t_motif_id +
                            10] ? 1
                                : 0;

                  if (term_dist_leading == 4)
#endif
                    terminator_found = true;
                  else {
                    cur_pos = cur_pos + 1 >= dna_size[indiv_id] ? cur_pos + 1 -
                                                                  dna_size[indiv_id]
                                                                :
                              cur_pos + 1;

//            if (indiv_id == 152) printf("Next cur value %d (prev %d)\n",cur_pos,term_dist_leading);
                    term_dist_leading = 0;
                    if (cur_pos == start_pos) {
                      no_terminator = true;
                      terminator_found = true;
                      //wno_terminator[0] = true;
                    }
                  }
                }

//        if (indiv_id == 152) printf("LOOP SIZE %d : start %d length %ld\n",loop_size,start_pos,dna_size[indiv_id]);

                if (!no_terminator) {

                  int32_t rna_end =
                      cur_pos + 10 >= dna_size[indiv_id] ?
                      cur_pos + 10 - dna_size[indiv_id] :
                      cur_pos + 10;

//        if (indiv_id == 152) printf("Adding new RNA %d (%d)\n",cur_pos,rna_end);
                  /*if (indiv_id == 309 && AeTime::time() == 105) {
                    printf("Looking for term from %d (start rna %d) : %d Computed end %d\n",k,internal_simd_struct[indiv_id]->promoters[rna_idx]->pos,
                           *it_rna_end,rna_end);
                  }*/
                  int32_t rna_length = 0;

                  if (prom_pos
                      > rna_end)
                    rna_length = dna_size[indiv_id] -
                                 prom_pos
                                 + rna_end;
                  else
                    rna_length = rna_end - prom_pos;

                  rna_length -= 21;

                  if (rna_length > 0) {
                    int glob_rna_idx = -1;
#pragma omp atomic capture
                    {
                      glob_rna_idx = internal_simd_struct[indiv_id]->rna_count_;
                      internal_simd_struct[indiv_id]->rna_count_ =
                          internal_simd_struct[indiv_id]->rna_count_ + 1;
                    }

                    internal_simd_struct[indiv_id]->rnas[glob_rna_idx] = new pRNA(
                        prom_pos,
                        rna_end,
                        !lead_lag,
                        1.0 -
                        prom_error /
                        5.0, rna_length);
                  }
                }
//        if (indiv_id == 152) printf("Hop to next\n");
              } else {
                /* Search for terminator */
                int cur_pos =
                    prom_pos - 22;
                cur_pos =
                    cur_pos < 0 ? dna_size[indiv_id] + (cur_pos) : cur_pos;
                int start_pos = cur_pos;
                bool terminator_found = false;
                bool no_terminator = false;
                int term_dist_lagging = 0;

                //if (indiv_id == 180) printf("Searching for RNA (OPT) for indiv %d RNA %d start at %d\n",indiv_id,rna_idx,start_pos);
                int loop_size = 0;

                while (!terminator_found) {
#ifdef WITH_BITSET
                  if (cur_pos >= dna_size[indiv_id] || cur_pos < 0) {
                      printf("---LAG------------> Position %d Length %d\n",cur_pos,dna_size[indiv_id]);
                  }

                  bool is_term = internal_simd_struct[indiv_id]->dna_->bitset_->is_terminator(
                      false, cur_pos);

                  if (is_term)
#else
                  //#pragma omp simd aligned(internal_simd_struct[indiv_id]->dna_->data_:64)
                  for (int t_motif_id = 0; t_motif_id < 4; t_motif_id++) {
                    term_dist_lagging +=
                        internal_simd_struct[indiv_id]->dna_->data_[
                            cur_pos - t_motif_id < 0 ? cur_pos -
                                                       t_motif_id +
                                                       dna_size[indiv_id]
                                                     : cur_pos -
                                                       t_motif_id] !=
                        internal_simd_struct[indiv_id]->dna_->data_[
                            cur_pos + t_motif_id - 10 < 0 ? cur_pos +
                                                            t_motif_id -
                                                            10 +
                                                            dna_size[indiv_id]
                                                          :
                            cur_pos + t_motif_id - 10] ? 1 : 0;
                  }

                  if (term_dist_lagging == 4)
#endif
                    terminator_found = true;
                  else {
                    //if (indiv_id == 180 && cur_pos - 1 < 0)
                    //printf("WHAT BEFORE ??? %d SIZE %d PREDIRECT %d\n",cur_pos,dna_size[indiv_id],
                    //       dna_size[indiv_id] + (cur_pos - 1));

                    cur_pos =
                        cur_pos - 1 < 0 ? dna_size[indiv_id] + (cur_pos - 1)
                                        : cur_pos - 1;

//            if (indiv_id == 180 && cur_pos > dna_size[indiv_id] - 1) {
//              printf("WHAT AFTER ??? %d SIZE %d\n", cur_pos,
//                     dna_size[indiv_id]);
//              exit(1654);
//            }

                    term_dist_lagging = 0;
                    if (cur_pos == start_pos) {
                      no_terminator = true;
                      terminator_found = true;
                      //wno_terminator[1] = true;
                    }
                  }
                  loop_size++;
                }

//        if (indiv_id == 152) printf("LOOP SIZE %d : start %d length %ld\n",loop_size,start_pos,dna_size[indiv_id]);

                if (!no_terminator) {


                  int32_t rna_end =
                      cur_pos - 10 < 0 ? dna_size[indiv_id] + (cur_pos - 10) :
                      cur_pos -
                      10;
//        if (indiv_id == 180) printf("Adding new RNA %d (%d) -- %d\n",cur_pos,rna_end,dna_size[indiv_id]);
                  /*if (indiv_id == 969 && AeTime::time() == 137) {
                    auto it_rn = it_rna_end;
                    it_rn++;
                    printf("Looking for term from %d (start rna %d) : %d Computed end %d (next end %d)\n",k,internal_simd_struct[indiv_id]->promoters[rna_idx]->pos,
                           *it_rna_end,rna_end,*it_rn);
                  }*/
                  int32_t rna_length = 0;

                  if (prom_pos <
                      rna_end)
                    rna_length =
                        prom_pos +
                        dna_size[indiv_id] - rna_end;
                  else
                    rna_length =
                        prom_pos -
                        rna_end;

                  rna_length -= 21;

                  if (rna_length >= 0) {
                    int glob_rna_idx = -1;
#pragma omp atomic capture
                    {
                      glob_rna_idx = internal_simd_struct[indiv_id]->rna_count_;
                      internal_simd_struct[indiv_id]->rna_count_ =
                          internal_simd_struct[indiv_id]->rna_count_ + 1;
                    }

                    internal_simd_struct[indiv_id]->rnas[glob_rna_idx] = new pRNA(
                        prom_pos,
                        rna_end,
                        !lead_lag,
                        1.0 -
                        prom_error /
                        5.0, rna_length);
                  }
//        if (indiv_id == 152) printf("Hop to next\n");

                }
              }
            }
          }
        }
//    if (indiv_id == 152) printf("--------> %d -- With no terminator %d %d\n",indiv_id,
//           wno_terminator[0],wno_terminator[1]);
      }
    //}
  //}
}



    void SIMD_Individual::opt_prom_compute_RNA(int indiv_id) {

        int nb_indiv = exp_m_->nb_indivs();


            if (exp_m_->dna_mutator_array_[indiv_id]->hasMutate()) {
                internal_simd_struct[indiv_id]->proteins.clear();
                internal_simd_struct[indiv_id]->rnas.clear();
                internal_simd_struct[indiv_id]->terminator_lead.clear();
                internal_simd_struct[indiv_id]->terminator_lag.clear();
            }



//#pragma omp parallel
//#pragma omp single
        //{
//#pragma omp parallel for schedule(dynamic)
        //for (int indiv_id = 0; indiv_id < nb_indiv; indiv_id++) {
            if (exp_m_->dna_mutator_array_[indiv_id]->hasMutate()) {
                internal_simd_struct[indiv_id]->rnas.resize(
                        internal_simd_struct[indiv_id]->promoters.size());
                /*internal_simd_struct[indiv_id]->proteins.clear();
                internal_simd_struct[indiv_id]->rnas.clear();
                internal_simd_struct[indiv_id]->terminator_lead.clear();
                internal_simd_struct[indiv_id]->terminator_lag.clear();*/
//#pragma omp parallel for schedule(dynamic)
//#pragma omp task firstprivate(indiv_id)
//          {

#ifdef WITH_FINETASKLOOP
#pragma omp taskloop grainsize(rna_grain_size)
#endif
                for (int prom_idx = 0; prom_idx< internal_simd_struct[indiv_id]->promoters.size(); prom_idx++) {

                    if (internal_simd_struct[indiv_id]->promoters[prom_idx] != nullptr) {
                        int rna_idx = prom_idx;
                        promoterStruct *prom;
                        prom = internal_simd_struct[indiv_id]->promoters[rna_idx];
                        //internal_simd_struct[indiv_id]->rnas[rna_idx] = nullptr;


                        if (prom != nullptr) {
                            int prom_pos;
                            bool lead_lag;
                            double prom_error;
                            prom_pos = internal_simd_struct[indiv_id]->promoters[rna_idx]->pos;
                            lead_lag = internal_simd_struct[indiv_id]->promoters[rna_idx]->leading_or_lagging;
                            prom_error = fabs(
                                    ((float) internal_simd_struct[indiv_id]->promoters[rna_idx]->error));

                            //if (indiv_id == 6) printf("Searching for RNA (OPT) for indiv %d RNA %d starting %d\n",indiv_id,rna_idx,prom_pos);

                            if (lead_lag) {

                                /* Search for terminators */
                                int cur_pos =
                                        prom_pos + 22;
                                cur_pos = cur_pos >= dna_size[indiv_id] ? cur_pos -
                                                                          dna_size[indiv_id] :
                                          cur_pos;
                                int start_pos = cur_pos;

                                bool terminator_found = false;
                                bool no_terminator = false;
                                int term_dist_leading = 0;

                                int loop_size = 0;

                                while (!terminator_found) {
                                    loop_size++;
#ifdef WITH_BITSET
                                    if (cur_pos >= dna_size[indiv_id] || cur_pos < 0) {
                          printf("-----LEAD-----> Position %d Length %d\n",cur_pos,dna_size[indiv_id]);
                      }

                      bool is_term = internal_simd_struct[indiv_id]->dna_->bitset_->is_terminator(
                          true, cur_pos);

                      if (is_term)
#else
                                    //#pragma omp simd aligned(internal_simd_struct[indiv_id]->dna_->data_:64)
                                    for (int t_motif_id = 0; t_motif_id < 4; t_motif_id++)
                                        term_dist_leading +=
                                                internal_simd_struct[indiv_id]->dna_->data_[
                                                        cur_pos + t_motif_id >= dna_size[indiv_id] ?
                                                        cur_pos +
                                                        t_motif_id -
                                                        dna_size[indiv_id] :
                                                        cur_pos + t_motif_id] !=
                                                internal_simd_struct[indiv_id]->dna_->data_[
                                                        cur_pos - t_motif_id + 10 >= dna_size[indiv_id] ?
                                                        cur_pos - t_motif_id + 10 - dna_size[indiv_id] :
                                                        cur_pos -
                                                        t_motif_id +
                                                        10] ? 1
                                                            : 0;

                                    if (term_dist_leading == 4)
#endif
                                        terminator_found = true;
                                    else {
                                        cur_pos = cur_pos + 1 >= dna_size[indiv_id] ? cur_pos + 1 -
                                                                                      dna_size[indiv_id]
                                                                                    :
                                                  cur_pos + 1;

//            if (indiv_id == 152) printf("Next cur value %d (prev %d)\n",cur_pos,term_dist_leading);
                                        term_dist_leading = 0;
                                        if (cur_pos == start_pos) {
                                            no_terminator = true;
                                            terminator_found = true;
                                            //wno_terminator[0] = true;
                                        }
                                    }
                                }

//        if (indiv_id == 152) printf("LOOP SIZE %d : start %d length %ld\n",loop_size,start_pos,dna_size[indiv_id]);

                                if (!no_terminator) {

                                    int32_t rna_end =
                                            cur_pos + 10 >= dna_size[indiv_id] ?
                                            cur_pos + 10 - dna_size[indiv_id] :
                                            cur_pos + 10;

//        if (indiv_id == 152) printf("Adding new RNA %d (%d)\n",cur_pos,rna_end);
                                    /*if (indiv_id == 309 && AeTime::time() == 105) {
                                      printf("Looking for term from %d (start rna %d) : %d Computed end %d\n",k,internal_simd_struct[indiv_id]->promoters[rna_idx]->pos,
                                             *it_rna_end,rna_end);
                                    }*/
                                    int32_t rna_length = 0;

                                    if (prom_pos
                                        > rna_end)
                                        rna_length = dna_size[indiv_id] -
                                                     prom_pos
                                                     + rna_end;
                                    else
                                        rna_length = rna_end - prom_pos;

                                    rna_length -= 21;

                                    if (rna_length > 0) {
                                        int glob_rna_idx = -1;
#pragma omp atomic capture
                                        {
                                            glob_rna_idx = internal_simd_struct[indiv_id]->rna_count_;
                                            internal_simd_struct[indiv_id]->rna_count_ =
                                                    internal_simd_struct[indiv_id]->rna_count_ + 1;
                                        }

                                        internal_simd_struct[indiv_id]->rnas[glob_rna_idx] = new pRNA(
                                                prom_pos,
                                                rna_end,
                                                !lead_lag,
                                                1.0 -
                                                prom_error /
                                                5.0, rna_length);
                                    }
                                }
//        if (indiv_id == 152) printf("Hop to next\n");
                            } else {
                                /* Search for terminator */
                                int cur_pos =
                                        prom_pos - 22;
                                cur_pos =
                                        cur_pos < 0 ? dna_size[indiv_id] + (cur_pos) : cur_pos;
                                int start_pos = cur_pos;
                                bool terminator_found = false;
                                bool no_terminator = false;
                                int term_dist_lagging = 0;

                                //if (indiv_id == 180) printf("Searching for RNA (OPT) for indiv %d RNA %d start at %d\n",indiv_id,rna_idx,start_pos);
                                int loop_size = 0;

                                while (!terminator_found) {
#ifdef WITH_BITSET
                                    if (cur_pos >= dna_size[indiv_id] || cur_pos < 0) {
                          printf("---LAG------------> Position %d Length %d\n",cur_pos,dna_size[indiv_id]);
                      }

                      bool is_term = internal_simd_struct[indiv_id]->dna_->bitset_->is_terminator(
                          false, cur_pos);

                      if (is_term)
#else
                                    //#pragma omp simd aligned(internal_simd_struct[indiv_id]->dna_->data_:64)
                                    for (int t_motif_id = 0; t_motif_id < 4; t_motif_id++) {
                                        term_dist_lagging +=
                                                internal_simd_struct[indiv_id]->dna_->data_[
                                                        cur_pos - t_motif_id < 0 ? cur_pos -
                                                                                   t_motif_id +
                                                                                   dna_size[indiv_id]
                                                                                 : cur_pos -
                                                                                   t_motif_id] !=
                                                internal_simd_struct[indiv_id]->dna_->data_[
                                                        cur_pos + t_motif_id - 10 < 0 ? cur_pos +
                                                                                        t_motif_id -
                                                                                        10 +
                                                                                        dna_size[indiv_id]
                                                                                      :
                                                        cur_pos + t_motif_id - 10] ? 1 : 0;
                                    }

                                    if (term_dist_lagging == 4)
#endif
                                        terminator_found = true;
                                    else {
                                        //if (indiv_id == 180 && cur_pos - 1 < 0)
                                        //printf("WHAT BEFORE ??? %d SIZE %d PREDIRECT %d\n",cur_pos,dna_size[indiv_id],
                                        //       dna_size[indiv_id] + (cur_pos - 1));

                                        cur_pos =
                                                cur_pos - 1 < 0 ? dna_size[indiv_id] + (cur_pos - 1)
                                                                : cur_pos - 1;

//            if (indiv_id == 180 && cur_pos > dna_size[indiv_id] - 1) {
//              printf("WHAT AFTER ??? %d SIZE %d\n", cur_pos,
//                     dna_size[indiv_id]);
//              exit(1654);
//            }

                                        term_dist_lagging = 0;
                                        if (cur_pos == start_pos) {
                                            no_terminator = true;
                                            terminator_found = true;
                                            //wno_terminator[1] = true;
                                        }
                                    }
                                    loop_size++;
                                }

//        if (indiv_id == 152) printf("LOOP SIZE %d : start %d length %ld\n",loop_size,start_pos,dna_size[indiv_id]);

                                if (!no_terminator) {


                                    int32_t rna_end =
                                            cur_pos - 10 < 0 ? dna_size[indiv_id] + (cur_pos - 10) :
                                            cur_pos -
                                            10;
//        if (indiv_id == 180) printf("Adding new RNA %d (%d) -- %d\n",cur_pos,rna_end,dna_size[indiv_id]);
                                    /*if (indiv_id == 969 && AeTime::time() == 137) {
                                      auto it_rn = it_rna_end;
                                      it_rn++;
                                      printf("Looking for term from %d (start rna %d) : %d Computed end %d (next end %d)\n",k,internal_simd_struct[indiv_id]->promoters[rna_idx]->pos,
                                             *it_rna_end,rna_end,*it_rn);
                                    }*/
                                    int32_t rna_length = 0;

                                    if (prom_pos <
                                        rna_end)
                                        rna_length =
                                                prom_pos +
                                                dna_size[indiv_id] - rna_end;
                                    else
                                        rna_length =
                                                prom_pos -
                                                rna_end;

                                    rna_length -= 21;

                                    if (rna_length >= 0) {
                                        int glob_rna_idx = -1;
#pragma omp atomic capture
                                        {
                                            glob_rna_idx = internal_simd_struct[indiv_id]->rna_count_;
                                            internal_simd_struct[indiv_id]->rna_count_ =
                                                    internal_simd_struct[indiv_id]->rna_count_ + 1;
                                        }

                                        internal_simd_struct[indiv_id]->rnas[glob_rna_idx] = new pRNA(
                                                prom_pos,
                                                rna_end,
                                                !lead_lag,
                                                1.0 -
                                                prom_error /
                                                5.0, rna_length);
                                    }
//        if (indiv_id == 152) printf("Hop to next\n");

                                }
                            }
                        }
                    }
                }
            }
//    if (indiv_id == 152) printf("--------> %d -- With no terminator %d %d\n",indiv_id,
//           wno_terminator[0],wno_terminator[1]);
        //}
        //}
    }

void SIMD_Individual::compute_RNA() {

  int nb_indiv = exp_m_->nb_indivs();
//#pragma omp parallel for schedule(dynamic)
#pragma omp parallel
#pragma omp single
    {
        for (int indiv_id = 0; indiv_id < nb_indiv; indiv_id++) {
            internal_simd_struct[indiv_id]->rnas.resize(
                    internal_simd_struct[indiv_id]->promoters.size());
//#pragma omp parallel for firstprivate(indiv_id) schedule(dynamic)
            for (int rna_idx = 0; rna_idx <
                                  (int) internal_simd_struct[indiv_id]->promoters.size(); rna_idx++) {
#pragma omp task firstprivate(indiv_id, rna_idx)
                {
                    /*if (indiv_id == 345 && AeTime::time() == 47 &&
                        internal_simd_struct[indiv_id]->promoters[rna_idx]->pos == 4744) {
                      printf("Searching for an end with start pos %d LorL %d\n",
                             internal_simd_struct[indiv_id]->promoters[rna_idx]->pos,
                             internal_simd_struct[indiv_id]->promoters[rna_idx]->leading_or_lagging);
                    }*/

                    if (internal_simd_struct[indiv_id]->promoters[rna_idx] != nullptr) {
                        if (internal_simd_struct[indiv_id]->promoters[rna_idx]->leading_or_lagging) {
                            if (internal_simd_struct[indiv_id]->terminator_lead.size() != 0) {


                                int k =
                                        internal_simd_struct[indiv_id]->promoters[rna_idx]->pos + 22;
                                k = k >= dna_size[indiv_id] ? k - dna_size[indiv_id] : k;

/*        if ((indiv_id == 309 && AeTime::time() == 105) ||
            (indiv_id == 915 && AeTime::time() == 26 && rna_idx == 11)) {
          printf("Looking at %d\n",k);
        }*/

                                auto it_rna_end = internal_simd_struct[indiv_id]->terminator_lead.lower_bound(
                                        k);

                                /*if ((indiv_id == 309 && AeTime::time() == 105) ||
                                    (indiv_id == 915 && AeTime::time() == 26 && rna_idx == 11)) {
                                  printf("Looking at %d\n",*it_rna_end);
                                }*/


                                if (it_rna_end ==
                                    internal_simd_struct[indiv_id]->terminator_lead.end()) {
                                    it_rna_end = internal_simd_struct[indiv_id]->terminator_lead.begin();

                                    /* if ((indiv_id == 309 && AeTime::time() == 105) ||
                                         (indiv_id == 915 && AeTime::time() == 26 && rna_idx == 11)) {
                                       printf("Found at %d (restart because end)\n",*it_rna_end);
                                     }*/
                                }

                                /*if ((indiv_id == 309 && AeTime::time() == 105) ||
                                    (indiv_id == 915 && AeTime::time() == 26 && rna_idx == 11)) {
                                  auto tmp_rna_end = it_rna_end;
                                  printf("Terminators for %d current %d : ",rna_idx,*it_rna_end);
                                  while (tmp_rna_end != internal_simd_struct[indiv_id]->terminator_lead.end()) {
                                    printf("%d ",(*tmp_rna_end));
                                    tmp_rna_end++;
                                  }
                                  printf("\n");
                                }*/

                                int32_t rna_end =
                                        *it_rna_end + 10 >= dna_size[indiv_id] ?
                                        *it_rna_end + 10 - dna_size[indiv_id] :
                                        *it_rna_end + 10;

                                /*if (indiv_id == 309 && AeTime::time() == 105) {
                                  printf("Looking for term from %d (start rna %d) : %d Computed end %d\n",k,internal_simd_struct[indiv_id]->promoters[rna_idx]->pos,
                                         *it_rna_end,rna_end);
                                }*/
                                int32_t rna_length = 0;

                                if (internal_simd_struct[indiv_id]->promoters[rna_idx]->pos
                                    > rna_end)
                                    rna_length = dna_size[indiv_id] -
                                                 internal_simd_struct[indiv_id]->promoters[rna_idx]->pos
                                                 + rna_end;
                                else
                                    rna_length = rna_end - internal_simd_struct[indiv_id]->
                                            promoters[rna_idx]->pos;

                                rna_length -= 21;

                                if (rna_length >= 0) {


                                    int glob_rna_idx = -1;
#pragma omp atomic capture
                                    {
                                        glob_rna_idx = internal_simd_struct[indiv_id]->rna_count_;
                                        internal_simd_struct[indiv_id]->rna_count_ =
                                                internal_simd_struct[indiv_id]->rna_count_ + 1;
                                    }

                                    internal_simd_struct[indiv_id]->rnas[glob_rna_idx] = new pRNA(
                                            internal_simd_struct[indiv_id]->promoters[rna_idx]->pos,
                                            rna_end,
                                            !internal_simd_struct[indiv_id]->promoters[rna_idx]->leading_or_lagging,
                                            1.0 -
                                            fabs(
                                                    ((float) internal_simd_struct[indiv_id]->promoters[rna_idx]->error)) /
                                            5.0, rna_length);
                                }
                            }
                        } else {
                            // LAGGING
                            if (internal_simd_struct[indiv_id]->terminator_lag.size() != 0) {




                                // Search for terminator
                                int k =
                                        internal_simd_struct[indiv_id]->promoters[rna_idx]->pos - 22;
                                k = k < 0 ? dna_size[indiv_id] + k : k;


                                auto it_rna_end = internal_simd_struct[indiv_id]->terminator_lag.upper_bound(
                                        k);


                                if (it_rna_end ==
                                    internal_simd_struct[indiv_id]->terminator_lag.begin()) {
                                    it_rna_end = internal_simd_struct[indiv_id]->terminator_lag.end();
                                    it_rna_end--;
                                } else if ((*it_rna_end) != k)
                                    it_rna_end--;


                                /* {
                                    it_rna_end--;
                                }*/


                                int32_t rna_end =
                                        *it_rna_end - 10 < 0 ? dna_size[indiv_id] + (*it_rna_end - 10)
                                                             :
                                        *it_rna_end - 10;

                                /*if (indiv_id == 969 && AeTime::time() == 137) {
                                  auto it_rn = it_rna_end;
                                  it_rn++;
                                  printf("Looking for term from %d (start rna %d) : %d Computed end %d (next end %d)\n",k,internal_simd_struct[indiv_id]->promoters[rna_idx]->pos,
                                         *it_rna_end,rna_end,*it_rn);
                                }*/
                                int32_t rna_length = 0;

                                if (internal_simd_struct[indiv_id]->promoters[rna_idx]->pos <
                                    rna_end)
                                    rna_length =
                                            internal_simd_struct[indiv_id]->promoters[rna_idx]->pos +
                                            dna_size[indiv_id] - rna_end;
                                else
                                    rna_length =
                                            internal_simd_struct[indiv_id]->promoters[rna_idx]->pos -
                                            rna_end;

                                rna_length -= 21;

                                if (rna_length >= 0) {


                                    int glob_rna_idx = -1;
#pragma omp atomic capture
                                    {
                                        glob_rna_idx = internal_simd_struct[indiv_id]->rna_count_;
                                        internal_simd_struct[indiv_id]->rna_count_ =
                                                internal_simd_struct[indiv_id]->rna_count_ + 1;
                                    }

                                    internal_simd_struct[indiv_id]->rnas[glob_rna_idx] = new pRNA(
                                            internal_simd_struct[indiv_id]->promoters[rna_idx]->pos,
                                            rna_end,
                                            !internal_simd_struct[indiv_id]->promoters[rna_idx]->leading_or_lagging,
                                            1.0 -
                                            fabs(
                                                    ((float) internal_simd_struct[indiv_id]->promoters[rna_idx]->error)) /
                                            5.0, rna_length);

                                }
                            }
                        }
                    }
                }
            }
        }
    }
}


    void SIMD_Individual::compute_RNA(int indiv_id) {

//#pragma omp parallel for schedule(dynamic)
        {
                internal_simd_struct[indiv_id]->rnas.resize(
                        internal_simd_struct[indiv_id]->promoters.size());
//#pragma omp parallel for firstprivate(indiv_id) schedule(dynamic)
#ifdef WITH_FINETASKLOOP
#pragma omp taskloop grainsize(rna_grain_size)
#endif
                for (int rna_idx = 0; rna_idx <
                                      (int) internal_simd_struct[indiv_id]->promoters.size(); rna_idx++) {
//#pragma omp task firstprivate(indiv_id, rna_idx)
                    {
                        /*if (indiv_id == 345 && AeTime::time() == 47 &&
                            internal_simd_struct[indiv_id]->promoters[rna_idx]->pos == 4744) {
                          printf("Searching for an end with start pos %d LorL %d\n",
                                 internal_simd_struct[indiv_id]->promoters[rna_idx]->pos,
                                 internal_simd_struct[indiv_id]->promoters[rna_idx]->leading_or_lagging);
                        }*/

                        if (internal_simd_struct[indiv_id]->promoters[rna_idx] != nullptr) {
                            if (internal_simd_struct[indiv_id]->promoters[rna_idx]->leading_or_lagging) {
                                if (internal_simd_struct[indiv_id]->terminator_lead.size() != 0) {


                                    int k =
                                            internal_simd_struct[indiv_id]->promoters[rna_idx]->pos + 22;
                                    k = k >= dna_size[indiv_id] ? k - dna_size[indiv_id] : k;

/*        if ((indiv_id == 309 && AeTime::time() == 105) ||
            (indiv_id == 915 && AeTime::time() == 26 && rna_idx == 11)) {
          printf("Looking at %d\n",k);
        }*/

                                    auto it_rna_end = internal_simd_struct[indiv_id]->terminator_lead.lower_bound(
                                            k);

                                    /*if ((indiv_id == 309 && AeTime::time() == 105) ||
                                        (indiv_id == 915 && AeTime::time() == 26 && rna_idx == 11)) {
                                      printf("Looking at %d\n",*it_rna_end);
                                    }*/


                                    if (it_rna_end ==
                                        internal_simd_struct[indiv_id]->terminator_lead.end()) {
                                        it_rna_end = internal_simd_struct[indiv_id]->terminator_lead.begin();

                                        /* if ((indiv_id == 309 && AeTime::time() == 105) ||
                                             (indiv_id == 915 && AeTime::time() == 26 && rna_idx == 11)) {
                                           printf("Found at %d (restart because end)\n",*it_rna_end);
                                         }*/
                                    }

                                    /*if ((indiv_id == 309 && AeTime::time() == 105) ||
                                        (indiv_id == 915 && AeTime::time() == 26 && rna_idx == 11)) {
                                      auto tmp_rna_end = it_rna_end;
                                      printf("Terminators for %d current %d : ",rna_idx,*it_rna_end);
                                      while (tmp_rna_end != internal_simd_struct[indiv_id]->terminator_lead.end()) {
                                        printf("%d ",(*tmp_rna_end));
                                        tmp_rna_end++;
                                      }
                                      printf("\n");
                                    }*/

                                    int32_t rna_end =
                                            *it_rna_end + 10 >= dna_size[indiv_id] ?
                                            *it_rna_end + 10 - dna_size[indiv_id] :
                                            *it_rna_end + 10;

                                    /*if (indiv_id == 309 && AeTime::time() == 105) {
                                      printf("Looking for term from %d (start rna %d) : %d Computed end %d\n",k,internal_simd_struct[indiv_id]->promoters[rna_idx]->pos,
                                             *it_rna_end,rna_end);
                                    }*/
                                    int32_t rna_length = 0;

                                    if (internal_simd_struct[indiv_id]->promoters[rna_idx]->pos
                                        > rna_end)
                                        rna_length = dna_size[indiv_id] -
                                                     internal_simd_struct[indiv_id]->promoters[rna_idx]->pos
                                                     + rna_end;
                                    else
                                        rna_length = rna_end - internal_simd_struct[indiv_id]->
                                                promoters[rna_idx]->pos;

                                    rna_length -= 21;

                                    if (rna_length >= 0) {


                                        int glob_rna_idx = -1;
#pragma omp atomic capture
                                        {
                                            glob_rna_idx = internal_simd_struct[indiv_id]->rna_count_;
                                            internal_simd_struct[indiv_id]->rna_count_ =
                                                    internal_simd_struct[indiv_id]->rna_count_ + 1;
                                        }

                                        internal_simd_struct[indiv_id]->rnas[glob_rna_idx] = new pRNA(
                                                internal_simd_struct[indiv_id]->promoters[rna_idx]->pos,
                                                rna_end,
                                                !internal_simd_struct[indiv_id]->promoters[rna_idx]->leading_or_lagging,
                                                1.0 -
                                                fabs(
                                                        ((float) internal_simd_struct[indiv_id]->promoters[rna_idx]->error)) /
                                                5.0, rna_length);
                                    }
                                }
                            } else {
                                // LAGGING
                                if (internal_simd_struct[indiv_id]->terminator_lag.size() != 0) {




                                    // Search for terminator
                                    int k =
                                            internal_simd_struct[indiv_id]->promoters[rna_idx]->pos - 22;
                                    k = k < 0 ? dna_size[indiv_id] + k : k;


                                    auto it_rna_end = internal_simd_struct[indiv_id]->terminator_lag.upper_bound(
                                            k);


                                    if (it_rna_end ==
                                        internal_simd_struct[indiv_id]->terminator_lag.begin()) {
                                        it_rna_end = internal_simd_struct[indiv_id]->terminator_lag.end();
                                        it_rna_end--;
                                    } else if ((*it_rna_end) != k)
                                        it_rna_end--;


                                    /* {
                                        it_rna_end--;
                                    }*/


                                    int32_t rna_end =
                                            *it_rna_end - 10 < 0 ? dna_size[indiv_id] + (*it_rna_end - 10)
                                                                 :
                                            *it_rna_end - 10;

                                    /*if (indiv_id == 969 && AeTime::time() == 137) {
                                      auto it_rn = it_rna_end;
                                      it_rn++;
                                      printf("Looking for term from %d (start rna %d) : %d Computed end %d (next end %d)\n",k,internal_simd_struct[indiv_id]->promoters[rna_idx]->pos,
                                             *it_rna_end,rna_end,*it_rn);
                                    }*/
                                    int32_t rna_length = 0;

                                    if (internal_simd_struct[indiv_id]->promoters[rna_idx]->pos <
                                        rna_end)
                                        rna_length =
                                                internal_simd_struct[indiv_id]->promoters[rna_idx]->pos +
                                                dna_size[indiv_id] - rna_end;
                                    else
                                        rna_length =
                                                internal_simd_struct[indiv_id]->promoters[rna_idx]->pos -
                                                rna_end;

                                    rna_length -= 21;

                                    if (rna_length >= 0) {


                                        int glob_rna_idx = -1;
#pragma omp atomic capture
                                        {
                                            glob_rna_idx = internal_simd_struct[indiv_id]->rna_count_;
                                            internal_simd_struct[indiv_id]->rna_count_ =
                                                    internal_simd_struct[indiv_id]->rna_count_ + 1;
                                        }

                                        internal_simd_struct[indiv_id]->rnas[glob_rna_idx] = new pRNA(
                                                internal_simd_struct[indiv_id]->promoters[rna_idx]->pos,
                                                rna_end,
                                                !internal_simd_struct[indiv_id]->promoters[rna_idx]->leading_or_lagging,
                                                1.0 -
                                                fabs(
                                                        ((float) internal_simd_struct[indiv_id]->promoters[rna_idx]->error)) /
                                                5.0, rna_length);

                                    }
                                }
                            }
                        }
                    }
                }
            }

    }

    void SIMD_Individual::start_protein(int indiv_id) {
#ifdef WITH_FINETASKLOOP
                    #pragma omp taskloop grainsize(rna_grain_size)
#endif
        for (int rna_idx = 0; rna_idx <
                              (int) internal_simd_struct[indiv_id]->rna_count_; rna_idx++) {
//#pragma omp task firstprivate(indiv_id, rna_idx) depend(out: internal_simd_struct[indiv_id])
            {
                if (internal_simd_struct[indiv_id]->rnas[rna_idx]->is_init_) {
                    int x = indiv_id / exp_m_->world()->height();
                    int y = indiv_id % exp_m_->world()->height();

                    int c_pos = internal_simd_struct[indiv_id]->rnas[rna_idx]->begin;

//      printf("Searching for proteins in %d of indiv %d\n",rna_idx,indiv_id);

                    if (internal_simd_struct[indiv_id]->rnas[rna_idx]->length >= 22) {
                        if (internal_simd_struct[indiv_id]->rnas[rna_idx]->leading_lagging ==
                            0) {
                            c_pos += 22;
                            c_pos =
                                    c_pos >= dna_size[indiv_id] ? c_pos - dna_size[indiv_id]
                                                                : c_pos;
                        } else {
                            c_pos -= 22;
                            c_pos =
                                    c_pos < 0 ? ((int) dna_size[indiv_id]) + c_pos : c_pos;
                        }

/*        if (indiv_id == 601 && AeTime::time() >= 322)
          printf("Search for protein at %d -> %d -- RNA IDX %d (LEADING/LAGGING %d)\n",c_pos,
                 internal_simd_struct[indiv_id]->rnas[rna_idx].end,
                 rna_idx,
                 internal_simd_struct[indiv_id]->rnas[rna_idx].leading_lagging);*/


                        while (c_pos !=
                               internal_simd_struct[indiv_id]->rnas[rna_idx]->end) {
                            bool start = false;
                            int t_pos, k_t;

                            if (internal_simd_struct[indiv_id]->rnas[rna_idx]->leading_lagging ==
                                0) {
                                // Search for Shine Dalgarro + START codon on LEADING
#ifdef WITH_BITSET
                                start = internal_simd_struct[indiv_id]->dna_->bitset_->is_shine_dalgarno_protein_start(
                true, c_pos);
#else

//#pragma omp simd aligned(internal_simd_struct[indiv_id]->dna_->data_, SHINE_DAL_SEQ_LEAD:64)
                                for (int k = 0; k < 9; k++) {
                                    k_t = k >= 6 ? k + 4 : k;
                                    t_pos = c_pos + k_t >= dna_size[indiv_id] ? c_pos + k_t -
                                                                                dna_size[indiv_id]
                                                                              :
                                            c_pos + k_t;

                                    if (internal_simd_struct[indiv_id]->dna_->data_[t_pos] ==
                                        SHINE_DAL_SEQ_LEAD[k]) {
                                        start = true;
                                    } else {
                                        start = false;
                                        break;
                                    }
                                }
#endif
                            } else {

                                // Search for Shine Dalgarro + START codon on LAGGING
#ifdef WITH_BITSET
                                start = internal_simd_struct[indiv_id]->dna_->bitset_->is_shine_dalgarno_protein_start(
                false, c_pos);
#else
//#pragma omp simd aligned(internal_simd_struct[indiv_id]->dna_->data_, SHINE_DAL_SEQ_LAG:64)
                                for (int k = 0; k < 9; k++) {
                                    k_t = k >= 6 ? k + 4 : k;
                                    t_pos =
                                            c_pos - k_t < 0 ? dna_size[indiv_id] + (c_pos - k_t) :
                                            c_pos - k_t;

                                    if (internal_simd_struct[indiv_id]->dna_->data_[t_pos] ==
                                        SHINE_DAL_SEQ_LAG[k]) {
                                        start = true;
                                    } else {
                                        start = false;
                                        break;
                                    }
                                }
#endif
                            }

//                  if (indiv_id == 905)
//                    printf("Search for start prot at %d : %d (%d)\n",c_pos,start,k_t);

                            if (start) {
                                /*if (indiv_id == 601 && AeTime::time() == 323 && internal_simd_struct[indiv_id]->rnas[rna_idx].leading_lagging == true)
                                  printf("Found Start LAG POS %d\n",c_pos);*/

                                //printf("Start protein %d\n",c_pos);

                                //if (indiv_id == 905)
                                //printf("Found start prot at %d : %d (%d)\n",c_pos,start,k_t);
#pragma omp critical
                                {
                                    internal_simd_struct[indiv_id]->rnas[rna_idx]->start_prot.
                                            push_back(c_pos);
                                }
                            }

                            if (internal_simd_struct[indiv_id]->rnas[rna_idx]->leading_lagging ==
                                0) {
                                c_pos++;
                                c_pos =
                                        c_pos >= dna_size[indiv_id] ? c_pos - dna_size[indiv_id]
                                                                    : c_pos;
                            } else {
                                c_pos--;
                                c_pos = c_pos < 0 ? dna_size[indiv_id] + c_pos : c_pos;
                            }

//          if (indiv_id == 180 && AeTime::time() == 2) {exit(-1);printf("Searching at %d (END %d LENGTH %ld DNASIZE %ld)\n",
//                                      c_pos,internal_simd_struct[indiv_id]->rnas[rna_idx].end,
//                                      internal_simd_struct[indiv_id]->rnas[rna_idx].length,
//                                      dna_size[indiv_id]);}
                        }
                    }
                }
            }
        }
    }


void SIMD_Individual::start_protein() {
  //int x, y;
//#pragma omp parallel for schedule(dynamic)

    for (int indiv_id = 0; indiv_id < (int) exp_m_->nb_indivs(); indiv_id++) {
      if (exp_m_->dna_mutator_array_[indiv_id]->hasMutate()) {
//#pragma omp parallel for firstprivate(indiv_id) schedule(dynamic)
        for (int rna_idx = 0; rna_idx <
                              (int) internal_simd_struct[indiv_id]->rna_count_; rna_idx++) {
//#pragma omp task firstprivate(indiv_id, rna_idx) depend(inout: internal_simd_struct[indiv_id])
          {
            if (internal_simd_struct[indiv_id]->rnas[rna_idx]->is_init_) {
              int x = indiv_id / exp_m_->world()->height();
              int y = indiv_id % exp_m_->world()->height();

              int c_pos = internal_simd_struct[indiv_id]->rnas[rna_idx]->begin;

//      printf("Searching for proteins in %d of indiv %d\n",rna_idx,indiv_id);

              if (internal_simd_struct[indiv_id]->rnas[rna_idx]->length >= 22) {
                if (internal_simd_struct[indiv_id]->rnas[rna_idx]->leading_lagging ==
                    0) {
                  c_pos += 22;
                  c_pos =
                      c_pos >= dna_size[indiv_id] ? c_pos - dna_size[indiv_id]
                                                  : c_pos;
                } else {
                  c_pos -= 22;
                  c_pos =
                      c_pos < 0 ? ((int) dna_size[indiv_id]) + c_pos : c_pos;
                }

/*        if (indiv_id == 601 && AeTime::time() >= 322)
          printf("Search for protein at %d -> %d -- RNA IDX %d (LEADING/LAGGING %d)\n",c_pos,
                 internal_simd_struct[indiv_id]->rnas[rna_idx].end,
                 rna_idx,
                 internal_simd_struct[indiv_id]->rnas[rna_idx].leading_lagging);*/


                while (c_pos !=
                       internal_simd_struct[indiv_id]->rnas[rna_idx]->end) {
                  bool start = false;
                  int t_pos, k_t;

                  if (internal_simd_struct[indiv_id]->rnas[rna_idx]->leading_lagging ==
                      0) {
                    // Search for Shine Dalgarro + START codon on LEADING
#ifdef WITH_BITSET
                    start = internal_simd_struct[indiv_id]->dna_->bitset_->is_shine_dalgarno_protein_start(
                        true, c_pos);
#else

//#pragma omp simd aligned(internal_simd_struct[indiv_id]->dna_->data_, SHINE_DAL_SEQ_LEAD:64)
                    for (int k = 0; k < 9; k++) {
                      k_t = k >= 6 ? k + 4 : k;
                      t_pos = c_pos + k_t >= dna_size[indiv_id] ? c_pos + k_t -
                                                                  dna_size[indiv_id]
                                                                :
                              c_pos + k_t;

                      if (internal_simd_struct[indiv_id]->dna_->data_[t_pos] ==
                          SHINE_DAL_SEQ_LEAD[k]) {
                        start = true;
                      } else {
                        start = false;
                        break;
                      }
                    }
#endif
                  } else {

                    // Search for Shine Dalgarro + START codon on LAGGING
#ifdef WITH_BITSET
                    start = internal_simd_struct[indiv_id]->dna_->bitset_->is_shine_dalgarno_protein_start(
                        false, c_pos);
#else
//#pragma omp simd aligned(internal_simd_struct[indiv_id]->dna_->data_, SHINE_DAL_SEQ_LAG:64)
                    for (int k = 0; k < 9; k++) {
                      k_t = k >= 6 ? k + 4 : k;
                      t_pos =
                          c_pos - k_t < 0 ? dna_size[indiv_id] + (c_pos - k_t) :
                          c_pos - k_t;

                      if (internal_simd_struct[indiv_id]->dna_->data_[t_pos] ==
                          SHINE_DAL_SEQ_LAG[k]) {
                        start = true;
                      } else {
                        start = false;
                        break;
                      }
                    }
#endif
                  }

//                  if (indiv_id == 905)
//                    printf("Search for start prot at %d : %d (%d)\n",c_pos,start,k_t);

                  if (start) {
                    /*if (indiv_id == 601 && AeTime::time() == 323 && internal_simd_struct[indiv_id]->rnas[rna_idx].leading_lagging == true)
                      printf("Found Start LAG POS %d\n",c_pos);*/

                    //printf("Start protein %d\n",c_pos);

                    //if (indiv_id == 905)
                      //printf("Found start prot at %d : %d (%d)\n",c_pos,start,k_t);

                    internal_simd_struct[indiv_id]->rnas[rna_idx]->start_prot.
                        push_back(c_pos);
                  }

                  if (internal_simd_struct[indiv_id]->rnas[rna_idx]->leading_lagging ==
                      0) {
                    c_pos++;
                    c_pos =
                        c_pos >= dna_size[indiv_id] ? c_pos - dna_size[indiv_id]
                                                    : c_pos;
                  } else {
                    c_pos--;
                    c_pos = c_pos < 0 ? dna_size[indiv_id] + c_pos : c_pos;
                  }

//          if (indiv_id == 180 && AeTime::time() == 2) {exit(-1);printf("Searching at %d (END %d LENGTH %ld DNASIZE %ld)\n",
//                                      c_pos,internal_simd_struct[indiv_id]->rnas[rna_idx].end,
//                                      internal_simd_struct[indiv_id]->rnas[rna_idx].length,
//                                      dna_size[indiv_id]);}
                }
              }
            }
          }
        }
      }
    }
 // }
}


void SIMD_Individual::compute_protein(int indiv_id) {
//#pragma omp task firstprivate(indiv_id) depend(inout: internal_simd_struct[indiv_id])
    {
        int resize_to = 0;
//#pragma omp taskloop
        for (int rna_idx = 0; rna_idx <
                              (int) internal_simd_struct[indiv_id]->rna_count_; rna_idx++) {
            if (internal_simd_struct[indiv_id]->rnas[rna_idx]->is_init_)
                resize_to += internal_simd_struct[indiv_id]->
                        rnas[rna_idx]->start_prot.size();
        }
        internal_simd_struct[indiv_id]->
                proteins.resize(resize_to);
    }
//#pragma omp parallel for firstprivate(indiv_id) schedule(dynamic)
#ifdef WITH_FINETASKLOOP
#pragma omp taskloop grainsize(rna_grain_size)
#endif
                for (int rna_idx = 0; rna_idx <
                                      (int) internal_simd_struct[indiv_id]->rna_count_; rna_idx++) {
                    if (internal_simd_struct[indiv_id]->rnas[rna_idx]->is_init_) {
//#pragma omp parallel for firstprivate(indiv_id, rna_idx) schedule(dynamic)
                        for (int protein_idx = 0;
                             protein_idx < (int) internal_simd_struct[indiv_id]->
                                     rnas[rna_idx]->start_prot.size(); protein_idx++) {
//#pragma omp task firstprivate(indiv_id, rna_idx, protein_idx) depend(inout: internal_simd_struct[indiv_id])
                            {
                                int x = indiv_id / exp_m_->world()->height();
                                int y = indiv_id % exp_m_->world()->height();

                                int start_protein_pos = internal_simd_struct[indiv_id]->
                                        rnas[rna_idx]->leading_lagging == 0 ?
                                                        internal_simd_struct[indiv_id]->
                                                                rnas[rna_idx]->start_prot[protein_idx] +
                                                        13 :
                                                        internal_simd_struct[indiv_id]->
                                                                rnas[rna_idx]->start_prot[protein_idx] -
                                                        13;
                                int length;

                                if (internal_simd_struct[indiv_id]->
                                        rnas[rna_idx]->leading_lagging == 0) {
                                    start_protein_pos = start_protein_pos >= dna_size[indiv_id] ?
                                                        start_protein_pos - dna_size[indiv_id]
                                                                                                : start_protein_pos;

                                    if (internal_simd_struct[indiv_id]->
                                            rnas[rna_idx]->start_prot[protein_idx] <
                                        internal_simd_struct[indiv_id]->
                                                rnas[rna_idx]->end) {
                                        length = internal_simd_struct[indiv_id]->
                                                rnas[rna_idx]->end -
                                                 internal_simd_struct[indiv_id]->
                                                         rnas[rna_idx]->start_prot[protein_idx];
                                    } else {
                                        length = dna_size[indiv_id] -
                                                 internal_simd_struct[indiv_id]->
                                                         rnas[rna_idx]->start_prot[protein_idx] +
                                                 internal_simd_struct[indiv_id]->
                                                         rnas[rna_idx]->end;

                                    }

                                    length -= 13;
                                } else {


                                    start_protein_pos = start_protein_pos < 0 ?
                                                        dna_size[indiv_id] + start_protein_pos
                                                                              : start_protein_pos;

                                    if (internal_simd_struct[indiv_id]->
                                            rnas[rna_idx]->start_prot[protein_idx] >
                                        internal_simd_struct[indiv_id]->
                                                rnas[rna_idx]->end) {
                                        length = internal_simd_struct[indiv_id]->
                                                rnas[rna_idx]->start_prot[protein_idx] -
                                                 internal_simd_struct[indiv_id]->
                                                         rnas[rna_idx]->end;
                                    } else {
                                        length = internal_simd_struct[indiv_id]->
                                                rnas[rna_idx]->start_prot[protein_idx] +
                                                 dna_size[indiv_id] -
                                                 internal_simd_struct[indiv_id]->
                                                         rnas[rna_idx]->end;
                                    }

                                    length -= 13;

                                    /*if (indiv_id == 107 && AeTime::time() == 6 && internal_simd_struct[indiv_id]->rnas[rna_idx].leading_lagging == true)
                                      printf("Found for Start Prot for LAG RNA IDX %d (%d) at %d SEARCHING LENGHT %d RNA END %d\n",rna_idx,protein_idx,internal_simd_struct[indiv_id]->
                                          rnas[rna_idx].start_prot[protein_idx],length,internal_simd_struct[indiv_id]->
                                          rnas[rna_idx].end);*/
                                }

                                bool is_protein = false;

                                //length -= 2;
                                length += 1;
                                length = length - (length % 3);

                                int j = 0;
                                int32_t transcribed_start = 0;

                                if (internal_simd_struct[indiv_id]->
                                        rnas[rna_idx]->leading_lagging == 0) {
                                    transcribed_start = internal_simd_struct[indiv_id]->
                                            rnas[rna_idx]->begin + 22;
                                    transcribed_start = transcribed_start >= dna_size[indiv_id] ?
                                                        transcribed_start - dna_size[indiv_id]
                                                                                                : transcribed_start;

                                    if (transcribed_start <= internal_simd_struct[indiv_id]->
                                            rnas[rna_idx]->start_prot[protein_idx]) {
                                        j = internal_simd_struct[indiv_id]->
                                                rnas[rna_idx]->start_prot[protein_idx] -
                                            transcribed_start;
                                    } else {
                                        j = dna_size[indiv_id] -
                                            transcribed_start +
                                            internal_simd_struct[indiv_id]->
                                                    rnas[rna_idx]->start_prot[protein_idx];

                                    }
                                } else {
                                    transcribed_start = internal_simd_struct[indiv_id]->
                                            rnas[rna_idx]->begin - 22;
                                    transcribed_start = transcribed_start < 0 ?
                                                        dna_size[indiv_id] + transcribed_start
                                                                              : transcribed_start;

                                    if (transcribed_start >=
                                        internal_simd_struct[indiv_id]->
                                                rnas[rna_idx]->start_prot[protein_idx]) {
                                        j = transcribed_start -
                                            internal_simd_struct[indiv_id]->
                                                    rnas[rna_idx]->start_prot[protein_idx];
                                    } else {
                                        j = transcribed_start +
                                            dna_size[indiv_id] - internal_simd_struct[indiv_id]->
                                                rnas[rna_idx]->start_prot[protein_idx];
                                    }
                                }

                                j += 13;

                                /*if (indiv_id == 906 && AeTime::time() == 69)
                                  printf("Length %d j %d DNA Size %d start prot %d start %d stop %d transcribed_start %d\n",internal_simd_struct[indiv_id]->
                                      rnas[rna_idx].length,j,dna_size[indiv_id],internal_simd_struct[indiv_id]->
                                      rnas[rna_idx].start_prot[protein_idx],internal_simd_struct[indiv_id]->
                                      rnas[rna_idx].begin,internal_simd_struct[indiv_id]->
                                      rnas[rna_idx].end,transcribed_start);*/

                                while (internal_simd_struct[indiv_id]->
                                        rnas[rna_idx]->length - j >= 3) {

                                    int t_k;

                                    /*if (indiv_id == 906 && AeTime::time() == 69)
                                      printf("Length %d j %d DNA Size %d start prot %d (%d) start %d stop %d\n",internal_simd_struct[indiv_id]->
                                          rnas[rna_idx].length,j,dna_size[indiv_id],internal_simd_struct[indiv_id]->
                                          rnas[rna_idx].start_prot[protein_idx],start_protein_pos,internal_simd_struct[indiv_id]->
                                          rnas[rna_idx].begin,internal_simd_struct[indiv_id]->
                                          rnas[rna_idx].end);*/

                                    if (internal_simd_struct[indiv_id]->
                                            rnas[rna_idx]->leading_lagging == 0) {
                                        start_protein_pos =
                                                start_protein_pos >= dna_size[indiv_id] ?
                                                start_protein_pos - dna_size[indiv_id]
                                                                                        : start_protein_pos;
#ifdef WITH_BITSET
                                        is_protein = internal_simd_struct[indiv_id]->dna_->bitset_->is_protein_stop(
                        true, start_protein_pos);
                    t_k = start_protein_pos + 2 >= dna_size[indiv_id] ?
                          start_protein_pos - dna_size[indiv_id] + 2 :
                          start_protein_pos + 2;
#else
                                        is_protein = false;

                                        for (int k = 0; k < 3; k++) {
                                            t_k = start_protein_pos + k >= dna_size[indiv_id] ?
                                                  start_protein_pos - dna_size[indiv_id] + k :
                                                  start_protein_pos + k;

                                            if (internal_simd_struct[indiv_id]->dna_->data_[t_k] ==
                                                PROTEIN_END_LEAD[k]) {
                                                is_protein = true;
                                            } else {
                                                is_protein = false;
                                                break;
                                            }
                                        }
#endif

                                        if (is_protein) {
                                            int prot_length = -1;

                                            if (internal_simd_struct[indiv_id]->
                                                    rnas[rna_idx]->start_prot[protein_idx] + 13 < t_k) {
                                                prot_length = t_k -
                                                              (internal_simd_struct[indiv_id]->
                                                                      rnas[rna_idx]->start_prot[protein_idx] +
                                                               13);
                                            } else {
                                                prot_length = dna_size[indiv_id] -
                                                              (internal_simd_struct[indiv_id]->
                                                                      rnas[rna_idx]->start_prot[protein_idx] +
                                                               13) + t_k;
                                            }

                                            if (prot_length >= 3) {
                                                int32_t glob_prot_idx = -1;
#pragma omp atomic capture
                                                {
                                                    glob_prot_idx = internal_simd_struct[indiv_id]->protein_count_;
                                                    internal_simd_struct[indiv_id]->protein_count_ =
                                                            internal_simd_struct[indiv_id]->protein_count_ +
                                                            1;
                                                }

//#pragma omp critical
                                                {
                                                    internal_simd_struct[indiv_id]->
                                                            proteins[glob_prot_idx] = new pProtein(
                                                            internal_simd_struct[indiv_id]->
                                                                    rnas[rna_idx]->start_prot[protein_idx], t_k,
                                                            prot_length,
                                                            internal_simd_struct[indiv_id]->rnas[rna_idx]->leading_lagging,
                                                            internal_simd_struct[indiv_id]->rnas[rna_idx]->e
                                                    );
                                                }

                                                internal_simd_struct[indiv_id]->
                                                        rnas[rna_idx]->is_coding_ = true;
                                            }/* else if (indiv_id == 906)
                  printf("Length %d j %d DNA Size %d start prot %d end prot %d start %d stop %d\n",internal_simd_struct[indiv_id]->
                      rnas[rna_idx].length,j,dna_size[indiv_id],internal_simd_struct[indiv_id]->
                      rnas[rna_idx].start_prot[protein_idx],t_k,
                         internal_simd_struct[indiv_id]->
                      rnas[rna_idx].begin,internal_simd_struct[indiv_id]->
                      rnas[rna_idx].end);*/
                                            break;
                                        }

                                        start_protein_pos += 3;
                                        start_protein_pos =
                                                start_protein_pos >= dna_size[indiv_id] ?
                                                start_protein_pos - dna_size[indiv_id]
                                                                                        : start_protein_pos;
                                    } else {


                                        start_protein_pos = start_protein_pos < 0 ?
                                                            dna_size[indiv_id] + start_protein_pos
                                                                                  : start_protein_pos;

#ifdef WITH_BITSET
                                        is_protein = internal_simd_struct[indiv_id]->dna_->bitset_->is_protein_stop(
                        false, start_protein_pos);
                    t_k = start_protein_pos - 2 < 0 ?
                          dna_size[indiv_id] + (start_protein_pos - 2) :
                          start_protein_pos - 2;
#else
                                        is_protein = false;
                                        for (int k = 0; k < 3; k++) {
                                            t_k = start_protein_pos - k < 0 ?
                                                  dna_size[indiv_id] + (start_protein_pos - k) :
                                                  start_protein_pos - k;

                                            if (internal_simd_struct[indiv_id]->dna_->data_[t_k] ==
                                                PROTEIN_END_LAG[k]) {
                                                is_protein = true;
                                            } else {
                                                is_protein = false;
                                                break;
                                            }
                                        }
#endif

                                        if (is_protein) {
                                            int prot_length = -1;
                                            if (internal_simd_struct[indiv_id]->
                                                    rnas[rna_idx]->start_prot[protein_idx] - 13 > t_k) {
                                                prot_length =
                                                        (internal_simd_struct[indiv_id]->
                                                                rnas[rna_idx]->start_prot[protein_idx] - 13) -
                                                        t_k;
                                            } else {
                                                prot_length =
                                                        (internal_simd_struct[indiv_id]->
                                                                rnas[rna_idx]->start_prot[protein_idx] - 13) +
                                                        dna_size[indiv_id] - t_k;
                                            }
                                            if (prot_length >= 3) {
                                                int32_t glob_prot_idx = -1;
#pragma omp atomic capture
                                                {
                                                    glob_prot_idx = internal_simd_struct[indiv_id]->protein_count_;
                                                    internal_simd_struct[indiv_id]->protein_count_ =
                                                            internal_simd_struct[indiv_id]->protein_count_ +
                                                            1;
                                                }

//#pragma omp critical
                                                {
                                                    internal_simd_struct[indiv_id]->
                                                            proteins[glob_prot_idx] = new pProtein(
                                                            internal_simd_struct[indiv_id]->
                                                                    rnas[rna_idx]->start_prot[protein_idx], t_k,
                                                            prot_length,
                                                            internal_simd_struct[indiv_id]->rnas[rna_idx]->leading_lagging,
                                                            internal_simd_struct[indiv_id]->rnas[rna_idx]->e
                                                    );
                                                }
//#pragma omp atomic
                                                internal_simd_struct[indiv_id]->
                                                        rnas[rna_idx]->is_coding_ = true;
                                            } /*else if (indiv_id == 906)
                printf("Length %d j %d DNA Size %d start prot %d end prot %d start %d stop %d\n",internal_simd_struct[indiv_id]->
                           rnas[rna_idx].length,j,dna_size[indiv_id],internal_simd_struct[indiv_id]->
                           rnas[rna_idx].start_prot[protein_idx],t_k,
                       internal_simd_struct[indiv_id]->
                           rnas[rna_idx].begin,internal_simd_struct[indiv_id]->
                        rnas[rna_idx].end);*/

                                            break;
                                        }
                                        start_protein_pos = start_protein_pos - 3;
                                        start_protein_pos = start_protein_pos < 0 ?
                                                            dna_size[indiv_id] + start_protein_pos
                                                                                  : start_protein_pos;
                                    }
                                    j += 3;
                                }
                            }
                        }
                    }
                }
            }


void SIMD_Individual::compute_protein() {
//#pragma omp parallel for schedule(dynamic)
/*#pragma omp parallel
#pragma omp single
  {*/
    for (int indiv_id = 0; indiv_id < (int) exp_m_->nb_indivs(); indiv_id++) {
      if (exp_m_->dna_mutator_array_[indiv_id]->hasMutate()) {
        int resize_to = 0;
        for (int rna_idx = 0; rna_idx <
                              (int) internal_simd_struct[indiv_id]->rna_count_; rna_idx++) {
          if (internal_simd_struct[indiv_id]->rnas[rna_idx]->is_init_)
            resize_to += internal_simd_struct[indiv_id]->
                rnas[rna_idx]->start_prot.size();
        }
        internal_simd_struct[indiv_id]->
            proteins.resize(resize_to);

//#pragma omp parallel for firstprivate(indiv_id) schedule(dynamic)
        for (int rna_idx = 0; rna_idx <
                              (int) internal_simd_struct[indiv_id]->rna_count_; rna_idx++) {
          if (internal_simd_struct[indiv_id]->rnas[rna_idx]->is_init_) {
//#pragma omp parallel for firstprivate(indiv_id, rna_idx) schedule(dynamic)
            for (int protein_idx = 0;
                 protein_idx < (int) internal_simd_struct[indiv_id]->
                     rnas[rna_idx]->start_prot.size(); protein_idx++) {
//#pragma omp task firstprivate(indiv_id, rna_idx, protein_idx) depend(in: internal_simd_struct[indiv_id])
              {
                int x = indiv_id / exp_m_->world()->height();
                int y = indiv_id % exp_m_->world()->height();

                int start_protein_pos = internal_simd_struct[indiv_id]->
                    rnas[rna_idx]->leading_lagging == 0 ?
                                        internal_simd_struct[indiv_id]->
                                            rnas[rna_idx]->start_prot[protein_idx] +
                                        13 :
                                        internal_simd_struct[indiv_id]->
                                            rnas[rna_idx]->start_prot[protein_idx] -
                                        13;
                int length;

                if (internal_simd_struct[indiv_id]->
                    rnas[rna_idx]->leading_lagging == 0) {
                  start_protein_pos = start_protein_pos >= dna_size[indiv_id] ?
                                      start_protein_pos - dna_size[indiv_id]
                                                                              : start_protein_pos;

                  if (internal_simd_struct[indiv_id]->
                      rnas[rna_idx]->start_prot[protein_idx] <
                      internal_simd_struct[indiv_id]->
                          rnas[rna_idx]->end) {
                    length = internal_simd_struct[indiv_id]->
                        rnas[rna_idx]->end -
                             internal_simd_struct[indiv_id]->
                                 rnas[rna_idx]->start_prot[protein_idx];
                  } else {
                    length = dna_size[indiv_id] -
                             internal_simd_struct[indiv_id]->
                                 rnas[rna_idx]->start_prot[protein_idx] +
                             internal_simd_struct[indiv_id]->
                                 rnas[rna_idx]->end;

                  }

                  length -= 13;
                } else {


                  start_protein_pos = start_protein_pos < 0 ?
                                      dna_size[indiv_id] + start_protein_pos
                                                            : start_protein_pos;

                  if (internal_simd_struct[indiv_id]->
                      rnas[rna_idx]->start_prot[protein_idx] >
                      internal_simd_struct[indiv_id]->
                          rnas[rna_idx]->end) {
                    length = internal_simd_struct[indiv_id]->
                        rnas[rna_idx]->start_prot[protein_idx] -
                             internal_simd_struct[indiv_id]->
                                 rnas[rna_idx]->end;
                  } else {
                    length = internal_simd_struct[indiv_id]->
                        rnas[rna_idx]->start_prot[protein_idx] +
                             dna_size[indiv_id] -
                             internal_simd_struct[indiv_id]->
                                 rnas[rna_idx]->end;
                  }

                  length -= 13;

                  /*if (indiv_id == 107 && AeTime::time() == 6 && internal_simd_struct[indiv_id]->rnas[rna_idx].leading_lagging == true)
                    printf("Found for Start Prot for LAG RNA IDX %d (%d) at %d SEARCHING LENGHT %d RNA END %d\n",rna_idx,protein_idx,internal_simd_struct[indiv_id]->
                        rnas[rna_idx].start_prot[protein_idx],length,internal_simd_struct[indiv_id]->
                        rnas[rna_idx].end);*/
                }

                bool is_protein = false;

                //length -= 2;
                length += 1;
                length = length - (length % 3);

                int j = 0;
                int32_t transcribed_start = 0;

                if (internal_simd_struct[indiv_id]->
                    rnas[rna_idx]->leading_lagging == 0) {
                  transcribed_start = internal_simd_struct[indiv_id]->
                      rnas[rna_idx]->begin + 22;
                  transcribed_start = transcribed_start >= dna_size[indiv_id] ?
                                      transcribed_start - dna_size[indiv_id]
                                                                              : transcribed_start;

                  if (transcribed_start <= internal_simd_struct[indiv_id]->
                      rnas[rna_idx]->start_prot[protein_idx]) {
                    j = internal_simd_struct[indiv_id]->
                        rnas[rna_idx]->start_prot[protein_idx] -
                        transcribed_start;
                  } else {
                    j = dna_size[indiv_id] -
                        transcribed_start +
                        internal_simd_struct[indiv_id]->
                            rnas[rna_idx]->start_prot[protein_idx];

                  }
                } else {
                  transcribed_start = internal_simd_struct[indiv_id]->
                      rnas[rna_idx]->begin - 22;
                  transcribed_start = transcribed_start < 0 ?
                                      dna_size[indiv_id] + transcribed_start
                                                            : transcribed_start;

                  if (transcribed_start >=
                      internal_simd_struct[indiv_id]->
                          rnas[rna_idx]->start_prot[protein_idx]) {
                    j = transcribed_start -
                        internal_simd_struct[indiv_id]->
                            rnas[rna_idx]->start_prot[protein_idx];
                  } else {
                    j = transcribed_start +
                        dna_size[indiv_id] - internal_simd_struct[indiv_id]->
                        rnas[rna_idx]->start_prot[protein_idx];
                  }
                }

                j += 13;

                /*if (indiv_id == 906 && AeTime::time() == 69)
                  printf("Length %d j %d DNA Size %d start prot %d start %d stop %d transcribed_start %d\n",internal_simd_struct[indiv_id]->
                      rnas[rna_idx].length,j,dna_size[indiv_id],internal_simd_struct[indiv_id]->
                      rnas[rna_idx].start_prot[protein_idx],internal_simd_struct[indiv_id]->
                      rnas[rna_idx].begin,internal_simd_struct[indiv_id]->
                      rnas[rna_idx].end,transcribed_start);*/

                while (internal_simd_struct[indiv_id]->
                    rnas[rna_idx]->length - j >= 3) {

                  int t_k;

                  /*if (indiv_id == 906 && AeTime::time() == 69)
                    printf("Length %d j %d DNA Size %d start prot %d (%d) start %d stop %d\n",internal_simd_struct[indiv_id]->
                        rnas[rna_idx].length,j,dna_size[indiv_id],internal_simd_struct[indiv_id]->
                        rnas[rna_idx].start_prot[protein_idx],start_protein_pos,internal_simd_struct[indiv_id]->
                        rnas[rna_idx].begin,internal_simd_struct[indiv_id]->
                        rnas[rna_idx].end);*/

                  if (internal_simd_struct[indiv_id]->
                      rnas[rna_idx]->leading_lagging == 0) {
                    start_protein_pos =
                        start_protein_pos >= dna_size[indiv_id] ?
                        start_protein_pos - dna_size[indiv_id]
                                                                : start_protein_pos;
#ifdef WITH_BITSET
                    is_protein = internal_simd_struct[indiv_id]->dna_->bitset_->is_protein_stop(
                        true, start_protein_pos);
                    t_k = start_protein_pos + 2 >= dna_size[indiv_id] ?
                          start_protein_pos - dna_size[indiv_id] + 2 :
                          start_protein_pos + 2;
#else
                    is_protein = false;

                    for (int k = 0; k < 3; k++) {
                      t_k = start_protein_pos + k >= dna_size[indiv_id] ?
                            start_protein_pos - dna_size[indiv_id] + k :
                            start_protein_pos + k;

                      if (internal_simd_struct[indiv_id]->dna_->data_[t_k] ==
                          PROTEIN_END_LEAD[k]) {
                        is_protein = true;
                      } else {
                        is_protein = false;
                        break;
                      }
                    }
#endif

                    if (is_protein) {
                      int prot_length = -1;

                      if (internal_simd_struct[indiv_id]->
                          rnas[rna_idx]->start_prot[protein_idx] + 13 < t_k) {
                        prot_length = t_k -
                                      (internal_simd_struct[indiv_id]->
                                          rnas[rna_idx]->start_prot[protein_idx] +
                                       13);
                      } else {
                        prot_length = dna_size[indiv_id] -
                                      (internal_simd_struct[indiv_id]->
                                          rnas[rna_idx]->start_prot[protein_idx] +
                                       13) + t_k;
                      }

                      if (prot_length >= 3) {
                        int32_t glob_prot_idx = -1;
#pragma omp atomic capture
                        {
                          glob_prot_idx = internal_simd_struct[indiv_id]->protein_count_;
                          internal_simd_struct[indiv_id]->protein_count_ =
                              internal_simd_struct[indiv_id]->protein_count_ +
                              1;
                        }

                        internal_simd_struct[indiv_id]->
                            proteins[glob_prot_idx] = new pProtein(
                            internal_simd_struct[indiv_id]->
                                rnas[rna_idx]->start_prot[protein_idx], t_k,
                            prot_length,
                            internal_simd_struct[indiv_id]->rnas[rna_idx]->leading_lagging,
                            internal_simd_struct[indiv_id]->rnas[rna_idx]->e
                        );


                        internal_simd_struct[indiv_id]->
                            rnas[rna_idx]->is_coding_ = true;
                      }/* else if (indiv_id == 906)
                  printf("Length %d j %d DNA Size %d start prot %d end prot %d start %d stop %d\n",internal_simd_struct[indiv_id]->
                      rnas[rna_idx].length,j,dna_size[indiv_id],internal_simd_struct[indiv_id]->
                      rnas[rna_idx].start_prot[protein_idx],t_k,
                         internal_simd_struct[indiv_id]->
                      rnas[rna_idx].begin,internal_simd_struct[indiv_id]->
                      rnas[rna_idx].end);*/
                      break;
                    }

                    start_protein_pos += 3;
                    start_protein_pos =
                        start_protein_pos >= dna_size[indiv_id] ?
                        start_protein_pos - dna_size[indiv_id]
                                                                : start_protein_pos;
                  } else {


                    start_protein_pos = start_protein_pos < 0 ?
                                        dna_size[indiv_id] + start_protein_pos
                                                              : start_protein_pos;

#ifdef WITH_BITSET
                    is_protein = internal_simd_struct[indiv_id]->dna_->bitset_->is_protein_stop(
                        false, start_protein_pos);
                    t_k = start_protein_pos - 2 < 0 ?
                          dna_size[indiv_id] + (start_protein_pos - 2) :
                          start_protein_pos - 2;
#else
                    is_protein = false;
                    for (int k = 0; k < 3; k++) {
                      t_k = start_protein_pos - k < 0 ?
                            dna_size[indiv_id] + (start_protein_pos - k) :
                            start_protein_pos - k;

                      if (internal_simd_struct[indiv_id]->dna_->data_[t_k] ==
                          PROTEIN_END_LAG[k]) {
                        is_protein = true;
                      } else {
                        is_protein = false;
                        break;
                      }
                    }
#endif

                    if (is_protein) {
                      int prot_length = -1;
                      if (internal_simd_struct[indiv_id]->
                          rnas[rna_idx]->start_prot[protein_idx] - 13 > t_k) {
                        prot_length =
                            (internal_simd_struct[indiv_id]->
                                rnas[rna_idx]->start_prot[protein_idx] - 13) -
                            t_k;
                      } else {
                        prot_length =
                            (internal_simd_struct[indiv_id]->
                                rnas[rna_idx]->start_prot[protein_idx] - 13) +
                            dna_size[indiv_id] - t_k;
                      }
                      if (prot_length >= 3) {
                        int32_t glob_prot_idx = -1;
#pragma omp atomic capture
                        {
                          glob_prot_idx = internal_simd_struct[indiv_id]->protein_count_;
                          internal_simd_struct[indiv_id]->protein_count_ =
                              internal_simd_struct[indiv_id]->protein_count_ +
                              1;
                        }

                        internal_simd_struct[indiv_id]->
                            proteins[glob_prot_idx] = new pProtein(
                            internal_simd_struct[indiv_id]->
                                rnas[rna_idx]->start_prot[protein_idx], t_k,
                            prot_length,
                            internal_simd_struct[indiv_id]->rnas[rna_idx]->leading_lagging,
                            internal_simd_struct[indiv_id]->rnas[rna_idx]->e
                        );
//#pragma omp atomic
                        internal_simd_struct[indiv_id]->
                            rnas[rna_idx]->is_coding_ = true;
                      } /*else if (indiv_id == 906)
                printf("Length %d j %d DNA Size %d start prot %d end prot %d start %d stop %d\n",internal_simd_struct[indiv_id]->
                           rnas[rna_idx].length,j,dna_size[indiv_id],internal_simd_struct[indiv_id]->
                           rnas[rna_idx].start_prot[protein_idx],t_k,
                       internal_simd_struct[indiv_id]->
                           rnas[rna_idx].begin,internal_simd_struct[indiv_id]->
                        rnas[rna_idx].end);*/

                      break;
                    }
                    start_protein_pos = start_protein_pos - 3;
                    start_protein_pos = start_protein_pos < 0 ?
                                        dna_size[indiv_id] + start_protein_pos
                                                              : start_protein_pos;
                  }
                  j += 3;
                }
              }
            }
          }
        }
      }
    }
//  }
}

void SIMD_Individual::translate_protein(double w_max) {
//#pragma omp parallel for schedule(dynamic)
/*#pragma omp parallel
#pragma omp single
  {*/
    for (int indiv_id = 0; indiv_id < (int) exp_m_->nb_indivs(); indiv_id++) {
      if (exp_m_->dna_mutator_array_[indiv_id]->hasMutate()) {
        /*printf("Translate protein for indiv %d : %ld (%d) -- %ld (%d)\n",
               indiv_id,internal_simd_struct[indiv_id]->proteins.size(),internal_simd_struct[indiv_id]->protein_count_,
               internal_simd_struct[indiv_id]->rnas.size(),internal_simd_struct[indiv_id]->rna_count_);*/

//#pragma omp parallel for firstprivate(indiv_id) schedule(dynamic)
        for (int protein_idx = 0; protein_idx <
                                  (int) internal_simd_struct[indiv_id]->protein_count_; protein_idx++) {
//#pragma omp task firstprivate(indiv_id, protein_idx) depend(inout: internal_simd_struct[indiv_id])
          {
            if (internal_simd_struct[indiv_id]->proteins[protein_idx]->is_init_) {
              int x = indiv_id / exp_m_->world()->height();
              int y = indiv_id % exp_m_->world()->height();

              int c_pos = internal_simd_struct[indiv_id]->proteins[protein_idx]->protein_start, t_pos;
              int end_pos = internal_simd_struct[indiv_id]->proteins[protein_idx]->protein_end;
              if (internal_simd_struct[indiv_id]->proteins[protein_idx]->leading_lagging ==
                  0) {
                c_pos += 13;
                end_pos -= 3;

                c_pos =
                    c_pos >= dna_size[indiv_id] ? c_pos - dna_size[indiv_id]
                                                : c_pos;
                end_pos = end_pos < 0 ? dna_size[indiv_id] + end_pos : end_pos;
              } else {
                c_pos -= 13;
                end_pos += 3;

                end_pos =
                    end_pos >= dna_size[indiv_id] ? end_pos - dna_size[indiv_id]
                                                  : end_pos;
                c_pos = c_pos < 0 ? dna_size[indiv_id] + c_pos : c_pos;
              }

              int8_t value = 0;
              int8_t codon_list[64] = {};
              int8_t codon_idx = 0;
              int32_t count_loop = 0;

              if (internal_simd_struct[indiv_id]->proteins[protein_idx]->leading_lagging ==
                  0) {
                // LEADING

                while (count_loop <
                       internal_simd_struct[indiv_id]->proteins[protein_idx]->protein_length /
                       3 &&
                       codon_idx < 64) {
#ifdef WITH_BITSET
                  codon_list[codon_idx] = internal_simd_struct[indiv_id]->dna_->bitset_->extract_codon(
                      true, c_pos);
#else
                  value = 0;
                  for (int8_t i = 0; i < 3; i++) {
                    t_pos =
                        c_pos + i >= dna_size[indiv_id] ? c_pos + i -
                                                          dna_size[indiv_id]
                                                        : c_pos + i;
                    if (internal_simd_struct[indiv_id]->dna_->data_[t_pos] ==
                        '1')
                      value += 1 << (CODON_SIZE - i - 1);
                  }
                  codon_list[codon_idx] = value;
#endif
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
                       internal_simd_struct[indiv_id]->proteins[protein_idx]->protein_length /
                       3 &&
                       codon_idx < 64) {
#ifdef WITH_BITSET
                  codon_list[codon_idx] = internal_simd_struct[indiv_id]->dna_->bitset_->extract_codon(
                      false, c_pos);
#else
                  value = 0;
                  for (int8_t i = 0; i < 3; i++) {
                    t_pos =
                        c_pos - i < 0 ? dna_size[indiv_id] + (c_pos - i) :
                        c_pos -
                        i;
                    if (internal_simd_struct[indiv_id]->dna_->data_[t_pos] !=
                        '1')
                      value += 1 << (CODON_SIZE - i - 1);
                  }
                  codon_list[codon_idx] = value;
#endif
                  codon_idx++;

                  count_loop++;

                  c_pos -= 3;
                  c_pos = c_pos < 0 ? c_pos + dna_size[indiv_id] : c_pos;
                }
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


              for (int i = 0; i < codon_idx; i++) {
                switch (codon_list[i]) {
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
              internal_simd_struct[indiv_id]->proteins[protein_idx]->m =
                  nb_m != 0 ? M / (pow(2, nb_m) - 1) : 0.5;
              internal_simd_struct[indiv_id]->proteins[protein_idx]->w =
                  nb_w != 0 ? W / (pow(2, nb_w) - 1) : 0.0;
              internal_simd_struct[indiv_id]->proteins[protein_idx]->h =
                  nb_h != 0 ? H / (pow(2, nb_h) - 1) : 0.5;

              //  ------------------------------------------------------------------------------------
              //  3) Normalize M, W and H values according to the allowed ranges (defined in macros.h)
              //  ------------------------------------------------------------------------------------
              // x_min <= M <= x_max
              // w_min <= W <= w_max
              // h_min <= H <= h_max
              internal_simd_struct[indiv_id]->proteins[protein_idx]->m =
                  (X_MAX - X_MIN) *
                  internal_simd_struct[indiv_id]->proteins[protein_idx]->m +
                  X_MIN;
              internal_simd_struct[indiv_id]->proteins[protein_idx]->w =
                  (w_max - W_MIN) *
                  internal_simd_struct[indiv_id]->proteins[protein_idx]->w +
                  W_MIN;
              internal_simd_struct[indiv_id]->proteins[protein_idx]->h =
                  (H_MAX - H_MIN) *
                  internal_simd_struct[indiv_id]->proteins[protein_idx]->h +
                  H_MIN;

              if (nb_m == 0 || nb_w == 0 || nb_h == 0 ||
                  internal_simd_struct[indiv_id]->proteins[protein_idx]->w ==
                  0.0 ||
                  internal_simd_struct[indiv_id]->proteins[protein_idx]->h ==
                  0.0) {
                internal_simd_struct[indiv_id]->proteins[protein_idx]->is_functional = false;
              } else {
                internal_simd_struct[indiv_id]->proteins[protein_idx]->is_functional = true;
              }
            }
          }
        }
      }
    }
//  }
}


    void SIMD_Individual::translate_protein(int indiv_id, double w_max) {
#ifdef WITH_FINETASKLOOP
                    #pragma omp taskloop grainsize(protein_grain_size)
#endif
        for (int protein_idx = 0; protein_idx <
                                  (int) internal_simd_struct[indiv_id]->protein_count_; protein_idx++) {
//#pragma omp task firstprivate(indiv_id, protein_idx) depend(inout: internal_simd_struct[indiv_id])
            {
                if (internal_simd_struct[indiv_id]->proteins[protein_idx]->is_init_) {
                    int x = indiv_id / exp_m_->world()->height();
                    int y = indiv_id % exp_m_->world()->height();

                    int c_pos = internal_simd_struct[indiv_id]->proteins[protein_idx]->protein_start, t_pos;
                    int end_pos = internal_simd_struct[indiv_id]->proteins[protein_idx]->protein_end;
                    if (internal_simd_struct[indiv_id]->proteins[protein_idx]->leading_lagging ==
                        0) {
                        c_pos += 13;
                        end_pos -= 3;

                        c_pos =
                                c_pos >= dna_size[indiv_id] ? c_pos - dna_size[indiv_id]
                                                            : c_pos;
                        end_pos = end_pos < 0 ? dna_size[indiv_id] + end_pos : end_pos;
                    } else {
                        c_pos -= 13;
                        end_pos += 3;

                        end_pos =
                                end_pos >= dna_size[indiv_id] ? end_pos - dna_size[indiv_id]
                                                              : end_pos;
                        c_pos = c_pos < 0 ? dna_size[indiv_id] + c_pos : c_pos;
                    }

                    int8_t value = 0;
                    int8_t codon_list[64] = {};
                    int8_t codon_idx = 0;
                    int32_t count_loop = 0;

                    if (internal_simd_struct[indiv_id]->proteins[protein_idx]->leading_lagging ==
                        0) {
                        // LEADING

                        while (count_loop <
                               internal_simd_struct[indiv_id]->proteins[protein_idx]->protein_length /
                               3 &&
                               codon_idx < 64) {
#ifdef WITH_BITSET
                            codon_list[codon_idx] = internal_simd_struct[indiv_id]->dna_->bitset_->extract_codon(
              true, c_pos);
#else
                            value = 0;
                            for (int8_t i = 0; i < 3; i++) {
                                t_pos =
                                        c_pos + i >= dna_size[indiv_id] ? c_pos + i -
                                                                          dna_size[indiv_id]
                                                                        : c_pos + i;
                                if (internal_simd_struct[indiv_id]->dna_->data_[t_pos] ==
                                    '1')
                                    value += 1 << (CODON_SIZE - i - 1);
                            }
                            codon_list[codon_idx] = value;
#endif
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
                               internal_simd_struct[indiv_id]->proteins[protein_idx]->protein_length /
                               3 &&
                               codon_idx < 64) {
#ifdef WITH_BITSET
                            codon_list[codon_idx] = internal_simd_struct[indiv_id]->dna_->bitset_->extract_codon(
              false, c_pos);
#else
                            value = 0;
                            for (int8_t i = 0; i < 3; i++) {
                                t_pos =
                                        c_pos - i < 0 ? dna_size[indiv_id] + (c_pos - i) :
                                        c_pos -
                                        i;
                                if (internal_simd_struct[indiv_id]->dna_->data_[t_pos] !=
                                    '1')
                                    value += 1 << (CODON_SIZE - i - 1);
                            }
                            codon_list[codon_idx] = value;
#endif
                            codon_idx++;

                            count_loop++;

                            c_pos -= 3;
                            c_pos = c_pos < 0 ? c_pos + dna_size[indiv_id] : c_pos;
                        }
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


                    for (int i = 0; i < codon_idx; i++) {
                        switch (codon_list[i]) {
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
                    internal_simd_struct[indiv_id]->proteins[protein_idx]->m =
                            nb_m != 0 ? M / (pow(2, nb_m) - 1) : 0.5;
                    internal_simd_struct[indiv_id]->proteins[protein_idx]->w =
                            nb_w != 0 ? W / (pow(2, nb_w) - 1) : 0.0;
                    internal_simd_struct[indiv_id]->proteins[protein_idx]->h =
                            nb_h != 0 ? H / (pow(2, nb_h) - 1) : 0.5;

                    //  ------------------------------------------------------------------------------------
                    //  3) Normalize M, W and H values according to the allowed ranges (defined in macros.h)
                    //  ------------------------------------------------------------------------------------
                    // x_min <= M <= x_max
                    // w_min <= W <= w_max
                    // h_min <= H <= h_max
                    internal_simd_struct[indiv_id]->proteins[protein_idx]->m =
                            (X_MAX - X_MIN) *
                            internal_simd_struct[indiv_id]->proteins[protein_idx]->m +
                            X_MIN;
                    internal_simd_struct[indiv_id]->proteins[protein_idx]->w =
                            (w_max - W_MIN) *
                            internal_simd_struct[indiv_id]->proteins[protein_idx]->w +
                            W_MIN;
                    internal_simd_struct[indiv_id]->proteins[protein_idx]->h =
                            (H_MAX - H_MIN) *
                            internal_simd_struct[indiv_id]->proteins[protein_idx]->h +
                            H_MIN;

                    if (nb_m == 0 || nb_w == 0 || nb_h == 0 ||
                        internal_simd_struct[indiv_id]->proteins[protein_idx]->w ==
                        0.0 ||
                        internal_simd_struct[indiv_id]->proteins[protein_idx]->h ==
                        0.0) {
                        internal_simd_struct[indiv_id]->proteins[protein_idx]->is_functional = false;
                    } else {
                        internal_simd_struct[indiv_id]->proteins[protein_idx]->is_functional = true;
                    }
                }
            }
        }


        std::map<int32_t,pProtein*> lookup;

        for (int protein_idx = 0; protein_idx <
                                  (int) internal_simd_struct[indiv_id]->protein_count_; protein_idx++) {
//#pragma omp task firstprivate(indiv_id, protein_idx) depend(inout: internal_simd_struct[indiv_id])
            {
                if (internal_simd_struct[indiv_id]->proteins[protein_idx]->is_init_) {
                    if (lookup.find(internal_simd_struct[indiv_id]->proteins[protein_idx]->protein_start) == lookup.end()) {
                        lookup[internal_simd_struct[indiv_id]->proteins[protein_idx]->protein_start] = internal_simd_struct[indiv_id]->proteins[protein_idx];
                    } else {
                        lookup[internal_simd_struct[indiv_id]->proteins[protein_idx]->protein_start]->e+=internal_simd_struct[indiv_id]->proteins[protein_idx]->e;
                        internal_simd_struct[indiv_id]->proteins[protein_idx]->is_init_ = false;
                    }
                }
            }
        }
    }


void SIMD_Individual::compute_phenotype() {

  int nb_indiv = exp_m_->nb_indivs();
//#pragma omp parallel for collapse(2)
  for (int indiv_id = 0; indiv_id < (int) nb_indiv; indiv_id++) {
    if (exp_m_->dna_mutator_array_[indiv_id]->hasMutate()) {
      for (int fuzzy_idx = 0; fuzzy_idx < 300; fuzzy_idx++) {

        internal_simd_struct[indiv_id]->phenotype[fuzzy_idx] = 0;
      }
    }
  }


/*#pragma omp parallel
#pragma omp single
  {*/
    for (int indiv_id = 0; indiv_id < (int) exp_m_->nb_indivs(); indiv_id++) {
      if (exp_m_->dna_mutator_array_[indiv_id]->hasMutate()) {
        //printf("%d -- Protein to phenotype for %ld\n",i,internal_simd_struct[i]->proteins.size());
//#pragma omp parallel for firstprivate(indiv_id)
        for (int protein_idx = 0; protein_idx <
                                  (int) internal_simd_struct[indiv_id]->protein_count_; protein_idx++) {

            if (internal_simd_struct[indiv_id]->proteins[protein_idx]->is_init_) {
              /*if (indiv_id == 908 && AeTime::time() == 87) {
                printf("Computing phenotype for 908 at 64\n");
              }*/

              if (fabs(
                  internal_simd_struct[indiv_id]->proteins[protein_idx]->w) <
                  1e-15 ||
                  fabs(
                      internal_simd_struct[indiv_id]->proteins[protein_idx]->h) <
                  1e-15)
                continue;

//#pragma omp task firstprivate(indiv_id, protein_idx) depend(inout: internal_simd_struct[indiv_id])
              {
              if (internal_simd_struct[indiv_id]->proteins[protein_idx]->is_functional) {

                // Compute triangle points' coordinates
                double x0 =
                    internal_simd_struct[indiv_id]->proteins[protein_idx]->m -
                    internal_simd_struct[indiv_id]->proteins[protein_idx]->w;
                double x1 = internal_simd_struct[indiv_id]->proteins[protein_idx]->m;
                double x2 =
                    internal_simd_struct[indiv_id]->proteins[protein_idx]->m +
                    internal_simd_struct[indiv_id]->proteins[protein_idx]->w;

                int ix0 = (int) (x0 * 300);
                int ix1 = (int) (x1 * 300);
                int ix2 = (int) (x2 * 300);

                if (ix0 < 0) ix0 = 0; else if (ix0 > (299)) ix0 = 299;
                if (ix1 < 0) ix1 = 0; else if (ix1 > (299)) ix1 = 299;
                if (ix2 < 0) ix2 = 0; else if (ix2 > (299)) ix2 = 299;

/*                if (indiv_id == 894) {
                  printf("Height %f %f %f %f -- %f\n",internal_simd_struct[indiv_id]->proteins[protein_idx]->m,
                         internal_simd_struct[indiv_id]->proteins[protein_idx]->w,
                         internal_simd_struct[indiv_id]->proteins[protein_idx]->h,
                         internal_simd_struct[indiv_id]->proteins[protein_idx]->e,
                         (internal_simd_struct[indiv_id]->proteins[protein_idx]->h *
                         internal_simd_struct[indiv_id]->proteins[protein_idx]->e)
                  );

                  printf("X0-1-2 %d %d %d (%f %f %f)\n", ix0, ix1, ix2, x0, x1,
                         x2);
                }*/

                /*if (indiv_id == 908 && AeTime::time() == 87) {
                  printf("Prot %d (%f %f %f) PHEN %d %d %d\n",protein_idx,internal_simd_struct[indiv_id]->proteins[protein_idx].m,
                  internal_simd_struct[indiv_id]->proteins[protein_idx].w,internal_simd_struct[indiv_id]->proteins[protein_idx].h*
                                                                          internal_simd_struct[indiv_id]->proteins[protein_idx].e,
                  ix0,ix1,ix2);
                }*/

                // Compute the first equation of the triangle
                double incY =
                    (internal_simd_struct[indiv_id]->proteins[protein_idx]->h *
                     internal_simd_struct[indiv_id]->proteins[protein_idx]->e) /
                    (ix1 - ix0);
                int count = 1;
                // Updating value between x0 and x1

                /*if (indiv_id == 894)
                  printf("incY first loop %f\n",incY);*/
                for (int i = ix0 + 1; i < ix1; i++) {
//#pragma omp critical
                  {
                    internal_simd_struct[indiv_id]->phenotype[i] =
                        internal_simd_struct[indiv_id]->phenotype[i] +
                        (incY * (count++));
                  }
                }

//#pragma omp critical
                {
                  internal_simd_struct[indiv_id]->phenotype[ix1] =
                      internal_simd_struct[indiv_id]->phenotype[ix1] +
                      (internal_simd_struct[indiv_id]->proteins[protein_idx]->h *
                       internal_simd_struct[indiv_id]->proteins[protein_idx]->e);
                }

                // Compute the second equation of the triangle
                incY =
                    (internal_simd_struct[indiv_id]->proteins[protein_idx]->h *
                     internal_simd_struct[indiv_id]->proteins[protein_idx]->e) /
                    (ix2 - ix1);
                count = 1;

/*                if (indiv_id == 894)
                  printf("incY second loop %f\n",incY);*/
                // Updating value between x1 and x2
                for (int i = ix1 + 1; i < ix2; i++) {
//#pragma omp atomic
                  {
                    internal_simd_struct[indiv_id]->phenotype[i] =
                        internal_simd_struct[indiv_id]->phenotype[i] +
                        ((internal_simd_struct[indiv_id]->proteins[protein_idx]->h *
                          internal_simd_struct[indiv_id]->proteins[protein_idx]->e) -
                         (incY * (count++)));
                  }
                }
              }
            }

/*              if (indiv_id == 894)
                for (int i = 213; i <= 227; i++) {
                  printf("SIMD -- X[%d] = %f\n",i,internal_simd_struct[indiv_id]->phenotype[i]);
                }*/
          }
        }
      }
    }
//  }
}


    void SIMD_Individual::compute_phenotype(int indiv_id) {

//#pragma omp task firstprivate(indiv_id) depend(inout: internal_simd_struct[indiv_id])

        double activ_phenotype[300];
        double inhib_phenotype[300];
        {
            for (int fuzzy_idx = 0; fuzzy_idx < 300; fuzzy_idx++) {
                activ_phenotype[fuzzy_idx] = 0;
                inhib_phenotype[fuzzy_idx] = 0;
            }
        }

        for (int protein_idx = 0; protein_idx <
                                  (int) internal_simd_struct[indiv_id]->protein_count_; protein_idx++) {

//#pragma omp task firstprivate(indiv_id, protein_idx) depend(inout: internal_simd_struct[indiv_id])
            {
                if (internal_simd_struct[indiv_id]->proteins[protein_idx]->is_init_) {
                    /*if (indiv_id == 908 && AeTime::time() == 87) {
                      printf("Computing phenotype for 908 at 64\n");
                    }*/

                    if (fabs(
                            internal_simd_struct[indiv_id]->proteins[protein_idx]->w) >=
                        1e-15 &&
                        fabs(
                                internal_simd_struct[indiv_id]->proteins[protein_idx]->h) >=
                        1e-15) {

                        if (internal_simd_struct[indiv_id]->proteins[protein_idx]->is_functional) {
                           /* if (indiv_id == 268) printf("Compute next prot\n");
                            if (indiv_id == 268)
                                for (int i = 200; i <= 216; i++) {
                                    printf("SIMD -- X[%d] = %f (%e %e %e)\n",i,internal_simd_struct[indiv_id]->phenotype[i],internal_simd_struct[indiv_id]->proteins[protein_idx]->m,
                                           internal_simd_struct[indiv_id]->proteins[protein_idx]->w,internal_simd_struct[indiv_id]->proteins[protein_idx]->h);
                                }*/

                            // Compute triangle points' coordinates
                            double x0 =
                                    internal_simd_struct[indiv_id]->proteins[protein_idx]->m -
                                    internal_simd_struct[indiv_id]->proteins[protein_idx]->w;
                            double x1 = internal_simd_struct[indiv_id]->proteins[protein_idx]->m;
                            double x2 =
                                    internal_simd_struct[indiv_id]->proteins[protein_idx]->m +
                                    internal_simd_struct[indiv_id]->proteins[protein_idx]->w;

                            int ix0 = (int) (x0 * 300);
                            int ix1 = (int) (x1 * 300);
                            int ix2 = (int) (x2 * 300);

                            if (ix0 < 0) ix0 = 0; else if (ix0 > (299)) ix0 = 299;
                            if (ix1 < 0) ix1 = 0; else if (ix1 > (299)) ix1 = 299;
                            if (ix2 < 0) ix2 = 0; else if (ix2 > (299)) ix2 = 299;

/*                if (indiv_id == 894) {
                  printf("Height %f %f %f %f -- %f\n",internal_simd_struct[indiv_id]->proteins[protein_idx]->m,
                         internal_simd_struct[indiv_id]->proteins[protein_idx]->w,
                         internal_simd_struct[indiv_id]->proteins[protein_idx]->h,
                         internal_simd_struct[indiv_id]->proteins[protein_idx]->e,
                         (internal_simd_struct[indiv_id]->proteins[protein_idx]->h *
                         internal_simd_struct[indiv_id]->proteins[protein_idx]->e)
                  );

                  printf("X0-1-2 %d %d %d (%f %f %f)\n", ix0, ix1, ix2, x0, x1,
                         x2);
                }*/

                            /*if (indiv_id == 908 && AeTime::time() == 87) {
                              printf("Prot %d (%f %f %f) PHEN %d %d %d\n",protein_idx,internal_simd_struct[indiv_id]->proteins[protein_idx].m,
                              internal_simd_struct[indiv_id]->proteins[protein_idx].w,internal_simd_struct[indiv_id]->proteins[protein_idx].h*
                                                                                      internal_simd_struct[indiv_id]->proteins[protein_idx].e,
                              ix0,ix1,ix2);
                            }*/

                            // Compute the first equation of the triangle
                            double incY =
                                    (internal_simd_struct[indiv_id]->proteins[protein_idx]->h *
                                     internal_simd_struct[indiv_id]->proteins[protein_idx]->e) /
                                    (ix1 - ix0);
                            int count = 1;
                            // Updating value between x0 and x1

                            /*if (indiv_id == 894)
                              printf("incY first loop %f\n",incY);*/
                            for (int i = ix0 + 1; i < ix1; i++) {
//#pragma omp critical
                                {
                                    /*if (indiv_id == 268) printf("add to %d : %f (%f %f)\n",i,(incY * (count)));*/
                                    if (internal_simd_struct[indiv_id]->proteins[protein_idx]->h > 0)
                                        activ_phenotype[i] += (incY * (count++));
                                    else
                                        inhib_phenotype[i] += (incY * (count++));
                                }
                            }

//#pragma omp critical
                            {
                               /* if (indiv_id == 268) printf("add to %d : %f\n",ix1,(internal_simd_struct[indiv_id]->proteins[protein_idx]->h *
                                                                                    internal_simd_struct[indiv_id]->proteins[protein_idx]->e));*/

                                if (internal_simd_struct[indiv_id]->proteins[protein_idx]->h > 0)
                                    activ_phenotype[ix1] +=
                                        (internal_simd_struct[indiv_id]->proteins[protein_idx]->h *
                                         internal_simd_struct[indiv_id]->proteins[protein_idx]->e);
                                else
                                    inhib_phenotype[ix1] +=
                                            (internal_simd_struct[indiv_id]->proteins[protein_idx]->h *
                                             internal_simd_struct[indiv_id]->proteins[protein_idx]->e);
                            }

                            // Compute the second equation of the triangle
                            incY =
                                    (internal_simd_struct[indiv_id]->proteins[protein_idx]->h *
                                     internal_simd_struct[indiv_id]->proteins[protein_idx]->e) /
                                    (ix2 - ix1);
                            count = 1;

/*                if (indiv_id == 894)
                  printf("incY second loop %f\n",incY);*/
                            // Updating value between x1 and x2
                            for (int i = ix1 + 1; i < ix2; i++) {
//#pragma omp atomic
                                {
                                    /*if (indiv_id == 268) printf("add to %d : %f\n",i,((internal_simd_struct[indiv_id]->proteins[protein_idx]->h *
                                                                                       internal_simd_struct[indiv_id]->proteins[protein_idx]->e) -
                                                                                      (incY * (count))));*/
                                    if (internal_simd_struct[indiv_id]->proteins[protein_idx]->h > 0)
                                        activ_phenotype[i] +=
                                            ((internal_simd_struct[indiv_id]->proteins[protein_idx]->h *
                                              internal_simd_struct[indiv_id]->proteins[protein_idx]->e) -
                                             (incY * (count++)));
                                    else
                                        inhib_phenotype[i] +=
                                                ((internal_simd_struct[indiv_id]->proteins[protein_idx]->h *
                                                  internal_simd_struct[indiv_id]->proteins[protein_idx]->e) -
                                                 (incY * (count++)));
                                }
                            }
                        }
                    }
/*
              if (indiv_id == 101)
                  for (int i = 0; i <= 1; i++) {
                  printf("SIMD -- X[%d] = %f (%e %e %e)\n",i,internal_simd_struct[indiv_id]->phenotype[i],internal_simd_struct[indiv_id]->proteins[protein_idx]->m,
                         internal_simd_struct[indiv_id]->proteins[protein_idx]->w,internal_simd_struct[indiv_id]->proteins[protein_idx]->h);
                }*/
                }
            }
        }

        for (int fuzzy_idx = 0; fuzzy_idx < 300; fuzzy_idx++) {

            if (activ_phenotype[fuzzy_idx] > 1)
                activ_phenotype[fuzzy_idx] = 1;
            if (inhib_phenotype[fuzzy_idx] < -1)
                inhib_phenotype[fuzzy_idx] = -1;

        }

        for (int fuzzy_idx = 0; fuzzy_idx < 300; fuzzy_idx++) {
            internal_simd_struct[indiv_id]->phenotype[fuzzy_idx] = activ_phenotype[fuzzy_idx] + inhib_phenotype[fuzzy_idx];
            if (internal_simd_struct[indiv_id]->phenotype[fuzzy_idx] < 0)
                internal_simd_struct[indiv_id]->phenotype[fuzzy_idx] = 0;
        }
    }

void SIMD_Individual::build_phenotypic_target(PhenotypicTargetHandler* phenotypic_target_handler) {
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
}

void SIMD_Individual::compute_fitness(double selection_pressure) {
//#pragma omp parallel for collapse(2)
  for (int indiv_id = 0; indiv_id < (int) exp_m_->nb_indivs(); indiv_id++) {
    if (exp_m_->dna_mutator_array_[indiv_id]->hasMutate()) {
      for (int fuzzy_idx = 0; fuzzy_idx < 300; fuzzy_idx++) {

        if (internal_simd_struct[indiv_id]->phenotype[fuzzy_idx] > 1)
          internal_simd_struct[indiv_id]->phenotype[fuzzy_idx] = 1;
        if (internal_simd_struct[indiv_id]->phenotype[fuzzy_idx] < 0)
          internal_simd_struct[indiv_id]->phenotype[fuzzy_idx] = 0;

        internal_simd_struct[indiv_id]->delta[fuzzy_idx] =
            internal_simd_struct[indiv_id]->phenotype[fuzzy_idx] -
            target[fuzzy_idx];
      }

      internal_simd_struct[indiv_id]->metaerror = 0;

      for (int fuzzy_idx = 0; fuzzy_idx < 299; fuzzy_idx++) {
        internal_simd_struct[indiv_id]->metaerror +=
            ((std::fabs(internal_simd_struct[indiv_id]->delta[fuzzy_idx]) +
              std::fabs(internal_simd_struct[indiv_id]->delta[fuzzy_idx + 1])) /
             (600.0));
      }

      internal_simd_struct[indiv_id]->fitness = exp(
          -selection_pressure *
          ((double) internal_simd_struct[indiv_id]->metaerror));
    }
  }
}


    void SIMD_Individual::compute_fitness(int indiv_id, double selection_pressure) {
        for (int fuzzy_idx = 0; fuzzy_idx < 300; fuzzy_idx++) {

            if (internal_simd_struct[indiv_id]->phenotype[fuzzy_idx] > 1)
                internal_simd_struct[indiv_id]->phenotype[fuzzy_idx] = 1;
            if (internal_simd_struct[indiv_id]->phenotype[fuzzy_idx] < 0)
                internal_simd_struct[indiv_id]->phenotype[fuzzy_idx] = 0;

            internal_simd_struct[indiv_id]->delta[fuzzy_idx] =
                    internal_simd_struct[indiv_id]->phenotype[fuzzy_idx] -
                    target[fuzzy_idx];
        }

        internal_simd_struct[indiv_id]->metaerror = 0;

        for (int fuzzy_idx = 0; fuzzy_idx < 299; fuzzy_idx++) {
            internal_simd_struct[indiv_id]->metaerror +=
                    ((std::fabs(internal_simd_struct[indiv_id]->delta[fuzzy_idx]) +
                      std::fabs(internal_simd_struct[indiv_id]->delta[fuzzy_idx + 1])) /
                     (600.0));
        }

        internal_simd_struct[indiv_id]->fitness = exp(
                -selection_pressure *
                ((double) internal_simd_struct[indiv_id]->metaerror));
    }


void SIMD_Individual::run_a_step(double w_max, double selection_pressure,bool optim_prom) {
#pragma omp parallel
#pragma omp single nowait
    {
        nb_clones_ = 0;
#ifndef WITH_L1TASK
        #pragma omp taskloop
#endif
        for (int indiv_id = 0; indiv_id < exp_m_->nb_indivs(); indiv_id++) {
#ifdef WITH_L1TASK
            #pragma omp task firstprivate(indiv_id)
#endif
            {

                //printf("Manage %d\n",indiv_id);

            //if (AeTime::time() > 0 && optim_prom) check_selection(indiv_id);
            if (standalone_ && optim_prom) {

                selection(indiv_id);
            }

//#pragma omp taskwait


  if (standalone_ && optim_prom) {

//    selection(indiv_id);

/*    for (int indiv_id = 0; indiv_id < exp_m_->nb_indivs(); indiv_id++) {
      int x = indiv_id / exp_m_->world()->height();
      int y = indiv_id % exp_m_->world()->height();
      delete exp_m_->dna_mutator_array_[indiv_id];

      exp_m_->dna_mutator_array_[indiv_id] = new DnaMutator(
          exp_m_->world()->grid(x, y)->mut_prng(),
          internal_simd_struct[indiv_id]->dna_->length(),
          exp_m_->exp_s()->mut_params()->duplication_rate(),
          exp_m_->exp_s()->mut_params()->deletion_rate(),
          exp_m_->exp_s()->mut_params()->translocation_rate(),
          exp_m_->exp_s()->mut_params()->inversion_rate(),
          exp_m_->exp_s()->mut_params()->point_mutation_rate(),
          exp_m_->exp_s()->mut_params()->small_insertion_rate(),
          exp_m_->exp_s()->mut_params()->small_deletion_rate(),
          exp_m_->exp_s()->mut_params()->max_indel_size(),
          exp_m_->exp_s()->min_genome_length(),
          exp_m_->exp_s()->max_genome_length());
      //exp_m_->dna_mutator_array_[indiv_id]->setMutate(true);
    }*/

    do_mutation(indiv_id);
  } else if (standalone_ && (!optim_prom)) {

      //for (int indiv_id = 0; indiv_id < exp_m_->nb_indivs(); indiv_id++) {
          int x = indiv_id / exp_m_->world()->height();
          int y = indiv_id % exp_m_->world()->height();
          delete exp_m_->dna_mutator_array_[indiv_id];

          exp_m_->dna_mutator_array_[indiv_id] = new DnaMutator(
                  exp_m_->world()->grid(x, y)->mut_prng(),
                  internal_simd_struct[indiv_id]->dna_->length(),
                  exp_m_->exp_s()->mut_params()->duplication_rate(),
                  exp_m_->exp_s()->mut_params()->deletion_rate(),
                  exp_m_->exp_s()->mut_params()->translocation_rate(),
                  exp_m_->exp_s()->mut_params()->inversion_rate(),
                  exp_m_->exp_s()->mut_params()->point_mutation_rate(),
                  exp_m_->exp_s()->mut_params()->small_insertion_rate(),
                  exp_m_->exp_s()->mut_params()->small_deletion_rate(),
                  exp_m_->exp_s()->mut_params()->max_indel_size(),
                  exp_m_->exp_s()->min_genome_length(),
                  exp_m_->exp_s()->max_genome_length());
          exp_m_->dna_mutator_array_[indiv_id]->setMutate(true);
      //}
  } else if (!standalone_ && optim_prom) {
      do_mutation(indiv_id);
  }


/*#pragma omp task firstprivate(indiv_id)
            {*/
  if (optim_prom) {
    //printf("Clear structures\n");
    //clear_struct_before_next_step();

    //printf("Apply mutation + Optimized search promoters\n");
    /*printf("Check DNA:  Optimized search promoters\n");
    */
    //check_dna();
    //printf("Optimized search stop RNA and Compute RNA\n");
    opt_prom_compute_RNA(indiv_id);
  } else {
    //if (standalone_) {
      /*double dupl_rate = exp_m_->exp_s()->mut_params()->duplication_rate();
      double del_rate = exp_m_->exp_s()->mut_params()->deletion_rate();
      double trans_rate = exp_m_->exp_s()->mut_params()->translocation_rate();
      double inv_rate = exp_m_->exp_s()->mut_params()->inversion_rate();
      double point_mut_rate = exp_m_->exp_s()->mut_params()->point_mutation_rate();
      double small_ins_rate = exp_m_->exp_s()->mut_params()->small_insertion_rate();
      double small_del_rate = exp_m_->exp_s()->mut_params()->small_deletion_rate();
      double max_indel_size = exp_m_->exp_s()->mut_params()->max_indel_size();
      double min_genome = exp_m_->exp_s()->min_genome_length();
      double max_genome = exp_m_->exp_s()->max_genome_length();


      for (int indiv_id = 0; indiv_id < exp_m_->nb_indivs(); indiv_id++) {
        int x = indiv_id / exp_m_->world()->height();
        int y = indiv_id % exp_m_->world()->height();
        delete exp_m_->dna_mutator_array_[indiv_id];

        exp_m_->dna_mutator_array_[indiv_id] = new DnaMutator(
            exp_m_->world()->grid(x, y)->mut_prng(),
            internal_simd_struct[indiv_id]->dna_->length(),
            dupl_rate,
            del_rate,
            trans_rate,
            inv_rate,
            point_mut_rate,
            small_ins_rate,
            small_del_rate,
            max_indel_size,
            min_genome,
            max_genome);
        exp_m_->dna_mutator_array_[indiv_id]->setMutate(true);
      }*/


/*      for (int indiv_id = 0; indiv_id < (int) exp_m_->nb_indivs(); indiv_id++) {
        int x = indiv_id / exp_m_->world()->height();
        int y = indiv_id % exp_m_->world()->height();

        delete exp_m_->world()->grid(x,y)->individual();
      }
    } else {

    }*/
      //for (int indiv_id = 0; indiv_id < exp_m_->nb_indivs(); indiv_id++) {
          int x = indiv_id / exp_m_->world()->height();
          int y = indiv_id % exp_m_->world()->height();
          delete exp_m_->dna_mutator_array_[indiv_id];

          exp_m_->dna_mutator_array_[indiv_id] = new DnaMutator(
                  exp_m_->world()->grid(x, y)->mut_prng(),
                  internal_simd_struct[indiv_id]->dna_->length(),
                  exp_m_->exp_s()->mut_params()->duplication_rate(),
                  exp_m_->exp_s()->mut_params()->deletion_rate(),
                  exp_m_->exp_s()->mut_params()->translocation_rate(),
                  exp_m_->exp_s()->mut_params()->inversion_rate(),
                  exp_m_->exp_s()->mut_params()->point_mutation_rate(),
                  exp_m_->exp_s()->mut_params()->small_insertion_rate(),
                  exp_m_->exp_s()->mut_params()->small_deletion_rate(),
                  exp_m_->exp_s()->mut_params()->max_indel_size(),
                  exp_m_->exp_s()->min_genome_length(),
                  exp_m_->exp_s()->max_genome_length());
          exp_m_->dna_mutator_array_[indiv_id]->setMutate(true);
      //}

    start_stop_RNA(indiv_id);

    compute_RNA(indiv_id);

      //if (indiv_id == 381) printf("Compute RNAs %d %d\n",internal_simd_struct[381]->rnas.size(),internal_simd_struct[381]->promoters.size());
  }

  /*
  for (int indiv_id = 0; indiv_id < (int) exp_m_->nb_indivs(); indiv_id++) {
      if (internal_simd_struct[indiv_id]->promoters.size() != internal_simd_struct[indiv_id]->leading_prom_pos.size() + internal_simd_struct[indiv_id]->lagging_prom_pos.size()) {
          printf("Error cache is not synchronized !!! -- %d\n",indiv_id);
      }

  }
*/
/*#pragma omp parallel
#pragma omp single
    {
        //printf("Search Protein start motifs\n");
        start_protein();
//#pragma omp taskwait
        //printf("Compute Proteins\n");
        compute_protein();
//#pragma omp taskwait
        //printf("Translate protein\n");
        translate_protein(w_max);
//#pragma omp taskwait
        //printf("Compute phenotype\n");
        compute_phenotype();
//#pragma omp taskwait
        //printf("Compute fitness\n");
        compute_fitness(selection_pressure);
    }*/

                if (exp_m_->dna_mutator_array_[indiv_id]->hasMutate()) {
//#pragma omp taskgroup
                    {
                        start_protein(indiv_id);
//#pragma omp taskwait
                        compute_protein(indiv_id);
//#pragma omp taskwait
                        translate_protein(indiv_id, w_max);
//#pragma omp taskwait
                        compute_phenotype(indiv_id);
//#pragma omp taskwait
                        compute_fitness(indiv_id, selection_pressure);
                    }
                }
                //if (indiv_id == 381) printf("Compute IndiS %d %d\n",internal_simd_struct[381]->rnas.size(),internal_simd_struct[381]->promoters.size());

                if (standalone_ && optim_prom) {
#pragma omp critical
                    {
                        int x = indiv_id / exp_m_->world()->height();
                        int y = indiv_id % exp_m_->world()->height();

                        EndReplicationEvent *eindiv = new EndReplicationEvent(
                                internal_simd_struct[indiv_id], x, y);
                        // Tell observers the replication is finished
                        internal_simd_struct[indiv_id]->notifyObservers(END_REPLICATION, eindiv);
                        delete eindiv;
                    }
                }
                //printf("Manage END %d\n",indiv_id);
            }
        }

//#pragma omp taskwait
    }


    if (optim_prom)
        notifyObservers(END_GENERATION);

    //printf("Compute BCLEAN %d %d\n",prev_internal_simd_struct[381]->rnas.size(),prev_internal_simd_struct[381]->promoters.size());
  //printf("Check results\n");
  //check_result();
  //exit(-44);
  if (optim_prom) {
    //printf("OPT -- Copy to old generation struct\n");
//#pragma omp parallel
//#pragma omp single
//#pragma omp taskloop
      for (int indiv_id = 0; indiv_id < (int) exp_m_->nb_indivs(); indiv_id++) {
          //if (indiv_id == 0)
          //    printf("%d -- usage %d -- \n", indiv_id, internal_simd_struct[indiv_id]->usage_count_);
          int usage_cpt = 0;

#pragma omp critical
          {

                  if (prev_internal_simd_struct[indiv_id]->usage_count_ == 1) {
                      prev_internal_simd_struct[indiv_id]->usage_count_ = -1;

                      delete prev_internal_simd_struct[indiv_id];
                  } else
                      prev_internal_simd_struct[indiv_id]->usage_count_--;
                  /*else {
                      printf("Still alive %d : %d\n",prev_internal_simd_struct[indiv_id]->global_id,
                             prev_internal_simd_struct[indiv_id]->usage_count_);
                  }*/


          }

          prev_internal_simd_struct[indiv_id] = internal_simd_struct[indiv_id];
          internal_simd_struct[indiv_id] = nullptr;
      }
  } else if (standalone_ && (!optim_prom)) {

      } else {
      //printf("Copy to old generation struct\n");
//#pragma omp parallel
//#pragma omp single
//#pragma omp taskloop
      for (int indiv_id = 0; indiv_id < (int) exp_m_->nb_indivs(); indiv_id++) {
          int usage_cpt = 0;

          #pragma omp critical
          {
              if (prev_internal_simd_struct[indiv_id] != internal_simd_struct[indiv_id]) {
                  if (prev_internal_simd_struct[indiv_id]->usage_count_ == 1) {
                      prev_internal_simd_struct[indiv_id]->usage_count_ = -1;

                      delete prev_internal_simd_struct[indiv_id];
                  } else
                      prev_internal_simd_struct[indiv_id]->usage_count_--;
                  //usage_cpt = prev_internal_simd_struct[indiv_id]->usage_count_;
              }

          }

          prev_internal_simd_struct[indiv_id] = internal_simd_struct[indiv_id];
          prev_internal_simd_struct[indiv_id]->clearAllObserver();
          internal_simd_struct[indiv_id] = nullptr;
    }
  }

  //printf("Compute XXX %d %d\n",prev_internal_simd_struct[381]->rnas.size(),prev_internal_simd_struct[381]->promoters.size());


    // Search for the best
  double best_fitness = prev_internal_simd_struct[0]->fitness;
  int idx_best = 0;
  for (int indiv_id = 1; indiv_id < (int) exp_m_->nb_indivs(); indiv_id++) {
    if (prev_internal_simd_struct[indiv_id]->fitness > best_fitness) {
      idx_best = indiv_id;
      best_fitness = prev_internal_simd_struct[indiv_id]->fitness;
    }
  }
/*
  printf("IDX BEST %d\n",idx_best);*/

  best_indiv = prev_internal_simd_struct[idx_best];


  // Stats
  if (!optim_prom) {
    stats_best = new Stats_SIMD(this, AeTime::time(), true);
    stats_mean = new Stats_SIMD(this, AeTime::time(), false);
  } else {
    stats_best->reinit(AeTime::time());
    stats_mean->reinit(AeTime::time());
  }

#pragma omp taskloop
    for (int indiv_id = 0; indiv_id < (int) exp_m_->nb_indivs(); indiv_id++) {
        for (int i = 0; i < prev_internal_simd_struct[indiv_id]->rna_count_; i++) {
            if (prev_internal_simd_struct[indiv_id]->rnas[i] != nullptr) {
                if (prev_internal_simd_struct[indiv_id]->rnas[i]->is_coding_)
                    prev_internal_simd_struct[indiv_id]->nb_coding_RNAs++;
                else
                    prev_internal_simd_struct[indiv_id]->nb_non_coding_RNAs++;
            }
        }

        for (int i = 0; i < prev_internal_simd_struct[indiv_id]->protein_count_; i++) {
            if (prev_internal_simd_struct[indiv_id]->rnas[i] != nullptr) {
                if (prev_internal_simd_struct[indiv_id]->proteins[i]->is_functional) {
                    prev_internal_simd_struct[indiv_id]->nb_func_genes++;
                } else {
                    prev_internal_simd_struct[indiv_id]->nb_non_func_genes++;
                }
                if (prev_internal_simd_struct[indiv_id]->proteins[i]->h > 0) {
                    prev_internal_simd_struct[indiv_id]->nb_genes_activ++;
                } else {
                    prev_internal_simd_struct[indiv_id]->nb_genes_inhib++;
                }
            }
        }
    }

  stats_best->write_best();
  stats_mean->write_average();

  if (standalone_ && AeTime::time() % exp_m_->output_m()->tree_step() == 0 && AeTime::time() > 0) {
      exp_m_->output_m()->write_tree();
  }


  if (standalone_ && AeTime::time() % exp_m_->backup_step() == 0) {

    for (int indiv_id = 0; indiv_id < (int) exp_m_->nb_indivs(); indiv_id++) {
      int x = indiv_id / exp_m_->world()->height();
      int y = indiv_id % exp_m_->world()->height();

        exp_m_->world()->grid(x, y)->individual()->clear_everything_except_dna_and_promoters();
        exp_m_->world()->grid(x, y)->individual()->genetic_unit_list_nonconst().clear();
        delete exp_m_->world()->grid(x, y)->individual();

      Individual * indiv = new Individual(exp_m_,
                                          exp_m_->world()->grid(x,y)->mut_prng(),
                                          exp_m_->world()->grid(x,y)->stoch_prng(),
                                          exp_m_->exp_s()->mut_params(),
                                          w_max,
                                          exp_m_->exp_s()->min_genome_length(),
                                          exp_m_->exp_s()->max_genome_length(),
                                          false,
                                          indiv_id,
                                          "",
                                          0);
      int32_t nb_blocks_ = prev_internal_simd_struct[indiv_id]->dna_->length()/BLOCK_SIZE + 1;
      char* dna_string = new char[nb_blocks_ * BLOCK_SIZE];
      memset(dna_string,0,
             (prev_internal_simd_struct[indiv_id]->dna_->length()+1) * sizeof(char));


      char* to_copy = prev_internal_simd_struct[indiv_id]->dna_->to_char();


      //printf("Copy DNA for indiv %d size %d (%d x %d)\n",indiv_id,prev_internal_simd_struct[indiv_id]->dna_->length(),nb_blocks_,BLOCK_SIZE);
      memcpy(dna_string, to_copy,
             (prev_internal_simd_struct[indiv_id]->dna_->length()+1) * sizeof(char));


      indiv->genetic_unit_list_.clear();
      indiv->add_GU(dna_string, prev_internal_simd_struct[indiv_id]->dna_->length());
      indiv->genetic_unit_nonconst(0).set_min_gu_length(exp_m_->exp_s()->min_genome_length());
      indiv->genetic_unit_nonconst(0).set_max_gu_length(exp_m_->exp_s()->max_genome_length());
      indiv->EvaluateInContext(exp_m_->world()->grid(x,y)->habitat());
      indiv->compute_statistical_data();

      exp_m_->world()->grid(x,y)->set_individual(indiv);

#ifdef WITH_BITSET
      delete [] to_copy;
#endif

    }

    // Create missing directories
    exp_m_->WriteDynamicFiles();

    std::ofstream last_gener_file(LAST_GENER_FNAME,
                                  std::ofstream::out);

    last_gener_file << AeTime::time() << std::endl;
    last_gener_file.close();

    if (AeTime::time() == exp_m_->end_step()) {
      for (int indiv_id = 0; indiv_id < (int) exp_m_->nb_indivs(); indiv_id++) {
        int x = indiv_id / exp_m_->world()->height();
        int y = indiv_id % exp_m_->world()->height();

          exp_m_->world()->grid(x, y)->individual()->clear_everything_except_dna_and_promoters();
          exp_m_->world()->grid(x, y)->individual()->genetic_unit_list_nonconst().clear();
        delete exp_m_->world()->grid(x, y)->individual();
      }
    }


  }

  //printf("Start to next gen\n");
//  for (int32_t index = 0; index < exp_m_->world()->width() * exp_m_->world()->height(); index++) {
//    int32_t x = index / exp_m_->world()->height();
//    int32_t y = index % exp_m_->world()->height();
//    auto indiv = exp_m_->world()->grid(x,y)->individual();
//
//    if (indiv->genetic_unit(0).dna()->length() == 8099)
//      printf("%d (%d %d : %d) (%d %d) ",indiv->genetic_unit(0).dna()->length(),
//          indiv->grid_cell()->x(),indiv->grid_cell()->y(),
//           indiv->grid_cell()->x()*exp_m_->world()->
//               height()+indiv->grid_cell()->y(),x,y);
//  }
//  printf("\n");

}

void SIMD_Individual::check_dna() {
#ifndef WITH_BITSET
  int x, y;
  for (int i = 0; i < (int) exp_m_->nb_indivs(); i++) {
    x = i / exp_m_->world()->height();
    y = i % exp_m_->world()->height();

    for (int dna_pos = 0; dna_pos < dna_size[i]; dna_pos++) {
      if (exp_m_->world()->grid(x, y)->individual()->genetic_unit(
          0).dna()->data()[dna_pos] != internal_simd_struct[i]->dna_->data_[dna_pos]) {

        printf("Check DNA indiv %d %ld %ld --- NB Mutation %ld\n",i,dna_size[i],exp_m_->world()->grid(x, y)->individual()->genetic_unit(
            0).dna()->length(),exp_m_->dna_mutator_array_[i]->mutation_list_.size());

        printf("Divergence between classic DNA and SIMD DNA %d %d at pos %d\n",
               exp_m_->world()->grid(x, y)->individual()->genetic_unit(
                   0).dna()->data()[dna_pos],
               internal_simd_struct[i]->dna_->data_[dna_pos],
               dna_pos);
        for (auto mute : exp_m_->dna_mutator_array_[i]->mutation_list_) {
          printf("Mutation type %d\n",mute->type());
        }

        break;
      }
    }
  }
#else
  printf("NOT YET IMPLEMENTED WITH BITSET !\n");
  exit(-1);
#endif
}

void SIMD_Individual::check_individual(int i, int x, int y) {
    exp_m_->world()->grid(x, y)->set_individual(exp_m_->world()->grid(x, y)->old_one);
    exp_m_->world()->grid(x, y)->old_one->Reevaluate();

    printf("%d %d %d -- ",i,x,y);

    printf(
            "Nb RNA SIMD/CPU %ld/%ld Protein %ld/%ld Metaerror %f/%f Fitness %e/%e DNA Size %d/%d\n",
            prev_internal_simd_struct[i]->rnas.size(),
            exp_m_->world()->grid(x, y)->individual()->rna_list().size(),
            prev_internal_simd_struct[i]->proteins.size(),
            exp_m_->world()->grid(x, y)->individual()->protein_list().size(),
            prev_internal_simd_struct[i]->metaerror,
            exp_m_->world()->grid(x, y)->individual()->dist_to_target_by_feature(
                    METABOLISM), prev_internal_simd_struct[i]->fitness,
            exp_m_->world()->grid(x, y)->individual()->fitness(), prev_internal_simd_struct[i]->dna_->length(),
            exp_m_->world()->grid(x, y)->individual()->genetic_unit(
                    0).seq_length());

    int idx = 0;

    for (auto rna : exp_m_->world()->grid(x, y)->old_one->rna_list()) {
        printf("RNA CPU %d Start %d Stop %d Leading/Lagging %d Length %d\n", idx,
               rna->promoter_pos(), rna->last_transcribed_pos(), rna->strand(), rna->transcript_length());
        idx++;
    }
    idx = 0;
    for (idx = 0; idx < (int) (prev_internal_simd_struct[i]->promoters.size()); idx++) {
        if (prev_internal_simd_struct[i]->promoters[idx] != nullptr)
            printf("Promoters found at %d\n",prev_internal_simd_struct[i]->promoters[idx]->pos);
    }

    idx = 0;
    for (idx = 0; idx < (int) (prev_internal_simd_struct[i]->rnas.size()); idx++) {
        if (prev_internal_simd_struct[i]->rnas[idx] != nullptr) {
            printf("RNA SIMD %d Start %d Stop %d Leading/Lagging %d Length %d\n", idx,
                   prev_internal_simd_struct[i]->rnas[idx]->begin,
                   prev_internal_simd_struct[i]->rnas[idx]->end,
                   prev_internal_simd_struct[i]->rnas[idx]->leading_lagging,
                   prev_internal_simd_struct[i]->rnas[idx]->length);
        }
    }

    int prot_cpt_a=0,prot_cpt_b=0;
    idx = 0;
    for (auto prot : exp_m_->world()->grid(x, y)->old_one->protein_list()) {
        bool found = false;

        for (int pidx = 0; pidx <
                           (int) prev_internal_simd_struct[i]->protein_count_; pidx++) {
            if (prev_internal_simd_struct[i]->proteins[pidx]->is_init_) {
                if ((prev_internal_simd_struct[i]->proteins[pidx]->e == prot->concentration()) &&
                    (prev_internal_simd_struct[i]->proteins[pidx]->protein_end == prot->last_STOP_base_pos())) {
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
                      (int) prev_internal_simd_struct[i]->protein_count_; idx++) {
        if (prev_internal_simd_struct[i]->proteins[idx]->is_init_) {


            bool found = false;

            for (auto prot : exp_m_->world()->grid(x, y)->old_one->protein_list()) {
                if (( prev_internal_simd_struct[i]->proteins[idx]->e ==  prot->concentration()) && ( prev_internal_simd_struct[i]->proteins[idx]->protein_end ==  prot->last_STOP_base_pos())) {
                    found = true;
                    break;
                }
            }

            //for (idx = 0; idx < (int) (internal_simd_struct[i]->proteins.size()); idx++) {
            if (!found)
                printf("Proteins SIMD %d Start %d (end %d) Length %d Leading/Lagging %d M/W/H %f/%f/%f Func %d -- Concentration %f\n", idx,
                       prev_internal_simd_struct[i]->proteins[idx]->protein_start,
                       prev_internal_simd_struct[i]->proteins[idx]->protein_end,
                       prev_internal_simd_struct[i]->proteins[idx]->protein_length,
                       prev_internal_simd_struct[i]->proteins[idx]->leading_lagging,
                       prev_internal_simd_struct[i]->proteins[idx]->m,
                       prev_internal_simd_struct[i]->proteins[idx]->w,
                       prev_internal_simd_struct[i]->proteins[idx]->h,
                       prev_internal_simd_struct[i]->proteins[idx]->is_functional,
                       prev_internal_simd_struct[i]->proteins[idx]->e
                );
            prot_cpt_b++;
        }
    }

}

void SIMD_Individual::check_result() {
  int x, y;

  for (int i = 0; i < (int) exp_m_->nb_indivs(); i++) {
    //if (i != 905) continue;

    x = i / exp_m_->world()->height();
    y = i % exp_m_->world()->height();


    double fit_1 = exp_m_->world()->grid(x,
                                         y)->individual()->dist_to_target_by_feature(
        METABOLISM);
    double fit_2 = internal_simd_struct[i]->metaerror;
    float i_fit_1 = roundf(fit_1 * 100);
    float i_fit_2 = roundf(fit_2 * 100);


    int prot_size = 0;


    for (auto prot : exp_m_->world()->grid(x, y)->individual()->protein_list()) {

        prot_size++;
      //}
    }

    if (i_fit_1 != i_fit_2 && dna_size[i] > 300)
    //if (i == 268) {
      if ((internal_simd_struct[i]->rnas.size() != exp_m_->world()->grid(x, y)->individual()->rna_list().size()) ||
        (internal_simd_struct[i]->protein_count_ != prot_size)) {
      printf(
          "ERROR -- Individual %d : Metaerror (CPU/GPU) : %e/%e || Fitness (CPU/GPU) : %e/%e \n",
          i,
          exp_m_->world()->grid(x, y)->individual()->dist_to_target_by_feature(
              METABOLISM),
          internal_simd_struct[i]->metaerror,
          exp_m_->world()->grid(x, y)->individual()->fitness(),
          internal_simd_struct[i]->fitness);

      printf(
          "Nb RNA SIMD/CPU %ld/%ld Protein %ld/%ld Metaerror %f/%f Fitness %e/%e DNA Size %d/%d\n",
          internal_simd_struct[i]->rnas.size(),
          exp_m_->world()->grid(x, y)->individual()->rna_list().size(),
          internal_simd_struct[i]->proteins.size(),
          exp_m_->world()->grid(x, y)->individual()->protein_list().size(),
          internal_simd_struct[i]->metaerror,
          exp_m_->world()->grid(x, y)->individual()->dist_to_target_by_feature(
              METABOLISM), internal_simd_struct[i]->fitness,
          exp_m_->world()->grid(x, y)->individual()->fitness(), dna_size[i],
          exp_m_->world()->grid(x, y)->individual()->genetic_unit(
              0).seq_length());

      int idx = 0;

      /*for (auto rna : exp_m_->world()->grid(x, y)->individual()->rna_list()) {
        printf("RNA CPU %d Start %d Stop %d Leading/Lagging %d Length %d\n", idx,
               rna->promoter_pos(), rna->last_transcribed_pos(), rna->strand(), rna->transcript_length());
        idx++;
      }

      idx = 0;
      for (idx = 0; idx < (int) (internal_simd_struct[i]->rnas.size()); idx++) {
        printf("RNA SIMD %d Start %d Stop %d Leading/Lagging %d Length %d\n", idx,
               internal_simd_struct[i]->rnas[idx]->begin,
               internal_simd_struct[i]->rnas[idx]->end,
               internal_simd_struct[i]->rnas[idx]->leading_lagging,
               internal_simd_struct[i]->rnas[idx]->length);
      }

*/

      idx = 0;
        int prot_cpt_a=0,prot_cpt_b=0;

      for (auto prot : exp_m_->world()->grid(x, y)->individual()->protein_list()) {
          bool found = false;

          for (int pidx = 0; pidx <
                            (int) internal_simd_struct[i]->protein_count_; pidx++) {
              if (internal_simd_struct[i]->proteins[pidx]->is_init_) {
                  if ((internal_simd_struct[i]->proteins[pidx]->e == prot->concentration()) &&
                      (internal_simd_struct[i]->proteins[pidx]->protein_end == prot->last_STOP_base_pos())) {
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
                                  (int) internal_simd_struct[i]->protein_count_; idx++) {
          if (internal_simd_struct[i]->proteins[idx]->is_init_) {


          bool found = false;

          for (auto prot : exp_m_->world()->grid(x, y)->individual()->protein_list()) {
              if (( internal_simd_struct[i]->proteins[idx]->e ==  prot->concentration()) && ( internal_simd_struct[i]->proteins[idx]->protein_end ==  prot->last_STOP_base_pos())) {
                  found = true;
                  break;
              }
          }

      //for (idx = 0; idx < (int) (internal_simd_struct[i]->proteins.size()); idx++) {
        if (!found)
          printf("Proteins SIMD %d Start %d (end %d) Length %d Leading/Lagging %d M/W/H %f/%f/%f Func %d -- Concentration %f\n", idx,
               internal_simd_struct[i]->proteins[idx]->protein_start,
               internal_simd_struct[i]->proteins[idx]->protein_end,
               internal_simd_struct[i]->proteins[idx]->protein_length,
               internal_simd_struct[i]->proteins[idx]->leading_lagging,
               internal_simd_struct[i]->proteins[idx]->m,
               internal_simd_struct[i]->proteins[idx]->w,
               internal_simd_struct[i]->proteins[idx]->h,internal_simd_struct[i]->proteins[idx]->is_functional,
               internal_simd_struct[i]->proteins[idx]->e
                );
        prot_cpt_b++;
          }
      }

      if (i==101) {
          compute_phenotype(i);

          for (auto &gen_unit: exp_m_->world()->grid(x, y)->individual()->genetic_unit_list_) {
              gen_unit.phenotypic_contributions_computed_ = false;
              gen_unit.activ_contribution()->clear();
              gen_unit.inhib_contribution()->clear();


              printf("Compute CPU phen (%p)\n", &gen_unit);
              gen_unit.compute_phenotypic_contribution(i);
          }

          for (int j = 0; j < 300; j++) {
              printf("%d : %f -- %f\n", j, internal_simd_struct[i]->phenotype[j],
                     ((HybridFuzzy *) exp_m_->world()->grid(x, y)->individual()->phenotype())->points()[j]);
          }
      }

/*
        idx = 0;
        for (idx = 0; idx < (int) (internal_simd_struct[i]->rnas.size()); idx++) {
          for (int idxb = 0; idxb <
                             (int) (internal_simd_struct[i]->rnas[idx].start_prot.size()); idxb++) {
            printf("Protein %d Start %d\n", idx,
                   internal_simd_struct[i]->rnas[idx].start_prot[idxb]);
          }
        }*/

  //    for (int j = 0; j < 300; j++) {
  //      printf("PHENOTYPE [%d] : %f/%f\n",j,internal_simd_struct[i]->phenotype[j],((HybridFuzzy*) exp_m_->world()->indiv_at(x, y)->phenotype())->points()[j]);
  //    }

/*      char c = getchar();
      if (c=='q')
        exit(-1);
      c = getchar();*/

    }
/*
    if (i == 0) {
      for (int j = 0; j < 300; j++) {
        printf("%d : %f/%f DELTA %f FABS %f\n",j,internal_simd_struct[i]->phenotype[j],target[j],internal_simd_struct[i]->delta[j],
               std::fabs(internal_simd_struct[i]->delta[j]));
      }
      printf("Nb RNA SIMD/CPU %ld/%ld Protein %ld/%ld Metaerror %f/%f Fitness %e/%e\n",internal_simd_struct[i]->rnas.size(),
             exp_m_->indiv_by_id(i)->rna_list().size(),internal_simd_struct[i]->proteins.size(),
             exp_m_->indiv_by_id(i)->protein_list().size(),internal_simd_struct[i]->metaerror,
             exp_m_->indiv_by_id(i)->dist_to_target_by_feature(METABOLISM),internal_simd_struct[i]->fitness,
             exp_m_->indiv_by_id(i)->fitness());

    }*/

  }
}

/** Internal_SIMD_Struct Constructor and Destructor **/
Internal_SIMD_Struct::Internal_SIMD_Struct(ExpManager* exp_m, Internal_SIMD_Struct* clone, bool copy_dna) {
  count_prom = 0;
  exp_m_ = exp_m;

  usage_count_ = 1;
  dna_ = new Dna_SIMD(clone->dna_,this,copy_dna);

  //promoters.resize(clone->promoters.size());

  for (const auto& prom : clone->promoters) {
    if (prom.second != nullptr) {
      auto prom_copy = new promoterStruct(prom.second->pos, prom.second->error,
                                          prom.second->leading_or_lagging);
      promoters[count_prom] = prom_copy;


      if (prom.second->leading_or_lagging) {
        leading_prom_pos[prom_copy->pos] = count_prom;
      } else {
        lagging_prom_pos[prom_copy->pos] = count_prom;
      }

      count_prom++;
    }
  }

  fitness = clone->fitness;
  metaerror = clone->metaerror;
  //leading_prom_pos = clone->leading_prom_pos;
  //lagging_prom_pos = clone->lagging_prom_pos;

}

Internal_SIMD_Struct::~Internal_SIMD_Struct() {
    if (promoters.size() > 0) {
        for (auto element = promoters.begin();
             element != promoters.end(); ++element) {
            if (element->second != nullptr) delete element->second;
        }
    }

  promoters.clear();

  leading_prom_pos.clear();
  lagging_prom_pos.clear();

  delete dna_;

  for (auto rn : rnas) {
    delete rn;
  }

  for (auto prot : proteins) {
    delete prot;
  }

  rnas.clear();
  proteins.clear();
}

/**
 * We need some index for the promoter optimization
 */
void Internal_SIMD_Struct::rebuild_index() {
  if (count_prom > promoters.size()/2) {
    /**
     * Do the reindexation process
     */
    auto old_promoters = promoters;
    promoters.clear();
    leading_prom_pos.clear();
    lagging_prom_pos.clear();
    count_prom = 0;
    for (auto prom : old_promoters) {
      promoters[count_prom] = prom.second;
      if (prom.second->leading_or_lagging) {
        leading_prom_pos[prom.second->pos] = count_prom;
      } else {
        lagging_prom_pos[prom.second->pos] = count_prom;
      }
      count_prom++;
    }
  }
}

/**
 * Optimized promoters search
 */


void Internal_SIMD_Struct::remove_promoters_around(int32_t pos) {
  if (dna_->length() >= PROM_SIZE) {
    remove_leading_promoters_starting_between(Utils::mod(pos - PROM_SIZE + 1,
                                                         dna_->length()),
                                              pos);
    remove_lagging_promoters_starting_between(pos,
                                              Utils::mod(pos + PROM_SIZE - 1,
                                                         dna_->length()));
  }
  else {
    remove_all_promoters();
  }
}

void Internal_SIMD_Struct::remove_promoters_around(int32_t pos_1, int32_t pos_2) {
  if (Utils::mod(pos_1 - pos_2, dna_->length()) >= PROM_SIZE) {
    remove_leading_promoters_starting_between(Utils::mod(pos_1 - PROM_SIZE + 1,
                                                         dna_->length()),
                                              pos_2);
    remove_lagging_promoters_starting_between(pos_1,
                                              Utils::mod(pos_2 + PROM_SIZE - 1,
                                                         dna_->length()));
  }
  else {
    remove_all_promoters();
  }
}


void Internal_SIMD_Struct::move_all_promoters_after(int32_t pos, int32_t delta_pos) {
  move_all_leading_promoters_after(pos, delta_pos);
  move_all_lagging_promoters_after(pos, delta_pos);
}

void Internal_SIMD_Struct::look_for_new_promoters_around(int32_t pos_1, int32_t pos_2) {
  if (dna_->length() >= PROM_SIZE) {
    look_for_new_leading_promoters_starting_between(
        Utils::mod(pos_1 - PROM_SIZE + 1,
                   dna_->length()), pos_2);

    look_for_new_lagging_promoters_starting_between(pos_1, Utils::mod(
        pos_2 + PROM_SIZE - 1,
        dna_->length()));
  }
}

void Internal_SIMD_Struct::look_for_new_promoters_around(int32_t pos) {
  if (dna_->length() >= PROM_SIZE) {
    look_for_new_leading_promoters_starting_between(
        Utils::mod(pos - PROM_SIZE + 1, dna_->length()),
        pos);
    look_for_new_lagging_promoters_starting_between(
        pos,
        Utils::mod(pos + PROM_SIZE - 1, dna_->length()));
  }
}

void Internal_SIMD_Struct::insert_promoters_at(std::vector<std::list<promoterStruct*>>& promoters_to_insert,
                                      int32_t pos) {

  for (auto strand: {LEADING, LAGGING}) {
    if (promoters_to_insert[strand].size() <= 0) {
      continue;
    }

    // Insert the promoters in the individual's RNA list
    for (auto& to_insert: promoters_to_insert[strand]) {
        int prev_pos = to_insert->pos;
      // Update promoter position
      to_insert->pos = Utils::mod(to_insert->pos + pos, dna_->length());
      if (strand == LEADING) {
          if ( leading_prom_pos.find(to_insert->pos) == leading_prom_pos.end() ) {
/*              if (indiv_id == 433) {
                printf("Before cache state %d %d %d\n",promoters.size(),leading_prom_pos.size(),lagging_prom_pos.size());
              }*/

              int prom_idx;
#pragma omp atomic capture
              {
                  prom_idx = count_prom;
                  count_prom = count_prom + 1;
              }

              promoters[prom_idx] = to_insert;


/*              if (indiv_id == 433) printf("Add new promoters %d at %d (from %d)\n", prom_idx, to_insert->pos, prev_pos);
              if (indiv_id == 433) printf("Will update leading prom %d : %d\n", to_insert->pos, prom_idx);
              if (indiv_id == 433) {
                  if (leading_prom_pos.find(to_insert->pos) == leading_prom_pos.end()) {
                      printf("Is available\n");
                  } else {
                      printf("A promoter is already at this position ?!? %d %d\n");
                  }


              }*/


              leading_prom_pos[to_insert->pos] = prom_idx;

/*              if (indiv_id == 433) {
                  printf("AFTER cache state %d %d %d\n",promoters.size(),leading_prom_pos.size(),lagging_prom_pos.size());
              }*/
          }
      } else {
              if (lagging_prom_pos.find(to_insert->pos) == lagging_prom_pos.end()) {
                  int prom_idx;
#pragma omp atomic capture
                  {
                      prom_idx = count_prom;
                      count_prom = count_prom + 1;
                  }

                  promoters[prom_idx] = to_insert;
                  lagging_prom_pos[to_insert->pos] = prom_idx;
              }
          }

    }
  }
}


void Internal_SIMD_Struct::duplicate_promoters_included_in(int32_t pos_1,
                                                  int32_t pos_2,
                                                  std::vector<std::list<promoterStruct*>>& duplicated_promoters) {
  // 1) Get promoters to be duplicated
  std::vector<std::list<promoterStruct*>> retrieved_promoters = {{},
                                           {}};

/*  if (indiv_id == 6) {
    printf("INSIDE_DUPLICATE 1 : Leading promoters lists : ");
    for (auto it : leading_prom_pos) {
      printf("%d (%d) || ", it.first, it.second);
    }
    printf("\n");

    printf("INSIDE_DUPLICATE 1 : Leading promoters lists (promoters): ");
    for (auto it : promoters) {
      printf("%d (%d) -- ", it.second->pos, it.first);
    }

    printf("\n");
  }*/

  promoters_included_in(pos_1, pos_2, retrieved_promoters);

/*  if (indiv_id == 6) {
    printf("INSIDE_DUPLICATE 2 : Leading promoters lists : ");
    for (auto it : leading_prom_pos) {
      printf("%d (%d) || ", it.first, it.second);
    }
    printf("\n");

    printf("INSIDE_DUPLICATE 2 : Leading promoters lists (promoters): ");
    for (auto it : promoters) {
      printf("%d (%d) -- ", it.second->pos, it.first);
    }

    printf("\n");
  }*/

  // 2) Set RNAs' position as their position on the duplicated segment
  for (auto& strand: {LEADING, LAGGING}) {
    for (auto& prom : retrieved_promoters[strand]) {
      // Make a copy of current RNA inside container
      duplicated_promoters[strand].push_back(new promoterStruct(prom));

      // Set RNA's position as it's position on the duplicated segment
      duplicated_promoters[strand].back()->pos = Utils::mod(duplicated_promoters[strand].back()->pos -pos_1,
                                                         dna_->length());
    }
  }

  /*if (indiv_id == 6) {
    printf("INSIDE_DUPLICATE 3 : Leading promoters lists : ");
    for (auto it : leading_prom_pos) {
      printf("%d (%d) || ", it.first, it.second);
    }
    printf("\n");

    printf("INSIDE_DUPLICATE 3 : Leading promoters lists (promoters): ");
    for (auto it : promoters) {
      printf("%d (%d) -- ", it.second->pos, it.first);
    }

    printf("\n");
  }*/
}

void Internal_SIMD_Struct::invert_promoters_included_in(int32_t pos1,
                                               int32_t pos2) {
  int32_t segment_length = pos2 - pos1;

  if (segment_length < PROM_SIZE) {
    return;
  }

  std::vector<std::list<promoterStruct*>> inverted_promoters = {{},
                                          {}};

  // 1) Extract the promoters completely included on the segment to be inverted
  extract_promoters_included_in(pos1, pos2, inverted_promoters);

  // 2) Invert segment's promoters
  Internal_SIMD_Struct::invert_promoters(inverted_promoters, pos1, pos2);

  // 3) Reinsert the inverted promoters
  insert_promoters(inverted_promoters);
}


void Internal_SIMD_Struct::extract_promoters_included_in(int32_t pos_1,
                                                int32_t pos_2,
                                                std::vector<std::list<promoterStruct*>>& extracted_promoters) {
  if (pos_2 - pos_1 < PROM_SIZE) {
    return;
  }

  extract_leading_promoters_starting_between(pos_1, pos_2 - PROM_SIZE + 1,
                                             extracted_promoters[LEADING]);
  extract_lagging_promoters_starting_between(pos_1 + PROM_SIZE - 1, pos_2,
                                             extracted_promoters[LAGGING]);
}

void Internal_SIMD_Struct::insert_promoters(std::vector<std::list<promoterStruct*>>& promoters_to_insert) {
  //promoters.resize(promoters.size()+promoters_to_insert.size()+10);

  for (auto strand: {LEADING, LAGGING}) {
    if (promoters_to_insert[strand].size() <= 0) {
      continue;
    }
    // Insert the promoters in the individual's RNA list
    for (auto& to_insert: promoters_to_insert[strand]) {
      if (strand == LEADING) {
        if (leading_prom_pos.find(to_insert->pos) == leading_prom_pos.end()) {

          int prom_idx;
#pragma omp atomic capture
          {
            prom_idx = count_prom;
            count_prom = count_prom + 1;
          }

          promoters[prom_idx] = to_insert;
          leading_prom_pos[to_insert->pos] = prom_idx;
        }
      } else {
        if (lagging_prom_pos.find(to_insert->pos) == lagging_prom_pos.end()) {
          int prom_idx;
#pragma omp atomic capture
          {
            prom_idx = count_prom;
            count_prom = count_prom + 1;
          }

          promoters[prom_idx] = to_insert;

          lagging_prom_pos[to_insert->pos] = prom_idx;
        }
      }
    }
  }
}

/*staticvoid Internal_SIMD_Struct::invert_promoters(std::vector<std::list<promoterStruct*>>& promoter_lists,
                                              int32_t pos1,
                                              int32_t pos2) {
  // Exchange LEADING and LAGGING lists
  promoter_lists[LEADING].swap(promoter_lists[LAGGING]);

  // Update the position and strand of each promoter to be inverted...
  for (auto& strand: {LEADING, LAGGING})
    for (auto& prom: promoter_lists[strand]) {
      prom->pos = pos1 + pos2 - prom->pos - 1;
      prom->leading_or_lagging = strand == LEADING ? true : false;
    }
}*/



void Internal_SIMD_Struct::shift_promoters(
    std::vector<std::list<promoterStruct*>>& promoters_to_shift,
    int32_t delta_pos, int32_t seq_length) {

  for (auto& strand: {LEADING, LAGGING})
    for (auto& prom: promoters_to_shift[strand])
      prom->pos = Utils::mod(prom->pos + delta_pos, seq_length);
}


/** Both (leading+lagging) **/
void Internal_SIMD_Struct::remove_all_promoters() {
  leading_prom_pos.clear();
  lagging_prom_pos.clear();

  //for (int prom_idx = 0; prom_idx < promoters.size(); prom_idx++) {
    //delete promoters[prom_idx];
  //}

  for (auto it = promoters.begin(),
           nextit = it;
       it != promoters.end();
       it = nextit) {
    delete it->second;
    nextit = next(it);
    promoters.erase(it);
  }

  promoters.clear();
  count_prom = 0;
}



void Internal_SIMD_Struct::invert_promoters(std::vector<std::list<promoterStruct*>>& promoter_lists,
                                   int32_t pos1,
                                   int32_t pos2) {
  // Exchange LEADING and LAGGING lists
  promoter_lists[LEADING].swap(promoter_lists[LAGGING]);

  // Update the position and strand of each promoter to be inverted...
  for (auto& strand: {LEADING, LAGGING})
    for (auto& prom: promoter_lists[strand]) {
      prom->pos = pos1 + pos2 - prom->pos - 1;
      prom->leading_or_lagging = !prom->leading_or_lagging;
    }
}

/** LEADING promoters **/
/** REMOVE **/
void Internal_SIMD_Struct::remove_leading_promoters_starting_between(int32_t pos_1,
                                                            int32_t pos_2) {
  if (pos_1 > pos_2) {
    remove_leading_promoters_starting_after(pos_1);
    remove_leading_promoters_starting_before(pos_2);
  }
  else {
    /*if (indiv_id == 128) {

      printf("Leading promoters lists : ");
      for (auto it : leading_prom_pos) {
        printf("%d (%d) || ", it.first, it.second);
      }
      printf("\n");

      printf("Leading promoters lists (promoters): ");
      for (auto it : promoters) {
        if (it.second->leading_or_lagging)
          printf("%d (%d) -- ", it.second->pos, it.first);
      }
      printf("\n");

      printf("Remove leading promoters between %d %d (starting at it %d %d)\n",
             pos_1, pos_2, leading_prom_pos.lower_bound(pos_1)->first,
             leading_prom_pos.lower_bound(pos_1)->second);
    }*/

    // STL Warning: don't erase the current iterator in the for-loop!
    for (auto it = leading_prom_pos.lower_bound(pos_1),
             nextit = it;
         it != leading_prom_pos.end() and it->first < pos_2;
         it = nextit) {

      int pidx = it->second;
      auto it_p = promoters[pidx];

/*      if (indiv_id == 6)
        printf("Deleting promoters at %d (%d) %d\n",it->first,promoters[it->second]->pos,it_p->first);*/

      delete it_p;
      promoters.erase(pidx);
      nextit = next(it);
      leading_prom_pos.erase(it);
    }

/*    if (indiv_id == 6) {
      printf("AFTER Leading promoters lists : ");
      for (auto it : leading_prom_pos) {
        printf("%d ", it.first);
      }
      printf("\n");

      printf("AFTER Leading promoters lists (promoters): ");
      for (auto it : promoters) {
        if (it.second->leading_or_lagging)
          printf("%d ", it.second->pos);
      }
      printf("\n");
    }*/
  }
}

void Internal_SIMD_Struct::remove_leading_promoters_starting_after(int32_t pos) {
  auto init_it = leading_prom_pos.lower_bound(pos);
  if (init_it == leading_prom_pos.end())
    return;



  /*if (indiv_id == 626) {
      printf("--------------------------------------------> Remove everything after %d: %d -- %d\n",pos, init_it->first,this->indiv_id);

  printf("BEFORE Leading promoters lists : ");
  for (auto it : leading_prom_pos) {
    printf("%d ", it.first);
  }
      printf("\n");
      printf("BEFORE Leading promoters lists (promoters): ");
      for (auto it : promoters) {
          if (it.second->leading_or_lagging)
              printf("%d ", it.second->pos);
      }
  printf("\n");
  }*/

  for (auto it = init_it,
           nextit = it;
       it != leading_prom_pos.end();
       it = nextit) {
   /* if (indiv_id == 626) printf("--------------------------------------------> Remove everything after %d : delete %d %d\n",pos,
           it->first,it->second,promoters[it->second]->pos);*/

    delete promoters[it->second];
    promoters.erase(it->second);
    nextit = next(it);
    leading_prom_pos.erase(it);
  }

   /* if (indiv_id == 626) {
    printf("AFTER Leading promoters lists : ");
    for (auto it : leading_prom_pos) {
        printf("%d ", it.first);
    }
        printf("\n");
    printf("AFTER Leading promoters lists (promoters): ");
    for (auto it : promoters) {
        if (it.second->leading_or_lagging)
            printf("%d ", it.second->pos);
    }
    printf("\n");
}*/
}

void Internal_SIMD_Struct::remove_leading_promoters_starting_before(int32_t pos) {

    /*if (indiv_id == 626) {
        printf("--------------------------------------------> Remove everything before %d: %d -- %d\n",pos, this->indiv_id);

        printf("BEFORE Leading promoters lists : ");
        for (auto it : leading_prom_pos) {
            printf("%d ", it.first);
        }
        printf("\n");
        printf("BEFORE Leading promoters lists (promoters): ");
        for (auto it : promoters) {
            if (it.second->leading_or_lagging)
                printf("%d ", it.second->pos);
        }
        printf("\n");
    }*/
  // Delete RNAs until we reach pos (or we reach the end of the list)
  for (auto it = leading_prom_pos.begin(),
           nextit = it;
       it != leading_prom_pos.end() and it->first < pos;
       it = nextit) {
    delete promoters[it->second];
    promoters.erase(it->second);
    nextit = next(it);
    leading_prom_pos.erase(it);
  }

    /*if (indiv_id == 626) {
        printf("AFTER Leading promoters lists : ");
        for (auto it : leading_prom_pos) {
            printf("%d ", it.first);
        }
        printf("\n");
        printf("AFTER Leading promoters lists (promoters): ");
        for (auto it : promoters) {
            if (it.second->leading_or_lagging)
                printf("%d ", it.second->pos);
        }
        printf("\n");
    }*/
}


/** MOVE **/
/// Shift (by delta_post) the positions of the promoters from the
/// LEADING strand whose starting positions are >= pos.
void Internal_SIMD_Struct::move_all_leading_promoters_after(int32_t pos,
                                                   int32_t delta_pos) {
  std::map<int32_t,int32_t> tmp_prom;

/*        if (indiv_id == 388) {
            printf("Cache state AA %d %d %d\n",promoters.size(),leading_prom_pos.size(),lagging_prom_pos.size());
        }*/
  for (auto it = leading_prom_pos.lower_bound(pos), nextit=it;
       it != leading_prom_pos.end();
       it = nextit) {

    int32_t new_pos = Utils::mod(it->first + delta_pos, dna_->length());
    int32_t prom_idx = it->second;

/*    if (indiv_id == 6) printf("Moving promoters %d : %d to %d DNA Size %d Delta %d\n",prom_idx,
                              it->first,new_pos,dna_->length(),delta_pos);*/

    promoters[it->second]->pos = new_pos;
    nextit = next(it);

      if (tmp_prom.find(new_pos) == tmp_prom.end()) {
          tmp_prom[new_pos] = prom_idx;
      } else {
          promoters.erase(it->second);
      }
      leading_prom_pos.erase(it);

  }

/*        if (indiv_id == 388) {
            printf("Cache state BB %d %d %d %d\n",promoters.size(),leading_prom_pos.size(),lagging_prom_pos.size(),tmp_prom.size());
        }*/
  for (auto to_insert : tmp_prom) {
    if (leading_prom_pos.find(to_insert.first) == leading_prom_pos.end()) {
        leading_prom_pos[to_insert.first] = to_insert.second;
    } else {
        promoters.erase(to_insert.second);
    }
  }

/*        if (indiv_id == 388) {
            printf("Cache state CC %d %d %d\n",promoters.size(),leading_prom_pos.size(),lagging_prom_pos.size());
        }*/
}

/** LOOK **/

void Internal_SIMD_Struct::locate_promoters() {
  look_for_new_leading_promoters_starting_between(0,dna_->length());
  look_for_new_lagging_promoters_starting_between(0,dna_->length());
}

void Internal_SIMD_Struct::look_for_new_leading_promoters_starting_between(int32_t pos_1,
                                                                  int32_t pos_2) {
   // When pos_1 > pos_2, we will perform the search in 2 steps.
  // As positions  0 and dna_->length() are equivalent, it's preferable to
  // keep 0 for pos_1 and dna_->length() for pos_2.

  if (pos_1 >= pos_2) {
    look_for_new_leading_promoters_starting_after(pos_1);
    look_for_new_leading_promoters_starting_before(pos_2);
    return;
  }
  // Hamming distance of the sequence from the promoter consensus

  for (int32_t i = pos_1; i < pos_2; i++) {
#ifdef WITH_BITSET
    int8_t dist = dna_->bitset_->is_promoter(true,i);
#else
    int8_t dist = is_promoter_leading(i);
#endif

    if (dist <= 4) {
      if (leading_prom_pos.find(i) == leading_prom_pos.end()) {
        promoterStruct* nprom = new promoterStruct(i, dist, true);
        {
          int prom_idx;
#pragma omp atomic capture
          {
            prom_idx = count_prom;
            count_prom = count_prom + 1;
          }

          promoters[prom_idx] = nprom;
          leading_prom_pos[i] = prom_idx;
        }
      }
    }
  }
}

void Internal_SIMD_Struct::look_for_new_leading_promoters_starting_after(int32_t pos) {
  // Hamming distance of the sequence from the promoter consensus
  //int8_t dist = 8;

  for (int32_t i = pos; i < dna_->length(); i++) {
#ifdef WITH_BITSET
    int8_t dist = dna_->bitset_->is_promoter(true,i);
#else
    int8_t dist = is_promoter_leading(i);
#endif
    if (dist <= 4) { // dist takes the hamming distance of the sequence from the consensus
      if (leading_prom_pos.find(i) == leading_prom_pos.end()) {
        promoterStruct* nprom = new promoterStruct(i, dist, true);
        {
          int prom_idx;
#pragma omp atomic capture
          {
            prom_idx = count_prom;
            count_prom = count_prom + 1;
          }

          promoters[prom_idx] = nprom;
          leading_prom_pos[i] = prom_idx;
        }
      }
    }
  }
}

void Internal_SIMD_Struct::look_for_new_leading_promoters_starting_before(int32_t pos) {
  // Hamming distance of the sequence from the promoter consensus

  for (int32_t i = 0; i < pos; i++) {
#ifdef WITH_BITSET
    int8_t dist = dna_->bitset_->is_promoter(true,i);
#else
    int8_t dist = is_promoter_leading(i);
#endif
    if (dist <= 4) { // dist takes the hamming distance of the sequence from the consensus
      if (leading_prom_pos.find(i) == leading_prom_pos.end()) {
        promoterStruct* nprom = new promoterStruct(i, dist, true);
        {
          int prom_idx;
#pragma omp atomic capture
          {
            prom_idx = count_prom;
            count_prom = count_prom + 1;
          }

          promoters[prom_idx] = nprom;
          leading_prom_pos[i] = prom_idx;
        }
      }
    }
  }
}

/** EXTRACT **/
void Internal_SIMD_Struct::extract_leading_promoters_starting_between(int32_t pos_1,
                                                             int32_t pos_2, std::list<promoterStruct*>& extracted_promoters) {
  if (pos_2 < pos_1) {

    auto first = leading_prom_pos.lower_bound(pos_1);

    if (first == leading_prom_pos.end() or first->first >= pos_2) {
      return;
    }

    // Extract the promoters (remove them from the individual's list and put them in extracted_promoters)

    for (auto it = first;
         it != leading_prom_pos.end();
         it++) {
      extracted_promoters.push_back(promoters[it->second]);
      promoters.erase(it->second);
    }

    leading_prom_pos.erase(first, leading_prom_pos.end());

    // Find the last promoters in the interval
    auto end = leading_prom_pos.lower_bound(pos_2);


    // Extract the promoters (remove them from the individual's list and put them in extracted_promoters)
    for (auto it = leading_prom_pos.begin();
         it != end;
         it++) {
      extracted_promoters.push_back(promoters[it->second]);
      promoters.erase(it->second);
    }

    leading_prom_pos.erase(leading_prom_pos.begin(),end);

  } else {

    auto first = leading_prom_pos.lower_bound(pos_1);

    if (first == leading_prom_pos.end() or first->first >= pos_2) {
      return;
    }

    // Find the last promoters in the interval
    auto end = leading_prom_pos.lower_bound(pos_2);


    // Extract the promoters (remove them from the individual's list and put them in extracted_promoters)
    for (auto it = first;
         it != end;
         it++) {
      extracted_promoters.push_back(promoters[it->second]);
      promoters.erase(it->second);
    }

    leading_prom_pos.erase(first, end);
  }
}


/** LAGGING Promoters **/
/** REMOVE **/
void Internal_SIMD_Struct::remove_lagging_promoters_starting_between(int32_t pos_1,
                                                            int32_t pos_2) {
  if (pos_1 == dna_->length()) pos_1 = 0;
  if (pos_2 == 0) pos_2 = dna_->length();

//  if (indiv_id == 77) printf("remove lagging between %d and %d\n",pos_1,pos_2);
  if (pos_1 >
      pos_2) { // vld: that's a weird case... really do this? used from remove_promoters_around()
//    if (indiv_id == 77) printf("remove lagging between POS1 > POS2 %d and %d\n",pos_1,pos_2);
    remove_lagging_promoters_starting_after(pos_1);
    remove_lagging_promoters_starting_before(pos_2);
  }

  else {
//    if (indiv_id == 77) {
//      printf("Lagging promoters lists : ");
//      for (auto it : lagging_prom_pos) {
//        printf("%d ", it.first);
//      }
//      printf("\n");
//
//      printf("Lagging promoters lists (promoters): ");
//      for (auto it : promoters) {
//        if (it.second->leading_or_lagging)
//          printf("%d ", it.first);
//      }
//      printf("\n");
//    }

    // Delete RNAs until we pass pos_1 (or we reach the end of the list)
    auto init_loop = lagging_prom_pos.lower_bound(pos_1);
    //init_loop--;

//    if (indiv_id == 77) {
//      auto ilopp = init_loop;
//      ilopp++;
//      printf("Remove Lagging promoters between %d %d (starting at it %d %d) %d\n",
//             pos_1, pos_2, init_loop->first,
//             init_loop->second,ilopp->first);
//    }

    for (auto it = init_loop,
             nextit = it;
         it != lagging_prom_pos.end() and it->first < pos_2;
         it = nextit) {

//      if (indiv_id == 77)
//        printf("----------> Deleting promoters at %d (%d)\n",it->first,promoters[it->second]->pos);
      /*if (it->first != promoters[it->second]->pos) {
        printf("errrrrrrrrrrrrrror\n");
        exit(-1);
      }*/

      delete promoters[it->second];
      promoters.erase(it->second);
      nextit = next(it);
      lagging_prom_pos.erase(it);
    }

//    if (indiv_id == 77) {
//      printf("AFTER Lagging promoters lists : ");
//
//      for (auto it : lagging_prom_pos) {
//        printf("%d ", it.first);
//      }
//      printf("\n");
//
//      printf("AFTER Lagging promoters lists (promoters): ");
//      for (auto it : promoters) {
//        if (it.second->leading_or_lagging)
//          printf("%d ", it.first);
//      }
//      printf("\n");
//    }

  }

/*  if (indiv_id == 55) {
    printf("Lagging promoters lists : ");
    for (auto it : lagging_prom_pos) {
      printf("%d ", it.first);
    }
    printf("\n");

    printf("Lagging promoters lists (promoters): ");
    for (auto it : promoters) {
      printf("READ %d ", it.first);
      if (!it.second->leading_or_lagging)
        printf("%d ", it.second->pos);
    }
    printf("\n");
  }*/
}

void Internal_SIMD_Struct::remove_lagging_promoters_starting_after(int32_t pos) {
  auto init_loop = lagging_prom_pos.lower_bound(pos);

  if (init_loop == lagging_prom_pos.end())
    return;

  /*printf("--------------------------------------------> Remove everything after %d: %d -- %d\n",pos, init_loop->first,this->indiv_id);
//  if (indiv_id == 63) {
    printf("Lagging promoters lists : ");
    for (auto it : lagging_prom_pos) {
      printf("%d ", it.first);
    }
    printf("\n");
//  }
*/
  // Delete RNAs until we pass pos (or we reach the end of the list)
  for (auto it = init_loop,
           nextit = it;
       it != lagging_prom_pos.end();
       it = nextit) {
  //    printf("--------------------------------------------> Remove everything after %d : delete %d %d\n",pos,
  //           it->first,promoters[it->second]->pos);

    delete promoters[it->second];
    promoters.erase(it->second);
    nextit = next(it);
    lagging_prom_pos.erase(it);
  }
}


void Internal_SIMD_Struct::remove_lagging_promoters_starting_before(int32_t pos) {
  // Delete RNAs until we reach pos (or we reach the end of the list)
  // TODO: optimize by starting from the end (with reverse iterators)
  auto init_loop = lagging_prom_pos.lower_bound(pos);
  if (init_loop == lagging_prom_pos.begin())
    return;

  //init_loop--;

/*  if (indiv_id == 30)
    printf("Remove everything before %d: %d (%d)\n",pos,init_loop->first,
          lagging_prom_pos.begin()->first);*/

  for (auto it = lagging_prom_pos.begin(),
           nextit = it;
       it != init_loop;
       it = nextit) {
    delete promoters[it->second];
    promoters.erase(it->second);
    nextit = next(it);
    lagging_prom_pos.erase(it);
  }
}

/** MOVE **/
void Internal_SIMD_Struct::move_all_lagging_promoters_after(int32_t pos,
                                                   int32_t delta_pos) {
  std::map<int32_t,int32_t> tmp_prom;

/*        if (indiv_id == 433) {
            printf("Cache state %d %d %d\n",promoters.size(),leading_prom_pos.size(),lagging_prom_pos.size());
        }*/
 for (auto it = lagging_prom_pos.lower_bound(pos), nextit = it;
       it != lagging_prom_pos.end() and it->first >= pos;
       it=nextit) {
   int32_t new_pos = Utils::mod(it->first + delta_pos, dna_->length());
   int32_t prom_idx = it->second;

/*   if (indiv_id == 6) printf("Moving LAGGING promoters %d : %d to %d DNA Size %d Delta %d\n",prom_idx,
                             it->first,new_pos,dna_->length(),delta_pos);*/

   promoters[it->second]->pos = new_pos;
   nextit = next(it);


     if (tmp_prom.find(new_pos) == tmp_prom.end()) {
         tmp_prom[new_pos] = prom_idx;
     } else {
         promoters.erase(it->second);
     }

     lagging_prom_pos.erase(it);
 }

// if (indiv_id == 433) {
//    printf("Cache state %d %d %d\n",promoters.size(),leading_prom_pos.size(),lagging_prom_pos.size());
// }


        for (auto to_insert : tmp_prom) {
            if (lagging_prom_pos.find(to_insert.first) == lagging_prom_pos.end()) {
                lagging_prom_pos[to_insert.first] = to_insert.second;
            }else {
                promoters.erase(to_insert.second);
            }
        }
}

/** MOVE **/

void Internal_SIMD_Struct::look_for_new_lagging_promoters_starting_between(int32_t pos_1,
                                                                  int32_t pos_2) {
  // When pos_1 > pos_2, we will perform the search in 2 steps.
  // As positions  0 and dna_->length() are equivalent, it's preferable to
  // keep 0 for pos_1 and dna_->length() for pos_2.

  if (pos_1 >= pos_2) {
    look_for_new_lagging_promoters_starting_after(pos_1);
    look_for_new_lagging_promoters_starting_before(pos_2);
    return;
  }

   // Hamming distance of the sequence from the promoter consensus
  for (int32_t i = pos_2 - 1; i >= pos_1; i--) {
#ifdef WITH_BITSET
    int8_t dist = dna_->bitset_->is_promoter(false,i);
#else
    int8_t dist = is_promoter_lagging(i);
#endif
    if (dist <= 4) { // dist takes the hamming distance of the sequence from the consensus
      if (lagging_prom_pos.find(i) == lagging_prom_pos.end()) {
        promoterStruct* nprom = new promoterStruct(i, dist, false);
        {
          int prom_idx;
#pragma omp atomic capture
          {
            prom_idx = count_prom;
            count_prom = count_prom + 1;
          }

          promoters[prom_idx] = nprom;
          lagging_prom_pos[i] = prom_idx;
        }
      }
    }
  }
}


void Internal_SIMD_Struct::look_for_new_lagging_promoters_starting_after(int32_t pos) {

  // Hamming distance of the sequence from the promoter consensus

  for (int32_t i = dna_->length() - 1; i >= pos; i--) {
#ifdef WITH_BITSET
    int8_t dist = dna_->bitset_->is_promoter(false,i);
#else
    int8_t dist = is_promoter_lagging(i);
#endif
    if (dist <= 4) { // dist takes the hamming distance of the sequence from the consensus
      if (lagging_prom_pos.find(i) == lagging_prom_pos.end()) {
        promoterStruct* nprom = new promoterStruct(i, dist, false);
        {
          int prom_idx;
#pragma omp atomic capture
          {
            prom_idx = count_prom;
            count_prom = count_prom + 1;
          }

          promoters[prom_idx] = nprom;
          lagging_prom_pos[i] = prom_idx;
        }
      }
    }
  }
}



void Internal_SIMD_Struct::look_for_new_lagging_promoters_starting_before(int32_t pos) {
  // Hamming distance of the sequence from the promoter consensus

  for (int32_t i = pos - 1; i >= 0; i--) {
#ifdef WITH_BITSET
    int8_t dist = dna_->bitset_->is_promoter(false,i);
#else
    int8_t dist = is_promoter_lagging(i);
#endif
    if (dist <= 4) { // dist takes the hamming distance of the sequence from the consensus
      if (lagging_prom_pos.find(i) == lagging_prom_pos.end()) {
        promoterStruct* nprom = new promoterStruct(i, dist, false);
        {
          int prom_idx;
#pragma omp atomic capture
          {
            prom_idx = count_prom;
            count_prom = count_prom + 1;
          }

          promoters[prom_idx] = nprom;
          lagging_prom_pos[i] = prom_idx;
        }
      }
    }
  }
}

/** EXTRACT **/
void Internal_SIMD_Struct::extract_lagging_promoters_starting_between(int32_t pos_1,
                                                             int32_t pos_2,
                                                                      std::list<promoterStruct*>& extracted_promoters) {

//  printf("extract between %d and %d\n",pos_1,pos_2);
  if (pos_1 > pos_2) {
    // From pos_1 to start


    // Find the last promoters in the interval
    auto end = lagging_prom_pos.lower_bound(pos_1);

    if (end != lagging_prom_pos.begin()) end--;

    // Extract the promoters (remove them from the individual's list and put the in extracted_promoters)

    /*printf("Extract promoters in between AXC %d and %d\n",lagging_prom_pos.begin()->first,end->first);*/


    for (auto it = lagging_prom_pos.begin();
         it != end;
         it++) {
      //printf("Adding to extract %d\n",it->second);

      extracted_promoters.push_back(promoters[it->second]);
      promoters.erase(it->second);
    }

    lagging_prom_pos.erase(lagging_prom_pos.begin(), end);

    // From end to pos_2

    auto first = lagging_prom_pos.lower_bound(pos_2);

    if (first != lagging_prom_pos.begin()) first--;

    if (first == lagging_prom_pos.end() or first->first < pos_2) {
      return;
    }

//    printf("Extract promoters in between AXD %d and %d\n",first->first,lagging_prom_pos.end());

    for (auto it = first;
         it != lagging_prom_pos.end();
         it++) {
      //printf("Adding to extract %d\n",it->second);

      extracted_promoters.push_back(promoters[it->second]);
      promoters.erase(it->second);
    }

    lagging_prom_pos.erase(first,lagging_prom_pos.end());

  } else {

    auto first = lagging_prom_pos.lower_bound(pos_2);

    //if (first != lagging_prom_pos.end()) first++;

    // Find the last promoters in the interval
    auto end = lagging_prom_pos.lower_bound(pos_1);

    //if (end != lagging_prom_pos.begin()) end--;

    //printf("BBB Extract promoters in between %d and %d\n",first->first,end->first);

    if (end == lagging_prom_pos.end() or end->first < pos_1)
      return;

    // Extract the promoters (remove them from the individual's list and put the in extracted_promoters)

//    printf("Extract promoters in between %d and %d\n",first->first,end->first);

    for (auto it = end;
         it != first;
         it++) {
      //printf("Adding to extract %d\n",it->second);

      extracted_promoters.push_back(promoters[it->second]);
      //printf("Erase %d (ou of %d)\n",it->second,promoters.size());
      promoters.erase(it->second);
    }

    lagging_prom_pos.erase(end, first);
  }
}

/** Generic function of SIMD_Individual **/
#ifndef WITH_BITSET
int8_t Internal_SIMD_Struct::is_promoter_leading(int pos) {
  int8_t prom_dist_leading[26];
  int len = dna_->length();

  for (int motif_id = 0; motif_id < 22; motif_id++) {
    prom_dist_leading[motif_id] =
        PROM_SEQ_LEAD[motif_id] ==
        dna_->data_[pos + motif_id >= len ? pos + motif_id - len : pos + motif_id]
        ? 0 : 1;
  }

  return prom_dist_leading[0] +
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

}


int8_t Internal_SIMD_Struct::is_promoter_lagging(int pos) {
  int8_t prom_dist[26];
  int len = dna_->length();

  for (int motif_id = 0; motif_id < 22; motif_id++) {
    prom_dist[motif_id] =
        PROM_SEQ_LAG[motif_id] ==
        dna_->data_[pos - motif_id < 0 ? len + pos - motif_id : pos - motif_id]
        ? 0 : 1;
  }

  return prom_dist[0] +
         prom_dist[1] +
         prom_dist[2] +
         prom_dist[3] +
         prom_dist[4] +
         prom_dist[5] +
         prom_dist[6] +
         prom_dist[7] +
         prom_dist[8] +
         prom_dist[9] +
         prom_dist[10] +
         prom_dist[11] +
         prom_dist[12] +
         prom_dist[13] +
         prom_dist[14] +
         prom_dist[15] +
         prom_dist[16] +
         prom_dist[17] +
         prom_dist[18] +
         prom_dist[19] +
         prom_dist[20] +
         prom_dist[21];
}
#endif

void Internal_SIMD_Struct::promoters_included_in(int32_t pos_1,
                                        int32_t pos_2,
                                        std::vector<std::list<promoterStruct*>>& promoters_list) {
  if (pos_1 < pos_2) {
    int32_t seg_length = pos_2 - pos_1;

    if (seg_length >= PROM_SIZE) {
      lst_promoters(LEADING, BETWEEN, pos_1, pos_2 - PROM_SIZE + 1,
                promoters_list[LEADING]);
      lst_promoters(LAGGING, BETWEEN, pos_2, pos_1 + PROM_SIZE - 1,
                promoters_list[LAGGING]);
    }
  }
  else {
    int32_t seg_length = dna_->length() + pos_2 - pos_1;

    /*printf("promoters included in %d and %d\n",pos_1,pos_2);*/

    if (seg_length >= PROM_SIZE) {
      bool is_near_end_of_genome = (pos_1 + PROM_SIZE > dna_->length());
      bool is_near_beginning_of_genome = (pos_2 - PROM_SIZE < 0);

      if (!is_near_end_of_genome && !is_near_beginning_of_genome) {
        /*printf("-----------------> leading promoters after %d (till end)\n",pos_1);*/
        lst_promoters(LEADING, AFTER, pos_1, -1, promoters_list[LEADING]);
        /*printf("-----------------> leading promoters before %d (till end)\n",pos_2 - PROM_SIZE + 1);*/
        lst_promoters(LEADING, BEFORE, -1, pos_2 - PROM_SIZE + 1,
                  promoters_list[LEADING]);
        /*printf("-----------------> lagging promoters after %d (till end)\n",pos_2);*/
        lst_promoters(LAGGING, AFTER, pos_2, -1, promoters_list[LAGGING]);
        /*printf("-----------------> lagging promoters before %d (till end)\n",pos_1 + PROM_SIZE - 1);*/
        lst_promoters(LAGGING, BEFORE, -1, pos_1 + PROM_SIZE - 1,
                  promoters_list[LAGGING]);
      }
      else if (!is_near_end_of_genome) // => && is_near_beginning_of_genome
      {
        // promoters(leading, between, pos_1, pos_2 + dna_->length() - PROM_SIZE + 1,
        //                                         promoters_list[LEADING]);
        lst_promoters(LEADING, BETWEEN, pos_1, pos_2 - PROM_SIZE + 1 +
                                           dna_->length(),
                  promoters_list[LEADING]);
        lst_promoters(LAGGING, AFTER, pos_2, -1, promoters_list[LAGGING]);
        lst_promoters(LAGGING, BEFORE, -1, pos_1 + PROM_SIZE - 1,
                  promoters_list[LAGGING]);
      }
      else if (!is_near_beginning_of_genome) // => && is_near_end_of_genome
      {
        lst_promoters(LEADING, AFTER, pos_1, -1, promoters_list[LEADING]);
        lst_promoters(LEADING, BEFORE, -1, pos_2 - PROM_SIZE + 1,
                  promoters_list[LEADING]);
        lst_promoters(LAGGING, BETWEEN, pos_2, pos_1 + PROM_SIZE - 1 -
                                           dna_->length(),
                  promoters_list[LAGGING]);
      }
      else // is_near_end_of_genome && is_near_beginning_of_genome
      {
        // promoters(leading, between, pos_1, pos_2 + dna_->length() - PROM_SIZE + 1,
        //                                         promoters_list[LEADING]);
        lst_promoters(LEADING, BETWEEN, pos_1, pos_2 - PROM_SIZE + 1 +
                                           dna_->length(),
                  promoters_list[LEADING]);
        lst_promoters(LAGGING, BETWEEN, pos_2, pos_1 + PROM_SIZE - 1 -
                                           dna_->length(),
                  promoters_list[LAGGING]);
      }
    }
  }
}


void Internal_SIMD_Struct::lst_promoters(bool lorl,
                            Position before_after_btw, // with regard to the strand's reading direction
                            int32_t pos1,
                            int32_t pos2,
                            std::list<promoterStruct*>& promoters_list) {
  auto it_begin = lagging_prom_pos.begin();
  auto it_end = lagging_prom_pos.end();


  if (lorl == LEADING) {
    it_begin = leading_prom_pos.begin();
    it_end = leading_prom_pos.end();
  }

  //if (indiv_id == 40)
//  printf("Pos 1 %d Pos 2 %d LorL %d Position %d\n",pos1,pos2,lorl,before_after_btw);

  if (before_after_btw != BEFORE && pos1 != -1) {
    if (lorl == LEADING) {
      /*printf("Leading prom list : ");
      for (auto prom : leading_prom_pos) {
        printf("%d (%d) ",prom.first,prom.second);
      }
      printf("\n");*/

      auto tmp_it = leading_prom_pos.lower_bound(pos1);
      if (tmp_it == leading_prom_pos.end())
        return;

      /*printf("Search begin (AFTER LEADING) %d : %d (%d)\n",pos1,tmp_it->first,tmp_it->second);*/

      if (tmp_it!=leading_prom_pos.end()) it_begin = tmp_it;

      /*printf("UPDATED Search begin (AFTER LEADING) %d : %d\n",pos1,tmp_it->first);*/
    } else {
      //if (indiv_id == 40) {
      //printf("Lagging prom list : ");

        //for (auto prom : lagging_prom_pos) {
          // printf("%d (%d) ",prom.first,prom.second);
          //}
        //printf("\n");
        //}
      auto tmp_it = lagging_prom_pos.lower_bound(pos1);

      //if (indiv_id == 40) {
        //printf("Search begin (AFTER LAGGING) %d : %d (%d)\n", pos1,
        //       tmp_it->first, tmp_it->second);

        //if (tmp_it != lagging_prom_pos.begin()) {
        //  tmp_it--;
        //}

        //printf("UPDATED Search begin (AFTER LAGGING) %d : %d (%d)\n", pos1,
        //       tmp_it->first, tmp_it->second);
      //}
      it_end = tmp_it;
    }
  }

  if (before_after_btw != AFTER && pos2 != -1) {
    if (lorl == LEADING) {
      auto tmp_it = leading_prom_pos.lower_bound(pos2);

      /*printf("Search begin (BEFORE LEADING) %d : %d\n",pos2,tmp_it->first);*/

      if (tmp_it!=leading_prom_pos.end()) it_end = tmp_it;

      /*printf("UPDATED Search begin (BEFORE LEADING) %d : %d\n",pos2,tmp_it->first);*/
    } else {
      //if (indiv_id == 40) {printf("Lagging prom list : ");
        //for (auto prom : lagging_prom_pos) {
        //printf("%d (%d) ",prom.first,prom.second);
        //}
        //printf("\n");}

      auto tmp_it = lagging_prom_pos.lower_bound(pos2);

      //if (indiv_id == 40) {printf("Search begin (BEFORE LAGGING) %d : %d\n",pos2,tmp_it->first);}

      /*if (tmp_it!=lagging_prom_pos.begin()) {
        tmp_it--;
      }*/
      //if (indiv_id == 40) {printf("UPDATED Search begin (BEFORE LAGGING) %d : %d\n",pos2,tmp_it->first);}

      it_begin = tmp_it;
    }
  }

  //if (indiv_id == 40) {printf("Searching for promoters in %d :  %d and %d : ",
  //       lorl,it_begin->first,it_end->first);}
//  printf("Searching for promoters in %d (pos_1 %d -- pos_2 %d):  %d and %d : ",
//         lorl,pos1,pos2,
//        it_begin->first,it_end->first);

  for (auto it = it_begin; it!=it_end; it++) {
    //if (indiv_id == 40) {}
//    printf("%d (%d) ",it->first,it->second);
    promoters_list.push_back(promoters[it->second]);
  }
  //if (indiv_id == 40) {}
//  printf("\n");
}

}
