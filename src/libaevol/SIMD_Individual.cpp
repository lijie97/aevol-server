//
// Created by arrouan on 27/07/17.
//

#include "SIMD_Individual.h"
#include "Dna_SIMD.h"
#include "DnaMutator.h"
#include "HybridFuzzy.h"
#include "Stats_SIMD.h"
#include "ExpManager.h"
//#include "SIMD_Abstract_Metadata.h"
#include "SIMD_Map_Metadata.h"
#include "SIMD_DynTab_Metadata.h"
#include "SIMD_List_Metadata.h"
#include "../../../../../../../../../usr/include/stdio.h"
#include "../../../../../../../../../usr/include/stdint.h"


#include <omp.h>
#include <chrono>

namespace aevol {

#ifndef WITH_STANDALONE_SIMD
     bool SIMD_Individual::standalone_simd = false;
#else
     bool SIMD_Individual::standalone_simd = true;
#endif




SIMD_Individual::SIMD_Individual(ExpManager* exp_m) {

    printf("  Loading SIMD Controller...");

    standalone_ = standalone_simd;
    exp_m_ = exp_m;

    nb_indivs_ = exp_m_->nb_indivs();

    internal_simd_struct = new Internal_SIMD_Struct *[exp_m_->nb_indivs()];
    prev_internal_simd_struct = new Internal_SIMD_Struct *[exp_m_->nb_indivs()];

    next_generation_reproducer_ = new int32_t[exp_m_->nb_indivs()];

    // Allocate Dna_SIMD
    int max_size_dna = -1;
    for (int indiv_id = 0; indiv_id < exp_m_->nb_indivs(); indiv_id++) {
        int x = indiv_id / exp_m_->world()->height();
        int y = indiv_id % exp_m_->world()->height();

        max_size_dna = max_size_dna < exp_m->world()->grid(x, y)->individual()->genetic_unit_seq_length(0) ?
                exp_m->world()->grid(x, y)->individual()->genetic_unit_seq_length(0) : max_size_dna;
    }

    dna_factory_ = new SIMD_DnaFactory(DnaFactory_Policy::FIRSTFIT,exp_m_->nb_indivs()*3,max_size_dna);


    for (int indiv_id = 0; indiv_id < exp_m_->nb_indivs(); indiv_id++) {
        int x = indiv_id / exp_m_->world()->height();
        int y = indiv_id % exp_m_->world()->height();

        internal_simd_struct[indiv_id] = new Internal_SIMD_Struct(exp_m_, exp_m_->best_indiv()->w_max(),dna_factory_);
        internal_simd_struct[indiv_id]->dna_ = dna_factory_->get_dna(exp_m->world()->grid(x, y)->individual()->genetic_unit_seq_length(0));
        internal_simd_struct[indiv_id]->dna_->set_indiv(exp_m->world()->grid(x, y)->individual()->genetic_unit(0).dna(),dna_factory_);
        internal_simd_struct[indiv_id]->indiv_id = indiv_id;
        internal_simd_struct[indiv_id]->parent_id = indiv_id;
        prev_internal_simd_struct[indiv_id] = internal_simd_struct[indiv_id];
        internal_simd_struct[indiv_id]->global_id = AeTime::time() * 1024 + indiv_id;
    }

    dna_size = new int[exp_m_->nb_indivs()];
    for (int indiv_id = 0; indiv_id < exp_m_->nb_indivs(); indiv_id++) {
        dna_size[indiv_id] = internal_simd_struct[indiv_id]->dna_->length();
    }


    target = new double[300];
    for (int i = 0; i < 300; i++) {
        double tmp = ((HybridFuzzy *) exp_m->world()->phenotypic_target_handler()->phenotypic_target().fuzzy())->points()[i];


        target[i] = tmp;

#ifdef WITH_PERF_TRACES
        std::ofstream perf_traces_file_;
        perf_traces_file_.open("simd_perf_traces.csv", std::ofstream::trunc);
        perf_traces_file_ << "Generation,Indiv_ID,Runtime" << std::endl;
        perf_traces_file_.close();
#endif
    }

    printf(" OK\n");

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

        //int * indiv_index = new int[neighborhood_size];

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


                //indiv_index[count] = cur_x * grid_height + cur_y;
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


        /*if (indiv_id==30 || indiv_id==61) {

            printf("New indiv %d parent %d\n",indiv_id,next_generation_reproducer_[indiv_id]);
        }*/
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
                                prev_internal_simd_struct[indiv_index[i]]->metadata_->rna_count(),
                                exp_m_->world()->grid(v_x, v_y)->individual()->rna_list().size(),
                                prev_internal_simd_struct[indiv_index[i]]->metadata_->proteins_count(),
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


    void SIMD_Individual::do_mutation(int indiv_id) {
    //printf("DO MUTATION for %d -- Begin\n",indiv_id);

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


//        if (indiv_id == 17 || indiv_id == 50 || indiv_id == 51 || indiv_id == 18) {
//            int parent_id = next_generation_reproducer_[indiv_id];
//
//            printf("FROM MUTATE -- %d -- %d -- Parent %d \n",time(),indiv_id,next_generation_reproducer_[indiv_id]);
//            printf("FROM MUTATE -- %d -- %d -- BEFORE -- Prom list LEAD : ",time(),next_generation_reproducer_[indiv_id]);
//            for (int prom_idx = 0; prom_idx < prev_internal_simd_struct[parent_id]->metadata_->promoter_count(); prom_idx++) {
//                if (prev_internal_simd_struct[parent_id]->metadata_->promoters(prom_idx) != nullptr)
//                    if (prev_internal_simd_struct[parent_id]->metadata_->promoters(prom_idx)->leading_or_lagging)
//                        printf("%d ",prev_internal_simd_struct[parent_id]->metadata_->promoters(prom_idx)->pos);
//            }
//            printf("\n");
//            printf("FROM MUTATE -- %d -- %d -- BEFORE -- Prom list LAG : ",time(),next_generation_reproducer_[indiv_id]);
//            for (int prom_idx = 0; prom_idx < prev_internal_simd_struct[parent_id]->metadata_->promoter_count(); prom_idx++) {
//                if (prev_internal_simd_struct[parent_id]->metadata_->promoters(prom_idx) != nullptr)
//                    if (!prev_internal_simd_struct[parent_id]->metadata_->promoters(prom_idx)->leading_or_lagging)
//                        printf("%d ",prev_internal_simd_struct[parent_id]->metadata_->promoters(prom_idx)->pos);
//            }
//            printf("\n");
//        }

        if (standalone_) {

                int x = indiv_id / exp_m_->world()->height();
                int y = indiv_id % exp_m_->world()->height();
                delete exp_m_->dna_mutator_array_[indiv_id];

                exp_m_->dna_mutator_array_[indiv_id] = new DnaMutator(
                        exp_m_->world()->grid(x, y)->mut_prng(),
                        prev_internal_simd_struct[next_generation_reproducer_[indiv_id]]->dna_->length(),
                        //dna_size[indiv_id],
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
                //printf("%d - has mutate or not %d\n",indiv_id,exp_m_->dna_mutator_array_[indiv_id]->hasMutate());
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
                        [next_generation_reproducer_[indiv_id]],dna_factory_);

                internal_simd_struct[indiv_id]->global_id = AeTime::time()*1024+indiv_id;
                internal_simd_struct[indiv_id]->indiv_id = indiv_id;
                internal_simd_struct[indiv_id]->parent_id =
                        next_generation_reproducer_[indiv_id];

  //              printf("%d -- %d -- A -- Number of RNAs %d (%d)\n",time(),indiv_id,internal_simd_struct[indiv_id]->metadata_->rna_count(),
   //                    internal_simd_struct[indiv_id]->metadata_->promoter_count());

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
#ifdef WITH_PERF_TRACES
                auto t_start = std::chrono::steady_clock::now();
#endif
                if (standalone_)
                    internal_simd_struct[indiv_id]->dna_->apply_mutations_standalone();
                else
                    internal_simd_struct[indiv_id]->dna_->apply_mutations();
#ifdef WITH_PERF_TRACES
                auto t_end = std::chrono::steady_clock::now();
                apply_mutation[indiv_id] = t_end.time_since_epoch().count() - t_start.time_since_epoch().count();
#endif
//                printf("%d -- %d -- B -- Number of RNAs %d (%d)\n",time(),indiv_id,internal_simd_struct[indiv_id]->metadata_->rna_count(),
//                       internal_simd_struct[indiv_id]->metadata_->promoter_count());
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
#ifdef WITH_PERF_TRACES
                apply_mutation[indiv_id] = -1;
#endif
            }

//        if (indiv_id == 17 || indiv_id == 50 || indiv_id == 51 || indiv_id == 18) {
//            printf("FROM MUTATE -- %d -- %d -- Parent %d (Mutate %d)\n",time(),indiv_id,internal_simd_struct[indiv_id]->parent_id,
//                   exp_m_->dna_mutator_array_[indiv_id]->hasMutate());
//            printf("FROM MUTATE -- %d -- %d -- AFTER -- Prom list LEAD : ",time(),internal_simd_struct[indiv_id]->indiv_id);
//            for (int prom_idx = 0; prom_idx < internal_simd_struct[indiv_id]->metadata_->promoter_count(); prom_idx++) {
//                if (internal_simd_struct[indiv_id]->metadata_->promoters(prom_idx) != nullptr)
//                    if (internal_simd_struct[indiv_id]->metadata_->promoters(prom_idx)->leading_or_lagging)
//                        printf("%d ",internal_simd_struct[indiv_id]->metadata_->promoters(prom_idx)->pos);
//            }
//            printf("\n");
//            printf("FROM MUTATE -- %d -- %d -- AFTER -- Prom list LAG : ",time(),internal_simd_struct[indiv_id]->indiv_id);
//            for (int prom_idx = 0; prom_idx < internal_simd_struct[indiv_id]->metadata_->promoter_count(); prom_idx++) {
//                if (internal_simd_struct[indiv_id]->metadata_->promoters(prom_idx) != nullptr)
//                    if (!internal_simd_struct[indiv_id]->metadata_->promoters(prom_idx)->leading_or_lagging)
//                        printf("%d ",internal_simd_struct[indiv_id]->metadata_->promoters(prom_idx)->pos);
//            }
//            printf("\n");
//        }
/*
        std::set<int> leading;
            printf("FROM MUTATE -- %d -- %d -- AFTER -- Prom list LEAD : ",time(),indiv_id);
            for (int prom_idx = 0; prom_idx < internal_simd_struct[indiv_id]->metadata_->promoter_count(); prom_idx++) {
                if (internal_simd_struct[indiv_id]->metadata_->promoters(prom_idx) != nullptr)
                    if (internal_simd_struct[indiv_id]->metadata_->promoters(prom_idx)->leading_or_lagging)
                        leading.insert(internal_simd_struct[indiv_id]->metadata_->promoters(prom_idx)->pos);
            }
            for (auto lead : leading) {
                printf("%d ",lead);
            }


            printf("\n");

        std::set<int> lagging;
            printf("FROM MUTATE -- %d -- %d -- AFTER -- Prom list LAG : ",time(),indiv_id);
            for (int prom_idx = 0; prom_idx < internal_simd_struct[indiv_id]->metadata_->promoter_count(); prom_idx++) {
                if (internal_simd_struct[indiv_id]->metadata_->promoters(prom_idx) != nullptr)
                    if (!internal_simd_struct[indiv_id]->metadata_->promoters(prom_idx)->leading_or_lagging)
                        lagging.insert(internal_simd_struct[indiv_id]->metadata_->promoters(prom_idx)->pos);
            }

        for (auto lag : lagging) {
            printf("%d ",lag);
        }
            printf("\n");
*/
            //printf("%d ",internal_simd_struct[indiv_id]->metadata_->promoters(prom_idx)->pos);
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

        //printf("DO MUTATION for %d -- End\n",indiv_id);
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

                  for (int rn = 0; rn < internal_simd_struct[indiv_id]->metadata_->rna_count(); rn++) {
                      delete internal_simd_struct[indiv_id]->metadata_->rnas(rn);
                  }

                  internal_simd_struct[indiv_id]->metadata_->rnas_clear();
                  for (int rn = 0; rn < internal_simd_struct[indiv_id]->metadata_->proteins_count(); rn++) {
                      delete internal_simd_struct[indiv_id]->metadata_->proteins(rn);
                  }
                  internal_simd_struct[indiv_id]->metadata_->proteins_clear();

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

              for (int rn = 0; rn < prev_internal_simd_struct[indiv_id]->metadata_->rna_count(); rn++) {
                  delete prev_internal_simd_struct[indiv_id]->metadata_->rnas(rn);
              }

              prev_internal_simd_struct[indiv_id]->metadata_->rnas_clear();
              for (int rn = 0; rn < prev_internal_simd_struct[indiv_id]->metadata_->proteins_count(); rn++) {
                  delete prev_internal_simd_struct[indiv_id]->metadata_->proteins(rn);
              }
              prev_internal_simd_struct[indiv_id]->metadata_->proteins_clear();

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

                for (int rn = 0; rn < prev_internal_simd_struct[indiv_id]->metadata_->rna_count(); rn++) {
                    delete prev_internal_simd_struct[indiv_id]->metadata_->rnas(rn);
                }

                prev_internal_simd_struct[indiv_id]->metadata_->rnas_clear();
                for (int rn = 0; rn < prev_internal_simd_struct[indiv_id]->metadata_->proteins_count(); rn++) {
                    delete prev_internal_simd_struct[indiv_id]->metadata_->proteins(rn);
                }
                prev_internal_simd_struct[indiv_id]->metadata_->proteins_clear();

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

  delete dna_factory_;



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


    void SIMD_Individual::start_stop_RNA(int indiv_id) {
        //int nb_indiv = exp_m_->nb_indivs();
        //int x, y;

        //
        //ExpManager* exp_m = exp_m_;

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
                //int x = indiv_id / exp_m->world()->height();
                //int y = indiv_id % exp_m->world()->height();

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

                            prom_idx = internal_simd_struct[indiv_id]->metadata_->promoter_count();
                            internal_simd_struct[indiv_id]->metadata_->set_promoters_count(
                                    internal_simd_struct[indiv_id]->metadata_->promoter_count()+ 1);

              /*if (indiv_id == 30)
                printf("Adding promoters LEADING %d at %d\n",dna_pos,prom_idx);*/

                            internal_simd_struct[indiv_id]->metadata_->promoter_add(prom_idx, nprom);
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
                        internal_simd_struct[indiv_id]->metadata_->terminator_add(LEADING, dna_pos);
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
                            prom_idx = internal_simd_struct[indiv_id]->metadata_->promoter_count();
                            internal_simd_struct[indiv_id]->metadata_->set_promoters_count(
                                    internal_simd_struct[indiv_id]->metadata_->promoter_count() + 1);


                            /*if (indiv_id == 30)
                                printf("Adding promoters LAGGING %d at %d\n",dna_pos,prom_idx);*/


                            internal_simd_struct[indiv_id]->metadata_->promoter_add(prom_idx, nprom);
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
                        internal_simd_struct[indiv_id]->metadata_->terminator_add(LAGGING, dna_pos);
                    }
                }
            }

        /*if (indiv_id==30) {
            printf("%d -- %d -- TOKEEP -- Prom list LEAD : ", time(), indiv_id);
            for (int prom_idx = 0; prom_idx < internal_simd_struct[indiv_id]->metadata_->promoter_count(); prom_idx++) {
                if (internal_simd_struct[indiv_id]->metadata_->promoters(prom_idx) != nullptr)
                    if (internal_simd_struct[indiv_id]->metadata_->promoters(prom_idx)->leading_or_lagging)
                        printf("%d ", internal_simd_struct[indiv_id]->metadata_->promoters(prom_idx)->pos);
            }
            printf("\n");
            printf("%d -- %d -- TOKEEP -- Prom list LAG : ", time(), internal_simd_struct[indiv_id]->indiv_id);
            for (int prom_idx = 0; prom_idx < internal_simd_struct[indiv_id]->metadata_->promoter_count(); prom_idx++) {
                if (internal_simd_struct[indiv_id]->metadata_->promoters(prom_idx) != nullptr)
                    if (!internal_simd_struct[indiv_id]->metadata_->promoters(prom_idx)->leading_or_lagging)
                        printf("%d ", internal_simd_struct[indiv_id]->metadata_->promoters(prom_idx)->pos);
            }
            printf("\n");


            printf("%d -- %d -- DIRECT -- Prom list LEAD : ", time(), internal_simd_struct[indiv_id]->indiv_id);
            for (auto prom : ((SIMD_List_Metadata*)internal_simd_struct[indiv_id]->metadata_)->promoters_list_[LEADING]) {
                printf("%d ", prom.pos);
            }
            printf("\n");
            printf("%d -- %d -- DIRECT -- Prom list LAG : ", time(), internal_simd_struct[indiv_id]->indiv_id);
            for (auto prom : ((SIMD_List_Metadata*)internal_simd_struct[indiv_id]->metadata_)->promoters_list_[LAGGING]) {
                printf("%d ", prom.pos);
            }
            printf("\n");


            printf("%d -- %d -- ADV-DIRECT -- Prom list LEAD : ", time(), internal_simd_struct[indiv_id]->indiv_id);
            for (int prom_idx = 0; prom_idx < ((SIMD_List_Metadata*)internal_simd_struct[indiv_id]->metadata_)->promoters_list_[LEADING].size(); prom_idx++) {
                auto it =  ((SIMD_List_Metadata*)internal_simd_struct[indiv_id]->metadata_)->promoters_list_[LEADING].begin();
                std::advance(it, prom_idx);
                printf("%d ", (*it).pos);
            }
            printf("\n");
            printf("%d -- %d -- ADV-DIRECT -- Prom list LAG : ", time(), internal_simd_struct[indiv_id]->indiv_id);
            for (int prom_idx = 0; prom_idx < ((SIMD_List_Metadata*)internal_simd_struct[indiv_id]->metadata_)->promoters_list_[LAGGING].size(); prom_idx++) {
                auto it =  ((SIMD_List_Metadata*)internal_simd_struct[indiv_id]->metadata_)->promoters_list_[LAGGING].begin();
                std::advance(it, prom_idx);
                printf("%d ", (*it).pos);
            }
            printf("\n");
        }*/

//}
//}
//#pragma omp taskwait
    }


    void SIMD_Individual::opt_prom_compute_RNA(int indiv_id) {
        if (exp_m_->dna_mutator_array_[indiv_id]->hasMutate()) {
            internal_simd_struct[indiv_id]->metadata_->proteins_clear();
            internal_simd_struct[indiv_id]->metadata_->rnas_clear();
            internal_simd_struct[indiv_id]->metadata_->terminators_clear();
        }

        if (exp_m_->dna_mutator_array_[indiv_id]->hasMutate()) {

            internal_simd_struct[indiv_id]->metadata_->rnas_resize(
                    internal_simd_struct[indiv_id]->metadata_->promoter_count());

            internal_simd_struct[indiv_id]->metadata_->promoter_begin();

            for (int prom_idx = 0;
                 prom_idx < (int) internal_simd_struct[indiv_id]->metadata_->promoter_count(); prom_idx++) {
                promoterStruct *prom = internal_simd_struct[indiv_id]->metadata_->promoter_next();

                if (prom != nullptr) {
                    Dna_SIMD *dna = internal_simd_struct[indiv_id]->dna_;
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
                                prom_pos + 22;
                        cur_pos = cur_pos >= dna_length ? cur_pos -
                                                          dna_length :
                                  cur_pos;
                        int start_pos = cur_pos;

                        bool terminator_found = false;
                        bool no_terminator = false;
                        int term_dist_leading = 0;

                        int loop_size = 0;

                        while (!terminator_found) {
                            loop_size++;
                            for (int t_motif_id = 0; t_motif_id < 4; t_motif_id++)
                                term_dist_leading +=
                                        dna->get_lead(cur_pos + t_motif_id) !=
                                        dna->get_lead(cur_pos - t_motif_id + 10)
                                        ? 1 : 0;

                            if (term_dist_leading == 4)
                                terminator_found = true;
                            else {
                                cur_pos = cur_pos + 1 >= dna_length ? cur_pos + 1 -
                                                                      dna_length
                                                                    :
                                          cur_pos + 1;

                                term_dist_leading = 0;
                                if (cur_pos == start_pos) {
                                    no_terminator = true;
                                    terminator_found = true;
                                }
                            }
                        }
                        if (!no_terminator) {

                            int32_t rna_end =
                                    cur_pos + 10 >= dna_length ?
                                    cur_pos + 10 - dna_length :
                                    cur_pos + 10;

                            int32_t rna_length = 0;

                            if (prom_pos
                                > rna_end)
                                rna_length = dna_length -
                                             prom_pos
                                             + rna_end;
                            else
                                rna_length = rna_end - prom_pos;

                            rna_length -= 21;

                            if (rna_length > 0) {
                                int glob_rna_idx = -1;
                                glob_rna_idx = internal_simd_struct[indiv_id]->metadata_->rna_count_++;
                                internal_simd_struct[indiv_id]->metadata_->rna_add(glob_rna_idx, new pRNA(
                                        prom_pos,
                                        rna_end,
                                        !lead_lag,
                                        1.0 -
                                        prom_error /
                                        5.0, rna_length));
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

                            if (term_dist_lagging == 4)
                                terminator_found = true;
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


                            int32_t rna_end =
                                    cur_pos - 10 < 0 ? dna_length + (cur_pos - 10) :
                                    cur_pos -
                                    10;

                            int32_t rna_length = 0;

                            if (prom_pos <
                                rna_end)
                                rna_length =
                                        prom_pos +
                                        dna_length - rna_end;
                            else
                                rna_length =
                                        prom_pos -
                                        rna_end;

                            rna_length -= 21;

                            if (rna_length >= 0) {
                                int glob_rna_idx = -1;
                                glob_rna_idx = internal_simd_struct[indiv_id]->metadata_->rna_count_++;


                                internal_simd_struct[indiv_id]->metadata_->rna_add(glob_rna_idx, new pRNA(
                                        prom_pos,
                                        rna_end,
                                        !lead_lag,
                                        1.0 -
                                        prom_error /
                                        5.0, rna_length));
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
                internal_simd_struct[indiv_id]->metadata_->rnas_resize(
                        internal_simd_struct[indiv_id]->metadata_->promoter_count());
//#pragma omp parallel for firstprivate(indiv_id) schedule(dynamic)
#ifdef WITH_FINETASKLOOP
#pragma omp taskloop grainsize(rna_grain_size)
#endif
                for (int rna_idx = 0; rna_idx <
                                      (int) internal_simd_struct[indiv_id]->metadata_->promoter_count(); rna_idx++) {
//#pragma omp task firstprivate(indiv_id, rna_idx)
                    {
                        /*if (indiv_id == 345 && AeTime::time() == 47 &&
                            internal_simd_struct[indiv_id]->promoters[rna_idx]->pos == 4744) {
                          printf("Searching for an end with start pos %d LorL %d\n",
                                 internal_simd_struct[indiv_id]->promoters[rna_idx]->pos,
                                 internal_simd_struct[indiv_id]->promoters[rna_idx]->leading_or_lagging);
                        }*/

                        if (internal_simd_struct[indiv_id]->metadata_->promoters(rna_idx) != nullptr) {
                            if (internal_simd_struct[indiv_id]->metadata_->promoters(rna_idx)->leading_or_lagging) {
                                if (internal_simd_struct[indiv_id]->metadata_->terminator_count(LEADING) != 0) {


                                    int k =
                                            internal_simd_struct[indiv_id]->metadata_->promoters(rna_idx)->pos + 22;
                                    k = k >= dna_size[indiv_id] ? k - dna_size[indiv_id] : k;

/*        if ((indiv_id == 309 && AeTime::time() == 105) ||
            (indiv_id == 915 && AeTime::time() == 26 && rna_idx == 11)) {
          printf("Looking at %d\n",k);
        }*/

                                    int32_t next_rna_end = internal_simd_struct[indiv_id]->metadata_->next_terminator(LEADING,k);

                                    int32_t rna_end =
                                            next_rna_end + 10 >= dna_size[indiv_id] ?
                                            next_rna_end + 10 - dna_size[indiv_id] :
                                            next_rna_end + 10;

                                    /*if (indiv_id == 309 && AeTime::time() == 105) {
                                      printf("Looking for term from %d (start rna %d) : %d Computed end %d\n",k,internal_simd_struct[indiv_id]->promoters[rna_idx]->pos,
                                             *it_rna_end,rna_end);
                                    }*/
                                    int32_t rna_length = 0;

                                    if (internal_simd_struct[indiv_id]->metadata_->promoters(rna_idx)->pos
                                        > rna_end)
                                        rna_length = dna_size[indiv_id] -
                                                     internal_simd_struct[indiv_id]->metadata_->promoters(rna_idx)->pos
                                                     + rna_end;
                                    else
                                        rna_length = rna_end - internal_simd_struct[indiv_id]->
                                                metadata_->promoters(rna_idx)->pos;

                                    rna_length -= 21;

                                    if (rna_length >= 0) {


                                        int glob_rna_idx = -1;
#pragma omp critical
                                        {
                                            glob_rna_idx = internal_simd_struct[indiv_id]->metadata_->rna_count();
                                            internal_simd_struct[indiv_id]->metadata_->set_rna_count(
                                                    internal_simd_struct[indiv_id]->metadata_->rna_count() + 1);
                                        }

                                        internal_simd_struct[indiv_id]->metadata_->rna_add(glob_rna_idx, new pRNA(
                                                internal_simd_struct[indiv_id]->metadata_->promoters(rna_idx)->pos,
                                                rna_end,
                                                !internal_simd_struct[indiv_id]->metadata_->promoters(rna_idx)->leading_or_lagging,
                                                1.0 -
                                                fabs(
                                                        ((float) internal_simd_struct[indiv_id]->metadata_->promoters(rna_idx)->error)) /
                                                5.0, rna_length));
                                    }
                                }
                            } else {
                                // LAGGING
                                if (internal_simd_struct[indiv_id]->metadata_->terminator_count(LAGGING) != 0) {




                                    // Search for terminator
                                    int k =
                                            internal_simd_struct[indiv_id]->metadata_->promoters(rna_idx)->pos - 22;
                                    k = k < 0 ? dna_size[indiv_id] + k : k;

                                    int32_t next_rna_end = internal_simd_struct[indiv_id]->metadata_->next_terminator(LAGGING,k);


                                    int32_t rna_end =
                                            next_rna_end - 10 < 0 ? dna_size[indiv_id] + (next_rna_end - 10)
                                                                 :
                                            next_rna_end - 10;

                                    /*if (indiv_id == 969 && AeTime::time() == 137) {
                                      auto it_rn = it_rna_end;
                                      it_rn++;
                                      printf("Looking for term from %d (start rna %d) : %d Computed end %d (next end %d)\n",k,internal_simd_struct[indiv_id]->promoters[rna_idx]->pos,
                                             *it_rna_end,rna_end,*it_rn);
                                    }*/
                                    int32_t rna_length = 0;

                                    if (internal_simd_struct[indiv_id]->metadata_->promoters(rna_idx)->pos <
                                        rna_end)
                                        rna_length =
                                                internal_simd_struct[indiv_id]->metadata_->promoters(rna_idx)->pos +
                                                dna_size[indiv_id] - rna_end;
                                    else
                                        rna_length =
                                                internal_simd_struct[indiv_id]->metadata_->promoters(rna_idx)->pos -
                                                rna_end;

                                    rna_length -= 21;

                                    if (rna_length >= 0) {


                                        int glob_rna_idx = -1;
#pragma omp critical
                                        {
                                            glob_rna_idx = internal_simd_struct[indiv_id]->metadata_->rna_count();
                                            internal_simd_struct[indiv_id]->metadata_->set_rna_count(
                                                    internal_simd_struct[indiv_id]->metadata_->rna_count() + 1);
                                        }

                                        internal_simd_struct[indiv_id]->metadata_->rna_add(rna_idx, new pRNA(
                                                internal_simd_struct[indiv_id]->metadata_->promoters(rna_idx)->pos,
                                                rna_end,
                                                !internal_simd_struct[indiv_id]->metadata_->promoters(rna_idx)->leading_or_lagging,
                                                1.0 -
                                                fabs(
                                                        ((float) internal_simd_struct[indiv_id]->metadata_->promoters(rna_idx)->error)) /
                                                5.0, rna_length));

                                    }
                                }
                            }
                        }
                    }
                }
            }

    }

    void SIMD_Individual::start_protein(int indiv_id) {
        internal_simd_struct[indiv_id]->metadata_->rna_begin();

        for (int rna_idx = 0; rna_idx <
                              (int) internal_simd_struct[indiv_id]->metadata_->rna_count(); rna_idx++) {
            {

                pRNA* rna = internal_simd_struct[indiv_id]->metadata_->rna_next();
                int32_t dna_length = dna_size[indiv_id];

                if (rna->is_init_) {
                    int c_pos = rna->begin;

                    if (rna->length >= 22) {
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

                        while (c_pos !=
                               rna->end) {
                            bool start = false;
                            int t_pos, k_t;

                            if (rna->leading_lagging ==
                                0) {
                                // Search for Shine Dalgarro + START codon on LEADING
                                for (int k = 0; k < 9; k++) {
                                    k_t = k >= 6 ? k + 4 : k;
                                    t_pos = c_pos + k_t;

                                    if (internal_simd_struct[indiv_id]->dna_->get_lead(t_pos) ==
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
                                    t_pos = c_pos - k_t;

                                    if (internal_simd_struct[indiv_id]->dna_->get_lag(t_pos) ==
                                        SHINE_DAL_SEQ_LAG[k]) {
                                        start = true;
                                    } else {
                                        start = false;
                                        break;
                                    }
                                }
                            }

                            if (start) {
                                int start_prot_idx = -1;
                                start_prot_idx = rna->start_prot_count_++;

                                rna->start_prot[start_prot_idx] = c_pos;
                            }

                            if (rna->leading_lagging ==
                                0) {
                                c_pos++;
                                c_pos =
                                        c_pos >= dna_length ? c_pos - dna_length
                                                                    : c_pos;
                            } else {
                                c_pos--;
                                c_pos = c_pos < 0 ? dna_length + c_pos : c_pos;
                            }
                        }
                    }
                }
            }
        }
    }


void SIMD_Individual::compute_protein(int indiv_id) {
    int resize_to = 0;

    internal_simd_struct[indiv_id]->metadata_->rna_begin();
    for (int rna_idx = 0; rna_idx <
                          (int) internal_simd_struct[indiv_id]->metadata_->rna_count(); rna_idx++) {
        pRNA* rna = internal_simd_struct[indiv_id]->metadata_->rna_next();

        if (rna->is_init_)
            resize_to += rna->start_prot_count_;
    }

    internal_simd_struct[indiv_id]->
            metadata_->proteins_resize(resize_to);

    Dna_SIMD* dna = internal_simd_struct[indiv_id]->dna_;
    int32_t dna_length = dna->length();

    internal_simd_struct[indiv_id]->metadata_->rna_begin();
    for (int rna_idx = 0; rna_idx <
                          (int) internal_simd_struct[indiv_id]->metadata_->rna_count(); rna_idx++) {

        pRNA* rna = internal_simd_struct[indiv_id]->metadata_->rna_next();

        if (rna->is_init_) {
            for (int protein_idx = 0;
                 protein_idx < (int) rna->start_prot_count_; protein_idx++) {

                int start_prot = rna->start_prot[protein_idx];
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
                        length = rna->start_prot[protein_idx] - rna->end;
                    } else {
                        length = rna->start_prot[protein_idx] +
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
                                glob_prot_idx = internal_simd_struct[indiv_id]->metadata_->proteins_count();
                                internal_simd_struct[indiv_id]->metadata_->set_proteins_count(
                                        internal_simd_struct[indiv_id]->metadata_->proteins_count() +
                                        1);


                                    internal_simd_struct[indiv_id]->
                                            metadata_->protein_add(glob_prot_idx, new pProtein(
                                            start_prot, t_k,
                                            prot_length,
                                            rna->leading_lagging,
                                            rna->e
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
                                glob_prot_idx = internal_simd_struct[indiv_id]->metadata_->proteins_count();
                                internal_simd_struct[indiv_id]->metadata_->set_proteins_count(
                                        internal_simd_struct[indiv_id]->metadata_->proteins_count() +
                                        1);
                                internal_simd_struct[indiv_id]->metadata_->protein_add(
                                        glob_prot_idx, new pProtein(
                                                start_prot, t_k,
                                                prot_length,
                                                rna->leading_lagging,
                                                rna->e
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
}


    void SIMD_Individual::translate_protein(int indiv_id, double w_max) {
        internal_simd_struct[indiv_id]->metadata_->protein_begin();
        for (int protein_idx = 0; protein_idx <
                                  (int) internal_simd_struct[indiv_id]->metadata_->proteins_count(); protein_idx++) {
            pProtein* prot = internal_simd_struct[indiv_id]->metadata_->protein_next();

            if (prot->is_init_) {
                int c_pos = prot->protein_start, t_pos;
                int end_pos = prot->protein_end;
                if (prot->leading_lagging ==
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

                if (prot->leading_lagging ==
                    0) {
                    // LEADING

                    while (count_loop <
                           prot->protein_length /
                           3 &&
                           codon_idx < 64) {
                        value = 0;
                        for (int8_t i = 0; i < 3; i++) {
                            t_pos = c_pos + i;
                            if (internal_simd_struct[indiv_id]->dna_->get_lead(t_pos) ==
                                '1')
                                value += 1 << (CODON_SIZE - i - 1);
                        }
                        codon_list[codon_idx] = value;

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
                           prot->protein_length /
                           3 &&
                           codon_idx < 64) {
                        value = 0;
                        for (int8_t i = 0; i < 3; i++) {
                            t_pos = c_pos - i;
                            if (internal_simd_struct[indiv_id]->dna_->get_lag(t_pos) !=
                                '1')
                                value += 1 << (CODON_SIZE - i - 1);
                        }
                        codon_list[codon_idx] = value;
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


        std::map<int32_t, pProtein *> lookup;

        internal_simd_struct[indiv_id]->metadata_->protein_begin();
        for (int protein_idx = 0; protein_idx <
                                  (int) internal_simd_struct[indiv_id]->metadata_->proteins_count(); protein_idx++) {
            {
                pProtein* prot  = internal_simd_struct[indiv_id]->metadata_->protein_next();

                if (prot->is_init_) {
                    if (lookup.find(prot->protein_start) ==
                        lookup.end()) {
                        lookup[prot->protein_start] = prot;
                    } else {
                        lookup[prot->protein_start]->e += prot->e;
                        prot->is_init_ = false;
                    }
                }
            }
        }
    }

    void SIMD_Individual::compute_phenotype(int indiv_id) {
        double activ_phenotype[300];
        double inhib_phenotype[300];
        {
            for (int fuzzy_idx = 0; fuzzy_idx < 300; fuzzy_idx++) {
                activ_phenotype[fuzzy_idx] = 0;
                inhib_phenotype[fuzzy_idx] = 0;
            }
        }

        internal_simd_struct[indiv_id]->metadata_->protein_begin();
        for (int protein_idx = 0; protein_idx < internal_simd_struct[indiv_id]->metadata_->proteins_count(); protein_idx++) {
            pProtein* prot = internal_simd_struct[indiv_id]->metadata_->protein_next();
            if (prot->is_init_) {
                if (fabs(
                        prot->w) >=
                    1e-15 &&
                    fabs(
                            prot->h) >=
                    1e-15) {

                    if (prot->is_functional) {
                        // Compute triangle points' coordinates
                        double x0 =
                                prot->m -
                                prot->w;
                        double x1 = prot->m;
                        double x2 =
                                prot->m +
                                prot->w;

                        double height = (prot->h *
                                         prot->e);

                        int loop_A_start = (int) std::ceil(x0 * 299.0);
                        loop_A_start = loop_A_start < 0 ? 0 : loop_A_start;
                        loop_A_start = loop_A_start > 299 ? 299 : loop_A_start;

                        int loop_A_end = (int) std::ceil(x1 * 299.0);
                        loop_A_end = loop_A_end < 0 ? 0 : loop_A_end;
                        loop_A_end = loop_A_end > 299 ? 299 : loop_A_end;

                        for (int i = loop_A_start; i < loop_A_end; i++) {
                            if (prot->h > 0)
                                activ_phenotype[i] += (((i / 299.0) - x0) / (x1 - x0)) * height;
                            else
                                inhib_phenotype[i] += (((i / 299.0) - x0) / (x1 - x0)) * height;
                        }

                        // Compute the second equation of the triangle
                        // Updating value between x1 and x2
                        int loop_B_start = (int) std::ceil(x1 * 299.0);
                        loop_B_start = loop_B_start < 0 ? 0 : loop_B_start;
                        loop_B_start = loop_B_start > 299 ? 299 : loop_B_start;

                        int loop_B_end = (int) std::ceil(x2 * 299.0);
                        if (loop_B_end > 299) {
                            if (prot->h > 0)
                                activ_phenotype[299] += height * ((x2 - 1.0) / (x2 - x1));
                            else
                                inhib_phenotype[299] += height * ((x2 - 1.0) / (x2 - x1));
                        }

                        loop_B_end = loop_B_end < 0 ? 0 : loop_B_end;
                        loop_B_end = loop_B_end > 299 ? 299 : loop_B_end;

                        for (int i = loop_B_start; i < loop_B_end; i++) {
                            if (prot->h > 0)
                                activ_phenotype[i] += height * ((x2 - (i / 299.0)) / (x2 - x1));
                            else
                                inhib_phenotype[i] += height * ((x2 - (i / 299.0)) / (x2 - x1));
                        }

                    }
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
#pragma omp single
    nb_clones_ = 0;
//
#pragma omp for schedule(dynamic)
    for (int g_indiv_id = 0; g_indiv_id < exp_m_->nb_indivs(); g_indiv_id+=1) {
        {
            for (int indiv_id = g_indiv_id; indiv_id < g_indiv_id + 1; indiv_id++) {
                //printf("COMPUTE INDIV %d -- Begin\n",indiv_id);
                if (standalone_ && optim_prom) {
                    selection(indiv_id);
                }

                if (standalone_ && optim_prom) {
                    do_mutation(indiv_id);
                } else if (standalone_ && (!optim_prom)) {
                    int x = indiv_id / exp_m_->world()->height();
                    int y = indiv_id % exp_m_->world()->height();
                    delete exp_m_->dna_mutator_array_[indiv_id];

                    exp_m_->dna_mutator_array_[indiv_id] = new DnaMutator(
                            exp_m_->world()->grid(x, y)->mut_prng(),
                            prev_internal_simd_struct[next_generation_reproducer_[indiv_id]]->dna_->length(),
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
                } else if (!standalone_ && optim_prom) {
                    do_mutation(indiv_id);
                }

                if (optim_prom) {
                    opt_prom_compute_RNA(indiv_id);
                } else {
                    int x = indiv_id / exp_m_->world()->height();
                    int y = indiv_id % exp_m_->world()->height();
                    delete exp_m_->dna_mutator_array_[indiv_id];

                    exp_m_->dna_mutator_array_[indiv_id] = new DnaMutator(
                            exp_m_->world()->grid(x, y)->mut_prng(),
                            prev_internal_simd_struct[next_generation_reproducer_[indiv_id]]->dna_->length(),
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
                    compute_phenotype(indiv_id);
                    compute_fitness(indiv_id, selection_pressure);
                }

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
                //printf("COMPUTE INDIV %d -- End\n",indiv_id);
            }
        }
    }

#pragma omp single
    {
        if (optim_prom)
            notifyObservers(END_GENERATION);

        if (optim_prom) {
            for (int indiv_id = 0; indiv_id < (int) exp_m_->nb_indivs(); indiv_id++) {
                {

                    if (prev_internal_simd_struct[indiv_id]->usage_count_ == 1) {
                        prev_internal_simd_struct[indiv_id]->usage_count_ = -1;
                        //printf("DESTRUCTOR SIMD INDIV -- Main Factory %p\n",dna_factory_);
                        //printf("DESTRUCTOR SIMD INDIV -- DNA Factory %p",
                        //       prev_internal_simd_struct[indiv_id]->dna_factory_);

                        //printf(" --- DNA %p\n",
                        //       prev_internal_simd_struct[indiv_id]->dna_);

                        delete prev_internal_simd_struct[indiv_id];
                    } else
                        prev_internal_simd_struct[indiv_id]->usage_count_--;
                }

                prev_internal_simd_struct[indiv_id] = internal_simd_struct[indiv_id];
                internal_simd_struct[indiv_id] = nullptr;
            }
        } else if (standalone_ && (!optim_prom)) {

        } else {
            for (int indiv_id = 0; indiv_id < (int) exp_m_->nb_indivs(); indiv_id++) {
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

        // Search for the best
        double best_fitness = prev_internal_simd_struct[0]->fitness;
        int idx_best = 0;
        for (int indiv_id = 1; indiv_id < (int) exp_m_->nb_indivs(); indiv_id++) {
            if (prev_internal_simd_struct[indiv_id]->fitness > best_fitness) {
                idx_best = indiv_id;
                best_fitness = prev_internal_simd_struct[indiv_id]->fitness;
            }
        }

        best_indiv = prev_internal_simd_struct[idx_best];

        // Traces
#ifdef WITH_PERF_TRACES
        std::ofstream perf_traces_file_;
        perf_traces_file_.open("simd_perf_traces.csv", std::ofstream::app);
        for (int indiv_id = 0; indiv_id < (int) exp_m_->nb_indivs(); indiv_id++) {
            perf_traces_file_ << AeTime::time() << "," << indiv_id << "," << apply_mutation[indiv_id] << std::endl;
        }
        perf_traces_file_.close();
#endif

        // Stats
        if (!optim_prom) {
            stats_best = new Stats_SIMD(this, AeTime::time(), true);
            stats_mean = new Stats_SIMD(this, AeTime::time(), false);
        } else {
            stats_best->reinit(AeTime::time());
            stats_mean->reinit(AeTime::time());
        }

        for (int indiv_id = 0; indiv_id < (int) exp_m_->nb_indivs(); indiv_id++) {
            prev_internal_simd_struct[indiv_id]->reset_stats();
            prev_internal_simd_struct[indiv_id]->metadata_->rna_begin();
            for (int i = 0; i < prev_internal_simd_struct[indiv_id]->metadata_->rna_count(); i++) {
                pRNA* rna = prev_internal_simd_struct[indiv_id]->metadata_->rna_next();
                if (rna != nullptr) {
                    if (rna->is_coding_)
                        prev_internal_simd_struct[indiv_id]->nb_coding_RNAs++;
                    else
                        prev_internal_simd_struct[indiv_id]->nb_non_coding_RNAs++;
                }
            }


            for (int i = 0; i < prev_internal_simd_struct[indiv_id]->metadata_->proteins_count(); i++) {
                pProtein* prot = prev_internal_simd_struct[indiv_id]->metadata_->proteins(i);
                if (prot != nullptr) {
                    if (prot->is_functional) {
                        prev_internal_simd_struct[indiv_id]->nb_func_genes++;
                    } else {
                        prev_internal_simd_struct[indiv_id]->nb_non_func_genes++;
                    }
                    if (prot->h > 0) {
                        prev_internal_simd_struct[indiv_id]->nb_genes_activ++;
                    } else {
                        prev_internal_simd_struct[indiv_id]->nb_genes_inhib++;
                    }
                }
            }
        }


        stats_best->write_best();
        stats_mean->write_average();


        if (standalone_ && exp_m_->record_light_tree()) {
            if (standalone_ && exp_m_->record_light_tree() && AeTime::time() % exp_m_->backup_step() == 0 &&
                AeTime::time() > 0) {
            }

            if (standalone_ && exp_m_->record_light_tree() && AeTime::time() > 0) {
                exp_m_->output_m()->light_tree()->update_tree(AeTime::time(), prev_internal_simd_struct);

                if (AeTime::time() % exp_m_->backup_step() == 0) {
                    std::cout << "writing light tree for gen : " << AeTime::time() << '\n';
                    exp_m_->output_m()->write_light_tree(AeTime::time());
                }
            }

            if (standalone_ && exp_m_->record_light_tree() && AeTime::time() % exp_m_->backup_step() == 0 &&
                AeTime::time() > 0) {
            }

        }


        if (standalone_ && exp_m_->record_tree() && AeTime::time() % exp_m_->output_m()->tree_step() == 0 &&
            AeTime::time() > 0) {
            printf("Tree SIMD backup\n");

            exp_m_->output_m()->write_tree(AeTime::time());
        }

        if (standalone_ && AeTime::time() % exp_m_->backup_step() == 0) {
            printf("Backup... OK\n");
            for (int indiv_id = 0; indiv_id < (int) exp_m_->nb_indivs(); indiv_id++) {
                int x = indiv_id / exp_m_->world()->height();
                int y = indiv_id % exp_m_->world()->height();

                exp_m_->world()->grid(x, y)->individual()->clear_everything_except_dna_and_promoters();
                exp_m_->world()->grid(x, y)->individual()->genetic_unit_list_nonconst().clear();
                delete exp_m_->world()->grid(x, y)->individual();

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
                int32_t nb_blocks_ = prev_internal_simd_struct[indiv_id]->dna_->length() / BLOCK_SIZE + 1;
                char *dna_string = new char[nb_blocks_ * BLOCK_SIZE];
                memset(dna_string, 0,
                       (prev_internal_simd_struct[indiv_id]->dna_->length() + 1) * sizeof(char));


                char *to_copy = prev_internal_simd_struct[indiv_id]->dna_->to_char();


                memcpy(dna_string, to_copy,
                       (prev_internal_simd_struct[indiv_id]->dna_->length() + 1) * sizeof(char));


                indiv->genetic_unit_list_.clear();
                indiv->add_GU(dna_string, prev_internal_simd_struct[indiv_id]->dna_->length());
                indiv->genetic_unit_nonconst(0).set_min_gu_length(exp_m_->exp_s()->min_genome_length());
                indiv->genetic_unit_nonconst(0).set_max_gu_length(exp_m_->exp_s()->max_genome_length());
                indiv->EvaluateInContext(exp_m_->world()->grid(x, y)->habitat());
                indiv->compute_statistical_data();

                exp_m_->world()->grid(x, y)->set_individual(indiv);
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
    }
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

        printf("Check DNA indiv %d %d %d --- NB Mutation %ld\n",i,dna_size[i],exp_m_->world()->grid(x, y)->individual()->genetic_unit(
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
            prev_internal_simd_struct[i]->metadata_->rna_count(),
            exp_m_->world()->grid(x, y)->individual()->rna_list().size(),
            prev_internal_simd_struct[i]->metadata_->proteins_count(),
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
    for (idx = 0; idx < (prev_internal_simd_struct[i]->metadata_->promoter_count()); idx++) {
        if (prev_internal_simd_struct[i]->metadata_->promoters(idx) != nullptr)
            printf("Promoters found at %d\n",prev_internal_simd_struct[i]->metadata_->promoters(idx)->pos);
    }

    idx = 0;
    for (idx = 0; idx < (prev_internal_simd_struct[i]->metadata_->rna_count()); idx++) {
        if (prev_internal_simd_struct[i]->metadata_->rnas(idx) != nullptr) {
            printf("RNA SIMD %d Start %d Stop %d Leading/Lagging %d Length %d\n", idx,
                   prev_internal_simd_struct[i]->metadata_->rnas(idx)->begin,
                   prev_internal_simd_struct[i]->metadata_->rnas(idx)->end,
                   prev_internal_simd_struct[i]->metadata_->rnas(idx)->leading_lagging,
                   prev_internal_simd_struct[i]->metadata_->rnas(idx)->length);
        }
    }

    int prot_cpt_b=0;
    idx = 0;
    for (auto prot : exp_m_->world()->grid(x, y)->old_one->protein_list()) {
        bool found = false;

        for (int pidx = 0; pidx <
                           (int) prev_internal_simd_struct[i]->metadata_->proteins_count(); pidx++) {
            if (prev_internal_simd_struct[i]->metadata_->proteins(pidx)->is_init_) {
                if ((prev_internal_simd_struct[i]->metadata_->proteins(pidx)->e == prot->concentration()) &&
                    (prev_internal_simd_struct[i]->metadata_->proteins(pidx)->protein_end == prot->last_STOP_base_pos())) {
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
                      (int) prev_internal_simd_struct[i]->metadata_->proteins_count(); idx++) {
        if (prev_internal_simd_struct[i]->metadata_->proteins(idx)->is_init_) {


            bool found = false;

            for (auto prot : exp_m_->world()->grid(x, y)->old_one->protein_list()) {
                if (( prev_internal_simd_struct[i]->metadata_->proteins(idx)->e ==  prot->concentration()) &&
                    ( prev_internal_simd_struct[i]->metadata_->proteins(idx)->protein_end ==  prot->last_STOP_base_pos())) {
                    found = true;
                    break;
                }
            }

            //for (idx = 0; idx < (int) (internal_simd_struct[i]->proteins.size()); idx++) {
            if (!found)
                printf("Proteins SIMD %d Start %d (end %d) Length %d Leading/Lagging %d M/W/H %f/%f/%f Func %d -- Concentration %f\n", idx,
                       prev_internal_simd_struct[i]->metadata_->proteins(idx)->protein_start,
                       prev_internal_simd_struct[i]->metadata_->proteins(idx)->protein_end,
                       prev_internal_simd_struct[i]->metadata_->proteins(idx)->protein_length,
                       prev_internal_simd_struct[i]->metadata_->proteins(idx)->leading_lagging,
                       prev_internal_simd_struct[i]->metadata_->proteins(idx)->m,
                       prev_internal_simd_struct[i]->metadata_->proteins(idx)->w,
                       prev_internal_simd_struct[i]->metadata_->proteins(idx)->h,
                       prev_internal_simd_struct[i]->metadata_->proteins(idx)->is_functional,
                       prev_internal_simd_struct[i]->metadata_->proteins(idx)->e
                );
            prot_cpt_b++;
        }
    }

}

void SIMD_Individual::check_struct() {
    for (int i = 0; i < (int) exp_m_->nb_indivs(); i++) {

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


    int prot_size = (int) exp_m_->world()->grid(x, y)->individual()->protein_list().size();


    if (i_fit_1 != i_fit_2 && dna_size[i] > 300)
    //if (i == 268) {
      if ((internal_simd_struct[i]->metadata_->rna_count() != exp_m_->world()->grid(x, y)->individual()->rna_list().size()) ||
        (internal_simd_struct[i]->metadata_->proteins_count() != prot_size)) {
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
          internal_simd_struct[i]->metadata_->rna_count(),
          exp_m_->world()->grid(x, y)->individual()->rna_list().size(),
          internal_simd_struct[i]->metadata_->proteins_count(),
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
        int prot_cpt_b=0;

      for (auto prot : exp_m_->world()->grid(x, y)->individual()->protein_list()) {
          bool found = false;

          for (int pidx = 0; pidx < internal_simd_struct[i]->metadata_->proteins_count(); pidx++) {
              if (internal_simd_struct[i]->metadata_->proteins(pidx)->is_init_) {
                  if ((internal_simd_struct[i]->metadata_->proteins(pidx)->e == prot->concentration()) &&
                      (internal_simd_struct[i]->metadata_->proteins(pidx)->protein_end == prot->last_STOP_base_pos())) {
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

      for (int idx = 0; idx < internal_simd_struct[i]->metadata_->proteins_count(); idx++) {
          if (internal_simd_struct[i]->metadata_->proteins(idx)->is_init_) {


          bool found = false;

          for (auto prot : exp_m_->world()->grid(x, y)->individual()->protein_list()) {
              if (( internal_simd_struct[i]->metadata_->proteins(idx)->e ==  prot->concentration()) &&
                  ( internal_simd_struct[i]->metadata_->proteins(idx)->protein_end ==  prot->last_STOP_base_pos())) {
                  found = true;
                  break;
              }
          }

      //for (idx = 0; idx < (int) (internal_simd_struct[i]->proteins.size()); idx++) {
        if (!found)
          printf("Proteins SIMD %d Start %d (end %d) Length %d Leading/Lagging %d M/W/H %f/%f/%f Func %d -- Concentration %f\n", idx,
               internal_simd_struct[i]->metadata_->proteins(idx)->protein_start,
                 internal_simd_struct[i]->metadata_->proteins(idx)->protein_end,
                 internal_simd_struct[i]->metadata_->proteins(idx)->protein_length,
                 internal_simd_struct[i]->metadata_->proteins(idx)->leading_lagging,
                 internal_simd_struct[i]->metadata_->proteins(idx)->m,
                 internal_simd_struct[i]->metadata_->proteins(idx)->w,
                 internal_simd_struct[i]->metadata_->proteins(idx)->h,
                 internal_simd_struct[i]->metadata_->proteins(idx)->is_functional,
                 internal_simd_struct[i]->metadata_->proteins(idx)->e
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
    Internal_SIMD_Struct::Internal_SIMD_Struct(ExpManager* exp_m, double w_max, SIMD_DnaFactory* dna_factory) {
        exp_m_ = exp_m;
        w_max_ = w_max;

        if (exp_m_->exp_s()->get_simd_metadata_flavor() == SIMDMetadataFlavor::STD_MAP)
            metadata_ = new SIMD_Map_Metadata(this);
        else if (exp_m_->exp_s()->get_simd_metadata_flavor() == SIMDMetadataFlavor::DYN_TAB)
            metadata_ = new SIMD_DynTab_Metadata(this);
        else if (exp_m_->exp_s()->get_simd_metadata_flavor() == SIMDMetadataFlavor::STD_LIST)
            metadata_ = new SIMD_List_Metadata(this);

        dna_factory_ = dna_factory;
    }


Internal_SIMD_Struct::Internal_SIMD_Struct(ExpManager* exp_m, Internal_SIMD_Struct* clone, SIMD_DnaFactory* dna_factory) {
    w_max_ = clone->w_max_;

  exp_m_ = exp_m;

  usage_count_ = 1;
  dna_ = dna_factory->get_dna(clone->dna_->length());
  //printf("DNA Factory -- %p %p\n",dna_,dna_->data_);
  dna_->set_indiv(clone->dna_,this);

  dna_factory_ = dna_factory;


  //promoters.resize(clone->promoters.size());
  if (exp_m_->exp_s()->get_simd_metadata_flavor() == SIMDMetadataFlavor::STD_MAP)
      metadata_ = new SIMD_Map_Metadata(this,dynamic_cast<SIMD_Map_Metadata*>(clone->metadata_));
  else if (exp_m_->exp_s()->get_simd_metadata_flavor() == SIMDMetadataFlavor::DYN_TAB)
      metadata_ = new SIMD_DynTab_Metadata(this,dynamic_cast<SIMD_DynTab_Metadata*>(clone->metadata_));
  else if (exp_m_->exp_s()->get_simd_metadata_flavor() == SIMDMetadataFlavor::STD_LIST)
      metadata_ = new SIMD_List_Metadata(this,dynamic_cast<SIMD_List_Metadata*>(clone->metadata_));

  fitness = clone->fitness;
  metaerror = clone->metaerror;
  //leading_prom_pos = clone->leading_prom_pos;
  //lagging_prom_pos = clone->lagging_prom_pos;

}

Internal_SIMD_Struct::~Internal_SIMD_Struct() {
        //printf("DESTRUCTOR SIMD INDIV -- DNA Factory %p --- DNA %p\n",dna_factory_,dna_);

  dna_factory_->give_back(dna_);

  delete metadata_;
}

/**
 * We need some index for the promoter optimization
 */
void Internal_SIMD_Struct::rebuild_index() {
        if (exp_m_->exp_s()->get_simd_metadata_flavor() == SIMDMetadataFlavor::STD_MAP)
            dynamic_cast<SIMD_Map_Metadata*>(metadata_)->rebuild_index();
}
}
