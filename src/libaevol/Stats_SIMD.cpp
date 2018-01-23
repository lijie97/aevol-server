//
// Created by arrouan on 19/01/18.
//

#include "Stats_SIMD.h"
#include "Dna_SIMD.h"
#include "SIMD_Individual.h"

namespace aevol {

Stats_SIMD::Stats_SIMD(SIMD_Individual* simd_individual, int64_t generation, bool best_or_not) {
  simd_individual_ = simd_individual;
  is_indiv_ = best_or_not;
  generation_ = generation;


  if (generation_==1) {
    statfile_best_.open("stats/stats_simd_best.csv",std::ofstream::trunc);
    statfile_mean_.open("stats/stats_simd_mean.csv",std::ofstream::trunc);
    if (is_indiv_)
      statfile_best_<<"Generation"<<","<<"fitness"<<","<<"metabolic_error"<<","<<
                  "amount_of_dna"<<","<<"nb_coding_rnas"<<","<<"nb_non_coding_rnas"<<","<<
                  "nb_functional_genes"<<","<<"nb_non_functional_genes"<<","<<"nb_mut"
                  <<","<<"nb_switch"<<","<<"nb_indels"<<","<<"nb_rear"<<","<<"nb_dupl"<<","<<
                  "nb_del"<<","<<"nb_trans"<<","<<"nb_inv"<<","<<"dupl_rate"<<","<<"del_rate"
                  <<","<<"trans_rate"<<","<<"inv_rate"
                  <<std::endl;
    else
      statfile_mean_<<"Generation"<<","<<"fitness"<<","<<"metabolic_error"<<","<<
                    "amount_of_dna"<<","<<"nb_coding_rnas"<<","<<"nb_non_coding_rnas"<<","<<
                    "nb_functional_genes"<<","<<"nb_non_functional_genes"<<","<<"nb_mut"
                    <<","<<"nb_switch"<<","<<"nb_indels"<<","<<"nb_rear"<<","<<"nb_dupl"<<","<<
                    "nb_del"<<","<<"nb_trans"<<","<<"nb_inv"<<","<<"dupl_rate"<<","<<"del_rate"
                    <<","<<"trans_rate"<<","<<"inv_rate"
                    <<std::endl;
  } else {
    statfile_best_.open("stats/stats_simd_best.csv",std::ofstream::app);
    statfile_mean_.open("stats/stats_simd_mean.csv",std::ofstream::app);
  }

}

void Stats_SIMD::compute_best() {
  is_indiv_ = true;

  fitness_ = simd_individual_->best_indiv->fitness;
  metabolic_error_  = simd_individual_->best_indiv->metaerror;

  amount_of_dna_ = simd_individual_->best_indiv->dna_->length();

  for (int i = 0; i < simd_individual_->best_indiv->rnas.size(); i++) {
    if (simd_individual_->best_indiv->rnas[i].is_coding_)
      nb_coding_rnas_++;
    else
      nb_non_coding_rnas_++;
  }

  for (int i = 0; i < simd_individual_->best_indiv->proteins.size(); i++) {
    if (simd_individual_->best_indiv->proteins[i].is_functional) {
      nb_functional_genes_++;
    } else {
      nb_non_functional_genes_++;
    }
  }

  nb_mut_ = simd_individual_->best_indiv->dna_->nb_mut_;
  nb_rear_ = simd_individual_->best_indiv->dna_->nb_rear_;
  nb_switch_ = simd_individual_->best_indiv->dna_->nb_swi_;
  nb_indels_ = simd_individual_->best_indiv->dna_->nb_indels_;
  nb_dupl_ = simd_individual_->best_indiv->dna_->nb_large_dupl_;
  nb_del_ = simd_individual_->best_indiv->dna_->nb_large_del_;
  nb_trans_ = simd_individual_->best_indiv->dna_->nb_large_trans_;
  nb_inv_ = simd_individual_->best_indiv->dna_->nb_large_inv_;

  dupl_rate_  = nb_dupl_  / simd_individual_->best_indiv->dna_->parent_length();
  del_rate_   = nb_del_   / simd_individual_->best_indiv->dna_->parent_length();
  trans_rate_ = nb_trans_ / simd_individual_->best_indiv->dna_->parent_length();
  inv_rate_   = nb_inv_   / simd_individual_->best_indiv->dna_->parent_length();

/*  nb_bases_in_0_CDS_;
  nb_bases_in_0_functional_CDS_;
  nb_bases_in_0_non_functional_CDS_;
  nb_bases_in_0_RNA_;
  nb_bases_in_0_coding_RNA_;
  nb_bases_in_0_non_coding_RNA_;

  nb_bases_non_essential_;
  nb_bases_non_essential_including_nf_genes_;*/
  is_computed_ = true;
}

void Stats_SIMD::compute_average() {
  is_indiv_ = false;
  pop_size_ = simd_individual_->nb_indivs_;

  for (int indiv_id = 0; indiv_id < pop_size_; indiv_id++) {
    fitness_ += simd_individual_->prev_internal_simd_struct[indiv_id]->fitness;
    metabolic_error_ += simd_individual_->prev_internal_simd_struct[indiv_id]->metaerror;

    amount_of_dna_ += simd_individual_->prev_internal_simd_struct[indiv_id]->dna_->length();

    for (int i = 0; i < simd_individual_->prev_internal_simd_struct[indiv_id]->rnas.size(); i++) {
      if (simd_individual_->prev_internal_simd_struct[indiv_id]->rnas[i].is_coding_)
        nb_coding_rnas_++;
      else
        nb_non_coding_rnas_++;
    }

    for (int i = 0; i < simd_individual_->prev_internal_simd_struct[indiv_id]->proteins.size(); i++) {
      if (simd_individual_->prev_internal_simd_struct[indiv_id]->proteins[i].is_functional) {
        nb_functional_genes_++;
      } else {
        nb_non_functional_genes_++;
      }
    }

    nb_mut_ += simd_individual_->prev_internal_simd_struct[indiv_id]->dna_->nb_mut_;
    nb_rear_ += simd_individual_->prev_internal_simd_struct[indiv_id]->dna_->nb_rear_;
    nb_switch_ += simd_individual_->prev_internal_simd_struct[indiv_id]->dna_->nb_swi_;
    nb_indels_ += simd_individual_->prev_internal_simd_struct[indiv_id]->dna_->nb_indels_;
    nb_dupl_ += simd_individual_->prev_internal_simd_struct[indiv_id]->dna_->nb_large_dupl_;
    nb_del_ += simd_individual_->prev_internal_simd_struct[indiv_id]->dna_->nb_large_del_;
    nb_trans_ += simd_individual_->prev_internal_simd_struct[indiv_id]->dna_->nb_large_trans_;
    nb_inv_ += simd_individual_->prev_internal_simd_struct[indiv_id]->dna_->nb_large_inv_;

    dupl_rate_ += nb_dupl_ / simd_individual_->prev_internal_simd_struct[indiv_id]->dna_->parent_length();
    del_rate_ += nb_del_ / simd_individual_->prev_internal_simd_struct[indiv_id]->dna_->parent_length();
    trans_rate_ +=
        nb_trans_ / simd_individual_->prev_internal_simd_struct[indiv_id]->dna_->parent_length();
    inv_rate_ += nb_inv_ / simd_individual_->prev_internal_simd_struct[indiv_id]->dna_->parent_length();

/*  nb_bases_in_0_CDS_;
  nb_bases_in_0_functional_CDS_;
  nb_bases_in_0_non_functional_CDS_;
  nb_bases_in_0_RNA_;
  nb_bases_in_0_coding_RNA_;
  nb_bases_in_0_non_coding_RNA_;

  nb_bases_non_essential_;
  nb_bases_non_essential_including_nf_genes_;*/
  }

  fitness_ /= pop_size_;
  metabolic_error_ /= pop_size_;

  amount_of_dna_ /= pop_size_;
  nb_coding_rnas_ /= pop_size_;
  nb_non_coding_rnas_ /= pop_size_;

  nb_functional_genes_ /= pop_size_;
  nb_non_functional_genes_ /= pop_size_;

  nb_mut_ /= pop_size_;
  nb_rear_ /= pop_size_;
  nb_switch_ /= pop_size_;
  nb_indels_ /= pop_size_;
  nb_dupl_ /= pop_size_;
  nb_del_ /= pop_size_;
  nb_trans_ /= pop_size_;
  nb_inv_ /= pop_size_;

  dupl_rate_ /= pop_size_;
  del_rate_ /= pop_size_;
  trans_rate_ /= pop_size_;
  inv_rate_ /= pop_size_;

  is_computed_ = true;
}


void Stats_SIMD::write_best() {
  if (is_indiv_ && !is_computed_)
    compute_best();

  if (is_indiv_ && is_computed_) {
    // Write best stats
    statfile_best_<<generation_<<","<<fitness_<<","<<metabolic_error_<<","<<
        amount_of_dna_<<","<<nb_coding_rnas_<<","<<nb_non_coding_rnas_<<","<<
        nb_functional_genes_<<","<<nb_non_functional_genes_<<","<<nb_mut_
        <<","<<nb_switch_<<","<<nb_indels_<<","<<nb_rear_<<","<<nb_dupl_<<","<<
        nb_del_<<","<<nb_trans_<<","<<nb_inv_<<","<<dupl_rate_<<","<<del_rate_
        <<","<<trans_rate_<<","<<inv_rate_
                  <<std::endl;
    statfile_best_.flush();
  }
}

void Stats_SIMD::write_average() {
  if (!is_indiv_ && !is_computed_)
    compute_average();

  if (!is_indiv_ && is_computed_) {
    // Write average stats
    statfile_mean_<<generation_<<","<<fitness_<<","<<metabolic_error_<<","<<
                  amount_of_dna_<<","<<nb_coding_rnas_<<","<<nb_non_coding_rnas_<<","<<
                  nb_functional_genes_<<","<<nb_non_functional_genes_<<","<<nb_mut_
                  <<","<<nb_switch_<<","<<nb_indels_<<","<<nb_rear_<<","<<nb_dupl_<<","<<
                  nb_del_<<","<<nb_trans_<<","<<nb_inv_<<","<<dupl_rate_<<","<<del_rate_
                  <<","<<trans_rate_<<","<<inv_rate_
                  <<std::endl;
  }
}

}
