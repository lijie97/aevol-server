//
// Created by arrouan on 11/01/18.
//

#include "DnaMutator.h"

namespace aevol {
DnaMutator::DnaMutator(std::shared_ptr<JumpingMT> mut_prng,
      int32_t length,
      double duplication_rate,
      double deletion_rate,
      double translocation_rate,
      double inversion_rate,
      double point_mutation_rate,
      double small_insertion_rate,
      double small_deletion_rate,
      int16_t max_indel_size) {
  mut_prng_ = mut_prng;
  length_ = length;
  duplication_rate_ = duplication_rate;
  deletion_rate_ = deletion_rate;
  translocation_rate_ = translocation_rate;
  inversion_rate_ = inversion_rate;
  point_mutation_rate_ = point_mutation_rate;
  small_insertion_rate_ = small_insertion_rate;
  small_deletion_rate_ = small_deletion_rate;
  max_indel_size_ = max_indel_size;
}

DnaMutator::DnaMutator(Individual * indiv) {
  mut_prng_ = indiv->mut_prng();
  length_ = indiv->amount_of_dna();
  duplication_rate_ = indiv->duplication_rate();
  deletion_rate_ = indiv->deletion_rate();
  translocation_rate_ = indiv->translocation_rate();
  inversion_rate_ = indiv->inversion_rate();
  point_mutation_rate_ = indiv->point_mutation_rate();
  small_insertion_rate_ = indiv->small_insertion_rate();
  small_deletion_rate_ = indiv->small_deletion_rate();
  max_indel_size_ = indiv->max_indel_size();
}

void DnaMutator::generate_mutations() {
  generate_rearrangements();
  generate_small_mutations();
}


void DnaMutator::generate_small_mutations() {
  int32_t nb_swi = mut_prng_->
      binomial_random(length_, point_mutation_rate_);
  int32_t nb_ins = mut_prng_->
      binomial_random(length_, small_insertion_rate_);
  int32_t nb_del = mut_prng_->
      binomial_random(length_, small_deletion_rate_);
  int32_t nb_mut = nb_swi + nb_ins + nb_del;

  if (nb_mut > 0) {
    if (!hasMutate_) {
      hasMutate_ = true;
    }
  }

  int32_t random_value;

  for (int32_t i = nb_mut; i >= 1; i--) {
    random_value = mut_prng_->random(i);

    if (random_value < nb_swi) {
      int32_t pos = mut_prng_->random(length_);

      MutationEvent* mevent = new MutationEvent();
      mevent->switch_pos(pos);
      mutation_list_.push_back(mevent);

      nb_swi--;
    }
    else if (random_value < nb_swi + nb_ins) {
      int32_t pos = mut_prng_->random(length_);
      int16_t nb_insert;

      if (max_indel_size_ == 1) {
        nb_insert = 1;
      }
      else {
        nb_insert = 1 + mut_prng_->random(max_indel_size_);
      }

      MutationEvent* mevent = new MutationEvent();
      mevent->small_insertion(pos,nb_insert);

      mutation_list_.push_back(mevent);

      nb_ins--;
    }
    else { // (random_value >= nb_swi + nb_ins) => del
      int32_t pos = mut_prng_->random(length_);
      int16_t nb_del;
      if (max_indel_size_ == 1) {
        nb_del = 1;
      }
      else {
        nb_del = 1 + mut_prng_->random(max_indel_size_);
      }

      MutationEvent* mevent = new MutationEvent();
      mevent->small_deletion(pos,nb_del);

      mutation_list_.push_back(mevent);

      nb_del--;
    }
  }
}

void DnaMutator::generate_rearrangements() {
  int32_t nb_dupl = mut_prng_->
      binomial_random(length_, duplication_rate_);
  int32_t nb_del = mut_prng_->
      binomial_random(length_, deletion_rate_);
  int32_t nb_trans = mut_prng_->
      binomial_random(length_, translocation_rate_);
  int32_t nb_inv = mut_prng_->
      binomial_random(length_, inversion_rate_);
  int32_t nb_rear = nb_dupl + nb_del + nb_trans + nb_inv;


  if (nb_rear > 0) {
    if (!hasMutate_) {
      hasMutate_ = true;
    }
  }

  int32_t random_value;

  for (int32_t i = nb_rear; i >= 1; i--) {
    random_value = mut_prng_->random(i);

    if (random_value < nb_dupl) {
      int32_t pos_1, pos_2, pos_3;
      pos_1 = mut_prng_->random(length_);
      pos_2 = mut_prng_->random(length_);
      pos_3 = mut_prng_->random(length_);

      MutationEvent* mevent = new MutationEvent();
      mevent->duplication(pos_1,pos_2,pos_3);

      mutation_list_.push_back(mevent);

      nb_dupl--;  // Updating the urn (no replacement!)...
    }
    else if (random_value < nb_dupl + nb_del) {
      int32_t pos_1, pos_2;
      pos_1 = mut_prng_->random(length_);
      pos_2 = mut_prng_->random(length_);

      MutationEvent* mevent = new MutationEvent();
      mevent->deletion(pos_1,pos_2);

      mutation_list_.push_back(mevent);

      nb_del--;
    }
    else if (random_value < nb_dupl + nb_del + nb_trans) {
      int32_t pos_1, pos_2, pos_3, pos_4;
      int32_t segment_length;
      bool invert;

      pos_1 = mut_prng_->random(length_);
      pos_2 = mut_prng_->random(length_);
      if (pos_1 == pos_2) continue;

      // As it is commented in do_translocation(int32_t pos_1, int32_t pos_2,
      // int32_t pos_3, int32_t pos_4, bool invert), translocating segment
      // [pos_1, pos_2] is the same as translocating segment [pos_2, pos_1]
      // Since OriC must be at position 0, we will always translocate segment
      // [pos_1, pos_2] with pos_1 < pos_2
      if (pos_1 > pos_2) Utils::exchange(pos_1, pos_2);

      segment_length = pos_2 - pos_1;

      // Generate a position between pos_1 and pos_2
      pos_3 = pos_1 + mut_prng_->random(segment_length);

      // Generate a position that is NOT between pos_1 and pos_2
      pos_4 = mut_prng_->random(length_ - segment_length);
      if (pos_4 >= pos_1) pos_4 += segment_length;

      invert = (mut_prng_->random(2) == 0);

      MutationEvent* mevent = new MutationEvent();
      mevent->translocation(pos_1,pos_2,pos_3,pos_4,invert);

      mutation_list_.push_back(mevent);

      nb_trans--;
    }
    else {
      int32_t pos_1, pos_2;
      int32_t segment_length;
      pos_1 = mut_prng_->random(length_);
      pos_2 = mut_prng_->random(length_);

      if (pos_1 == pos_2) continue; // Invert everything <=> Invert nothing!

      // Invert the segment that don't contain OriC
      if (pos_1 > pos_2) Utils::exchange(pos_1, pos_2);

      MutationEvent* mevent = new MutationEvent();
      mevent->inversion(pos_1,pos_2);

      mutation_list_.push_back(mevent);

      nb_inv--;
    }
  }
}

}
