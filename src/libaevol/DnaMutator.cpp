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
      int16_t max_indel_size,
      int32_t min_genome_length, int32_t max_genome_length) {
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

  min_genome_length_ = min_genome_length;
  max_genome_length_ = max_genome_length;
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

  min_genome_length_ = indiv->min_genome_length();
  max_genome_length_ = indiv->max_genome_length();

  id_ = indiv->id();
}


void DnaMutator::generate_mutations() {
  generate_rearrangements();
  generate_small_mutations();

/*  if (nb_rear_+nb_mut_ > 0)
    printf("%d -- DNAMutator -- Dupli %d Large_Del %d Trans %d Inv %d Swi %d Ins %d Del %d\n",
      id_,
      nb_large_dupl_,nb_large_del_,nb_large_trans_,nb_large_inv_,nb_swi_,
      nb_ins_,nb_del_);*/
}

void DnaMutator::generate_rearrangements() {
  nb_large_dupl_ = mut_prng_->
      binomial_random(length_, duplication_rate_);
  nb_large_del_ = mut_prng_->
      binomial_random(length_, deletion_rate_);
  nb_large_trans_ = mut_prng_->
      binomial_random(length_, translocation_rate_);
  nb_large_inv_ = mut_prng_->
      binomial_random(length_, inversion_rate_);
  nb_rear_ = nb_large_dupl_ + nb_large_del_ + nb_large_trans_ + nb_large_inv_;
  cpt_rear_ = nb_rear_;

  if (nb_rear_ > 0) {
    if (!hasMutate_) {
      hasMutate_ = true;
    }
  }
}

void DnaMutator::generate_small_mutations() {
  nb_swi_ = mut_prng_->
      binomial_random(length_, point_mutation_rate_);
  nb_ins_ = mut_prng_->
      binomial_random(length_, small_insertion_rate_);
  nb_del_ = mut_prng_->
      binomial_random(length_, small_deletion_rate_);
  nb_mut_ = nb_swi_ + nb_ins_ + nb_del_;
  cpt_mut_ = nb_mut_;

  if (nb_mut_ > 0) {
    if (!hasMutate_) {
      hasMutate_ = true;
    }
  }
}

MutationEvent* DnaMutator::generate_next_mutation(int32_t length) {
  length_ = length;
  int32_t random_value;
  MutationEvent* mevent = nullptr;

  if (cpt_rear_>0) {
    random_value = mut_prng_->random(nb_rear_);
    cpt_rear_--;

    if (random_value < nb_large_dupl_) {
      nb_large_dupl_--;  // Updating the urn (no replacement!)...

      int32_t pos_1, pos_2, pos_3;
      pos_1 = mut_prng_->random(length_);
      pos_2 = mut_prng_->random(length_);
      pos_3 = mut_prng_->random(length_);

      int32_t genome_size_after = length_ + Utils::mod(pos_2 - pos_1 - 1, length_) + 1;
      if (genome_size_after > max_genome_length_)
        return nullptr;

      mevent = new MutationEvent();
      mevent->duplication(pos_1,pos_2,pos_3);

      mutation_list_.push_back(mevent);

    }
    else if (random_value < nb_large_dupl_ + nb_large_del_) {
      nb_large_del_--;

      int32_t pos_1, pos_2;
      pos_1 = mut_prng_->random(length_);
      pos_2 = mut_prng_->random(length_);


      int32_t genome_size_after = length_ - Utils::mod(pos_2 - pos_1 - 1, length_) + 1;

      if (genome_size_after < min_genome_length_)
        return nullptr;

      mevent = new MutationEvent();
      mevent->deletion(pos_1,pos_2);

      mutation_list_.push_back(mevent);

    }
    else if (random_value < nb_large_dupl_ + nb_large_del_ + nb_large_trans_) {
      nb_large_trans_--;

      int32_t pos_1, pos_2, pos_3, pos_4;
      int32_t segment_length;
      bool invert;

      pos_1 = mut_prng_->random(length_);
      pos_2 = mut_prng_->random(length_);
      if (pos_1 == pos_2) return nullptr;

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

      mevent = new MutationEvent();
      mevent->translocation(pos_1,pos_2,pos_3,pos_4,invert);

      mutation_list_.push_back(mevent);

    }
    else {
      nb_large_inv_--;

      int32_t pos_1, pos_2;
      pos_1 = mut_prng_->random(length_);
      pos_2 = mut_prng_->random(length_);

      if (pos_1 == pos_2) return nullptr; // Invert everything <=> Invert nothing!

      // Invert the segment that don't contain OriC
      if (pos_1 > pos_2) Utils::exchange(pos_1, pos_2);

      mevent = new MutationEvent();
      mevent->inversion(pos_1,pos_2);

      mutation_list_.push_back(mevent);

    }
  } else if (cpt_mut_>0) {
    random_value = mut_prng_->random(nb_mut_);

    cpt_mut_--;

    if (random_value < nb_swi_) {
      nb_swi_--;

      int32_t pos = mut_prng_->random(length_);

      MutationEvent* mevent = new MutationEvent();
      mevent->switch_pos(pos);
      mutation_list_.push_back(mevent);

    }
    else if (random_value < nb_swi_ + nb_ins_) {
      nb_ins_--;

      int32_t pos = mut_prng_->random(length_);
      int16_t nb_insert;

      if (max_indel_size_ == 1) {
        nb_insert = 1;
      }
      else {
        nb_insert = 1 + mut_prng_->random(max_indel_size_);
      }

      if (length_ + nb_insert > max_genome_length_) {
        return nullptr;
      }

      // Insert the sequence
      // Prepare the sequence to be inserted
#ifdef WITH_BITSET
      BitSet_SIMD* inserted_seq = new BitSet_SIMD(nb_insert);
#else
      char* inserted_seq = new char[nb_insert + 1];
#endif
      char inserted_char;
      for (int16_t j = 0; j < nb_insert; j++) {
        inserted_char = static_cast<char>('0' + mut_prng_->random(NB_BASE));
#ifdef WITH_BITSET
        if (inserted_char == '0')
#ifdef _DYNAMIC_CUSTOM_BITSET
          inserted_seq->set_to_0(j);
#elif _STATIC_BITSET
          inserted_seq->data_.set(j,false);
#endif
        else
#ifdef _DYNAMIC_CUSTOM_BITSET
          inserted_seq->set_to_1(j);
#elif _STATIC_BITSET
          inserted_seq->data_.set(j,true);
#endif
#else
        inserted_seq[j] = inserted_char;
#endif
      }
#ifndef WITH_BITSET
      inserted_seq[nb_insert] = '\0';
#endif

      mevent = new MutationEvent();
#ifdef WITH_BITSET
      mevent->small_insertion(pos,inserted_seq);
#else
      mevent->small_insertion(pos,nb_insert,inserted_seq);
#endif

      mutation_list_.push_back(mevent);

    }
    else { // (random_value >= nb_swi + nb_ins) => del
      nb_del_--;

      int32_t pos = mut_prng_->random(length_);
      int16_t nb_del;
      if (max_indel_size_ == 1) {
        nb_del_ = 1;
      }
      else {
        nb_del_ = 1 + mut_prng_->random(max_indel_size_);
      }

      if (length_ - nb_del < min_genome_length_)
        return nullptr;

      mevent = new MutationEvent();
      mevent->small_deletion(pos,nb_del);

      mutation_list_.push_back(mevent);

    }
  }

  return mevent;
}


}
