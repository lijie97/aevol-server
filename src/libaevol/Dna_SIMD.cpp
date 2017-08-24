//
// Created by Jonathan Rouzaud-Cornabas on 24/08/17.
//

#include "Dna_SIMD.h"

namespace aevol {
Dna_SIMD::Dna_SIMD(Dna* dna) {
  length_ = dna->length();
  nb_blocks_ = nb_blocks(length_);
  data_ = new char[nb_blocks_ * BLOCK_SIZE];
  memcpy(data_, dna->data(), (length_+1) * sizeof(char));
}

Dna_SIMD::~Dna_SIMD() {
  delete [] data;
}

void Dna_SIMD::remove(int32_t pos_1, int32_t pos_2) {
  assert(pos_1 >= 0 && pos_2 >= pos_1 && pos_2 <= length_);

  // Compute size of new genome
  int32_t new_length    = length_ - (pos_2 - pos_1);
  int32_t new_nb_blocks = nb_blocks(new_length);
  char*   new_genome    = new char[new_nb_blocks * BLOCK_SIZE * sizeof(char)];

  // Copy the remaining of the genome in tmp (preceeding and following parts)
  memcpy(new_genome, data_, pos_1 * sizeof(char));
  memcpy(&new_genome[pos_1], &data_[pos_2],
         (new_length - pos_1) * sizeof(char));
  new_genome[new_length] = '\0';

  // Replace previous genome with the new one
  delete [] data_;
  data_ = new_genome;

  // Update length data
  length_ = new_length;
  nb_blocks_ = new_nb_blocks;
}

void Dna_SIMD::insert(int32_t pos, const char* seq, int32_t seq_length) {
// Insert sequence 'seq' at position 'pos'
  assert(pos >= 0 && pos < length_);

  // If the sequence's length was not provided, compute it
  if (seq_length == -1) {
    seq_length = strlen(seq);
  }

  // Compute size of new genome
  int32_t new_length    = length_ + seq_length;
  int32_t new_nb_blocks = nb_blocks(new_length);
  char*   new_genome    = new char[new_nb_blocks * BLOCK_SIZE * sizeof(char)];

  // Build new genome from previous genome and sequence to insert
  memcpy(new_genome, data_, pos * sizeof(char));
  memcpy(&new_genome[pos], seq, seq_length * sizeof(char));
  memcpy(&new_genome[pos+seq_length], &data_[pos],
         (length_ - pos) * sizeof(char));
  new_genome[new_length] = '\0';

  // Replace the previous genome with the new one
  delete [] data_;
  data_ = new_genome;

  // Update length-related data
  length_ = new_length;
  nb_blocks_ = new_nb_blocks;
}

void Dna_SIMD::replace(int32_t pos, char* seq, int32_t seq_length) {
// Invert the sequence between positions 'first' and 'last'
  // Check pos value
  assert(pos >= 0 && pos < length_);

  // If the sequence's length was not provided, compute it
  if (seq_length == -1) {
    seq_length = strlen(seq);
  }

  // Check that the sequence is contiguous
  assert(pos + seq_length <= length_);

  // Perform the replacement
  memcpy(&data_[pos], seq, seq_length * sizeof(char));
}

}
