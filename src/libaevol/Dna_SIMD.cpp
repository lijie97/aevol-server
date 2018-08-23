//
// Created by Jonathan Rouzaud-Cornabas on 24/08/17.
//


#include "SIMD_Individual.h"
#include "ExpManager.h"
#include "Dna_SIMD.h"
#include "DnaMutator.h"
#include "BitSet_SIMD.h"
#include "Individual.h"


namespace aevol {
Dna_SIMD::Dna_SIMD(Dna* dna, Internal_SIMD_Struct* indiv) {
#ifdef WITH_BITSET
  //printf("Building BitSet\n");
  bitset_ = new BitSet_SIMD(dna->data(),dna->length());
  //printf("BitSet is BUILT\n");
#else
  length_ = dna->length();

  nb_blocks_ = nb_blocks(length_);// _mm_free(ptr)
  posix_memalign((void **)&data_ ,64, nb_blocks_ * BLOCK_SIZE* sizeof(char)); //new char[nb_blocks_ * BLOCK_SIZE];
  memset(data_,0,(length_+1) * sizeof(char));
  memcpy(data_, dna->data(), (length_+1) * sizeof(char));
#endif
  parent_length_ = dna->length();
  indiv_ = indiv;
}


Dna_SIMD::Dna_SIMD(Dna_SIMD* dna, Internal_SIMD_Struct* indiv, bool copy_dna) {
#ifdef WITH_BITSET
  if (copy_dna)
    bitset_ = new BitSet_SIMD(dna->bitset_);
  else
    bitset_ = nullptr;
  parent_length_ = dna->bitset_->length_;
#else
  length_ = dna->length_;

  nb_blocks_ = nb_blocks(length_);
  posix_memalign((void **)&data_ ,64, nb_blocks_ * BLOCK_SIZE* sizeof(char));//new char[nb_blocks_ * BLOCK_SIZE];
  memset(data_,0,(length_+1) * sizeof(char));
  memcpy(data_, dna->data_, (length_+1) * sizeof(char));

  parent_length_ = dna->length_;
#endif
  indiv_ = indiv;
}

Dna_SIMD::Dna_SIMD(Dna* dna) {
#ifdef WITH_BITSET
 // printf("Building BitSet\n");
  bitset_ = new BitSet_SIMD(dna->data(),dna->length());
  //printf("BitSet is BUILT\n");
  //int a=1;
  //printf("BitSet is BUILT...\n");
#else
  length_ = dna->length();

  nb_blocks_ = nb_blocks(length_);
  posix_memalign((void **)&data_ ,64, nb_blocks_ * BLOCK_SIZE* sizeof(char));//new char[nb_blocks_ * BLOCK_SIZE];
  memset(data_,0,(length_+1) * sizeof(char));
  memcpy(data_, dna->data(), (length_+1) * sizeof(char));
#endif

  //printf("Calling LENGTH\n");
  parent_length_ = dna->length();
  //printf("Dna_SIMD out of constructor\n");
}

Dna_SIMD::~Dna_SIMD() {
#ifdef WITH_BITSET
  delete bitset_;
#else
  if (data_ != nullptr) { free(data_); data_=nullptr;}
#endif

if (mutation_list.size() > 0) {
    for (auto repl : mutation_list)
        delete repl;
    mutation_list.clear();
}


}

#ifndef WITH_BITSET
void Dna_SIMD::remove(int32_t pos_1, int32_t pos_2) {
//  printf("%d -- REMOVE POS 1 %d POS 2 %d LENGTH %d\n",indiv_->indiv_id,pos_1,pos_2,length());

  assert(pos_1 >= 0 && pos_2 >= pos_1 && pos_2 <= length_);

  // Compute size of new genome
  int32_t new_length    = length_ - (pos_2 - pos_1);
  int32_t new_nb_blocks = nb_blocks(new_length);

//  printf("POS 1 %d POS 2 %d  LENGTH %d NEW Length %d Nb Block %d\n",
//         pos_1,pos_2,length(),new_length,new_nb_blocks);

  char*   new_genome;
  posix_memalign((void **)&new_genome,64,nb_blocks_ * BLOCK_SIZE* sizeof(char));//new char[new_nb_blocks * BLOCK_SIZE];

  memset(new_genome,0,(new_length) * sizeof(char));

  // Copy the remaining of the genome in tmp (preceeding and following parts)
  memcpy(new_genome, data_, pos_1 * sizeof(char));
  memcpy(&new_genome[pos_1], &data_[pos_2],
         (new_length - pos_1) * sizeof(char));
  new_genome[new_length] = '\0';

  // Replace previous genome with the new one
  free(data_);
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
  char*   new_genome;
  posix_memalign((void **)&new_genome ,64, new_nb_blocks * BLOCK_SIZE* sizeof(char));//new char[new_nb_blocks * BLOCK_SIZE * sizeof(char)];

  memset(new_genome,0,(new_length) * sizeof(char));

  // Build new genome from previous genome and sequence to insert
  memcpy(new_genome, data_, pos * sizeof(char));
  memcpy(&new_genome[pos], seq, seq_length * sizeof(char));
  memcpy(&new_genome[pos+seq_length], &data_[pos],
         (length_ - pos) * sizeof(char));
  new_genome[new_length] = '\0';

  // Replace the previous genome with the new one
  free(data_);
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
#endif

bool Dna_SIMD::do_switch(int32_t pos) {
#ifdef WITH_BITSET
  bitset_->flip(pos);
#else
  // Perform the mutation
  if (data_[pos] == '0') data_[pos] = '1';
  else data_[pos] = '0';
#endif

  // Remove promoters containing the switched base
  indiv_->remove_promoters_around(pos, Utils::mod(pos + 1, length()));

  // Look for potential new promoters containing the switched base
  if (length() >= PROM_SIZE)
    indiv_->look_for_new_promoters_around(pos, Utils::mod(pos + 1, length()));

  if (SIMD_Individual::standalone_simd) {
      PointMutation *mut = new PointMutation(pos);
/*      printf("SEND MUTATION EVENT (%ld)\n",indiv_->observers_[MUTATION].size());*/
      indiv_->notifyObservers(MUTATION, mut);
      delete mut;
  }

  return true;
}

#ifdef WITH_BITSET
bool Dna_SIMD::do_small_insertion(int32_t pos, BitSet_SIMD* seq) {
  // Remove the promoters that will be broken
  indiv_->remove_promoters_around(pos);

  bitset_->insert_at(seq,pos);

  // Look for new promoters
  if (length() >= PROM_SIZE) {
    if (length() - seq->length_ < PROM_SIZE) {
      // Special case where the genome was smaller than a promoter before the
      // insertion and greater than (or as big as) a promoter after the
      // insertion.
      // In that case, we must look for new promoters thoroughly on the whole
      // genome using locate_promoters
      indiv_->locate_promoters();
    }
    else {
      indiv_->move_all_promoters_after(pos, seq->length_);
      indiv_->look_for_new_promoters_around(pos, Utils::mod(pos + seq->length_,
                                                            length()));
    }
  }

  if (SIMD_Individual::standalone_simd) {
    SmallInsertion* mut = new SmallInsertion(pos, seq->length_, seq->to_char());
    indiv_->notifyObservers(MUTATION, mut);
    delete mut;
  }

  return true;
}
#else
bool Dna_SIMD::do_small_insertion(int32_t pos, int16_t nb_insert, char* seq) {
  // Remove the promoters that will be broken
  indiv_->remove_promoters_around(pos);

  insert(pos, seq, nb_insert);

  // Look for new promoters
  if (length_ >= PROM_SIZE) {
    if (length_ - nb_insert < PROM_SIZE) {
      // Special case where the genome was smaller than a promoter before the
      // insertion and greater than (or as big as) a promoter after the
      // insertion.
      // In that case, we must look for new promoters thoroughly on the whole
      // genome using locate_promoters
      indiv_->locate_promoters();
    }
    else {
      indiv_->move_all_promoters_after(pos, nb_insert);
      indiv_->look_for_new_promoters_around(pos, Utils::mod(pos + nb_insert,
                                                               length_));
    }
  }

    if (SIMD_Individual::standalone_simd) {
        SmallInsertion *mut = new SmallInsertion(pos, nb_insert, seq);
        indiv_->notifyObservers(MUTATION, mut);
        delete mut;
    }

  return true;
}
#endif

bool Dna_SIMD::do_small_deletion(int32_t pos, int16_t nb_del) {
  // Remove promoters containing at least one nucleotide from the sequence to
  // delete
  indiv_->remove_promoters_around(pos, Utils::mod(pos + nb_del, length()));

  // Do the deletion and update promoter list
  if (pos + nb_del <= length()) { // the deletion does not contain the origin of
    // replication
    // Do the deletion
#ifdef WITH_BITSET
    bitset_->remove(pos,pos+nb_del);
#else
    remove(pos, pos + nb_del);
#endif

    // Update promoter list
    if (length() >= PROM_SIZE) {
      indiv_->move_all_promoters_after(pos, -nb_del);
      indiv_->look_for_new_promoters_around(Utils::mod(pos, length()));
    }
  }
  else { // the deletion contains the origin of replication
    // Do the deletion
    int32_t nb_del_at_pos_0 = nb_del - length() + pos;

#ifdef WITH_BITSET
    bitset_->remove(pos, length());
    bitset_->remove(0, nb_del_at_pos_0);
#else
    remove(pos, length_);
    remove(0, nb_del_at_pos_0);
#endif

    pos -= nb_del_at_pos_0;

    // Update promoter list
    if (length() >= PROM_SIZE) {
      indiv_->move_all_promoters_after(0, -nb_del_at_pos_0);
      indiv_->look_for_new_promoters_around(0);
    }
  }

    if (SIMD_Individual::standalone_simd) {
        SmallDeletion *mut = new SmallDeletion(pos, nb_del);
        indiv_->notifyObservers(MUTATION, mut);
        delete mut;
    }

  return true;
}


bool Dna_SIMD::do_duplication(int32_t pos_1, int32_t pos_2, int32_t pos_3) {
// Duplicate segment [pos_1; pos_2[ and insert the duplicate before pos_3
#ifdef WITH_BITSET
  BitSet_SIMD *duplicate_segment = nullptr;
#else
  char* duplicate_segment = NULL;
#endif
  int32_t seg_length;
    //printf("Mutation is %d %d %d -- %d\n",pos_1,pos_2,pos_3,length());

  if (pos_1 < pos_2) {
    //
    //       pos_1         pos_2                   -> 0-
    //         |             |                   -       -
    // 0--------------------------------->      -         -
    //         ===============                  -         - pos_1
    //           tmp (copy)                      -       -
    //                                             -----      |
    //                                             pos_2    <-'
    //
    seg_length = pos_2 - pos_1;
#ifdef WITH_BITSET
    duplicate_segment = bitset_->duplicate(pos_1,pos_2);
#else
    posix_memalign((void **)&duplicate_segment,64,(seg_length+1)* sizeof(char));//new char[seg_length + 1];
    memcpy(duplicate_segment, &data_[pos_1], seg_length);
    duplicate_segment[seg_length] = '\0';
#endif
  }
  else { // if (pos_1 >= pos_2)
    // The segment to duplicate includes the origin of replication.
    // The copying process will be done in two steps.
    //
    //                                            ,->
    //    pos_2                 pos_1            |      -> 0-
    //      |                     |                   -       - pos_2
    // 0--------------------------------->     pos_1 -         -
    // ======                     =======            -         -
    //  tmp2                        tmp1              -       -
    //                                                  -----
    //
    //
    int32_t tmp1_len = length() - pos_1;
    int32_t tmp2_len = pos_2;
    seg_length = tmp1_len + tmp2_len;

#ifdef WITH_BITSET
    duplicate_segment = bitset_->duplicate(pos_1,pos_1+tmp1_len);
    BitSet_SIMD* append_seg = bitset_->duplicate(0,pos_2);
    duplicate_segment->append(append_seg);
    delete append_seg;
#else
    posix_memalign((void **)&duplicate_segment,64,(seg_length+1)* sizeof(char));//new char[seg_length + 1];
    memcpy(duplicate_segment, &data_[pos_1], tmp1_len);     // Copy tmp1
    memcpy(&duplicate_segment[tmp1_len], data_, tmp2_len);  // Copy tmp2
    duplicate_segment[seg_length] = '\0';
#endif
  }

/*  if (indiv_->indiv_id==49) {
    printf(
        "Duplication _lengh %d\n", seg_length);
  }*/

  if (seg_length <= 0) {
#ifdef WITH_BITSET
    delete duplicate_segment;
#else
    free(duplicate_segment);
#endif
    return true;
  }

  // Create a copy of the promoters beared by the segment to be duplicated
  // (they will be inserted in the individual's RNA list later)
  std::vector<std::list<promoterStruct*>> duplicated_promoters = {{},
                                            {}};
  indiv_->duplicate_promoters_included_in(pos_1, pos_2,
                                             duplicated_promoters);

 /* if (indiv_->indiv_id == 433) {
    printf("DUPLICATE : Leading promoters lists : ");
    for (auto it : indiv_->leading_prom_pos) {
      printf("%d (%d) || ", it.first, it.second);
    }
    printf("\n");

    printf("DUPLICATE : Lagging promoters lists : ");
    for (auto it : indiv_->lagging_prom_pos) {
      printf("%d (%d) || ", it.first, it.second);
    }
    printf("\n");

    printf("DUPLICATE : Leading promoters lists (promoters): ");
    for (auto it : indiv_->promoters) {
      printf("%d (%d) -- ", it.second->pos, it.first);
    }

    printf("\n");

    printf("DUPLICATE : Leading DUPLICATED promoters lists : ");
    for (auto it : duplicated_promoters[LEADING]) {
      printf("%d || ", it->pos);
    }
    printf("\n");

    printf("DUPLICATE : Lagging DUPLICATED promoters lists: ");
    for (auto it : duplicated_promoters[LAGGING]) {
      printf("%d -- ", it->pos);
    }

    printf("\n");

    if (indiv_->promoters.size() != indiv_->leading_prom_pos.size() + indiv_->lagging_prom_pos.size()) {
        printf("Error unsynchronized cache !\n");
    }
  }*/

  indiv_->remove_promoters_around(pos_3);

/*  if (indiv_->indiv_id == 433) {
    printf("REMOVE : Leading promoters lists : ");
    for (auto it : indiv_->leading_prom_pos) {
      printf("%d (%d) || ", it.first, it.second);
    }
    printf("\n");
    printf("REMOVE : Lagging promoters lists : ");
    for (auto it : indiv_->lagging_prom_pos) {
      printf("%d (%d) || ", it.first, it.second);
    }
    printf("\n");
    printf("REMOVE : Leading promoters lists (promoters): ");
    for (auto it : indiv_->promoters) {
      printf("%d (%d) -- ", it.second->pos, it.first);
    }

    printf("\n");
      if (indiv_->promoters.size() != indiv_->leading_prom_pos.size() + indiv_->lagging_prom_pos.size()) {
          printf("Error unsynchronized cache !\n");
      }
  }*/

#ifdef WITH_BITSET
  bitset_->insert_at(duplicate_segment,pos_3);
#else
  insert(pos_3, duplicate_segment, seg_length);
#endif

  if (length() >= PROM_SIZE) {
    if (length() - seg_length < PROM_SIZE) {
      // Special case where the genome was smaller than a promoter before
      // the insertion and greater than (or as big as) a promoter after the
      // insertion.
      // In that case, we must look for new promoters thoroughly on the whole
      // genome using locate_promoters
      indiv_->locate_promoters();
    }
    else {
      indiv_->move_all_promoters_after(pos_3, seg_length);
      /*if (indiv_->indiv_id == 433) {
        printf("MOVE : Leading promoters lists : ");
        for (auto it : indiv_->leading_prom_pos) {
          printf("%d (%d) || ", it.first, it.second);
        }
        printf("\n");
        printf("MOVE : Lagging promoters lists : ");
        for (auto it : indiv_->lagging_prom_pos) {
          printf("%d (%d) || ", it.first, it.second);
        }
        printf("\n");
        printf("MOVE : Leading promoters lists (promoters): ");
        for (auto it : indiv_->promoters) {
          printf("%d (%d) -- ", it.second->pos, it.first);
        }

        printf("\n");

          printf("MOVE : Leading duplicated promoters lists : ");
          for (auto it : duplicated_promoters[LEADING]) {
            printf("%d || ", it->pos);
          }
          printf("\n");

          printf("MOVE : Lagging duplicated promoters lists (promoters): ");
          for (auto it : duplicated_promoters[LAGGING]) {
            printf("%d -- ", it->pos);
          }

          printf("\n");
          if (indiv_->promoters.size() != indiv_->leading_prom_pos.size() + indiv_->lagging_prom_pos.size()) {
              printf("Error unsynchronized cache !\n");
          }
      }*/

      indiv_->insert_promoters_at(duplicated_promoters, pos_3);

      /*if (indiv_->indiv_id == 433) {
        printf("INSERT : Leading promoters lists : ");
        for (auto it : indiv_->leading_prom_pos) {
          printf("%d (%d) || ", it.first, it.second);
        }
        printf("\n");

        printf("INSERT : Lagging promoters lists : ");
        for (auto it : indiv_->lagging_prom_pos) {
          printf("%d (%d) || ", it.first, it.second);
        }
        printf("\n");

        printf("INSERT : Leading promoters lists (promoters): ");
        for (auto it : indiv_->promoters) {
          printf("%d (%d) -- ", it.second->pos, it.first);
        }

        printf("\n");
          if (indiv_->promoters.size() != indiv_->leading_prom_pos.size() + indiv_->lagging_prom_pos.size()) {
              printf("Error unsynchronized cache !\n");
          }
      }*/

      indiv_->look_for_new_promoters_around(pos_3);
      indiv_->look_for_new_promoters_around(pos_3 + seg_length);
    }
  }


#ifdef WITH_BITSET
  delete duplicate_segment;
#else
  free(duplicate_segment);
#endif

    if (SIMD_Individual::standalone_simd) {
        Duplication *mut = new Duplication(pos_1, pos_2, pos_3, seg_length);

        indiv_->notifyObservers(MUTATION, mut);
        delete mut;
    }
  return true;
}

bool Dna_SIMD::do_translocation(int32_t pos_1, int32_t pos_2, int32_t pos_3,
                           int32_t pos_4, bool invert) {
  int32_t pos_min = Utils::min(pos_1,
                               Utils::min(pos_2, Utils::min(pos_3, pos_4)));

  if (not invert) {
    if (pos_min == pos_1) {
//      printf("A\n");
      ABCDE_to_ADCBE(pos_1, pos_3, pos_2, pos_4);
    }
    else if (pos_min == pos_2) {
//      printf("B\n");
      ABCDE_to_ADCBE(pos_2, pos_4, pos_1, pos_3);
    }
    else if (pos_min == pos_3) {
//      printf("C\n");
      ABCDE_to_ADCBE(pos_3, pos_2, pos_4, pos_1);
    }
    else { // if (pos_min == pos_4)
//      printf("D\n");
      ABCDE_to_ADCBE(pos_4, pos_1, pos_3, pos_2);
    }
  }
  else { // invert
    if (pos_min == pos_1) {
//      printf("E\n");
      ABCDE_to_ADBpCpE(pos_1, pos_3, pos_2, pos_4);
    }
    else if (pos_min == pos_2) {
//      printf("F\n");
      ABCDE_to_ADBpCpE(pos_2, pos_4, pos_1, pos_3);
    }
    else if (pos_min == pos_3) {
//      printf("G\n");
      ABCDE_to_ACpDpBE(pos_3, pos_2, pos_4, pos_1);
    }
    else { // if (pos_min == pos_4)
//      printf("H\n");
      ABCDE_to_ACpDpBE(pos_4, pos_1, pos_3, pos_2);
    }
  }

  int32_t segment_length = pos_2 - pos_1;

    if (SIMD_Individual::standalone_simd) {
        Translocation *mut = new Translocation(pos_1, pos_2, pos_3, pos_4,
                                               segment_length, invert);
        indiv_->notifyObservers(MUTATION, mut);
        delete mut;
    }

  return true;
}


/*
bool Dna_SIMD::do_translocation(int32_t pos_1, int32_t pos_2, int32_t pos_3,
                           int32_t pos_4, bool invert) {
  // Determine which position comes first and do the corresponding rearrangement
  // TODO(dpa) use min from std
  int32_t pos_min = Utils::min(pos_1,
                               Utils::min(pos_2, Utils::min(pos_3, pos_4)));

  if (not invert) {
    if (pos_min == pos_1) {
      ABCDE_to_ADCBE(pos_1, pos_3, pos_2, pos_4);
    }
    else if (pos_min == pos_2) {
      ABCDE_to_ADCBE(pos_2, pos_4, pos_1, pos_3);
    }
    else if (pos_min == pos_3) {
      ABCDE_to_ADCBE(pos_3, pos_2, pos_4, pos_1);
    }
    else { // if (pos_min == pos_4)
      ABCDE_to_ADCBE(pos_4, pos_1, pos_3, pos_2);
    }
  }
  else { // invert
    if (pos_min == pos_1) {
      ABCDE_to_ADBpCpE(pos_1, pos_3, pos_2, pos_4);
    }
    else if (pos_min == pos_2) {
      ABCDE_to_ADBpCpE(pos_2, pos_4, pos_1, pos_3);
    }
    else if (pos_min == pos_3) {
      ABCDE_to_ACpDpBE(pos_3, pos_2, pos_4, pos_1);
    }
    else { // if (pos_min == pos_4)
      ABCDE_to_ACpDpBE(pos_4, pos_1, pos_3, pos_2);
    }
  }

  return true;
}
*/

bool Dna_SIMD::do_inversion(int32_t pos_1, int32_t pos_2) {
// Invert segment going from pos_1 (included) to pos_2 (excluded)
// Exemple : sequence 011101001100 => 110011010001
  if (pos_1 == pos_2) return false; // Invert everything <=> Invert nothing!
  //
  //       pos_1         pos_2                   -> 0-
  //         |             |                   -       -
  // 0--------------------------------->      -         -
  //         ===============                  -         - pos_1
  //           tmp (copy)                      -       -
  //                                             -----      |
  //                                             pos_2    <-'
  //

  int32_t seg_length = pos_2 - pos_1;

  // Create the inverted sequence
#ifdef WITH_BITSET
  BitSet_SIMD* inverted_segment = bitset_->duplicate(pos_1,pos_2,true);

  inverted_segment->flip_all();
#else
  char* inverted_segment = NULL;
  posix_memalign((void **)&inverted_segment,64,(seg_length+1)* sizeof(char));//new char[seg_length + 1];


#ifdef __SIMD
#pragma omp simd
#endif
  for (int32_t i = 0, j = pos_2 - 1; i < seg_length; i++, j--) {
    if (data_[j] == '0') inverted_segment[i] = '1';
    else inverted_segment[i] = '0';
  }
  inverted_segment[seg_length] = '\0';
#endif
//  if (indiv_->indiv_id == 37) {
//    printf("DO_INVERSION 1 : Leading promoters lists : ");
//    for (auto it : indiv_->leading_prom_pos) {
//      printf("%d (%d) || ", it.first, it.second);
//    }
//    printf("\n");
//
//    printf("DO_INVERSION 1 : Lagging promoters lists : ");
//    for (auto it : indiv_->lagging_prom_pos) {
//      printf("%d (%d) || ", it.first, it.second);
//    }
//    printf("\n");
//
//    printf("DO_INVERSION 1 : Leading promoters lists (promoters): ");
//    for (auto it : indiv_->promoters) {
//      printf("%d (%d) -- ", it.second->pos, it.first);
//    }
//
//    printf("\n");
//  }
  // Remove promoters that included a breakpoint
  indiv_->remove_promoters_around(pos_1);

//  if (indiv_->indiv_id == 37) {
//    printf("DO_INVERSION 2.1 : Leading promoters lists : ");
//    for (auto it : indiv_->leading_prom_pos) {
//      printf("%d (%d) || ", it.first, it.second);
//    }
//    printf("\n");
//
//    printf("DO_INVERSION 2.1 : Lagging promoters lists : ");
//    for (auto it : indiv_->lagging_prom_pos) {
//      printf("%d (%d) || ", it.first, it.second);
//    }
//    printf("\n");
//
//    printf("DO_INVERSION 2.1 : Leading promoters lists (promoters): ");
//    for (auto it : indiv_->promoters) {
//      printf("%d (%d) -- ", it.second->pos, it.first);
//    }
//
//    printf("\n");
//  }

  indiv_->remove_promoters_around(pos_2);

//  if (indiv_->indiv_id == 37) {
//    printf("DO_INVERSION 2 : Leading promoters lists : ");
//    for (auto it : indiv_->leading_prom_pos) {
//      printf("%d (%d) || ", it.first, it.second);
//    }
//    printf("\n");
//
//    printf("DO_INVERSION 2 : Lagging promoters lists : ");
//    for (auto it : indiv_->lagging_prom_pos) {
//      printf("%d (%d) || ", it.first, it.second);
//    }
//    printf("\n");
//
//    printf("DO_INVERSION 2 : Leading promoters lists (promoters): ");
//    for (auto it : indiv_->promoters) {
//      printf("%d (%d) -- ", it.second->pos, it.first);
//    }
//
//    printf("\n");
//  }
  // Invert the sequence
#ifdef WITH_BITSET
  bitset_->replace(inverted_segment,pos_1);
#else
  replace(pos_1, inverted_segment, seg_length);
#endif

  // Update promoter list
  if (length() >= PROM_SIZE) {
    indiv_->invert_promoters_included_in(pos_1, pos_2);


//    if (indiv_->indiv_id == 37) {
//      printf("DO_INVERSION 3 : Leading promoters lists : ");
//      for (auto it : indiv_->leading_prom_pos) {
//        printf("%d (%d) || ", it.first, it.second);
//      }
//      printf("\n");
//
//      printf("DO_INVERSION 3 : Lagging promoters lists : ");
//      for (auto it : indiv_->lagging_prom_pos) {
//        printf("%d (%d) || ", it.first, it.second);
//      }
//      printf("\n");
//
//      printf("DO_INVERSION 3 : Leading promoters lists (promoters): ");
//      for (auto it : indiv_->promoters) {
//        printf("%d (%d) -- ", it.second->pos, it.first);
//      }
//
//      printf("\n");
//    }

    indiv_->look_for_new_promoters_around(pos_1);
    indiv_->look_for_new_promoters_around(pos_2);

//    if (indiv_->indiv_id == 37) {
//      printf("DO_INVERSION 4 : Leading promoters lists : ");
//      for (auto it : indiv_->leading_prom_pos) {
//        printf("%d (%d) || ", it.first, it.second);
//      }
//      printf("\n");
//
//      printf("DO_INVERSION 4 : Lagging promoters lists : ");
//      for (auto it : indiv_->lagging_prom_pos) {
//        printf("%d (%d) || ", it.first, it.second);
//      }
//      printf("\n");
//
//      printf("DO_INVERSION 4 : Leading promoters lists (promoters): ");
//      for (auto it : indiv_->promoters) {
//        printf("%d (%d) -- ", it.second->pos, it.first);
//      }
//
//      printf("\n");
//    }
  }

#ifdef WITH_BITSET
  delete inverted_segment;
#else
  free(inverted_segment);
#endif

  if (SIMD_Individual::standalone_simd) {
      Inversion *mut = new Inversion(pos_1, pos_2, seg_length);
      indiv_->notifyObservers(MUTATION, mut);
      delete mut;
  }

  return true;
}


bool Dna_SIMD::do_deletion(int32_t pos_1, int32_t pos_2) {
    //printf("DO DELETION is %d %d -- %d\n",pos_1,pos_2,length());
  // Delete segment going from pos_1 (included) to pos_2 (excluded)
  if (pos_1 < pos_2) {
    //
    //       pos_1         pos_2                   -> 0-
    //         |             |                   -       -
    // 0--------------------------------->      -         -
    //         ===============                  -         - pos_1
    //           tmp (copy)                      -       -
    //                                             -----      |
    //                                             pos_2    <-'
    //

    int32_t segment_length = pos_2 - pos_1;

    /*if (indiv_->indiv_id == 128) {
      printf("DO_DELETION 1 : Leading promoters lists : ");
      for (auto it : indiv_->leading_prom_pos) {
        printf("%d (%d) || ", it.first, it.second);
      }
      printf("\n");

      printf("DO_DELETION 1 : Lagging promoters lists : ");
      for (auto it : indiv_->lagging_prom_pos) {
        printf("%d (%d) || ", it.first, it.second);
      }
      printf("\n");

      printf("DO_DELETION 1 : Leading promoters lists (promoters): ");
      for (auto it : indiv_->promoters) {
        printf("%d (%d) -- ", it.second->pos, it.first);
      }

      printf("\n");
    }*/
    // Remove promoters containing at least one nucleotide from the sequence
    // to delete
    indiv_->remove_promoters_around(pos_1, pos_2);

//    if (indiv_->indiv_id == 30) {
//      printf("DO_DELETION 2 : Leading promoters lists : ");
//      for (auto it : indiv_->leading_prom_pos) {
//        printf("%d (%d) || ", it.first, it.second);
//      }
//      printf("\n");
//
//      printf("DO_DELETION 2 : Lagging promoters lists : ");
//      for (auto it : indiv_->lagging_prom_pos) {
//        printf("%d (%d) || ", it.first, it.second);
//      }
//      printf("\n");
//
//      printf("DO_DELETION 2 : Leading promoters lists (promoters): ");
//      for (auto it : indiv_->promoters) {
//        printf("%d (%d) -- ", it.second->pos, it.first);
//      }
//
//      printf("\n");
//    }

    // Delete the sequence between pos_1 and pos_2
#ifdef WITH_BITSET
    bitset_->remove(pos_1,pos_2);
#else
    remove(pos_1, pos_2);
#endif


    // Update promoter list
    if (length() >= PROM_SIZE) {
      indiv_->move_all_promoters_after(pos_1, -segment_length);

      indiv_->look_for_new_promoters_around(pos_1);
    }
  }
  else { // if (pos_1 >= pos_2)
    // The segment to delete includes the origin of replication.
    // The deletion process will be done in two steps.
    //
    //                                            ,->
    //    pos_2                 pos_1            |      -> 0-
    //      |                     |                   -       - pos_2
    // 0--------------------------------->     pos_1 -         -
    // =====                      =======            -         -
    //  tmp2                        tmp1              -       -
    //                                                  -----
    //
    //

    // int32_t segment_length = length_ + pos_2 - pos_1; //useless variable

    // Remove promoters containing at least one nucleotide from the sequence
    // to delete
/*    if (indiv_->indiv_id == 626) {
      printf("DO_DELETION 1 : Leading promoters lists : ");
      for (auto it : indiv_->leading_prom_pos) {
        printf("%d (%d) || ", it.first, it.second);
      }
      printf("\n");

      printf("DO_DELETION 1 : Lagging promoters lists : ");
      for (auto it : indiv_->lagging_prom_pos) {
        printf("%d (%d) || ", it.first, it.second);
      }
      printf("\n");

      printf("DO_DELETION 1 : Leading promoters lists (promoters): ");
      for (auto it : indiv_->promoters) {
        printf("%d (%d) -- ", it.second->pos, it.first);
      }

      printf("\n");
    }*/
    indiv_->remove_promoters_around(pos_1, pos_2);

/*    if (indiv_->indiv_id == 626) {
      printf("DO_DELETION 2 : Leading promoters lists : ");
      for (auto it : indiv_->leading_prom_pos) {
        printf("%d (%d) || ", it.first, it.second);
      }
      printf("\n");

      printf("DO_DELETION 2 : Lagging promoters lists : ");
      for (auto it : indiv_->lagging_prom_pos) {
        printf("%d (%d) || ", it.first, it.second);
      }
      printf("\n");

      printf("DO_DELETION 2 : Leading promoters lists (promoters): ");
      for (auto it : indiv_->promoters) {
        printf("%d (%d) -- ", it.second->pos, it.first);
      }

      printf("\n");
    }*/

    // Delete the sequence between pos_1 and pos_2
#ifdef WITH_BITSET
    bitset_->remove(pos_1,bitset_->length_);
    bitset_->remove(0,pos_2);
#else
    remove(pos_1, length_); // delete tmp1 from genome
    remove(0, pos_2);       // delete tmp2 from genome
#endif

    // Update promoter list
    if (length() >= PROM_SIZE) {
/*      if (indiv_->indiv_id == 626) {
        printf("DO_DELETION 3 : Leading promoters lists : ");
        for (auto it : indiv_->leading_prom_pos) {
          printf("%d (%d) || ", it.first, it.second);
        }
        printf("\n");

        printf("DO_DELETION 3 : Lagging promoters lists : ");
        for (auto it : indiv_->lagging_prom_pos) {
          printf("%d (%d) || ", it.first, it.second);
        }
        printf("\n");

        printf("DO_DELETION 3 : Leading promoters lists (promoters): ");
        for (auto it : indiv_->promoters) {
          printf("%d (%d) -- ", it.second->pos, it.first);
        }

        printf("\n");
      }*/
      indiv_->move_all_promoters_after(0, -pos_2);

/*      if (indiv_->indiv_id == 626) {
        printf("DO_DELETION 4 : Leading promoters lists : ");
        for (auto it : indiv_->leading_prom_pos) {
          printf("%d (%d) || ", it.first, it.second);
        }
        printf("\n");

        printf("DO_DELETION 4 : Lagging promoters lists : ");
        for (auto it : indiv_->lagging_prom_pos) {
          printf("%d (%d) || ", it.first, it.second);
        }
        printf("\n");

        printf("DO_DELETION 4 : Leading promoters lists (promoters): ");
        for (auto it : indiv_->promoters) {
          printf("%d (%d) -- ", it.second->pos, it.first);
        }

        printf("\n");
      }*/
      indiv_->look_for_new_promoters_around(0);
    }
  }
    if (SIMD_Individual::standalone_simd) {
        int32_t segment_length = Utils::mod(pos_2 - pos_1 - 1, length()) + 1;
        Deletion *mut = new Deletion(pos_1, pos_2, segment_length);
        indiv_->notifyObservers(MUTATION, mut);
        delete mut;
    }

  return true;
}

void Dna_SIMD::apply_mutations_standalone() {
  MutationEvent* repl = indiv_->exp_m_->
      dna_mutator_array_[indiv_->indiv_id]->generate_next_mutation(length());

  while (indiv_->exp_m_->dna_mutator_array_[indiv_->indiv_id]->mutation_available()) {
    if (repl != nullptr)
      switch (repl->type()) {
        case DO_SWITCH:
//            if (indiv_->indiv_id == 433) printf("Start switch at %d\n",repl->pos_1());
          do_switch(repl->pos_1());
          nb_swi_++;
          nb_mut_++;
//        printf("End switch at %d\n",repl->pos_1());
          break;
        case SMALL_INSERTION:
//            if (indiv_->indiv_id == 433)         printf("Start insertion at %d (%d %s)\n",repl->pos_1(),repl->number(),repl->seq());
#ifdef WITH_BITSET
          do_small_insertion(repl->pos_1(), repl->seq());
#else
          do_small_insertion(repl->pos_1(), repl->number(), repl->seq());
#endif
          nb_indels_++;
          nb_mut_++;
//        printf("End insertion at %d (%d)\n",repl->pos_1(),repl->number(),repl->seq());
          break;
        case SMALL_DELETION:
          //if (indiv_->indiv_id == 626 && AeTime::time() == 21)
//              if (indiv_->indiv_id == 433) printf("Start deletion at %d (%d)\n",repl->pos_1(),repl->number());
          do_small_deletion(repl->pos_1(), repl->number());
          nb_indels_++;
          nb_mut_++;
//        printf("End deletion at %d (%d)\n",repl->pos_1(),repl->number());
          break;
        case DUPLICATION:
//            if (indiv_->indiv_id == 433)         printf("Start duplication at %d (%d %d)\n",repl->pos_1(),repl->pos_2(),repl->pos_3());
          do_duplication(repl->pos_1(), repl->pos_2(), repl->pos_3());
          nb_large_dupl_++;
          nb_rear_++;
//        printf("End duplication at %d (%d %d)\n",repl->pos_1(),repl->pos_2(),repl->pos_3());
          break;
        case TRANSLOCATION:
//            if (indiv_->indiv_id == 433)         printf("Start translocation at %d (%d %d %d %d)\n",repl->pos_1(),repl->pos_2(),repl->pos_3(),repl->pos_4(),repl->invert());
          do_translocation(repl->pos_1(), repl->pos_2(), repl->pos_3(),
                           repl->pos_4(), repl->invert());
          nb_large_trans_++;
          nb_rear_++;
//        printf("End translocation at %d (%d %d %d %d)\n",repl->pos_1(),repl->pos_2(),repl->pos_3(),repl->pos_4(),repl->invert());
          break;
        case INVERSION:
//            if (indiv_->indiv_id == 433)         printf("Start invertion at %d (%d)\n",repl->pos_1(),repl->pos_2());
          do_inversion(repl->pos_1(), repl->pos_2());
          nb_large_inv_++;
          nb_rear_++;
//        printf("End invertion at %d (%d)\n",repl->pos_1(),repl->pos_2());
          break;
        case DELETION:
            //if (indiv_->indiv_id == 626 && AeTime::time() == 21)
//                if (indiv_->indiv_id == 433) printf("Start LARGE deletion at %d (%d)\n",repl->pos_1(),repl->pos_2());
          do_deletion(repl->pos_1(), repl->pos_2());
          nb_large_del_++;
          nb_rear_++;
//        printf("End LARGE deletion at %d (%d)\n",repl->pos_1(),repl->pos_2());
          break;
      }

    repl = indiv_->exp_m_->
        dna_mutator_array_[indiv_->indiv_id]->generate_next_mutation(length());
  }
}


void Dna_SIMD::apply_mutations() {
//
//  if (indiv_->indiv_id == 30) {
//    printf("APPLY_MUTATION : Leading promoters lists : ");
//    for (auto it : indiv_->leading_prom_pos) {
//      printf("%d (%d) || ", it.first, it.second);
//    }
//    printf("\n");
//
//    printf("APPLY_MUTATION : Leading promoters lists (promoters): ");
//    for (auto it : indiv_->promoters) {
//      printf("%d (%d) -- ", it.second->pos, it.first);
//    }
//
//    printf("\n");
//  }


  for (auto repl : indiv_->exp_m_->dna_mutator_array_[indiv_->indiv_id]->mutation_list_) {
/*    if (indiv_->indiv_id == 49) {
      printf("Mutation type %d\n",repl->type());
    }*/

    switch(repl->type()) {
      case DO_SWITCH:
        //printf("Start switch at %d\n",repl->pos_1());
        do_switch(repl->pos_1());
//        printf("End switch at %d\n",repl->pos_1());
        break;
      case SMALL_INSERTION:
//        printf("Start insertion at %d (%d %s)\n",repl->pos_1(),repl->number(),repl->seq());
#ifdef WITH_BITSET
        do_small_insertion(repl->pos_1(),repl->seq());
#else
        do_small_insertion(repl->pos_1(),repl->number(),repl->seq());
#endif
//        printf("End insertion at %d (%d)\n",repl->pos_1(),repl->number(),repl->seq());
        break;
      case SMALL_DELETION:
//        printf("Start deletion at %d (%d)\n",repl->pos_1(),repl->number());
        do_small_deletion(repl->pos_1(),repl->number());
//        printf("End deletion at %d (%d)\n",repl->pos_1(),repl->number());
        break;
      case DUPLICATION:
        /*if (indiv_->indiv_id == 49) printf("Start duplication at %d (%d %d) size %d\n",repl->pos_1(),repl->pos_2(),repl->pos_3(),length_);*/
        do_duplication(repl->pos_1(),repl->pos_2(),repl->pos_3());
//        printf("End duplication at %d (%d %d)\n",repl->pos_1(),repl->pos_2(),repl->pos_3());
        break;
      case TRANSLOCATION:
//        printf("Start translocation at %d (%d %d %d %d)\n",repl->pos_1(),repl->pos_2(),repl->pos_3(),repl->pos_4(),repl->invert());
        do_translocation(repl->pos_1(), repl->pos_2(), repl->pos_3(),
                         repl->pos_4(), repl->invert());
//        printf("End translocation at %d (%d %d %d %d)\n",repl->pos_1(),repl->pos_2(),repl->pos_3(),repl->pos_4(),repl->invert());
        break;
      case INVERSION:
//        printf("Start invertion at %d (%d)\n",repl->pos_1(),repl->pos_2());
        do_inversion(repl->pos_1(), repl->pos_2());
//        printf("End invertion at %d (%d)\n",repl->pos_1(),repl->pos_2());
        break;
      case DELETION:
//        printf("Start LARGE deletion at %d (%d)\n",repl->pos_1(),repl->pos_2());
        do_deletion(repl->pos_1(), repl->pos_2());
//        printf("End LARGE deletion at %d (%d)\n",repl->pos_1(),repl->pos_2());
        break;
    }
//
//    if (indiv_->indiv_id == 30) {
//      printf("APPLY_MUTATION : Leading promoters lists : ");
//      for (auto it : indiv_->leading_prom_pos) {
//        printf("%d (%d) || ", it.first, it.second);
//      }
//      printf("\n");
//
//      printf("APPLY_MUTATION : Leading promoters lists (promoters): ");
//      for (auto it : indiv_->promoters) {
//        printf("%d (%d) -- ", it.second->pos, it.first);
//      }
//
//      printf("\n");
//    }

  }

}


void Dna_SIMD::ABCDE_to_ADCBE(int32_t pos_B, int32_t pos_C, int32_t pos_D,
                         int32_t pos_E) {
  // Rearrange the sequence from ABCDE to ADCBE (complex translocation
  // of segment defined between positions pos_B and pos_D)
  //
  // Segments are identified by pos_x values as shown below.
  //
  // TODO(dpa) CHECK THIS !!!
  // WARNING : Segment C includes nucleotide at pos_D // NOTE : WTF???
  //
  //         A      B        C       D       E
  //      |----->=======[>=======>-------[>-------|        =>
  //          pos_B   pos_C    pos_D   pos_E
  //
  //                         |
  //                         V
  //
  //         A      D        C       B        E
  //      |----->-------[>=======>=======[>-------|
  // Compute segment lengths
  int32_t len_A = pos_B;
  int32_t len_B = pos_C - pos_B;
  int32_t len_C = pos_D - pos_C;
  int32_t len_D = pos_E - pos_D;
  int32_t len_E = length() - pos_E;
  int32_t len_AD = len_A + len_D;
  int32_t len_ADC = len_AD + len_C;
  int32_t len_ADCB = len_ADC + len_B;

  // Create new sequence
#ifdef WITH_BITSET
  BitSet_SIMD* bitset = new BitSet_SIMD(bitset_->length_);
  bitset->replace(bitset_,0,len_A,0);
  bitset->replace(bitset_,pos_D,pos_D+len_D,len_A);
  bitset->replace(bitset_,pos_C,pos_C+len_C,len_AD);
  bitset->replace(bitset_,pos_B,pos_B+len_B,len_ADC);
  bitset->replace(bitset_,pos_E,pos_E+len_E,len_ADCB);

  delete bitset_;
  bitset_ = bitset;
#else
  char* new_genome;
  posix_memalign((void **)&new_genome,64,nb_blocks_ * BLOCK_SIZE* sizeof(char));//new char[nb_blocks_ * BLOCK_SIZE];

  memset(new_genome,0,(length_) * sizeof(char));

  memcpy(new_genome, data_, len_A * sizeof(char));
  memcpy(&new_genome[len_A], &data_[pos_D], len_D * sizeof(char));
  memcpy(&new_genome[len_AD], &data_[pos_C], len_C * sizeof(char));
  memcpy(&new_genome[len_ADC], &data_[pos_B], len_B * sizeof(char));
  memcpy(&new_genome[len_ADCB], &data_[pos_E], len_E * sizeof(char));
  new_genome[length_] = '\0';

  // Replace sequence
  // NB : The size of the genome doesn't change. Therefore, we don't nee
  // to update length_ and nb_blocks_
  free(data_);
  data_ = new_genome;
#endif


  // ========== Update promoter list ==========
  if (length() >= PROM_SIZE) {
    //printf("WHO AMI ??? %d\n",indiv_->indiv_id);
//    if (indiv_->indiv_id == 152) {
//      printf("BEFORE_REMOVE : Leading promoters lists : ");
//      for (auto it : indiv_->leading_prom_pos) {
//        printf("%d (%d) || ", it.first, it.second);
//      }
//      printf("\n");
//
//      printf("BEFORE_REMOVE : Lagging promoters lists : ");
//      for (auto it : indiv_->lagging_prom_pos) {
//        printf("%d (%d) || ", it.first, it.second);
//      }
//      printf("\n");
//
//      printf("BEFORE_REMOVE : Promoters lists (promoters): ");
//      for (auto it : indiv_->promoters) {
//        printf("%d (%d) -- ", it.second->pos, it.first);
//      }
//
//      printf("\n");
//    }
    // Remove promoters that include a breakpoint
    indiv_->remove_promoters_around(pos_B);
    indiv_->remove_promoters_around(pos_C);
    indiv_->remove_promoters_around(pos_D);
    indiv_->remove_promoters_around(pos_E);
//    if (indiv_->indiv_id == 152) {
//      printf("AFTER_REMOVE : Leading promoters lists : ");
//      for (auto it : indiv_->leading_prom_pos) {
//        printf("%d (%d) || ", it.first, it.second);
//      }
//      printf("\n");
//
//      printf("AFTER_REMOVE : Lagging promoters lists : ");
//      for (auto it : indiv_->lagging_prom_pos) {
//        printf("%d (%d) || ", it.first, it.second);
//      }
//      printf("\n");
//
//      printf("AFTER_REMOVE : Promoters lists (promoters): ");
//      for (auto it : indiv_->promoters) {
//        printf("%d (%d) -- ", it.second->pos, it.first);
//      }
//
//      printf("\n");
//    }
    // Create temporary lists for promoters to move and/or invert
    std::vector<std::list<promoterStruct*>> promoters_B = {{},
                                     {}};
    std::vector<std::list<promoterStruct*>> promoters_C = {{},
                                     {}};
    std::vector<std::list<promoterStruct*>> promoters_D = {{},
                                     {}};
    // Extract promoters that are totally included in each segment to be moved
    // and shift them to their new positions
    if (len_B >= PROM_SIZE) {
      indiv_->extract_promoters_included_in(pos_B, pos_C, promoters_B);

      Internal_SIMD_Struct::shift_promoters(promoters_B, len_D + len_C,
                                   length());
    }
    if (len_C >= PROM_SIZE) {
      indiv_->extract_promoters_included_in(pos_C, pos_D, promoters_C);

      Internal_SIMD_Struct::shift_promoters(promoters_C, len_D - len_B,
                                   length());
    }
    if (len_D >= PROM_SIZE) {
      indiv_->extract_promoters_included_in(pos_D, pos_E, promoters_D);

      Internal_SIMD_Struct::shift_promoters(promoters_D, -len_B - len_C,
                                   length());
    }
//    if (indiv_->indiv_id == 152) {
//      printf("AFTER_EXTRACT : Leading promoters lists : ");
//      for (auto it : indiv_->leading_prom_pos) {
//        printf("%d (%d) || ", it.first, it.second);
//      }
//      printf("\n");
//
//      printf("AFTER_EXTRACT : Lagging promoters lists : ");
//      for (auto it : indiv_->lagging_prom_pos) {
//        printf("%d (%d) || ", it.first, it.second);
//      }
//      printf("\n");
//
//      printf("AFTER_EXTRACT : Promoters lists (promoters): ");
//      for (auto it : indiv_->promoters) {
//        printf("%d (%d) -- ", it.second->pos, it.first);
//      }
//
//      printf("\n");
//    }
    // Reinsert the shifted promoters
    indiv_->insert_promoters(promoters_B);

    indiv_->insert_promoters(promoters_C);

    indiv_->insert_promoters(promoters_D);
//    if (indiv_->indiv_id == 152) {
//      printf("AFTER_INSERT : Leading promoters lists : ");
//      for (auto it : indiv_->leading_prom_pos) {
//        printf("%d (%d) || ", it.first, it.second);
//      }
//      printf("\n");
//
//      printf("AFTER_INSERT : Lagging promoters lists : ");
//      for (auto it : indiv_->lagging_prom_pos) {
//        printf("%d (%d) || ", it.first, it.second);
//      }
//      printf("\n");
//
//      printf("AFTER_INSERT : Promoters lists (promoters): ");
//      for (auto it : indiv_->promoters) {
//        printf("%d (%d) -- ", it.second->pos, it.first);
//      }
//
//      printf("\n");
//    }
    // 5) Look for new promoters including a breakpoint
    indiv_->look_for_new_promoters_around(len_A);
    indiv_->look_for_new_promoters_around(len_AD);
    indiv_->look_for_new_promoters_around(len_ADC);
    indiv_->look_for_new_promoters_around(len_ADCB);
//    if (indiv_->indiv_id == 152) {
//      printf("AFTER_LOOK : Leading promoters lists : ");
//      for (auto it : indiv_->leading_prom_pos) {
//        printf("%d (%d) || ", it.first, it.second);
//      }
//      printf("\n");
//
//      printf("AFTER_LOOK : Lagging promoters lists : ");
//      for (auto it : indiv_->lagging_prom_pos) {
//        printf("%d (%d) || ", it.first, it.second);
//      }
//      printf("\n");
//
//      printf("AFTER_LOOK : Promoters lists (promoters): ");
//      for (auto it : indiv_->promoters) {
//        printf("%d (%d) -- ", it.second->pos, it.first);
//      }
//
//      printf("\n");
//    }
  }
}

void Dna_SIMD::ABCDE_to_ADBpCpE(int32_t pos_B, int32_t pos_C, int32_t pos_D,
                           int32_t pos_E) {
  // Rearrange the sequence from ABCDE to ADBpCpE (complex translocation
  // with inversion of segment defined between positions pos_B and pos_D)
  // Bp (resp Cp) stands for inverted B (resp C)
  //
  // Segments are identified by pos_x values as shown below.
  //
  // TODO(dpa) CHECK THIS !!!
  // WARNING : Segment C includes nucleotide at pos_D // NOTE : WTF???
  //
  //         A      B        C       D        E
  //      |----->=======[>=======>-------<]-------|
  //          pos_B   pos_C    pos_D   pos_E
  //
  //                         |
  //                         V
  //
  //         A      D        Bp      Cp       E
  //      |----->-------<]=======<=======<]-------|
  // Compute segment lengths
  int32_t len_A = pos_B;
  int32_t len_B = pos_C - pos_B;
  int32_t len_C = pos_D - pos_C;
  int32_t len_D = pos_E - pos_D;
  int32_t len_E = length() - pos_E;
  int32_t len_AD = len_A + len_D;
  int32_t len_ADB = len_AD + len_B;
  int32_t len_ADBC = len_ADB + len_C;

#ifdef WITH_BITSET
  // Create new sequence
  BitSet_SIMD* bitset = new BitSet_SIMD(bitset_->length_);

  // Copy segments A and D
  bitset->replace(bitset_,0,len_A,0);
  bitset->replace(bitset_,pos_D,pos_D+len_D,len_A);

  // Build Bp and put it in the new genome
  BitSet_SIMD* inverted_segment = bitset_->duplicate(pos_C-len_B,pos_C,true);
  inverted_segment->flip_all();

  bitset->replace(inverted_segment,0,len_B,len_AD);

  delete inverted_segment;


  // Build Cp and put it in the new genome
  inverted_segment = bitset_->duplicate(pos_D-len_C,pos_D,true);
  inverted_segment->flip_all();

  bitset->replace(inverted_segment,0,len_C,len_ADB);

  delete inverted_segment;

  // Copy segment E into the new genome
  bitset->replace(bitset_,pos_E,pos_E+len_E,len_ADBC);

  // Replace sequence
  delete bitset_;
  bitset_ = bitset;
#else
  // Create new sequence
  char* new_genome;
  posix_memalign((void **)&new_genome,64,nb_blocks_ * BLOCK_SIZE* sizeof(char));//new char[nb_blocks_ * BLOCK_SIZE];
  memset(new_genome,0,(length_) * sizeof(char));

  // Copy segments A and D
  memcpy(new_genome, data_, len_A * sizeof(char));
  memcpy(&new_genome[len_A], &data_[pos_D], len_D * sizeof(char));


  // Build Bp and put it in the new genome
  char* inverted_segment;
  posix_memalign((void **)&inverted_segment,64,(len_B+1)* sizeof(char));//new char[len_B + 1];

//#pragma simd
//#pragma distribute_point
#ifdef __SIMD
#pragma omp simd
#endif
  for (int32_t i = 0, j = pos_C - 1; i < len_B; i++, j--) {
    if (data_[j] == '0') inverted_segment[i] = '1';
    else inverted_segment[i] = '0';
  }
  inverted_segment[len_B] = '\0';

  memcpy(&new_genome[len_AD], inverted_segment, len_B * sizeof(char));

  free(inverted_segment);


  // Build Cp and put it in the new genome
  posix_memalign((void **)&inverted_segment,64,(len_C+1)* sizeof(char));//new char[len_C + 1];

//#pragma simd
//#pragma distribute_point
#ifdef __SIMD
#pragma omp simd
#endif
  for (int32_t i = 0, j = pos_D - 1; i < len_C; i++, j--) {
    if (data_[j] == '0') inverted_segment[i] = '1';
    else inverted_segment[i] = '0';
  }
  inverted_segment[len_C] = '\0';

  memcpy(&new_genome[len_ADB], inverted_segment, len_C * sizeof(char));

  free(inverted_segment);

  // Copy segment E into the new genome
  memcpy(&new_genome[len_ADBC], &data_[pos_E], len_E * sizeof(char));
  new_genome[length_] = '\0';


  // Replace sequence
  free(data_);
  data_ = new_genome;
#endif

  // ========== Update promoter list ==========
  if (length() >= PROM_SIZE) {
    /*if (indiv_->indiv_id == 26) {
      printf("BEFORE_REMOVE : Leading promoters lists : ");
      for (auto it : indiv_->leading_prom_pos) {
        printf("%d (%d) || ", it.first, it.second);
      }
      printf("\n");

      printf("BEFORE_REMOVE : Lagging promoters lists : ");
      for (auto it : indiv_->lagging_prom_pos) {
        printf("%d (%d) || ", it.first, it.second);
      }
      printf("\n");

      printf("BEFORE_REMOVE : Promoters lists (promoters): ");
      for (auto it : indiv_->promoters) {
        printf("%d (%d) -- ", it.second->pos, it.first);
      }

      printf("\n");
    }*/
    // Remove promoters that include a breakpoint
    indiv_->remove_promoters_around(pos_B);
    indiv_->remove_promoters_around(pos_C);
    indiv_->remove_promoters_around(pos_D);
    indiv_->remove_promoters_around(pos_E);

    /*if (indiv_->indiv_id == 26) {
      printf("AFTER_REMOVE : Leading promoters lists : ");
      for (auto it : indiv_->leading_prom_pos) {
        printf("%d (%d) || ", it.first, it.second);
      }
      printf("\n");

      printf("AFTER_REMOVE : Lagging promoters lists : ");
      for (auto it : indiv_->lagging_prom_pos) {
        printf("%d (%d) || ", it.first, it.second);
      }
      printf("\n");

      printf("AFTER_REMOVE : Promoters lists (promoters): ");
      for (auto it : indiv_->promoters) {
        printf("%d (%d) -- ", it.second->pos, it.first);
      }

      printf("\n");
    }*/

    // Create temporary lists for promoters to move and/or invert
    std::vector<std::list<promoterStruct*>> promoters_B = {{},
                                     {}};
    std::vector<std::list<promoterStruct*>> promoters_C = {{},
                                     {}};
    std::vector<std::list<promoterStruct*>> promoters_D = {{},
                                     {}};

    // 2) Extract promoters that are totally included in each segment to be
    //    moved (B, C and D)
    if (len_B >= PROM_SIZE) {
      indiv_->extract_promoters_included_in(pos_B, pos_C, promoters_B);
    }
    if (len_C >= PROM_SIZE) {
      indiv_->extract_promoters_included_in(pos_C, pos_D, promoters_C);
    }
    if (len_D >= PROM_SIZE) {
      indiv_->extract_promoters_included_in(pos_D, pos_E, promoters_D);
    }

    /*if (indiv_->indiv_id == 26) {
      printf("AFTER_EXTRACT : Leading promoters lists : ");
      for (auto it : indiv_->leading_prom_pos) {
        printf("%d (%d) || ", it.first, it.second);
      }
      printf("\n");

      printf("AFTER_EXTRACT : Lagging promoters lists : ");
      for (auto it : indiv_->lagging_prom_pos) {
        printf("%d (%d) || ", it.first, it.second);
      }
      printf("\n");

      printf("AFTER_EXTRACT : Promoters lists (promoters): ");
      for (auto it : indiv_->promoters) {
        printf("%d (%d) -- ", it.second->pos, it.first);
      }

      printf("\n");

      printf("PROMOTER_LIST BEFORE_INVERT B: ");
      for (auto strand: {LEADING, LAGGING}) {
        for (auto it : promoters_B[strand]) {
          printf("%d || ", it->pos);
        }
      }
      printf("\n");

      printf("PROMOTER_LIST BEFORE_INVERT C: ");
      for (auto strand: {LEADING, LAGGING}) {
        for (auto it : promoters_C[strand]) {
          printf("%d || ", it->pos);
        }
      }
      printf("\n");

      printf("PROMOTER_LIST BEFORE_INVERT D: ");
      for (auto strand: {LEADING, LAGGING}) {
        for (auto it : promoters_D[strand]) {
          printf("%d || ", it->pos);
        }
      }
      printf("\n");
    }*/


    // 3a) Invert promoters of segments B and C
    Internal_SIMD_Struct::invert_promoters(promoters_B, pos_B, pos_C);

    Internal_SIMD_Struct::invert_promoters(promoters_C, pos_C, pos_D);

    /*if (indiv_->indiv_id == 12) {
      printf("PROMOTER_LIST AFTER_INVERT B: ");
      for (auto strand: {LEADING, LAGGING}) {
        for (auto it : promoters_B[strand]) {
          printf("%d || ", it->pos);
        }
      }
      printf("\n");

      printf("PROMOTER_LIST AFTER_INVERT C: ");
      for (auto strand: {LEADING, LAGGING}) {
        for (auto it : promoters_C[strand]) {
          printf("%d || ", it->pos);
        }
      }
      printf("\n");
    }*/

    // 3b) Shift these promoters positions
    Internal_SIMD_Struct::shift_promoters(promoters_B, len_D,
                                 length());

    Internal_SIMD_Struct::shift_promoters(promoters_C, len_D,
                                 length());

    Internal_SIMD_Struct::shift_promoters(promoters_D, -len_B - len_C,
                                 length());

    /*if (indiv_->indiv_id == 12) {
      printf("PROMOTER_LIST AFTER_SHIFT B: ");
      for (auto strand: {LEADING, LAGGING}) {
        for (auto it : promoters_B[strand]) {
          printf("%d || ", it->pos);
        }
      }
      printf("\n");

      printf("PROMOTER_LIST AFTER_SHIFT C: ");
      for (auto strand: {LEADING, LAGGING}) {
        for (auto it : promoters_C[strand]) {
          printf("%d || ", it->pos);
        }
      }
      printf("\n");

      printf("PROMOTER_LIST AFTER_SHIFT D: ");
      for (auto strand: {LEADING, LAGGING}) {
        for (auto it : promoters_D[strand]) {
          printf("%d || ", it->pos);
        }
      }
      printf("\n");
    }*/

    // 4) Reinsert the shifted promoters
    indiv_->insert_promoters(promoters_C);

    indiv_->insert_promoters(promoters_B);

    indiv_->insert_promoters(promoters_D);

    /*if (indiv_->indiv_id == 12) {
      printf("AFTER_INSERT : Leading promoters lists : ");
      for (auto it : indiv_->leading_prom_pos) {
        printf("%d (%d) || ", it.first, it.second);
      }
      printf("\n");

      printf("AFTER_INSERT : Lagging promoters lists : ");
      for (auto it : indiv_->lagging_prom_pos) {
        printf("%d (%d) || ", it.first, it.second);
      }
      printf("\n");

      printf("AFTER_INSERT : Promoters lists (promoters): ");
      for (auto it : indiv_->promoters) {
        printf("%d (%d) -- ", it.second->pos, it.first);
      }

      printf("\n");
    }*/

    // 5) Look for new promoters including a breakpoint
    indiv_->look_for_new_promoters_around(len_A);
    indiv_->look_for_new_promoters_around(len_AD);
    indiv_->look_for_new_promoters_around(len_ADB);
    indiv_->look_for_new_promoters_around(len_ADBC);
  }
}

void Dna_SIMD::ABCDE_to_ACpDpBE(int32_t pos_B, int32_t pos_C, int32_t pos_D,
                           int32_t pos_E) {
  // Rearrange the sequence from ABCDE to ACpDpBE (complex translocation with
  // inversion of segment defined between positions pos_C and pos_E)
  // Cp (resp Dp) stands for inverted C (resp D)
  //
  // Segments are identified by pos_x values as shown below.
  //
  // TODO(dpa) CHECK THIS !!!
  // WARNING : Segment D includes nucleotide at pos_E // NOTE : WTF???
  //
  //         A      B        C       D       E
  //      |----<]-------->=======[>=======>-------|
  //          pos_B    pos_C    pos_D   pos_E
  //
  //                         |
  //                         V
  //
  //          A       C'      D'       B       E
  //       |-----<]=======>=======<]------->-------|



  // Compute segment lengths
  int32_t len_A = pos_B;
  int32_t len_B = pos_C - pos_B;
  int32_t len_C = pos_D - pos_C;
  int32_t len_D = pos_E - pos_D;
  int32_t len_E = length() - pos_E;
  int32_t len_AC = len_A + len_C;
  int32_t len_ACD = len_AC + len_D;
  int32_t len_ACDB = len_ACD + len_B;

  //printf("%d %d %d %d -- %d %d %d %d %d\n",pos_B,pos_C,pos_D,pos_E,len_A,len_B,len_C,len_D,len_E);

#ifdef WITH_BITSET
  // Create new sequence
  BitSet_SIMD* bitset = new BitSet_SIMD(bitset_->length_);

  // Copy segment A
  bitset->replace(bitset_,0);

  // Build Cp and put it in the new genome
  BitSet_SIMD* inverted_segment = bitset_->duplicate(pos_D-len_C,pos_D,true);

  bitset->replace(inverted_segment,len_A);

  delete inverted_segment;

  // Build Dp and put it in the new genome
  inverted_segment = bitset_->duplicate(pos_E-len_D,pos_E,true);

  bitset->replace(inverted_segment,len_AC);

  delete inverted_segment;


  // Copy segments B and E
  bitset->replace(bitset_,pos_B,pos_B+len_B,len_ACD);
  bitset->replace(bitset_,pos_E,pos_E+len_E,len_ACDB);

  // Replace sequence
  delete bitset_;
  bitset_ = bitset;
#else
  // Create new sequence
  char* new_genome;
  posix_memalign((void **)&new_genome,64,nb_blocks_ * BLOCK_SIZE* sizeof(char));//new char[nb_blocks_ * BLOCK_SIZE];
  memset(new_genome,0,(length_) * sizeof(char));

  // Copy segment A
  memcpy(new_genome, data_, len_A * sizeof(char));


  // Build Cp and put it in the new genome
  char* inverted_segment;
  posix_memalign((void **)&inverted_segment,64,(len_C+1)* sizeof(char));//new char[len_C + 1];

#ifdef __SIMD
  #pragma simd
#pragma distribute_point
#endif
  for (int32_t i = 0, j = pos_D - 1; i < len_C; i++, j--) {
    if (data_[j] == '0') inverted_segment[i] = '1';
    else inverted_segment[i] = '0';
  }
  inverted_segment[len_C] = '\0';

  memcpy(&new_genome[len_A], inverted_segment, len_C * sizeof(char));

  free(inverted_segment);


  // Build Dp and put it in the new genome
  posix_memalign((void **)&inverted_segment,64,(len_D+1)* sizeof(char));//new char[len_D + 1];

#ifdef __SIMD
  #pragma simd
#pragma distribute_point
#endif
  for (int32_t i = 0, j = pos_E - 1; i < len_D; i++, j--) {
    if (data_[j] == '0') inverted_segment[i] = '1';
    else inverted_segment[i] = '0';
  }
  inverted_segment[len_D] = '\0';

  memcpy(&new_genome[len_AC], inverted_segment, len_D * sizeof(char));

  free(inverted_segment);

  // Copy segments B and E
  memcpy(&new_genome[len_ACD], &data_[pos_B], len_B * sizeof(char));
  memcpy(&new_genome[len_ACDB], &data_[pos_E], len_E * sizeof(char));
  new_genome[length_] = '\0';


  // Replace sequence
  free(data_);
  data_ = new_genome;
#endif


  // ========== Update promoter list ==========
  // 1) Remove promoters that include a breakpoint
  // 2) Extract promoters that are totally included in each segment to be
  //    moved (B, C and D)
  // 3) Shift (and invert when needed) these promoters positions
  // 4) Reinsert the shifted promoters
  // 5) Look for new promoters including a breakpoint
  if (length() >= PROM_SIZE) {

//    if (indiv_->indiv_id == 77) {
//      printf("BEFORE_REMOVE : Leading promoters lists : ");
//      for (auto it : indiv_->leading_prom_pos) {
//        printf("%d (%d) || ", it.first, it.second);
//      }
//      printf("\n");
//
//      printf("BEFORE_REMOVE : Lagging promoters lists : ");
//      for (auto it : indiv_->lagging_prom_pos) {
//        printf("%d (%d) || ", it.first, it.second);
//      }
//      printf("\n");
//
//      printf("BEFORE_REMOVE : Promoters lists (promoters): ");
//      for (auto it : indiv_->promoters) {
//        printf("%d (%d) -- ", it.second->pos, it.first);
//      }
//
//      printf("\n");
//    }

    // 1) Remove promoters that include a breakpoint
    indiv_->remove_promoters_around(pos_B);
    indiv_->remove_promoters_around(pos_C);
    indiv_->remove_promoters_around(pos_D);
    indiv_->remove_promoters_around(pos_E);

//    if (indiv_->indiv_id == 77) {
//      printf("AFTER_REMOVE : Leading promoters lists : ");
//      for (auto it : indiv_->leading_prom_pos) {
//        printf("%d (%d) || ", it.first, it.second);
//      }
//      printf("\n");
//
//      printf("AFTER_REMOVE : Lagging promoters lists : ");
//      for (auto it : indiv_->lagging_prom_pos) {
//        printf("%d (%d) || ", it.first, it.second);
//      }
//      printf("\n");
//
//      printf("AFTER_REMOVE : Promoters lists (promoters): ");
//      for (auto it : indiv_->promoters) {
//        printf("%d (%d) -- ", it.second->pos, it.first);
//      }
//
//      printf("\n");
//    }


    // Create temporary lists for promoters to move and/or invert
    std::vector<std::list<promoterStruct*>> promoters_B = {{},
                                     {}};
    std::vector<std::list<promoterStruct*>> promoters_C = {{},
                                     {}};
    std::vector<std::list<promoterStruct*>> promoters_D = {{},
                                     {}};

    // 2) Extract promoters that are totally included in each segment to be
    //    moved (B, C and D)
    if (len_B >= PROM_SIZE) {
      indiv_->extract_promoters_included_in(pos_B, pos_C, promoters_B);
    }
    if (len_C >= PROM_SIZE) {
      indiv_->extract_promoters_included_in(pos_C, pos_D, promoters_C);
    }
    if (len_D >= PROM_SIZE) {
      indiv_->extract_promoters_included_in(pos_D, pos_E, promoters_D);
    }

//    if (indiv_->indiv_id == 77) {
//      printf("AFTER_EXTRACT : Leading promoters lists : ");
//      for (auto it : indiv_->leading_prom_pos) {
//        printf("%d (%d) || ", it.first, it.second);
//      }
//      printf("\n");
//
//      printf("AFTER_EXTRACT : Lagging promoters lists : ");
//      for (auto it : indiv_->lagging_prom_pos) {
//        printf("%d (%d) || ", it.first, it.second);
//      }
//      printf("\n");
//
//      printf("AFTER_EXTRACT : Promoters lists (promoters): ");
//      for (auto it : indiv_->promoters) {
//        printf("%d (%d) -- ", it.second->pos, it.first);
//      }
//
//      printf("\n");
//
//      printf("PROMOTER_LIST BEFORE_INVERT B: ");
//      for (auto strand: {LEADING, LAGGING}) {
//        for (auto it : promoters_B[strand]) {
//          printf("%d || ", it->pos);
//        }
//      }
//      printf("\n");
//
//      printf("PROMOTER_LIST BEFORE_INVERT C: ");
//      for (auto strand: {LEADING, LAGGING}) {
//        for (auto it : promoters_C[strand]) {
//          printf("%d || ", it->pos);
//        }
//      }
//      printf("\n");
//
//      printf("PROMOTER_LIST BEFORE_INVERT D: ");
//      for (auto strand: {LEADING, LAGGING}) {
//        for (auto it : promoters_D[strand]) {
//          printf("%d || ", it->pos);
//        }
//      }
//      printf("\n");
//    }

    // 3a) Invert promoters of segments C and D
    Internal_SIMD_Struct::invert_promoters(promoters_C, pos_C, pos_D);

    Internal_SIMD_Struct::invert_promoters(promoters_D, pos_D, pos_E);

//    if (indiv_->indiv_id == 77) {
//      printf("PROMOTER_LIST AFTER_INVERT C: ");
//      for (auto strand: {LEADING, LAGGING}) {
//        for (auto it : promoters_C[strand]) {
//          printf("%d || ", it->pos);
//        }
//      }
//      printf("\n");
//
//      printf("PROMOTER_LIST AFTER_INVERT D: ");
//      for (auto strand: {LEADING, LAGGING}) {
//        for (auto it : promoters_D[strand]) {
//          printf("%d || ", it->pos);
//        }
//      }
//      printf("\n");
//    }

    // 3b) Shift these promoters positions
    Internal_SIMD_Struct::shift_promoters(promoters_B, len_C + len_D,
                                 length());

    Internal_SIMD_Struct::shift_promoters(promoters_C, -len_B,
                                 length());

    Internal_SIMD_Struct::shift_promoters(promoters_D, -len_B,
                                 length());

//    if (indiv_->indiv_id == 77) {
//      printf("PROMOTER_LIST AFTER_SHIFT B: ");
//      for (auto strand: {LEADING, LAGGING}) {
//        for (auto it : promoters_B[strand]) {
//          printf("%d || ", it->pos);
//        }
//      }
//      printf("\n");
//
//      printf("PROMOTER_LIST AFTER_SHIFT C: ");
//      for (auto strand: {LEADING, LAGGING}) {
//        for (auto it : promoters_C[strand]) {
//          printf("%d || ", it->pos);
//        }
//      }
//      printf("\n");
//
//      printf("PROMOTER_LIST AFTER_SHIFT D: ");
//      for (auto strand: {LEADING, LAGGING}) {
//        for (auto it : promoters_D[strand]) {
//          printf("%d || ", it->pos);
//        }
//      }
//      printf("\n");
//    }
    // 4) Reinsert the shifted promoters
    indiv_->insert_promoters(promoters_B);

    indiv_->insert_promoters(promoters_D);

    indiv_->insert_promoters(promoters_C);

//    if (indiv_->indiv_id == 77) {
//      printf("AFTER_INSERT : Leading promoters lists : ");
//      for (auto it : indiv_->leading_prom_pos) {
//        printf("%d (%d) || ", it.first, it.second);
//      }
//      printf("\n");
//
//      printf("AFTER_INSERT : Lagging promoters lists : ");
//      for (auto it : indiv_->lagging_prom_pos) {
//        printf("%d (%d) || ", it.first, it.second);
//      }
//      printf("\n");
//
//      printf("AFTER_INSERT : Promoters lists (promoters): ");
//      for (auto it : indiv_->promoters) {
//        printf("%d (%d) -- ", it.second->pos, it.first);
//      }
//
//      printf("\n");
//    }

    // 5) Look for new promoters including a breakpoint
    indiv_->look_for_new_promoters_around(len_A);
    indiv_->look_for_new_promoters_around(len_AC);
    indiv_->look_for_new_promoters_around(len_ACD);
    indiv_->look_for_new_promoters_around(len_ACDB);
  }
}



}
