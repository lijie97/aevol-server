// ****************************************************************************
//
//          Aevol - An in silico experimental evolution platform
//
// ****************************************************************************
//
// Copyright: See the AUTHORS file provided with the package or <www.aevol.fr>
// Web: http://www.aevol.fr/
// E-mail: See <http://www.aevol.fr/contact/>
// Original Authors : Guillaume Beslon, Carole Knibbe, David Parsons
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ****************************************************************************

#include "Dna_7.h"

#include "DnaMutator.h"
#include "DynTab_Metadata.h"
#include "ExpManager.h"
#include "Individual.h"
#include "Individual_7.h"
#include "List_Metadata.h"
#include "Map_Metadata.h"

#define REDUCTION_FACTOR 16

namespace aevol {
Dna_7::Dna_7(Dna* dna, Individual_7 * indiv, DnaFactory * dna_factory) {
  length_ = dna->length();

  nb_blocks_ = nb_blocks(length_);// _mm_free(ptr)
  posix_memalign((void **)&data_ ,64, nb_blocks_ * BLOCK_SIZE* sizeof(char)); //new char[nb_blocks_ * BLOCK_SIZE];
  memset(data_,0,(length_+1) * sizeof(char));
  memcpy(data_, dna->data(), (length_+1) * sizeof(char));

  parent_length_ = dna->length();
  indiv_ = indiv;
  dna_factory_ = dna_factory;
}

Dna_7::Dna_7(Dna* dna, DnaFactory * dna_factory) {
        length_ = dna->length();

        nb_blocks_ = nb_blocks(length_);
        posix_memalign((void **)&data_ ,64, nb_blocks_ * BLOCK_SIZE* sizeof(char));//new char[nb_blocks_ * BLOCK_SIZE];
        memset(data_,0,(length_+1) * sizeof(char));
        memcpy(data_, dna->data(), (length_+1) * sizeof(char));

        //printf("Calling LENGTH\n");
        parent_length_ = dna->length();
        //printf("Dna_7 out of constructor\n");
    dna_factory_ = dna_factory;

}

Dna_7::Dna_7(Dna_7 * dna, Individual_7 * indiv, DnaFactory * dna_factory) {
  length_ = dna->length_;

  nb_blocks_ = nb_blocks(length_);
  posix_memalign((void **)&data_ ,64, nb_blocks_ * BLOCK_SIZE* sizeof(char));//new char[nb_blocks_ * BLOCK_SIZE];
  memset(data_,0,(length_+1) * sizeof(char));
  memcpy(data_, dna->data_, (length_+1) * sizeof(char));

  parent_length_ = dna->length_;

  indiv_ = indiv;
    dna_factory_ = dna_factory;

}

Dna_7::Dna_7(int length, DnaFactory * dna_factory) {
    length_ = length;

    nb_blocks_ = nb_blocks(length_);
    posix_memalign((void **) &data_, 64, nb_blocks_ * BLOCK_SIZE * sizeof(char));//new char[nb_blocks_ * BLOCK_SIZE];
    memset(data_, 0, (length_ + 1) * sizeof(char));
    dna_factory_ = dna_factory;

}

void Dna_7::set_indiv(Dna_7 * dna, Individual_7 * indiv) {
    int new_nb_blocks = nb_blocks(dna->length_);

    //printf("Length %d %d (%d %d)\n",length_,dna->length(),nb_blocks_,dna->nb_blocks_);
    length_ = dna->length();
    if ((new_nb_blocks > nb_blocks_) || (new_nb_blocks<nb_blocks_/ REDUCTION_FACTOR)) {
        //printf("REALLOCATE DNA -- Old Block %d, New Block %d, Length %d\n",nb_blocks_,new_nb_blocks,length());
        nb_blocks_ = new_nb_blocks;
        if (data_ != nullptr) { free(data_); data_=nullptr;}

        posix_memalign((void **)&data_ ,64, nb_blocks_ * BLOCK_SIZE* sizeof(char));
        memset(data_,0,(length_+1) * sizeof(char));
    }


    memcpy(data_, dna->data_, (length_+1) * sizeof(char));

    parent_length_ = dna->length_;
    indiv_ = indiv;
    dna_factory_ = indiv->dna_factory_;


}

void Dna_7::set_indiv(Dna* dna, DnaFactory * dna_factory) {
    int new_nb_blocks = nb_blocks(dna->length());

    length_ = dna->length();

    if ((new_nb_blocks > nb_blocks_) || (new_nb_blocks<nb_blocks_/REDUCTION_FACTOR)) {
        //printf("REALLOCATE DNA -- Old Block %d, New Block %d, Length %d\n",nb_blocks_,new_nb_blocks,length());
        nb_blocks_ = new_nb_blocks;
        if (data_ != nullptr) { free(data_); data_=nullptr;}

        posix_memalign((void **)&data_ ,64, nb_blocks_ * BLOCK_SIZE* sizeof(char));
        memset(data_,0,(length_+1) * sizeof(char));
    }

    memcpy(data_, dna->data(), (length_+1) * sizeof(char));

    parent_length_ = dna->length();

    indiv_ = nullptr;

    dna_factory_ = dna_factory;

}

void Dna_7::set_indiv(int req_length, int parent_length, Individual_7 * indiv) {
    int new_nb_blocks = nb_blocks(req_length);
    length_ = req_length;

    if ((new_nb_blocks > nb_blocks_) || (new_nb_blocks < nb_blocks_ / REDUCTION_FACTOR)) {
        //printf("REALLOCATE DNA -- Old Block %d, New Block %d, Length %d\n", nb_blocks_, new_nb_blocks, length());
        nb_blocks_ = new_nb_blocks;
        if (data_ != nullptr) {
            free(data_);
            data_ = nullptr;
        }

        posix_memalign((void **) &data_, 64, nb_blocks_ * BLOCK_SIZE * sizeof(char));
        memset(data_, 0, (length_ + 1) * sizeof(char));
    }

    parent_length_ = parent_length;
    indiv_ = indiv;
    dna_factory_ = indiv->dna_factory_;

}

Dna_7::~Dna_7() {
  if (data_ != nullptr) { free(data_); data_=nullptr;}

if (mutation_list.size() > 0) {
    for (auto repl : mutation_list)
        delete repl;
    mutation_list.clear();
}


}

void Dna_7::remove(int32_t pos_1, int32_t pos_2) {
// if (indiv_->indiv_id == 586) printf("%d -- BEGIN -- REMOVE POS 1 %d POS 2 %d LENGTH %d\n",indiv_->indiv_id,pos_1,pos_2,length());

  assert(pos_1 >= 0 && pos_2 >= pos_1 && pos_2 <= length_);

  // Compute size of new genome
  int32_t new_length    = length_ - (pos_2 - pos_1);

//  printf("POS 1 %d POS 2 %d  LENGTH %d NEW Length %d Nb Block %d\n",
//         pos_1,pos_2,length(),new_length,new_nb_blocks);

  Dna_7 * new_genome = indiv_->dna_factory_->get_dna(new_length);
  //printf("Dna_7 -- REMOVE -- Start\n");
    new_genome->set_indiv(new_length,parent_length_,indiv_);

  // Copy the remaining of the genome in tmp (preceeding and following parts)
  memcpy(new_genome->data_, data_, pos_1 * sizeof(char));
  memcpy(&(new_genome->data_[pos_1]), &data_[pos_2],
         (new_length - pos_1) * sizeof(char));
  new_genome->data_[new_length] = '\0';

  // Replace previous genome with the new one
  char* old_data = data_;
    int new_nb_block = new_genome->nb_blocks_;
    new_genome->nb_blocks_ = nb_blocks_;
    data_ = new_genome->data_;
  new_genome->data_ = old_data;
  new_genome->length_ = length();

  indiv_->dna_factory_->give_back(new_genome);
   // printf("Dna_7 -- REMOVE -- End\n");

  // Update length data
  length_ = new_length;
  nb_blocks_ = new_nb_block;
//    if (indiv_->indiv_id == 586) printf("%d -- END -- REMOVE POS 1 %d POS 2 %d LENGTH %d\n",indiv_->indiv_id,pos_1,pos_2,length());
}

void Dna_7::insert(int32_t pos, const char* seq, int32_t seq_length) {
// Insert sequence 'seq' at position 'pos'
  assert(pos >= 0 && pos < length_);

  // If the sequence's length was not provided, compute it
  if (seq_length == -1) {
    seq_length = strlen(seq);
  }

  // Compute size of new genome
  int32_t new_length    = length_ + seq_length;

  Dna_7 * new_genome = indiv_->dna_factory_->get_dna(new_length);
    //printf("Dna_7 -- INSERT -- Start\n");

    new_genome->set_indiv(new_length,parent_length_,indiv_);

  // Build new genome from previous genome and sequence to insert
  memcpy(new_genome->data_, data_, pos * sizeof(char));
  memcpy(&(new_genome->data_[pos]), seq, seq_length * sizeof(char));
  memcpy(&(new_genome->data_[pos+seq_length]), &data_[pos],
         (length_ - pos) * sizeof(char));
  new_genome->data_[new_length] = '\0';

  // Replace the previous genome with the new one
    char* old_data = data_;
    data_ = new_genome->data_;
    int new_nb_block = new_genome->nb_blocks_;
    new_genome->data_ = old_data;
    new_genome->nb_blocks_ = nb_blocks_;
    new_genome->length_ = length();

    indiv_->dna_factory_->give_back(new_genome);

    //printf("Dna_7 -- INSERT -- End\n");


    // Update length-related data
  length_ = new_length;
  nb_blocks_ = new_nb_block;
}

void Dna_7::replace(int32_t pos, char* seq, int32_t seq_length) {
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


bool Dna_7::do_switch(int32_t pos) {
  // Perform the mutation
  if (data_[pos] == '0') data_[pos] = '1';
  else data_[pos] = '0';

  // Remove promoters containing the switched base
  indiv_->metadata_->remove_promoters_around(pos, Utils::mod(pos + 1, length()));

  // Look for potential new promoters containing the switched base
  if (length() >= PROM_SIZE)
      indiv_->metadata_->look_for_new_promoters_around(pos, Utils::mod(pos + 1, length()));

  if (ExpManager_7::standalone_simd && indiv_->exp_m_->record_tree()) {
      PointMutation *mut = new PointMutation(pos);
/*      printf("SEND MUTATION EVENT (%ld)\n",indiv_->observers_[MUTATION].size());*/
      //indiv_->notifyObservers(MUTATION, mut);
      indiv_->exp_m_->tree()->report_by_index(AeTime::time(),indiv_->indiv_id)->dna_replic_report().add_mut(mut);
      delete mut;
  }

  return true;
}

bool Dna_7::do_small_insertion(int32_t pos, int16_t nb_insert, char* seq) {
  // Remove the promoters that will be broken
    indiv_->metadata_->remove_promoters_around(pos);

  insert(pos, seq, nb_insert);

  // Look for new promoters
  if (length_ >= PROM_SIZE) {
    if (length_ - nb_insert < PROM_SIZE) {
      // Special case where the genome was smaller than a promoter before the
      // insertion and greater than (or as big as) a promoter after the
      // insertion.
      // In that case, we must look for new promoters thoroughly on the whole
      // genome using locate_promoters
        indiv_->metadata_->locate_promoters();
    }
    else {
        indiv_->metadata_->move_all_promoters_after(pos, nb_insert);
        indiv_->metadata_->look_for_new_promoters_around(pos, Utils::mod(pos + nb_insert,
                                                               length_));
    }
  }

    if (ExpManager_7::standalone_simd && indiv_->exp_m_->record_tree()) {
        SmallInsertion *mut = new SmallInsertion(pos, nb_insert, seq);
        //indiv_->notifyObservers(MUTATION, mut);
        indiv_->exp_m_->tree()->report_by_index(AeTime::time(),indiv_->indiv_id)->dna_replic_report().add_mut(mut);
        delete mut;
    }

  return true;
}

bool Dna_7::do_small_deletion(int32_t pos, int16_t nb_del) {
    int32_t old_pos = pos;
  // Remove promoters containing at least one nucleotide from the sequence to
  // delete
    indiv_->metadata_->remove_promoters_around(pos, Utils::mod(pos + nb_del, length()));

  // Do the deletion and update promoter list
  if (pos + nb_del <= length()) { // the deletion does not contain the origin of
    // replication
    // Do the deletion
    remove(pos, pos + nb_del);

    // Update promoter list
    if (length() >= PROM_SIZE) {
        indiv_->metadata_->move_all_promoters_after(pos, -nb_del);
        indiv_->metadata_->look_for_new_promoters_around(Utils::mod(pos, length()));
    }
  }
  else { // the deletion contains the origin of replication
    // Do the deletion
    int32_t nb_del_at_pos_0 = nb_del - length() + pos;

    remove(pos, length_);
    remove(0, nb_del_at_pos_0);


    pos -= nb_del_at_pos_0;

    // Update promoter list
    if (length() >= PROM_SIZE) {
        indiv_->metadata_->move_all_promoters_after(0, -nb_del_at_pos_0);
        indiv_->metadata_->look_for_new_promoters_around(0);
    }
  }

    if (ExpManager_7::standalone_simd && indiv_->exp_m_->record_tree()) {
        SmallDeletion *mut = new SmallDeletion(old_pos, nb_del);
        //indiv_->notifyObservers(MUTATION, mut);
        indiv_->exp_m_->tree()->report_by_index(AeTime::time(),indiv_->indiv_id)->dna_replic_report().add_mut(mut);
        delete mut;
    }

  return true;
}


bool Dna_7::do_duplication(int32_t pos_1, int32_t pos_2, int32_t pos_3) {
// Duplicate segment [pos_1; pos_2[ and insert the duplicate before pos_3
  char* duplicate_segment = NULL;

    if (length_ == 1)
    {
        printf("*** genome of size 1 ; duplication not done *** \n");
        return false;
    }

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

    posix_memalign((void **)&duplicate_segment,64,(seg_length+1)* sizeof(char));//new char[seg_length + 1];
    memcpy(duplicate_segment, &data_[pos_1], seg_length);
    duplicate_segment[seg_length] = '\0';

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

    posix_memalign((void **)&duplicate_segment,64,(seg_length+1)* sizeof(char));//new char[seg_length + 1];
    memcpy(duplicate_segment, &data_[pos_1], tmp1_len);     // Copy tmp1
    memcpy(&duplicate_segment[tmp1_len], data_, tmp2_len);  // Copy tmp2
    duplicate_segment[seg_length] = '\0';
  }

/*  if (indiv_->indiv_id==49) {
    printf(
        "Duplication _lengh %d\n", seg_length);
  }*/

  if (seg_length <= 0) {
    free(duplicate_segment);
    return true;
  }

  // Create a copy of the promoters beared by the segment to be duplicated
  // (they will be inserted in the individual's RNA list later)
  std::vector<std::list<PromoterStruct *>> duplicated_promoters = {{},
                                            {}};
    indiv_->metadata_->duplicate_promoters_included_in(pos_1, pos_2,
                                             duplicated_promoters);


//   printf("Duplicate %d %d %d -- %d %d\n",pos_1,pos_2,pos_3,duplicated_promoters[0].size(),duplicated_promoters[1].size());

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
//    printf("Prom list LEAD : ");
//    for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//        if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//            if (indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                printf("%d ",indiv_->metadata_->promoters(prom_idx)->pos);
//    }
//    printf("\n");
//    printf("Prom list LAG : ");
//    for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//        if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//            if (!indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                printf("%d ",indiv_->metadata_->promoters(prom_idx)->pos);
//    }
//    printf("\n");

//  printf("%d -- %d -- DUP-1 -- Number of RNAs %d (%d)\n",time(),indiv_->indiv_id,indiv_->metadata_->rna_count(),
//           indiv_->metadata_->promoter_count());

    indiv_->metadata_->remove_promoters_around(pos_3);

//    printf("%d -- %d -- DUP-2 -- Number of RNAs %d (%d)\n",time(),indiv_->indiv_id,indiv_->metadata_->rna_count(),
//           indiv_->metadata_->promoter_count());

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

  insert(pos_3, duplicate_segment, seg_length);

  if (length() >= PROM_SIZE) {
    if (length() - seg_length < PROM_SIZE) {
        for (auto strand: {LEADING, LAGGING}) {
            for (std::list<PromoterStruct *>::iterator it_rna = duplicated_promoters[strand].begin();
                    it_rna != duplicated_promoters[strand].end(); it_rna++) {
                delete (*(it_rna));
            }
        }

      // Special case where the genome was smaller than a promoter before
      // the insertion and greater than (or as big as) a promoter after the
      // insertion.
      // In that case, we must look for new promoters thoroughly on the whole
      // genome using locate_promoters
        indiv_->metadata_->locate_promoters();
    }
    else {
//        printf("Prom list LEAD : ");
//        for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//            if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                if (indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                    printf("%d ",indiv_->metadata_->promoters(prom_idx)->pos);
//        }
//        printf("\n");
//        printf("Prom list LAG : ");
//        for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//            if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                if (!indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                    printf("%d ",indiv_->metadata_->promoters(prom_idx)->pos);
//        }
//        printf("\n");

        indiv_->metadata_->move_all_promoters_after(pos_3, seg_length);
//        printf("Prom list LEAD : ");
//        for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//            if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                if (indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                    printf("%d ",indiv_->metadata_->promoters(prom_idx)->pos);
//        }
//        printf("\n");
//        printf("Prom list LAG : ");
//        for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//            if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                if (!indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                    printf("%d ",indiv_->metadata_->promoters(prom_idx)->pos);
//        }
//        printf("\n");


//        printf("Duplicated PROMs LEAD : ");
//        for (auto prom :duplicated_promoters[LEADING]) {
//            if (prom != nullptr)
//                if (prom->leading_or_lagging)
//                    printf("%d ",prom->pos);
//        }
//        printf("\n");
//        printf("Duplicated PROMs LAG : ");
//        for (auto prom :duplicated_promoters[LAGGING]) {
//            if (prom != nullptr)
//                if (!prom->leading_or_lagging)
//                    printf("%d ",prom->pos);
//        }
//        printf("\n");
//
//
//        printf("%d -- %d -- DUP-3 -- Number of RNAs %d (%d)\n",time(),indiv_->indiv_id,indiv_->metadata_->rna_count(),
//               indiv_->metadata_->promoter_count());

        indiv_->metadata_->insert_promoters_at(duplicated_promoters, pos_3);
        /*for (auto strand: {LEADING, LAGGING}) {
            for (std::list<PromoterStruct *>::iterator it_rna = duplicated_promoters[strand].begin();
                 it_rna != duplicated_promoters[strand].end(); it_rna++) {
                delete (*(it_rna));
            }
        }*/

//        printf("%d -- %d -- DUP-4 -- Number of RNAs %d (%d)\n",time(),indiv_->indiv_id,indiv_->metadata_->rna_count(),
//               indiv_->metadata_->promoter_count());
//        printf("Prom list LEAD : ");
//        for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//            if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                if (indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                    printf("%d ",indiv_->metadata_->promoters(prom_idx)->pos);
//        }
//        printf("\n");
//        printf("Prom list LAG : ");
//        for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//            if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                if (!indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                    printf("%d ",indiv_->metadata_->promoters(prom_idx)->pos);
//        }
//        printf("\n");
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
//
//        printf("%d -- %d -- DUP-5 -- Number of RNAs %d (%d)\n",time(),indiv_->indiv_id,indiv_->metadata_->rna_count(),
//               indiv_->metadata_->promoter_count());

        indiv_->metadata_->look_for_new_promoters_around(pos_3);

//
//        printf("%d -- %d -- DUP-6 -- Number of RNAs %d (%d)\n",time(),indiv_->indiv_id,indiv_->metadata_->rna_count(),
//               indiv_->metadata_->promoter_count());


        indiv_->metadata_->look_for_new_promoters_around(pos_3 + seg_length);


//        printf("%d -- %d -- DUP-7 -- Number of RNAs %d (%d)\n",time(),indiv_->indiv_id,indiv_->metadata_->rna_count(),
//               indiv_->metadata_->promoter_count());

//        printf("Prom list LEAD : ");
//        for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//            if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                if (indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                    printf("%d ",indiv_->metadata_->promoters(prom_idx)->pos);
//        }
//        printf("\n");
//        printf("Prom list LAG : ");
//        for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//            if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                if (!indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                    printf("%d ",indiv_->metadata_->promoters(prom_idx)->pos);
//        }
//        printf("\n");
        for (auto strand: {LEADING, LAGGING}) {
            for (std::list<PromoterStruct *>::iterator it_rna = duplicated_promoters[strand].begin();
                 it_rna != duplicated_promoters[strand].end(); it_rna++) {
                delete (*(it_rna));
            }
        }
    }
  } else {
      for (auto strand: {LEADING, LAGGING}) {
          for (std::list<PromoterStruct *>::iterator it_rna = duplicated_promoters[strand].begin();
               it_rna != duplicated_promoters[strand].end(); it_rna++) {
              delete (*(it_rna));
          }
      }
  }

  free(duplicate_segment);

    if (ExpManager_7::standalone_simd && indiv_->exp_m_->record_tree()) {
        Duplication *mut = new Duplication(pos_1, pos_2, pos_3, seg_length);
        indiv_->exp_m_->tree()->report_by_index(AeTime::time(),indiv_->indiv_id)->dna_replic_report().add_mut(mut);
        //indiv_->notifyObservers(MUTATION, mut);
        delete mut;
    }
  return true;
}

bool Dna_7::do_translocation(int32_t pos_1, int32_t pos_2, int32_t pos_3,
                           int32_t pos_4, bool invert) {
  int32_t pos_min = Utils::min(pos_1,
                               Utils::min(pos_2, Utils::min(pos_3, pos_4)));

  if (not invert) {
    if (pos_min == pos_1) {
      //if (indiv_->indiv_id==799) printf("A\n");
      ABCDE_to_ADCBE(pos_1, pos_3, pos_2, pos_4);
    }
    else if (pos_min == pos_2) {
        //if (indiv_->indiv_id==799) printf("B\n");
      ABCDE_to_ADCBE(pos_2, pos_4, pos_1, pos_3);
    }
    else if (pos_min == pos_3) {
       // if (indiv_->indiv_id==799) printf("C\n");
      ABCDE_to_ADCBE(pos_3, pos_2, pos_4, pos_1);
    }
    else { // if (pos_min == pos_4)
        //if (indiv_->indiv_id==799) printf("D\n");
      ABCDE_to_ADCBE(pos_4, pos_1, pos_3, pos_2);
    }
  }
  else { // invert
    if (pos_min == pos_1) {
        //if (indiv_->indiv_id==799) printf("E\n");
      ABCDE_to_ADBpCpE(pos_1, pos_3, pos_2, pos_4);
    }
    else if (pos_min == pos_2) {
        //if (indiv_->indiv_id==799) printf("F\n");
      ABCDE_to_ADBpCpE(pos_2, pos_4, pos_1, pos_3);
    }
    else if (pos_min == pos_3) {
        //if (indiv_->indiv_id==799) printf("G\n");
      ABCDE_to_ACpDpBE(pos_3, pos_2, pos_4, pos_1);
    }
    else { // if (pos_min == pos_4)
        //if (indiv_->indiv_id==799) printf("H\n");
      ABCDE_to_ACpDpBE(pos_4, pos_1, pos_3, pos_2);
    }
  }

  int32_t segment_length = pos_2 - pos_1;

    if (ExpManager_7::standalone_simd && indiv_->exp_m_->record_tree()) {
        Translocation *mut = new Translocation(pos_1, pos_2, pos_3, pos_4,
                                               segment_length, invert);
        //indiv_->notifyObservers(MUTATION, mut);
        indiv_->exp_m_->tree()->report_by_index(AeTime::time(),indiv_->indiv_id)->dna_replic_report().add_mut(mut);
        delete mut;
    }

  return true;
}


/*
bool Dna_7::do_translocation(int32_t pos_1, int32_t pos_2, int32_t pos_3,
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

bool Dna_7::do_inversion(int32_t pos_1, int32_t pos_2) {
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

        if (length_ == 1)
        {
            printf("*** genome of size 1 ; inversion not done *** \n");
            return false;
        }

  int32_t seg_length = pos_2 - pos_1;

  // Create the inverted sequence
  char* inverted_segment = NULL;
  posix_memalign((void **)&inverted_segment,64,(seg_length+1)* sizeof(char));//new char[seg_length + 1];

  for (int32_t i = 0, j = pos_2 - 1; i < seg_length; i++, j--) {
    if (data_[j] == '0') inverted_segment[i] = '1';
    else inverted_segment[i] = '0';
  }
  inverted_segment[seg_length] = '\0';


        /*if (indiv_->indiv_id == 37) {
            printf("Before remove %d\n",pos_1,pos_2);
            printf("%d -- %d -- TOKEEP -- Prom list LEAD : ", time(), indiv_->indiv_id);
            for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
                if (indiv_->metadata_->promoters(prom_idx) != nullptr)
                    if (indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
                        printf("%d ", indiv_->metadata_->promoters(prom_idx)->pos);
            }
            printf("\n");
            printf("%d -- %d -- TOKEEP -- Prom list LAG : ", time(), indiv_->indiv_id);
            for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
                if (indiv_->metadata_->promoters(prom_idx) != nullptr)
                    if (!indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
                        printf("%d ", indiv_->metadata_->promoters(prom_idx)->pos);
            }
            printf("\n");

            printf("%d -- %d -- ADV-DIRECT -- Prom list LEAD : ", time(), indiv_->indiv_id);
            for (int prom_idx = 0; prom_idx < ((List_Metadata*)indiv_->metadata_)->promoters_list_[LEADING].size(); prom_idx++) {
                auto it =  ((List_Metadata*)indiv_->metadata_)->promoters_list_[LEADING].begin();
                std::advance(it, prom_idx);
                printf("%d ", (*it).pos);
            }
            printf("\n");
            printf("%d -- %d -- ADV-DIRECT -- Prom list LAG : ", time(), indiv_->indiv_id);
            for (int prom_idx = 0; prom_idx < ((List_Metadata*)indiv_->metadata_)->promoters_list_[LAGGING].size(); prom_idx++) {
                auto it =  ((List_Metadata*)indiv_->metadata_)->promoters_list_[LAGGING].begin();
                std::advance(it, prom_idx);
                printf("%d ", (*it).pos);
            }
            printf("\n");
        }*/

  // Remove promoters that included a breakpoint
        indiv_->metadata_->remove_promoters_around(pos_1);

       /* if (indiv_->indiv_id == 37) {

            printf("Before remove 2\n");
        printf("%d -- %d -- TOKEEP -- Prom list LEAD : ", time(), indiv_->indiv_id);
        for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
            if (indiv_->metadata_->promoters(prom_idx) != nullptr)
                if (indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
                    printf("%d ", indiv_->metadata_->promoters(prom_idx)->pos);
        }
        printf("\n");
        printf("%d -- %d -- TOKEEP -- Prom list LAG : ", time(), indiv_->indiv_id);
        for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
            if (indiv_->metadata_->promoters(prom_idx) != nullptr)
                if (!indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
                    printf("%d ", indiv_->metadata_->promoters(prom_idx)->pos);
        }
        printf("\n");

        printf("%d -- %d -- ADV-DIRECT -- Prom list LEAD : ", time(), indiv_->indiv_id);
        for (int prom_idx = 0; prom_idx < ((List_Metadata*)indiv_->metadata_)->promoters_list_[LEADING].size(); prom_idx++) {
            auto it =  ((List_Metadata*)indiv_->metadata_)->promoters_list_[LEADING].begin();
            std::advance(it, prom_idx);
            printf("%d ", (*it).pos);
        }
        printf("\n");
        printf("%d -- %d -- ADV-DIRECT -- Prom list LAG : ", time(), indiv_->indiv_id);
        for (int prom_idx = 0; prom_idx < ((List_Metadata*)indiv_->metadata_)->promoters_list_[LAGGING].size(); prom_idx++) {
            auto it =  ((List_Metadata*)indiv_->metadata_)->promoters_list_[LAGGING].begin();
            std::advance(it, prom_idx);
            printf("%d ", (*it).pos);
        }
        printf("\n");
    }*/


        indiv_->metadata_->remove_promoters_around(pos_2);
        /*if (indiv_->indiv_id == 37) {

            printf("Before remove 3\n");
            printf("%d -- %d -- TOKEEP -- Prom list LEAD : ", time(), indiv_->indiv_id);
            for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
                if (indiv_->metadata_->promoters(prom_idx) != nullptr)
                    if (indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
                        printf("%d ", indiv_->metadata_->promoters(prom_idx)->pos);
            }
            printf("\n");
            printf("%d -- %d -- TOKEEP -- Prom list LAG : ", time(), indiv_->indiv_id);
            for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
                if (indiv_->metadata_->promoters(prom_idx) != nullptr)
                    if (!indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
                        printf("%d ", indiv_->metadata_->promoters(prom_idx)->pos);
            }
            printf("\n");

            printf("%d -- %d -- ADV-DIRECT -- Prom list LEAD : ", time(), indiv_->indiv_id);
            for (int prom_idx = 0; prom_idx < ((List_Metadata*)indiv_->metadata_)->promoters_list_[LEADING].size(); prom_idx++) {
                auto it =  ((List_Metadata*)indiv_->metadata_)->promoters_list_[LEADING].begin();
                std::advance(it, prom_idx);
                printf("%d ", (*it).pos);
            }
            printf("\n");
            printf("%d -- %d -- ADV-DIRECT -- Prom list LAG : ", time(), indiv_->indiv_id);
            for (int prom_idx = 0; prom_idx < ((List_Metadata*)indiv_->metadata_)->promoters_list_[LAGGING].size(); prom_idx++) {
                auto it =  ((List_Metadata*)indiv_->metadata_)->promoters_list_[LAGGING].begin();
                std::advance(it, prom_idx);
                printf("%d ", (*it).pos);
            }
            printf("\n");
        }*/
  // Invert the sequence
  replace(pos_1, inverted_segment, seg_length);

  // Update promoter list
  if (length() >= PROM_SIZE) {
      indiv_->metadata_->invert_promoters_included_in(pos_1, pos_2);

      /*if (indiv_->indiv_id == 37) {

      printf("Before look for 1+2\n");

          printf("%d -- %d -- TOKEEP -- Prom list LEAD : ", time(), indiv_->indiv_id);
          for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
              if (indiv_->metadata_->promoters(prom_idx) != nullptr)
                  if (indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
                      printf("%d ", indiv_->metadata_->promoters(prom_idx)->pos);
          }
          printf("\n");
          printf("%d -- %d -- TOKEEP -- Prom list LAG : ", time(), indiv_->indiv_id);
          for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
              if (indiv_->metadata_->promoters(prom_idx) != nullptr)
                  if (!indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
                      printf("%d ", indiv_->metadata_->promoters(prom_idx)->pos);
          }
          printf("\n");

          printf("%d -- %d -- ADV-DIRECT -- Prom list LEAD : ", time(), indiv_->indiv_id);
          for (int prom_idx = 0; prom_idx < ((List_Metadata*)indiv_->metadata_)->promoters_list_[LEADING].size(); prom_idx++) {
              auto it =  ((List_Metadata*)indiv_->metadata_)->promoters_list_[LEADING].begin();
              std::advance(it, prom_idx);
              printf("%d ", (*it).pos);
          }
          printf("\n");
          printf("%d -- %d -- ADV-DIRECT -- Prom list LAG : ", time(), indiv_->indiv_id);
          for (int prom_idx = 0; prom_idx < ((List_Metadata*)indiv_->metadata_)->promoters_list_[LAGGING].size(); prom_idx++) {
              auto it =  ((List_Metadata*)indiv_->metadata_)->promoters_list_[LAGGING].begin();
              std::advance(it, prom_idx);
              printf("%d ", (*it).pos);
          }
          printf("\n");
      }*/

      indiv_->metadata_->look_for_new_promoters_around(pos_1);
      indiv_->metadata_->look_for_new_promoters_around(pos_2);
      /*if (indiv_->indiv_id == 37) {

      printf("Finish\n");

          printf("%d -- %d -- TOKEEP -- Prom list LEAD : ", time(), indiv_->indiv_id);
          for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
              if (indiv_->metadata_->promoters(prom_idx) != nullptr)
                  if (indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
                      printf("%d ", indiv_->metadata_->promoters(prom_idx)->pos);
          }
          printf("\n");
          printf("%d -- %d -- TOKEEP -- Prom list LAG : ", time(), indiv_->indiv_id);
          for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
              if (indiv_->metadata_->promoters(prom_idx) != nullptr)
                  if (!indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
                      printf("%d ", indiv_->metadata_->promoters(prom_idx)->pos);
          }
          printf("\n");

          printf("%d -- %d -- ADV-DIRECT -- Prom list LEAD : ", time(), indiv_->indiv_id);
          for (int prom_idx = 0; prom_idx < ((List_Metadata*)indiv_->metadata_)->promoters_list_[LEADING].size(); prom_idx++) {
              auto it =  ((List_Metadata*)indiv_->metadata_)->promoters_list_[LEADING].begin();
              std::advance(it, prom_idx);
              printf("%d ", (*it).pos);
          }
          printf("\n");
          printf("%d -- %d -- ADV-DIRECT -- Prom list LAG : ", time(), indiv_->indiv_id);
          for (int prom_idx = 0; prom_idx < ((List_Metadata*)indiv_->metadata_)->promoters_list_[LAGGING].size(); prom_idx++) {
              auto it =  ((List_Metadata*)indiv_->metadata_)->promoters_list_[LAGGING].begin();
              std::advance(it, prom_idx);
              printf("%d ", (*it).pos);
          }
          printf("\n");
      }*/
  }

  free(inverted_segment);

  if (ExpManager_7::standalone_simd && indiv_->exp_m_->record_tree()) {
      Inversion *mut = new Inversion(pos_1, pos_2, seg_length);
      //indiv_->notifyObservers(MUTATION, mut);
      indiv_->exp_m_->tree()->report_by_index(AeTime::time(),indiv_->indiv_id)->dna_replic_report().add_mut(mut);
      delete mut;
  }

  return true;
}


bool Dna_7::do_deletion(int32_t pos_1, int32_t pos_2) {
    //if (indiv_->indiv_id==30) printf("DO DELETION is %d %d -- %d\n",pos_1,pos_2,length());
  // Delete segment going from pos_1 (included) to pos_2 (excluded)
    if (length_ == 1)
    {
        printf("*** genome of size 1 ; deletion not done *** \n");
        return false;
    }

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

    /*if (indiv_->indiv_id == 30) {
      printf("Before remove %d %d\n",pos_1,pos_2);
        printf("%d -- %d -- ADV-DIRECT -- Prom list LEAD : ", time(), indiv_->indiv_id);
        for (int prom_idx = 0; prom_idx < ((List_Metadata*)indiv_->metadata_)->promoters_list_[LEADING].size(); prom_idx++) {
            auto it =  ((List_Metadata*)indiv_->metadata_)->promoters_list_[LEADING].begin();
            std::advance(it, prom_idx);
            printf("%d ", (*it).pos);
        }
        printf("\n");
        printf("%d -- %d -- ADV-DIRECT -- Prom list LAG : ", time(), indiv_->indiv_id);
        for (int prom_idx = 0; prom_idx < ((List_Metadata*)indiv_->metadata_)->promoters_list_[LAGGING].size(); prom_idx++) {
            auto it =  ((List_Metadata*)indiv_->metadata_)->promoters_list_[LAGGING].begin();
            std::advance(it, prom_idx);
            printf("%d ", (*it).pos);
        }
        printf("\n");
    }*/
    // Remove promoters containing at least one nucleotide from the sequence
    // to delete
      indiv_->metadata_->remove_promoters_around(pos_1, pos_2);

      /*if (indiv_->indiv_id == 30)  {
          printf("After remove\n");
          printf("%d -- %d -- ADV-DIRECT -- Prom list LEAD : ", time(), indiv_->indiv_id);
          for (int prom_idx = 0; prom_idx < ((List_Metadata*)indiv_->metadata_)->promoters_list_[LEADING].size(); prom_idx++) {
              auto it =  ((List_Metadata*)indiv_->metadata_)->promoters_list_[LEADING].begin();
              std::advance(it, prom_idx);
              printf("%d ", (*it).pos);
          }
          printf("\n");
          printf("%d -- %d -- ADV-DIRECT -- Prom list LAG : ", time(), indiv_->indiv_id);
          for (int prom_idx = 0; prom_idx < ((List_Metadata*)indiv_->metadata_)->promoters_list_[LAGGING].size(); prom_idx++) {
              auto it =  ((List_Metadata*)indiv_->metadata_)->promoters_list_[LAGGING].begin();
              std::advance(it, prom_idx);
              printf("%d ", (*it).pos);
          }
          printf("\n");
      }*/

    // Delete the sequence between pos_1 and pos_2
    remove(pos_1, pos_2);


    // Update promoter list
    if (length() >= PROM_SIZE) {
        indiv_->metadata_->move_all_promoters_after(pos_1, -segment_length);
        /*if (indiv_->indiv_id == 30) {
            printf("After move %d %d\n",pos_1, -segment_length);
            printf("%d -- %d -- ADV-DIRECT -- Prom list LEAD : ", time(), indiv_->indiv_id);
            for (int prom_idx = 0; prom_idx < ((List_Metadata*)indiv_->metadata_)->promoters_list_[LEADING].size(); prom_idx++) {
                auto it =  ((List_Metadata*)indiv_->metadata_)->promoters_list_[LEADING].begin();
                std::advance(it, prom_idx);
                printf("%d ", (*it).pos);
            }
            printf("\n");
            printf("%d -- %d -- ADV-DIRECT -- Prom list LAG : ", time(), indiv_->indiv_id);
            for (int prom_idx = 0; prom_idx < ((List_Metadata*)indiv_->metadata_)->promoters_list_[LAGGING].size(); prom_idx++) {
                auto it =  ((List_Metadata*)indiv_->metadata_)->promoters_list_[LAGGING].begin();
                std::advance(it, prom_idx);
                printf("%d ", (*it).pos);
            }
            printf("\n");
        }*/
        indiv_->metadata_->look_for_new_promoters_around(pos_1);
        /*if (indiv_->indiv_id == 30) {
            printf("After look\n");
            printf("%d -- %d -- ADV-DIRECT -- Prom list LEAD : ", time(), indiv_->indiv_id);
            for (int prom_idx = 0; prom_idx < ((List_Metadata*)indiv_->metadata_)->promoters_list_[LEADING].size(); prom_idx++) {
                auto it =  ((List_Metadata*)indiv_->metadata_)->promoters_list_[LEADING].begin();
                std::advance(it, prom_idx);
                printf("%d ", (*it).pos);
            }
            printf("\n");
            printf("%d -- %d -- ADV-DIRECT -- Prom list LAG : ", time(), indiv_->indiv_id);
            for (int prom_idx = 0; prom_idx < ((List_Metadata*)indiv_->metadata_)->promoters_list_[LAGGING].size(); prom_idx++) {
                auto it =  ((List_Metadata*)indiv_->metadata_)->promoters_list_[LAGGING].begin();
                std::advance(it, prom_idx);
                printf("%d ", (*it).pos);
            }
            printf("\n");
        }*/
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
     /*if (indiv_->indiv_id == 30) {
          printf("Before remove %d %d\n",pos_1,pos_2);
         printf("%d -- %d -- ADV-DIRECT -- Prom list LEAD : ", time(), indiv_->indiv_id);
         for (int prom_idx = 0; prom_idx < ((List_Metadata*)indiv_->metadata_)->promoters_list_[LEADING].size(); prom_idx++) {
             auto it =  ((List_Metadata*)indiv_->metadata_)->promoters_list_[LEADING].begin();
             std::advance(it, prom_idx);
             printf("%d ", (*it).pos);
         }
         printf("\n");
         printf("%d -- %d -- ADV-DIRECT -- Prom list LAG : ", time(), indiv_->indiv_id);
         for (int prom_idx = 0; prom_idx < ((List_Metadata*)indiv_->metadata_)->promoters_list_[LAGGING].size(); prom_idx++) {
             auto it =  ((List_Metadata*)indiv_->metadata_)->promoters_list_[LAGGING].begin();
             std::advance(it, prom_idx);
             printf("%d ", (*it).pos);
         }
         printf("\n");
      }*/
      indiv_->metadata_->remove_promoters_around(pos_1, pos_2);

         /*if (indiv_->indiv_id == 30) {
          printf("After remove\n");
             printf("%d -- %d -- ADV-DIRECT -- Prom list LEAD : ", time(), indiv_->indiv_id);
             for (int prom_idx = 0; prom_idx < ((List_Metadata*)indiv_->metadata_)->promoters_list_[LEADING].size(); prom_idx++) {
                 auto it =  ((List_Metadata*)indiv_->metadata_)->promoters_list_[LEADING].begin();
                 std::advance(it, prom_idx);
                 printf("%d ", (*it).pos);
             }
             printf("\n");
             printf("%d -- %d -- ADV-DIRECT -- Prom list LAG : ", time(), indiv_->indiv_id);
             for (int prom_idx = 0; prom_idx < ((List_Metadata*)indiv_->metadata_)->promoters_list_[LAGGING].size(); prom_idx++) {
                 auto it =  ((List_Metadata*)indiv_->metadata_)->promoters_list_[LAGGING].begin();
                 std::advance(it, prom_idx);
                 printf("%d ", (*it).pos);
             }
             printf("\n");
      }*/

    // Delete the sequence between pos_1 and pos_2
    remove(pos_1, length_); // delete tmp1 from genome
    remove(0, pos_2);       // delete tmp2 from genome

    // Update promoter list
    if (length() >= PROM_SIZE) {
        /*if (indiv_->indiv_id == 30) {
            printf("Before move XA %d %d\n",0,-pos_2);
            printf("%d -- %d -- ADV-DIRECT -- Prom list LEAD : ", time(), indiv_->indiv_id);
            for (int prom_idx = 0; prom_idx < ((List_Metadata*)indiv_->metadata_)->promoters_list_[LEADING].size(); prom_idx++) {
                auto it =  ((List_Metadata*)indiv_->metadata_)->promoters_list_[LEADING].begin();
                std::advance(it, prom_idx);
                printf("%d ", (*it).pos);
            }
            printf("\n");
            printf("%d -- %d -- ADV-DIRECT -- Prom list LAG : ", time(), indiv_->indiv_id);
            for (int prom_idx = 0; prom_idx < ((List_Metadata*)indiv_->metadata_)->promoters_list_[LAGGING].size(); prom_idx++) {
                auto it =  ((List_Metadata*)indiv_->metadata_)->promoters_list_[LAGGING].begin();
                std::advance(it, prom_idx);
                printf("%d ", (*it).pos);
            }
            printf("\n");
        }*/
        indiv_->metadata_->move_all_promoters_after(0, -pos_2);
        /*if (indiv_->indiv_id == 30) {
            printf("After move XA %d %d\n",0,-pos_2);
            printf("%d -- %d -- ADV-DIRECT -- Prom list LEAD : ", time(), indiv_->indiv_id);
            for (int prom_idx = 0; prom_idx < ((List_Metadata*)indiv_->metadata_)->promoters_list_[LEADING].size(); prom_idx++) {
                auto it =  ((List_Metadata*)indiv_->metadata_)->promoters_list_[LEADING].begin();
                std::advance(it, prom_idx);
                printf("%d ", (*it).pos);
            }
            printf("\n");
            printf("%d -- %d -- ADV-DIRECT -- Prom list LAG : ", time(), indiv_->indiv_id);
            for (int prom_idx = 0; prom_idx < ((List_Metadata*)indiv_->metadata_)->promoters_list_[LAGGING].size(); prom_idx++) {
                auto it =  ((List_Metadata*)indiv_->metadata_)->promoters_list_[LAGGING].begin();
                std::advance(it, prom_idx);
                printf("%d ", (*it).pos);
            }
            printf("\n");
        }*/

        indiv_->metadata_->look_for_new_promoters_around(0);
    }
  }
    if (ExpManager_7::standalone_simd && indiv_->exp_m_->record_tree()) {
        int32_t segment_length = Utils::mod(pos_2 - pos_1 - 1, length()) + 1;
        Deletion *mut = new Deletion(pos_1, pos_2, segment_length);
        indiv_->exp_m_->tree()->report_by_index(AeTime::time(),indiv_->indiv_id)->dna_replic_report().add_mut(mut);//notifyObservers(MUTATION, mut);
        delete mut;
    }

  return true;
}

void Dna_7::apply_mutations_standalone() {
  MutationEvent* repl;
//    printf("%d -- %d -- AMS-1 -- Number of RNAs %d (%d)\n",time(),indiv_->indiv_id,indiv_->metadata_->rna_count(),
//           indiv_->metadata_->promoter_count());
//    printf("%d -- %d -- BM -- Prom list LEAD : ",time(),indiv_->indiv_id);
//    for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//        if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//            if (indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                printf("%d ",indiv_->metadata_->promoters(prom_idx)->pos);
//    }
//    printf("\n");
//    printf("%d -- %d -- BM -- Prom list LAG : ",time(),indiv_->indiv_id);
//    for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//        if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//            if (!indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                printf("%d ",indiv_->metadata_->promoters(prom_idx)->pos);
//    }
//    printf("\n");
//    if (indiv_->indiv_id == 6) {
//        printf("%d -- %d -- Parent %d\n",time(),indiv_->indiv_id,indiv_->parent_id);
//        printf("%d -- %d -- BEFORE -- Prom list LEAD : ",time(),indiv_->indiv_id);
//        for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//            if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                if (indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                    printf("%d ",indiv_->metadata_->promoters(prom_idx)->pos);
//        }
//        printf("\n");
//        printf("%d -- %d -- BEFORE -- Prom list LAG : ",time(),indiv_->indiv_id);
//        for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//            if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                if (!indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                    printf("%d ",indiv_->metadata_->promoters(prom_idx)->pos);
//        }
//        printf("\n");
//    }

    /*
    if (AeTime::time() > 182) {
        //printf("Parent ID %d\n",indiv_->parent_id);
        /*printf("Promoter size %ld %ld %ld\n", indiv_->metadata_->promoter_count(),
               ((List_Metadata*)indiv_->metadata_)->promoters_list_[LEADING].size(),
                       ((List_Metadata*)indiv_->metadata_)->promoters_list_[LAGGING].size());*/
/*
        std::set<int> leading;
        printf("FROM_MUTATE -- %d -- %d -- AFTER -- Prom list LEAD : ",time(),indiv_->indiv_id);
        for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
            if (indiv_->metadata_->promoters(prom_idx) != nullptr)
                if (indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
                    leading.insert(indiv_->metadata_->promoters(prom_idx)->pos);
        }
        for (auto lead : leading) {
            printf("%d ",lead);
        }


        printf("\n");

        std::set<int> lagging;
        printf("FROM_MUTATE -- %d -- %d -- AFTER -- Prom list LAG : ",time(),indiv_->indiv_id);
        for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
            if (indiv_->metadata_->promoters(prom_idx) != nullptr)
                if (!indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
                    lagging.insert(indiv_->metadata_->promoters(prom_idx)->pos);
        }

        for (auto lag : lagging) {
            printf("%d ",lag);
        }
        printf("\n");
*/
/*        printf("%d -- %d -- ADV-DIRECT -- Prom list LEAD : ", time(), indiv_->indiv_id);
        for (int prom_idx = 0; prom_idx < ((List_Metadata*)indiv_->metadata_)->promoters_list_[LEADING].size(); prom_idx++) {
            auto it =  ((List_Metadata*)indiv_->metadata_)->promoters_list_[LEADING].begin();
            std::advance(it, prom_idx);
            printf("%d ", (*it).pos);
        }
        printf("\n");
        printf("%d -- %d -- ADV-DIRECT -- Prom list LAG : ", time(), indiv_->indiv_id);
        for (int prom_idx = 0; prom_idx < ((List_Metadata*)indiv_->metadata_)->promoters_list_[LAGGING].size(); prom_idx++) {
            auto it =  ((List_Metadata*)indiv_->metadata_)->promoters_list_[LAGGING].begin();
            std::advance(it, prom_idx);
            printf("%d ",(*it).pos);
        }
        printf("\n");*/
    //}
  do {
      repl = indiv_->exp_m_->
              dna_mutator_array_[indiv_->indiv_id]->generate_next_mutation(length());

      if (repl != nullptr) {
          //if (indiv_->indiv_id == 30 || indiv_->indiv_id==61 || indiv_->indiv_id==62)
          //printf("%d -- %d -- Mutation type %d\n", time(), indiv_->indiv_id, repl->type());

          switch (repl->type()) {
              case DO_SWITCH:
//                  printf("%d -- %d -- Start switch at %d\n", AeTime::time(),indiv_->indiv_id, repl->pos_1());
                  do_switch(repl->pos_1());
                  nb_swi_++;
                  nb_mut_++;
//        printf("End switch at %d\n",repl->pos_1());
                  break;
              case SMALL_INSERTION:
//                  printf("%d -- %d -- Start insertion at %d (%d %s)\n", AeTime::time(),indiv_->indiv_id, repl->pos_1(), repl->number(), repl->seq());
                  do_small_insertion(repl->pos_1(), repl->number(), repl->seq());

                  nb_indels_++;
                  nb_mut_++;
        //printf("End insertion at %d (%d)\n",repl->pos_1(),repl->number(),repl->seq());
                  break;
              case SMALL_DELETION:
                  //if (indiv_->indiv_id == 626 && AeTime::time() == 21)
//                  printf("%d -- %d -- Start deletion at %d (%d)\n",AeTime::time(),indiv_->indiv_id,  repl->pos_1(), repl->number());
                  do_small_deletion(repl->pos_1(), repl->number());
                  nb_indels_++;
                  nb_mut_++;
//        printf("End deletion at %d (%d)\n",repl->pos_1(),repl->number());
                  break;
              case DUPLICATION:
//                  printf("%d -- %d -- Start duplication at %d (%d %d)\n",AeTime::time(),indiv_->indiv_id,  repl->pos_1(), repl->pos_2(), repl->pos_3());
                  do_duplication(repl->pos_1(), repl->pos_2(), repl->pos_3());
                  nb_large_dupl_++;
                  nb_rear_++;
//        printf("End duplication at %d (%d %d)\n",repl->pos_1(),repl->pos_2(),repl->pos_3());
                  break;
              case TRANSLOCATION:
//                  printf("%d -- %d -- Start translocation at %d (%d %d %d %d)\n",AeTime::time(),indiv_->indiv_id,  repl->pos_1(), repl->pos_2(), repl->pos_3(),
//                         repl->pos_4(), repl->invert());
                  do_translocation(repl->pos_1(), repl->pos_2(), repl->pos_3(),
                                   repl->pos_4(), repl->invert());
                  nb_large_trans_++;
                  nb_rear_++;
//        printf("End translocation at %d (%d %d %d %d)\n",repl->pos_1(),repl->pos_2(),repl->pos_3(),repl->pos_4(),repl->invert());
                  break;
              case INVERSION:
//                  printf("%d -- %d -- Start invertion at %d (%d)\n", AeTime::time(),indiv_->indiv_id, repl->pos_1(), repl->pos_2());
                  do_inversion(repl->pos_1(), repl->pos_2());
                  nb_large_inv_++;
                  nb_rear_++;
//        printf("End invertion at %d (%d)\n",repl->pos_1(),repl->pos_2());
                  break;
              case DELETION:
                  //if (indiv_->indiv_id == 626 && AeTime::time() == 21)
//                  printf("%d -- %d -- Start LARGE deletion at %d (%d)\n", AeTime::time(),indiv_->indiv_id, repl->pos_1(), repl->pos_2());
                  do_deletion(repl->pos_1(), repl->pos_2());
                  nb_large_del_++;
                  nb_rear_++;
//        printf("End LARGE deletion at %d (%d)\n",repl->pos_1(),repl->pos_2());
                  break;
          }

//          if (indiv_->indiv_id == 799 && time() == 52) {
//              printf("%d -- %d -- TOKEEP -- Prom list LEAD : ", time(), indiv_->indiv_id);
//              for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//                  if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                      if (indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                          printf("%d ", indiv_->metadata_->promoters(prom_idx)->pos);
//              }
//              printf("\n");
//              printf("%d -- %d -- TOKEEP -- Prom list LAG : ", time(), indiv_->indiv_id);
//              for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//                  if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                      if (!indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                          printf("%d ", indiv_->metadata_->promoters(prom_idx)->pos);
//              }
//              printf("\n");
//          }
      }

  } while (indiv_->exp_m_->dna_mutator_array_[indiv_->indiv_id]->mutation_available() > 0);
   // printf("End mutation for %d\n",indiv_->indiv_id);
//    if (indiv_->indiv_id == 6) {
//        printf("%d -- %d -- AFTER -- Prom list LEAD : ",time(),indiv_->indiv_id);
//        for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//            if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                if (indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                    printf("%d ",indiv_->metadata_->promoters(prom_idx)->pos);
//        }
//        printf("\n");
//        printf("%d -- %d -- AFTER -- Prom list LAG : ",time(),indiv_->indiv_id);
//        for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//            if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                if (!indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                    printf("%d ",indiv_->metadata_->promoters(prom_idx)->pos);
//        }
//        printf("\n");
//    }


    //if (indiv_->indiv_id==30 || indiv_->indiv_id==61 || indiv_->indiv_id==62) {
//
//        printf("%d -- %d -- AFTER -- TOKEEP -- Prom list LEAD : ", time(), indiv_->indiv_id);
//        for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//            if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                if (indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                    printf("%d ", indiv_->metadata_->promoters(prom_idx)->pos);
//        }
//        printf("\n");
//        printf("%d -- %d -- AFTER -- TOKEEP -- Prom list LAG : ", time(), indiv_->indiv_id);
//        for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//            if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                if (!indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                    printf("%d ", indiv_->metadata_->promoters(prom_idx)->pos);
//        }
//        printf("\n");
//
//        printf("%d -- %d -- AFTER -- ADV-DIRECT -- Prom list LEAD : ", time(), indiv_->indiv_id);
//        for (int prom_idx = 0; prom_idx < ((List_Metadata*)indiv_->metadata_)->promoters_list_[LEADING].size(); prom_idx++) {
//            auto it =  ((List_Metadata*)indiv_->metadata_)->promoters_list_[LEADING].begin();
//            std::advance(it, prom_idx);
//            printf("%d ", (*it).pos);
//        }
//        printf("\n");
//        printf("%d -- %d -- AFTER -- ADV-DIRECT -- Prom list LAG : ", time(), indiv_->indiv_id);
//        for (int prom_idx = 0; prom_idx < ((List_Metadata*)indiv_->metadata_)->promoters_list_[LAGGING].size(); prom_idx++) {
//            auto it =  ((List_Metadata*)indiv_->metadata_)->promoters_list_[LAGGING].begin();
//            std::advance(it, prom_idx);
//            printf("%d ", (*it).pos);
//        }
//        printf("\n");
    //}
}


void Dna_7::apply_mutations() {
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

/*    if (indiv_->indiv_id == 632) {
          printf("%d -- %d -- BEFORE -- Prom list LEAD : ",time(),indiv_->indiv_id);
          for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
              if (indiv_->metadata_->promoters(prom_idx) != nullptr)
                  if (indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
                      printf("%d ",indiv_->metadata_->promoters(prom_idx)->pos);
          }
          printf("\n");
          printf("%d -- %d -- BEFORE -- Prom list LAG : ",time(),indiv_->indiv_id);
          for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
              if (indiv_->metadata_->promoters(prom_idx) != nullptr)
                  if (!indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
                      printf("%d ",indiv_->metadata_->promoters(prom_idx)->pos);
          }
          printf("\n");
    }*/
//    if (indiv_->indiv_id == 6) {
//        printf("%d -- %d -- Parent %d\n",time(),indiv_->indiv_id,indiv_->parent_id);
//        printf("%d -- %d -- BEFORE -- Prom list LEAD : ",time(),indiv_->indiv_id);
//        for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//            if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                if (indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                    printf("%d ",indiv_->metadata_->promoters(prom_idx)->pos);
//        }
//        printf("\n");
//        printf("%d -- %d -- BEFORE -- Prom list LAG : ",time(),indiv_->indiv_id);
//        for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//            if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                if (!indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                    printf("%d ",indiv_->metadata_->promoters(prom_idx)->pos);
//        }
//        printf("\n");
//    }
  for (auto repl : indiv_->exp_m_->dna_mutator_array_[indiv_->indiv_id]->mutation_list_) {
//    if (indiv_->indiv_id == 179) {
//      printf("Mutation type %d\n",repl->type());
//    }
      //printf("Mutation type %d\n",repl->type());
//    if (indiv_->indiv_id == 686) printf("%d -- ",indiv_->indiv_id);
    switch(repl->type()) {
        case DO_SWITCH:
//            printf("Start switch at %d\n",repl->pos_1());
            do_switch(repl->pos_1());
            nb_swi_++;
            nb_mut_++;
//        printf("End switch at %d\n",repl->pos_1());
            break;
        case SMALL_INSERTION:
//           printf("Start insertion at %d (%d %s)\n",repl->pos_1(),repl->number(),repl->seq());
            do_small_insertion(repl->pos_1(), repl->number(), repl->seq());
            nb_indels_++;
            nb_mut_++;
//        printf("End insertion at %d (%d)\n",repl->pos_1(),repl->number(),repl->seq());
            break;
        case SMALL_DELETION:
            //if (indiv_->indiv_id == 626 && AeTime::time() == 21)
//            printf("Start deletion at %d (%d)\n",repl->pos_1(),repl->number());
            do_small_deletion(repl->pos_1(), repl->number());
            nb_indels_++;
            nb_mut_++;
//        printf("End deletion at %d (%d)\n",repl->pos_1(),repl->number());
            break;
        case DUPLICATION:
//            printf("Start duplication at %d (%d %d)\n",repl->pos_1(),repl->pos_2(),repl->pos_3());
            do_duplication(repl->pos_1(), repl->pos_2(), repl->pos_3());
            nb_large_dupl_++;
            nb_rear_++;
//        printf("End duplication at %d (%d %d)\n",repl->pos_1(),repl->pos_2(),repl->pos_3());
            break;
        case TRANSLOCATION:
//            printf("Start translocation at %d (%d %d %d %d)\n",repl->pos_1(),repl->pos_2(),repl->pos_3(),repl->pos_4(),repl->invert());
            do_translocation(repl->pos_1(), repl->pos_2(), repl->pos_3(),
                             repl->pos_4(), repl->invert());
            nb_large_trans_++;
            nb_rear_++;
//        printf("End translocation at %d (%d %d %d %d)\n",repl->pos_1(),repl->pos_2(),repl->pos_3(),repl->pos_4(),repl->invert());
            break;
        case INVERSION:
//            printf("Start invertion at %d (%d)\n",repl->pos_1(),repl->pos_2());
            do_inversion(repl->pos_1(), repl->pos_2());
            nb_large_inv_++;
            nb_rear_++;
//        printf("End invertion at %d (%d)\n",repl->pos_1(),repl->pos_2());
            break;
        case DELETION:
//            if (indiv_->indiv_id == 179)
//            printf("%d -- Start LARGE deletion at %d (%d) %ld\n",indiv_->indiv_id,repl->pos_1(),repl->pos_2(),length());
            do_deletion(repl->pos_1(), repl->pos_2());
            nb_large_del_++;
            nb_rear_++;
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

//    if (indiv_->indiv_id == 6) {
//        printf("%d -- %d -- AFTER -- Prom list LEAD : ",time(),indiv_->indiv_id);
//        for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//            if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                if (indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                    printf("%d ",indiv_->metadata_->promoters(prom_idx)->pos);
//        }
//        printf("\n");
//        printf("%d -- %d -- AFTER -- Prom list LAG : ",time(),indiv_->indiv_id);
//        for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//            if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                if (!indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                    printf("%d ",indiv_->metadata_->promoters(prom_idx)->pos);
//        }
//        printf("\n");
//    }

/*    if (indiv_->indiv_id == 632) {
        printf("%d -- %d -- AFTER -- Prom list LEAD : ",time(),indiv_->indiv_id);
        for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
            if (indiv_->metadata_->promoters(prom_idx) != nullptr)
                if (indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
                    printf("%d ",indiv_->metadata_->promoters(prom_idx)->pos);
        }
        printf("\n");
        printf("%d -- %d -- AFTER -- Prom list LAG : ",time(),indiv_->indiv_id);
        for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
            if (indiv_->metadata_->promoters(prom_idx) != nullptr)
                if (!indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
                    printf("%d ",indiv_->metadata_->promoters(prom_idx)->pos);
        }
        printf("\n");
    }*/

}


void Dna_7::ABCDE_to_ADCBE(int32_t pos_B, int32_t pos_C, int32_t pos_D,
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
  Dna_7 * new_genome = indiv_->dna_factory_->get_dna(length_);
    //printf("BEFORE :: DNA %d || New DNA %d || Block %d || New Block %d\n",length_,new_genome->length_,nb_blocks_,new_genome->nb_blocks_);

    new_genome->set_indiv(length_,parent_length_,indiv_);
    //printf("DNA %d || New DNA %d || Block %d || New Block %d\n",length_,new_genome->length_,nb_blocks_,new_genome->nb_blocks_);
    //printf("LEN %d %d %d\n",len_A,len_D,pos_D);
  memcpy(new_genome->data_, data_, len_A * sizeof(char));
  memcpy(&(new_genome->data_[len_A]), &data_[pos_D], len_D * sizeof(char));
  memcpy(&(new_genome->data_[len_AD]), &data_[pos_C], len_C * sizeof(char));
  memcpy(&(new_genome->data_[len_ADC]), &data_[pos_B], len_B * sizeof(char));
  memcpy(&(new_genome->data_[len_ADCB]), &data_[pos_E], len_E * sizeof(char));
  new_genome->data_[length_] = '\0';

  // Replace sequence
  // NB : The size of the genome doesn't change. Therefore, we don't nee
  // to update length_ and nb_blocks_
    char* old_data = data_;
    data_ = new_genome->data_;
    int new_nb_block = new_genome->nb_blocks_;
    new_genome->nb_blocks_ = nb_blocks_;

    new_genome->data_ = old_data;
    new_genome->length_ = length();

    nb_blocks_ = new_nb_block;
    indiv_->dna_factory_->give_back(new_genome);


  // ========== Update promoter list ==========
  if (length() >= PROM_SIZE) {
//      if (indiv_->indiv_id == 799 && time() == 52) {
//          printf("Before remove XA %d %d %d %d\n",0,pos_B,pos_C,pos_D,pos_E);
//          std::set<int> leading;
//          printf("FROM MUTATE -- %d -- %d -- AFTER -- Prom list LEAD : ",time(),indiv_->indiv_id);
//          for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//              if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                  if (indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                      leading.insert(indiv_->metadata_->promoters(prom_idx)->pos);
//          }
//          for (auto lead : leading) {
//              printf("%d ",lead);
//          }
//
//      printf("BEFORE_REMOVE : Lagging promoters lists : ");
//      for (auto it : indiv_->lagging_prom_pos) {
//        printf("%d (%d) || ", it.first, it.second);
//      }
//      printf("\n");
//
//          printf("\n");
//
//          std::set<int> lagging;
//          printf("FROM MUTATE -- %d -- %d -- AFTER -- Prom list LAG : ",time(),indiv_->indiv_id);
//          for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//              if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                  if (!indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                      lagging.insert(indiv_->metadata_->promoters(prom_idx)->pos);
//          }
//
//          for (auto lag : lagging) {
//              printf("%d ",lag);
//          }
//          printf("\n");
//      }
    // Remove promoters that include a breakpoint
      indiv_->metadata_->remove_promoters_around(pos_B);
      indiv_->metadata_->remove_promoters_around(pos_C);
      indiv_->metadata_->remove_promoters_around(pos_D);
      indiv_->metadata_->remove_promoters_around(pos_E);

//      if (indiv_->indiv_id == 799 && time() == 52) {
//          printf("After remove XA %d %d %d %d\n",0,pos_B,pos_C,pos_D,pos_E);
//          std::set<int> leading;
//          printf("FROM MUTATE -- %d -- %d -- AFTER -- Prom list LEAD : ",time(),indiv_->indiv_id);
//          for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//              if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                  if (indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                      leading.insert(indiv_->metadata_->promoters(prom_idx)->pos);
//          }
//          for (auto lead : leading) {
//              printf("%d ",lead);
//          }
//
//
//          printf("\n");
//
//          std::set<int> lagging;
//          printf("FROM MUTATE -- %d -- %d -- AFTER -- Prom list LAG : ",time(),indiv_->indiv_id);
//          for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//              if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                  if (!indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                      lagging.insert(indiv_->metadata_->promoters(prom_idx)->pos);
//          }
//
//          for (auto lag : lagging) {
//              printf("%d ",lag);
//          }
//          printf("\n");
//      }

    // Create temporary lists for promoters to move and/or invert
    std::vector<std::list<PromoterStruct *>> promoters_B = {{},
                                     {}};
    std::vector<std::list<PromoterStruct *>> promoters_C = {{},
                                     {}};
    std::vector<std::list<PromoterStruct *>> promoters_D = {{},
                                     {}};
    // Extract promoters that are totally included in each segment to be moved
    // and shift them to their new positions
    if (len_B >= PROM_SIZE) {
        indiv_->metadata_->extract_promoters_included_in(pos_B, pos_C, promoters_B);

//        if (indiv_->indiv_id == 799 && time() == 52) {
//            printf("Extracted promoters\n");
//            printf("%d -- %d -- TOKEEP -- Prom list B LEAD : ",time(),indiv_->indiv_id);
//            for (auto prom : promoters_B[LEADING]) {
//                if (prom != nullptr)
//                        printf("%d ",prom->pos);
//            }
//            printf("\n");
//            printf("%d -- %d -- TOKEEP -- Prom list B LAG : ",time(),indiv_->indiv_id);
//            for (auto prom : promoters_B[LAGGING]) {
//                if (prom != nullptr)
//                    printf("%d ",prom->pos);
//            }
//            printf("\n");
//        }



        if (indiv_->exp_m_->exp_s()->get_simd_metadata_flavor() == SIMDMetadataFlavor::STD_MAP)
          Map_Metadata::shift_promoters(promoters_B, len_D + len_C,
                                   length());
        else if (indiv_->exp_m_->exp_s()->get_simd_metadata_flavor() == SIMDMetadataFlavor::DYN_TAB)
          DynTab_Metadata::shift_promoters(promoters_B, len_D + len_C,
                                                   length());
        else if (indiv_->exp_m_->exp_s()->get_simd_metadata_flavor() == SIMDMetadataFlavor::STD_LIST)
          List_Metadata::shift_promoters(promoters_B, len_D + len_C,
                                                  length());
//        if (indiv_->indiv_id == 799 && time() == 52) {
//            printf("After shift Extracted promoters\n");
//            printf("%d -- %d -- TOKEEP -- Prom list B LEAD : ",time(),indiv_->indiv_id);
//            for (auto prom : promoters_B[LEADING]) {
//                if (prom != nullptr)
//                    printf("%d ",prom->pos);
//            }
//            printf("\n");
//            printf("%d -- %d -- TOKEEP -- Prom list B LAG : ",time(),indiv_->indiv_id);
//            for (auto prom : promoters_B[LAGGING]) {
//                if (prom != nullptr)
//                    printf("%d ",prom->pos);
//            }
//            printf("\n");
//        }
    }
    if (len_C >= PROM_SIZE) {
        indiv_->metadata_->extract_promoters_included_in(pos_C, pos_D, promoters_C);

//        if (indiv_->indiv_id == 799 && time() == 52) {
//            printf("Extracted promoters\n");
//            printf("%d -- %d -- TOKEEP -- Prom list C LEAD : ",time(),indiv_->indiv_id);
//            for (auto prom : promoters_C[LEADING]) {
//                if (prom != nullptr)
//                    printf("%d ",prom->pos);
//            }
//            printf("\n");
//            printf("%d -- %d -- TOKEEP -- Prom list C LAG : ",time(),indiv_->indiv_id);
//            for (auto prom : promoters_C[LAGGING]) {
//                if (prom != nullptr)
//                    printf("%d ",prom->pos);
//            }
//            printf("\n");
//        }



        if (indiv_->exp_m_->exp_s()->get_simd_metadata_flavor() == SIMDMetadataFlavor::STD_MAP)
          Map_Metadata::shift_promoters(promoters_C, len_D - len_B,
                                   length());
        else if (indiv_->exp_m_->exp_s()->get_simd_metadata_flavor() == SIMDMetadataFlavor::DYN_TAB)
          DynTab_Metadata::shift_promoters(promoters_C, len_D - len_B,
                                               length());
        else if (indiv_->exp_m_->exp_s()->get_simd_metadata_flavor() == SIMDMetadataFlavor::STD_LIST)
          List_Metadata::shift_promoters(promoters_C, len_D - len_B,
                                                  length());
//        if (indiv_->indiv_id == 799 && time() == 52) {
//            printf("After shift extracted promoters\n");
//            printf("%d -- %d -- TOKEEP -- Prom list C LEAD : ",time(),indiv_->indiv_id);
//            for (auto prom : promoters_C[LEADING]) {
//                if (prom != nullptr)
//                    printf("%d ",prom->pos);
//            }
//            printf("\n");
//            printf("%d -- %d -- TOKEEP -- Prom list C LAG : ",time(),indiv_->indiv_id);
//            for (auto prom : promoters_C[LAGGING]) {
//                if (prom != nullptr)
//                    printf("%d ",prom->pos);
//            }
//            printf("\n");
//        }
    }
    if (len_D >= PROM_SIZE) {
        indiv_->metadata_->extract_promoters_included_in(pos_D, pos_E, promoters_D);

//        if (indiv_->indiv_id == 799 && time() == 52) {
//            printf("Extracted promoters\n");
//            printf("%d -- %d -- TOKEEP -- Prom list D LEAD : ",time(),indiv_->indiv_id);
//            for (auto prom : promoters_D[LEADING]) {
//                if (prom != nullptr)
//                    printf("%d ",prom->pos);
//            }
//            printf("\n");
//            printf("%d -- %d -- TOKEEP -- Prom list D LAG : ",time(),indiv_->indiv_id);
//            for (auto prom : promoters_D[LAGGING]) {
//                if (prom != nullptr)
//                    printf("%d ",prom->pos);
//            }
//            printf("\n");
//        }

        if (indiv_->exp_m_->exp_s()->get_simd_metadata_flavor() == SIMDMetadataFlavor::STD_MAP)
          Map_Metadata::shift_promoters(promoters_D, -len_B - len_C,
                                   length());
        else if (indiv_->exp_m_->exp_s()->get_simd_metadata_flavor() == SIMDMetadataFlavor::DYN_TAB)
          DynTab_Metadata::shift_promoters(promoters_D, -len_B - len_C,
                                               length());
        else if (indiv_->exp_m_->exp_s()->get_simd_metadata_flavor() == SIMDMetadataFlavor::STD_LIST)
          List_Metadata::shift_promoters(promoters_D, -len_B - len_C,
                                                  length());

//        if (indiv_->indiv_id == 799 && time() == 52) {
//            printf("After shift extracted promoters\n");
//            printf("%d -- %d -- TOKEEP -- Prom list D LEAD : ",time(),indiv_->indiv_id);
//            for (auto prom : promoters_D[LEADING]) {
//                if (prom != nullptr)
//                    printf("%d ",prom->pos);
//            }
//            printf("\n");
//            printf("%d -- %d -- TOKEEP -- Prom list D LAG : ",time(),indiv_->indiv_id);
//            for (auto prom : promoters_D[LAGGING]) {
//                if (prom != nullptr)
//                    printf("%d ",prom->pos);
//            }
//            printf("\n");
//        }
    }

//    printf("INsert\n");
//      printf("%d -- %d -- Prom list LEAD : ",time(),indiv_->indiv_id);
//      for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//          if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//              if (indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                  printf("%d ",indiv_->metadata_->promoters(prom_idx)->pos);
//      }
//      printf("\n");
//      printf("%d -- %d -- Prom list LAG : ",time(),indiv_->indiv_id);
//      for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//          if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//              if (!indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                  printf("%d ",indiv_->metadata_->promoters(prom_idx)->pos);
//      }
//      printf("\n");
//
//      printf("%d -- %d -- Prom list promoter_B : ",time(),indiv_->indiv_id);
//      for (auto pb : promoters_B)
//          for (auto prom : pb)
//            if (prom != nullptr)
//                printf("%d ",prom->pos);
//
//      printf("\n");


//        printf("%d -- %d -- Prom list promoter_C : ",time(),indiv_->indiv_id);
//        for (auto pb : promoters_C)
//            for (auto prom : pb)
//                if (prom != nullptr)
//                    printf("%d ",prom->pos);
//      printf("\n");
//
//
//      printf("%d -- %d -- Prom list promoter_D : ",time(),indiv_->indiv_id);
//    for (auto pb : promoters_D)
//        for (auto prom : pb)
//            if (prom != nullptr)
//                printf("%d ",prom->pos);
//
//    printf("\n");

//    }
    // Reinsert the shifted promoters
      indiv_->metadata_->insert_promoters(promoters_B);
      for (auto strand: {LEADING, LAGGING}) {
          for (std::list<PromoterStruct *>::iterator it_rna = promoters_B[strand].begin();
               it_rna != promoters_B[strand].end(); it_rna++) {
              delete (*(it_rna));
          }
      }

      indiv_->metadata_->insert_promoters(promoters_C);
      for (auto strand: {LEADING, LAGGING}) {
          for (std::list<PromoterStruct *>::iterator it_rna = promoters_C[strand].begin();
               it_rna != promoters_C[strand].end(); it_rna++) {
              delete (*(it_rna));
          }
      }

      indiv_->metadata_->insert_promoters(promoters_D);
      for (auto strand: {LEADING, LAGGING}) {
          for (std::list<PromoterStruct *>::iterator it_rna = promoters_D[strand].begin();
               it_rna != promoters_D[strand].end(); it_rna++) {
              delete (*(it_rna));
          }
      }

//      if (indiv_->indiv_id == 799 && time() == 52) {
//          printf("After insert\n");
//          std::set<int> leading;
//          printf("FROM MUTATE -- %d -- %d -- AFTER -- Prom list LEAD : ",time(),indiv_->indiv_id);
//          for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//              if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                  if (indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                      leading.insert(indiv_->metadata_->promoters(prom_idx)->pos);
//          }
//          for (auto lead : leading) {
//              printf("%d ",lead);
//          }
//
//
//          printf("\n");
//
//          std::set<int> lagging;
//          printf("FROM MUTATE -- %d -- %d -- AFTER -- Prom list LAG : ",time(),indiv_->indiv_id);
//          for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//              if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                  if (!indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                      lagging.insert(indiv_->metadata_->promoters(prom_idx)->pos);
//          }
//
//          for (auto lag : lagging) {
//              printf("%d ",lag);
//          }
//          printf("\n");
//      }

//
//      printf("Look for\n");
//      printf("%d -- %d -- Prom list LEAD : ",time(),indiv_->indiv_id);
//      for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//          if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//              if (indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                  printf("%d ",indiv_->metadata_->promoters(prom_idx)->pos);
//      }
//      printf("\n");
//      printf("%d -- %d -- Prom list LAG : ",time(),indiv_->indiv_id);
//      for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//          if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//              if (!indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                  printf("%d ",indiv_->metadata_->promoters(prom_idx)->pos);
//      }
//
//      printf("\n");
//    }
    // 5) Look for new promoters including a breakpoint
      indiv_->metadata_->look_for_new_promoters_around(len_A);
      indiv_->metadata_->look_for_new_promoters_around(len_AD);
      indiv_->metadata_->look_for_new_promoters_around(len_ADC);
      indiv_->metadata_->look_for_new_promoters_around(len_ADCB);

//      if (indiv_->indiv_id == 799 && time() == 52) {
//          printf("After look\n");
//          std::set<int> leading;
//          printf("FROM MUTATE -- %d -- %d -- AFTER -- Prom list LEAD : ",time(),indiv_->indiv_id);
//          for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//              if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                  if (indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                      leading.insert(indiv_->metadata_->promoters(prom_idx)->pos);
//          }
//          for (auto lead : leading) {
//              printf("%d ",lead);
//          }
//
//
//          printf("\n");
//
//          std::set<int> lagging;
//          printf("FROM MUTATE -- %d -- %d -- AFTER -- Prom list LAG : ",time(),indiv_->indiv_id);
//          for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//              if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                  if (!indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                      lagging.insert(indiv_->metadata_->promoters(prom_idx)->pos);
//          }
//
//          for (auto lag : lagging) {
//              printf("%d ",lag);
//          }
//          printf("\n");
//      }
//      printf("Finished\n");
//      printf("%d -- %d -- Prom list LEAD : ",time(),indiv_->indiv_id);
//      for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//          if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//              if (indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                  printf("%d ",indiv_->metadata_->promoters(prom_idx)->pos);
//      }
//      printf("\n");
//      printf("%d -- %d -- Prom list LAG : ",time(),indiv_->indiv_id);
//      for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//          if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//              if (!indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                  printf("%d ",indiv_->metadata_->promoters(prom_idx)->pos);
//      }
//
//      printf("\n");
//    }
  }
}

void Dna_7::ABCDE_to_ADBpCpE(int32_t pos_B, int32_t pos_C, int32_t pos_D,
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

  // Create new sequence
  Dna_7 * new_genome = indiv_->dna_factory_->get_dna(length_);
    new_genome->set_indiv(length_,parent_length_,indiv_);

  // Copy segments A and D
  memcpy(new_genome->data_, data_, len_A * sizeof(char));
  memcpy(&(new_genome->data_[len_A]), &data_[pos_D], len_D * sizeof(char));


  // Build Bp and put it in the new genome
  char* inverted_segment;
  posix_memalign((void **)&inverted_segment,64,(len_B+1)* sizeof(char));//new char[len_B + 1];

//#pragma simd
//#pragma distribute_point
  for (int32_t i = 0, j = pos_C - 1; i < len_B; i++, j--) {
    if (data_[j] == '0') inverted_segment[i] = '1';
    else inverted_segment[i] = '0';
  }
  inverted_segment[len_B] = '\0';

  memcpy(&(new_genome->data_[len_AD]), inverted_segment, len_B * sizeof(char));

  free(inverted_segment);


  // Build Cp and put it in the new genome
  posix_memalign((void **)&inverted_segment,64,(len_C+1)* sizeof(char));//new char[len_C + 1];

//#pragma simd
//#pragma distribute_point
  for (int32_t i = 0, j = pos_D - 1; i < len_C; i++, j--) {
    if (data_[j] == '0') inverted_segment[i] = '1';
    else inverted_segment[i] = '0';
  }
  inverted_segment[len_C] = '\0';

  memcpy(&(new_genome->data_[len_ADB]), inverted_segment, len_C * sizeof(char));

  free(inverted_segment);

  // Copy segment E into the new genome
  memcpy(&(new_genome->data_[len_ADBC]), &data_[pos_E], len_E * sizeof(char));
  new_genome->data_[length_] = '\0';


  // Replace sequence
    char* old_data = data_;
    data_ = new_genome->data_;
    int new_nb_block = new_genome->nb_blocks_;
    new_genome->nb_blocks_ = nb_blocks_;

    new_genome->data_ = old_data;
    new_genome->length_ = length();

    nb_blocks_ = new_nb_block;

    indiv_->dna_factory_->give_back(new_genome);


  // ========== Update promoter list ==========
  if (length() >= PROM_SIZE) {
//    if (indiv_->indiv_id == 18) {
//        printf("BEFORE Remove\n");
//        printf("%d -- %d -- Prom list LEAD : ", time(), indiv_->indiv_id);
//        for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//            if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                if (indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                    printf("%d ", indiv_->metadata_->promoters(prom_idx)->pos);
//        }
//        printf("\n");
//        printf("%d -- %d -- Prom list LAG : ", time(), indiv_->indiv_id);
//        for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//            if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                if (!indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                    printf("%d ", indiv_->metadata_->promoters(prom_idx)->pos);
//        }
//        printf("\n");
//    }
    // Remove promoters that include a breakpoint
      indiv_->metadata_->remove_promoters_around(pos_B);
      indiv_->metadata_->remove_promoters_around(pos_C);
      indiv_->metadata_->remove_promoters_around(pos_D);
      indiv_->metadata_->remove_promoters_around(pos_E);

//    if (indiv_->indiv_id == 18) {
//        printf("After remove\n");
//        printf("%d -- %d -- Prom list LEAD : ", time(), indiv_->indiv_id);
//        for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//            if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                if (indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                    printf("%d ", indiv_->metadata_->promoters(prom_idx)->pos);
//        }
//        printf("\n");
//        printf("%d -- %d -- Prom list LAG : ", time(), indiv_->indiv_id);
//        for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//            if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                if (!indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                    printf("%d ", indiv_->metadata_->promoters(prom_idx)->pos);
//        }
//        printf("\n");
//    }

    // Create temporary lists for promoters to move and/or invert
    std::vector<std::list<PromoterStruct *>> promoters_B = {{},
                                     {}};
    std::vector<std::list<PromoterStruct *>> promoters_C = {{},
                                     {}};
    std::vector<std::list<PromoterStruct *>> promoters_D = {{},
                                     {}};

    // 2) Extract promoters that are totally included in each segment to be
    //    moved (B, C and D)
    if (len_B >= PROM_SIZE) {
        indiv_->metadata_->extract_promoters_included_in(pos_B, pos_C, promoters_B);
    }
    if (len_C >= PROM_SIZE) {
        indiv_->metadata_->extract_promoters_included_in(pos_C, pos_D, promoters_C);
    }
    if (len_D >= PROM_SIZE) {
        indiv_->metadata_->extract_promoters_included_in(pos_D, pos_E, promoters_D);
    }

//            if (indiv_->indiv_id == 18) {
//          printf("AFTER EXTRACT \n");
//          printf("%d -- %d -- Prom list LEAD : ", time(), indiv_->indiv_id);
//          for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//              if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                  if (indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                      printf("%d ", indiv_->metadata_->promoters(prom_idx)->pos);
//          }
//          printf("\n");
//          printf("%d -- %d -- Prom list LAG : ", time(), indiv_->indiv_id);
//          for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//              if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                  if (!indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                      printf("%d ", indiv_->metadata_->promoters(prom_idx)->pos);
//          }
//          printf("\n");
//      printf("%d -- %d -- Prom list promoter_B : ",time(),indiv_->indiv_id);
//      for (auto pb : promoters_B)
//          for (auto prom : pb)
//            if (prom != nullptr)
//                printf("%d ",prom->pos);
//
//      printf("\n");
//
//
//        printf("%d -- %d -- Prom list promoter_C : ",time(),indiv_->indiv_id);
//        for (auto pb : promoters_C)
//            for (auto prom : pb)
//                if (prom != nullptr)
//                    printf("%d ",prom->pos);
//      printf("\n");
//
//
//      printf("%d -- %d -- Prom list promoter_D : ",time(),indiv_->indiv_id);
//    for (auto pb : promoters_D)
//        for (auto prom : pb)
//            if (prom != nullptr)
//                printf("%d ",prom->pos);
//
//    printf("\n");
//
//
//      }

//      printf("AFTER_EXTRACT : Lagging promoters lists : ");
//      for (auto it : indiv_->lagging_prom_pos) {
//        printf("%d (%d) || ", it.first, it.second);
//      }
//      printf("\n");

    // 3a) Invert promoters of segments B and C
      if (indiv_->exp_m_->exp_s()->get_simd_metadata_flavor() == SIMDMetadataFlavor::STD_MAP) {
        Map_Metadata::invert_promoters(promoters_B, pos_B, pos_C);

        Map_Metadata::invert_promoters(promoters_C, pos_C, pos_D);
      } else  if (indiv_->exp_m_->exp_s()->get_simd_metadata_flavor() == SIMDMetadataFlavor::DYN_TAB) {
        DynTab_Metadata::invert_promoters(promoters_B, pos_B, pos_C);

        DynTab_Metadata::invert_promoters(promoters_C, pos_C, pos_D);
      } else  if (indiv_->exp_m_->exp_s()->get_simd_metadata_flavor() == SIMDMetadataFlavor::STD_LIST) {
        List_Metadata::invert_promoters(promoters_B, pos_B, pos_C);

        List_Metadata::invert_promoters(promoters_C, pos_C, pos_D);
      }

//      if (indiv_->indiv_id == 18) {
//          printf("AFTER INVERT \n");
//          printf("%d -- %d -- Prom list LEAD : ", time(), indiv_->indiv_id);
//          for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//              if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                  if (indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                      printf("%d ", indiv_->metadata_->promoters(prom_idx)->pos);
//          }
//          printf("\n");
//          printf("%d -- %d -- Prom list LAG : ", time(), indiv_->indiv_id);
//          for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//              if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                  if (!indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                      printf("%d ", indiv_->metadata_->promoters(prom_idx)->pos);
//          }
//          printf("\n");
//          printf("%d -- %d -- Prom list promoter_B : ",time(),indiv_->indiv_id);
//          for (auto pb : promoters_B)
//              for (auto prom : pb)
//                  if (prom != nullptr)
//                      printf("%d ",prom->pos);
//
//          printf("\n");
//
//
//          printf("%d -- %d -- Prom list promoter_C : ",time(),indiv_->indiv_id);
//          for (auto pb : promoters_C)
//              for (auto prom : pb)
//                  if (prom != nullptr)
//                      printf("%d ",prom->pos);
//          printf("\n");
//
//
//          printf("%d -- %d -- Prom list promoter_D : ",time(),indiv_->indiv_id);
//          for (auto pb : promoters_D)
//              for (auto prom : pb)
//                  if (prom != nullptr)
//                      printf("%d ",prom->pos);
//
//          printf("\n");
//
//
//      }

    // 3b) Shift these promoters positions
      if (indiv_->exp_m_->exp_s()->get_simd_metadata_flavor() == SIMDMetadataFlavor::STD_MAP) {
        Map_Metadata::shift_promoters(promoters_B, len_D,
                                               length());

        Map_Metadata::shift_promoters(promoters_C, len_D,
                                                length());

        Map_Metadata::shift_promoters(promoters_D, -len_B - len_C,
                                                length());
      } else if (indiv_->exp_m_->exp_s()->get_simd_metadata_flavor() == SIMDMetadataFlavor::DYN_TAB) {
        DynTab_Metadata::shift_promoters(promoters_B, len_D,
                                             length());

        DynTab_Metadata::shift_promoters(promoters_C, len_D,
                                             length());

        DynTab_Metadata::shift_promoters(promoters_D, -len_B - len_C,
                                             length());
      } else if (indiv_->exp_m_->exp_s()->get_simd_metadata_flavor() == SIMDMetadataFlavor::STD_LIST) {
        List_Metadata::shift_promoters(promoters_B, len_D,
                                                length());

        List_Metadata::shift_promoters(promoters_C, len_D,
                                                length());

        List_Metadata::shift_promoters(promoters_D, -len_B - len_C,
                                                length());
      }

//      if (indiv_->indiv_id == 18) {
//          printf("AFTER SHIFT \n");
//          printf("%d -- %d -- Prom list LEAD : ", time(), indiv_->indiv_id);
//          for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//              if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                  if (indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                      printf("%d ", indiv_->metadata_->promoters(prom_idx)->pos);
//          }
//          printf("\n");
//          printf("%d -- %d -- Prom list LAG : ", time(), indiv_->indiv_id);
//          for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//              if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                  if (!indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                      printf("%d ", indiv_->metadata_->promoters(prom_idx)->pos);
//          }
//          printf("\n");
//          printf("%d -- %d -- Prom list promoter_B : ",time(),indiv_->indiv_id);
//          for (auto pb : promoters_B)
//              for (auto prom : pb)
//                  if (prom != nullptr)
//                      printf("%d ",prom->pos);
//
//          printf("\n");
//
//
//          printf("%d -- %d -- Prom list promoter_C : ",time(),indiv_->indiv_id);
//          for (auto pb : promoters_C)
//              for (auto prom : pb)
//                  if (prom != nullptr)
//                      printf("%d ",prom->pos);
//          printf("\n");
//
//
//          printf("%d -- %d -- Prom list promoter_D : ",time(),indiv_->indiv_id);
//          for (auto pb : promoters_D)
//              for (auto prom : pb)
//                  if (prom != nullptr)
//                      printf("%d ",prom->pos);
//
//          printf("\n");
//
//
//      }

    // 4) Reinsert the shifted promoters
      indiv_->metadata_->insert_promoters(promoters_C);
      for (auto strand: {LEADING, LAGGING}) {
          for (std::list<PromoterStruct *>::iterator it_rna = promoters_C[strand].begin();
               it_rna != promoters_C[strand].end(); it_rna++) {
              delete (*(it_rna));
          }
      }

      indiv_->metadata_->insert_promoters(promoters_B);
      for (auto strand: {LEADING, LAGGING}) {
          for (std::list<PromoterStruct *>::iterator it_rna = promoters_B[strand].begin();
               it_rna != promoters_B[strand].end(); it_rna++) {
              delete (*(it_rna));
          }
      }

      indiv_->metadata_->insert_promoters(promoters_D);
      for (auto strand: {LEADING, LAGGING}) {
          for (std::list<PromoterStruct *>::iterator it_rna = promoters_D[strand].begin();
               it_rna != promoters_D[strand].end(); it_rna++) {
              delete (*(it_rna));
          }
      }


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
      indiv_->metadata_->look_for_new_promoters_around(len_A);
      indiv_->metadata_->look_for_new_promoters_around(len_AD);
      indiv_->metadata_->look_for_new_promoters_around(len_ADB);
      indiv_->metadata_->look_for_new_promoters_around(len_ADBC);
  }
}

void Dna_7::ABCDE_to_ACpDpBE(int32_t pos_B, int32_t pos_C, int32_t pos_D,
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

  // Create new sequence
  Dna_7 * new_genome = indiv_->dna_factory_->get_dna(length_);
    new_genome->set_indiv(length_,parent_length_,indiv_);

  // Copy segment A
  memcpy(new_genome->data_, data_, len_A * sizeof(char));


  // Build Cp and put it in the new genome
  char* inverted_segment;
  posix_memalign((void **)&inverted_segment,64,(len_C+1)* sizeof(char));//new char[len_C + 1];

  for (int32_t i = 0, j = pos_D - 1; i < len_C; i++, j--) {
    if (data_[j] == '0') inverted_segment[i] = '1';
    else inverted_segment[i] = '0';
  }
  inverted_segment[len_C] = '\0';

  memcpy(&(new_genome->data_[len_A]), inverted_segment, len_C * sizeof(char));

  free(inverted_segment);


  // Build Dp and put it in the new genome
  posix_memalign((void **)&inverted_segment,64,(len_D+1)* sizeof(char));//new char[len_D + 1];

  for (int32_t i = 0, j = pos_E - 1; i < len_D; i++, j--) {
    if (data_[j] == '0') inverted_segment[i] = '1';
    else inverted_segment[i] = '0';
  }
  inverted_segment[len_D] = '\0';

  memcpy(&(new_genome->data_[len_AC]), inverted_segment, len_D * sizeof(char));

  free(inverted_segment);

  // Copy segments B and E
  memcpy(&(new_genome->data_[len_ACD]), &data_[pos_B], len_B * sizeof(char));
  memcpy(&(new_genome->data_[len_ACDB]), &data_[pos_E], len_E * sizeof(char));
    new_genome->data_[length_] = '\0';


  // Replace sequence
    char* old_data = data_;
    data_ = new_genome->data_;

    int new_nb_block = new_genome->nb_blocks_;
    new_genome->nb_blocks_ = nb_blocks_;

    new_genome->data_ = old_data;
    new_genome->length_ = length();

    nb_blocks_ = new_nb_block;
    indiv_->dna_factory_->give_back(new_genome);


  // ========== Update promoter list ==========
  // 1) Remove promoters that include a breakpoint
  // 2) Extract promoters that are totally included in each segment to be
  //    moved (B, C and D)
  // 3) Shift (and invert when needed) these promoters positions
  // 4) Reinsert the shifted promoters
  // 5) Look for new promoters including a breakpoint
  if (length() >= PROM_SIZE) {

//    if (indiv_->indiv_id == 17) {
//        printf("Remove\n");
//        printf("%d -- %d -- Prom list LEAD : ", time(), indiv_->indiv_id);
//        for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//            if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                if (indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                    printf("%d ", indiv_->metadata_->promoters(prom_idx)->pos);
//        }
//        printf("\n");
//        printf("%d -- %d -- Prom list LAG : ", time(), indiv_->indiv_id);
//        for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//            if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                if (!indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                    printf("%d ", indiv_->metadata_->promoters(prom_idx)->pos);
//        }
//        printf("\n");
//    }

    // 1) Remove promoters that include a breakpoint
      indiv_->metadata_->remove_promoters_around(pos_B);
      indiv_->metadata_->remove_promoters_around(pos_C);
      indiv_->metadata_->remove_promoters_around(pos_D);
      indiv_->metadata_->remove_promoters_around(pos_E);

//      if (indiv_->indiv_id == 17) {
//          printf("AFTER Remove\n");
//          printf("%d -- %d -- Prom list LEAD : ", time(), indiv_->indiv_id);
//          for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//              if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                  if (indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                      printf("%d ", indiv_->metadata_->promoters(prom_idx)->pos);
//          }
//          printf("\n");
//          printf("%d -- %d -- Prom list LAG : ", time(), indiv_->indiv_id);
//          for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//              if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                  if (!indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                      printf("%d ", indiv_->metadata_->promoters(prom_idx)->pos);
//          }
//          printf("\n");
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
    std::vector<std::list<PromoterStruct *>> promoters_B = {{},
                                     {}};
    std::vector<std::list<PromoterStruct *>> promoters_C = {{},
                                     {}};
    std::vector<std::list<PromoterStruct *>> promoters_D = {{},
                                     {}};

    // 2) Extract promoters that are totally included in each segment to be
    //    moved (B, C and D)
    if (len_B >= PROM_SIZE) {
        indiv_->metadata_->extract_promoters_included_in(pos_B, pos_C, promoters_B);
    }
    if (len_C >= PROM_SIZE) {
        indiv_->metadata_->extract_promoters_included_in(pos_C, pos_D, promoters_C);
    }
    if (len_D >= PROM_SIZE) {
        indiv_->metadata_->extract_promoters_included_in(pos_D, pos_E, promoters_D);
    }
//      if (indiv_->indiv_id == 17) {
//          printf("AFTER EXTRACT \n");
//          printf("%d -- %d -- Prom list LEAD : ", time(), indiv_->indiv_id);
//          for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//              if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                  if (indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                      printf("%d ", indiv_->metadata_->promoters(prom_idx)->pos);
//          }
//          printf("\n");
//          printf("%d -- %d -- Prom list LAG : ", time(), indiv_->indiv_id);
//          for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//              if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                  if (!indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                      printf("%d ", indiv_->metadata_->promoters(prom_idx)->pos);
//          }
//          printf("\n");
//      printf("%d -- %d -- Prom list promoter_B : ",time(),indiv_->indiv_id);
//      for (auto pb : promoters_B)
//          for (auto prom : pb)
//            if (prom != nullptr)
//                printf("%d ",prom->pos);
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
//        printf("%d -- %d -- Prom list promoter_C : ",time(),indiv_->indiv_id);
//        for (auto pb : promoters_C)
//            for (auto prom : pb)
//                if (prom != nullptr)
//                    printf("%d ",prom->pos);
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
//      printf("%d -- %d -- Prom list promoter_D : ",time(),indiv_->indiv_id);
//    for (auto pb : promoters_D)
//        for (auto prom : pb)
//            if (prom != nullptr)
//                printf("%d ",prom->pos);
//
//    printf("\n");
//
//
//      }
//      printf("\n");
//    }

    // 3a) Invert promoters of segments C and D
      if (indiv_->exp_m_->exp_s()->get_simd_metadata_flavor() == SIMDMetadataFlavor::STD_MAP) {
        Map_Metadata::invert_promoters(promoters_C, pos_C, pos_D);

        Map_Metadata::invert_promoters(promoters_D, pos_D, pos_E);
      } else if (indiv_->exp_m_->exp_s()->get_simd_metadata_flavor() == SIMDMetadataFlavor::DYN_TAB) {
        DynTab_Metadata::invert_promoters(promoters_C, pos_C, pos_D);

        DynTab_Metadata::invert_promoters(promoters_D, pos_D, pos_E);
      } else if (indiv_->exp_m_->exp_s()->get_simd_metadata_flavor() == SIMDMetadataFlavor::STD_LIST) {
        List_Metadata::invert_promoters(promoters_C, pos_C, pos_D);

        List_Metadata::invert_promoters(promoters_D, pos_D, pos_E);
      }

//      if (indiv_->indiv_id == 17) {
//          printf("AFTER Invert\n");
//          printf("%d -- %d -- Prom list LEAD : ", time(), indiv_->indiv_id);
//          for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//              if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                  if (indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                      printf("%d ", indiv_->metadata_->promoters(prom_idx)->pos);
//          }
//          printf("\n");
//          printf("%d -- %d -- Prom list LAG : ", time(), indiv_->indiv_id);
//          for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//              if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                  if (!indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                      printf("%d ", indiv_->metadata_->promoters(prom_idx)->pos);
//          }
//          printf("\n");
//
//          printf("%d -- %d -- Prom list promoter_B : ",time(),indiv_->indiv_id);
//          for (auto pb : promoters_B)
//              for (auto prom : pb)
//                  if (prom != nullptr)
//                      printf("%d ",prom->pos);
//
//          printf("\n");
//
//
//          printf("%d -- %d -- Prom list promoter_C : ",time(),indiv_->indiv_id);
//          for (auto pb : promoters_C)
//              for (auto prom : pb)
//                  if (prom != nullptr)
//                      printf("%d ",prom->pos);
//          printf("\n");
//
//
//          printf("%d -- %d -- Prom list promoter_D : ",time(),indiv_->indiv_id);
//          for (auto pb : promoters_D)
//              for (auto prom : pb)
//                  if (prom != nullptr)
//                      printf("%d ",prom->pos);
//
//          printf("\n");
//      }
//      printf("\n");
//    }

    // 3b) Shift these promoters positions
      if (indiv_->exp_m_->exp_s()->get_simd_metadata_flavor() == SIMDMetadataFlavor::STD_MAP) {
        Map_Metadata::shift_promoters(promoters_B, len_C + len_D,
                                                length());

        Map_Metadata::shift_promoters(promoters_C, -len_B,
                                                length());

        Map_Metadata::shift_promoters(promoters_D, -len_B,
                                                length());
      } else if (indiv_->exp_m_->exp_s()->get_simd_metadata_flavor() == SIMDMetadataFlavor::DYN_TAB) {
        DynTab_Metadata::shift_promoters(promoters_B, len_C + len_D,
                                             length());

        DynTab_Metadata::shift_promoters(promoters_C, -len_B,
                                             length());

        DynTab_Metadata::shift_promoters(promoters_D, -len_B,
                                             length());
      } else if (indiv_->exp_m_->exp_s()->get_simd_metadata_flavor() == SIMDMetadataFlavor::STD_LIST) {
        List_Metadata::shift_promoters(promoters_B, len_C + len_D,
                                                length());

        List_Metadata::shift_promoters(promoters_C, -len_B,
                                                length());

        List_Metadata::shift_promoters(promoters_D, -len_B,
                                                length());
      }

//      if (indiv_->indiv_id == 17) {
//          printf("AFTER Shift\n");
//          printf("%d -- %d -- Prom list LEAD : ", time(), indiv_->indiv_id);
//          for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//              if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                  if (indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                      printf("%d ", indiv_->metadata_->promoters(prom_idx)->pos);
//          }
//          printf("\n");
//          printf("%d -- %d -- Prom list LAG : ", time(), indiv_->indiv_id);
//          for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//              if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                  if (!indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                      printf("%d ", indiv_->metadata_->promoters(prom_idx)->pos);
//          }
//          printf("\n");
//
//          printf("%d -- %d -- Prom list promoter_B : ",time(),indiv_->indiv_id);
//          for (auto pb : promoters_B)
//              for (auto prom : pb)
//                  if (prom != nullptr)
//                      printf("%d ",prom->pos);
//
//          printf("\n");
//
//
//          printf("%d -- %d -- Prom list promoter_C : ",time(),indiv_->indiv_id);
//          for (auto pb : promoters_C)
//              for (auto prom : pb)
//                  if (prom != nullptr)
//                      printf("%d ",prom->pos);
//          printf("\n");
//
//
//          printf("%d -- %d -- Prom list promoter_D : ",time(),indiv_->indiv_id);
//          for (auto pb : promoters_D)
//              for (auto prom : pb)
//                  if (prom != nullptr)
//                      printf("%d ",prom->pos);
//
//          printf("\n");
//      }

    // 4) Reinsert the shifted promoters
      indiv_->metadata_->insert_promoters(promoters_B);
      for (auto strand: {LEADING, LAGGING}) {
          for (std::list<PromoterStruct *>::iterator it_rna = promoters_B[strand].begin();
               it_rna != promoters_B[strand].end(); it_rna++) {
              delete (*(it_rna));
          }
      }

      indiv_->metadata_->insert_promoters(promoters_D);
      for (auto strand: {LEADING, LAGGING}) {
          for (std::list<PromoterStruct *>::iterator it_rna = promoters_D[strand].begin();
               it_rna != promoters_D[strand].end(); it_rna++) {
              delete (*(it_rna));
          }
      }

      indiv_->metadata_->insert_promoters(promoters_C);
      for (auto strand: {LEADING, LAGGING}) {
          for (std::list<PromoterStruct *>::iterator it_rna = promoters_C[strand].begin();
               it_rna != promoters_C[strand].end(); it_rna++) {
              delete (*(it_rna));
          }
      }

//      if (indiv_->indiv_id == 17) {
//          printf("AFTER Insert\n");
//          printf("%d -- %d -- Prom list LEAD : ", time(), indiv_->indiv_id);
//          for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//              if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                  if (indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                      printf("%d ", indiv_->metadata_->promoters(prom_idx)->pos);
//          }
//          printf("\n");
//          printf("%d -- %d -- Prom list LAG : ", time(), indiv_->indiv_id);
//          for (int prom_idx = 0; prom_idx < indiv_->metadata_->promoter_count(); prom_idx++) {
//              if (indiv_->metadata_->promoters(prom_idx) != nullptr)
//                  if (!indiv_->metadata_->promoters(prom_idx)->leading_or_lagging)
//                      printf("%d ", indiv_->metadata_->promoters(prom_idx)->pos);
//          }
//          printf("\n");
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
      indiv_->metadata_->look_for_new_promoters_around(len_A);
      indiv_->metadata_->look_for_new_promoters_around(len_AC);
      indiv_->metadata_->look_for_new_promoters_around(len_ACD);
      indiv_->metadata_->look_for_new_promoters_around(len_ACDB);
  }
}



}