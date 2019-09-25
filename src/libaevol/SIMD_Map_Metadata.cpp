//
// Created by arrouan on 18/07/19.
//

#include "SIMD_Map_Metadata.h"
#include "AeTime.h"

namespace aevol {
    int8_t SIMD_Map_Metadata::is_promoter_leading(int pos) {
        int8_t prom_dist_leading[26];
        int len = indiv_->dna_->length();

        for (int motif_id = 0; motif_id < 22; motif_id++) {
            prom_dist_leading[motif_id] =
                    PROM_SEQ_LEAD[motif_id] ==
                            indiv_->dna_->data_[pos + motif_id >= len ? pos + motif_id - len : pos + motif_id]
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

    int8_t SIMD_Map_Metadata::is_promoter_lagging(int pos) {
        int8_t prom_dist[26];
        int len = indiv_->dna_->length();

        for (int motif_id = 0; motif_id < 22; motif_id++) {
            prom_dist[motif_id] =
                    PROM_SEQ_LAG[motif_id] ==
                            indiv_->dna_->data_[pos - motif_id < 0 ? len + pos - motif_id : pos - motif_id]
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

    void SIMD_Map_Metadata::lst_promoters(bool lorl, Position before_after_btw, int32_t pos1, int32_t pos2,
                                          std::list<promoterStruct*>& motif_list) {
        auto it_begin = lagging_prom_pos_.begin();
        auto it_end = lagging_prom_pos_.end();


        if (lorl == LEADING) {
            it_begin = leading_prom_pos_.begin();
            it_end = leading_prom_pos_.end();
        }

        if (before_after_btw != BEFORE && pos1 != -1) {
            if (lorl == LEADING) {
                auto tmp_it = leading_prom_pos_.lower_bound(pos1);
                if (tmp_it == leading_prom_pos_.end())
                    return;
                if (tmp_it!=leading_prom_pos_.end()) it_begin = tmp_it;
            } else {
                auto tmp_it = lagging_prom_pos_.lower_bound(pos1);
                it_end = tmp_it;
            }
        }

        if (before_after_btw != AFTER && pos2 != -1) {
            if (lorl == LEADING) {
                auto tmp_it = leading_prom_pos_.lower_bound(pos2);
                if (tmp_it!=leading_prom_pos_.end()) it_end = tmp_it;
                if (tmp_it!=leading_prom_pos_.end()) it_end = tmp_it;
            } else {
                auto tmp_it = lagging_prom_pos_.lower_bound(pos2);
                it_begin = tmp_it;
            }
        }

        for (auto it = it_begin; it!=it_end; it++) {
            motif_list.push_back(promoters_[it->second]);
        }
    }

    void SIMD_Map_Metadata::remove_promoters_around(int32_t pos) {
        if (length() >= PROM_SIZE) {
            remove_leading_promoters_starting_between(Utils::mod(pos - PROM_SIZE + 1,
                                                                 length()),
                                                      pos);
            remove_lagging_promoters_starting_between(pos,
                                                      Utils::mod(pos + PROM_SIZE - 1,
                                                                 length()));
        }
        else {
            remove_all_promoters();
        }
    }

    void SIMD_Map_Metadata::remove_promoters_around(int32_t pos_1, int32_t pos_2) {
        if (Utils::mod(pos_1 - pos_2, length()) >= PROM_SIZE) {
//            printf("Remove LEADING between %d %d\n",Utils::mod(pos_1 - PROM_SIZE + 1,
//                                                               length()),
//                   pos_2);
            remove_leading_promoters_starting_between(Utils::mod(pos_1 - PROM_SIZE + 1,
                                                                 length()),
                                                      pos_2);
//            printf("Remove LAGGING between %d %d\n",pos_1,
//                   Utils::mod(pos_2 + PROM_SIZE - 1,
//                              length()));
            remove_lagging_promoters_starting_between(pos_1,
                                                      Utils::mod(pos_2 + PROM_SIZE - 1,
                                                                 length()));
        }
        else {
            remove_all_promoters();
        }
    }

    void SIMD_Map_Metadata::remove_all_promoters() {
        leading_prom_pos_.clear();
        lagging_prom_pos_.clear();

        for (auto it = promoters_.begin(),
                     nextit = it;
             it != promoters_.end();
             it = nextit) {
            delete it->second;
            nextit = next(it);
            promoters_.erase(it);
        }

        promoters_.clear();
        count_promoters_ = 0;
    }

    void SIMD_Map_Metadata::look_for_new_promoters_around(int32_t pos_1, int32_t pos_2) {
        if (length() >= PROM_SIZE) {
            look_for_new_leading_promoters_starting_between(
                    Utils::mod(pos_1 - PROM_SIZE + 1,
                               length()), pos_2);

            look_for_new_lagging_promoters_starting_between(pos_1, Utils::mod(
                    pos_2 + PROM_SIZE - 1,
                    length()));
        }
    }

    void SIMD_Map_Metadata::look_for_new_promoters_around(int32_t pos) {
        if (length() >= PROM_SIZE) {
            //printf("%d -- %d -- LNPA-1 -- Number of RNAs %d (%d)\n",AeTime::time(),indiv_->indiv_id,indiv_->metadata_->rna_count(),
            //       indiv_->metadata_->promoter_count());
            look_for_new_leading_promoters_starting_between(
                    Utils::mod(pos - PROM_SIZE + 1, length()),
                    pos);
            //printf("%d -- %d -- LNPA-2 -- Number of RNAs %d (%d)\n",AeTime::time(),indiv_->indiv_id,indiv_->metadata_->rna_count(),
            //       indiv_->metadata_->promoter_count());
            look_for_new_lagging_promoters_starting_between(
                    pos,
                    Utils::mod(pos + PROM_SIZE - 1, length()));
            //printf("%d -- %d -- LNPA-3 -- Number of RNAs %d (%d)\n",AeTime::time(),indiv_->indiv_id,indiv_->metadata_->rna_count(),
            //       indiv_->metadata_->promoter_count());
        }
    }

    void SIMD_Map_Metadata::locate_promoters() {
        look_for_new_leading_promoters_starting_between(0,length());
        look_for_new_lagging_promoters_starting_between(0,length());
    }

    void SIMD_Map_Metadata::move_all_promoters_after(int32_t pos, int32_t delta_pos) {
        move_all_leading_promoters_after(pos, delta_pos);
        move_all_lagging_promoters_after(pos, delta_pos);
    }

    void SIMD_Map_Metadata::duplicate_promoters_included_in(int32_t pos_1, int32_t pos_2,
                                                            std::vector<std::list<promoterStruct *>> &duplicated_promoters) {
        // 1) Get promoters to be duplicated
        std::vector<std::list<promoterStruct*>> retrieved_promoters = {{},
                                                                       {}};

        promoters_included_in(pos_1, pos_2, retrieved_promoters);

//        printf("RETRIEVED PROMs LEAD : ");
//        for (auto prom :retrieved_promoters[LEADING]) {
//            if (prom != nullptr)
//                if (prom->leading_or_lagging)
//                    printf("%d ",prom->pos);
//        }
//        printf("\n");
//        printf("RETRIEVED PROMs LAG : ");
//        for (auto prom :retrieved_promoters[LAGGING]) {
//            if (prom != nullptr)
//                if (!prom->leading_or_lagging)
//                    printf("%d ",prom->pos);
//        }
//        printf("\n");

        // 2) Set RNAs' position as their position on the duplicated segment
        for (auto& strand: {LEADING, LAGGING}) {
            for (auto& prom : retrieved_promoters[strand]) {
                // Make a copy of current RNA inside container
                duplicated_promoters[strand].push_back(new promoterStruct(prom));

                // Set RNA's position as it's position on the duplicated segment
                duplicated_promoters[strand].back()->pos = Utils::mod(duplicated_promoters[strand].back()->pos -pos_1,
                                                                      length());
            }
        }

    }

    void SIMD_Map_Metadata::extract_promoters_included_in(int32_t pos_1, int32_t pos_2,
                                                          std::vector<std::list<promoterStruct *>> &extracted_promoters) {
        if (pos_2 - pos_1 < PROM_SIZE) {
            return;
        }

        extract_leading_promoters_starting_between(pos_1, pos_2 - PROM_SIZE + 1,
                                                   extracted_promoters[LEADING]);
        extract_lagging_promoters_starting_between(pos_1 + PROM_SIZE - 1, pos_2,
                                                   extracted_promoters[LAGGING]);
    }

    void SIMD_Map_Metadata::insert_promoters(std::vector<std::list<promoterStruct *>> &promoters_to_insert) {
        for (auto strand: {LEADING, LAGGING}) {
            if (promoters_to_insert[strand].size() <= 0) {
                continue;
            }
            // Insert the promoters in the individual's RNA list
            for (auto& to_insert: promoters_to_insert[strand]) {
                if (strand == LEADING) {
                    if (leading_prom_pos_.find(to_insert->pos) == leading_prom_pos_.end()) {

                        int prom_idx;
#pragma omp atomic capture
                        {
                            prom_idx = count_promoters_;
                            count_promoters_ = count_promoters_ + 1;
                        }

                        promoters_[prom_idx] = to_insert;
                        leading_prom_pos_[to_insert->pos] = prom_idx;
                    }
                } else {
                    if (lagging_prom_pos_.find(to_insert->pos) == lagging_prom_pos_.end()) {
                        int prom_idx;
#pragma omp atomic capture
                        {
                            prom_idx = count_promoters_;
                            count_promoters_ = count_promoters_ + 1;
                        }

                        promoters_[prom_idx] = to_insert;

                        lagging_prom_pos_[to_insert->pos] = prom_idx;
                    }
                }
            }
        }
    }

    void SIMD_Map_Metadata::insert_promoters_at(std::vector<std::list<promoterStruct*>>& promoters_to_insert,
                             int32_t pos) {
        for (auto strand: {LEADING, LAGGING}) {
            if (promoters_to_insert[strand].size() <= 0) {
                continue;
            }

            // Insert the promoters in the individual's RNA list
            for (auto &to_insert: promoters_to_insert[strand]) {
                // Update promoter position
                to_insert->pos = Utils::mod(to_insert->pos + pos, length());
                if (strand == LEADING) {
                    if (leading_prom_pos_.find(to_insert->pos) == leading_prom_pos_.end()) {

                        int prom_idx;
#pragma omp atomic capture
                        {
                            prom_idx = count_promoters_;
                            count_promoters_ = count_promoters_ + 1;
                        }

                        promoters_[prom_idx] = to_insert;

                        leading_prom_pos_[to_insert->pos] = prom_idx;

                    }
                } else {
                    if (lagging_prom_pos_.find(to_insert->pos) == lagging_prom_pos_.end()) {
                        int prom_idx;
#pragma omp atomic capture
                        {
                            prom_idx = count_promoters_;
                            count_promoters_ = count_promoters_ + 1;
                        }

                        promoters_[prom_idx] = to_insert;
                        lagging_prom_pos_[to_insert->pos] = prom_idx;
                    }
                }

            }
        }
    }

    void SIMD_Map_Metadata::invert_promoters_included_in(int32_t pos1,
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
        SIMD_Map_Metadata::invert_promoters(inverted_promoters, pos1, pos2);

        // 3) Reinsert the inverted promoters
        insert_promoters(inverted_promoters);
    }

    void SIMD_Map_Metadata::shift_promoters(
            std::vector<std::list<promoterStruct*>>& promoters_to_shift,
            int32_t delta_pos,
            int32_t seq_length) {
        for (auto& strand: {LEADING, LAGGING})
            for (auto& prom: promoters_to_shift[strand])
                prom->pos = Utils::mod(prom->pos + delta_pos, seq_length);
    }

    void SIMD_Map_Metadata::invert_promoters(std::vector<std::list<promoterStruct*>>& promoter_lists,
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

    void SIMD_Map_Metadata::remove_leading_promoters_starting_between(int32_t pos_1,
                                                   int32_t pos_2) {
        if (pos_1 > pos_2) {
            remove_leading_promoters_starting_after(pos_1);
            remove_leading_promoters_starting_before(pos_2);
        }
        else {
            // STL Warning: don't erase the current iterator in the for-loop!
            for (auto it = leading_prom_pos_.lower_bound(pos_1),
                         nextit = it;
                 it != leading_prom_pos_.end() and it->first < pos_2;
                 it = nextit) {

                int pidx = it->second;
                auto it_p = promoters_[pidx];

                delete it_p;
                promoters_.erase(pidx);
                nextit = next(it);
                leading_prom_pos_.erase(it);
            }
        }

    }
    void SIMD_Map_Metadata::remove_leading_promoters_starting_after(int32_t pos) {
        auto init_it = leading_prom_pos_.lower_bound(pos);
        if (init_it == leading_prom_pos_.end())
            return;


        for (auto it = init_it,
                     nextit = it;
             it != leading_prom_pos_.end();
             it = nextit) {

            delete promoters_[it->second];
            promoters_.erase(it->second);
            nextit = next(it);
            leading_prom_pos_.erase(it);
        }
    }

    void SIMD_Map_Metadata::remove_leading_promoters_starting_before(int32_t pos) {
        // Delete RNAs until we reach pos (or we reach the end of the list)
        for (auto it = leading_prom_pos_.begin(),
                     nextit = it;
             it != leading_prom_pos_.end() and it->first < pos;
             it = nextit) {
            delete promoters_[it->second];
            promoters_.erase(it->second);
            nextit = next(it);
            leading_prom_pos_.erase(it);
        }
    }

    void SIMD_Map_Metadata::remove_lagging_promoters_starting_between(int32_t pos_1,
                                                   int32_t pos_2) {
        if (pos_1 == length()) pos_1 = 0;
        if (pos_2 == 0) pos_2 = length();

        if (pos_1 >
            pos_2) { // vld: that's a weird case... really do this? used from remove_promoters_around()
            remove_lagging_promoters_starting_after(pos_1);
            remove_lagging_promoters_starting_before(pos_2);
        } else {

            // Delete RNAs until we pass pos_1 (or we reach the end of the list)
            auto init_loop = lagging_prom_pos_.lower_bound(pos_1);
            for (auto it = init_loop,
                         nextit = it;
                 it != lagging_prom_pos_.end() and it->first < pos_2;
                 it = nextit) {

                delete promoters_[it->second];
                promoters_.erase(it->second);
                nextit = next(it);
                lagging_prom_pos_.erase(it);
            }
        }
    }

    void SIMD_Map_Metadata::remove_lagging_promoters_starting_after(int32_t pos) {
        auto init_loop = lagging_prom_pos_.lower_bound(pos);

        if (init_loop == lagging_prom_pos_.end())
            return;

        // Delete RNAs until we pass pos (or we reach the end of the list)
        for (auto it = init_loop,
                     nextit = it;
             it != lagging_prom_pos_.end();
             it = nextit) {
            delete promoters_[it->second];
            promoters_.erase(it->second);
            nextit = next(it);
            lagging_prom_pos_.erase(it);
        }
    }

    void SIMD_Map_Metadata::remove_lagging_promoters_starting_before(int32_t pos) {
        // Delete RNAs until we reach pos (or we reach the end of the list)
        // TODO: optimize by starting from the end (with reverse iterators)
        auto init_loop = lagging_prom_pos_.lower_bound(pos);
        if (init_loop == lagging_prom_pos_.begin())
            return;
  /*      if (indiv_->indiv_id == 797 && time() == 9) {
            printf("Deleting from start to %d (%d)\n", pos, init_loop->first);
        }*/
        for (auto it = lagging_prom_pos_.begin(),
                     nextit = it;
             it != init_loop;
             it = nextit) {
            delete promoters_[it->second];
            promoters_.erase(it->second);
            nextit = next(it);
            lagging_prom_pos_.erase(it);
/*            if (indiv_->indiv_id == 797 && time() == 9) {
                printf("Delete %d\n", it->first);
            }*/
        }
    }

    void SIMD_Map_Metadata::move_all_leading_promoters_after(int32_t pos, int32_t delta_pos) {
        std::map<int32_t,int32_t> tmp_prom;

        for (auto it = leading_prom_pos_.lower_bound(pos), nextit=it;
             it != leading_prom_pos_.end();
             it = nextit) {

            int32_t new_pos = Utils::mod(it->first + delta_pos, length());
            int32_t prom_idx = it->second;

/*            printf("LEAD -- Moving %d to %d (%d)\n",promoters_[it->second]->pos,
                   Utils::mod(promoters_[it->second]->pos + delta_pos, length()),length());*/

            promoters_[it->second]->pos = new_pos;
            nextit = next(it);


            if (tmp_prom.find(new_pos) == tmp_prom.end()) {
                tmp_prom[new_pos] = prom_idx;
            } else {
                promoters_.erase(it->second);
            }
            leading_prom_pos_.erase(it);

        }
        for (auto to_insert : tmp_prom) {
            if (leading_prom_pos_.find(to_insert.first) == leading_prom_pos_.end()) {
                leading_prom_pos_[to_insert.first] = to_insert.second;
            } else {
                promoters_.erase(to_insert.second);
            }
        }
    }

    void SIMD_Map_Metadata::move_all_lagging_promoters_after(int32_t pos,int32_t delta_pos) {
        std::map<int32_t,int32_t> tmp_prom;

        for (auto it = lagging_prom_pos_.lower_bound(pos), nextit = it;
             it != lagging_prom_pos_.end() and it->first >= pos;
             it=nextit) {
            int32_t new_pos = Utils::mod(it->first + delta_pos, length());
            int32_t prom_idx = it->second;

/*            printf("LAG -- Moving %d to %d (%d)\n",promoters_[it->second]->pos,
                   Utils::mod(promoters_[it->second]->pos + delta_pos, length()),length());*/

            promoters_[it->second]->pos = new_pos;
            nextit = next(it);

            if (tmp_prom.find(new_pos) == tmp_prom.end()) {
                tmp_prom[new_pos] = prom_idx;
            } else {
                promoters_.erase(it->second);
            }

            lagging_prom_pos_.erase(it);
        }

        for (auto to_insert : tmp_prom) {
            if (lagging_prom_pos_.find(to_insert.first) == lagging_prom_pos_.end()) {
                lagging_prom_pos_[to_insert.first] = to_insert.second;
            }else {
                promoters_.erase(to_insert.second);
            }
        }
    }

    void SIMD_Map_Metadata::look_for_new_leading_promoters_starting_between(int32_t pos_1, int32_t pos_2) {
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
                if (leading_prom_pos_.find(i) == leading_prom_pos_.end()) {
                    promoterStruct* nprom = new promoterStruct(i, dist, true);
                    {
                        int prom_idx;
#pragma omp atomic capture
                        {
                            prom_idx = count_promoters_;
                            count_promoters_ = count_promoters_ + 1;
                        }

                        promoters_[prom_idx] = nprom;
                        leading_prom_pos_[i] = prom_idx;
                    }
                }
            }
        }
    }

    void SIMD_Map_Metadata::look_for_new_leading_promoters_starting_after(int32_t pos) {
        // Hamming distance of the sequence from the promoter consensus
        //int8_t dist = 8;

        for (int32_t i = pos; i < length(); i++) {
#ifdef WITH_BITSET
            int8_t dist = dna_->bitset_->is_promoter(true,i);
#else
            int8_t dist = is_promoter_leading(i);
#endif
            if (dist <= 4) { // dist takes the hamming distance of the sequence from the consensus
                if (leading_prom_pos_.find(i) == leading_prom_pos_.end()) {
                    promoterStruct* nprom = new promoterStruct(i, dist, true);
                    {
                        int prom_idx;
#pragma omp atomic capture
                        {
                            prom_idx = count_promoters_;
                            count_promoters_ = count_promoters_ + 1;
                        }

                        promoters_[prom_idx] = nprom;
                        leading_prom_pos_[i] = prom_idx;
                    }
                }
            }
        }
    }

    void SIMD_Map_Metadata::look_for_new_leading_promoters_starting_before(int32_t pos) {
        // Hamming distance of the sequence from the promoter consensus

        for (int32_t i = 0; i < pos; i++) {
#ifdef WITH_BITSET
            int8_t dist = dna_->bitset_->is_promoter(true,i);
#else
            int8_t dist = is_promoter_leading(i);
#endif
            if (dist <= 4) { // dist takes the hamming distance of the sequence from the consensus
                if (leading_prom_pos_.find(i) == leading_prom_pos_.end()) {
                    promoterStruct* nprom = new promoterStruct(i, dist, true);
                    {
                        int prom_idx;
#pragma omp atomic capture
                        {
                            prom_idx = count_promoters_;
                            count_promoters_ = count_promoters_ + 1;
                        }

                        promoters_[prom_idx] = nprom;
                        leading_prom_pos_[i] = prom_idx;
                    }
                }
            }
        }
    }

    void SIMD_Map_Metadata::look_for_new_lagging_promoters_starting_between(int32_t pos_1,int32_t pos_2) {
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
                if (lagging_prom_pos_.find(i) == lagging_prom_pos_.end()) {
                    promoterStruct* nprom = new promoterStruct(i, dist, false);
                    {
                        int prom_idx;
#pragma omp atomic capture
                        {
                            prom_idx = count_promoters_;
                            count_promoters_ = count_promoters_ + 1;
                        }

                        promoters_[prom_idx] = nprom;
                        lagging_prom_pos_[i] = prom_idx;
                    }
                }
            }
        }
    }

    void SIMD_Map_Metadata::look_for_new_lagging_promoters_starting_after(int32_t pos) {
        // Hamming distance of the sequence from the promoter consensus

        for (int32_t i = length() - 1; i >= pos; i--) {
#ifdef WITH_BITSET
            int8_t dist = dna_->bitset_->is_promoter(false,i);
#else
            int8_t dist = is_promoter_lagging(i);
#endif
            if (dist <= 4) { // dist takes the hamming distance of the sequence from the consensus
                if (lagging_prom_pos_.find(i) == lagging_prom_pos_.end()) {
                    promoterStruct* nprom = new promoterStruct(i, dist, false);
                    {
                        int prom_idx;
#pragma omp atomic capture
                        {
                            prom_idx = count_promoters_;
                            count_promoters_ = count_promoters_ + 1;
                        }

                        promoters_[prom_idx] = nprom;
                        lagging_prom_pos_[i] = prom_idx;
                    }
                }
            }
        }
    }

    void SIMD_Map_Metadata::look_for_new_lagging_promoters_starting_before(int32_t pos) {
        // Hamming distance of the sequence from the promoter consensus

        for (int32_t i = pos - 1; i >= 0; i--) {
#ifdef WITH_BITSET
            int8_t dist = dna_->bitset_->is_promoter(false,i);
#else
            int8_t dist = is_promoter_lagging(i);
#endif
            if (dist <= 4) { // dist takes the hamming distance of the sequence from the consensus
                if (lagging_prom_pos_.find(i) == lagging_prom_pos_.end()) {
                    promoterStruct* nprom = new promoterStruct(i, dist, false);
                    {
                        int prom_idx;
#pragma omp atomic capture
                        {
                            prom_idx = count_promoters_;
                            count_promoters_ = count_promoters_ + 1;
                        }

                        promoters_[prom_idx] = nprom;
                        lagging_prom_pos_[i] = prom_idx;
                    }
                }
            }
        }
    }

    void SIMD_Map_Metadata::promoters_included_in(int32_t pos_1,
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
            int32_t seg_length = length() + pos_2 - pos_1;

            /*printf("promoters included in %d and %d\n",pos_1,pos_2);*/

            if (seg_length >= PROM_SIZE) {
                bool is_near_end_of_genome = (pos_1 + PROM_SIZE > length());
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
                                                           length(),
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
                                                           length(),
                                  promoters_list[LAGGING]);
                }
                else // is_near_end_of_genome && is_near_beginning_of_genome
                {
                    // promoters(leading, between, pos_1, pos_2 + dna_->length() - PROM_SIZE + 1,
                    //                                         promoters_list[LEADING]);
                    lst_promoters(LEADING, BETWEEN, pos_1, pos_2 - PROM_SIZE + 1 +
                                                           length(),
                                  promoters_list[LEADING]);
                    lst_promoters(LAGGING, BETWEEN, pos_2, pos_1 + PROM_SIZE - 1 -
                                                           length(),
                                  promoters_list[LAGGING]);
                }
            }
        }
    }

    void SIMD_Map_Metadata::extract_leading_promoters_starting_between(int32_t pos_1,
                                                    int32_t pos_2, std::list<promoterStruct*>& extracted_promoters) {
        if (pos_2 < pos_1) {

            auto first = leading_prom_pos_.lower_bound(pos_1);

            if (first == leading_prom_pos_.end() or first->first >= pos_2) {
                return;
            }

            // Extract the promoters (remove them from the individual's list and put them in extracted_promoters)

            for (auto it = first;
                 it != leading_prom_pos_.end();
                 it++) {
                extracted_promoters.push_back(promoters_[it->second]);
                promoters_.erase(it->second);
            }

            leading_prom_pos_.erase(first, leading_prom_pos_.end());

            // Find the last promoters in the interval
            auto end = leading_prom_pos_.lower_bound(pos_2);


            // Extract the promoters (remove them from the individual's list and put them in extracted_promoters)
            for (auto it = leading_prom_pos_.begin();
                 it != end;
                 it++) {
                extracted_promoters.push_back(promoters_[it->second]);
                promoters_.erase(it->second);
            }

            leading_prom_pos_.erase(leading_prom_pos_.begin(),end);

        } else {

            auto first = leading_prom_pos_.lower_bound(pos_1);

            if (first == leading_prom_pos_.end() or first->first >= pos_2) {
                return;
            }

            // Find the last promoters in the interval
            auto end = leading_prom_pos_.lower_bound(pos_2);


            // Extract the promoters (remove them from the individual's list and put them in extracted_promoters)
            for (auto it = first;
                 it != end;
                 it++) {
                extracted_promoters.push_back(promoters_[it->second]);
                promoters_.erase(it->second);
            }

            leading_prom_pos_.erase(first, end);
        }
    }

    void SIMD_Map_Metadata::extract_lagging_promoters_starting_between(int32_t pos_1,
                                                    int32_t pos_2,
                                                    std::list<promoterStruct*>& extracted_promoters) {
        if (pos_1 > pos_2) {
            // From pos_1 to start
            // Find the last promoters in the interval
            auto end = lagging_prom_pos_.lower_bound(pos_1);

            if (end != lagging_prom_pos_.begin()) end--;

            // Extract the promoters (remove them from the individual's list and put the in extracted_promoters)
            for (auto it = lagging_prom_pos_.begin();
                 it != end;
                 it++) {
                extracted_promoters.push_back(promoters_[it->second]);
                promoters_.erase(it->second);
            }

            lagging_prom_pos_.erase(lagging_prom_pos_.begin(), end);

            // From end to pos_2

            auto first = lagging_prom_pos_.lower_bound(pos_2);

            if (first != lagging_prom_pos_.begin()) first--;

            if (first == lagging_prom_pos_.end() or first->first < pos_2) {
                return;
            }

            for (auto it = first;
                 it != lagging_prom_pos_.end();
                 it++) {
                extracted_promoters.push_back(promoters_[it->second]);
                promoters_.erase(it->second);
            }

            lagging_prom_pos_.erase(first,lagging_prom_pos_.end());

        } else {

            auto first = lagging_prom_pos_.lower_bound(pos_2);

            // Find the last promoters in the interval
            auto end = lagging_prom_pos_.lower_bound(pos_1);

            if (end == lagging_prom_pos_.end() or end->first < pos_1)
                return;

            // Extract the promoters (remove them from the individual's list and put the in extracted_promoters)
            for (auto it = end;
                 it != first;
                 it++) {
                extracted_promoters.push_back(promoters_[it->second]);
                promoters_.erase(it->second);
            }

            lagging_prom_pos_.erase(end, first);
        }
    }

    promoterStruct *SIMD_Map_Metadata::promoters(int idx) {
        return promoters_[idx];
    }

    void SIMD_Map_Metadata::promoter_add(int idx, promoterStruct *prom) {
        promoters_[idx] = prom;

        if (prom->leading_or_lagging)
            leading_prom_pos_[prom->pos] = idx;
        else
            lagging_prom_pos_[prom->pos] = idx;
    }

    int SIMD_Map_Metadata::promoter_count() {
        return count_promoters_;
    }

    void SIMD_Map_Metadata::set_promoters_count(int pcount) {
        count_promoters_ = pcount;
    }

    int SIMD_Map_Metadata::terminator_count(int LoL) {
        if (LoL == LEADING)
            return (int) terminator_lead_.size();
        else
            return (int) terminator_lag_.size();
    }

    int SIMD_Map_Metadata::terminator_add(int LoL, int dna_pos) {
        if (LoL == LEADING)
            terminator_lead_.insert(dna_pos);
        else
            terminator_lag_.insert(dna_pos);
    }

    int SIMD_Map_Metadata::next_terminator(int LoL, int dna_pos) {
        if (LoL == LEADING) {
            auto it_rna_end = terminator_lead_.lower_bound(dna_pos);

            if (it_rna_end == terminator_lead_.end()) {
                it_rna_end = terminator_lead_.begin();
            }

            return *it_rna_end;
        } else {
            auto it_rna_end = terminator_lag_.upper_bound(dna_pos);


            if (it_rna_end == terminator_lag_.begin()) {
                it_rna_end = terminator_lag_.end();
                it_rna_end--;
            } else if ((*it_rna_end) != dna_pos)
                it_rna_end--;

            return *it_rna_end;
        }

    }

    void SIMD_Map_Metadata::terminators_clear() {
        terminator_lead_.clear();
        terminator_lag_.clear();
    }

    pRNA *SIMD_Map_Metadata::rnas(int idx) {
        return rnas_[idx];
    }

    void SIMD_Map_Metadata::rna_add(int idx, pRNA *rna) {
        rnas_[idx] = rna;
    }

    int SIMD_Map_Metadata::rna_count() {
        return rna_count_;
    }

    void SIMD_Map_Metadata::set_rna_count(int rcount) {
        rna_count_ = rcount;
    }

    void SIMD_Map_Metadata::rnas_resize(int resize) {
        rnas_.resize(resize);
    }

    void SIMD_Map_Metadata::rnas_clear() {
        rnas_.clear();
    }

    pProtein *SIMD_Map_Metadata::proteins(int idx) {
        return proteins_[idx];
    }

    void SIMD_Map_Metadata::protein_add(int idx, pProtein *prot) {
        proteins_[idx] = prot;
    }

    int SIMD_Map_Metadata::proteins_count() {
        return protein_count_;
    }

    void SIMD_Map_Metadata::set_proteins_count(int pcount) {
        protein_count_ = pcount;
    }

    void SIMD_Map_Metadata::proteins_resize(int resize) {
        proteins_.resize(resize);
    }

    void SIMD_Map_Metadata::proteins_clear() {
        proteins_.clear();
    }
}
