//
// Created by arrouan on 18/07/19.
//

#include "SIMD_DynTab_Metadata.h"
#include "AeTime.h"

namespace aevol {




    void SIMD_DynTab_Metadata::lst_promoters(bool lorl, Position before_after_btw, int32_t pos1, int32_t pos2,
                                          std::list<promoterStruct*>& motif_list) {

        /*
         * pos_1    >=      LEADING
         * pos_1    <       LAGGING
         *
         * pos_2    >=      LEADING
         * pos_2    <       LAGGING
         */

        for (int prom_idx = 0; prom_idx < count_promoters_; prom_idx++) {
            if (promoters_[prom_idx] != nullptr) {


                if (before_after_btw == BEFORE) {
                    //  begin -> pos_2
                    if (promoters_[prom_idx]->leading_or_lagging && lorl == LEADING) {
                        if (promoters_[prom_idx]->pos < pos2) {
//                            printf("Adding -- Looking on %d (%d) for %d in %d %d : %d\n",lorl,promoters_[prom_idx]->leading_or_lagging,
//                                   before_after_btw,pos1,pos2,promoters_[prom_idx]->pos);
                            motif_list.push_back(promoters_[prom_idx]);
                        }
                    } else if (!promoters_[prom_idx]->leading_or_lagging && lorl == LAGGING)
                        if (promoters_[prom_idx]->pos >= pos2) {
//                            printf("Adding -- Looking on %d (%d) for %d in %d %d : %d\n", lorl,
//                                   promoters_[prom_idx]->leading_or_lagging,
//                                   before_after_btw, pos1, pos2, promoters_[prom_idx]->pos);
                            motif_list.push_back(promoters_[prom_idx]);
                        }
                } else if (before_after_btw == AFTER) {
                    // pos_1 -> end
                    if (promoters_[prom_idx]->leading_or_lagging && lorl == LEADING) {
                        if (promoters_[prom_idx]->pos >= pos1) {
//                            printf("Adding -- Looking on %d (%d) for %d in %d %d : %d\n",lorl,promoters_[prom_idx]->leading_or_lagging,
//                                   before_after_btw,pos1,pos2,promoters_[prom_idx]->pos);
                            motif_list.push_back(promoters_[prom_idx]);
                        }
                    } else if (!promoters_[prom_idx]->leading_or_lagging && lorl == LAGGING)
                        if (promoters_[prom_idx]->pos < pos1) {
//                            printf("Adding -- Looking on %d (%d) for %d in %d %d : %d\n", lorl,
//                                   promoters_[prom_idx]->leading_or_lagging,
//                                   before_after_btw, pos1, pos2, promoters_[prom_idx]->pos);
                            motif_list.push_back(promoters_[prom_idx]);
                        }
                } else {
                    // pos_1 -> pos_2
                    if (promoters_[prom_idx]->leading_or_lagging && lorl == LEADING) {
                        if ((promoters_[prom_idx]->pos >= pos1) && (promoters_[prom_idx]->pos < pos2)) {
//                            printf("Adding -- Looking on %d (%d) for %d in %d %d : %d\n",lorl,promoters_[prom_idx]->leading_or_lagging,
//                                   before_after_btw,pos1,pos2,promoters_[prom_idx]->pos);
                            motif_list.push_back(promoters_[prom_idx]);
                        }
                    } else if (!promoters_[prom_idx]->leading_or_lagging && lorl == LAGGING)
                        if ((promoters_[prom_idx]->pos < pos1) && (promoters_[prom_idx]->pos >= pos2)) {
//                            printf("Adding -- Looking on %d (%d) for %d in %d %d : %d\n", lorl,
//                                   promoters_[prom_idx]->leading_or_lagging,
//                                   before_after_btw, pos1, pos2, promoters_[prom_idx]->pos);
                            motif_list.push_back(promoters_[prom_idx]);
                        }
                }
            }
        }

    }

    void SIMD_DynTab_Metadata::remove_promoters_around(int32_t pos) {
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

    void SIMD_DynTab_Metadata::remove_promoters_around(int32_t pos_1, int32_t pos_2) {
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

    void SIMD_DynTab_Metadata::remove_all_promoters() {

        for (int prom_idx = 0; prom_idx < count_promoters_; prom_idx++) {
            delete promoters_[prom_idx];
        }

        count_promoters_ = 0;
    }

    void SIMD_DynTab_Metadata::look_for_new_promoters_around(int32_t pos_1, int32_t pos_2) {
        if (length() >= PROM_SIZE) {
            look_for_new_leading_promoters_starting_between(
                    Utils::mod(pos_1 - PROM_SIZE + 1,
                               length()), pos_2);

            look_for_new_lagging_promoters_starting_between(pos_1, Utils::mod(
                    pos_2 + PROM_SIZE - 1,
                    length()));
        }
    }

    void SIMD_DynTab_Metadata::look_for_new_promoters_around(int32_t pos) {
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

    void SIMD_DynTab_Metadata::locate_promoters() {
        look_for_new_leading_promoters_starting_between(0,length());
        look_for_new_lagging_promoters_starting_between(0,length());
    }

    void SIMD_DynTab_Metadata::move_all_promoters_after(int32_t pos, int32_t delta_pos) {
        move_all_leading_promoters_after(pos, delta_pos);
        move_all_lagging_promoters_after(pos, delta_pos);
    }

    void SIMD_DynTab_Metadata::duplicate_promoters_included_in(int32_t pos_1, int32_t pos_2,
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

    void SIMD_DynTab_Metadata::extract_promoters_included_in(int32_t pos_1, int32_t pos_2,
                                                          std::vector<std::list<promoterStruct *>> &extracted_promoters) {
        if (pos_2 - pos_1 < PROM_SIZE) {
            return;
        }

        extract_leading_promoters_starting_between(pos_1, pos_2 - PROM_SIZE + 1,
                                                   extracted_promoters[LEADING]);
        extract_lagging_promoters_starting_between(pos_1 + PROM_SIZE - 1, pos_2,
                                                   extracted_promoters[LAGGING]);
    }

    void SIMD_DynTab_Metadata::insert_promoters(std::vector<std::list<promoterStruct *>> &promoters_to_insert) {
        if (count_promoters_ + promoters_to_insert[LEADING].size() + promoters_to_insert[LAGGING].size() >= dyntab_size_)
            reallocate_promoters();

        for (auto strand: {LEADING, LAGGING}) {
            if (promoters_to_insert[strand].empty()) {
                continue;
            }
            // Insert the promoters in the individual's RNA list
            for (auto& to_insert: promoters_to_insert[strand]) {

                bool to_add = true;
                if (strand == LEADING) {
                    for (int prom_idx = 0; prom_idx < count_promoters_; prom_idx++) {
                        if (promoters_[prom_idx] != nullptr) {
                            if ((promoters_[prom_idx]->pos == to_insert->pos) &&
                                (promoters_[prom_idx]->leading_or_lagging == to_insert->leading_or_lagging)) {
                                to_add = false;
                                break;
                            }
                        }
                    }
                } else {
                    for (int prom_idx = 0; prom_idx < count_promoters_; prom_idx++) {
                        if (promoters_[prom_idx] != nullptr) {
                            if ((promoters_[prom_idx]->pos == to_insert->pos) &&
                                (promoters_[prom_idx]->leading_or_lagging == to_insert->leading_or_lagging)) {
                                to_add = false;
                                break;
                            }
                        }
                    }
                }
                if (to_add) {
                    int prom_idx;
//#pragma omp atomic capture
                    {
                        prom_idx = count_promoters_;
                        count_promoters_ = count_promoters_ + 1;
                    }
                    //printf("%d -- Adding a promoters %d (ouf of %d) --> %ld \n",AeTime::time(),prom_idx,dyntab_size_,indiv_->dna_->length());
                    promoters_[prom_idx] = to_insert;
                }

            }
        }
    }

    void SIMD_DynTab_Metadata::insert_promoters_at(std::vector<std::list<promoterStruct*>>& promoters_to_insert,
                             int32_t pos) {
        if (count_promoters_ + promoters_to_insert[LEADING].size() + promoters_to_insert[LAGGING].size() >= dyntab_size_)
            reallocate_promoters();

        for (auto strand: {LEADING, LAGGING}) {
            if (promoters_to_insert[strand].size() <= 0) {
                continue;
            }

            // Insert the promoters in the individual's RNA list
            for (auto &to_insert: promoters_to_insert[strand]) {
                // Update promoter position
                to_insert->pos = Utils::mod(to_insert->pos + pos, length());
                bool to_add = true;
                if (strand == LEADING) {
                    for (int prom_idx = 0; prom_idx < count_promoters_; prom_idx++) {
                        if (promoters_[prom_idx] != nullptr) {
                            if ((promoters_[prom_idx]->pos == to_insert->pos) &&
                                (promoters_[prom_idx]->leading_or_lagging == to_insert->leading_or_lagging)) {
                                to_add = false;
                                break;
                            }
                        }
                    }
                } else {
                    for (int prom_idx = 0; prom_idx < count_promoters_; prom_idx++) {
                        if (promoters_[prom_idx] != nullptr) {
                            if ((promoters_[prom_idx]->pos == to_insert->pos) &&
                                (promoters_[prom_idx]->leading_or_lagging == to_insert->leading_or_lagging)) {
                                to_add = false;
                                break;
                            }
                        }
                    }
                }
                if (to_add) {
                    int prom_idx;
//#pragma omp atomic capture
                    {
                        prom_idx = count_promoters_;
                        count_promoters_ = count_promoters_ + 1;
                    }
                    promoters_[prom_idx] = to_insert;
                }

            }
        }
    }

    void SIMD_DynTab_Metadata::invert_promoters_included_in(int32_t pos1,
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
        SIMD_DynTab_Metadata::invert_promoters(inverted_promoters, pos1, pos2);

        // 3) Reinsert the inverted promoters
        insert_promoters(inverted_promoters);
    }

    void SIMD_DynTab_Metadata::shift_promoters(
            std::vector<std::list<promoterStruct*>>& promoters_to_shift,
            int32_t delta_pos,
            int32_t seq_length) {
        for (auto& strand: {LEADING, LAGGING})
            for (auto& prom: promoters_to_shift[strand])
                prom->pos = Utils::mod(prom->pos + delta_pos, seq_length);
    }

    void SIMD_DynTab_Metadata::invert_promoters(std::vector<std::list<promoterStruct*>>& promoter_lists,
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

    void SIMD_DynTab_Metadata::remove_leading_promoters_starting_between(int32_t pos_1,
                                                   int32_t pos_2) {
        if (pos_1 > pos_2) {
            remove_leading_promoters_starting_after(pos_1);
            remove_leading_promoters_starting_before(pos_2);
        }
        else {
            for (int prom_idx = 0; prom_idx < count_promoters_; prom_idx++) {
                if (promoters_[prom_idx] != nullptr) {
                    if (promoters_[prom_idx]->leading_or_lagging)
                        if ((promoters_[prom_idx]->pos >= pos_1) && (promoters_[prom_idx]->pos < pos_2)) {
                            delete promoters_[prom_idx];
                            promoters_[prom_idx] = nullptr;
                        }
                }
            }
        }

    }
    void SIMD_DynTab_Metadata::remove_leading_promoters_starting_after(int32_t pos) {
        for (int prom_idx = 0; prom_idx < count_promoters_; prom_idx++) {
            if (promoters_[prom_idx] != nullptr) {
                if (promoters_[prom_idx]->leading_or_lagging)
                    if (promoters_[prom_idx]->pos >= pos) {
                        delete promoters_[prom_idx];
                        promoters_[prom_idx] = nullptr;
                    }
            }
        }
    }

    void SIMD_DynTab_Metadata::remove_leading_promoters_starting_before(int32_t pos) {
        for (int prom_idx = 0; prom_idx < count_promoters_; prom_idx++) {
            if (promoters_[prom_idx] != nullptr) {
                if (promoters_[prom_idx]->leading_or_lagging)
                    if (promoters_[prom_idx]->pos < pos) {
                        delete promoters_[prom_idx];
                        promoters_[prom_idx] = nullptr;
                    }
            }
        }
    }

    void SIMD_DynTab_Metadata::remove_lagging_promoters_starting_between(int32_t pos_1,
                                                   int32_t pos_2) {
        if (pos_1 == length()) pos_1 = 0;
        if (pos_2 == 0) pos_2 = length();

        if (pos_1 >
            pos_2) { // vld: that's a weird case... really do this? used from remove_promoters_around()
            remove_lagging_promoters_starting_after(pos_1);
            remove_lagging_promoters_starting_before(pos_2);
        } else {

            // Delete RNAs until we pass pos_1 (or we reach the end of the list)
            for (int prom_idx = 0; prom_idx < count_promoters_; prom_idx++) {
                if (promoters_[prom_idx] != nullptr) {
                    if (!promoters_[prom_idx]->leading_or_lagging)
                        if ((promoters_[prom_idx]->pos >= pos_1) && (promoters_[prom_idx]->pos < pos_2)) {
                            delete promoters_[prom_idx];
                            promoters_[prom_idx] = nullptr;
                        }
                }
            }
        }
    }

    void SIMD_DynTab_Metadata::remove_lagging_promoters_starting_after(int32_t pos) {
        for (int prom_idx = 0; prom_idx < count_promoters_; prom_idx++) {
            if (promoters_[prom_idx] != nullptr) {
                if (!promoters_[prom_idx]->leading_or_lagging)
                    if (promoters_[prom_idx]->pos >= pos) {
                        delete promoters_[prom_idx];
                        promoters_[prom_idx] = nullptr;
                    }
            }
        }
    }

    void SIMD_DynTab_Metadata::remove_lagging_promoters_starting_before(int32_t pos) {
        // Delete RNAs until we reach pos (or we reach the end of the list)
        // TODO: optimize by starting from the end (with reverse iterators)
        for (int prom_idx = 0; prom_idx < count_promoters_; prom_idx++) {
            if (promoters_[prom_idx] != nullptr) {
                if (!promoters_[prom_idx]->leading_or_lagging)
                    if (promoters_[prom_idx]->pos < pos) {
                        delete promoters_[prom_idx];
                        promoters_[prom_idx] = nullptr;
                    }
            }
        }
    }

    void SIMD_DynTab_Metadata::move_all_leading_promoters_after(int32_t pos, int32_t delta_pos) {
        for (int prom_idx = 0; prom_idx < count_promoters_; prom_idx++) {
            if (promoters_[prom_idx] != nullptr) {
                if (promoters_[prom_idx]->leading_or_lagging)
                    if (promoters_[prom_idx]->pos >= pos) {
                        /*printf("LEAD -- Moving %d to %d (%d)\n",promoters_[prom_idx]->pos,
                               Utils::mod(promoters_[prom_idx]->pos + delta_pos, length()),length());*/
                        promoters_[prom_idx]->pos = Utils::mod(promoters_[prom_idx]->pos + delta_pos, length());
                    }
            }
        }
    }

    void SIMD_DynTab_Metadata::move_all_lagging_promoters_after(int32_t pos,int32_t delta_pos) {
        for (int prom_idx = 0; prom_idx < count_promoters_; prom_idx++) {
            if (promoters_[prom_idx] != nullptr) {
                if (!promoters_[prom_idx]->leading_or_lagging)
                    if (promoters_[prom_idx]->pos >= pos) {
/*                        printf("LAG -- Moving %d to %d (%d)\n",promoters_[prom_idx]->pos,
                               Utils::mod(promoters_[prom_idx]->pos + delta_pos, length()),length());*/
                        promoters_[prom_idx]->pos = Utils::mod(promoters_[prom_idx]->pos + delta_pos, length());
                    }
            }
        }
    }

    void SIMD_DynTab_Metadata::look_for_new_leading_promoters_starting_between(int32_t pos_1, int32_t pos_2) {
        // When pos_1 > pos_2, we will perform the search in 2 steps.
        // As positions  0 and dna_->length() are equivalent, it's preferable to
        // keep 0 for pos_1 and dna_->length() for pos_2.


        if (pos_1 >= pos_2) {
            if (count_promoters_ + (pos_1 + (length() - pos_2))/PROM_SIZE >= dyntab_size_)
                reallocate_promoters();

            look_for_new_leading_promoters_starting_after(pos_1);
            look_for_new_leading_promoters_starting_before(pos_2);
            return;
        }
        // Hamming distance of the sequence from the promoter consensus

        if (count_promoters_ + (pos_2 - pos_1) / PROM_SIZE >= dyntab_size_)
            reallocate_promoters();

        for (int32_t i = pos_1; i < pos_2; i++) {
#ifdef WITH_BITSET
            int8_t dist = dna_->bitset_->is_promoter(true,i);
#else
            int8_t dist = is_promoter_leading(i);
#endif

            if (dist <= 4) {
                bool notFound = true;
                for (int prom_idx = 0; prom_idx < count_promoters_; prom_idx++) {
                    if (promoters_[prom_idx] != nullptr) {
                        if (promoters_[prom_idx]->leading_or_lagging)
                            if (promoters_[prom_idx]->pos == i) {
                                notFound = false;
                            }
                    }
                }

                if (notFound) {
                    int prom_idx;
//#pragma omp atomic capture
                    {
                        prom_idx = count_promoters_;
                        count_promoters_ = count_promoters_ + 1;
                    }

                    promoters_[prom_idx] = new promoterStruct(i, dist, true);
                }

            }
        }
    }

    void SIMD_DynTab_Metadata::look_for_new_leading_promoters_starting_after(int32_t pos) {
        // Hamming distance of the sequence from the promoter consensus
        //int8_t dist = 8;

        if (count_promoters_ + (length() - pos)/PROM_SIZE >= dyntab_size_)
            reallocate_promoters();

        for (int32_t i = pos; i < length(); i++) {
#ifdef WITH_BITSET
            int8_t dist = dna_->bitset_->is_promoter(true,i);
#else
            int8_t dist = is_promoter_leading(i);
#endif
            if (dist <= 4) { // dist takes the hamming distance of the sequence from the consensus
                bool notFound = true;
                for (int prom_idx = 0; prom_idx < count_promoters_; prom_idx++) {
                    if (promoters_[prom_idx] != nullptr) {
                        if (promoters_[prom_idx]->leading_or_lagging)
                            if (promoters_[prom_idx]->pos == i) {
                                notFound = false;
                            }
                    }
                }

                if (notFound) {
                    int prom_idx;
//#pragma omp atomic capture
                    {
                        prom_idx = count_promoters_;
                        count_promoters_ = count_promoters_ + 1;
                    }

                    promoters_[prom_idx] = new promoterStruct(i, dist, true);;
                }
            }
        }
    }

    void SIMD_DynTab_Metadata::look_for_new_leading_promoters_starting_before(int32_t pos) {
        // Hamming distance of the sequence from the promoter consensus
        if (count_promoters_ + (pos)/PROM_SIZE >= dyntab_size_)
            reallocate_promoters();

        for (int32_t i = 0; i < pos; i++) {
#ifdef WITH_BITSET
            int8_t dist = dna_->bitset_->is_promoter(true,i);
#else
            int8_t dist = is_promoter_leading(i);
#endif
            if (dist <= 4) { // dist takes the hamming distance of the sequence from the consensus
                bool notFound = true;
                for (int prom_idx = 0; prom_idx < count_promoters_; prom_idx++) {
                    if (promoters_[prom_idx] != nullptr) {
                        if (promoters_[prom_idx]->leading_or_lagging)
                            if (promoters_[prom_idx]->pos == i) {
                                notFound = false;
                            }
                    }
                }

                if (notFound) {
                    int prom_idx;
//#pragma omp atomic capture
                    {
                        prom_idx = count_promoters_;
                        count_promoters_ = count_promoters_ + 1;
                    }

                    promoters_[prom_idx] = new promoterStruct(i, dist, true);;
                }
            }
        }
    }

    void SIMD_DynTab_Metadata::look_for_new_lagging_promoters_starting_between(int32_t pos_1,int32_t pos_2) {
        // When pos_1 > pos_2, we will perform the search in 2 steps.
        // As positions  0 and dna_->length() are equivalent, it's preferable to
        // keep 0 for pos_1 and dna_->length() for pos_2.

        if (pos_1 >= pos_2) {
            if (count_promoters_ + (pos_1 + (length() - pos_2))/PROM_SIZE >= dyntab_size_)
                reallocate_promoters();

            look_for_new_lagging_promoters_starting_after(pos_1);
            look_for_new_lagging_promoters_starting_before(pos_2);
            return;
        }

        if (count_promoters_ + (pos_2 - pos_1) / PROM_SIZE >= dyntab_size_)
            reallocate_promoters();

        // Hamming distance of the sequence from the promoter consensus
        for (int32_t i = pos_2 - 1; i >= pos_1; i--) {
#ifdef WITH_BITSET
            int8_t dist = dna_->bitset_->is_promoter(false,i);
#else
            int8_t dist = is_promoter_lagging(i);
#endif
            if (dist <= 4) { // dist takes the hamming distance of the sequence from the consensus
                bool notFound = true;
                for (int prom_idx = 0; prom_idx < count_promoters_; prom_idx++) {
                    if (promoters_[prom_idx] != nullptr) {
                        if (!promoters_[prom_idx]->leading_or_lagging)
                            if (promoters_[prom_idx]->pos == i) {
                                notFound = false;
                            }
                    }
                }

                if (notFound) {
                    int prom_idx;
//#pragma omp atomic capture
                    {
                        prom_idx = count_promoters_;
                        count_promoters_ = count_promoters_ + 1;
                    }

                    promoters_[prom_idx] = new promoterStruct(i, dist, false);;
                }
            }
        }
    }

    void SIMD_DynTab_Metadata::look_for_new_lagging_promoters_starting_after(int32_t pos) {
        // Hamming distance of the sequence from the promoter consensus
        if (count_promoters_ + (length() - pos)/PROM_SIZE >= dyntab_size_)
            reallocate_promoters();

        for (int32_t i = length() - 1; i >= pos; i--) {
#ifdef WITH_BITSET
            int8_t dist = dna_->bitset_->is_promoter(false,i);
#else
            int8_t dist = is_promoter_lagging(i);
#endif
            if (dist <= 4) { // dist takes the hamming distance of the sequence from the consensus
                bool notFound = true;
                for (int prom_idx = 0; prom_idx < count_promoters_; prom_idx++) {
                    if (promoters_[prom_idx] != nullptr) {
                        if (!promoters_[prom_idx]->leading_or_lagging)
                            if (promoters_[prom_idx]->pos == i) {
                                notFound = false;
                            }
                    }
                }

                if (notFound) {
                    int prom_idx;
//#pragma omp atomic capture
                    {
                        prom_idx = count_promoters_;
                        count_promoters_ = count_promoters_ + 1;
                    }

                    promoters_[prom_idx] = new promoterStruct(i, dist, false);;
                }
            }
        }
    }

    void SIMD_DynTab_Metadata::look_for_new_lagging_promoters_starting_before(int32_t pos) {
        // Hamming distance of the sequence from the promoter consensus
        if (count_promoters_ + (pos)/PROM_SIZE >= dyntab_size_)
            reallocate_promoters();

        for (int32_t i = pos - 1; i >= 0; i--) {
#ifdef WITH_BITSET
            int8_t dist = dna_->bitset_->is_promoter(false,i);
#else
            int8_t dist = is_promoter_lagging(i);
#endif
            if (dist <= 4) { // dist takes the hamming distance of the sequence from the consensus
                bool notFound = true;
                for (int prom_idx = 0; prom_idx < count_promoters_; prom_idx++) {
                    if (promoters_[prom_idx] != nullptr) {
                        if (!promoters_[prom_idx]->leading_or_lagging)
                            if (promoters_[prom_idx]->pos == i) {
                                notFound = false;
                            }
                    }
                }

                if (notFound) {
                    int prom_idx;
//#pragma omp atomic capture
                    {
                        prom_idx = count_promoters_;
                        count_promoters_ = count_promoters_ + 1;
                    }

                    promoters_[prom_idx] = new promoterStruct(i, dist, false);;
                }
            }
        }
    }

    void SIMD_DynTab_Metadata::promoters_included_in(int32_t pos_1,
                               int32_t pos_2,
                               std::vector<std::list<promoterStruct*>>& promoters_list) {
        if (pos_1 < pos_2) {
//            printf("promoters << included in %d and %d\n",pos_1,pos_2);
            int32_t seg_length = pos_2 - pos_1;

            if (seg_length >= PROM_SIZE) {
//                printf("LOOKING FOR PROM in %d -> %d\n", pos_1, pos_2 - PROM_SIZE + 1);
                lst_promoters(LEADING, BETWEEN, pos_1, pos_2 - PROM_SIZE + 1,
                              promoters_list[LEADING]);
//                printf("PROM INCL IN PROMs LEAD : ");
//                for (auto prom :promoters_list[LEADING]) {
//                    if (prom != nullptr)
//                        printf("%d ",prom->pos);
//                }
//                printf("\n");
//
//                printf("LOOKING FOR PROM in %d -> %d\n",pos_2, pos_1 + PROM_SIZE - 1);
                lst_promoters(LAGGING, BETWEEN, pos_2, pos_1 + PROM_SIZE - 1,
                              promoters_list[LAGGING]);

//                printf("PROM INCL IN PROMs LAG : ");
//                for (auto prom :promoters_list[LAGGING]) {
//                    if (prom != nullptr)
//                        printf("%d ",prom->pos);
//                }
//                printf("\n");
            }
        }
        else {
            int32_t seg_length = length() + pos_2 - pos_1;

//            printf("promoters included in %d and %d\n",pos_1,pos_2);

            if (seg_length >= PROM_SIZE) {
                bool is_near_end_of_genome = (pos_1 + PROM_SIZE > length());
                bool is_near_beginning_of_genome = (pos_2 - PROM_SIZE < 0);

                if (!is_near_end_of_genome && !is_near_beginning_of_genome) {
//                    printf("-----------------> leading promoters after %d (till end)\n",pos_1);
                    lst_promoters(LEADING, AFTER, pos_1, -1, promoters_list[LEADING]);
//                    printf("-----------------> leading promoters before %d (till end)\n",pos_2 - PROM_SIZE + 1);
                    lst_promoters(LEADING, BEFORE, -1, pos_2 - PROM_SIZE + 1,
                                  promoters_list[LEADING]);
//                    printf("-----------------> lagging promoters after %d (till end)\n",pos_2);
                    lst_promoters(LAGGING, AFTER, pos_2, -1, promoters_list[LAGGING]);
//                    printf("-----------------> lagging promoters before %d (till end)\n",pos_1 + PROM_SIZE - 1);
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

    void SIMD_DynTab_Metadata::extract_leading_promoters_starting_between(int32_t pos_1,
                                                    int32_t pos_2, std::list<promoterStruct*>& extracted_promoters) {
        for (int prom_idx = 0; prom_idx < count_promoters_; prom_idx++) {
            if (promoters_[prom_idx] != nullptr) {
                if (promoters_[prom_idx]->leading_or_lagging)
                    if ((promoters_[prom_idx]->pos >= pos_1) && (promoters_[prom_idx]->pos < pos_2)) {
                        extracted_promoters.push_back(promoters_[prom_idx]);
                        promoters_[prom_idx] = nullptr;
                    }
            }
        }
    }


    void SIMD_DynTab_Metadata::extract_lagging_promoters_starting_between(int32_t pos_1,
                                                    int32_t pos_2,
                                                    std::list<promoterStruct*>& extracted_promoters) {

        for (int prom_idx = 0; prom_idx < count_promoters_; prom_idx++) {
            if (promoters_[prom_idx] != nullptr) {
                if (!promoters_[prom_idx]->leading_or_lagging)
                    if ((promoters_[prom_idx]->pos >= pos_1) && (promoters_[prom_idx]->pos < pos_2)) {
                        extracted_promoters.push_back(promoters_[prom_idx]);
                        promoters_[prom_idx] = nullptr;
                    }
            }
        }
    }

    promoterStruct *SIMD_DynTab_Metadata::promoters(int idx) {
        return promoters_[idx];
    }

    void SIMD_DynTab_Metadata::promoter_add(int idx, promoterStruct *prom) {
        promoters_[idx] = prom;
    }

    promoterStruct* SIMD_DynTab_Metadata::promoter_next() {
        promoterStruct* prom = promoters_[it_promoter_];
        it_promoter_++;
        return prom;
    }

    void SIMD_DynTab_Metadata::promoter_begin() {
        it_promoter_ = 0;
    }

    bool SIMD_DynTab_Metadata::promoter_end() {
        return it_promoter_ == promoter_count();
    }


    int SIMD_DynTab_Metadata::promoter_count() {
        return count_promoters_;
    }

    void SIMD_DynTab_Metadata::set_promoters_count(int pcount) {
        count_promoters_ = pcount;
    }

    int SIMD_DynTab_Metadata::terminator_count(int LoL) {
        if (LoL == LEADING)
            return (int) terminator_lead_.size();
        else
            return (int) terminator_lag_.size();
    }

    void SIMD_DynTab_Metadata::terminator_add(int LoL, int dna_pos) {
        if (LoL == LEADING)
            terminator_lead_.insert(dna_pos);
        else
            terminator_lag_.insert(dna_pos);
    }

    int SIMD_DynTab_Metadata::next_terminator(int LoL, int dna_pos) {
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

    void SIMD_DynTab_Metadata::terminators_clear() {
        terminator_lead_.clear();
        terminator_lag_.clear();
    }

    pRNA *SIMD_DynTab_Metadata::rnas(int idx) {
        return rnas_[idx];
    }

    void SIMD_DynTab_Metadata::rna_add(int idx, pRNA *rna) {
        rnas_[idx] = rna;
    }

    pRNA* SIMD_DynTab_Metadata::rna_next() {
        pRNA* rna = rnas_[it_rna_];
        it_rna_++;
        return rna;
    }

    void SIMD_DynTab_Metadata::rna_begin() {
        it_rna_ = 0;
    }

    bool SIMD_DynTab_Metadata::rna_end() {
        return it_rna_ == rna_count();
    }

    int SIMD_DynTab_Metadata::rna_count() {
        return rna_count_;
    }

    void SIMD_DynTab_Metadata::set_rna_count(int rcount) {
        rna_count_ = rcount;
    }

    void SIMD_DynTab_Metadata::rnas_resize(int resize) {
        rnas_.resize(resize);
    }

    void SIMD_DynTab_Metadata::rnas_clear() {
        rnas_.clear();
    }

    pProtein *SIMD_DynTab_Metadata::proteins(int idx) {
        return proteins_[idx];
    }

    void SIMD_DynTab_Metadata::protein_add(int idx, pProtein *prot) {
        proteins_[idx] = prot;
    }

    pProtein* SIMD_DynTab_Metadata::protein_next() {
        pProtein* prot = proteins_[it_protein_];
        it_protein_++;
        return prot;
    }

    void SIMD_DynTab_Metadata::protein_begin() {
        it_promoter_ = 0;
    }

    bool SIMD_DynTab_Metadata::protein_end() {
        return it_protein_ == proteins_count();
    }

    int SIMD_DynTab_Metadata::proteins_count() {
        return protein_count_;
    }

    void SIMD_DynTab_Metadata::set_proteins_count(int pcount) {
        protein_count_ = pcount;
    }

    void SIMD_DynTab_Metadata::proteins_resize(int resize) {
        proteins_.resize(resize);
    }

    void SIMD_DynTab_Metadata::proteins_clear() {
        proteins_.clear();
    }
}