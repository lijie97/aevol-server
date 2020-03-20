//
// Created by arrouan on 26/09/19.
//

#include "SIMD_List_Metadata.h"
#include "AeTime.h"
#include <list>
#include <algorithm>
#include <iterator>

namespace aevol {
    void SIMD_List_Metadata::lst_promoters(bool lorl,
                       Position before_after_btw, // with regard to the strand's reading direction
                       int32_t pos1,
                       int32_t pos2,
                       std::list<promoterStruct*>&  motif_list) {
        Strand strand_id;
        //if (lorl)
            strand_id = (lorl==LEADING) ? LEADING : LAGGING;
        //else
        //    strand_id = LAGGING;

        auto strand = promoters_list_[strand_id];
        auto it_begin = strand.begin();
        auto it_end = strand.end();

        if (before_after_btw != BEFORE) {
            it_begin = find_if(strand.begin(),
                               strand.end(),
                               [pos1, strand_id](promoterStruct& p) {
                                   if (strand_id == LEADING) {
                                       return p.pos >= pos1;
                                   }
                                   else {
                                       return p.pos < pos1;
                                   }
                               });
        }

        if (before_after_btw != AFTER) {
            it_end = find_if(it_begin,
                             strand.end(),
                             [pos2, strand_id](promoterStruct& p) {
                                 if (strand_id == LEADING) {
                                     return p.pos >= pos2;
                                 }
                                 else {
                                     return p.pos < pos2;
                                 }
                             });
        }

        std::list<promoterStruct*> promoters_1D;
        for (auto it = it_begin; it != it_end; it++) {
            promoterStruct* ptr = new promoterStruct(*it);
            promoters_1D.push_back(ptr);
        }

        motif_list.insert(motif_list.end(), promoters_1D.begin(), promoters_1D.end());
    }

    void SIMD_List_Metadata::remove_promoters_around(int32_t pos_1) {
        if (length() >= PROM_SIZE) {
            remove_leading_promoters_starting_between(Utils::mod(pos_1 - PROM_SIZE + 1,
                                                                 length()),
                                                      pos_1);

            remove_lagging_promoters_starting_between(pos_1,
                                                      Utils::mod(pos_1 + PROM_SIZE - 1,
                                                                 length()));
        }
        else {
            remove_all_promoters();
        }
    }

    void SIMD_List_Metadata::remove_promoters_around(int32_t pos_1, int32_t pos_2) {
        if (Utils::mod(pos_1 - pos_2, length()) >= PROM_SIZE) {
            remove_leading_promoters_starting_between(Utils::mod(pos_1 - PROM_SIZE + 1,
                                                                 length()),
                                                      pos_2);

            remove_lagging_promoters_starting_between(pos_1,
                                                      Utils::mod(pos_2 + PROM_SIZE - 1,
                                                                 length()));
        }
        else {
            remove_all_promoters();
        }
    }

    void SIMD_List_Metadata::remove_all_promoters() {
        promoters_list_[LEADING].clear();
        promoters_list_[LAGGING].clear();
    }

    void SIMD_List_Metadata::look_for_new_promoters_around(int32_t pos_1, int32_t pos_2) {
        if (length() >= PROM_SIZE) {
            look_for_new_leading_promoters_starting_between(
                    Utils::mod(pos_1 - PROM_SIZE + 1,
                               length()), pos_2);
            look_for_new_lagging_promoters_starting_between(pos_1, Utils::mod(
                    pos_2 + PROM_SIZE - 1,
                    length()));
        }
    }

    void SIMD_List_Metadata::look_for_new_promoters_around(int32_t pos) {
        if (length() >= PROM_SIZE) {
            look_for_new_leading_promoters_starting_between(
                    Utils::mod(pos - PROM_SIZE + 1, length()),
                    pos);
            look_for_new_lagging_promoters_starting_between(
                    pos,
                    Utils::mod(pos + PROM_SIZE - 1, length()));
        }
    }

    void SIMD_List_Metadata::locate_promoters() {
        int8_t dist; // Hamming distance of the sequence from the promoter consensus

        // Empty RNA list
        for (auto& strand: promoters_list_) {
            strand.clear();
        }

        if (length() < PROM_SIZE) {
            return;
        }

        for (int32_t i = 0; i < length(); i++) {
            int8_t dist = is_promoter_leading(i);
            if (dist <= 4) {// dist takes the hamming distance of the sequence from the consensus
                promoters_list_[LEADING].emplace_back(i, dist,true);

            }

            dist = is_promoter_lagging(i);
            if (dist <= 4) {
                promoters_list_[LAGGING].emplace_back(length() - i - 1,
                                                dist,false);
            }
        }
    }

    void SIMD_List_Metadata::move_all_promoters_after(int32_t pos, int32_t delta_pos) {
        move_all_leading_promoters_after(pos, delta_pos);
        move_all_lagging_promoters_after(pos, delta_pos);
    }

    void SIMD_List_Metadata::duplicate_promoters_included_in(int32_t pos_1,
                                         int32_t pos_2,
                                         std::vector<std::list<promoterStruct*>>& duplicated_promoters) {
        // 1) Get promoters to be duplicated
        std::vector<std::list<promoterStruct*>> retrieved_promoters = {{},
                                                 {}};
        promoters_included_in(pos_1, pos_2, retrieved_promoters);

        // 2) Set RNAs' position as their position on the duplicated segment
        for (auto& strand: {LEADING, LAGGING}) {
            for (auto& rna: retrieved_promoters[strand]) {
                // Make a copy of current RNA inside container
                duplicated_promoters[strand].emplace_back(rna);


                // Set RNA's position as it's position on the duplicated segment
                duplicated_promoters[strand].back()->pos = Utils::mod(duplicated_promoters[strand].back()->pos -pos_1, length());
            }
        }
    }

    void SIMD_List_Metadata::extract_promoters_included_in(int32_t pos_1,
                                       int32_t pos_2, std::vector<std::list<promoterStruct*>>& extracted_promoters) {
        if (pos_2 - pos_1 < PROM_SIZE) {
            return;
        }

        extract_leading_promoters_starting_between(pos_1, pos_2 - PROM_SIZE + 1,
                                                   extracted_promoters[LEADING]);

        extract_lagging_promoters_starting_between(pos_1 + PROM_SIZE - 1, pos_2,
                                                   extracted_promoters[LAGGING]);
    }

    void SIMD_List_Metadata::insert_promoters(std::vector<std::list<promoterStruct*>>& promoters_to_insert) {
        for (auto strand: {LEADING, LAGGING}) {
            if (promoters_to_insert[strand].size() <= 0) {
                continue;
            }
            // Get to the right position in individual's list (first promoter after the inserted segment)
            int32_t from_pos = promoters_to_insert[strand].back()->pos;

            auto pos = find_if(promoters_list_[strand].begin(),
                               promoters_list_[strand].end(),
                               [from_pos, strand](promoterStruct& r) {
                                   if (strand == LEADING) {
                                       return r.pos >= from_pos;
                                   }
                                   else {
                                       return r.pos < from_pos;
                                   }
                               });

            // Insert the promoters in the individual's RNA list
            for (auto& to_insert: promoters_to_insert[strand])
                // TODO vld: could be compacted in a unique emplace(pos, to_insert) ?
                if (pos != promoters_list_[strand].end()) {
                    promoters_list_[strand].insert(pos, *to_insert);
                }
                else {
                    promoters_list_[strand].push_back(*to_insert);
                }
        }
    }

    void SIMD_List_Metadata::insert_promoters_at(std::vector<std::list<promoterStruct*>>& promoters_to_insert,
                             int32_t pos) {
        for (auto strand: {LEADING, LAGGING}) {
            if (promoters_to_insert[strand].size() <= 0) {
                continue;
            }
            // Get to the right position in individual's list (first promoter after the inserted segment)
            auto first = find_if(promoters_list_[strand].begin(),
                                 promoters_list_[strand].end(),
                                 [pos, strand](promoterStruct& r) {
                                     if (strand == LEADING) {
                                         return r.pos >= pos;
                                     }
                                     else {
                                         return r.pos < pos;
                                     }
                                 });

            // Insert the promoters in the individual's RNA list
            for (auto& to_insert: promoters_to_insert[strand]) {
                // Update promoter position
                to_insert->pos = Utils::mod(to_insert->pos + pos, length()); //shift_position(pos, dna_->length());
                // Insert
                if (first != promoters_list_[strand].end()) {
                    promoters_list_[strand].insert(first, *to_insert);
                }
                else {
                    promoters_list_[strand].push_back(*to_insert);
                }
            }
        }
    }

    void SIMD_List_Metadata::invert_promoters_included_in(int32_t pos1,
                                      int32_t pos2) {
        int32_t segment_length = pos2 - pos1;

        if (segment_length < PROM_SIZE) {
            return;
        }

        std::vector<std::list<promoterStruct*>> inverted_promoters = {{},
                                                {}};

        // 1) Extract the promoters completely included on the segment to be inverted
        extract_promoters_included_in(pos1, pos2, inverted_promoters);

/*        if (indiv_->indiv_id == 37) {

            printf("Before look for 1+2\n");

            printf("%d -- %d -- EXTRACTED -- Prom list LEAD : ", AeTime::time(), indiv_->indiv_id);
            for (int prom_idx = 0; prom_idx < inverted_promoters[LEADING].size(); prom_idx++) {
                auto it =  inverted_promoters[LEADING].begin();
                std::advance(it, prom_idx);
                printf("%d ", (*it)->pos);
            }
            printf("\n");
            printf("%d -- %d -- EXTRACTED -- Prom list LAG : ", AeTime::time(), indiv_->indiv_id);
            for (int prom_idx = 0; prom_idx < inverted_promoters[LAGGING].size(); prom_idx++) {
                auto it =  inverted_promoters[LAGGING].begin();
                std::advance(it, prom_idx);
                printf("%d ", (*it)->pos);
            }
            printf("\n");
        }*/

        // 2) Invert segment's promoters
        SIMD_List_Metadata::invert_promoters(inverted_promoters, pos1, pos2);

/*        if (indiv_->indiv_id == 37) {

            printf("Before look for 1+2\n");

            printf("%d -- %d -- INVERTED -- Prom list LEAD : ", AeTime::time(), indiv_->indiv_id);
            for (int prom_idx = 0; prom_idx < inverted_promoters[LEADING].size(); prom_idx++) {
                auto it =  inverted_promoters[LEADING].begin();
                std::advance(it, prom_idx);
                printf("%d ", (*it)->pos);
            }
            printf("\n");
            printf("%d -- %d -- INVERTED -- Prom list LAG : ", AeTime::time(), indiv_->indiv_id);
            for (int prom_idx = 0; prom_idx < inverted_promoters[LAGGING].size(); prom_idx++) {
                auto it =  inverted_promoters[LAGGING].begin();
                std::advance(it, prom_idx);
                printf("%d ", (*it)->pos);
            }
            printf("\n");
        }*/


        // 3) Reinsert the inverted promoters
        insert_promoters(inverted_promoters);
    }


    void SIMD_List_Metadata::shift_promoters(
            std::vector<std::list<promoterStruct*>>& promoters_to_shift,
            int32_t delta_pos,
            int32_t seq_length) {
        for (auto& strand: {LEADING, LAGGING})
            for (auto& rna: promoters_to_shift[strand])
                rna->pos = Utils::mod(rna->pos + delta_pos, seq_length);
    }

    void SIMD_List_Metadata::invert_promoters(std::vector<std::list<promoterStruct*>>& promoter_lists,
                                 int32_t pos1,
                                 int32_t pos2) {
        // Exchange LEADING and LAGGING lists
        promoter_lists[LEADING].swap(promoter_lists[LAGGING]);

        //if (indiv_->indiv_id == 37) {



            /*printf("%d -- SWAP -- Prom list LEAD : ", AeTime::time());
            for (int prom_idx = 0; prom_idx < promoter_lists[LEADING].size(); prom_idx++) {
                auto it =  promoter_lists[LEADING].begin();
                std::advance(it, prom_idx);
                printf("%d ", (*it)->pos);
            }
            printf("\n");
            printf("%d -- SWAP -- Prom list LAG : ", AeTime::time());
            for (int prom_idx = 0; prom_idx < promoter_lists[LAGGING].size(); prom_idx++) {
                auto it =  promoter_lists[LAGGING].begin();
                std::advance(it, prom_idx);
                printf("%d ", (*it)->pos);
            }
            printf("\n");*/
        //}

        // Update the position and strand of each promoter to be inverted...
        for (auto strand: {LEADING, LAGGING})
            for (auto rna: promoter_lists[strand]) {
                //printf("Update %d : %d %d %d // %d\n",(pos1 + pos2 - rna->pos - 1),pos1,pos2,rna->pos,strand);
                rna->pos = pos1 + pos2 - rna->pos - 1;
                rna->leading_or_lagging = (strand==LEADING);
                //printf("Update DD %d : %d %d %d // %d\n",(pos1 + pos2 - rna->pos - 1),pos1,pos2,rna->pos,strand);
            }



/*        printf("%d -- UPDATED_LIST -- Prom list LEAD : ", AeTime::time());
        for (int prom_idx = 0; prom_idx < promoter_lists[LEADING].size(); prom_idx++) {
            auto it =  promoter_lists[LEADING].begin();
            std::advance(it, prom_idx);
            printf("%d ", (*it)->pos);
        }
        printf("\n");
        printf("%d -- UPDATED_LIST -- Prom list LAG : ", AeTime::time());
        for (int prom_idx = 0; prom_idx < promoter_lists[LAGGING].size(); prom_idx++) {
            auto it =  promoter_lists[LAGGING].begin();
            std::advance(it, prom_idx);
            printf("%d ", (*it)->pos);
        }
        printf("\n");*/
    }

    void SIMD_List_Metadata::remove_leading_promoters_starting_between(int32_t pos_1,
                                                   int32_t pos_2) {
        if (pos_1 > pos_2) {
            remove_leading_promoters_starting_after(pos_1);
            remove_leading_promoters_starting_before(pos_2);
        }
        else {
            auto& strand = promoters_list_[LEADING];
            // Delete RNAs until we pass pos_2 (or we reach the end of the list)
            // STL Warning: don't erase the current iterator in the for-loop!
            auto init_loop = find_if(strand.begin(),
                                     strand.end(),
                                     [pos_1](promoterStruct& r) {
                                         return r.pos >= pos_1;
                                     });

            for (auto it = init_loop,
                         nextit = it;
                 it != strand.end() and it->pos < pos_2;
                 it = nextit) {
                nextit = next(it);
                strand.erase(it);
            }
        }
    }

    void SIMD_List_Metadata::remove_leading_promoters_starting_after(int32_t pos) {
        auto& strand = promoters_list_[LEADING];

        auto init_loop = find_if(strand.begin(), strand.end(),
                                 [pos](promoterStruct& r) { return r.pos >= pos; });

        for (auto it = init_loop,
                     nextit = it;
             it != strand.end();
             it = nextit) {
            nextit = next(it);
            strand.erase(it);
        }
    }

    void SIMD_List_Metadata::remove_leading_promoters_starting_before(int32_t pos) {
        auto& strand = promoters_list_[LEADING];
        // Delete RNAs until we reach pos (or we reach the end of the list)
        for (auto it = strand.begin(),
                     nextit = it;
             it != strand.end() and it->pos < pos;
             it = nextit) {
            nextit = next(it);
            strand.erase(it);
        }
    }

    void SIMD_List_Metadata::remove_lagging_promoters_starting_between(int32_t pos_1,
                                                   int32_t pos_2) {
        if (pos_1 == length()) pos_1 = 0;
        if (pos_2 == 0) pos_2 = length();
        if (pos_1 >
            pos_2) { // vld: that's a weird case... really do this? used from remove_promoters_around()
            remove_lagging_promoters_starting_after(pos_1);
            remove_lagging_promoters_starting_before(pos_2);
        }
        else {
            auto& strand = promoters_list_[LAGGING];
            // Delete RNAs until we pass pos_1 (or we reach the end of the list)
            auto init_loop = find_if(strand.begin(),
                                     strand.end(),
                                     [pos_2](promoterStruct& r) {
                                         return r.pos < pos_2;
                                     });

            for (auto it = init_loop,
                         nextit = it;
                 it != strand.end() and it->pos >= pos_1;
                 it = nextit) {
                nextit = next(it);
                strand.erase(it);
            }
        }
    }

    void SIMD_List_Metadata::remove_lagging_promoters_starting_after(int32_t pos) {
        auto& strand = promoters_list_[LAGGING];
        // Delete RNAs until we pass pos (or we reach the end of the list)
        for (auto it = strand.begin(),
                     nextit = it;
             it != strand.end() and it->pos >= pos;
             it = nextit) {
            nextit = next(it);
            strand.erase(it);
        }
    }

    void SIMD_List_Metadata::remove_lagging_promoters_starting_before(int32_t pos) {
        auto& strand = promoters_list_[LAGGING];
        // Delete RNAs until we reach pos (or we reach the end of the list)

        auto init_loop = find_if(strand.begin(),
                                 strand.end(),
                                 [pos](promoterStruct& r) { return r.pos < pos; });

        for (auto it = init_loop,
                     nextit = it;
             it != strand.end();
             it = nextit) {
            nextit = next(it);
            strand.erase(it);
        }
    }

    void SIMD_List_Metadata::move_all_leading_promoters_after(int32_t pos, int32_t delta_pos) {
        auto& strand = promoters_list_[LEADING];
        auto init_loop = find_if(strand.begin(), strand.end(), [pos](promoterStruct& r) {
            return r.pos >= pos;
        });

        for (auto rna = init_loop;
             rna != strand.end();
             ++rna)
            (*rna).pos = Utils::mod((*rna).pos + delta_pos, length());
    }

    void SIMD_List_Metadata::move_all_lagging_promoters_after(int32_t pos,int32_t delta_pos) {
        auto& strand = promoters_list_[LAGGING];
        // Update RNAs until we pass pos (or we reach the end of the list)
        for (auto rna = strand.begin();
             rna != strand.end() and rna->pos >= pos;
             ++rna)
            (*rna).pos = Utils::mod((*rna).pos + delta_pos, length());    }

    void SIMD_List_Metadata::look_for_new_leading_promoters_starting_between(int32_t pos_1, int32_t pos_2) {
        if (pos_1 >= pos_2) {
            look_for_new_leading_promoters_starting_after(pos_1);
            look_for_new_leading_promoters_starting_before(pos_2);
            return;
        }
        int8_t dist; // Hamming distance of the sequence from the promoter consensus

        for (int32_t i = pos_1; i < pos_2; i++) {
            dist = is_promoter_leading(i);
            if (dist <= 4) // dist takes the hamming distance of the sequence from the consensus
            {

                // Look for the right place to insert the new promoter in the list
                auto& strand = promoters_list_[LEADING];

                auto first = find_if(strand.begin(),
                                     strand.end(),
                                     [i](promoterStruct& r) { return r.pos >= i; });


                if (first == strand.end() or first->pos != i) {
                    promoters_list_[LEADING].emplace(first, i, dist, true);
                }
            }
        }
    }

    void SIMD_List_Metadata::look_for_new_leading_promoters_starting_after(int32_t pos) {
// Hamming distance of the sequence from the promoter consensus
        int8_t dist;

        // rna list node used to find the new promoter's place in the list
        for (int32_t i = pos; i < length(); i++) {
            dist = is_promoter_leading(i);
            if (dist <= 4) { // dist takes the hamming distance of the sequence from the consensus
                // Look for the right place to insert the new promoter in the list
                auto& strand = promoters_list_[LEADING];
                auto first = find_if(strand.begin(),
                                     strand.end(),
                                     [i](promoterStruct& r) {
                                         return r.pos >= i;
                                     });

                if (first == strand.end() or first->pos != i) {
                    promoters_list_[LEADING].emplace(first, i, dist,true);
                }
            }
        }
    }

    void SIMD_List_Metadata::look_for_new_leading_promoters_starting_before(int32_t pos) {
        // Hamming distance of the sequence from the promoter consensus
        int8_t dist;

        auto& strand = promoters_list_[LEADING];
        auto first = strand.begin(); // TODO vld: should it not be reset at each loop step?

        for (int32_t i = 0; i < pos; i++) {
            dist = is_promoter_leading(i);
            if (dist <= 4) {
                // Look for the right place to insert the new promoter in the list

                first = find_if(first,
                                strand.end(),
                                [i](promoterStruct& r) { return r.pos >= i; });

                if (first == strand.end() or first->pos != i) {
                    promoters_list_[LEADING].emplace(first, i, dist, true);

                }
            }
        }
    }

    void SIMD_List_Metadata::look_for_new_lagging_promoters_starting_between(int32_t pos_1,int32_t pos_2) {
        if (pos_1 >= pos_2) {
            look_for_new_lagging_promoters_starting_after(pos_1);
            look_for_new_lagging_promoters_starting_before(pos_2);
            return;
        }

        int8_t dist; // Hamming distance of the sequence from the promoter consensus
        for (int32_t i = pos_2 - 1; i >= pos_1; i--) {
            dist = is_promoter_lagging(i);
            if (dist <= 4) {
                // Look for the right place to insert the new promoter in the list
                auto& strand = promoters_list_[LAGGING];

                auto first = find_if(strand.begin(),
                                     strand.end(),
                                     [i](promoterStruct& r) { return r.pos <= i; });

                if (first == strand.end() or first->pos != i) {
                    promoters_list_[LAGGING].emplace(first, i, dist, false);
                }
            }
        }
    }

    void SIMD_List_Metadata::look_for_new_lagging_promoters_starting_after(int32_t pos) {
// Hamming distance of the sequence from the promoter consensus
        int8_t dist;
        auto& strand = promoters_list_[LAGGING];
        auto first = strand.begin();

        for (int32_t i = length() - 1; i >= pos; i--) {
            dist = is_promoter_lagging(i);
            if (dist <= 4) {
                // Look for the right place to insert the new promoter in the list
                first = find_if(first,
                                strand.end(),
                                [i](promoterStruct& r) { return r.pos <= i; });

                if (first == strand.end() or first->pos != i) {
                    promoters_list_[LAGGING].emplace(first, i, dist, false);
                }
            }
        }
    }

    void SIMD_List_Metadata::look_for_new_lagging_promoters_starting_before(int32_t pos) {
        int8_t dist;

        // rna list node used to find the new promoter's place in the list
        auto& strand = promoters_list_[LAGGING];
        auto first = strand.begin();

        for (int32_t i = pos - 1; i >= 0; i--) {
            dist = is_promoter_lagging(i);
            if (dist <= 4) {
                assert (i >= 0 && i < length());
                // Look for the right place to insert the new promoter in the list
                first = find_if(first,
                                strand.end(),
                                [i](promoterStruct& r) {
                                    return r.pos <= i;
                                });

                if (first == strand.end() or first->pos != i) {
                    promoters_list_[LAGGING].emplace(first, i, dist, false);

                }
            }
        }
    }

    void SIMD_List_Metadata::promoters_included_in(int32_t pos_1,
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

            if (seg_length >= PROM_SIZE) {
                bool is_near_end_of_genome = (pos_1 + PROM_SIZE > length());
                bool is_near_beginning_of_genome = (pos_2 - PROM_SIZE < 0);

                if (!is_near_end_of_genome && !is_near_beginning_of_genome) {
                    lst_promoters(LEADING, AFTER, pos_1, -1, promoters_list[LEADING]);
                    lst_promoters(LEADING, BEFORE, -1, pos_2 - PROM_SIZE + 1,
                              promoters_list[LEADING]);
                    lst_promoters(LAGGING, AFTER, pos_2, -1, promoters_list[LAGGING]);
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

    void SIMD_List_Metadata::extract_leading_promoters_starting_between(int32_t pos_1,
                                                    int32_t pos_2, std::list<promoterStruct*>& extracted_promoters) {
        // Find the first promoters in the interval
        auto& strand = promoters_list_[LEADING];

        auto first = find_if(strand.begin(),
                             strand.end(),
                             [pos_1](promoterStruct& p) {
                                 return p.pos >= pos_1;
                             });

        if (first == strand.end() or first->pos >= pos_2) {
            return;
        }

        // Find the last promoters in the interval

        auto end = find_if(first,
                           strand.end(),
                           [pos_2](promoterStruct& p) { return p.pos >= pos_2; });

        // Extract the promoters (remove them from the individual's list and put them in extracted_promoters)
        std::list<promoterStruct*> promoters_1D;
        for (auto it = first; it != end; it++) {
            promoterStruct* ptr = new promoterStruct(*it);
            promoters_1D.push_back(ptr);
        }

        extracted_promoters.insert(extracted_promoters.end(), promoters_1D.begin(), promoters_1D.end());

        strand.erase(first, end);
    }

    void SIMD_List_Metadata::extract_lagging_promoters_starting_between(int32_t pos_1,
                                                    int32_t pos_2,
                                                    std::list<promoterStruct*>& extracted_promoters) {
        // Find the first promoters in the interval (if any)
        auto& strand = promoters_list_[LAGGING];

        auto first = find_if(strand.begin(),
                             strand.end(),
                             [pos_2](promoterStruct& r) {
                                 return r.pos < pos_2;
                             });

        if (first == strand.end() or first->pos < pos_1) {
            return;
        }

        // Find the last promoters in the interval
        auto end = find_if(first,
                           strand.end(),
                           [pos_1](promoterStruct& r) { return r.pos < pos_1; });

        // Extract the promoters (remove them from the individual's list and put the in extracted_promoters)
        std::list<promoterStruct*> promoters_1D;
        for (auto it = first; it != end; it++) {
            promoterStruct* ptr = new promoterStruct(*it);
            promoters_1D.push_back(ptr);
        }

        extracted_promoters.insert(extracted_promoters.end(), promoters_1D.begin(), promoters_1D.end());
        strand.erase(first, end);
    }

    promoterStruct *SIMD_List_Metadata::promoters(int idx) {
        //for (int i = 0; i <= idx; i++)
        if (idx >= promoters_list_[LEADING].size()) {
            auto it = promoters_list_[LAGGING].begin();
            std::advance(it, idx-promoters_list_[LEADING].size());
            return &*(it);
        } else {
            auto it = promoters_list_[LEADING].begin();
            std::advance(it,idx);
            //printf(" [POS %d] ",idx);
            return &*(it);
        }
    }

    void SIMD_List_Metadata::promoter_add(int idx, promoterStruct *prom) {
        if (prom->leading_or_lagging) {
            promoters_list_[LEADING].push_back(prom);
        } else {
            promoters_list_[LAGGING].push_front(prom);
        }
    }

    promoterStruct *SIMD_List_Metadata::promoter_next() {
        promoterStruct* prom = &*(it_promoter_);
        it_promoter_++;
        if (it_promoter_ == promoters_list_[LEADING].end())
            it_promoter_ = promoters_list_[LAGGING].begin();
        return prom;
    }

    void SIMD_List_Metadata::promoter_begin() {
        it_promoter_ = promoters_list_[LEADING].begin();
    }

    bool SIMD_List_Metadata::promoter_end() {
        return it_promoter_ == promoters_list_[LAGGING].end();
    }


    int SIMD_List_Metadata::promoter_count() {
        return promoters_list_[LEADING].size()+ promoters_list_[LAGGING].size();
    }

    void SIMD_List_Metadata::set_promoters_count(int pcount) {
        //exit(-10);
    }


    int SIMD_List_Metadata::terminator_count(int LoL) {
        if (LoL == LEADING)
            return (int) terminator_lead_.size();
        else
            return (int) terminator_lag_.size();
    }

    void SIMD_List_Metadata::terminator_add(int LoL, int dna_pos) {
        if (LoL == LEADING)
            terminator_lead_.insert(dna_pos);
        else
            terminator_lag_.insert(dna_pos);
    }

    int SIMD_List_Metadata::next_terminator(int LoL, int dna_pos) {
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

    void SIMD_List_Metadata::terminators_clear() {
        terminator_lead_.clear();
        terminator_lag_.clear();
    }

    pRNA *SIMD_List_Metadata::rnas(int idx) {
        auto it = rnas_.begin();
        std::advance(it, idx);
        return *(it);
    }

    void SIMD_List_Metadata::rna_add(int idx, pRNA *rna) {
        rnas_.push_front(rna);
        it_rna_ = rnas_.begin();
    }


    pRNA *SIMD_List_Metadata::rna_next() {
        pRNA* rna = *(it_rna_);
        it_rna_++;
        cmp_rna++;
        return rna;
    }

    void SIMD_List_Metadata::rna_begin() {
        it_rna_ = rnas_.begin();
    }

    bool SIMD_List_Metadata::rna_end() {
        return it_rna_ == rnas_.end();
    }

    int SIMD_List_Metadata::rna_count() {
        return rnas_.size();
    }

    void SIMD_List_Metadata::set_rna_count(int rcount) {
        //rna_count_ = rcount;
    }

    void SIMD_List_Metadata::rnas_resize(int resize) {
        //rnas_.resize(resize);
    }

    void SIMD_List_Metadata::rnas_clear() {
        rnas_.clear();
    }

    pProtein *SIMD_List_Metadata::proteins(int idx) {
        auto it = proteins_.begin();
        std::advance(it, idx);
        return *(it);
    }

    void SIMD_List_Metadata::protein_add(int idx, pProtein *prot) {
        proteins_.push_front(prot);
    }


    pProtein *SIMD_List_Metadata::protein_next() {
        pProtein* prot = *(it_protein_);
        it_protein_++;
        return prot;
    }

    void SIMD_List_Metadata::protein_begin() {
        it_protein_ = proteins_.begin();
    }

    bool SIMD_List_Metadata::protein_end() {
        return it_protein_ == proteins_.end();
    }

    int SIMD_List_Metadata::proteins_count() {
        return proteins_.size();
    }

    void SIMD_List_Metadata::set_proteins_count(int pcount) {
        //protein_count_ = pcount;
    }

    void SIMD_List_Metadata::proteins_resize(int resize) {
        //proteins_.resize(resize);
    }

    void SIMD_List_Metadata::proteins_clear() {
        proteins_.clear();
    }
}