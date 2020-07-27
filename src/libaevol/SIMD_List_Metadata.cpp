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

                //delete (rna);
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

        for (std::list<promoterStruct*>::iterator it_prom = inverted_promoters[LEADING].begin();
             it_prom != inverted_promoters[LEADING].end(); it_prom++) {
            delete (*(it_prom));
        }
        inverted_promoters[LEADING].clear();
        for (std::list<promoterStruct*>::iterator it_prom = inverted_promoters[LAGGING].begin();
             it_prom != inverted_promoters[LAGGING].end(); it_prom++) {
            delete (*(it_prom));
        }
        inverted_promoters[LAGGING].clear();
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
        auto first = strand.begin();

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

    void SIMD_List_Metadata::update_metadata_before_new_generation() {
        // Set to_delete, to_recompute, to_translate, to_compute (Protein/Rna) to false and to_compute to false for promoters
        for (auto &strand: promoters_list_)
            for (auto &prom: strand) {
                prom.to_compute = false;
                prom.to_delete = false;
            }

        for (auto rna : rnas_) {
            rna->to_recompute = false;
            rna->to_delete = false;
            rna->to_search_for_end = false;
        }

        for (auto prot : proteins_) {
            prot->to_recompute = false;
            prot->to_retranslate = false;
            prot->to_delete = false;
        }

        recompute_phenotype = false;
    }

    void SIMD_List_Metadata::remove_range(int32_t pos) {
        remove_range(pos,pos);
    }

    void SIMD_List_Metadata::remove_range(int32_t pos_1, int32_t pos_2) {

//        printf("======= BEGIN :: REMOVE RANGE [%d => %d]============\n",pos_1,pos_2);
//        display(false);
        // 0. Remove RNA with promoter
        // LEAD : within [pos_1-21,pos_2]
        for (int32_t cpos = pos_1-21; cpos <= pos_2; cpos++) {
            int32_t rcpos = Utils::mod(cpos,length());
            for (auto& prom: promoters_list_[LEADING]) {
                if (prom.pos == rcpos) {
//                    printf("Promoters to delete LEAD %d\n",prom.pos);
                    prom.to_delete = true;
                }
            }
        }

        // LAG  : within [pos_1,pos_2+21]
        for (int32_t cpos = pos_1; cpos <= pos_2+21; cpos++) {
            int32_t rcpos = Utils::mod(cpos,length());
            for (auto& prom: promoters_list_[LAGGING]) {
                if (prom.pos == rcpos) {
//                    printf("Promoters to delete LAG %d\n",prom.pos);
                    prom.to_delete = true;
                }
            }
        }

        // 1. Remove RNA with promoter
        // LEAD : within [pos_1-21,pos_2]
//        printf("Search for RNA between %d and %d\n",pos_1-21,pos_2);
        for (int32_t cpos = pos_1-21; cpos <= pos_2; cpos++) {
            int32_t rcpos = Utils::mod(cpos,length());
            for (auto rna : rnas_) {
//                printf("Process RNA %d %d\n",rna->begin,rna->leading_lagging);
                if (rna->leading_lagging == 0) {
//                    printf("To delete or NOT LEAD : %d == %d\n",rna->begin,rcpos);

                    if (rna->begin == rcpos) {
                        rna->to_delete = true;

//                        printf("RNA to delete LEAD %d\n",rna->begin);

                        for (auto prot : proteins_) {
                            for (auto rna2 = prot->rna_list_.begin(); rna2 != prot->rna_list_.end();) {
                                if (((*rna2)->begin == rna->begin) & ((*rna2)->leading_lagging == rna->leading_lagging)) {
                                    rna2 = prot->rna_list_.erase(rna2);

                                    prot->e -= rna->e;
                                }
                                else ++ rna2;
                            }

                            if (prot->rna_list_.empty()) {
//                                printf("Protein to delete %d : RNA %d\n",rna->begin,prot->protein_start);
                                prot->to_delete = true;
                            }
                        }
                    }
                }
            }
        }

        // LAG  : within [pos_1,pos_2+21]
//                printf("Remove RNA with %d to %d\n",pos_1,pos_2+21);
//                printf("Remove RNA with %d to %d\n",pos_1,pos_2+21);

        for (int32_t cpos = pos_1; cpos <= pos_2+21; cpos++) {
            int32_t rcpos = Utils::mod(cpos,length());
            for (auto rna : rnas_) {
//                printf("Process RNA %d %d\n",rna->begin,rna->leading_lagging);
                if (rna->leading_lagging == 1) {
//                    printf("To delete or NOT LAG : %d == %d\n",rna->begin,rcpos);

                    if (rna->begin == rcpos) {
                        rna->to_delete = true;

//                        printf("RNA to delete LAG %d\n",rna->begin);

                        for (auto prot : proteins_) {
                            for (auto rna2 = prot->rna_list_.begin(); rna2 != prot->rna_list_.end();) {
                                if (((*rna2)->begin == rna->begin) & ((*rna2)->leading_lagging == rna->leading_lagging)) {
                                    rna2 = prot->rna_list_.erase(rna2);
                                    prot->e -= rna->e;
                                }
                                else ++ rna2;
                            }

                            if (prot->rna_list_.empty()) {
//                                printf("Protein to delete %d : RNA %d\n", rna->begin, prot->protein_start);
                                prot->to_delete = true;
                            }
                        }
                    }
                }
            }
        }

        // 2. Recompute RNA with terminator
        // LEAD : within [pos_1-10,pos_2]
        for (int32_t cpos = pos_1-10; cpos <= pos_2; cpos++) {
            int32_t rcpos = Utils::mod(cpos,length());
            for (auto rna : rnas_) {
                if (rna->leading_lagging == 0) {
                    if (rna->end == rcpos) {
//                        printf("RECOMPUTE RNA %d => %d\n",rna->begin,rna->end);
                        rna->to_search_for_end = true;
                        rna->to_recompute = true;

                        for (auto prot : proteins_) {
                            for (auto rna2 = prot->rna_list_.begin(); rna2 != prot->rna_list_.end();) {
                                if (((*rna2)->begin == rna->begin) & ((*rna2)->leading_lagging == rna->leading_lagging)) {
                                    rna2 = prot->rna_list_.erase(rna2);
                                    prot->e -= rna->e;
                                }
                                else ++ rna2;
                            }

                            if (prot->rna_list_.empty())
                                prot->to_delete = true;
                        }
                    }
                }
            }
        }

        // LAG  : within [pos_1,pos_2+10]
        for (int32_t cpos = pos_1; cpos <= pos_2+10; cpos++) {
            int32_t rcpos = Utils::mod(cpos,length());
            for (auto rna : rnas_) {
                if (rna->leading_lagging == 1) {
                    if (rna->end == rcpos) {
//                        printf("RECOMPUTE RNA %d => %d\n",rna->begin,rna->end);
                        rna->to_search_for_end = true;
                        rna->to_recompute = true;

                        for (auto prot : proteins_) {
                            for (auto rna2 = prot->rna_list_.begin(); rna2 != prot->rna_list_.end();) {
                                if (((*rna2)->begin == rna->begin) & ((*rna2)->leading_lagging == rna->leading_lagging)) {
                                    rna2 = prot->rna_list_.erase(rna2);
                                    prot->e -= rna->e;
                                }
                                else ++ rna2;
                            }

                            if (prot->rna_list_.empty())
                                prot->to_delete = true;
                        }
                    }
                }
            }
        }

        // 2bis. Recompute if mutation within RNA
        for (auto rna : rnas_) {
            if (rna->leading_lagging == 0) {
                if (rna->begin < pos_1 && rna->end > pos_2) {
//                    printf("RECOMPUTE RNA MUTATION WITHIN %d => %d\n",rna->begin,rna->end);
                    rna->to_recompute = true;
                    rna->to_search_for_end = true;

                    for (auto prot : proteins_) {
                        for (auto rna2 = prot->rna_list_.begin(); rna2 != prot->rna_list_.end();) {
                            if (((*rna2)->begin == rna->begin) & ((*rna2)->leading_lagging == rna->leading_lagging)) {
                                rna2 = prot->rna_list_.erase(rna2);
                                prot->e -= rna->e;
                            }
                            else ++ rna2;
                        }

                        if (prot->rna_list_.empty())
                            prot->to_delete = true;
                    }
                }
            } else {
                if (rna->begin > pos_2 && rna->end < pos_1) {
//                    printf("RECOMPUTE RNA MUTATION WITHIN %d => %d\n",rna->begin,rna->end);
                    rna->to_recompute = true;
                    rna->to_search_for_end = true;

                    for (auto prot : proteins_) {
                        for (auto rna2 = prot->rna_list_.begin(); rna2 != prot->rna_list_.end();) {
                            if (((*rna2)->begin == rna->begin) & ((*rna2)->leading_lagging == rna->leading_lagging)) {
                                rna2 = prot->rna_list_.erase(rna2);
                                prot->e -= rna->e;
                            }
                            else ++ rna2;
                        }

                        if (prot->rna_list_.empty())
                            prot->to_delete = true;
                    }
                }
            }
        }

        // 3. Remove Protein with ShineDal+Start
        // LEAD : within [pos_1-12,pos_2]
//        printf("Search for Protein_STart (LEAD) between %d and %d\n",pos_1-12,pos_2);

        for (int32_t cpos = pos_1-12; cpos <= pos_2; cpos++) {
            int32_t rcpos = Utils::mod(cpos,length());
            for (auto prot : proteins_) {
                if (prot->leading_lagging == 0) {
                    if (prot->protein_start == rcpos) {
                        prot->to_delete = true;
                    }
                }
            }
        }

        // LAG  : within [pos_1,pos_2+12]
//        printf("Search for Protein_STart (LAG) between %d and %d\n",pos_1,pos_2+12);
        for (int32_t cpos = pos_1; cpos <= pos_2+12; cpos++) {
            int32_t rcpos = Utils::mod(cpos,length());
            for (auto prot : proteins_) {
                if (prot->leading_lagging == 1) {
                    if (prot->protein_start == rcpos) {
                        prot->to_delete = true;
                    }
                }
            }
        }

        // 4. Recompute Protein with Stop
        // LEAD : within [pos_1-3,pos_2]
//        printf("Search for Protein_Stop (LEAD) between %d and %d\n",pos_1-3,pos_2);

        for (int32_t cpos = pos_1-3; cpos <= pos_2; cpos++) {
            int32_t rcpos = Utils::mod(cpos,length());
            for (auto prot : proteins_) {
                if (prot->leading_lagging == 0) {
                    if (prot->protein_end == rcpos) {
                        prot->to_recompute = true;
                    }
                }
            }
        }

        // LAG  : within [pos_1,pos_2+3]
//        printf("Search for Protein_Stop (LEAD) between %d and %d\n",pos_1+3,pos_2);


        for (int32_t cpos = pos_1; cpos <= pos_2+3; cpos++) {
            int32_t rcpos = Utils::mod(cpos,length());
            for (auto prot : proteins_) {
                if (prot->leading_lagging == 1) {
                    if (prot->protein_end == rcpos) {
                        prot->to_recompute = true;
                    }
                }
            }
        }

        // 5. Recompute if mutation within protein
        for (auto prot : proteins_) {
            if (prot->leading_lagging == 0) {
                if (prot->protein_start < pos_1 && prot->protein_end > pos_2) {
//                    printf("RECOMPUTE PROTEIN MUTATION WITHIN %d => %d\n",
//                           prot->protein_start,prot->protein_end);

                    prot->to_delete = true;
                }
            } else {
                if (prot->protein_start > pos_2 && prot->protein_end < pos_1) {
//                    printf("RECOMPUTE PROTEIN MUTATION WITHIN %d => %d\n",
//                           prot->protein_start,prot->protein_end);

                    prot->to_delete = true;
                }
            }
        }

//        printf("======= END :: REMOVE RANGE ============ %d => %d\n",pos_1,pos_2);
//        display(false);
    }

    void SIMD_List_Metadata::update_range(int32_t pos) {
        update_range(pos,pos);
    }

    void SIMD_List_Metadata::update_range(int32_t pos_1, int32_t pos_2) {

//        printf("======= BEGIN :: UPDATE RANGE [%d => %d] ============\n",pos_1,pos_2);
//        display(false);

        // 1. Search promoter + mark as to compute
        // LEAD : within [pos_1-21,pos_2]
        for (int32_t cpos = pos_1-21; cpos <= pos_2; cpos++) {
            int32_t rcpos = Utils::mod(cpos,length());
            int8_t dist = is_promoter_leading(rcpos);
            if (dist <= 4) {// dist takes the hamming distance of the sequence from the consensus
                auto &strand = promoters_list_[LEADING];
                auto first = find_if(strand.begin(),
                                     strand.end(),
                                     [rcpos](promoterStruct &r) {
                                         return r.pos >= rcpos;
                                     });

                bool toAdd = false;
//                if (first != strand.end())
//                    printf("UPDATE_PROMOTER_LEADING : ALREADY FOUND %d => %d (%d)\n", first->pos, rcpos,first->to_delete);

                if (first == strand.end())
                    toAdd = true;

                if ((!first->to_delete) and first->pos != rcpos) {
                    toAdd = true;
                }

                if (toAdd) {
                    promoters_list_[LEADING].emplace(first, rcpos, dist, true);
//                    printf("UPDATE_PROMOTER_LEADING : NEW %d => %d\n", rcpos, dist);
                }
            }
        }

        // LAG  : within [pos_1,pos_2+21]
//                printf("Looking for prom LAG %d %d\n",pos_1,pos_2+21);
        for (int32_t cpos = pos_1; cpos <= pos_2+21; cpos++) {
            int32_t rcpos = Utils::mod(cpos,length());
            int8_t dist = is_promoter_lagging(rcpos);
            if (dist <= 4) {
                auto& strand = promoters_list_[LAGGING];

                auto first = find_if(strand.begin(),
                                     strand.end(),
                                     [rcpos](promoterStruct& r) { return r.pos <= rcpos; });

                bool toAdd = false;
//                if (first != strand.end())
//                    printf("UPDATE_PROMOTER_LAGGING : ALREADY FOUND %d => %d\n", rcpos, dist);

                if (first == strand.end())
                    toAdd = true;

                if ((!first->to_delete) and first->pos != rcpos) {
                    toAdd = true;
                }

                if (toAdd) {
                    promoters_list_[LAGGING].emplace(first, rcpos, dist, false);
//                    printf("UPDATE_PROMOTER_LAGGING : NEW %d => %d\n", rcpos, dist);
                }
            }
        }

//        printf("======= AFTER_PROMOTER_UPDATE :: UPDATE RANGE ============\n");
//        display(false);

        // 2. Within RNA, search terminator + update terminator position + delete all Protein on RNA
        // LEAD : within [pos_1-10,pos_2]
        Dna_SIMD* dna = indiv_->dna_;
        for (int32_t cpos = pos_1-10; cpos <= pos_2; cpos++) {
            int32_t rcpos = Utils::mod(cpos,length());
            int8_t term_dist_leading = 0;

            for (int t_motif_id = 0; t_motif_id < 4; t_motif_id++) {

                term_dist_leading +=
                        dna->get_lead(rcpos + t_motif_id) !=
                        dna->get_lead(rcpos - t_motif_id + 10)
                        ? 1 : 0;
            }

            if (term_dist_leading == 4) {
                for (auto rna : rnas_) {
                    if (!rna->to_delete) {
                        if (rna->leading_lagging == 0) {
                            if (rna->begin < rna->end) {
                                if (rna->begin + 22 < rcpos && rna->end > rcpos) {
                                    int32_t old_new = rna->end;

                                    rna->end = rcpos;

                                    rna->length = rna->end - rna->begin - 11;

                                    rna->to_recompute = true;

//                                    printf("UPDATE_RNA_LEADING : NEW %d => %d (PREV %d => %d)\n", rna->begin, old_new,
//                                           rna->begin, rna->end);
                                }
                            } else {
                                if ((rna->begin + 22 < rcpos && rcpos < length()) && (rna->end > rcpos && rcpos >= 0)) {
                                    int32_t old_new = rna->end;

                                    rna->end = rcpos;

                                    rna->length = dna->length() - rna->begin + rna->end;

                                    rna->to_recompute = true;

//                                    printf("UPDATE_RNA_LAGGING : NEW %d => %d (PREV %d => %d)\n", rna->begin, old_new,
//                                           rna->begin, rna->end);
                                }
                            }

                            if (rna->to_recompute) {
                                for (auto& prot : proteins_) {
                                    for (auto rna2 = prot->rna_list_.begin(); rna2 != prot->rna_list_.end();) {
                                        if (((*rna2)->begin == rna->begin) & ((*rna2)->leading_lagging == rna->leading_lagging)) {
                                            rna2 = prot->rna_list_.erase(rna2);
                                            prot->e -= rna->e;
                                        }
                                        else ++rna2;
                                    }

                                    if (prot->rna_list_.empty()) {
//                                        printf("Mark protein %d to delete\n",prot->protein_start);
                                        prot->to_delete = true;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // LAG  : within [pos_1,pos_2+10]
        for (int32_t cpos = pos_1-10; cpos <= pos_2; cpos++) {
            int32_t rcpos = Utils::mod(cpos,length());
            int8_t term_dist_lagging = 0;

            for (int32_t t_motif_id = 0; t_motif_id < 4; t_motif_id++) {
                term_dist_lagging +=
                        dna->get_lag(rcpos - t_motif_id) !=
                        dna->get_lag(rcpos + t_motif_id - 10)
                        ? 1 : 0;
            }

            if (term_dist_lagging == 4) {
                for (auto rna : rnas_) {
                    if (!rna->to_delete) {
                        if (rna->leading_lagging == 1) {
                            if (rna->begin > rna->end) {
                                if (rna->begin - 22 > rcpos && rna->end < rcpos) {
                                    int32_t old_new = rna->end;

                                    rna->end = rcpos;

                                    rna->length = rna->begin - rna->end - 21;

                                    rna->to_recompute = true;

//                                    printf("UPDATE_RNA_LAGGING: NEW %d => %d (PREV %d => %d)\n", rna->begin, old_new,
//                                           rna->begin, rna->end);
                                }
                            } else {
                                if ((rna->begin - 22 > rcpos && rcpos >= 0) && (rna->end < rcpos && rcpos < length())) {
                                    int32_t old_new = rna->end;

                                    rna->end = rcpos;

                                    rna->length = dna->length() + rna->begin - rna->end;

                                    rna->to_recompute = true;

//                                    printf("UPDATE_RNA_LAGGING : NEW %d => %d (PREV %d => %d)\n", rna->begin, old_new,
//                                           rna->begin, rna->end);
                                }
                            }

                            if (rna->to_recompute) {
                                for (auto prot : proteins_) {
                                    for (auto rna2 = prot->rna_list_.begin(); rna2 != prot->rna_list_.end();) {
                                        if (((*rna2)->begin == rna->begin) & ((*rna2)->leading_lagging == rna->leading_lagging)) {
                                            rna2 = prot->rna_list_.erase(rna2);
                                            prot->e -= rna->e;
                                        }
                                        else ++rna2;
                                    }

                                    if (prot->rna_list_.empty()) {
//                                        printf("Mark protein %d to delete\n",prot->protein_start);
                                        prot->to_delete = true;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

//        printf("======= AFTER_UPDATE_RNA :: UPDATE RANGE ============\n");
//        display(false);



        // 4. Within Protein, search for STOP codon + update Protein STOP position + mark as to retranslate
        // LEAD : within [pos_1-3,pos_2]
        for (int32_t cpos = pos_1-3; cpos <= pos_2; cpos++) {
            int32_t rcpos = Utils::mod(cpos,length());
            bool is_protein = false;

            for (int k = 0; k < 3; k++) {
                if (dna->get_lead(rcpos + k) ==
                    PROTEIN_END_LEAD[k]) {
                    is_protein = true;
                } else {
                    is_protein = false;
                    break;
                }
            }

            if (is_protein) {
                for (auto prot : proteins_) {
                    if (!prot->to_delete && prot->protein_end!=-1) {
                        if (prot->leading_lagging == 0) {
                            if (prot->protein_start < prot->protein_end) {
                                if (prot->protein_start < rcpos && prot->protein_end > rcpos) {
//                                    printf("UPDATE PROTEIN END %d => %d\n",prot->protein_end,rcpos);
                                    prot->protein_end = rcpos;
                                    prot->protein_length = (prot->protein_end - (prot->protein_start))/3;
                                    prot->is_init_ = true;
                                    prot->to_retranslate = true;
                                }
                            } else {
                                if ((prot->protein_start < rcpos && rcpos < length()) &&
                                    (prot->protein_end > rcpos && rcpos >= 0)) {
//                                    printf("UPDATE PROTEIN END %d => %d\n",prot->protein_end,rcpos);
                                    prot->protein_end = rcpos;
                                    prot->protein_length =
                                            (dna->length() - (prot->protein_start) + prot->protein_end)/3;
                                    prot->is_init_ = true;
                                    prot->to_retranslate = true;
                                }
                            }
                        }
                    }
                }

            }
        }

        // LAG  : within [pos_1,pos_2+3]
        for (int32_t cpos = pos_1-3; cpos <= pos_2; cpos++) {
            int32_t rcpos = Utils::mod(cpos,length());
            bool is_protein = false;

            for (int k = 0; k < 3; k++) {
                if (dna->get_lag(rcpos - k) == PROTEIN_END_LAG[k]) {
                    is_protein = true;
                } else {
                    is_protein = false;
                    break;
                }
            }

            if (is_protein) {
                for (auto prot : proteins_) {
                    if (!prot->to_delete && prot->protein_end!=-1) {
                        if (prot->leading_lagging == 1) {
                            if (prot->protein_start > prot->protein_end) {
                                if (prot->protein_start > rcpos && prot->protein_end < rcpos) {
//                                    printf("UPDATE PROTEIN END %d => %d\n",prot->protein_end,rcpos);
                                    prot->protein_end = rcpos;
                                    prot->protein_length = ((prot->protein_start) - prot->protein_end)/3;
                                    prot->to_retranslate = true;
                                    prot->is_init_ = true;
                                }
                            } else {
                                if ((prot->protein_start < rcpos && rcpos >= 0) &&
                                    (prot->protein_end < rcpos && rcpos < length())) {
//                                    printf("UPDATE PROTEIN END %d => %d\n",prot->protein_end,rcpos);
                                    prot->protein_end = rcpos;
                                    prot->protein_length = ((prot->protein_start) +
                                                           dna->length() - prot->protein_end)/3;
                                    prot->to_retranslate = true;
                                    prot->is_init_ = true;
                                }
                            }
                        }
                    }
                }
            }
        }

        // 3. Within RNA, search for ShineDal+START codon + add Protein to recompute
        // LEAD : within [pos_1-12,pos_2]
        for (int32_t cpos = pos_1-12; cpos <= pos_2; cpos++) {
            int32_t rcpos = Utils::mod(cpos,length());
            bool start = false;
            int32_t k_t;

            for (int32_t k = 0; k < 9; k++) {
                k_t = k >= 6 ? k + 4 : k;

                if (dna->get_lead(rcpos + k_t) == SHINE_DAL_SEQ_LEAD[k]) {
                    start = true;
                } else {
                    start = false;
                    break;
                }
            }

            if (start) {
                int32_t glob_prot_idx = -1;
                for (auto rna : rnas_) {
                    if (!rna->to_delete)
                        if (rna->leading_lagging == 0) {
                                if (rna->begin + 22 < rcpos && rna->end > rcpos) {
                                    int mrcpos = Utils::mod(rcpos+13,length());
                                    auto prot = std::find_if(proteins_.begin(),proteins_.end(),
                                                             [mrcpos](const pProtein* q)
                                                             {return q->protein_start == mrcpos;});

                                    bool toAdd=false;
                                    if ((prot == proteins_.end()))
                                        toAdd = true;
                                    else if ((!(*prot)->to_delete) and ((*prot)->leading_lagging != rna->leading_lagging)) {
                                        toAdd = true;
                                    } else if ((!(*prot)->to_delete) and ((*prot)->protein_start == rcpos+13) and
                                            ((*prot)->leading_lagging != rna->leading_lagging)) {
                                        (*prot)->e += rna->e;
                                        (*prot)->rna_list_.push_back(rna);
                                    }

                                    if (toAdd) {
                                        glob_prot_idx = proteins_count();
                                        set_proteins_count(proteins_count() + 1);

//                                        printf("]]]]]]]]]]]] ADD == Protein (IF) %d\n", rcpos);

//                                        printf("]]]]]]]]]]]] ADD == Protein (IF) %d -- RNA %d %d\n",
//                                               Utils::mod(rcpos+13,length()),
//                                               rna->begin, rna->end);


                                        protein_add(glob_prot_idx, new pProtein(Utils::mod(rcpos+13,length()), -1, -1,
                                                                                rna->leading_lagging,
                                                                                rna->e, rna, false
                                        ));
                                    }
                                }
                            } else {
                                if ((rna->begin + 22 < rcpos && rcpos < length()) &&
                                    (rna->end > rcpos && rcpos >= 0)) {
                                    int mrcpos = Utils::mod(rcpos-13,length());

                                    auto prot = std::find_if(proteins_.begin(),proteins_.end(),
                                                             [mrcpos](const pProtein* q)
                                                             {return q->protein_start == mrcpos;});


                                    bool toAdd=false;
                                    if ((prot == proteins_.end()))
                                        toAdd = true;
                                    else if ((!(*prot)->to_delete) and  ((*prot)->leading_lagging != rna->leading_lagging)) {
                                        toAdd = true;
                                    } else if ((!(*prot)->to_delete) and ((*prot)->protein_start == rcpos-13) and
                                               ((*prot)->leading_lagging != rna->leading_lagging)) {
                                        (*prot)->e += rna->e;
                                        (*prot)->rna_list_.push_back(rna);
                                    }

                                    if (toAdd) {
                                        glob_prot_idx = proteins_count();
                                        set_proteins_count(proteins_count() + 1);

//                                        printf("]]]]]]]]]]]] ADD == Protein (ELSE) %d -- RNA %d %d\n", Utils::mod(rcpos-13,length()),
//                                               rna->begin, rna->end);

                                        protein_add(glob_prot_idx, new pProtein(Utils::mod(rcpos-13,length()), -1, -1,
                                                                                rna->leading_lagging,
                                                                                rna->e, rna, false
                                        ));
                                    }
                                }
                            }
                        }

            }
        }

        // LAG  : within [pos_1,pos_2+12]
        for (int32_t cpos = pos_1; cpos <= pos_2+12; cpos++) {
            int32_t rcpos = Utils::mod(cpos,length());
            bool start = false;
            int32_t k_t;

            for (int k = 0; k < 9; k++) {
                k_t = k >= 6 ? k + 4 : k;

                if (dna->get_lag(rcpos - k_t) == SHINE_DAL_SEQ_LAG[k]) {
                    start = true;
                } else {
                    start = false;
                    break;
                }
            }

            if (start) {
                int32_t glob_prot_idx = -1;
                for (auto rna : rnas_) {
                    if (!rna->to_delete)
                        if (rna->leading_lagging == 0) {
                            if (rna->begin - 22 > rcpos && rna->end < rcpos) {
                                int mrcpos = Utils::mod(rcpos+13,length());

                                auto prot = std::find_if(proteins_.begin(), proteins_.end(),
                                                         [mrcpos](const pProtein *q) {
                                                             return q->protein_start == mrcpos;
                                                         });
                                bool toAdd = false;
                                if ((prot == proteins_.end()))
                                    toAdd = true;
                                else if ((!(*prot)->to_delete) and ((*prot)->leading_lagging != rna->leading_lagging)) {
                                    toAdd = true;
                                } else if ((!(*prot)->to_delete) and ((*prot)->protein_start == rcpos+13) and
                                           ((*prot)->leading_lagging != rna->leading_lagging)) {
                                    (*prot)->e += rna->e;
                                    (*prot)->rna_list_.push_back(rna);
                                }

                                if (toAdd) {
                                    glob_prot_idx = proteins_count();
                                    set_proteins_count(proteins_count() + 1);
//                                    printf("]]]]]]]]]]]] ADD == Protein (LEAD 2) %d -- RNA %d %d\n", Utils::mod(rcpos+13,length()),
//                                           rna->begin, rna->end);
                                    protein_add(glob_prot_idx, new pProtein(Utils::mod(rcpos+13,length()), -1, -1,
                                                                            rna->leading_lagging,
                                                                            rna->e, rna, false
                                    ));
                                }
                            }
                        } else {
                            if ((rna->begin - 22 > rcpos && rcpos >= 0) &&
                                (rna->end < rcpos && rcpos < length())) {
                                int mrcpos = Utils::mod(rcpos-13,length());

                                auto prot = std::find_if(proteins_.begin(), proteins_.end(),
                                                         [mrcpos](const pProtein *q) {
                                                             return q->protein_start == mrcpos;
                                                         });
                                bool toAdd = false;
                                if ((prot == proteins_.end()))
                                    toAdd = true;
                                else if ((!(*prot)->to_delete) and ((*prot)->leading_lagging != rna->leading_lagging)) {
                                    toAdd = true;
                                } else if (((*prot)->protein_start == rcpos-13)
                                           and ((*prot)->leading_lagging == rna->leading_lagging)) {

                                } else if ((!(*prot)->to_delete) and ((*prot)->protein_start == rcpos-13) and
                                           ((*prot)->leading_lagging != rna->leading_lagging)) {
                                    (*prot)->e += rna->e;
                                    (*prot)->rna_list_.push_back(rna);
                                }

                                if (toAdd) {
                                    glob_prot_idx = proteins_count();
                                    set_proteins_count(proteins_count() + 1);
//                                    printf("]]]]]]]]]]]] ADD == Protein (LAG 2) %d -- RNA %d %d %d\n", Utils::mod(rcpos-13,length()),
//                                           rna->begin, rna->end, rna->leading_lagging);
                                    protein_add(glob_prot_idx, new pProtein(Utils::mod(rcpos-13,length()), -1, -1,
                                                                            rna->leading_lagging,
                                                                            rna->e, rna, false
                                    ));
                                }
                            }
                        }

                }

            }
        }

//        printf("======= END :: UPDATE RANGE ============ %d %d\n",pos_1,pos_2);
//        display(false);
    }

    void SIMD_List_Metadata::display(bool check_to_delete) {
        bool error = false;
        printf("[SIMD] -- %d -- [ BEGIN -------------------------------------------------------]\n",indiv_->indiv_id);
        printf("[SIMD] -- Length %d\n",length());
        printf("[SIMD] -- Promoters Start : ");

        for (auto &strand: promoters_list_)
            for (auto &prom: strand) {
            printf("[%d (%d) => %d %d]",prom.pos,prom.leading_or_lagging,prom.to_compute,prom.to_delete);
            if ((prom.pos >= length())) {
                if (prom.to_delete && check_to_delete) {
                    error = true;
                    printf("\nPromoters ERROR : %d %d %d\nPromoters Start : ", prom.pos, length(), prom.to_delete);
                }
            }

        }
        printf("\n[SIMD] -- RNA Start : ");

        for (auto rna : rnas_) {
            printf("[%d => %d (R %d D %d LoL %d SfT %d)]",rna->begin,rna->end,rna->to_recompute,rna->to_delete,
                   rna->leading_lagging,rna->to_search_for_end);
            if ((rna->begin >= length()))
                if ((rna->to_delete && check_to_delete)) {
                    error = true;
                    printf("\nRNA ERROR : %d %d %d\nRNA Start : ", rna->begin, length(), rna->to_delete);
                }

        }

        printf("\n[SIMD] -- Proteins : ");
        for (auto prot : proteins_) {
            printf("[%d (%d) => %d == %lf %lf %lf] %d : ",prot->protein_start,prot->to_delete,prot->protein_end,
            prot->w,prot->h,prot->m,prot->rna_list_.size());
            for (auto rna : prot->rna_list_) {
                printf("[%d => %d] ",rna->begin,rna->end);
            }
            printf("\n");
            if (prot->protein_start >= length())
                if ((prot->to_delete && check_to_delete)) {
                    error = true;
                    printf("\nProtein ERROR : %d %d %d\nProtein Start : ", prot->protein_start, length(), prot->to_delete);
                }
        }
        printf("\n[SIMD] -- %d -- [ END   -------------------------------------------------------]\n",indiv_->indiv_id);

        if (error) {
            printf("ERROR !!!!\n");exit(11);
        }

        if (check_to_delete) {
            for (auto rna : rnas_) {
                if (rna->to_delete)
                    printf("RNA -- Must has been deleted\n");
            }

            for (auto prot : proteins_) {
                if (prot->to_delete)
                    printf("Proteins -- Must has been deleted\n");
            }
        }
    }

    void SIMD_List_Metadata::update_metadata() {
//        printf("Start processing Individual %d\n",indiv_->indiv_id);
        // 0. Cleanup Promoter/Rna/Protein mark as to delete
        cleanup();
//        display(true);
//        proteins_print();


//        printf("%d -- Compute promoters\n",indiv_->indiv_id);
        // 1. Compute new promoter
        compute_promoters();

//        display(true);

        // 2. Recompute RNA mark as to_recompute
//        printf("%d -- Search RNAs Terminator\n",indiv_->indiv_id);
        search_for_rna_terminator();
        cleanup();
//        printf("%d -- Compute RNAs\n",indiv_->indiv_id);
        recompute_rna();

//        display(true);
//        printf("%d -- Compute proteins\n",indiv_->indiv_id);
        // 3. Recompute Protein mark as to_recompute
        recompute_proteins();
//        display(true);


//        proteins_print();

//        display(true);
//        printf("%d -- Translate proteins\n",indiv_->indiv_id);
        // 4. Retranslate Protein mark as to_retranslate
        retranslate_proteins();

//        proteins_print();
//        display(true);
        // 5. If the proteome change (at least one new protein or a protein has changed), recompute phenotype+fitness
        //if (recompute_phenotype)
//        display(true);
//        printf("%d -- Individual has been computed\n",indiv_->indiv_id);

    }


    void SIMD_List_Metadata::cleanup() {
        // 1. Promoters
        for (auto &strand: promoters_list_)
            for (auto it = strand.begin(); it != strand.end();) {
                if ((*it).to_delete) it = strand.erase(it);
                else ++it;
            }

        // 2. RNA
        for (std::list<pRNA*>::iterator it_rna = rnas_.begin(); it_rna != rnas_.end();) {
            if ((*(it_rna))->to_delete) {
                delete (*(it_rna));
                it_rna = rnas_.erase(it_rna);
            }
            else ++it_rna;
        }

        // 3. Proteins
        for (std::list<pProtein*>::iterator it_protein = proteins_.begin();
                it_protein != proteins_.end(); ) {
            if ((*(it_protein))->to_delete) {
                delete (*(it_protein));
                it_protein = proteins_.erase(it_protein);
            } else ++it_protein;
        }
    }

    void SIMD_List_Metadata::compute_promoters() {
        Dna_SIMD* dna = indiv_->dna_;
        for (auto &strand: promoters_list_)
            for (auto &prom: strand) {
                if (prom.to_compute) {
                    recompute_phenotype = true;

                    Dna_SIMD *dna = indiv_->dna_;
                    int32_t dna_length = dna->length();
                    int prom_pos;
                    bool lead_lag;
                    double prom_error;

                    prom_pos = prom.pos;
                    lead_lag = prom.leading_or_lagging;
                    prom_error = fabs(
                            ((float) prom.error));

                    if (lead_lag) {
                        /* Search for terminators */
                        int cur_pos =
                                Utils::mod(prom_pos + 22,dna_length);
                        int start_pos = cur_pos;

                        bool terminator_found = false;
                        bool no_terminator = false;
                        int term_dist_leading = 0;

                        int loop_size = 0;

                        while (!terminator_found) {
                            for (int t_motif_id = 0; t_motif_id < 4; t_motif_id++) {

                                term_dist_leading +=
                                        dna->get_lead(cur_pos + t_motif_id) !=
                                        dna->get_lead(cur_pos - t_motif_id + 10)
                                        ? 1 : 0;
                            }

                            if (term_dist_leading == 4) {
                                terminator_found = true;
                            }
                            else {
                                cur_pos = Utils::mod(cur_pos + 1,dna_length);

                                term_dist_leading = 0;
                                if (cur_pos == start_pos) {
                                    no_terminator = true;
                                    terminator_found = true;
                                }
                            }

                            loop_size++;

                        }

                        if (!no_terminator) {

                            int32_t rna_end =
                                    cur_pos;
                            int32_t rna_length = 0;

                            rna_length = loop_size + TERM_SIZE -1;

                            if (rna_length > 0) {


                                int glob_rna_idx = -1;
                                glob_rna_idx = rna_count_++;
                                rna_add(glob_rna_idx,
                                                                                   prom_pos,
                                                                                   rna_end,
                                                                                   !lead_lag,
                                                                                   1.0 -
                                                                                   prom_error /
                                                                                   5.0, rna_length);
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

                            if (term_dist_lagging == 4) {
                                terminator_found = true;
                            }
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


                            int32_t rna_end = Utils::mod(cur_pos - 10,dna_length);

                            int32_t rna_length = 0;

                            rna_length = loop_size + TERM_SIZE -1;

                            if (rna_length >= 0) {
                                int glob_rna_idx = -1;
                                glob_rna_idx = rna_count_++;


                                 rna_add(glob_rna_idx,
                                                                                   prom_pos,
                                                                                   rna_end,
                                                                                   !lead_lag,
                                                                                   1.0 -
                                                                                   prom_error /
                                                                                   5.0, rna_length);
                            }

                        }
                    }

                }
            }
    }

    void SIMD_List_Metadata::search_for_rna_terminator() {
//        display(true);
        Dna_SIMD* dna = indiv_->dna_;
        for (auto rna: rnas_)
            if (rna->to_search_for_end) {
                recompute_phenotype = true;
                Dna_SIMD *dna = indiv_->dna_;
                int32_t dna_length = dna->length();
                int prom_pos;
                bool lead_lag;
                double prom_error;

                prom_pos = rna->begin;
                lead_lag = rna->leading_lagging;

//                printf("Search for new terminator %d (%d)\n",prom_pos,lead_lag);

                if (lead_lag == 0) {
                    /* Search for terminators */
                    int cur_pos =
                            Utils::mod(prom_pos + 22,dna_length);
                    int start_pos = cur_pos;

                    bool terminator_found = false;
                    bool no_terminator = false;
                    int term_dist_leading = 0;

                    int loop_size = 0;

                    while (!terminator_found) {
                        for (int t_motif_id = 0; t_motif_id < 4; t_motif_id++) {

                            term_dist_leading +=
                                    dna->get_lead(cur_pos + t_motif_id) !=
                                    dna->get_lead(cur_pos - t_motif_id + 10)
                                    ? 1 : 0;
                        }

                        if (term_dist_leading == 4) {
                            terminator_found = true;
                        }
                        else {
                            cur_pos = Utils::mod(cur_pos + 1,dna_length);

                            term_dist_leading = 0;
                            if (cur_pos == start_pos) {
                                no_terminator = true;
                                terminator_found = true;
                            }
                        }

                        loop_size++;

                    }

                    if (!no_terminator) {

                        int32_t rna_end =
                                cur_pos;

                        int32_t rna_length = 0;

                        rna_length = loop_size + TERM_SIZE -1;

                        if (rna_length > 0) {
//                            printf("Update terminator %d OLD %d Begin %d\n",rna_end,rna->end,rna->begin);

                            rna->end = rna_end;
                            rna->length = rna_length;


                            for (auto prot : proteins_) {
                                for (auto rnax : prot->rna_list_) {
                                    if (rnax->begin == rna->begin && rna->leading_lagging == rnax->leading_lagging) {
//                                        printf("Protein still within RNA ? %d => %d RNA %d => %d\n",prot->protein_start,
//                                               prot->protein_end,rna->begin,rna->end);

                                        prot->to_delete = true;
                                    }
                                }
                            }

                        } else
                            rna->to_delete = true;
                    }else
                        rna->to_delete = true;
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

                        if (term_dist_lagging == 4) {
                            terminator_found = true;
                        }
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

//                    printf("New terminator found or not ? %d : %d (%d) %d\n",
//                           no_terminator,cur_pos,rna->begin,length());

                    if (!no_terminator) {


                        int32_t rna_end = Utils::mod(cur_pos - 10,dna_length);

                        int32_t rna_length = 0;
                        rna_length = loop_size + TERM_SIZE -1;


                        if (rna_length > 0) {

//                            printf("Update terminator %d OLD %d Begin %d\n",rna_end,rna->end,rna->begin);

                            rna->end = rna_end;
                            rna->length = rna_length;


                            for (auto prot : proteins_) {
                                for (auto rnax : prot->rna_list_) {
                                    if (rnax->begin == rna->begin && rna->leading_lagging == rnax->leading_lagging) {
//                                        printf("Protein still within RNA ? %d => %d RNA %d => %d\n",prot->protein_start,
//                                        prot->protein_end,rna->begin,rna->end);

                                        prot->to_delete = true;
                                    }
                                }
                            }
                        } else
                            rna->to_delete = true;

                    }else
                        rna->to_delete = true;
                }

            }
//        display(true);
    }

    void SIMD_List_Metadata::recompute_rna() {
        Dna_SIMD* dna = indiv_->dna_;
        int32_t dna_length = length();
//        display(true);
        int rna_idx = 0;
        for (auto rna : rnas_) {
//            printf("Compute RNA %d %d %d\n",rna_idx,rna->to_recompute,rna->to_delete);
            rna_idx++;

            if (rna->to_recompute) {
                rna->reset();

                recompute_phenotype = true;

                if (rna->is_init_) {
                    int c_pos = rna->begin;

                    if (rna->length >= 21) {
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

                        int loop_size = 0;
//                        printf("RNA Begin loop %d => %d :: %d\n",rna->begin,rna->end,loop_size);
                        while (loop_size < rna->length) {
                            bool start = false;
                            int k_t;

                            if (rna->leading_lagging ==
                                0) {
                                // Search for Shine Dalgarro + START codon on LEADING
                                for (int k = 0; k < 9; k++) {
                                    k_t = k >= 6 ? k + 4 : k;

                                    if (indiv_->dna_->data_[Utils::mod((c_pos + k_t),dna_length)] ==
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

                                    if (indiv_->dna_->data_[Utils::mod((c_pos - k_t),dna_length)] ==
                                        SHINE_DAL_SEQ_LAG[k]) {
                                        start = true;
                                    } else {
                                        start = false;
                                        break;
                                    }
                                }
                            }

                            if (start) {
                                int32_t glob_prot_idx = proteins_count();
                                set_proteins_count(proteins_count() + 1);

                                int start_protein_pos = rna->leading_lagging == 0 ?
                                                        c_pos +
                                                        13 :
                                                        c_pos -
                                                        13;
                                start_protein_pos = Utils::mod(start_protein_pos,dna_length);

                                auto prot = std::find_if(proteins_.begin(), proteins_.end(),
                                                         [start_protein_pos](const pProtein *q) {
                                                             return q->protein_start == start_protein_pos;
                                                         });
                                bool toAdd = false;
                                if ((prot == proteins_.end()))
                                    toAdd = true;
                                else if ((!(*prot)->to_delete) and ((*prot)->leading_lagging != rna->leading_lagging)) {
                                    toAdd = true;
                                } else if ((!(*prot)->to_delete) and ((*prot)->protein_start == start_protein_pos) and
                                           ((*prot)->leading_lagging != rna->leading_lagging)) {
                                    (*prot)->e += rna->e;
                                    (*prot)->rna_list_.push_back(rna);
                                }

                                if (toAdd) {
                                    protein_add(glob_prot_idx, new pProtein(start_protein_pos, -1, -1,
                                                                            rna->leading_lagging,
                                                                            rna->e, rna, false
                                    ));

//                                    printf("Add proteins %d -- %d %d ==> ",start_protein_pos,rna->begin,rna->end);
//                                    if (proteins_.front()->protein_start == start_protein_pos) {
//                                        for (auto rna : proteins_.front()->rna_list_) {
//                                            printf("[ %d %d ] ",
//                                                   proteins_.front()->rna_list_.front()->begin,
//                                                   proteins_.front()->rna_list_.front()->end);
//                                        }
//                                    }
//                                    printf("\n");
                                }

                            }

                            if (rna->leading_lagging == 0) {
                                c_pos++;
                                c_pos =
                                        c_pos >= dna_length ? c_pos - dna_length
                                                            : c_pos;
                            } else {
                                c_pos--;
                                c_pos = c_pos < 0 ? dna_length + c_pos : c_pos;
                            }
                            loop_size++;
                        }
                    }

                }
            }
        }
//        display(true);
    }

    void SIMD_List_Metadata::recompute_proteins() {
        Dna_SIMD* dna = indiv_->dna_;
        int32_t dna_length = length();

        for (auto& prot : proteins_) {
            if (prot->to_recompute || prot->protein_end == -1) {
                recompute_phenotype = true;

                pRNA* rna = prot->rna_list_.front();

                // Search for proteins STOP codon
                int start_prot = prot->protein_start;
                int start_protein_pos = rna->leading_lagging == 0 ?
                                        start_prot :
                                        start_prot ;
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
                        length =  prot->protein_start - rna->end;
                    } else {
                        length =  prot->protein_start +
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

//                    printf("RNA Begin %d Transcribed start %d Start Prot %d\n",
//                           rna->begin,transcribed_start,start_prot);

                    if (transcribed_start <= start_prot) {
                        j = start_prot - transcribed_start;
                    } else {
                        j = dna_length -
                            transcribed_start +
                            start_prot;

                    }
                } else {
                    transcribed_start = rna->begin - 22;
//                    printf("A -- RNA Begin %d Transcribed start %d Start Prot %d\n",
//                           rna->begin,transcribed_start,start_prot);

                    transcribed_start = transcribed_start < 0 ?
                                        dna_length + transcribed_start
                                                              : transcribed_start;

//                    printf("B -- RNA Begin %d Transcribed start %d Start Prot %d\n",
//                           rna->begin,transcribed_start,start_prot);

                    if (transcribed_start >=
                        start_prot) {
                        j = transcribed_start -
                            start_prot;
                    } else {
                        j = transcribed_start +
                            dna_length - start_prot;
                    }
                }

//                printf("Searching for protein END J %d Lenght %d TS %d ProtStart %d RNA %d => %d\n",
//                       j,rna->length,transcribed_start,prot->protein_start,rna->begin,rna->end);

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



//                        printf("LOOKing for start prot at %d : %d\n",start_protein_pos,is_protein);
//

                        if (is_protein && j > 1) {
                            int prot_length = -1;

                            if (start_prot < t_k) {
                                prot_length = t_k -
                                              (start_prot);
                            } else {
                                prot_length = dna_length -
                                              (start_prot) + t_k;
                            }

//                            printf("ADD PROT LEAD %d %d %d %d\n",start_prot,prot_length,t_k,j);
                            if (prot_length >= 3) {

                                prot->protein_length = prot_length/3;
                                prot->protein_end = t_k;

                                prot->is_init_ = true;
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

//                        printf("Search for prot end starting at %d : IsP %d J %d -- L %d\n",start_protein_pos,
//                               is_protein,j,
//                               rna->length);


                        if (is_protein && j > 1) {
                            int prot_length = -1;
                            if (start_prot > t_k) {
                                prot_length =
                                        (start_prot -
                                        t_k);
                            } else {
                                prot_length =
                                        (start_prot) +
                                        dna_length - t_k;
                            }

//                            printf("ADD PROT %d %d %d %d\n",start_prot,prot_length,t_k,j);
                            if (prot_length >= 3) {
                                prot->protein_length = prot_length/3;
                                prot->protein_end = t_k;

                                prot->is_init_ = true;
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

    void SIMD_List_Metadata::retranslate_proteins() {
        Dna_SIMD* dna = indiv_->dna_;
        int32_t dna_length = length();

        for (auto prot : proteins_) {
//            printf("TRANSLATING PROTEIN %d -- %d :: %d %d\n",prot->protein_start,prot->protein_length,
//                   prot->to_retranslate,prot->is_init_);
            if (prot->to_retranslate && prot->is_init_) {
                recompute_phenotype = true;

                if (prot->is_init_) {
                    int c_pos = prot->protein_start, t_pos;
                    int end_pos = prot->protein_end;
                    if (prot->leading_lagging ==
                        0) {
                        //c_pos += 13;
                        end_pos -= 3;

                        c_pos =
                                c_pos >= dna_length ? c_pos - dna_length
                                                            : c_pos;
                        end_pos = end_pos < 0 ? dna_length + end_pos : end_pos;
                    } else {
                        //c_pos -= 13;
                        end_pos += 3;

                        end_pos =
                                end_pos >= dna_length ? end_pos - dna_length
                                                              : end_pos;
                        c_pos = c_pos < 0 ? dna_length + c_pos : c_pos;
                    }

                    int8_t value = 0;
                    int8_t codon_list[64] = {};
                    int8_t codon_idx = 0;
                    int32_t count_loop = 0;

                    if (prot->leading_lagging ==
                        0) {
                        // LEADING

                        while (count_loop <
                               prot->protein_length &&
                               codon_idx < 64) {
                            value = 0;
                            for (int8_t i = 0; i < 3; i++) {
                                t_pos = c_pos + i;
                                if (indiv_->dna_->get_lead(t_pos) ==
                                    '1')
                                    value += 1 << (CODON_SIZE - i - 1);
                            }
//                            printf("Codon %d  : %d\n",codon_idx,value);
                            codon_list[codon_idx] = value;

                            codon_idx++;

                            count_loop++;
                            c_pos += 3;
                            c_pos =
                                    c_pos >= dna_length ? c_pos - dna_length
                                                                : c_pos;
                        }
                    } else {
                        // LAGGING
                        while (count_loop <
                               prot->protein_length &&
                               codon_idx < 64) {
                            value = 0;
                            for (int8_t i = 0; i < 3; i++) {
                                t_pos = c_pos - i;
                                if (indiv_->dna_->get_lag(t_pos) !=
                                    '1')
                                    value += 1 << (CODON_SIZE - i - 1);
                            }
                            codon_list[codon_idx] = value;
                            codon_idx++;

                            count_loop++;

                            c_pos -= 3;
                            c_pos = c_pos < 0 ? c_pos + dna_length : c_pos;
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
                            (indiv_->w_max_ - W_MIN) *
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
        }


        std::map<int32_t, pProtein *> lookup;

        protein_begin();
        //((SIMD_List_Metadata*)internal_simd_struct[indiv_id]->metadata_)->proteins_print();

        for (int protein_idx = 0; protein_idx <
                                  (int) proteins_count(); protein_idx++) {
            {
                //internal_simd_struct[indiv_id]->metadata_->protein_begin();

                pProtein* prot  = protein_next();
//                printf("Protein %d Indiv %d : %p %p -- %lf %lf %lf\n",protein_idx,indiv_->indiv_id,prot,
//                       proteins(protein_idx),prot->m,prot->w,prot->h*prot->e);
                if (prot->is_init_ && prot->leading_lagging==0) {
                    if (lookup.find(prot->protein_start) ==
                        lookup.end()) {
                        lookup[prot->protein_start] = prot;
                    } else {
                        lookup[prot->protein_start]->e += prot->e;
                        prot->is_init_ = false;
//                        printf("Protein %d is already there, fuz it with another one %f : %f\n",prot->protein_start,
//                               prot->e,lookup[prot->protein_start]->e);
                    }
                }
            }
        }

        lookup.clear();

        protein_begin();
        //proteins_print();
        protein_begin();

        for (int protein_idx = 0; protein_idx <
                                  (int) proteins_count(); protein_idx++) {
            {
                //internal_simd_struct[indiv_id]->metadata_->protein_begin();

                pProtein* prot  = protein_next();
//                printf("Protein %d Indiv %d : %p %p -- %lf %lf %lf\n",protein_idx,indiv_->indiv_id,prot,
//                       proteins(protein_idx),prot->m,prot->w,prot->h*prot->e);

                if (prot->is_init_ && prot->leading_lagging==1) {
                    if (lookup.find(prot->protein_start) ==
                        lookup.end()) {
                        lookup[prot->protein_start] = prot;
                    } else {
                        lookup[prot->protein_start]->e += prot->e;
                        prot->is_init_ = false;
//                        printf("Protein %d is already there, fuz it with another one %f : %f\n",prot->protein_start,
//                               prot->e,lookup[prot->protein_start]->e);
                    }
                }
            }
        }
    }

    void SIMD_List_Metadata::update_positions(int pos_after, bool insert_or_remove, int length_diff) {
//        printf("Update positions %d %d %d -- %d\n",pos_after,insert_or_remove,length_diff,length());
//        printf("%d -- UPDATE_POSITIONS : Begin\n",indiv_->indiv_id);


        if (insert_or_remove) {
            // INSERT
            // UPDATE Promoters
            for (auto &strand: promoters_list_) {
                for (auto &prom: strand) {
                    if (prom.pos >= pos_after) {
//                        printf("INSERT -- Update Promoters from %d to %d\n",prom.pos,prom.pos+length_diff);
                        prom.pos += length_diff;
                    }
                }
            }

            promoters_list_[LEADING].sort();
            promoters_list_[LAGGING].sort();
            promoters_list_[LAGGING].reverse();

            // UPDATE Rna
            for (auto rna : rnas_) {
                if (rna->begin >= pos_after) {
//                    printf("INSERT -- Update RNA from %d to %d\n",rna->begin,rna->begin+length_diff);
                    rna->begin += length_diff;
                }
                if (rna->end >= pos_after)
                    rna->end   += length_diff;
            }

            // UPDATE Proteins
            for (auto prot : proteins_) {
                if (prot->protein_start >= pos_after) {
//                    printf("INSERT -- Update Proteins from %d to %d\n",prot->protein_start,prot->protein_start+length_diff);
                    prot->protein_start += length_diff;
                }
                if (prot->protein_end >= pos_after)
                    prot->protein_end   += length_diff;
            }
        } else {
            // REMOVE
//
//            printf("Promoters list : ");
//            for (auto &strand: promoters_list_)
//                for (auto &prom: strand) {
//                    printf("%d ",prom.pos);
//                }
//            printf("\n");

            // UPDATE Promoters
            for (auto &strand: promoters_list_) {
                for (auto &prom: strand) {
                    if (prom.pos >= pos_after) {
//                        printf("DELETE -- Update Promoters from %d to %d\n",prom.pos,prom.pos-length_diff);
                        prom.pos -= length_diff;
                    }
                }
            }

            promoters_list_[LEADING].sort();
            promoters_list_[LAGGING].sort();
            promoters_list_[LAGGING].reverse();
            // UPDATE Rna
            for (auto rna : rnas_) {
                if (rna->begin >= pos_after) {
//                    printf("DELETE -- Update RNA from %d to %d\n",rna->begin,rna->begin-length_diff);
                    rna->begin -= length_diff;
                }
                if (rna->end >= pos_after)
                    rna->end   -= length_diff;
            }

            // UPDATE Proteins
            for (auto prot : proteins_) {
                if (prot->protein_start >= pos_after) {
//                    printf("DELETE -- Update Proteins from %d to %d\n",prot->protein_start,prot->protein_start-length_diff);
                    prot->protein_start -= length_diff;
                }
                if (prot->protein_end >= pos_after)
                    prot->protein_end   -= length_diff;
            }
        }

//        printf("%d -- UPDATE_POSITIONS : End\n",indiv_->indiv_id);

        //display(false);

    }

    void SIMD_List_Metadata::update_positions(int pos_after, int pos_before, int8_t add_or_reduce, int length_diff, int8_t mark) {
//        printf("Update positions %d %d %d -- %d -- %d\n",pos_after,pos_before,length_diff,length(),mark);
//        printf("%d -- UPDATE_POSITIONS : Begin\n",indiv_->indiv_id);

        if (add_or_reduce == 1) {
            // INCREASE
            // UPDATE Promoters
            for (auto &strand: promoters_list_) {
                for (auto &prom: strand) {
                    if (prom.pos >= pos_after && prom.pos < pos_before) {
                        if (mark == -1 || prom.transloc_move == mark) {
//                            printf("INSERT -- Update Promoters from %d to %d\n", prom.pos, prom.pos + length_diff);
                            prom.pos = prom.pos+length_diff;
                        }
                    }
                }
            }

            // UPDATE Rna
            for (auto rna : rnas_) {
                if (rna->begin >= pos_after && rna->begin < pos_before) {
                    if (mark == -1 || rna->transloc_move ==  mark) {
//                        printf("INSERT -- Update RNA from %d to %d\n", rna->begin, rna->begin + length_diff);
                        rna->begin = rna->begin+length_diff;
                    } else {

                    }
                }

                if (rna->end >= pos_after && rna->end < pos_before) {
                    if (mark == -1 || rna->transloc_move ==  mark) {
//                        printf("INSERT -- Update RNA from %d to %d\n", rna->begin, rna->begin + length_diff);

                        rna->end = rna->end+length_diff;
                    } else {

                    }
                }

            }

            // UPDATE Proteins
            for (auto prot : proteins_) {
                if (prot->protein_start >= pos_after && prot->protein_start < pos_before) {
                    if (mark == -1 ||prot->transloc_move ==  mark) {
//                        printf("INSERT -- Update Proteins from %d to %d\n", prot->protein_start,
//                               prot->protein_start + length_diff);
                        prot->protein_start = prot->protein_start+length_diff;
                    }
                }


                if (prot->protein_end >= pos_after && prot->protein_end < pos_before) {
                    if (mark == -1 ||prot->transloc_move ==  mark) {
//                        printf("INSERT -- Update Proteins from %d to %d\n", prot->protein_start,
//                               prot->protein_start + length_diff);
                        prot->protein_end = prot->protein_end+length_diff;
                    }
                }
            }
        } else if (add_or_reduce == 0){
            //REDUCE
            for (auto &strand: promoters_list_) {
                for (auto &prom: strand) {
                    if (prom.pos >= pos_after && prom.pos < pos_before) {
                        if (mark == -1 || prom.transloc_move ==  mark) {
//                            printf("DELETE -- Update Promoters from %d to %d\n", prom.pos, prom.pos - length_diff);
                            prom.pos = prom.pos-length_diff;
                        }
                    }
                }
            }

            // UPDATE Rna
            for (auto rna : rnas_) {
                if (rna->begin >= pos_after && rna->begin < pos_before) {
                    if (mark == -1 ||rna->transloc_move ==  mark) {
//                        printf("DELETE -- Update RNA from %d to %d\n", rna->begin, rna->begin - length_diff);
                        rna->begin = rna->begin-length_diff;
                    }
                }

                if (rna->end >= pos_after && rna->end < pos_before) {
                    if (mark == -1 ||rna->transloc_move ==  mark) {
                        rna->end = rna->end-length_diff;
                    }
                }
            }

            // UPDATE Proteins
            for (auto prot : proteins_) {
                if (prot->protein_start >= pos_after && prot->protein_start < pos_before) {
                    if (mark == -1 ||  prot->transloc_move ==  mark) {
//                        printf("DELETE -- Update Proteins from %d to %d\n", prot->protein_start,
//                               prot->protein_start - length_diff);
                        prot->protein_start = prot->protein_start-length_diff;
                    }
                }

                if (prot->protein_end >= pos_after && prot->protein_end < pos_before) {
                    if (mark == -1 ||  prot->transloc_move ==  mark) {
//                        printf("DELETE -- Update Proteins from %d to %d\n", prot->protein_start,
//                               prot->protein_start - length_diff);
                        prot->protein_end = prot->protein_end-length_diff;
                    }
                }
            }
        } else if (add_or_reduce == 2)  {
            // SHIFT
            for (auto &strand: promoters_list_) {
                for (auto &prom: strand) {
                    if (prom.pos >= pos_after && prom.pos < pos_before) {
                        if (mark == -1 ||prom.transloc_move == mark) {
//                            printf("SHIFT -- Update Promoters from %d to %d\n", prom.pos, prom.pos + length_diff);
                            prom.pos = prom.pos+length_diff;
                        }
                    } else {
//                        printf("NO SHIFT (INCORRECT MARK) -- Update Promoters from %d to %d __ %d != %d\n", prom.pos,
//                               prom.pos + length_diff,prom.transloc_move,mark);
                    }
                }
            }

            // UPDATE Rna
            for (auto rna : rnas_) {
                if (rna->begin >= pos_after && rna->begin < pos_before) {
                    if (mark == -1 ||  rna->transloc_move ==  mark) {
//                        printf("SHIFT -- Update RNA from %d to %d\n", rna->begin, rna->begin + length_diff);
                        rna->begin = rna->begin+length_diff;
                    }
                }

                if (rna->end >= pos_after && rna->end < pos_before) {
                    if (mark == -1 ||  rna->transloc_move ==  mark) {
//                        printf("SHIFT -- Update RNA from %d to %d\n", rna->begin, rna->begin + length_diff);
                        rna->end = rna->end+length_diff;
                    }
                }
            }

            // UPDATE Proteins
            for (auto prot : proteins_) {
                if (prot->protein_start >= pos_after && prot->protein_start < pos_before) {
                    if (mark == -1 || prot->transloc_move ==  mark ) {
//                        printf("%d -- SHIFT -- Update Proteins from %d (%d) to %d (%d)\n", AeTime::time(),
//                               prot->protein_start,
//                               prot->protein_end,
//                               prot->protein_start + length_diff,
//                               prot->protein_end + length_diff);
                        prot->protein_start = prot->protein_start+length_diff;
                    }
                }

                if (prot->protein_end >= pos_after && prot->protein_end < pos_before) {
                    if (mark == -1 || prot->transloc_move ==  mark ) {
//                        printf("%d -- SHIFT -- Update Proteins from %d (%d) to %d (%d)\n", AeTime::time(),
//                               prot->protein_start,
//                               prot->protein_end,
//                               prot->protein_start + length_diff,
//                               prot->protein_end + length_diff);
                        prot->protein_end = prot->protein_end+length_diff;
                    }
                }
            }
        }
//        printf("%d -- UPDATE_POSITIONS : End\n",indiv_->indiv_id);
    }

    void SIMD_List_Metadata::mark_positions(int pos_1, int pos_2, int8_t before) {
        for (auto &strand: promoters_list_) {
            for (auto &prom: strand) {
                if (prom.pos >= pos_1 && prom.pos < pos_2) {
                    prom.transloc_move = before;
                }
            }
        }

        // UPDATE Rna
        for (auto rna : rnas_) {
            if (rna->begin >= pos_1 && rna->begin < pos_2) {
                rna->transloc_move = before;
            }
        }

        // UPDATE Proteins
        for (auto prot : proteins_) {
            if (prot->protein_start >= pos_1 && prot->protein_start < pos_2) {
                prot->transloc_move = before;
            }
        }
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

        if (it_promoter_ == promoters_list_[LEADING].end()) {
            it_promoter_ = promoters_list_[LAGGING].begin();
        }

        return prom;
    }

    void SIMD_List_Metadata::promoter_begin() {
        it_promoter_ = promoters_list_[LEADING].begin();

        if (promoters_list_[LEADING].size() == 0) it_promoter_ = promoters_list_[LAGGING].begin();
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

    void SIMD_List_Metadata::rna_add(int idx, int32_t t_begin, int32_t t_end,
                 int8_t t_leading_lagging, double t_e,
                 int32_t t_length) {
        rnas_.push_back(new pRNA(t_begin, t_end,
                 t_leading_lagging, t_e,
                 t_length));
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

    void SIMD_List_Metadata::proteins_print() {
        int prot_idx=0;
        for (auto prot : proteins_) {
            printf("Protein %d : %lf %lf %lf %lf : ",prot_idx,prot->m,prot->w,prot->h,prot->e);
            for (auto rna : prot->rna_list_)
                printf("[%d => %d] ",rna->begin,rna->end);
            printf("\n");
            prot_idx++;
        }
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