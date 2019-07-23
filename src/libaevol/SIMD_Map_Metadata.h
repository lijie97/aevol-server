//
// Created by arrouan on 18/07/19.
//

#ifndef AEVOL_SIMD_MAP_METADATA_H
#define AEVOL_SIMD_MAP_METADATA_H

#include "SIMD_Individual.h"
#include "SIMD_Abstract_Metadata.h"
#include "Dna_SIMD.h"


namespace aevol {
    class SIMD_Map_Metadata : public SIMD_Abstract_Metadata {
    public:
        SIMD_Map_Metadata(Internal_SIMD_Struct* indiv) { indiv_ = indiv; count_promoters_ = 0; };

        SIMD_Map_Metadata(Internal_SIMD_Struct* indiv, SIMD_Map_Metadata* metadata) {
            indiv_ = indiv;
            count_promoters_ = 0;

            for (const auto& prom : metadata->promoters_) {
                if (prom.second != nullptr) {
                    auto prom_copy = new promoterStruct(prom.second->pos, prom.second->error,
                                                        prom.second->leading_or_lagging);
                    promoters_[count_promoters_] = prom_copy;


                    if (prom.second->leading_or_lagging) {
                        leading_prom_pos_[prom_copy->pos] = count_promoters_;
                    } else {
                        lagging_prom_pos_[prom_copy->pos] = count_promoters_;
                    }

                    count_promoters_++;
                }
            }
        };

        ~SIMD_Map_Metadata() override {
            if (! promoters_.empty()) {
                for (auto element = promoters_.begin();
                     element != promoters_.end(); ++element) {
                    if (element->second != nullptr) delete element->second;
                }
            }

            promoters_.clear();

            leading_prom_pos_.clear();
            lagging_prom_pos_.clear();


            for (auto rn : rnas_) {
                delete rn;
            }

            for (auto prot : proteins_) {
                delete prot;
            }

            rnas_.clear();
            proteins_.clear();
        };

        /** Getter **/

        /*** Promoters ***/
        promoterStruct* promoters(int idx) override;
        void promoter_add(int idx, promoterStruct* prom) override;

        int promoter_count() override;
        void set_promoters_count(int pcount) override;

        /*** Terminators ***/
        int terminator_count(int LoL) override;
        int terminator_add(int LoL, int dna_pos) override;

        int next_terminator(int LoL, int dna_pos) override;

        void terminators_clear() override;

        /*** RNAs ***/
        pRNA* rnas(int idx) override;
        void rna_add(int idx, pRNA* rna) override;

        int rna_count() override;
        void set_rna_count(int rcount) override;

        void rnas_resize(int resize) override;
        void rnas_clear() override;

        /*** Proteins ***/
        pProtein* proteins(int idx) override;
        void protein_add(int idx, pProtein* prot) override;

        int proteins_count() override;
        void set_proteins_count(int pcount) override;

        void proteins_resize(int resize) override;
        void proteins_clear() override;

        /*** Promoters ***/
        int8_t is_promoter_leading(int pos) override;
        int8_t is_promoter_lagging(int pos) override;
        void lst_promoters(bool lorl,
                                   Position before_after_btw, // with regard to the strand's reading direction
                                   int32_t pos1,
                                   int32_t pos2,
                                   std::list<promoterStruct*>&  motif_list) override;

        /*** Terminator ***/
        /*int8_t is_terminator_leading(int pos);
        int8_t is_terminator_lagging(int pos);
        void lst_terminator(bool lorl,
                                    Position before_after_btw, // with regard to the strand's reading direction
                                    int32_t pos1,
                                    int32_t pos2,
                                    std::list<promoterStruct*>& promoters_list);
*/
        /*** Shine Dal + Start Codon ***/
  /*      int8_t is_shine_dal_start_prot_leading(int pos);
        int8_t is_shine_dal_start_prot_lagging(int pos);
        void lst_shine_dal_start_prot(bool lorl,
                                              Position before_after_btw, // with regard to the strand's reading direction
                                              int32_t pos1,
                                              int32_t pos2,
                                              std::list<SIMD_Abstract_Metadata*>& promoters_list);
*/
        /*** Stop Codon ***/
/*        int8_t is_stop_prot_leading(int pos);
        int8_t is_stop_prot_lagging(int pos);
        void lst_stop_prot(bool lorl,
                                   Position before_after_btw, // with regard to the strand's reading direction
                                   int32_t pos1,
                                   int32_t pos2,
                                   std::list<SIMD_Abstract_Metadata*>& promoters_list);
*/
        /** Search and update **/
        void remove_promoters_around(int32_t pos_1) override;
        void remove_promoters_around(int32_t pos_1, int32_t pos_2) override;
        void remove_all_promoters() override;

        void look_for_new_promoters_around(int32_t pos_1, int32_t pos_2) override;
        void look_for_new_promoters_around(int32_t pos) override;

        void locate_promoters() override;

        void move_all_promoters_after(int32_t pos, int32_t delta_pos) override;

        void duplicate_promoters_included_in(int32_t pos_1,
                                                     int32_t pos_2,
                                                     std::vector<std::list<promoterStruct*>>& duplicated_promoters) override;
        void extract_promoters_included_in(int32_t pos_1,
                                                   int32_t pos_2, std::vector<std::list<promoterStruct*>>& extracted_promoters) override;
        void insert_promoters(std::vector<std::list<promoterStruct*>>& promoters_to_insert) override;
        void insert_promoters_at(std::vector<std::list<promoterStruct*>>& promoters_to_insert,
                                         int32_t pos) override;

        void invert_promoters_included_in(int32_t pos1,
                                                  int32_t pos2) override;


        static void shift_promoters(
                std::vector<std::list<promoterStruct*>>& promoters_to_shift,
                int32_t delta_pos,
                int32_t seq_length);
        static void invert_promoters(std::vector<std::list<promoterStruct*>>& promoter_lists,
                                             int32_t pos1,
                                             int32_t pos2);

        void remove_leading_promoters_starting_between(int32_t pos_1,
                                                               int32_t pos_2) override;
        void remove_leading_promoters_starting_after(int32_t pos) override;
        void remove_leading_promoters_starting_before(int32_t pos) override;

        void remove_lagging_promoters_starting_between(int32_t pos_1,
                                                               int32_t pos_2) override;
        void remove_lagging_promoters_starting_after(int32_t pos) override;
        void remove_lagging_promoters_starting_before(int32_t pos) override;

        void move_all_leading_promoters_after(int32_t pos, int32_t delta_pos) override;
        void move_all_lagging_promoters_after(int32_t pos,int32_t delta_pos) override;

        void look_for_new_leading_promoters_starting_between(int32_t pos_1, int32_t pos_2) override;
        void look_for_new_leading_promoters_starting_after(int32_t pos) override;
        void look_for_new_leading_promoters_starting_before(int32_t pos) override;

        void look_for_new_lagging_promoters_starting_between(int32_t pos_1,int32_t pos_2) override;
        void look_for_new_lagging_promoters_starting_after(int32_t pos) override;
        void look_for_new_lagging_promoters_starting_before(int32_t pos) override;

        void promoters_included_in(int32_t pos_1,
                                           int32_t pos_2,
                                           std::vector<std::list<promoterStruct*>>& promoters_list) override;

        void extract_leading_promoters_starting_between(int32_t pos_1,
                                                                int32_t pos_2, std::list<promoterStruct*>& extracted_promoters) override;

        void extract_lagging_promoters_starting_between(int32_t pos_1,
                                                                int32_t pos_2,
                                                                std::list<promoterStruct*>& extracted_promoters) override;


        int32_t length() override { indiv_->dna_->length(); };

        void rebuild_index() {
            if (count_promoters_ > (int)promoters_.size()/2) {
                /**
                 * Do the reindexation process
                 */
                auto old_promoters = promoters_;
                promoters_.clear();
                leading_prom_pos_.clear();
                lagging_prom_pos_.clear();
                count_promoters_ = 0;
                for (auto prom : old_promoters) {
                    promoters_[count_promoters_] = prom.second;
                    if (prom.second->leading_or_lagging) {
                        leading_prom_pos_[prom.second->pos] = count_promoters_;
                    } else {
                        lagging_prom_pos_[prom.second->pos] = count_promoters_;
                    }
                    count_promoters_++;
                }
            }
        }

    protected:
        std::map<int32_t,promoterStruct*> promoters_;
        std::map<int32_t,int32_t> leading_prom_pos_;
        std::map<int32_t,int32_t> lagging_prom_pos_;

        std::set<int> terminator_lag_;
        std::set<int> terminator_lead_;
        std::vector<pRNA*> rnas_;
        std::vector<pProtein*> proteins_;

        int32_t count_promoters_;
        int32_t protein_count_ = 0;
        int32_t rna_count_ = 0;

        Internal_SIMD_Struct* indiv_;
    };
}


#endif //AEVOL_SIMD_MAP_METADATA_H
