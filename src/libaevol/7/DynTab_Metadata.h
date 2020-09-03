//
// Created by arrouan on 18/07/19.
//

#ifndef AEVOL_DYNTAB_METADATA_H
#define AEVOL_DYNTAB_METADATA_H

#include "Abstract_Metadata.h"
#include "Dna_7.h"
#include "Individual_7.h"

#include <cstdint>
#include <list>

namespace aevol {
    constexpr int32_t DYNTAB_BLOCK_SIZE = 2000;
    constexpr int32_t DYNTAB_BLOCK_MUL = 2;

    class SIMD_DynTab_Metadata : public SIMD_Abstract_Metadata {
    public:
        SIMD_DynTab_Metadata(Internal_SIMD_Struct* indiv) : SIMD_Abstract_Metadata(indiv) {
            count_promoters_ = 0;

            nb_block_dyntab_ = 2;

            dyntab_size_ = nb_block_dyntab_ * DYNTAB_BLOCK_SIZE;

            promoters_ = new promoterStruct*[dyntab_size_];
            for (int prom_idx = 0; prom_idx < dyntab_size_; prom_idx++)
                promoters_[prom_idx] = nullptr;
        };

        SIMD_DynTab_Metadata(Internal_SIMD_Struct* indiv, SIMD_DynTab_Metadata* metadata) : SIMD_Abstract_Metadata(indiv,metadata) {
            count_promoters_ = 0;

            nb_block_dyntab_ = ((metadata->promoter_count() * DYNTAB_BLOCK_MUL) / DYNTAB_BLOCK_SIZE) + 1;
            dyntab_size_ = nb_block_dyntab_ * DYNTAB_BLOCK_SIZE;

//            printf("Nb block %d (past %d) for a dyntab size of %d (past %d) -- Genome size %d -- Nb promoters %d\n",nb_block_dyntab_,metadata->nb_block_dyntab_,
//                                                dyntab_size_,metadata->dyntab_size_,indiv->dna_->length(),metadata->count_promoters_);
            promoters_ = new promoterStruct*[dyntab_size_];

            for (int prom_idx = 0; prom_idx < dyntab_size_; prom_idx++)
                promoters_[prom_idx] = nullptr;

            for (int prom_idx = 0; prom_idx < metadata->count_promoters_; prom_idx++) {
                if (metadata->promoters_[prom_idx] != nullptr) {
                    auto prom_copy = new promoterStruct(metadata->promoters_[prom_idx]->pos,
                                                        metadata->promoters_[prom_idx]->error,
                                                        metadata->promoters_[prom_idx]->leading_or_lagging);
                    promoters_[count_promoters_] = prom_copy;

                    count_promoters_++;
                }

            }
        };

        ~SIMD_DynTab_Metadata() override {
            for (int prom_idx = 0; prom_idx < dyntab_size_; prom_idx++) {
                    delete promoters_[prom_idx];
                }


            delete [] promoters_;

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

        promoterStruct* promoter_next() override ;
        void promoter_begin() override ;
        bool promoter_end() override ;

        int promoter_count() override;
        void set_promoters_count(int pcount) override;

        /*** Terminators ***/
        int terminator_count(int LoL) override;
        void terminator_add(int LoL, int dna_pos) override;

        int next_terminator(int LoL, int dna_pos) override;

        void terminators_clear() override;

        /*** RNAs ***/
        pRNA* rnas(int idx) override;
        void rna_add(int idx, pRNA* rna) override;
        void rna_add(int idx, int32_t t_begin, int32_t t_end,
                     int8_t t_leading_lagging, double t_e,
                     int32_t t_length) override;

        pRNA* rna_next() override ;
        void rna_begin() override ;
        bool rna_end() override ;

        int rna_count() override;
        void set_rna_count(int rcount) override;

        void rnas_resize(int resize) override;
        void rnas_clear() override;

        /*** Proteins ***/
        pProtein* proteins(int idx) override;
        void protein_add(int idx, pProtein* prot) override;

        pProtein* protein_next() override ;
        void protein_begin() override ;
        bool protein_end() override ;

        int proteins_count() override;
        void set_proteins_count(int pcount) override;

        void proteins_resize(int resize) override;
        void proteins_clear() override;

        /*** Promoters ***/
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


        void rebuild_index() {

        }

        void reallocate_promoters() {
            int old_dyntab_size_ = dyntab_size_;
            int old_nb_block_dyntab_ = nb_block_dyntab_;

            nb_block_dyntab_ = ((promoter_count() * DYNTAB_BLOCK_MUL) / DYNTAB_BLOCK_SIZE) + 1;
            dyntab_size_ = nb_block_dyntab_ * DYNTAB_BLOCK_SIZE;

//            printf("Nb block %d (past %d) for a dyntab size of %d (past %d) -- Genome size %d -- Nb promoters %d\n",
//                   nb_block_dyntab_,old_nb_block_dyntab_,
//                   dyntab_size_,old_dyntab_size_,indiv_->dna_->length(),count_promoters_);

            promoterStruct** old_promoters_ = promoters_;
            promoters_ = new promoterStruct*[dyntab_size_];
            for (int prom_idx = 0; prom_idx < dyntab_size_; prom_idx++)
                promoters_[prom_idx] = nullptr;

            int old_count_promoters_ = count_promoters_;
            count_promoters_ = 0;

            for (int prom_idx = 0; prom_idx < old_count_promoters_; prom_idx++) {
                if (old_promoters_[prom_idx] != nullptr) {
                    promoters_[count_promoters_] = old_promoters_[prom_idx];

                    count_promoters_++;
                }

            }

            delete [] old_promoters_;

        }

    //protected:
        promoterStruct** promoters_ = nullptr;

        std::set<int> terminator_lag_;
        std::set<int> terminator_lead_;
        std::vector<pRNA*> rnas_;
        std::vector<pProtein*> proteins_;

        int it_promoter_ = 0;
        int it_promoter_count_ = 0;
        std::vector<pRNA*>::iterator it_rna_ = rnas_.begin();
        std::vector<pProtein*>::iterator it_protein_ = proteins_.begin();

        int32_t count_promoters_ = 0;

        int32_t protein_count_ = 0;


        int32_t dyntab_size_ = 5000;
        int32_t nb_block_dyntab_ = 5;
    };
}


#endif //AEVOL_MAP_METADATA_H
