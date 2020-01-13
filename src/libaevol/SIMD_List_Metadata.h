//
// Created by arrouan on 26/09/19.
//

#ifndef AEVOL_SIMD_LIST_METADATA_H
#define AEVOL_SIMD_LIST_METADATA_H

#include "SIMD_Individual.h"
#include "SIMD_Abstract_Metadata.h"
#include "Dna_SIMD.h"

using std::list;


namespace aevol {

    using SIMD_Promoters1Strand = std::list<promoterStruct>;
    using SIMD_Promoters2Strands = std::vector<SIMD_Promoters1Strand>;

    class SIMD_List_Metadata  : public SIMD_Abstract_Metadata {
    public:
        SIMD_List_Metadata(Internal_SIMD_Struct* indiv) : SIMD_Abstract_Metadata(indiv) { };

        SIMD_List_Metadata(Internal_SIMD_Struct* indiv, SIMD_List_Metadata* metadata) : SIMD_Abstract_Metadata(indiv,metadata) {

            for (auto& strand: {LEADING, LAGGING}) {
                for (auto& rna: metadata->promoters_list_[strand]) {
                    promoters_list_[strand].emplace_back(rna);

                }
            }
        };

        ~SIMD_List_Metadata() override {
           promoters_list_.clear();
        };

        /** Getter **/

        /*** Promoters ***/
        promoterStruct* promoters(int idx) override;
        void promoter_add(int idx, promoterStruct* prom) override;

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


        SIMD_Promoters2Strands promoters_list_ = {{},
                                                  {}};
    protected:



        std::set<int> terminator_lag_;
        std::set<int> terminator_lead_;
        std::vector<pRNA*> rnas_;
        std::vector<pProtein*> proteins_;

        int32_t protein_count_ = 0;

        Internal_SIMD_Struct* indiv_;
    };
}


#endif //AEVOL_SIMD_LIST_METADATA_H
