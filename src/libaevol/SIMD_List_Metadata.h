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
        SIMD_List_Metadata(Internal_SIMD_Struct* indiv) { indiv_ = indiv; set_iterators(); };

        SIMD_List_Metadata(Internal_SIMD_Struct* indiv, SIMD_List_Metadata* metadata)  {
            indiv_ = indiv;

            // Copy promoters
            for (auto& strand: {LEADING, LAGGING}) {
                for (auto& rna: metadata->promoters_list_[strand]) {
                    promoters_list_[strand].emplace_back(rna);
                }
            }

            // Copy RNAs
            for (auto rna : metadata->rnas_) {
                pRNA* nrna = new pRNA(rna);
                rnas_.push_back(nrna);
                nrna->to_recompute = false;
            }

            // Copy Proteins
            for (auto prot : metadata->proteins_) {
                pProtein* nprot = new pProtein(prot);
                proteins_.push_back(nprot);
                for (auto rna : prot->rna_list_) {
                    for (auto rna2 : rnas_) {
                        if (rna->begin == rna2->begin) {
//                            printf("Add RNA %d to Protein %d\n",rna2->begin,nprot->protein_start);
                            proteins_.back()->rna_list_.push_back(rna2);
                            proteins_.back()->to_retranslate = false;
                            proteins_.back()->to_recompute = false;
                        }
                    }
                }
            }

            set_iterators();

//            printf("End of constructor List Metadata\n");
//            display(false);
        };


        ~SIMD_List_Metadata() override {
           promoters_list_.clear();

           for (std::list<pRNA*>::iterator it_rna = rnas_.begin(); it_rna != rnas_.end(); it_rna++) {
               delete (*(it_rna));
           }


            for (std::list<pProtein*>::iterator it_protein = proteins_.begin(); it_protein != proteins_.end(); it_protein++) {
                delete (*(it_protein));
            }
        };

        /** Getter **/

        void set_iterators() {
            it_promoter_ = promoters_list_[LEADING].begin();
            it_rna_ = rnas_.begin();
            it_protein_ = proteins_.begin();
        };

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
                     int32_t t_length) override ;

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

        void proteins_print();

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

        void update_range(int32_t begin, int32_t end) override;
        void update_range(int32_t pos) override;
        void remove_range(int32_t begin, int32_t end) override;
        void remove_range(int32_t pos) override;

        void cleanup() override;
        void compute_promoters() override;
        void search_for_rna_terminator() ;
        void recompute_rna() override;
        void recompute_proteins() override;
        void retranslate_proteins() override;

        void update_metadata() override;
        void update_metadata_before_new_generation() override;

        void update_positions(int pos_after, bool insert_or_remove, int length_diff) override;
        void update_positions(int pos_after, int pos_before, int8_t add_or_reduce, int length_diff, int8_t mark = -1) override;
        void mark_positions(int pos_1, int pos_2, int8_t before) override;

        void display(bool check_to_delete) override ;

        SIMD_Promoters2Strands promoters_list_ = {{},
                                                  {}};

        std::list<pProtein*> proteins_;
    protected:
        SIMD_Promoters1Strand::iterator it_promoter_;
        int it_promoter_pos_;

        std::list<pRNA*>::iterator it_rna_;
        std::list<pProtein*>::iterator it_protein_;

        std::set<int> terminator_lag_;
        std::set<int> terminator_lead_;
        std::list<pRNA*> rnas_;

        int32_t protein_count_ = 0;

        int cmp_rna = 0;

    };
}


#endif //AEVOL_SIMD_LIST_METADATA_H
