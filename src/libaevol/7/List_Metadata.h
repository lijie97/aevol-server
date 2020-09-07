//
// Created by arrouan on 26/09/19.
//

#ifndef AEVOL_LIST_METADATA_H
#define AEVOL_LIST_METADATA_H


#include "Rna_7.h"
#include "Abstract_Metadata.h"
#include "Dna_7.h"
#include "Individual_7.h"

using std::list;


namespace aevol {

    using Promoters1Strand_7     = std::list<PromoterStruct>;
    using Promoters2Strands_7    = std::vector<Promoters1Strand_7>;

    class List_Metadata : public Abstract_Metadata {
    public:
     List_Metadata(Individual_7* indiv) : Abstract_Metadata(indiv) { set_iterators(); };

     List_Metadata(Individual_7* indiv, List_Metadata* metadata) : Abstract_Metadata(indiv,metadata) {

            for (auto& strand: {LEADING, LAGGING}) {
                for (auto& rna: metadata->promoters_list_[strand]) {
                    promoters_list_[strand].emplace_back(rna);

                }
            }

            set_iterators();
        };


        ~List_Metadata() override {
           promoters_list_.clear();

           for (std::list<Rna_7*>::iterator it_rna = rnas_.begin(); it_rna != rnas_.end(); it_rna++) {
               delete (*(it_rna));
           }


            for (std::list<Protein_7*>::iterator it_protein = proteins_.begin(); it_protein != proteins_.end(); it_protein++) {
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
        PromoterStruct* promoters(int idx) override;
        void promoter_add(int idx, PromoterStruct* prom) override;

        PromoterStruct* promoter_next() override ;
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
        Rna_7* rnas(int idx) override;
        void rna_add(int idx, Rna_7* rna) override;
        void rna_add(int idx, int32_t t_begin, int32_t t_end,
                     int8_t t_leading_lagging, double t_e,
                     int32_t t_length) override ;

        Rna_7* rna_next() override ;
        void rna_begin() override ;
        bool rna_end() override ;

        int rna_count() override;
        void set_rna_count(int rcount) override;

        void rnas_resize(int resize) override;
        void rnas_clear() override;

        /*** Proteins ***/
        Protein_7* proteins(int idx) override;
        void protein_add(int idx, Protein_7* prot) override;

        Protein_7* protein_next() override ;
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
                           std::list<PromoterStruct*>&  motif_list) override;

        /*** Terminator ***/
        /*int8_t is_terminator_leading(int pos);
        int8_t is_terminator_lagging(int pos);
        void lst_terminator(bool lorl,
                                    Position before_after_btw, // with regard to the strand's reading direction
                                    int32_t pos1,
                                    int32_t pos2,
                                    std::list<PromoterStruct*>& promoters_list);
*/
        /*** Shine Dal + Start Codon ***/
        /*      int8_t is_shine_dal_start_prot_leading(int pos);
              int8_t is_shine_dal_start_prot_lagging(int pos);
              void lst_shine_dal_start_prot(bool lorl,
                                                    Position before_after_btw, // with regard to the strand's reading direction
                                                    int32_t pos1,
                                                    int32_t pos2,
                                                    std::list<Abstract_Metadata*>& promoters_list);
      */
        /*** Stop Codon ***/
/*        int8_t is_stop_prot_leading(int pos);
        int8_t is_stop_prot_lagging(int pos);
        void lst_stop_prot(bool lorl,
                                   Position before_after_btw, // with regard to the strand's reading direction
                                   int32_t pos1,
                                   int32_t pos2,
                                   std::list<Abstract_Metadata*>& promoters_list);
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
                                             std::vector<std::list<PromoterStruct*>>& duplicated_promoters) override;
        void extract_promoters_included_in(int32_t pos_1,
                                           int32_t pos_2, std::vector<std::list<PromoterStruct*>>& extracted_promoters) override;
        void insert_promoters(std::vector<std::list<PromoterStruct*>>& promoters_to_insert) override;
        void insert_promoters_at(std::vector<std::list<PromoterStruct*>>& promoters_to_insert,
                                 int32_t pos) override;

        void invert_promoters_included_in(int32_t pos1,
                                          int32_t pos2) override;


        static void shift_promoters(
                std::vector<std::list<PromoterStruct*>>& promoters_to_shift,
                int32_t delta_pos,
                int32_t seq_length);
        static void invert_promoters(std::vector<std::list<PromoterStruct*>>& promoter_lists,
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
                                   std::vector<std::list<PromoterStruct*>>& promoters_list) override;

        void extract_leading_promoters_starting_between(int32_t pos_1,
                                                        int32_t pos_2, std::list<PromoterStruct*>& extracted_promoters) override;

        void extract_lagging_promoters_starting_between(int32_t pos_1,
                                                        int32_t pos_2,
                                                        std::list<PromoterStruct*>& extracted_promoters) override;

        Promoters2Strands_7 promoters_list_ = {{},
                                                  {}};

        std::list<Protein_7*> proteins_;
    protected:
     Promoters1Strand_7::iterator it_promoter_;
        int it_promoter_pos_;

        std::list<Rna_7*>::iterator it_rna_;
        std::list<Protein_7*>::iterator it_protein_;

        std::set<int> terminator_lag_;
        std::set<int> terminator_lead_;
        std::list<Rna_7*> rnas_;

        int32_t protein_count_ = 0;

        int cmp_rna = 0;

        Individual_7* indiv_;
    };
}


#endif //AEVOL_LIST_METADATA_H
