//
// Created by arrouan on 18/07/19.
//

#ifndef AEVOL_SIMD_METADATA_H
#define AEVOL_SIMD_METADATA_H

#include <cstdint>
#include <vector>
#include <list>
#include "ae_enums.h"
#include "SIMD_Individual.h"

namespace aevol {
    class promoterStruct;

    class SIMD_Abstract_Metadata {
    public:
        //SIMD_Abstract_Metadata(SIMD_Abstract_Metadata* metadata);
        //SIMD_Abstract_Metadata(Internal_SIMD_Struct* indiv);
        virtual ~SIMD_Abstract_Metadata() {};
        /** Getter **/
        /*** Promoters ***/
        virtual int8_t is_promoter_leading(int pos) = 0;
        virtual int8_t is_promoter_lagging(int pos) = 0;
        virtual void lst_promoters(bool lorl,
                           Position before_after_btw, // with regard to the strand's reading direction
                           int32_t pos1,
                           int32_t pos2, std::list<promoterStruct*>& motif_list) = 0;

        /*** Terminator ***/
        /*virtual int8_t is_terminator_leading(int pos);
        virtual int8_t is_terminator_lagging(int pos);
        virtual void lst_terminator(bool lorl,
                           Position before_after_btw, // with regard to the strand's reading direction
                           int32_t pos1,
                           int32_t pos2,
                           std::list<promoterStruct*>& promoters_list);
*/
        /*** Shine Dal + Start Codon ***/
/*        virtual int8_t is_shine_dal_start_prot_leading(int pos);
        virtual int8_t is_shine_dal_start_prot_lagging(int pos);
        virtual void lst_shine_dal_start_prot(bool lorl,
                           Position before_after_btw, // with regard to the strand's reading direction
                           int32_t pos1,
                           int32_t pos2,
                           std::list<promoterStruct*>& promoters_list);
*/
        /*** Stop Codon ***/
/*        virtual int8_t is_stop_prot_leading(int pos);
        virtual int8_t is_stop_prot_lagging(int pos);
        virtual void lst_stop_prot(bool lorl,
                                      Position before_after_btw, // with regard to the strand's reading direction
                                      int32_t pos1,
                                      int32_t pos2,
                                      std::list<promoterStruct*>& promoters_list);


*/
        /*** Promoters ***/
        virtual promoterStruct* promoters(int idx) = 0;
        virtual void promoter_add(int idx, promoterStruct* prom) = 0;

        virtual int promoter_count() = 0;
        virtual void set_promoters_count(int rcount) = 0;

        /*** Terminators ***/
        virtual int terminator_count(int LoL) = 0;
        virtual int terminator_add(int LoL, int dna_pos) = 0;

        virtual int next_terminator(int LoL, int dna_pos) = 0;

        virtual void terminators_clear() = 0;

        /*** RNAs ***/
        virtual pRNA* rnas(int idx) = 0;
        virtual void rna_add(int idx, pRNA* prot) = 0;

        virtual int rna_count() = 0;
        virtual void set_rna_count(int rcount) = 0;

        virtual void rnas_resize(int resize) = 0;
        virtual void rnas_clear() = 0;

        /*** Proteins ***/
        virtual pProtein* proteins(int idx) = 0;
        virtual void protein_add(int idx, pProtein* prot) = 0;

        virtual int proteins_count() = 0;
        virtual void set_proteins_count(int pcount) = 0;

        virtual void proteins_resize(int resize) = 0;
        virtual void proteins_clear() = 0;

        /** Search and update **/
        virtual void remove_promoters_around(int32_t pos_1) = 0;
        virtual void remove_promoters_around(int32_t pos_1, int32_t pos_2) = 0;
        virtual void remove_all_promoters() = 0;

        virtual void look_for_new_promoters_around(int32_t pos_1, int32_t pos_2) = 0;
        virtual void look_for_new_promoters_around(int32_t pos) = 0;

        virtual void locate_promoters() = 0;

        virtual void move_all_promoters_after(int32_t pos, int32_t delta_pos) = 0;

        virtual void duplicate_promoters_included_in(int32_t pos_1,
                                             int32_t pos_2,
                                             std::vector<std::list<promoterStruct*>>& duplicated_promoters) = 0;
        virtual void extract_promoters_included_in(int32_t pos_1,
                                           int32_t pos_2, std::vector<std::list<promoterStruct*>>& extracted_promoters) = 0;
        virtual void insert_promoters(std::vector<std::list<promoterStruct*>>& promoters_to_insert) = 0;
        virtual void insert_promoters_at(std::vector<std::list<promoterStruct*>>& promoters_to_insert,
                                 int32_t pos) = 0;

        virtual void invert_promoters_included_in(int32_t pos1,
                                          int32_t pos2) = 0;


        static void shift_promoters(
                std::vector<std::list<promoterStruct*>>& promoters_to_shift,
                int32_t delta_pos,
                int32_t seq_length);
        static void invert_promoters(std::vector<std::list<promoterStruct*>>& promoter_lists,
                                     int32_t pos1,
                                     int32_t pos2);

        virtual void remove_leading_promoters_starting_between(int32_t pos_1,
                                                       int32_t pos_2) = 0;
        virtual void remove_leading_promoters_starting_after(int32_t pos) = 0;
        virtual void remove_leading_promoters_starting_before(int32_t pos) = 0;

        virtual void remove_lagging_promoters_starting_between(int32_t pos_1,
                                                       int32_t pos_2) = 0;
        virtual void remove_lagging_promoters_starting_after(int32_t pos) = 0;
        virtual void remove_lagging_promoters_starting_before(int32_t pos) = 0;

        virtual void move_all_leading_promoters_after(int32_t pos, int32_t delta_pos) = 0;
        virtual void move_all_lagging_promoters_after(int32_t pos,int32_t delta_pos) = 0;

        virtual void look_for_new_leading_promoters_starting_between(int32_t pos_1, int32_t pos_2) = 0;
        virtual void look_for_new_leading_promoters_starting_after(int32_t pos) = 0;
        virtual void look_for_new_leading_promoters_starting_before(int32_t pos) = 0;

        virtual void look_for_new_lagging_promoters_starting_between(int32_t pos_1,int32_t pos_2) = 0;
        virtual void look_for_new_lagging_promoters_starting_after(int32_t pos) = 0;
        virtual void look_for_new_lagging_promoters_starting_before(int32_t pos) = 0;

        virtual void promoters_included_in(int32_t pos_1,
                                   int32_t pos_2,
                                   std::vector<std::list<promoterStruct*>>& promoters_list) = 0;

        virtual void extract_leading_promoters_starting_between(int32_t pos_1,
                                                        int32_t pos_2, std::list<promoterStruct*>& extracted_promoters) = 0;

        virtual void extract_lagging_promoters_starting_between(int32_t pos_1,
                                                        int32_t pos_2,
                                                        std::list<promoterStruct*>& extracted_promoters) = 0;


        virtual int32_t length() = 0;
    };
}

#endif //AEVOL_SIMD_METADATA_H
