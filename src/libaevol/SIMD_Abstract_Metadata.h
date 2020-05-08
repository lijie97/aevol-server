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
#include "Dna_SIMD.h"

namespace aevol {
    class promoterStruct;

    class SIMD_Abstract_Metadata {
    public:
        SIMD_Abstract_Metadata(Internal_SIMD_Struct* indiv, SIMD_Abstract_Metadata* metadata) {
            indiv_ = indiv;
        }

        SIMD_Abstract_Metadata(Internal_SIMD_Struct* indiv) {
            indiv_ = indiv;
        }


        virtual ~SIMD_Abstract_Metadata() {};
        /** Getter **/
        /*** Promoters ***/
        inline int8_t is_promoter_leading(int pos) {
#ifdef VANILLA_SEARCH
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
#elif CMOD_SEARCH
            int8_t prom_dist_leading[26];
                int len = indiv_->dna_->length();

                for (int motif_id = 0; motif_id < 22; motif_id++) {
                        prom_dist_leading[motif_id] =
                                PROM_SEQ_LEAD[motif_id] ==
                                indiv_->dna_->data_[Utils::mod(pos + motif_id, len)]
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

#elif AMOD_SEARCH
            int8_t prom_dist_leading[26];
                int len = indiv_->dna_->length();

                if (pos+22 >= len) {
                    for (int motif_id = 0; motif_id < 22; motif_id++) {
                            prom_dist_leading[motif_id] =
                                    PROM_SEQ_LEAD[motif_id] ==
                                    indiv_->dna_->data_[Utils::mod(pos + motif_id, len)]
                                    ? 0 : 1;
                    }
                } else {
                    for (int motif_id = 0; motif_id < 22; motif_id++) {
                            prom_dist_leading[motif_id] =
                                    PROM_SEQ_LEAD[motif_id] ==
                                    indiv_->dna_->data_[pos + motif_id]
                                    ? 0 : 1;
                    }
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

#elif B1MOD_SEARCH
            int8_t prom_dist_leading[26];
                int32_t len = indiv_->dna_->length();

                for (int motif_id = 0; motif_id < 22; motif_id++) {
                        prom_dist_leading[motif_id] =
                                PROM_SEQ_LEAD[motif_id] ==
                                indiv_->dna_->data_[(pos + motif_id) - ((pos + motif_id) >= len) * len]
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
#elif B2MOD_SEARCH
            int8_t prom_dist_leading[26];
                int32_t len = indiv_->dna_->length();

                for (int motif_id = 0; motif_id < 22; motif_id++) {
                        prom_dist_leading[motif_id] =
                                PROM_SEQ_LEAD[motif_id] ==
                                indiv_->dna_->data_[(pos + motif_id) + ((((unsigned int32_t)((pos + motif_id) - len)) >> 31) -1 )* len]
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
#else // OLD_SEARCH
            int8_t prom_dist_leading = 0;
            int len = indiv_->dna_->length();

            for (int motif_id = 0; motif_id < 22; motif_id++) {
                prom_dist_leading +=
                        PROM_SEQ_LEAD[motif_id] ==
                        indiv_->dna_->data_[pos + motif_id >= len ? pos + motif_id - len : pos + motif_id]
                        ? 0 : 1;
                if (prom_dist_leading>PROM_MAX_DIFF)
                    break;
            }

            return prom_dist_leading;
#endif
        }

        inline int8_t is_promoter_lagging(int pos) {
#ifdef VANILLA_SEARCH
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
#elif CMOD_SEARCH
            int8_t prom_dist[26];
            int len = indiv_->dna_->length();

            for (int motif_id = 0; motif_id < 22; motif_id++) {
                prom_dist[motif_id] =
                        PROM_SEQ_LAG[motif_id] ==
                        indiv_->dna_->data_[Utils::mod(pos - motif_id,len)]
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
#elif AMOD_SEARCH
            int8_t prom_dist[26];
            int len = indiv_->dna_->length();

            if (pos-22 < 0) {

                for (int motif_id = 0; motif_id < 22; motif_id++) {
                    prom_dist[motif_id] =
                            PROM_SEQ_LAG[motif_id] ==
                            indiv_->dna_->data_[Utils::mod(pos - motif_id,len)]
                            ? 0 : 1;
                }
            } else {
                for (int motif_id = 0; motif_id < 22; motif_id++) {
                    prom_dist[motif_id] =
                            PROM_SEQ_LAG[motif_id] ==
                            indiv_->dna_->data_[pos - motif_id]
                            ? 0 : 1;
                }
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
#elif B1MOD
            int8_t prom_dist[26];
                int32_t len = indiv_->dna_->length();

                for (int motif_id = 0; motif_id < 22; motif_id++) {
                        prom_dist[motif_id] =
                                PROM_SEQ_LAG[motif_id] ==
                                indiv_->dna_->data_[(pos - motif_id) + ((pos - motif_id) < 0) * len]
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

#elif B2MOD
            int8_t prom_dist[26];
                int32_t len = indiv_->dna_->length();

                for (int motif_id = 0; motif_id < 22; motif_id++) {
                        prom_dist[motif_id] =
                                PROM_SEQ_LAG[motif_id] ==
                                indiv_->dna_->data_[ (pos - motif_id) + (((unsigned int32_t)((pos - motif_id))) >> 31) * len]
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
#else
            int8_t prom_dist = 0;
            int len = indiv_->dna_->length();

            for (int motif_id = 0; motif_id < 22; motif_id++) {
                prom_dist+=
                        PROM_SEQ_LAG[motif_id] ==
                        indiv_->dna_->data_[Utils::mod(pos - motif_id,len)]
                        ? 0 : 1;
                if (prom_dist>PROM_MAX_DIFF)
                    break;
            }

            return prom_dist;
#endif
            }

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
        inline int8_t is_shine_dal_start_prot_leading(int pos) {
            int32_t len = indiv_->dna_->length();
            int8_t start[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

            for (int k = 0; k < 9; k++) {
                int k_t = k >= 6 ? k + 4 : k;

                if (indiv_->dna_->data_[(pos + k_t) + ((((unsigned int32_t)((pos + k_t) - len)) >> 31) -1 )* len] ==
                    SHINE_DAL_SEQ_LEAD[k]) {
                    start[k] = true;
                } else {
                    start[k] = false;
                    break;
                }
            }

            return start[0] + start[1] + start[2] + start[3] + start[4] + start[5] + start[6]
                   + start[7] + start[8];
        }

        inline int8_t is_shine_dal_start_prot_lagging(int pos) {
            int32_t len = indiv_->dna_->length();
            int8_t start = false;

            for (int k = 0; k < 9; k++) {
                int k_t = k >= 6 ? k + 4 : k;

                if (indiv_->dna_->data_[(pos - k_t) + (((unsigned int32_t)((pos - k_t))) >> 31) * len] ==
                    SHINE_DAL_SEQ_LAG[k]) {
                    start = true;
                } else {
                    start = false;
                    break;
                }
            }
        }
/*        virtual void lst_shine_dal_start_prot(bool lorl,
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

        virtual promoterStruct* promoter_next() = 0;
        virtual void promoter_begin() = 0;
        virtual bool promoter_end() = 0;

        /*** Terminators ***/
        virtual int terminator_count(int LoL) = 0;
        virtual void terminator_add(int LoL, int dna_pos) = 0;

        virtual int next_terminator(int LoL, int dna_pos) = 0;

        virtual void terminators_clear() = 0;

        /*** RNAs ***/
        virtual pRNA* rnas(int idx) = 0;
        virtual void rna_add(int idx, pRNA* prot) = 0;
        virtual void rna_add(int idx, int32_t t_begin, int32_t t_end, int8_t t_leading_lagging, double t_e,
                             int32_t t_length) = 0;

        virtual int rna_count() = 0;
        virtual void set_rna_count(int rcount) = 0;

        virtual void rnas_resize(int resize) = 0;
        virtual void rnas_clear() = 0;

        virtual pRNA* rna_next() = 0;
        virtual void rna_begin() = 0;
        virtual bool rna_end() = 0;

        /*** Proteins ***/
        virtual pProtein* proteins(int idx) = 0;
        virtual void protein_add(int idx, pProtein* prot) = 0;

        virtual int proteins_count() = 0;
        virtual void set_proteins_count(int pcount) = 0;

        virtual void proteins_resize(int resize) = 0;
        virtual void proteins_clear() = 0;


        virtual pProtein* protein_next() = 0;
        virtual void protein_begin() = 0;
        virtual bool protein_end() = 0;
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


        int32_t length() { return indiv_->dna_->length(); };


        int32_t rna_count_ = 0;
    protected:

        Internal_SIMD_Struct* indiv_;
    };
}

#endif //AEVOL_SIMD_METADATA_H
