// ****************************************************************************
//
//          Aevol - An in silico experimental evolution platform
//
// ****************************************************************************
//
// Copyright: See the AUTHORS file provided with the package or <www.aevol.fr>
// Web: http://www.aevol.fr/
// E-mail: See <http://www.aevol.fr/contact/>
// Original Authors : Guillaume Beslon, Carole Knibbe, David Parsons
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ****************************************************************************


#ifndef AEVOL_DNA_H__
#define AEVOL_DNA_H__


// =================================================================
//                              Includes
// =================================================================
#include <cinttypes>
#include <cstdlib>

#include <list>
#include <vector>
#include <memory>

#include <zlib.h>

#include "ae_enums.h"
#include "Mutation.h"
#include "ae_string.h"
#include "JumpingMT.h"


namespace aevol {


// =================================================================
//                          Class declarations
// =================================================================
class ExpManager;
class Individual;
class GeneticUnit;
class VisAVis;
class Rna;
class Rna_R;






class Dna : public ae_string
{
  public :
    // =================================================================
    //                             Constructors
    // =================================================================
    Dna(void) = delete; // TODO <david.parsons@inria.fr> merge unsure
    Dna(const Dna &) = delete; // TODO <david.parsons@inria.fr> merge unsure
    Dna(GeneticUnit* gen_unit,
           int32_t length,
           std::shared_ptr<JumpingMT> prng);
    Dna(GeneticUnit* gen_unit, const Dna &model);
    Dna(GeneticUnit* gen_unit, Dna * const parent_dna);
    Dna(GeneticUnit* gen_unit, char* seq, int32_t length);
    Dna(GeneticUnit* gen_unit, gzFile backup_file);
    Dna(GeneticUnit* gen_unit, char* organism_file_name);

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~Dna(void);

    // =================================================================
    //                              Accessors
    // =================================================================
    // From ae_string
    //   inline const char*   get_data(void) const;
    //   inline       void    set_data(char* data, int32_t length = -1);
    //   inline       int32_t get_length(void) const;

    inline GeneticUnit*   get_genetic_unit(void) const;
    inline Individual * get_indiv(void) const;

    char* get_subsequence(int32_t from, int32_t to, Strand strand) const; // WARNING : creates a new char[...] (up to you to delete it!)


    // =================================================================
    //                            Public Methods
    // =================================================================
    // Perform all the mutations (local mutations, rearrangements and transfer)
    void perform_mutations(int32_t parent_id);

    // Perform all the local mutations (point mutations and indels)
    void do_small_mutations(void);

    // Perform all the chromosomal rearrangements (duplications, deletions,
    // translocations and inversions)
    void do_rearrangements(void);
    void do_rearrangements_with_align(void);

    // Perform all transfer (with insertion and with replacement)
    void do_transfer(int32_t parent_id);

    // Perform a single local mutation at a random position
    Mutation* do_switch(void);
    Mutation* do_small_insertion(void);
    Mutation* do_small_deletion(void);

    // Perform a single local mutation at a specified position (useful to replay the evolution)
    bool do_switch(int32_t pos);
    bool do_small_insertion(int32_t pos, int16_t nb_insert, char * seq);
    bool do_small_deletion(int32_t pos, int16_t nb_del);

    // Perform a single rearrangement at random positions
    Mutation* do_duplication(void);
    Mutation* do_deletion(void);
    Mutation* do_translocation(void);
    Mutation* do_inter_GU_translocation(void);
    Mutation* do_inversion(void);
    Mutation* do_insertion(const char* seq_to_insert, int32_t seq_length = -1);

    // Perform a single rearrangement at specified positions
    bool do_duplication(int32_t pos_1, int32_t pos_2, int32_t pos_3);
    bool do_deletion(int32_t pos_1, int32_t pos_2);
    bool do_translocation(int32_t pos_1, int32_t pos_2, int32_t pos_3, int32_t pos_4, bool invert);
    bool do_inter_GU_translocation(int32_t pos_1, int32_t pos_2, int32_t pos_3, int32_t pos_4, bool invert);
    bool do_inversion(int32_t pos_1, int32_t pos_2);
    bool do_insertion(int32_t pos, const char* seq_to_insert, int32_t seq_length);

    // Perform transfer with the search of alignments
    Mutation* do_ins_HT(int32_t parent_id);
    Mutation* do_repl_HT(int32_t parent_id);

    // Perform a single transfer at specified positions
    bool do_ins_HT(int32_t pos, const char* seq_to_insert, int32_t seq_length);
    bool do_repl_HT(int32_t pos1, int32_t pos2, const char* seq_to_insert, int32_t seq_length);

    GeneticUnit*  extract_into_new_GU(int32_t pos_1, int32_t pos_2);
    GeneticUnit*  copy_into_new_GU(int32_t pos_1, int32_t pos_2) const;
    void insert_GU(GeneticUnit* GU_to_insert, int32_t pos_B, int32_t pos_D, bool invert);

    VisAVis* search_alignment(Dna * chrom2, int32_t& nb_pairs, AlignmentSense sense);
    VisAVis* search_alignment_around_positions(Dna * chrom2, int32_t chrom1_pos_1, int32_t chrom1_pos_2, AlignmentSense sense, int8_t& research_sense);

    // Undergo a specific mutation
    // useful when we replay the evolution
    void undergo_this_mutation(const Mutation& mut);

    void compute_statistical_data(void);

    // TODO vld rna_list should be passed by reference
    #ifndef __REGUL
    static void set_GU(std::vector<std::list<Rna>> rna_list, const GeneticUnit* GU);
    #else
    static void set_GU(std::vector<std::list<Rna_R>> rna_list, const GeneticUnit* GU);
    #endif

    // =================================================================
    //                           Public Attributes
    // =================================================================





  protected :
    // =================================================================
    //                           Protected Methods
    // =================================================================
    void ABCDE_to_ADCBE(  int32_t pos_B, int32_t pos_C, int32_t pos_D, int32_t pos_E);
    void ABCDE_to_ADBpCpE(int32_t pos_B, int32_t pos_C, int32_t pos_D, int32_t pos_E);
    void ABCDE_to_ACpDpBE(int32_t pos_B, int32_t pos_C, int32_t pos_D, int32_t pos_E);
    void inter_GU_ABCDE_to_ACDBE(int32_t pos_B, int32_t pos_C, int32_t pos_E);
    void inter_GU_ABCDE_to_BDCAE(int32_t pos_B, int32_t pos_C, int32_t pos_E);


    // =================================================================
    //                          Protected Attributes
    // =================================================================
    ExpManager* _exp_m;

    // From ae_string
    //   char*   _data;
    //   int32_t _length;
    //   int32_t _nb_blocks;

    Individual* _indiv;
    GeneticUnit* _gen_unit; // Genetic unit which the dna sequence belongs to
};


// =====================================================================
//                          Accessors definitions
// =====================================================================
inline Individual *Dna::get_indiv(void) const
{
  return _indiv;
}

inline GeneticUnit*Dna::get_genetic_unit(void) const
{
  return _gen_unit;
}


// =====================================================================
//                       Inline functions' definition
// =====================================================================

} // namespace aevol
#endif // AEVOL_DNA_H__
