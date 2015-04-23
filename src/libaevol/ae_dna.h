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


 #ifndef __AE_DNA_H__
#define  __AE_DNA_H__


// =================================================================
//                              Libraries
// =================================================================
#include <cinttypes>
#include <cstdlib>

#include <list>
#include <vector>
#include <memory>

#include <zlib.h>

// =================================================================
//                            Project Files
// =================================================================
#include "dna_replic_report.h"
#include "ae_enums.h"
#include "ae_mutation.h"
#include "ae_string.h"
#include "ae_jumping_mt.h"


namespace aevol {


// =================================================================
//                          Class declarations
// =================================================================
class ae_exp_manager;
class ae_individual;
class GeneticUnit;
class ae_vis_a_vis;
class ae_rna;






class ae_dna : public ae_string
{
  public :
    // =================================================================
    //                             Constructors
    // =================================================================
    ae_dna(void) = delete; // TODO <david.parsons@inria.fr> merge unsure
    ae_dna(const ae_dna&) = delete; // TODO <david.parsons@inria.fr> merge unsure
    ae_dna(GeneticUnit* gen_unit,
           int32_t length,
           std::shared_ptr<ae_jumping_mt> prng);
    ae_dna(GeneticUnit* gen_unit, const ae_dna &model);
    ae_dna(GeneticUnit* gen_unit, ae_dna* const parent_dna);
    ae_dna(GeneticUnit* gen_unit, char* seq, int32_t length);
    ae_dna(GeneticUnit* gen_unit, gzFile backup_file);
    ae_dna(GeneticUnit* gen_unit, char* organism_file_name);

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ae_dna(void);

    // =================================================================
    //                              Accessors
    // =================================================================
    // From ae_string
    //   inline const char*   get_data(void) const;
    //   inline       void    set_data(char* data, int32_t length = -1);
    //   inline       int32_t get_length(void) const;
    inline DnaReplicReport*  get_replic_report(void) const;
    inline void                   set_replic_report(DnaReplicReport * rep); // for post-treatment only

    inline GeneticUnit*      get_genetic_unit(void) const;
    inline ae_individual*         get_indiv(void) const;

    char* get_subsequence(int32_t from, int32_t to, ae_strand strand) const; // WARNING : creates a new char[...] (up to you to delete it!)


    // =================================================================
    //                            Public Methods
    // =================================================================
    // Perform all the mutations (local mutations, rearrangements and transfer)
    void perform_mutations(int32_t parent_id);

    // Perform all the local mutations (point mutations and indels) of the replication
    void do_small_mutations(void);

    // Perform all the chromosomic rearrangements (duplications, deletions, translocations and inversions)
    // of the replication
    void do_rearrangements(void);
    void do_rearrangements_with_align(void);

    // Perform all transfer (with insertion and with replacement)
    void do_transfer(int32_t parent_id);

    // Perform a single local mutation at a random position
    ae_mutation* do_switch(void);
    ae_mutation* do_small_insertion(void);
    ae_mutation* do_small_deletion(void);

    // Perform a single local mutation at a specified position (useful to replay the evolution)
    bool do_switch(int32_t pos);
    bool do_small_insertion(int32_t pos, int16_t nb_insert, char * seq);
    bool do_small_deletion(int32_t pos, int16_t nb_del);

    // Perform a single rearrangement at random positions
    ae_mutation* do_duplication(void);
    ae_mutation* do_deletion(void);
    ae_mutation* do_translocation(void);
    ae_mutation* do_inter_GU_translocation(void);
    ae_mutation* do_inversion(void);
    ae_mutation* do_insertion(const char* seq_to_insert, int32_t seq_length = -1);

    // Perform a single rearrangement at specified positions
    bool do_duplication(int32_t pos_1, int32_t pos_2, int32_t pos_3);
    bool do_deletion(int32_t pos_1, int32_t pos_2);
    bool do_translocation(int32_t pos_1, int32_t pos_2, int32_t pos_3, int32_t pos_4, bool invert);
    bool do_inter_GU_translocation(int32_t pos_1, int32_t pos_2, int32_t pos_3, int32_t pos_4, bool invert);
    bool do_inversion(int32_t pos_1, int32_t pos_2);
    bool do_insertion(int32_t pos, const char* seq_to_insert, int32_t seq_length);

    // Perform transfer with the search of alignments
    ae_mutation* do_ins_HT(int32_t parent_id);
    ae_mutation* do_repl_HT(int32_t parent_id);

    // Perform a single transfer at specified positions
    bool do_ins_HT(int32_t pos, const char* seq_to_insert, int32_t seq_length);
    bool do_repl_HT(int32_t pos1, int32_t pos2, const char* seq_to_insert, int32_t seq_length);

    GeneticUnit*  extract_into_new_GU(int32_t pos_1, int32_t pos_2);
    GeneticUnit*  copy_into_new_GU(int32_t pos_1, int32_t pos_2) const;
    void insert_GU(GeneticUnit* GU_to_insert, int32_t pos_B, int32_t pos_D, bool invert);

    ae_vis_a_vis* search_alignment(ae_dna* chrom2, int32_t& nb_pairs, ae_sense sense);
    ae_vis_a_vis* search_alignment_around_positions(ae_dna* chrom2, int32_t chrom1_pos_1, int32_t chrom1_pos_2, ae_sense sense, int8_t& research_sense);

    void undergo_this_mutation(const ae_mutation * mut); // useful when we replay the evolution

    void compute_statistical_data(void);

    // TODO rna_list should be passed by ref
    static void set_GU(std::vector<std::list<ae_rna*>> rna_list, const GeneticUnit* GU);


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
    ae_exp_manager* _exp_m;

    // From ae_string
    //   char*   _data;
    //   int32_t _length;
    //   int32_t _nb_blocks;

    ae_individual* _indiv;
    GeneticUnit* _gen_unit; // Genetic unit which the dna sequence belongs to
    DnaReplicReport* _replic_report;
};


// =====================================================================
//                          Accessors definitions
// =====================================================================
inline DnaReplicReport* ae_dna::get_replic_report(void) const
{
  return _replic_report;
}

inline ae_individual* ae_dna::get_indiv(void) const
{
  return _indiv;
}

 // for post-treatment only
inline void ae_dna::set_replic_report(DnaReplicReport * rep)
{
  _replic_report = rep;
}

inline GeneticUnit * ae_dna::get_genetic_unit(void) const
{
  return _gen_unit;
}


// =====================================================================
//                       Inline functions' definition
// =====================================================================

} // namespace aevol
#endif // __AE_DNA_H__
