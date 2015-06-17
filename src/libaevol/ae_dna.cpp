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
//*****************************************************************************




// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>
#include <stdio.h>
#include <math.h>

#include <list>
#include <vector>

// =================================================================
//                            Project Files
// =================================================================
#include "ExpManager.h"
#include "ExpSetup.h"
#include "ae_dna.h"
#include "genetic_unit.h"
#include "Individual.h"
#include "ae_rna.h"
#include "ae_utils.h"
#include "ae_vis_a_vis.h"
#include "ae_align.h"

namespace aevol {

//##############################################################################
//                                                                             #
//                                Class ae_dna                                 #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
/**
 * Create a random dna sequence of length <length> belonging to <gen_unit>.
 */
ae_dna::ae_dna(GeneticUnit* gen_unit,
               int32_t length,
               std::shared_ptr<ae_jumping_mt> prng) :
    ae_string(length, prng)
{
  _gen_unit = gen_unit;
  _exp_m    = gen_unit->get_exp_m();
  _indiv    = gen_unit->get_indiv();

  _replic_report  = NULL;
}

/**
 * Create a new piece of dna identical to the model but belonging to <gen_unit>
 * The replication report is copied if it exists
 */
ae_dna::ae_dna(GeneticUnit* gen_unit, const ae_dna &model) :
    ae_string(model)
{
  _gen_unit = gen_unit;
  _exp_m    = gen_unit->get_exp_m();
  _indiv    = gen_unit->get_indiv();

  _replic_report = model._replic_report ?
                   new DnaReplicReport(*(model._replic_report)) :
                   NULL;
}

/**
 * Creates a new piece of dna identical to the parent's but belonging to <gen_unit>
 * The replication report is set to NULL
 */
ae_dna::ae_dna(GeneticUnit* gen_unit, ae_dna* const parent_dna) :
ae_string(parent_dna->_data, parent_dna->_length)
{
_gen_unit = gen_unit;
_exp_m    = gen_unit->get_exp_m();
_indiv    = gen_unit->get_indiv();

_replic_report = NULL;
}

/**
 * Creates a new piece of dna with sequence <seq> (of length <length>).
 * WARNING : <seq> will be used directly as the new dna sequence (it will not be copied),
 *           which means the caller must not delete it.
 * The replication report is set to NULL
 */
ae_dna::ae_dna(GeneticUnit* gen_unit, char* seq, int32_t length) :
    ae_string(seq, length, true)
{
  _gen_unit = gen_unit;
  _exp_m    = gen_unit->get_exp_m();
  _indiv    = gen_unit->get_indiv();

  _replic_report = NULL;
}

/**
 * Loads a piece of dna from <backup_file>
 * The replication report is set to NULL
 */
ae_dna::ae_dna(GeneticUnit* gen_unit, gzFile backup_file) : ae_string(backup_file)
{
  _gen_unit = gen_unit;
  _exp_m    = gen_unit->get_exp_m();
  _indiv    = gen_unit->get_indiv();

  _replic_report = NULL;
}

/**
 * Creates a dna sequence from a text file
 * The replication report is set to NULL
 */
ae_dna::ae_dna(GeneticUnit* gen_unit, char* organism_file_name) : ae_string(organism_file_name)
{
  _gen_unit = gen_unit;
  _exp_m    = gen_unit->get_exp_m();
  _indiv    = gen_unit->get_indiv();

  _replic_report  = NULL;
}

// =================================================================
//                             Destructors
// =================================================================
ae_dna::~ae_dna(void)
{
}

// =================================================================
//                         Non inline Accessors
// =================================================================
char* ae_dna::get_subsequence(int32_t from, int32_t to, ae_strand strand) const
{
  char* subseq = NULL;

  from  = ae_utils::mod(from, _length);
  to    = ae_utils::mod(to, _length);

  if (strand == LEADING)
  {
    if (from < to)
    {
      subseq = new char[to-from+1];
      subseq[to-from] = '\0';
      strncpy(subseq, &(_data[from]), to-from);
    }
    else
    {
      subseq = new char[_length-from+to+1];
      subseq[_length-from+to] = '\0';
      strncpy(subseq, &(_data[from]), _length-from);
      strncpy(&subseq[_length-from], _data, to);
    }
  }
  else // if (strand == LAGGING)
  {
    if (from > to)
    {
      subseq = new char[from-to+1];
      subseq[from-to] = '\0';

      for (int32_t i = 0 ; i < from - to ; i++)
      {
        subseq[i] = (_data[from-1-i] == '1') ? '0' : '1';
      }
    }
    else
    {
      subseq = new char[from+_length-to+1];
      subseq[from+_length-to] = '\0';

      for (int32_t i = 0 ; i < from ; i++)
      {
        subseq[i] = (_data[from-1-i] == '1') ? '0' : '1';
      }
      for (int32_t i = 0 ; i < _length-to ; i++)
      {
        subseq[from+i] = (_data[_length-1-i] == '1') ? '0' : '1';
      }
    }
  }

  return subseq;
}

// =================================================================
//                            Public Methods
// =================================================================
void ae_dna::perform_mutations(int32_t parent_id)
{
  if (_exp_m->get_output_m()->get_record_tree() && _exp_m->get_output_m()->get_tree_mode() == NORMAL)
  {
    _replic_report = new DnaReplicReport();
  }

  if(_indiv->get_with_HT())
  {
    do_transfer(parent_id);
  }

  if (_indiv->get_with_alignments())
  {
    do_rearrangements_with_align();
  }
  else
  {
    do_rearrangements();
  }

  do_small_mutations();

  // Store mutation events in replication report (temporary lists are emptied)
  // TODO : next 2 lines
  //~ new_indiv->_replic_report->set_parent_genome_size((this->_genome)->get_length());
  // if crossover => please set also donor_index, donor_misadaptation_value and donor_genome_size

  // TODO!!!
  //~ new_indiv->_replic_report->compute_statistical_data();
}

void ae_dna::do_small_mutations(void)
{
  // ==============================================================
  //  1. Compute how many rearrangements this genome will undertake
  // ==============================================================
  //
  // Given the rate p (by nucl.) of insertion - for instance -, the number of
  // insertions we perform on the genome follows a binomial law B(n, p), with
  // n = genome length.

  int32_t nb_swi = _indiv->_mut_prng->binomial_random(_length, _indiv->get_point_mutation_rate());
  int32_t nb_ins = _indiv->_mut_prng->binomial_random(_length, _indiv->get_small_insertion_rate());
  int32_t nb_del = _indiv->_mut_prng->binomial_random(_length, _indiv->get_small_deletion_rate());
  int32_t nb_mut = nb_swi + nb_ins + nb_del;



  // ====================================================
  //  2. Perform those small mutations in a random order
  // ====================================================
  //
  // We put the '_nb_small_mutations' mutation events in an "urn". Then we repeat a random drawing
  // of one mutation event in this urn, without replacement, until no mutation
  // event is left in the urn. Here is the "urn" we use at the beginning:
  //
  //     -----------------------------------------------------------
  //    | swi | swi | swi | ins | ins | ins | del | del | del | del |
  //     -----------------------------------------------------------
  //                      ^                 ^                       ^
  //                  nb_swi             nb_swi                   nb_swi
  //                                    +nb_ins                  +nb_ins
  //                                                             +nb_del
  //
  // Random draw of one mutation = random draw of one position in this "urn".
  // Given this position, we know what kind of mutation we have drawn.


  int32_t random_value;
  ae_mutation* mut = NULL;

  for (int32_t i = nb_mut ; i >= 1 ; i--)
  {
    random_value = _indiv->_mut_prng->random(i);

    if (random_value < nb_swi)
    {
      mut = do_switch();
      assert(mut != NULL || !(_exp_m->get_output_m()->get_record_tree() && _exp_m->get_output_m()->get_tree_mode() == NORMAL));

      nb_swi--;  // updating the urn (no replacement!)...
    }
    else if (random_value < nb_swi + nb_ins)
    {
      mut = do_small_insertion();

      nb_ins--;
    }
    else // (random_value >= nb_swi + nb_ins) => del
    {
      mut = do_small_deletion();

      nb_del--;
    }

    // Record mutation in tree
    if (_exp_m->get_output_m()->get_record_tree() && _exp_m->get_output_m()->get_tree_mode() == NORMAL)
    {
      if (mut != NULL)
      {
        _replic_report->add_mut(*mut);
      }
    }
    else
    {
      if (mut != NULL)
      {
        delete mut;
        mut = NULL;
      }
    }
  }
}

void ae_dna::do_rearrangements(void)
{
  // ==============================================================
  //  1. Compute how many rearrangements this genome will undertake
  // ==============================================================
  //
  // Given the rate p (by nucl.) of duplication - for instance -, the number of
  // duplications we perform on the genome follows a binomial law B(n, p), with
  // n = genome length.


  int32_t nb_dupl  = _indiv->_mut_prng->binomial_random(_length, _indiv->get_duplication_rate());
  int32_t nb_del   = _indiv->_mut_prng->binomial_random(_length, _indiv->get_deletion_rate());
  int32_t nb_trans = _indiv->_mut_prng->binomial_random(_length, _indiv->get_translocation_rate());
  int32_t nb_inv   = _indiv->_mut_prng->binomial_random(_length, _indiv->get_inversion_rate());
  int32_t nb_rear  = nb_dupl + nb_del + nb_trans + nb_inv;


  // ===================================================
  //  2. Perform those rearrangements in a random order
  // ===================================================
  //
  // We put the nb_rea rearrangements in an "urn". Then we repeat a random draw
  // of one rearrangement in this urn, without replacement, until no rearrange-
  // -ment is left in the urn. Here is the "urn" we use at the beginning:
  //
  //     ------------------------------------------------------------------------------
  //    | Dupl | Dupl | Dupl | Del | Del | Del | Del | Trans | Trans | Inv | Inv | Inv |
  //     ------------------------------------------------------------------------------
  //                         ^                       ^               ^                 ^
  //                      nb_dupl                 nb_dupl         nb_dupl           nb_dupl
  //                                             +nb_del         +nb_del           +nb_del
  //                                                             +nb_trans         +nb_trans
  //                                                                               +nb_inv
  //
  // Random draw of one rearrangement = random draw of one position in this urn.
  // Given this position, we know what kind of rearrangement we have drawn.

  int32_t random_value;
  ae_mutation* mut = NULL;

  for (int32_t i = nb_rear ; i >= 1 ; i--)
  {
    random_value = _indiv->_mut_prng->random(i);

    if (random_value < nb_dupl)
    {
      mut = do_duplication();
      nb_dupl--;  // Updating the urn (no replacement!)...
    }
    else if (random_value < nb_dupl + nb_del)
    {
      mut = do_deletion();
      nb_del--;
    }
    else if (random_value < nb_dupl + nb_del + nb_trans)
    {
      mut = do_translocation();
      nb_trans--;
    }
    else
    {
      mut = do_inversion();
      nb_inv--;
    }

    // Record rearrangement in tree
    if (_exp_m->get_output_m()->get_record_tree() && _exp_m->get_output_m()->get_tree_mode() == NORMAL)
    {
      if (mut != NULL)
      {
        _replic_report->add_rear(*mut);
      }
    }
    else
    {
      if (mut != NULL)
      {
        delete mut;
        mut = NULL;
      }
    }
  }
}

void ae_dna::do_rearrangements_with_align(void)
{
  bool    direct_sense; // Whether we look for a direct or indirect alignment
  double  rand1 = 0.0;  // Determines the type of rearrangement that will be done if an alignment is found
  int16_t needed_score; // Minimum alignment score needed to recombine (stochastic)
  int32_t seed1, seed2; // Points defining the sequences between which we will look for an alignment

  double ttl = 1.0; // Indiv's Time To Live
  int32_t nb_pairs; // Number of pairs of sequences we will try to align
  int32_t genome_size = _length; // Keep trace of the original length of the genome

  ae_mutation* mut = NULL;
  ae_vis_a_vis* alignment = NULL;

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // For each pair of points to be tested (i.e. while the organism is still "alive"),
  // 1) Draw a random sense (direct or indirect).
  // 2) Determine the minimum alignment score needed for a rearrangement to occur.
  // 3) Test the existence of an alignment with a high enough score.
  // 4) If such an alignment was found, determine the type of rearrangement to be performed
  //    and proceed (WARNING : translocations require another alignment to be found between
  //    the sequence to be translocated and the rest of the chromosome).
  // 5) If there was a change in the chromosome's length, update the individual's TTL and nb_pairs
  //    according to new genome size.
  // 6) If there was a rearrangement, we either save its record in the tree or delete it.
  //////////////////////////////////////////////////////////////////////////////////////////////////
  for (nb_pairs = (int32_t)(ceil(ttl * _length * _indiv->get_neighbourhood_rate()))
      ; nb_pairs > 0
      ; nb_pairs--)
  {
    /////////////////////////////////////////////////
    // 1) Draw a random sense (direct or indirect) //
    /////////////////////////////////////////////////
    direct_sense  = (_indiv->_mut_prng->random() < 0.5); // Determine whether we look for a direct or indirect alignment
    rand1         = _indiv->_mut_prng->random();         // Determine the type of rearrangement to be done. This is an
    // anticipation on 4) for optimization purpose (save computation
    // time if the type of rear is "none"


    //////////////////////////////////////////////////////////////////////////////////
    // 2) Determine the minimum alignment score needed for a rearrangement to occur //
    //////////////////////////////////////////////////////////////////////////////////
    if (_indiv->get_align_fun_shape() == LINEAR)
    {
      needed_score = (int16_t) ceil(_indiv->get_align_lin_min() + _indiv->_mut_prng->random() * (_indiv->get_align_lin_max() - _indiv->get_align_lin_min()));
    }
    else
    {
      // I want the probability of rearrangement for an alignment of score <score> to be
      // prob = 1 / (1 + exp(-(score-mean)/lambda))
      // The score needed for a rearrangement to take place with a given random drawing is hence
      // needed_score = ceil(-lambda * log(1/rand - 1) + mean)
      needed_score = (int16_t) ceil(- _indiv->get_align_sigm_lambda() * log(1/_indiv->_mut_prng->random() - 1) + _indiv->get_align_sigm_mean());
      if (needed_score < 0) needed_score = 0;

      //~ <DEBUG>
      //~ FILE* tmp_file = fopen("scores.out", "a");
      //~ fprintf(tmp_file, "%"PRId16"\n", needed_score);
      //~ fclose(tmp_file);
      //~ </DEBUG>
    }

    // Determine where to look for an alignment (draw seeds)
    seed1 = _indiv->_mut_prng->random(_length);
    seed2 = _indiv->_mut_prng->random(_length);


    if (direct_sense)
    {
      if (rand1 >= _indiv->get_duplication_proportion() + _indiv->get_deletion_proportion() + _indiv->get_translocation_proportion())
      {
        // rand1 corresponds to "no rearrangement" => Nothing to do
        continue;
      }

      ////////////////////////////////////////////////////////////////////
      // 3) Test the existence of an alignment with a high enough score //
      ////////////////////////////////////////////////////////////////////
      alignment = ae_align::search_alignment_direct(this, seed1, this, seed2, needed_score);

      if (alignment == NULL)
      {
        // No alignment found
        continue;
      }

      //~ printf("direct   needed_score : %"PRId32"\n", needed_score);

      ////////////////////////////////////////////////////////////////////////
      // 4) Determine the type of rearrangement to be performed and proceed //
      ////////////////////////////////////////////////////////////////////////
      if (rand1 < _indiv->get_duplication_proportion())
      {
        // Remember the length of the segment to be duplicated and of the genome before the duplication
        int32_t segment_length  = ae_utils::mod(alignment->get_i_2() -
                                                alignment->get_i_1(),
                                                _length);
        int32_t gu_size_before  = _length;
        int32_t gu_size_after   = gu_size_before + segment_length;
        int32_t genome_size_before = _indiv->get_amount_of_dna();
        int32_t genome_size_after = genome_size_before + segment_length;

        if ((genome_size_after > _indiv->get_max_genome_length()) ||
            (gu_size_after > _gen_unit->get_max_gu_length()))
        {
          if (_exp_m->get_output_m()->is_logged(LOG_BARRIER))
          {
            // Write an entry in the barrier log file
            fprintf(_exp_m->get_output_m()->get_log(LOG_BARRIER),
                    "%" PRId64 " %" PRId32 " DUPLICATION %" PRId32 " %" PRId32
                        " %" PRId32 " %" PRId32 "\n",
                    Time::get_time(), _indiv->get_id(), segment_length, 0,
                    gu_size_before, genome_size_before);
          }
        }
        else
        {
          // Perform in situ (tandem) DUPLICATION
          do_duplication(alignment->get_i_1(),
                         alignment->get_i_2(),
                         alignment->get_i_2());

          // Report the duplication
          if (_exp_m->get_output_m()->get_record_tree() &&
              _exp_m->get_output_m()->get_tree_mode() == NORMAL)
          {
            // Report the insertion
            mut = new ae_mutation();
            mut->report_duplication(alignment->get_i_1(),
                                    alignment->get_i_2(),
                                    alignment->get_i_2(),
                                    segment_length, needed_score);
          }

          // Write a line in rearrangement logfile
          if (_exp_m->get_output_m()->is_logged(LOG_REAR))
          {
            fprintf(_exp_m->get_output_m()->get_log(LOG_REAR),
                    "%" PRId64 " %" PRId32 " %" PRId8 " %" PRId32 " %" PRId32
                        " %" PRId16 "\n",
                    Time::get_time(), _indiv->get_id(), int8_t(DUPL),
                    segment_length, genome_size_before, needed_score);
          }
        }
      }
      else if (rand1 < _indiv->get_duplication_proportion() +
                       _indiv->get_deletion_proportion())
      {
        // Remember the length of the segment to be duplicated and of the genome before the deletion
        int32_t segment_length  = ae_utils::mod(alignment->get_i_2() -
                                                alignment->get_i_1() - 1,
                                                _length) + 1;
        int32_t gu_size_before  = _length;
        int32_t gu_size_after   = gu_size_before - segment_length;
        int32_t genome_size_before = _indiv->get_amount_of_dna();
        int32_t genome_size_after = genome_size_before - _length;

        if ((genome_size_after < _indiv->get_min_genome_length()) ||
            (gu_size_after < _gen_unit->get_min_gu_length()))
        {
          if (_exp_m->get_output_m()->is_logged(LOG_BARRIER))
          {
            // Write an entry in the barrier log file
            fprintf(_exp_m->get_output_m()->get_log(LOG_BARRIER),
                    "%" PRId64 " %" PRId32 " DELETION %" PRId32 " %" PRId32
                        " %" PRId32 " %" PRId32 "\n",
                    Time::get_time(), _indiv->get_id(), segment_length, 0,
                    gu_size_before, genome_size_before);
          }
        }
        else
        {
          // Perform DELETION
          do_deletion(alignment->get_i_1(), alignment->get_i_2());

          // Report the deletion
          if (_exp_m->get_output_m()->get_record_tree() &&
              _exp_m->get_output_m()->get_tree_mode() == NORMAL)
          {
            // Report the insertion
            mut = new ae_mutation();
            mut->report_deletion(alignment->get_i_1(),
                                 alignment->get_i_2(),
                                 segment_length, needed_score);
          }

          // Write a line in rearrangement logfile
          if (_exp_m->get_output_m()->is_logged(LOG_REAR))
          {
            fprintf(_exp_m->get_output_m()->get_log(LOG_REAR),
                    "%" PRId64 " %" PRId32 " %" PRId8 " %" PRId32 " %" PRId32
                        " %" PRId16 "\n",
                    Time::get_time(), _indiv->get_id(), int8_t(DEL),
                    segment_length, genome_size_before, needed_score);
          }
        }
      }
      else
      {
        assert(rand1 < _indiv->get_duplication_proportion() +
                       _indiv->get_deletion_proportion() +
                       _indiv->get_translocation_proportion());

        // Perform TRANSLOCATION
        // Make sure the segment to be translocated doesn't contain OriC TODO : is that still necessary?
        if (alignment->get_i_1() > alignment->get_i_2())
        {
          alignment->swap();
        }

        // Remember the length of the segment to be translocated
        int32_t segment_length = ae_utils::mod(alignment->get_i_2() - alignment->get_i_1(), _length);

        // Extract the segment to be translocated
        GeneticUnit* translocated_segment = extract_into_new_GU(alignment->get_i_1(), alignment->get_i_2());

        // Look for an alignment between the segment to be translocated and the rest of the genome
        bool direct_sense;
        int16_t needed_score_2;
        ae_vis_a_vis* alignment_2 = NULL;
        int32_t seed1, seed2;
        for (nb_pairs = (int32_t)(ceil(ttl * _length * _indiv->get_neighbourhood_rate()))
            ; nb_pairs > 0
            ; nb_pairs--)
        {
          direct_sense = (_indiv->_mut_prng->random() < 0.5);

          if (_indiv->get_align_fun_shape() == LINEAR)
          {
            needed_score_2  = (int16_t) ceil(_indiv->get_align_lin_min() + _indiv->_mut_prng->random() * (_indiv->get_align_lin_max() - _indiv->get_align_lin_min()));
          }
          else
          {
            needed_score_2 = (int16_t) ceil(- _indiv->get_align_sigm_lambda() * log(1/_indiv->_mut_prng->random() - 1) + _indiv->get_align_sigm_mean());
            if (needed_score_2 < 0) needed_score_2 = 0;
          }

          seed1 = _indiv->_mut_prng->random(_length);
          seed2 = _indiv->_mut_prng->random(segment_length);

          if (direct_sense)
          {
            alignment_2 = ae_align::search_alignment_direct(this, seed1, translocated_segment->get_dna(), seed2, needed_score_2);
          }
          else // if indirect
          {
            alignment_2 = ae_align::search_alignment_indirect(this, seed1, translocated_segment->get_dna(), seed2, needed_score_2);
          }

          if (alignment_2 != NULL)
          {
            //~ printf("transloc needed_score : %"PRId32"\n", needed_score_2);
            break;
          }
        }


        // If an alignment was found between the segment to be translocated and the rest of the genome, proceed to the
        // translocation, otherwise, replace the extracted segment at its former position (cancel the translocation event).
        if (alignment_2 != NULL)
        {
          // Proceed to the translocation
          insert_GU(translocated_segment,
                    alignment_2->get_i_1(), alignment_2->get_i_2(),
                    (alignment_2->get_sense() == INDIRECT));

          // Report the translocation
          if (_exp_m->get_output_m()->get_record_tree() &&
              _exp_m->get_output_m()->get_tree_mode() == NORMAL)
          {
            mut = new ae_mutation();
            mut->report_translocation(alignment->get_i_1(),
                                      alignment->get_i_2(),
                                      alignment_2->get_i_1(),
                                      alignment_2->get_i_2(),
                                      segment_length,
                                      (alignment_2->get_sense() == INDIRECT),
                                      needed_score, needed_score_2);
          }

          // Write a line in rearrangement logfile
          if (_exp_m->get_output_m()->is_logged(LOG_REAR))
          {
            fprintf(_exp_m->get_output_m()->get_log(LOG_REAR),
                    "%" PRId64 " %" PRId32 " %" PRId8 " %" PRId32 " %" PRId32
                        " %" PRId16 "\n",
                    Time::get_time(), _indiv->get_id(), int8_t(TRANS),
                    segment_length, _length, needed_score_2);
          }
          delete alignment_2;
        }
        else
        {
          // Cancel the translocation (replace the extracted segment at its former position)
          insert_GU(translocated_segment, alignment->get_i_1(), 0, false);
        }

        delete translocated_segment;
      }

      delete alignment;
    }
    else // if indirect
    {
      if (rand1 >= _indiv->get_inversion_proportion())
      {
        // rand1 corresponds to no rearrangement => Nothing to do
        continue;
      }

      //~ printf("indirect needed_score : %"PRId32"\n", needed_score);

      ////////////////////////////////////////////////////////////////////
      // 3) Test the existence of an alignment with a high enough score //
      ////////////////////////////////////////////////////////////////////
      alignment = ae_align::search_alignment_indirect(this, seed1, this, seed2, needed_score);

      if (alignment == NULL)
      {
        // No alignment found
        continue;
      }

      /////////////////////////////
      // 4) Proceed to inversion //
      /////////////////////////////
      // Make sure the segment to be inverted doesn't contain OriC
      if (alignment->get_i_1() > alignment->get_i_2())
      {
        alignment->swap();
      }

      // Remember the length of the segment to be duplicated
      int32_t segment_length = ae_utils::mod(alignment->get_i_2() -
                                             alignment->get_i_1(),
                                             _length);

      // Proceed
      do_inversion(alignment->get_i_1(), alignment->get_i_2());

      // Report the inversion
      if (_exp_m->get_output_m()->get_record_tree() &&
          _exp_m->get_output_m()->get_tree_mode() == NORMAL)
      {
        mut = new ae_mutation();
        mut->report_inversion(alignment->get_i_1(),
                              alignment->get_i_2(),
                              segment_length, needed_score);
      }

      // Write a line in rearrangement logfile
      if (_exp_m->get_output_m()->is_logged(LOG_REAR))
      {
        fprintf(_exp_m->get_output_m()->get_log(LOG_REAR),
                "%" PRId64 " %" PRId32 " %" PRId8 " %" PRId32 " %" PRId32
                    " %" PRId16 "\n",
                Time::get_time(), _indiv->get_id(), int8_t(INV), segment_length,
                _length, needed_score);
      }

      delete alignment;
    }


    //////////////////////////////////////////////////////////////////////////////
    // 5) If there was a change in the chromosome's length,                     //
    //    update the individual's TTL and nb_pairs according to new genome size //
    //////////////////////////////////////////////////////////////////////////////
    if (genome_size != _length)
    {
      ttl = ((double)(nb_pairs-1)) / ((double)genome_size) / _indiv->get_neighbourhood_rate();
      genome_size = _length;
      nb_pairs = (int32_t)(ceil(ttl * _length * _indiv->get_neighbourhood_rate())) + 1;
    }

    //////////////////////////////////////////////////////////////////////////////////////////
    // 6) If there was a rearrangement, we either save its record in the tree or delete it. //
    //////////////////////////////////////////////////////////////////////////////////////////
    if (mut != NULL)
    {
      if (_exp_m->get_output_m()->get_record_tree() && _exp_m->get_output_m()->get_tree_mode() == NORMAL)
      {
        _replic_report->add_rear(*mut);
        mut = NULL;
      }
      else
      {
        delete mut;
        mut = NULL;
      }
    }
  }
}

void ae_dna::do_transfer(int32_t parent_id)
{
  ae_mutation* mut = NULL;
  if (_indiv->get_mut_prng()->random() < _indiv->get_HT_ins_rate())
  {
    mut = do_ins_HT(parent_id);
    if (_exp_m->get_output_m()->get_record_tree() && _exp_m->get_output_m()->get_tree_mode() == NORMAL)
    {
      if (mut != NULL)
      {
        _replic_report->add_HT(*mut);
      }
    }
    else
    {
      if (mut != NULL)
      {
        delete mut;
        mut = NULL;
      }
    }
  }

  if (_indiv->get_mut_prng()->random() < _indiv->get_HT_repl_rate())
  {
    mut = do_repl_HT(parent_id);
    if (_exp_m->get_output_m()->get_record_tree() && _exp_m->get_output_m()->get_tree_mode() == NORMAL)
    {
      if (mut != NULL)
      {
        _replic_report->add_HT(*mut);
      }
    }
    else
    {
      if (mut != NULL)
      {
        delete mut;
        mut = NULL;
      }
    }
  }
}

ae_mutation* ae_dna::do_switch(void)
{
  ae_mutation* mut = NULL;

  int32_t pos = _indiv->_mut_prng->random(_length);

  if (do_switch(pos))
  {
    if (_exp_m->get_output_m()->get_record_tree() && _exp_m->get_output_m()->get_tree_mode() == NORMAL)
    {
      // Report the mutation
      mut = new ae_mutation();
      mut->report_point_mutation(pos);
    }
  }

  return mut;
}

ae_mutation* ae_dna::do_small_insertion(void)
{
  ae_mutation* mut = NULL;

  // Determine the position and size of the small insertion
  int32_t pos = _indiv->_mut_prng->random(_length);
  int16_t nb_insert;
  if (_indiv->get_max_indel_size() == 1)
  {
    nb_insert = 1;
  }
  else
  {
    nb_insert = 1 + _indiv->_mut_prng->random(_indiv->get_max_indel_size());
    // <nb_insert> must be in [1 ; max_indel_size]
  }

  // Check that the insertion won't throw the genome size over the limit
  if ((_indiv->get_amount_of_dna() + nb_insert >
       _indiv->get_max_genome_length()) ||
      (_length + nb_insert > _gen_unit->get_max_gu_length()))
  {
    if (_exp_m->get_output_m()->is_logged(LOG_BARRIER))
    {
      // Write an entry in the barrier log file
      fprintf(_exp_m->get_output_m()->get_log(LOG_BARRIER),
              "%" PRId64 " %" PRId32 " S_INS %" PRId32 " %" PRId32 " %" PRId32
                  " %" PRId32 "\n",
              Time::get_time(), _indiv->get_id(), nb_insert, 0, _length,
              _indiv->get_amount_of_dna());
    }

    return NULL;
  }

  // Prepare the sequence to be inserted
  char* inserted_seq = new char[nb_insert + 1];
  char  inserted_char;
  for (int16_t j = 0 ; j < nb_insert ; j++)
  {
    inserted_char = (char) '0' + _indiv->_mut_prng->random(NB_BASE);
    inserted_seq[j] = inserted_char;
  }
  inserted_seq[nb_insert] = '\0';

  // Proceed to the insertion and report it
  if (do_small_insertion(pos, nb_insert, inserted_seq))
  {
    if (_exp_m->get_output_m()->get_record_tree() &&
        _exp_m->get_output_m()->get_tree_mode() == NORMAL)
    {
      // Report the insertion
      mut = new ae_mutation();
      mut->report_small_insertion(pos, nb_insert, inserted_seq);
    }
  }

  // Delete the sequence
  delete [] inserted_seq;

  return mut;
}

ae_mutation* ae_dna::do_small_deletion(void)
{
  ae_mutation* mut = NULL;

  // Determine the position and size of the small deletion
  int32_t pos = _indiv->_mut_prng->random(_length);
  int16_t nb_del;
  if (_indiv->get_max_indel_size() == 1)
  {
    nb_del = 1;
  }
  else
  {
    nb_del = 1 + _indiv->_mut_prng->random(_indiv->get_max_indel_size());
    // <nb_del> must be in [1 ; max_indel_size]
  }

  // Check that the insertion won't shrink the genome size under the limit nor to nothing
  if ((_indiv->get_amount_of_dna() - nb_del <
       _indiv->get_min_genome_length()) ||
      (_length - nb_del < _gen_unit->get_min_gu_length()))
  {
    if (_exp_m->get_output_m()->is_logged(LOG_BARRIER))
    {
      // Write an entry in the barrier log file
      fprintf(_exp_m->get_output_m()->get_log(LOG_BARRIER),
              "%" PRId64 " %" PRId32 " S_DEL %" PRId32 " %" PRId32
                  " %" PRId32 " %" PRId32 "\n",
              Time::get_time(), _indiv->get_id(), nb_del, 0, _length,
              _indiv->get_amount_of_dna());
    }

    return NULL;
  }

  if (do_small_deletion(pos, nb_del))
  {
    if (_exp_m->get_output_m()->get_record_tree() && _exp_m->get_output_m()->get_tree_mode() == NORMAL)
    {
      mut = new ae_mutation();
      mut->report_small_deletion(pos, nb_del);
    }
  }

  return mut;
}

bool ae_dna::do_switch(int32_t pos)
{
  // Perform the mutation
  if (_data[pos] == '0')  _data[pos] = '1';
  else                      _data[pos] = '0';

  // Remove promoters containing the switched base
  _gen_unit->remove_promoters_around(pos, ae_utils::mod(pos + 1, _length));

  // Look for potential new promoters containing the switched base
  if (_length >= PROM_SIZE)
  {
    _gen_unit->look_for_new_promoters_around(pos, ae_utils::mod(pos + 1, _length));
  }


  return true;
}

bool ae_dna::do_small_insertion(int32_t pos, int16_t nb_insert, char * seq)
{
  // Check genome size limit
  assert(_length + nb_insert <= _gen_unit->get_max_gu_length());
  assert(_indiv->get_amount_of_dna() + nb_insert <= _indiv->get_max_genome_length());

  // Remove the promoters that will be broken
  _gen_unit->remove_promoters_around(pos);

  // Insert the sequence
  insert(pos, seq, nb_insert);

  // Look for new promoters
  if (_length >= PROM_SIZE)
  {
    if (_length - nb_insert < PROM_SIZE)
    {
      // Special case where the genome was smaller than a promoter before the insertion and greater than (or as big as) a promoter after the insertion.
      // In that case, we must look for new promoters thoroughly on the whole genome using locate_promoters
      _gen_unit->locate_promoters();
    }
    else
    {
      _gen_unit->move_all_promoters_after(pos, nb_insert);
      _gen_unit->look_for_new_promoters_around(pos, ae_utils::mod(pos + nb_insert, _length));
    }
  }

  return true;
}

bool ae_dna::do_small_deletion(int32_t pos, int16_t nb_del)
{
  // Check genome size limit
  assert(_length - nb_del >= _gen_unit->get_min_gu_length());
  assert(_indiv->get_amount_of_dna() - nb_del >= _indiv->get_min_genome_length());

  // Remove promoters containing at least one nucleotide from the sequence to delete
  _gen_unit->remove_promoters_around(pos, ae_utils::mod(pos + nb_del, _length));

  // Do the deletion and update promoter list
  if (pos + nb_del <= _length) // the deletion does not contain the replication origin
  {
    // Do the deletion
    remove(pos, pos + nb_del);

    // Update promoter list
    if (_length >= PROM_SIZE)
    {
      _gen_unit->move_all_promoters_after(pos, -nb_del);
      _gen_unit->look_for_new_promoters_around(ae_utils::mod(pos, _length));
    }
  }
  else // the deletion contains the replication origin
  {
    // Do the deletion
    int32_t nb_del_at_pos_0 = nb_del - _length + pos;
    remove(pos, _length);
    remove(0, nb_del_at_pos_0);
    pos -= nb_del_at_pos_0;

    // Update promoter list
    if (_length >= PROM_SIZE)
    {
      _gen_unit->move_all_promoters_after(0, -nb_del_at_pos_0);
      _gen_unit->look_for_new_promoters_around(0);
    }
  }

  return true;
}

ae_mutation* ae_dna::do_duplication(void)
{
  ae_mutation* mut = NULL;

  int32_t pos_1, pos_2, pos_3;
  pos_1 = _indiv->_mut_prng->random(_length);
  pos_2 = _indiv->_mut_prng->random(_length);
  pos_3 = _indiv->_mut_prng->random(_length);

  // Remember the length of the segment to be duplicated and of the former genome
  int32_t segment_length      = ae_utils::mod(pos_2 - pos_1 - 1, _length) + 1;
  int32_t gu_size_before  = _length;
  int32_t gu_size_after   = gu_size_before + segment_length;
  int32_t genome_size_before = _indiv->get_amount_of_dna();
  int32_t genome_size_after = genome_size_before + segment_length;
  if ((gu_size_after > _gen_unit->get_max_gu_length()) ||
      (genome_size_after > _indiv->get_max_genome_length()))
  {
    if (_exp_m->get_output_m()->is_logged(LOG_BARRIER))
    {
      // Write an entry in the barrier log file
      fprintf(_exp_m->get_output_m()->get_log(LOG_BARRIER),
              "%" PRId64 " %" PRId32 " DUPLICATION %" PRId32 " %" PRId32
                  " %" PRId32 " %" PRId32 "\n",
              Time::get_time(), _indiv->get_id(), segment_length, 0,
              gu_size_before, genome_size_before);
    }
  }
  else
  {
    // Perform the duplication
    do_duplication(pos_1, pos_2, pos_3);

    // Report the duplication
    if (_exp_m->get_output_m()->get_record_tree() &&
        _exp_m->get_output_m()->get_tree_mode() == NORMAL)
    {
      mut = new ae_mutation();
      mut->report_duplication(pos_1, pos_2, pos_3, segment_length);
    }

    // Write a line in rearrangement logfile
    if (_exp_m->get_output_m()->is_logged(LOG_REAR))
    {
      fprintf(_exp_m->get_output_m()->get_log(LOG_REAR),
              "%" PRId64 " %" PRId32 " %" PRId8 " %" PRId32 " %" PRId32 "\n",
              Time::get_time(), _indiv->get_id(), int8_t(DUPL), segment_length,
              genome_size_before);
    }
  }

  return mut;
}

ae_mutation* ae_dna::do_deletion(void)
{
  ae_mutation* mut = NULL;

  int32_t pos_1, pos_2;
  pos_1 = _indiv->_mut_prng->random(_length);
  pos_2 = _indiv->_mut_prng->random(_length);

  // Remember the length of the segment to be deleted and of the genome before the deletion
  int32_t segment_length  = ae_utils::mod(pos_2 - pos_1 - 1, _length) + 1;
  int32_t gu_size_before  = _length;
  int32_t gu_size_after   = gu_size_before - segment_length;
  int32_t genome_size_before = _indiv->get_amount_of_dna();
  int32_t genome_size_after = genome_size_before - segment_length;


  if ((gu_size_after < _gen_unit->get_min_gu_length()) ||
      (genome_size_after < _indiv->get_min_genome_length()))
  {
    if (_exp_m->get_output_m()->is_logged(LOG_BARRIER))
    {
      // Write an entry in the barrier log file
      fprintf(_exp_m->get_output_m()->get_log(LOG_BARRIER),
              "%" PRId64 " %" PRId32 " DELETION %" PRId32 " %" PRId32
                  " %" PRId32 " %" PRId32 "\n",
              Time::get_time(), _indiv->get_id(), segment_length, 0,
              gu_size_before, genome_size_before);
    }
  }
  else
  {
    // Perform the deletion
    do_deletion(pos_1, pos_2);

    // Report the deletion
    if (_exp_m->get_output_m()->get_record_tree() &&
        _exp_m->get_output_m()->get_tree_mode() == NORMAL)
    {
      mut = new ae_mutation();
      mut->report_deletion(pos_1, pos_2, segment_length);
    }

    // Write a line in rearrangement logfile
    if (_exp_m->get_output_m()->is_logged(LOG_REAR))
    {
      fprintf(_exp_m->get_output_m()->get_log(LOG_REAR),
              "%" PRId64 " %" PRId32 " %" PRId8 " %" PRId32 " %" PRId32 "\n",
              Time::get_time(), _indiv->get_id(), int8_t(DEL), segment_length,
              genome_size_before);
    }
  }

  return mut;
}

ae_mutation* ae_dna::do_translocation(void)
{

  ae_mutation* mut = NULL;

  int32_t pos_1, pos_2, pos_3, pos_4;
  int32_t segment_length;
  bool invert;

  if (_indiv->get_allow_plasmids())
  {
    // -----------------------------------------------------------------
    // WARNING : This is only valid when there is only 1 plasmid allowed
    // -----------------------------------------------------------------
    int32_t pos_1_rel, pos_2_rel, pos_3_rel, pos_4_rel;

    Individual * indiv = _indiv;
    const GeneticUnit* chromosome = &indiv->get_genetic_unit_list().front();
    const GeneticUnit* plasmid    = &*std::next(indiv->get_genetic_unit_list().begin());
    int32_t chrom_length        = chromosome->get_dna()->get_length();
    int32_t total_amount_of_dna = indiv->get_amount_of_dna();

    // 1) What sequence are we translocating?
    pos_1_rel = _indiv->_mut_prng->random(_length);
    pos_2_rel = _indiv->_mut_prng->random(_length);

    int32_t segment_length = ae_utils::mod(pos_2_rel - pos_1_rel, _length);

    pos_3_rel = ae_utils::mod(pos_1_rel + _indiv->_mut_prng->random(segment_length), _length);

    if (_gen_unit == chromosome)
    {
      pos_1 = pos_1_rel;
      pos_2 = pos_2_rel;
      pos_3 = pos_3_rel;
    }
    else // (_gen_unit == plasmid)
    {
      pos_1 = pos_1_rel + chrom_length;
      pos_2 = pos_2_rel + chrom_length;
      pos_3 = pos_3_rel + chrom_length;
    }


    // 2) Where are we translocating it?
    pos_4 = _indiv->_mut_prng->random(total_amount_of_dna - segment_length);

    if (_gen_unit == chromosome)
    {
      if (pos_1 <= pos_2)
      {
        if (pos_4 >= pos_1)
        {
          pos_4 += segment_length;
        }
      }
      else
      {
        if (pos_4 >= chrom_length - segment_length)
        {
          pos_4 += segment_length;
        }
        else
        {
          pos_4 += pos_2;
        }
      }
      if (pos_4 >= chrom_length)
      {
        pos_4_rel = pos_4 - chrom_length;
      }
      else
      {
        pos_4_rel = pos_4;
      }
    }
    else // (_gen_unit == plasmid)
    {
      if (pos_1 <= pos_2)
      {
        if (pos_4 >= pos_1)
        {
          pos_4 += segment_length;
        }
      }
      else
      {
        if (pos_4 >= chrom_length)
        {
          pos_4 += pos_2_rel;
        }
      }


      if (pos_4 >= chrom_length)
      {
        pos_4_rel = pos_4 - chrom_length;
      }
      else
      {
        pos_4_rel = pos_4;
      }

    }

    invert = (_indiv->_mut_prng->random(2) == 0);


    // If inter GU translocation
    if ((_gen_unit == chromosome && pos_4 >= chrom_length) ||
        (_gen_unit == plasmid && pos_4 < chrom_length))
    {
      // printf("Translocation from the %s to the %s\n", _gen_unit==chromosome?"chromosome":"plasmid", _gen_unit==chromosome?"plasmid":"chromosome");

      //~ printf("  Chromosome length : %"PRId32"     Plasmid length : %"PRId32"     Total length : %"PRId32"\n",
      //~ chrom_length, plasmid->get_dna()->get_length(), chrom_length + plasmid->get_dna()->get_length());
      //~ printf("  %"PRId32" %"PRId32" %"Pthere are plasmidsRId32" %"PRId32" (length : %"PRId32")\n", pos_1, pos_2, pos_3, pos_4, segment_length);
      //~ printf("  former pos : %"PRId32" %"PRId32" %"PRId32" %"PRId32"\n", former_pos_1, former_pos_2, former_pos_3, former_pos_4);

      if (do_inter_GU_translocation(pos_1_rel, pos_2_rel,
                                    pos_3_rel, pos_4_rel, invert))
      {
        // Report the translocation
        if (_exp_m->get_output_m()->get_record_tree() &&
            _exp_m->get_output_m()->get_tree_mode() == NORMAL)
        {
          mut = new ae_mutation();
          mut->report_translocation(pos_1_rel, pos_2_rel,
                                    pos_3_rel, pos_4_rel,
                                    segment_length, invert);
        }

        // Write a line in rearrangement logfile
        if (_exp_m->get_output_m()->is_logged(LOG_REAR))
        {
          fprintf(_exp_m->get_output_m()->get_log(LOG_REAR),
                  "%" PRId64 " %" PRId32 " %" PRId8 " %" PRId32 " %" PRId32 "\n",
                  Time::get_time(), _indiv->get_id(), int8_t(TRANS),
                  segment_length, _length);
        }
      }
    }
    else
    {
      if (pos_1_rel > pos_2_rel)
      {

      }
      // printf("Translocation intra %s\n", _gen_unit==chromosome?"chromosome":"plasmid");
      //~ printf("  Chromosome length : %"PRId32"     Plasmid length : %"PRId32"     Total length : %"PRId32"\n",
      //~ chrom_length, plasmid->get_dna()->get_length(), chrom_length + plasmid->get_dna()->get_length());
      //~ printf("  %"PRId32" %"PRId32" %"PRId32" %"PRId32" (length : %"PRId32")\n", pos_1, pos_2, pos_3, pos_4, segment_length);
      //~ printf("  former pos : %"PRId32" %"PRId32" %"PRId32" %"PRId32"\n", former_pos_1, former_pos_2, former_pos_3, former_pos_4);
      if (do_translocation(pos_1_rel, pos_2_rel, pos_3_rel, pos_4_rel, invert))
      {
        // Report the translocation
        if (_exp_m->get_output_m()->get_record_tree() &&
            _exp_m->get_output_m()->get_tree_mode() == NORMAL)
        {
          mut = new ae_mutation();
          mut->report_translocation(pos_1_rel, pos_2_rel,
                                    pos_3_rel, pos_4_rel,
                                    segment_length, invert);
        }

        // Write a line in rearrangement logfile
        if (_exp_m->get_output_m()->is_logged(LOG_REAR))
        {
          fprintf(_exp_m->get_output_m()->get_log(LOG_REAR),
                  "%" PRId64 " %" PRId32 " %" PRId8 " %" PRId32 " %" PRId32 "\n",
                  Time::get_time(), _indiv->get_id(), int8_t(TRANS),
                  segment_length, _length);
        }
      }
    }
  }
  else // (! ae_common::params->get_allow_plasmids())
  {
    pos_1 = _indiv->_mut_prng->random(_length);
    pos_2 = _indiv->_mut_prng->random(_length);
    if (pos_1 == pos_2) return NULL;

    // As it will be seen in do_translocation(int32_t pos_1, int32_t pos_2, int32_t pos_3, int32_t pos_4, bool invert),
    // translocating segment [pos_1, pos_2] is the same as translocating segment  [pos_2, pos_1]
    // Since OriC must be at position 0, we will always translocate segment [pos_1, pos_2] with pos_1 < pos_2
    if (pos_1 > pos_2) ae_utils::exchange(pos_1, pos_2);

    segment_length = pos_2 - pos_1;

    // Generate a position between pos_1 and pos_2
    pos_3 = pos_1 + _indiv->_mut_prng->random(segment_length);

    // Generate a position that is NOT between pos_1 and pos_2
    pos_4 = _indiv->_mut_prng->random(_length - segment_length);
    if (pos_4 >= pos_1) pos_4 += segment_length;

    invert = (_indiv->_mut_prng->random(2) == 0);

    if (do_translocation(pos_1, pos_2, pos_3, pos_4, invert))
    {
      // Report the translocation
      if (_exp_m->get_output_m()->get_record_tree() &&
          _exp_m->get_output_m()->get_tree_mode() == NORMAL)
      {
        mut = new ae_mutation();
        mut->report_translocation(pos_1, pos_2, pos_3, pos_4,
                                  segment_length, invert);
      }

      // Write a line in rearrangement logfile
      if (_exp_m->get_output_m()->is_logged(LOG_REAR))
      {
        fprintf(_exp_m->get_output_m()->get_log(LOG_REAR),
                "%" PRId64 " %" PRId32 " %" PRId8 " %" PRId32 " %" PRId32 "\n",
                Time::get_time(), _indiv->get_id(), int8_t(TRANS),
                segment_length, _length);
      }
    }
  }

  return mut;
}

ae_mutation* ae_dna::do_inversion(void)
{
  ae_mutation* mut = NULL;

  int32_t pos_1, pos_2;
  int32_t segment_length;
  pos_1 = _indiv->_mut_prng->random(_length);
  pos_2 = _indiv->_mut_prng->random(_length);

  if (pos_1 == pos_2) return NULL; // Invert everything <=> Invert nothing!
  if (pos_1 >  pos_2) ae_utils::exchange(pos_1, pos_2); // Invert the segment that don't contain OriC

  segment_length = pos_2 - pos_1;

  if(do_inversion(pos_1, pos_2))
  {
    // Report the inversion
    if (_exp_m->get_output_m()->get_record_tree() &&
        _exp_m->get_output_m()->get_tree_mode() == NORMAL)
    {
      mut = new ae_mutation();
      mut->report_inversion(pos_1, pos_2, segment_length);
    }

    // Write a line in rearrangement logfile
    if (_exp_m->get_output_m()->is_logged(LOG_REAR))
    {
      fprintf(_exp_m->get_output_m()->get_log(LOG_REAR),
              "%" PRId64 " %" PRId32 " %" PRId8 " %" PRId32 " %" PRId32 "\n",
              Time::get_time(), _indiv->get_id(), int8_t(INV),
              segment_length, _length);
    }
  }

  return mut;
}

ae_mutation* ae_dna::do_insertion(const char* seq_to_insert, int32_t seq_length /*= -1*/)
{
  ae_mutation* mut = NULL;

  // Compute seq_length if not known
  if (seq_length == -1)
  {
    seq_length = strlen(seq_to_insert);
  }

  // Where to insert the sequence
  int32_t pos = _indiv->_mut_prng->random(_length);

  if (do_insertion(pos, seq_to_insert, seq_length))
  {
    if (_exp_m->get_output_m()->get_record_tree() && _exp_m->get_output_m()->get_tree_mode() == NORMAL)
    {
      // Report the insertion
      mut = new ae_mutation();
      mut->report_insertion(pos, seq_length, seq_to_insert);
    }
  }

  return mut;
}


bool ae_dna::do_duplication(int32_t pos_1, int32_t pos_2, int32_t pos_3)
// Duplicate segment [pos_1; pos_2[ and insert the duplicate before pos_3
{
  char* duplicate_segment = NULL;
  int32_t seg_length;

  if (pos_1 < pos_2)
  {
    //
    //       pos_1         pos_2                   -> 0-
    //         |             |                   -       -
    // 0--------------------------------->      -         -
    //         ===============                  -         - pos_1
    //           tmp (copy)                      -       -
    //                                             -----      |
    //                                             pos_2    <-'
    //

    seg_length = pos_2 - pos_1;
    duplicate_segment = new char[seg_length+1];
    memcpy(duplicate_segment, &_data[pos_1], seg_length);
    duplicate_segment[seg_length] = '\0';
  }
  else // if (pos_1 >= pos_2)
  {
    // The segment to duplicate includes the replication origin.
    // The copying process will be done in two steps.
    //
    //                                            ,->
    //    pos_2                 pos_1            |      -> 0-
    //      |                     |                   -       - pos_2
    // 0--------------------------------->     pos_1 -         -
    // ======                     =======            -         -
    //  tmp2                        tmp1              -       -
    //                                                  -----
    //
    //

    int32_t tmp1_len = _length - pos_1;
    int32_t tmp2_len = pos_2;
    seg_length = tmp1_len + tmp2_len;
    duplicate_segment = new char[seg_length+1];
    memcpy(duplicate_segment, &_data[pos_1], tmp1_len);     // Copy tmp1
    memcpy(&duplicate_segment[tmp1_len], _data, tmp2_len);  // Copy tmp2
    duplicate_segment[seg_length] = '\0';
  }

  // Create a copy of the promoters beared by the segment to be duplicated
  // (they will be inserted in the individual's RNA list later)

  Promoters2Strands duplicated_promoters = {{},{}};
  _gen_unit->duplicate_promoters_included_in(pos_1, pos_2, duplicated_promoters);

  _gen_unit->remove_promoters_around(pos_3);

  insert(pos_3, duplicate_segment, seg_length);

  if (_length >= PROM_SIZE)
  {
    if (_length - seg_length < PROM_SIZE)
    {
      // Special case where the genome was smaller than a promoter before the insertion and greater than (or as big as) a promoter after the insertion.
      // In that case, we must look for new promoters thoroughly on the whole genome using locate_promoters
      _gen_unit->locate_promoters();
    }
    else
    {
      _gen_unit->move_all_promoters_after(pos_3, seg_length);

      _gen_unit->insert_promoters_at(duplicated_promoters, pos_3);

      _gen_unit->look_for_new_promoters_around(pos_3);
      _gen_unit->look_for_new_promoters_around(pos_3 + seg_length);
    }
  }


  delete [] duplicate_segment;

  return true;
}

bool ae_dna::do_deletion(int32_t pos_1, int32_t pos_2)
// Delete segment going from pos_1 (included) to pos_2 (excluded)
{
  if (pos_1 < pos_2)
  {
    //
    //       pos_1         pos_2                   -> 0-
    //         |             |                   -       -
    // 0--------------------------------->      -         -
    //         ===============                  -         - pos_1
    //           tmp (copy)                      -       -
    //                                             -----      |
    //                                             pos_2    <-'
    //

    int32_t segment_length = pos_2 - pos_1;

    // Remove promoters containing at least one nucleotide from the sequence to delete
    _gen_unit->remove_promoters_around(pos_1, pos_2);

    // Delete the sequence between pos_1 and pos_2
    remove(pos_1, pos_2);

    // Update promoter list
    if (_length >= PROM_SIZE)
    {
      _gen_unit->move_all_promoters_after(pos_1, -segment_length);

      _gen_unit->look_for_new_promoters_around(pos_1);
    }
  }
  else // if (pos_1 >= pos_2)
  {
    // The segment to delete includes the replication origin.
    // The deletion process will be done in two steps.
    //
    //                                            ,->
    //    pos_2                 pos_1            |      -> 0-
    //      |                     |                   -       - pos_2
    // 0--------------------------------->     pos_1 -         -
    // =====                      =======            -         -
    //  tmp2                        tmp1              -       -
    //                                                  -----
    //
    //

    //int32_t segment_length = _length + pos_2 - pos_1; //useless variable

    // Remove promoters containing at least one nucleotide from the sequence to delete
    _gen_unit->remove_promoters_around(pos_1, pos_2);

    // Delete the sequence between pos_1 and pos_2
    remove(pos_1, _length); // delete tmp1 from genome
    remove(0, pos_2);       // delete tmp2 from genome

    // Update promoter list
    if (_length >= PROM_SIZE)
    {
      _gen_unit->move_all_promoters_after(0, -pos_2);

      _gen_unit->look_for_new_promoters_around(0);
    }
  }

  return true;
}

bool ae_dna::do_translocation(int32_t pos_1, int32_t pos_2, int32_t pos_3, int32_t pos_4, bool invert)
{
  // Provided that OriC must be at position 0
  //
  //    1) Note that in Case 1 (without inversion), whichever position comes first, translocating segment [pos_1->pos_2]
  //        to pos_4 through pos_3 is always equivalent to rearrange the sequences from an ABCDE order to ADCBE
  //    2) In Case 2, depending on which position comes first, we may have the following rearrangements :
  //        (ABCDE => ADB'C'E) or (ABCDE => AC'D'BE) where X' stands for "inverted X"
  //
  //  Case 1 : Without inversion
  //
  //         A      B        C       D       E                        A      D        C       B        E
  //      |----->=======[>=======>-------[>-------|        =>      |----->-------[>=======>=======[>-------|
  //          pos_1   pos_3    pos_2   pos_4
  //
  //         A      B        C       D       E                        A      D        C       B        E
  //      |=====>-------[>------->=======[>=======|        =>      |=====>=======[>------->-------[>=======|
  //          pos_2   pos_4    pos_1   pos_3
  //
  //         A      B        C       D       E                        A       D       C        B        E
  //      |====[>========>-------[>------->=======|        =>      |=====[>------->-------[>=======[>=======|
  //          pos_3    pos_2    pos_4   pos_1
  //
  //         A      B        C       D       E                        A       D       C        B        E
  //      |----[>-------->=======[>=======>-------|        =>      |-----[>=======>=======[>-------[>-------|
  //          pos_4    pos_1    pos_3   pos_2
  //
  //
  //  Case 2 : With inversion
  //
  //    Case 2.A
  //
  //         A      B        C       D        E                       A      D        B'      C'       E
  //      |----->=======[>=======>-------<]-------|        =>      |----->-------<]=======<=======<]-------|
  //          pos_1   pos_3    pos_2   pos_4
  //
  //         A      B        C       D       E                        A      D        B'      C'       E
  //      |=====>-------[>------->=======<]=======|        =>      |=====>=======<]-------<-------<]=======|
  //          pos_2   pos_4    pos_1   pos_3
  //
  //
  //    Case 2.B
  //
  //         A      B        C       D       E                        A       C'      D'       B       E
  //      |====[>========>-------<]------->=======|        =>      |=====[>-------<-------[>=======>=======|
  //          pos_3    pos_2    pos_4   pos_1
  //
  //         A      B        C       D       E                        A       C'      D'       B       E
  //      |----<]-------->=======[>=======>-------|        =>      |-----<]=======>=======<]------->-------|
  //          pos_4    pos_1    pos_3   pos_2
  //

  // Determine which position comes first and do the corresponding rearrangement
  int32_t pos_min = ae_utils::min(pos_1, ae_utils::min(pos_2, ae_utils::min(pos_3, pos_4)));

  if (! invert)
  {
    if (pos_min == pos_1)
    {
      ABCDE_to_ADCBE(pos_1, pos_3, pos_2, pos_4);
    }
    else if (pos_min == pos_2)
    {
      ABCDE_to_ADCBE(pos_2, pos_4, pos_1, pos_3);
    }
    else if (pos_min == pos_3)
    {
      ABCDE_to_ADCBE(pos_3, pos_2, pos_4, pos_1);
    }
    else // if (pos_min == pos_4)
    {
      ABCDE_to_ADCBE(pos_4, pos_1, pos_3, pos_2);
    }
  }
  else // invert
  {
    if (pos_min == pos_1)
    {
      ABCDE_to_ADBpCpE(pos_1, pos_3, pos_2, pos_4);
    }
    else if (pos_min == pos_2)
    {
      ABCDE_to_ADBpCpE(pos_2, pos_4, pos_1, pos_3);
    }
    else if (pos_min == pos_3)
    {
      ABCDE_to_ACpDpBE(pos_3, pos_2, pos_4, pos_1);
    }
    else // if (pos_min == pos_4)
    {
      ABCDE_to_ACpDpBE(pos_4, pos_1, pos_3, pos_2);
    }
  }

  return true;
}

bool ae_dna::do_inter_GU_translocation(int32_t pos_1_rel, int32_t pos_2_rel, int32_t pos_3_rel, int32_t pos_4_rel, bool invert)
{
  // TODO check GU lengths according to positions and size limit
  int32_t segment_length = ae_utils::mod(pos_2_rel - pos_1_rel, _length);

  if (pos_1_rel == pos_2_rel) // TODO : should'nt that raise an error?
  {
    return false;
  }

  // Do not allow translocation if it would decrease the size of the origin GU below a threshold
  if ((_length - segment_length) < _gen_unit->get_min_gu_length())
  {
    if (_exp_m->get_output_m()->is_logged(LOG_BARRIER))
    {
      // Write an entry in the barrier log file
      fprintf(_exp_m->get_output_m()->get_log(LOG_BARRIER),
              "%" PRId64 " %" PRId32 " TRANS %" PRId32 " %" PRId32 " %" PRId32
                  " %" PRId32 "\n",
              Time::get_time(), _indiv->get_id(), segment_length, 0, _length,
              _indiv->get_amount_of_dna());
    }
    return false;
  }

  //
  const GeneticUnit& chromosome = _indiv->get_genetic_unit(0);
  const GeneticUnit& plasmid    = _indiv->get_genetic_unit(1);
  // TODO vld (2015-02-23): check if this == is sound
  const GeneticUnit& destination_GU = (_gen_unit == &chromosome) ?
                                      plasmid :
                                      chromosome;

  int32_t dest_gu_size_before = destination_GU.get_seq_length();

  // Do not allow translocation if it would increase the size of the receiving GU above a threshold
  if (dest_gu_size_before + segment_length >
      destination_GU.get_max_gu_length())
  {
    if (_exp_m->get_output_m()->is_logged(LOG_BARRIER))
    {
      // Write an entry in the barrier log file
      fprintf(_exp_m->get_output_m()->get_log(LOG_BARRIER),
              "%" PRId64 " %" PRId32 " TRANS %" PRId32 " %" PRId32 " %" PRId32
                  " %" PRId32 "\n",
              Time::get_time(), _indiv->get_id(), segment_length, 0,
              dest_gu_size_before, _indiv->get_amount_of_dna());
    }
    return false;
  }

  //~ // Provided that OriC must be at position 0
  //~ //
  //~ //      Gen unit 1           Gen unit 2              Gen unit 1       Gen unit 2
  //~ //
  //~ //   Case 1: inter_GU_ABCDE_to_ACDBE
  //~ //
  //~ //         A    B     C        D     E                  A    C          D     B     E
  //~ //      |---->=====>----|    |----[>----|             |---->----|     |----[>====[>----|
  //~ //          p1r   p2r            p4r                      p1r             p4r   p4r+(p2r-p1r)
  //~ //
  //~ //
  //~ //   Case 2: inter_GU_ABCDE_to_BDCAE
  //~ //
  //~ //         A    B     C        D     E                   B               D     C     A   E
  //~ //      |====>----->====|    |----[>----|             |-----|         |----[>====[>====>----|
  //~ //          p2r   p1r             p4r                                     p4r
  //~ //                                                                            p4r+(_length-p1r)
  //~ //                                                                                  p4r+(_length-(p1r-p2r))

  // Determine which position comes first and do the corresponding rearrangement
  // int32_t pos_min = ae_utils::min(pos_1, pos_2);

  if (! invert)
  {
    if (pos_1_rel < pos_2_rel)
    {
      // if (_gen_unit == _indiv->get_genetic_unit(0) /*chromosome*/)
      //      {
      //           printf("do_inter_GU_translocation(%"PRId32", %"PRId32", %"PRId32", %"PRId32", %s) (sizes : %"PRId32" %"PRId32" %"PRId32")\n",
      //                   pos_1_rel, pos_2_rel, pos_3_rel, pos_4_rel + _length, invert?"invert":"plain",
      //                   _indiv->get_genetic_unit(0)->get_dna()->get_length(),
      //                   _indiv->get_genetic_unit(1)->get_dna()->get_length(),
      //                   _indiv->get_amount_of_dna());
      //       }
      //       else
      //       {
      //           printf("do_inter_GU_translocation(%"PRId32", %"PRId32", %"PRId32", %"PRId32", %s) (sizes : %"PRId32" %"PRId32" %"PRId32")\n",
      //                   pos_1_rel + _indiv->get_genetic_unit(0)->get_dna()->get_length(),
      //                   pos_2_rel + _indiv->get_genetic_unit(0)->get_dna()->get_length(),
      //                   pos_3_rel + _indiv->get_genetic_unit(0)->get_dna()->get_length(),
      //                   pos_4_rel, invert?"invert":"plain",
      //                   _indiv->get_genetic_unit(0)->get_dna()->get_length(),
      //                   _indiv->get_genetic_unit(1)->get_dna()->get_length(),
      //                   _indiv->get_amount_of_dna());
      //       }

      segment_length = ae_utils::mod(pos_2_rel - pos_1_rel, _length);
      inter_GU_ABCDE_to_ACDBE(pos_1_rel, pos_2_rel, pos_4_rel);
    }
    else
    {
      //       if (_gen_unit == _indiv->get_genetic_unit(0) /*chromosome*/)
      //       {
      //           printf("do_inter_GU_translocation(%"PRId32", %"PRId32", %"PRId32", %"PRId32", %s) (sizes : %"PRId32" %"PRId32" %"PRId32")\n",
      //                   pos_1_rel, pos_2_rel, pos_3_rel, pos_4_rel + _length, invert?"invert":"plain",
      //                   _indiv->get_genetic_unit(0)->get_dna()->get_length(),
      //                   _indiv->get_genetic_unit(1)->get_dna()->get_length(),
      //                   _indiv->get_amount_of_dna());
      //       }
      //       else
      //       {
      //           printf("do_inter_GU_translocation(%"PRId32", %"PRId32", %"PRId32", %"PRId32", %s) (sizes : %"PRId32" %"PRId32" %"PRId32")\n",
      //                   pos_1_rel + _indiv->get_genetic_unit(0)->get_dna()->get_length(),
      //                   pos_2_rel + _indiv->get_genetic_unit(0)->get_dna()->get_length(),
      //                   pos_3_rel + _indiv->get_genetic_unit(0)->get_dna()->get_length(),
      //                   pos_4_rel, invert?"invert":"plain",
      //                   _indiv->get_genetic_unit(0)->get_dna()->get_length(),
      //                   _indiv->get_genetic_unit(1)->get_dna()->get_length(),
      //                   _indiv->get_amount_of_dna());
      //       }

      segment_length = ae_utils::mod(pos_1_rel - pos_2_rel, _length);
      inter_GU_ABCDE_to_BDCAE(pos_2_rel, pos_1_rel, pos_4_rel);
    }
  }
  else // invert
  {
    if (pos_1_rel < pos_2_rel)
    {
      segment_length = ae_utils::mod(pos_2_rel - pos_1_rel, _length);
      do_inversion(pos_1_rel, pos_2_rel);
      inter_GU_ABCDE_to_ACDBE(pos_1_rel, pos_2_rel, pos_4_rel);
    }
    else // pos_1_rel > pos_2_rel
    {
      segment_length = ae_utils::mod(pos_1_rel - pos_2_rel, _length);
      if (pos_2_rel != 0)       { do_inversion(0, pos_2_rel); }
      if (pos_1_rel != _length) { do_inversion(pos_1_rel, _length); }
      inter_GU_ABCDE_to_BDCAE(pos_2_rel, pos_1_rel, pos_4_rel);
    }
  }

  return true;
}

bool ae_dna::do_inversion(int32_t pos_1, int32_t pos_2)
// Invert segment going from pos_1 (included) to pos_2 (excluded)
// Exemple : sequence 011101001100 => 110011010001
{
  if (pos_1 == pos_2) return false; // Invert everything <=> Invert nothing!
  assert(pos_1 < pos_2);

  //
  //       pos_1         pos_2                   -> 0-
  //         |             |                   -       -
  // 0--------------------------------->      -         -
  //         ===============                  -         - pos_1
  //           tmp (copy)                      -       -
  //                                             -----      |
  //                                             pos_2    <-'
  //

  int32_t seg_length = pos_2 - pos_1;

  // Create the inverted sequence
  char* inverted_segment = NULL;
  inverted_segment = new char[seg_length+1];
  for (int32_t i = 0, j = pos_2 - 1 ; i < seg_length ; i++, j--)
  {
    if (_data[j] == '0') inverted_segment[i] = '1';
    else                   inverted_segment[i] = '0';
  }
  inverted_segment[seg_length] = '\0';

  // Remove promoters that included a breakpoint
  _gen_unit->remove_promoters_around(pos_1);
  _gen_unit->remove_promoters_around(pos_2);

  // Invert the sequence
  replace(pos_1, inverted_segment, seg_length);

  // Update promoter list
  if (_length >= PROM_SIZE)
  {
    _gen_unit->invert_promoters_included_in(pos_1, pos_2);

    _gen_unit->look_for_new_promoters_around(pos_1);
    _gen_unit->look_for_new_promoters_around(pos_2);
  }

  delete [] inverted_segment;

  return true;
}

bool ae_dna::do_insertion(int32_t pos, const char* seq_to_insert, int32_t seq_length)
{
  // Remove the promoters that will be broken
  _gen_unit->remove_promoters_around(pos);

  // Insert the sequence
  insert(pos, seq_to_insert, seq_length);

  // Look for new promoters
  if (_length >= PROM_SIZE)
  {
    _gen_unit->move_all_promoters_after(pos, seq_length);
    _gen_unit->look_for_new_promoters_around(pos, pos + seq_length);
  }

  return true;
}


ae_mutation* ae_dna::do_ins_HT(int32_t parent_id)
{
  ae_mutation* mut = NULL;

  // TODO <david.parsons@inria.fr> disabled
//  int32_t nb_indivs = _exp_m->get_pop()->get_nb_indivs();
//
//  // Insertion transfer
//  // Requirements:
//  //    * A circular exogenote => an alignment on the donor chromosome
//  //    * An alignment between the exogenote and the endogenote
//
//  // 1) Draw a random donor (uniform drawing).
//  // We use the rank because indivs are sorted by rank (1 for the worst, POP_SIZE for the best).
//  Individual * donor = NULL;
//  do donor = _exp_m->get_pop()->get_indiv_by_rank(_exp_m->get_sel()->get_prng()->random(nb_indivs) + 1);
//  while (donor->get_id() == parent_id);
//
//  // 2) Look for an alignment within the donor genome
//  ae_vis_a_vis* alignment_1   = NULL;
//  ae_dna*       donor_dna     = donor->get_genetic_unit(0).get_dna();
//  int32_t       nb_pairs_1    = (int32_t)(ceil(donor_dna->get_length() * _indiv->get_neighbourhood_rate()));
//
//  alignment_1 = donor_dna->search_alignment(donor_dna, nb_pairs_1, DIRECT);
//
//    if (alignment_1 != NULL)
//    {
//      // 3) Make a copy of the sequence to be transferred (the exogenote)
//      GeneticUnit* exogenote = donor_dna->copy_into_new_GU(alignment_1->get_i_1(), alignment_1->get_i_2());
//
//      // 4) Look for an alignments between the exogenote and the endogenote
//      ae_vis_a_vis* alignment_2 = NULL;
//      int32_t       nb_pairs_2  = (int32_t)(ceil(get_length() * _indiv->get_neighbourhood_rate()));
//
//      alignment_2 = exogenote->get_dna()->search_alignment(this, nb_pairs_2, BOTH_SENSES);
//
//      if (alignment_2 != NULL)
//      {
//        int32_t gu_length_before  = _length;
//        int32_t gu_length_after   = gu_length_before + exogenote->get_dna()->get_length();
//        int32_t genome_length_before = _indiv->get_amount_of_dna();
//        int32_t genome_length_after = genome_length_before + exogenote->get_dna()->get_length();
//
//        if ((genome_length_after > _indiv->get_max_genome_length()) ||
//            (gu_length_after > _gen_unit->get_max_gu_length()))
//        {
//          if (_exp_m->get_output_m()->is_logged(LOG_BARRIER))
//          {
//            // Write an entry in the barrier log file
//            fprintf(_exp_m->get_output_m()->get_log(LOG_BARRIER),
//                "%" PRId64 " %" PRId32 " INS_TRANSFER %" PRId32 " %" PRId32
//                " %" PRId32 " %" PRId32 "\n",
//                Time::get_time(),
//                _indiv->get_id(),
//                exogenote->get_dna()->get_length(),
//                0,
//                gu_length_before,
//                genome_length_before);
//          }
//        }
//        else
//        {
//          insert_GU(exogenote,
//              alignment_2->get_i_2(), alignment_2->get_i_1(),
//              (alignment_2->get_sense() == INDIRECT));
//          //~ fprintf(logfile, "RESULT:\n%s\n\n\n", new_indiv_dna->get_data());
//          //~ fflush(logfile);
//
//          // Write a line in transfer logfile
//          if (_exp_m->get_output_m()->is_logged(LOG_TRANSFER))
//          {
//            fprintf(_exp_m->get_output_m()->get_log(LOG_TRANSFER),
//                    "%" PRId64 " %" PRId32 " %" PRId32 " 0 %" PRId32
//                    " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32
//                    " %" PRId32 " %" PRId16 " %" PRId32 " %" PRId32
//                    " %" PRId16 "\n",
//                    Time::get_time(),
//                    _indiv->get_id(),
//                    donor->get_id(),
//                    exogenote->get_dna()->get_length(),
//                    0,
//                    genome_length_before,
//                    get_length(),
//                    alignment_1->get_i_1(),
//                    alignment_1->get_i_2(),
//                    alignment_1->get_score(),
//                    alignment_2->get_i_1(),
//                    alignment_2->get_i_2(),
//                    alignment_2->get_score());
//          }
//
//          #ifdef BIG_DEBUG
//            ae_common::sim->get_logs()->flush();
//            _indiv->assert_promoters();
//            _indiv->assert_promoters_order();
//          #endif
//
//          if (_exp_m->get_output_m()->get_record_tree() && _exp_m->get_output_m()->get_tree_mode() == NORMAL)
//          {
//            char* donor_seq;
//            if(alignment_2->get_sense() == DIRECT)
//            {
//              donor_seq = exogenote->get_dna()->get_subsequence(alignment_2->get_i_1(),alignment_2->get_i_1(), LEADING);
//            }
//            else
//            {
//              donor_seq = exogenote->get_dna()->get_subsequence(alignment_2->get_i_1(),alignment_2->get_i_1(), LAGGING);
//            }
//            // Report the transfer
//            mut = new ae_mutation();
//            mut->report_ins_HT(alignment_1->get_i_1(), alignment_1->get_i_2(), alignment_2->get_i_1(), alignment_2->get_i_2(), exogenote->get_dna()->get_length(), alignment_1->get_score(),alignment_2->get_score(),  donor->get_id(), alignment_2->get_sense(), donor_seq);
//            delete [] donor_seq;
//          }
//        }
//
//        delete alignment_2;
//      }
//
//      delete exogenote;
//      delete alignment_1;
//    }

  return mut;
}

ae_mutation* ae_dna::do_repl_HT(int32_t parent_id)
{
  ae_mutation* mut = NULL;

  // TODO <david.parsons@inria.fr> disabled
//  int32_t nb_indivs = _exp_m->get_pop()->get_nb_indivs();
//
//  // Replacement transfer
//  // Requirements:
//  //    * 2 distinct alignments between the (linear) exogenote and the endogenote
//
//  //printf("Indiv: %d\n",new_indiv->get_id());
//  // 1) Draw a random donor (uniform drawing).
//  // We use the rank because indivs are sorted by rank (1 for the worst, POP_SIZE for the best).
//  Individual * donor = NULL;
//  do donor = _exp_m->get_pop()->get_indiv_by_rank(_exp_m->get_sel()->get_prng()->random(nb_indivs) + 1);
//  while (donor->get_id() == parent_id);
//
//  // 2) Look for an alignment between the parent genome and the donor genome
//  ae_vis_a_vis* alignment_1   = NULL;
//  ae_vis_a_vis* alignment_2   = NULL;
//  ae_dna*       donor_dna     = donor->get_genetic_unit(0).get_dna();
//  ae_sense      sense         = (_exp_m->get_sel()->get_prng()->random() < 0.5) ? DIRECT : INDIRECT;
//  int32_t       nb_pairs_1    = (int32_t)(ceil(get_length() * _indiv->get_neighbourhood_rate()));
//  int32_t       nb_pairs_2    = (int32_t)(ceil(get_length() * _indiv->get_neighbourhood_rate()));
//  int8_t        research_sense = 0;
//
//  alignment_1 = search_alignment(donor_dna, nb_pairs_1, sense);
//  if (alignment_1 != NULL)
//  {
//    //printf("\talignement 1: %d %d\n", alignment_1->get_i_1(), alignment_1->get_i_2());
//    if(_exp_m->get_repl_HT_with_close_points())
//    {
//      //printf("Id: %d\n",new_indiv->get_id());
//
//      alignment_2 = search_alignment_around_positions(donor_dna, alignment_1->get_i_1(), alignment_1->get_i_2(), alignment_1->get_sense(), research_sense);
//      if (alignment_2 != NULL && alignment_2->get_i_1() == alignment_1->get_i_1())
//      {
//        delete alignment_2;
//        alignment_2 = NULL;
//      }
//      if (alignment_2 != NULL)
//      {
//        // If the second alignment is found upstream of the first alignment, they are inverse to facilitate
//        if(research_sense == -1)
//        {
//          ae_vis_a_vis* tmp_alignment = new ae_vis_a_vis(*alignment_1);
//          alignment_1->copy(alignment_2);
//          alignment_2->copy(tmp_alignment);
//          delete tmp_alignment;
//        }
//      }
//    }
//    else
//    {
//      // Look for a second alignement between the parent and the donor (must be different from alignment_1)
//      while (alignment_2 == NULL && nb_pairs_2 > 0)
//      {
//        alignment_2 = search_alignment(donor_dna, nb_pairs_2, sense);
//
//        // Forbid the replacement of the whole genome of the parent
//        if (alignment_2 != NULL && alignment_2->get_i_1() == alignment_1->get_i_1())
//        {
//          delete alignment_2;
//          alignment_2 = NULL;
//        }
//      }
//    }
//
//    // If both alignments were found, proceed to the transfer
//    if (alignment_2 != NULL)
//    {
//      int32_t gu_length_before  = _length;
//      int32_t exogenote_length      = ae_utils::mod(alignment_2->get_i_2() - alignment_1->get_i_2() - 1, donor_dna->get_length()) + 1;
//      int32_t replaced_seq_length   = ae_utils::mod(alignment_2->get_i_1() - alignment_1->get_i_1() - 1, gu_length_before) + 1;
//      int32_t gu_length_after   = gu_length_before - replaced_seq_length + exogenote_length;
//
//      int32_t genome_length_before = _indiv->get_amount_of_dna();
//      int32_t genome_length_after = genome_length_before - replaced_seq_length + exogenote_length;
//
//      if (genome_length_after < _indiv->get_min_genome_length() ||
//          genome_length_after > _indiv->get_max_genome_length() ||
//          gu_length_after < _gen_unit->get_min_gu_length() ||
//          gu_length_after > _gen_unit->get_max_gu_length())
//      {
//        if (_exp_m->get_output_m()->is_logged(LOG_BARRIER))
//        {
//          // Write an entry in the barrier log file
//          fprintf(_exp_m->get_output_m()->get_log(LOG_BARRIER),
//              "%" PRId64 " %" PRId32 " REPL_TRANSFER %" PRId32 " %" PRId32
//              " %" PRId32 " %" PRId32 "\n",
//              Time::get_time(),
//              _indiv->get_id(),
//              exogenote_length,
//              replaced_seq_length,
//              gu_length_before,
//              genome_length_before);
//        }
//      }
//      else
//      {
//        // 3) Make a copy of the sequence to be transferred (the exogenote)
//        GeneticUnit* exogenote = NULL;
//        if (sense == DIRECT)
//        {
//          exogenote = donor_dna->copy_into_new_GU(alignment_1->get_i_2(), alignment_2->get_i_2());
//        }
//        else
//        {
//          exogenote = donor_dna->copy_into_new_GU(alignment_2->get_i_2(), alignment_1->get_i_2());
//        }
//
//        char* alignment1_parent_dna = nullptr;
//        char* alignment2_parent_dna = nullptr;
//        char* alignment1_donor_dna  = nullptr;
//        char* alignment2_donor_dna  = nullptr;
//        if (_exp_m->get_output_m()->is_logged(LOG_TRANSFER) == true)
//        {
//          if(sense  == DIRECT)
//          {
//            alignment1_parent_dna = get_subsequence(alignment_1->get_i_1(), alignment_1->get_i_1() + 2*_indiv->get_align_w_zone_h_len() + _indiv->get_align_max_shift(), LEADING);
//            alignment2_parent_dna = get_subsequence(alignment_2->get_i_1(), alignment_2->get_i_1() + 2*_indiv->get_align_w_zone_h_len() + _indiv->get_align_max_shift(), LEADING);
//            alignment1_donor_dna = donor_dna->get_subsequence(alignment_1->get_i_2(), alignment_1->get_i_2() + 2*_indiv->get_align_w_zone_h_len() + _indiv->get_align_max_shift(), LEADING);
//            alignment2_donor_dna = donor_dna->get_subsequence(alignment_2->get_i_2(), alignment_2->get_i_2() + 2*_indiv->get_align_w_zone_h_len() + _indiv->get_align_max_shift(), LEADING);
//          }
//          else
//          {
//            alignment1_parent_dna = get_subsequence(alignment_1->get_i_1(), alignment_1->get_i_1() + 2*_indiv->get_align_w_zone_h_len() + _indiv->get_align_max_shift(), LEADING);
//            alignment2_parent_dna = get_subsequence(alignment_2->get_i_1(), alignment_2->get_i_1() + 2*_indiv->get_align_w_zone_h_len() + _indiv->get_align_max_shift(), LEADING);
//            alignment1_donor_dna = donor_dna->get_subsequence(alignment_1->get_i_2(), alignment_1->get_i_2() - 2*_indiv->get_align_w_zone_h_len() - _indiv->get_align_max_shift(), LAGGING);
//            alignment2_donor_dna = donor_dna->get_subsequence(alignment_2->get_i_2(), alignment_2->get_i_2() - 2*_indiv->get_align_w_zone_h_len() - _indiv->get_align_max_shift(), LAGGING);
//          }
//        }
//
//        // Delete the sequence to be replaced
//        do_deletion(alignment_1->get_i_1(), alignment_2->get_i_1());
//        if (alignment_1->get_i_1() < alignment_2->get_i_1())
//        {
//          insert_GU(exogenote, alignment_1->get_i_1(), 0, sense == INDIRECT);
//        }
//        else
//        {
//          insert_GU(exogenote, 0, 0, sense == INDIRECT);
//        }
//
//        // Write a line in transfer logfile
//        if (_exp_m->get_output_m()->is_logged(LOG_TRANSFER))
//        {
//          fprintf(_exp_m->get_output_m()->get_log(LOG_TRANSFER),
//              "%" PRId64 " %" PRId32 " %" PRId32 " 1 %" PRId32 " %" PRId32
//              " %" PRId32 " %" PRId32 " %" PRId16 " %" PRId32 " %" PRId32
//              " %" PRId16 " %" PRId32 " %" PRId32 " %" PRId16 " %" PRId16 "\n",
//              Time::get_time(),
//              _indiv->get_id(),
//              donor->get_id(),
//              exogenote->get_dna()->get_length(),
//              replaced_seq_length,
//              genome_length_before,
//              get_length(),
//              (int16_t) alignment_1->get_sense(),
//              alignment_1->get_i_1(),
//              alignment_1->get_i_2(),
//              alignment_1->get_score(),
//              alignment_2->get_i_1(),
//              alignment_2->get_i_2(),
//              alignment_2->get_score() ,
//              (int16_t) research_sense);
//
//        fprintf(_exp_m->get_output_m()->get_log(LOG_TRANSFER),
//            "\tAlignment 1:\n\t\t%s\n\t\t%s\n\tAlignment 2:\n\t\t%s\n\t\t%s\n",
//            alignment1_parent_dna, alignment1_donor_dna,
//            alignment2_parent_dna, alignment2_donor_dna);
//
//            delete [] alignment1_parent_dna;
//            delete [] alignment2_parent_dna;
//            delete [] alignment1_donor_dna;
//            delete [] alignment2_donor_dna;
//        }
//
//        if (_exp_m->get_output_m()->get_record_tree() && _exp_m->get_output_m()->get_tree_mode() == NORMAL)
//        {
//          // Report the transfer
//          char* donor_seq;
//          if(alignment_2->get_sense() == DIRECT)
//          {
//            donor_seq = exogenote->get_dna()->get_subsequence(0,exogenote->get_dna()->get_length(), LEADING);
//          }
//          else
//          {
//            donor_seq = exogenote->get_dna()->get_subsequence(0,exogenote->get_dna()->get_length(), LAGGING);
//          }
//          mut = new ae_mutation();
//          mut->report_repl_HT(alignment_1->get_i_1(), alignment_1->get_i_2(), alignment_2->get_i_1(), alignment_2->get_i_2(), replaced_seq_length, exogenote->get_dna()->get_length(), alignment_1->get_score(),alignment_2->get_score(),  donor->get_id(), alignment_2->get_sense(),donor_seq);
//          delete [] donor_seq;
//        }
//
//        delete exogenote;
//      }
//      delete alignment_2;
//    }
//    delete alignment_1;
//  }
  return mut;
}

bool ae_dna::do_ins_HT(int32_t pos, const char* seq_to_insert, int32_t seq_length)
{
  // Remove the promoters that will be broken
  _gen_unit->remove_promoters_around(pos);

  // Insert the sequence
  insert(pos, seq_to_insert, seq_length);

  // Look for new promoters
  if (_length >= PROM_SIZE)
  {
    _gen_unit->move_all_promoters_after(pos, seq_length);
    _gen_unit->look_for_new_promoters_around(pos, pos + seq_length);
  }

  return true;
}

bool ae_dna::do_repl_HT(int32_t pos1, int32_t pos2, const char* seq_to_insert, int32_t seq_length)
{
  // Remove the promoters that will be broken
  _gen_unit->remove_promoters_around(pos1);

  // Delete the replaced segment
  do_deletion(pos1, pos2);

  // Insert the sequence
  int32_t insertion_position;
  if(pos1 < pos2)
  {
    insertion_position = pos1;
  }
  else
  {
    insertion_position = 0;
  }
  insert(insertion_position, seq_to_insert, seq_length);

  // Look for new promoters
  if (_length >= PROM_SIZE)
  {
    _gen_unit->move_all_promoters_after(insertion_position, seq_length);
    _gen_unit->look_for_new_promoters_around(insertion_position, insertion_position + seq_length);
  }

  return true;
}

void ae_dna::undergo_this_mutation(const ae_mutation * mut)
{
  if(mut == NULL) return;

  int32_t pos1, pos2, pos3, pos4;
  int32_t length;
  bool invert;
  char *seq = NULL;
  ae_sense sense;

  switch(mut->get_mut_type())
  {
    case SWITCH:
      mut->get_infos_point_mutation(&pos1);
      do_switch(pos1);
      break;
    case S_INS:
      mut->get_infos_small_insertion(&pos1, &length);
      seq = new char[length + 1];
      mut->get_sequence_small_insertion(seq);
      do_small_insertion(pos1, length, seq);
      delete [] seq;
      break;
    case S_DEL:
      mut->get_infos_small_deletion(&pos1, &length);
      do_small_deletion(pos1, length);
      break;
    case DUPL:
      mut->get_infos_duplication(&pos1, &pos2, &pos3);
      do_duplication(pos1, pos2, pos3);
      break;
    case DEL:
      mut->get_infos_deletion(&pos1, &pos2);
      do_deletion(pos1, pos2);
      break;
    case TRANS:
      mut->get_infos_translocation(&pos1, &pos2, &pos3, &pos4, &invert);
      if (_indiv->get_with_alignments())
      {
        // Extract the segment to be translocated
        GeneticUnit* translocated_segment = extract_into_new_GU(pos1, pos2);

        // Reinsert the segment
        insert_GU(translocated_segment, pos3, pos4, invert);
      }
      else
      {
        do_translocation(pos1, pos2, pos3, pos4, invert);
      }
      break;
    case INV:
      mut->get_infos_inversion(&pos1, &pos2);
      do_inversion(pos1, pos2);
      break;
    case INS_HT:
      mut->get_infos_ins_HT(&pos1, &pos2, &pos3, &pos4, &sense, &length);
      seq = new char[length + 1];
      mut->get_sequence_ins_HT(seq);
      do_ins_HT(pos4, seq, length);
      delete [] seq;
      break;
    case REPL_HT:
      mut->get_infos_repl_HT(&pos1, &pos2, &pos3, &pos4, &sense, &length);
      seq = new char[length + 1];
      mut->get_sequence_repl_HT(seq);
      do_repl_HT(pos1, pos3, seq, length);
      delete [] seq;
      break;
    default :
      fprintf(stderr, "ERROR, invalid mutation type in file %s:%d\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      break;
  }
}

void ae_dna::compute_statistical_data(void)
// *************************************************************************************************
//                                      WARNING     Deprecated
// *************************************************************************************************
// Values that are computed :
//    Number of each type of mutations and rearrangements undergone by the individual
//    Number of small mutations
//    Number of rearrangements
{
  assert(false);
  //~ ae_list_node<ae_mutation*>* mut_node  = _replic_report->_mutations->get_first();
  //~ ae_mutation*  mut;

  //~ while (mut_node != NULL)
  //~ {
  //~ mut = mut_node->get_obj();

  //~ switch(mut->get_mut_type())
  //~ {
  //~ case SWITCH :
  //~ _replic_report->_nb_switch++;
  //~ break;
  //~ case S_INS :
  //~ _replic_report->_nb_small_insertions++;
  //~ break;
  //~ case S_DEL :
  //~ _replic_report->_nb_small_deletions++;
  //~ break;
  //~ case DUPL :
  //~ _replic_report->_nb_duplications++;
  //~ break;
  //~ case DEL :
  //~ _replic_report->_nb_deletions++;
  //~ break;
  //~ case TRANS :
  //~ _replic_report->_nb_translocations++;
  //~ // TODO : if inter GU ...
  //~ break;
  //~ case INV :
  //~ _replic_report->_nb_inversions++;
  //~ break;
  //~ }

  //~ mut_node = mut_node->get_next();
  //~ }
}

void ae_dna::set_GU(std::vector<std::list<ae_rna>> rna_list, const GeneticUnit* GU) {
  for (auto& strand: {LEADING, LAGGING})
    for (auto& rna: rna_list[strand])
      rna.set_genetic_unit(GU);
}

GeneticUnit* ae_dna::extract_into_new_GU(int32_t pos_1, int32_t pos_2)
{
  assert(pos_1 < pos_2);
  int32_t seq_length = pos_2 - pos_1;

  // ==================== Remove/Extract promoters from old sequence ====================
  // Remove promoters around breakpoints
  _gen_unit->remove_promoters_around(pos_1);
  _gen_unit->remove_promoters_around(pos_2);

  // Remove promoters belonging to the sequence (to be extracted) from the "old" GU
  // and put them in a stand-alone promoter list (with indices ranging from 0 to seq_length-1)
  Promoters2Strands proms_GU_1 = {{},{}};
  _gen_unit->extract_promoters_included_in(pos_1, pos_2, proms_GU_1);
  GeneticUnit::shift_promoters(proms_GU_1, -pos_1, _length);

  // ==================== Manage sequences ====================
  // Copy the sequence in a stand-alone char* (size must be multiple of BLOCK_SIZE)
  int32_t length_GU_1     = seq_length;
  char*   sequence_GU_1   = new char[nb_blocks(length_GU_1) * BLOCK_SIZE];
  memcpy(sequence_GU_1, &_data[pos_1], length_GU_1 * sizeof(char));
  sequence_GU_1[length_GU_1] = '\0';

  // Remove the sequence from the "old" GU
  int32_t  length_GU_0     = _length - seq_length;
  char*    sequence_GU_0   = new char[nb_blocks(length_GU_0) * BLOCK_SIZE];
  memcpy(sequence_GU_0,          _data,          pos_1 * sizeof(char));
  memcpy(&sequence_GU_0[pos_1],  &_data[pos_2],  (_length - pos_2) * sizeof(char));
  sequence_GU_0[length_GU_0] = '\0';

  set_data(sequence_GU_0, length_GU_0);

  // ==================== Create the new genetic unit ====================
  GeneticUnit* GU_1 = new GeneticUnit(_indiv, sequence_GU_1, length_GU_1, proms_GU_1);

  // ==================== Update promoter lists ====================
  // Shift the position of the promoters of the "old" GU
  _gen_unit->move_all_promoters_after(pos_1, -seq_length);

  // Look for new promoters around breakpoints
  _gen_unit->look_for_new_promoters_around(pos_1);
  GU_1->look_for_new_promoters_around(0);

  return GU_1;
}


/*!
  \brief Copy the sequence going from pos_1 (included) to pos_2 (excluded) into a new standalone genetic unit.

  The new genetic unit's list of promoter is up-to-date.
  if (pos_1 == pos_2), the whole genome is copied
*/
GeneticUnit* ae_dna::copy_into_new_GU(int32_t pos_1, int32_t pos_2) const
{
  int32_t seq_length = ae_utils::mod(pos_2 - pos_1, _length);
  if (seq_length == 0) seq_length = _length;

  // ==================== Copy promoters from old sequence ====================
  // Copy the promoters belonging to the sequence to be copied from the "old" GU
  // into a stand-alone promoter list (with indices ranging from 0 to seq_length-1)
  Promoters2Strands proms_new_GU = {{},{}};
  _gen_unit->copy_promoters_included_in(pos_1, pos_2, proms_new_GU);
  GeneticUnit::shift_promoters(proms_new_GU, -pos_1, _length);


  // ==================== Manage sequences ====================
  // Copy the sequence in a stand-alone char* (size must be multiple of BLOCK_SIZE)
  int32_t length_new_GU     = seq_length;
  char*   sequence_new_GU   = new char[nb_blocks(length_new_GU) * BLOCK_SIZE];
  if (pos_1 < pos_2)
  {
    memcpy(sequence_new_GU, &_data[pos_1], length_new_GU * sizeof(char));
  }
  else
  {
    memcpy(sequence_new_GU, &_data[pos_1], (_length - pos_1) * sizeof(char));
    memcpy(&sequence_new_GU[_length - pos_1], &_data[0], pos_2 * sizeof(char));
  }
  sequence_new_GU[length_new_GU] = '\0';


  // ==================== Create the new genetic unit ====================
  GeneticUnit* new_GU = new GeneticUnit(_indiv, sequence_new_GU, length_new_GU, proms_new_GU);

  // ==================== Update new GU promoter list ====================
  // Look for new promoters around breakpoints
  new_GU->look_for_new_promoters_around(0);

  return new_GU;
}

/*!
  \brief Insert the genetic unit GU_to_insert at pos_B, through pos_D. (sequence is inverted if invert == true)

  GU to insert: segments C and D, breakpoint pos_D.
  Current GU:   segments A and B, breakpoint pos_B.

  \verbatim
  If invert is false, the insertion will render ADCB
           A        B                  C        D                          A       D       C        B
       |-------[>-------|     +     |=====[>========|        =>        |-------[>=====|========[>-------|
             pos_B                      pos_D


  If invert is true, the insertion will render ACpDpB, with Cp (resp. Dp) = inverted C (resp. D).
           A        B                  C        D                          A       Cp     Dp        B
       |-------[>-------|     +     |=====<]========|        =>        |-------[>=====|========[>-------|
             pos_B                      pos_D
  \endverbatim

  Sequence from GU_to_insert is untouched but its list of promoters is emptied
*/
void ae_dna::insert_GU(GeneticUnit* GU_to_insert, int32_t pos_B, int32_t pos_D, bool invert)
{
  // Compute segment lengths
  const char* GUti_data = GU_to_insert->get_dna()->get_data();
  int32_t len_A     = pos_B;
  int32_t len_B     = _length - pos_B;
  int32_t len_C     = pos_D;
  int32_t len_D     = GU_to_insert->get_dna()->get_length() - pos_D;
  int32_t len_AB    = _length;
  int32_t len_CD    = GU_to_insert->get_dna()->get_length();
  int32_t len_ABCD  = len_AB + len_CD;

  //~ printf("len_A : %"PRId32", len_B : %"PRId32", len_C : %"PRId32", len_D : %"PRId32",\nlen_AB : %"PRId32", len_CD : %"PRId32", len_ABCD : %"PRId32"\n",
  //~ len_A, len_B, len_C, len_D, len_AB, len_CD, len_ABCD);


  // ==================== Insert the sequence ====================
  // Create new genome
  char* new_seq = new char[ BLOCK_SIZE * nb_blocks(len_ABCD) ];

  // Insert A
  memcpy(new_seq, _data, len_A * sizeof(char));

  // Insert C and D (inverted?)
  if (invert)
  {
    // Build Cp and Dp
    char* seq_Cp = new char[pos_D+1];
    for (int32_t i = 0, j = pos_D - 1 ; i < len_C ; i++, j--)
    {
      if (GUti_data[j] == '0') seq_Cp[i] = '1';
      else                       seq_Cp[i] = '0';
    }
    seq_Cp[len_C] = '\0';

    char* seq_Dp = new char[len_D+1];
    for (int32_t i = 0, j = len_CD - 1 ; i < len_D ; i++, j--)
    {
      if (GUti_data[j] == '0') seq_Dp[i] = '1';
      else                       seq_Dp[i] = '0';
    }
    seq_Dp[len_D] = '\0';

    // Insert Cp and DP
    // TODO : Maybe we should construct Cp and Dp directly at their rightful place...
    memcpy(&new_seq[len_A],        seq_Cp, len_C * sizeof(char));
    memcpy(&new_seq[len_A+len_C],  seq_Dp, len_D * sizeof(char));

    delete [] seq_Cp;
    delete [] seq_Dp;
  }
  else // if (invert == false)
  {
    // Insert D and C
    memcpy(&new_seq[len_A],        &GUti_data[pos_D],  len_D * sizeof(char));
    memcpy(&new_seq[len_A+len_D],  GUti_data,          len_C * sizeof(char));
  }

  // Insert B
  memcpy(&new_seq[len_A+len_C+len_D], &_data[pos_B], len_B * sizeof(char));
  new_seq[len_ABCD] = '\0';


  // Remove promoters that are astride segments A and B : breakpoint pos_B
  _gen_unit->remove_promoters_around(pos_B);


  set_data(new_seq, len_ABCD);




  // ==================== Manage promoters ====================
  // Remove promoters that are astride segments C and D : breakpoint pos_D
  GU_to_insert->remove_promoters_around(pos_D);

  // Shift the position of the promoters of segment B
  _gen_unit->move_all_promoters_after(pos_B, len_CD);

  // Extract promoters of segments C and D.
  // NOTE : We want ALL THE PROMOTERS to be transfered, not only those that are completely included
  //        in segment C or D. Hence, we will use extract_promoters_starting_between() instead of
  //        extract_promoters_included_in().
  // NOTE : Once removed the promoters starting on sequence D, the remaining is precisely the promoters
  //        starting on sequence C (and they are at their rightful position). We can hence directly use
  //        the list of promoters from GU_to_insert.
  Promoters2Strands proms_C = {{},{}};
  Promoters2Strands proms_D = {{},{}};

  if (pos_D != 0) // TODO : Manage this in the different functions? with a parameter WholeGenomeEventHandling ?
  {
    GU_to_insert->extract_promoters_starting_between(0, pos_D, proms_C);
  }
  GU_to_insert->extract_promoters_starting_between(pos_D, len_CD, proms_D);
  assert(GU_to_insert->get_rna_list()[LEADING].empty());
  assert(GU_to_insert->get_rna_list()[LAGGING].empty());
  GeneticUnit::shift_promoters(proms_D, -len_C, len_D);

  if (invert)
  {
    //~ printf("+++++++++++++++++++++++++++++++++++++\n");
    //~ GeneticUnit::print_rnas(proms_C);
    //~ GeneticUnit::print_rnas(proms_D);
    //~ printf("//////////////////////////////////////\n");

    GeneticUnit::invert_promoters(proms_C, len_C);
    GeneticUnit::invert_promoters(proms_D, len_D);

    //~ GeneticUnit::print_rnas(proms_C);
    //~ GeneticUnit::print_rnas(proms_D);
    //~ printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");

    _gen_unit->insert_promoters_at(proms_C, len_A);
    _gen_unit->insert_promoters_at(proms_D, len_A + len_C);
  }
  else // if (invert == false)
  {
    _gen_unit->insert_promoters_at(proms_D, len_A);
    _gen_unit->insert_promoters_at(proms_C, len_A + len_D);
  }

  // Look for new promoters around breakpoints
  _gen_unit->look_for_new_promoters_around(pos_B);
  _gen_unit->look_for_new_promoters_around(pos_B + len_CD);

  _gen_unit->take_ownership_of_all_rnas();
}


/*!
  \brief Looks for an alignment between this and chrom2 in the given sense with max nb_pairs trials. nb_pairs is updated accordingly

  Performs local alignment searches between this and chrom2 around randomly drawn pairs of points.
  The minimum score will be generated according to align_fun_shape and associated parameters for each pair of points.
  The parameter nb_pairs will be updated according to how many trials were necessary for an alignment to be found.

  The sense of the searched alignment can be either DIRECT, INDIRECT or BOTH_SENSE. \
  In the latter case, the sense will be randomly drawn (uniformly between DIRECT and INDIRECT) for each pair of points.
*/
ae_vis_a_vis* ae_dna::search_alignment(ae_dna* chrom2, int32_t& nb_pairs, ae_sense sense)
{
  ae_vis_a_vis* alignment = NULL;
  ae_sense cur_sense = sense; // Which sense (direct or indirect)
  int16_t needed_score;       // Minimum alignement score needed to recombine (stochastic)

  for (; nb_pairs > 0 ; nb_pairs--)
  {
    ///////////////////////////////////////
    //  1) Draw random sense (if needed) //
    ///////////////////////////////////////
    if (sense == BOTH_SENSES)
    {
      cur_sense = (_indiv->_mut_prng->random() < 0.5) ? DIRECT : INDIRECT;
    }

    /////////////////////////////////////////////////////
    // 2) Determine the minimum alignment score needed //
    /////////////////////////////////////////////////////
    if (_indiv->get_align_fun_shape() == LINEAR)
    {
      needed_score = (int16_t) ceil(_indiv->get_align_lin_min() + _indiv->_mut_prng->random() * (_indiv->get_align_lin_max() - _indiv->get_align_lin_min()));
    }
    else
    {
      // I want the probability of rearrangement for an alignment of score <score> to be
      // prob = 1 / (1 + exp(-(score-mean_score)/lambda))
      // The score needed for a rearrangement to take place with a given random drawing is hence
      // needed_score = ceil(-lambda * log(1/rand - 1) + mean)
      needed_score = (int16_t) ceil(- _indiv->get_align_sigm_lambda() * log(1/_indiv->_mut_prng->random() - 1) + _indiv->get_align_sigm_mean());
      if (needed_score < 0) needed_score = 0;
    }

    ///////////////////////////////////////////////////////////////
    // 3) Determine where to look for an alignement (draw seeds) //
    ///////////////////////////////////////////////////////////////
    int32_t seed1 = _indiv->_mut_prng->random(_length);
    int32_t seed2 = _indiv->_mut_prng->random(chrom2->get_length());

    ////////////////////////////////////////////////////////////////////
    // 3) Test the existence of an alignment with a high enough score //
    ////////////////////////////////////////////////////////////////////
    if (cur_sense == DIRECT)
    {
      alignment = ae_align::search_alignment_direct(this, seed1, chrom2, seed2, needed_score);
      if (alignment != NULL)
      {
        return alignment;
      }
    }
    else // if (cur_sense = INDIRECT)
    {
      alignment = ae_align::search_alignment_indirect(this, seed1, chrom2, seed2, needed_score);
      if (alignment != NULL)
      {
        return alignment;
      }
    }
  }

  return NULL;
}

/*!
  \brief Looks for an alignment between this and chrom2 in the given sense around the given positions

  Performs local alignment searches between this and chrom2 around the given positions
  The minimum score will be generated according to align_fun_shape and associated parameters for each pair of points.
  The parameter nb_pairs will be updated according to how many trials were necessary for an alignment to be found.

  The sense of the searched alignment can be either DIRECT, INDIRECT or BOTH_SENSE. \
  In the latter case, the sense will be randomly drawn (uniformly between DIRECT and INDIRECT) for each pair of points.
*/
ae_vis_a_vis* ae_dna::search_alignment_around_positions(ae_dna* chrom2, int32_t chrom1_pos_1, int32_t chrom2_pos_1, ae_sense sense, int8_t& research_sense)
{
  ae_vis_a_vis* alignment = NULL;
  ae_vis_a_vis* tmp_alignment = NULL;
  ae_sense cur_sense = sense; // Which sense (direct or indirect)
  int16_t needed_score;       // Minimum alignement score needed to recombine (stochastic)
  int32_t chrom1_pos_for_research;
  int32_t chrom2_pos_for_research;
  int32_t size_between_two_alignments = 3 * _indiv->get_align_w_zone_h_len();

  ///////////////////////////////////////
  //  1) Draw random sense (if needed) //
  ///////////////////////////////////////
  if (sense == BOTH_SENSES)
  {
    printf("WARNING : Alignment could not be searched in both senses in %s:%d\n", __FILE__, __LINE__);
    return(NULL);
  }

  /////////////////////////////////////////////////////
  // 2) Determine the minimum alignment score needed //
  /////////////////////////////////////////////////////
  if (_indiv->get_align_fun_shape() == LINEAR)
  {
    needed_score = (int16_t) ceil(_indiv->get_align_lin_min() + _indiv->_mut_prng->random() * (_indiv->get_align_lin_max() - _indiv->get_align_lin_min()));
  }
  else
  {
    // I want the probability of rearrangement for an alignment of score <score> to be
    // prob = 1 / (1 + exp(-(score-mean_score)/lambda))
    // The score needed for a rearrangement to take place with a given random drawing is hence
    // needed_score = ceil(-lambda * log(1/rand - 1) + mean)
    needed_score = (int16_t) ceil(- _indiv->get_align_sigm_lambda() * log(1/_indiv->_mut_prng->random() - 1) + _indiv->get_align_sigm_mean());
    if (needed_score < 0) needed_score = 0;
  }

  /////////////////////////////////////////////////////////
  // 3) Determine the sense by which the research begins //
  /////////////////////////////////////////////////////////
  int16_t first_research_sense = (_indiv->_mut_prng->random() < 0.5) ? 1 : -1;
  int16_t second_research_sense = -1*first_research_sense;

  ////////////////////////////////////////////////////////////////////////////////////////
  // 4) Test the first sense for the existence of an alignment with a high enough score //
  ////////////////////////////////////////////////////////////////////////////////////////
  chrom1_pos_for_research = chrom1_pos_1;
  chrom2_pos_for_research = chrom2_pos_1;
  int16_t i = 0;

  while(_indiv->_mut_prng->random() < 1-_exp_m->get_repl_HT_detach_rate())
  {
    //printf("%d longueur genome %d\n",i, this->get_length());
    //printf("chrom1_pos_for_research : %d, chrom2_pos_for_research : %d\n", chrom1_pos_for_research, chrom1_pos_for_research);
    chrom1_pos_for_research = ae_utils::mod(chrom1_pos_for_research + first_research_sense * size_between_two_alignments, this->get_length()) ;
    if (cur_sense == DIRECT)
    {
      chrom2_pos_for_research = ae_utils::mod(chrom2_pos_for_research + first_research_sense * size_between_two_alignments, chrom2->get_length());
      tmp_alignment = ae_align::search_alignment_direct(this, chrom1_pos_for_research, chrom2, chrom2_pos_for_research, needed_score);
    }
    else // if (cur_sense = INDIRECT)
    {
      chrom2_pos_for_research = ae_utils::mod(chrom2_pos_for_research - first_research_sense * size_between_two_alignments, chrom2->get_length());
      tmp_alignment = ae_align::search_alignment_indirect(this, chrom1_pos_for_research, chrom2, chrom2_pos_for_research, needed_score);
    }
    //printf("chrom1_pos_for_research : %d, chrom2_pos_for_research : %d\n", chrom1_pos_for_research, chrom1_pos_for_research);

    if(tmp_alignment == NULL)
    {
      if(alignment != NULL)
      {
        //printf("tmp_alignment == NULL alignment != NULL\n",i);
        research_sense = first_research_sense;
        return alignment;
      }
      else
      {
        //printf("tmp_alignment == NULL alignment == NULL\n",i);
        break;
      }
    }
    else
    {
      if(alignment != NULL)
      {
        //printf("tmp_alignment != NULL alignment != NULL\n",i);
        alignment->copy(tmp_alignment);
      }
      else
      {
        //printf("tmp_alignment != NULL alignment == NULL\n",i);
        alignment = new ae_vis_a_vis(*tmp_alignment);
      }
      delete tmp_alignment;
      chrom1_pos_for_research = alignment->get_i_1();
      chrom2_pos_for_research = alignment->get_i_2();
    }
    i++;
  }

  if(alignment != NULL)
  {

    research_sense = first_research_sense;
    return alignment;
  }

  /////////////////////////////////////////////////////////////////////////////////////////
  // 5) Test the second sense for the existence of an alignment with a high enough score //
  /////////////////////////////////////////////////////////////////////////////////////////
  alignment = NULL;
  chrom1_pos_for_research = chrom1_pos_1;
  chrom2_pos_for_research = chrom2_pos_1;
  i = 0 ;
  while(_indiv->_mut_prng->random() < 1-_exp_m->get_repl_HT_detach_rate())
  {
    chrom1_pos_for_research = ae_utils::mod(chrom1_pos_for_research + second_research_sense * size_between_two_alignments, this->get_length());
    if (cur_sense == DIRECT)
    {
      chrom2_pos_for_research = ae_utils::mod(chrom2_pos_for_research + second_research_sense * size_between_two_alignments, chrom2->get_length());
      tmp_alignment = ae_align::search_alignment_direct(this, chrom1_pos_for_research, chrom2, chrom2_pos_for_research, needed_score);
    }
    else // if (cur_sense = INDIRECT)
    {
      chrom2_pos_for_research = ae_utils::mod(chrom2_pos_for_research - second_research_sense * size_between_two_alignments, chrom2->get_length());
      tmp_alignment = ae_align::search_alignment_indirect(this, chrom1_pos_for_research, chrom2, chrom2_pos_for_research, needed_score);
    }

    if(tmp_alignment == NULL)
    {
      if(alignment != NULL)
      {
        research_sense = second_research_sense;
        return alignment;
      }
      else
      {
        break;
      }
    }
    else
    {
      if(alignment != NULL)
      {
        alignment->copy(tmp_alignment);
      }
      else
      {
        alignment = new ae_vis_a_vis(*tmp_alignment);
      }
      delete tmp_alignment;
      chrom1_pos_for_research = alignment->get_i_1();
      chrom2_pos_for_research = alignment->get_i_2();
    }
    i++;
  }
  if(alignment != NULL)
  {
    research_sense = second_research_sense;
    return alignment;
  }

  return NULL;
}

// =================================================================
//                           Protected Methods
// =================================================================
/**
 * Extract the sequence going from pos_1 (included) to pos_2 (excluded) into a new standalone genetic unit.
 * Promoter lists are created / updated accordingly
 */
void ae_dna::ABCDE_to_ADCBE(int32_t pos_B, int32_t pos_C, int32_t pos_D, int32_t pos_E)
{
  // Rearrange the sequence from ABCDE to ADCBE (complex translocation of segment defined
  // between positions pos_B and pos_D)
  //
  // Segments are identified by pos_x values as shown below.
  //
  // WARNING : Segment C includes nucleotide at pos_D // NOTE : WTF???
  //
  //         A      B        C       D       E                        A      D        C       B        E
  //      |----->=======[>=======>-------[>-------|        =>      |----->-------[>=======>=======[>-------|
  //          pos_B   pos_C    pos_D   pos_E
  //

  // Check points' order and range
  assert(pos_B >= 0 && pos_C >= pos_B && pos_D >= pos_C && pos_E >= pos_D && pos_E <= _length);

  // Compute segment lengths
  int32_t len_A = pos_B;
  int32_t len_B = pos_C - pos_B;
  int32_t len_C = pos_D - pos_C;
  int32_t len_D = pos_E - pos_D;
  int32_t len_E = _length - pos_E;
  int32_t len_AD   = len_A + len_D;
  int32_t len_ADC  = len_AD + len_C;
  int32_t len_ADCB = len_ADC + len_B;

  // Create new sequence
  char* new_genome = new char[_nb_blocks * BLOCK_SIZE];

  memcpy(new_genome,             _data,          len_A * sizeof(char));
  memcpy(&new_genome[len_A],     &_data[pos_D],  len_D * sizeof(char));
  memcpy(&new_genome[len_AD],    &_data[pos_C],  len_C * sizeof(char));
  memcpy(&new_genome[len_ADC],   &_data[pos_B],  len_B * sizeof(char));
  memcpy(&new_genome[len_ADCB],  &_data[pos_E],  len_E * sizeof(char));
  new_genome[_length] = '\0';

  // Replace sequence
  // NB : The size of the genome doesn't change. Therefore, we don't need to update _length and _nb_blocks
  delete [] _data;
  _data = new_genome;


  // ========== Update promoter list ==========
  if (_length >= PROM_SIZE)
  {
    // Remove promoters that include a breakpoint
    _gen_unit->remove_promoters_around(pos_B);
    _gen_unit->remove_promoters_around(pos_C);
    _gen_unit->remove_promoters_around(pos_D);
    _gen_unit->remove_promoters_around(pos_E);

    // Create temporary lists for promoters to move and/or invert
    Promoters2Strands promoters_B = {{},{}};
    Promoters2Strands promoters_C = {{},{}};
    Promoters2Strands promoters_D = {{},{}};

    // Extract promoters that are totally included in each segment to be moved
    // and shift them to their new positions
    if (len_B >= PROM_SIZE)
    {
      _gen_unit->extract_promoters_included_in(pos_B, pos_C, promoters_B);
      GeneticUnit::shift_promoters(promoters_B, len_D + len_C,  _gen_unit->get_dna()->get_length());
    }
    if (len_C >= PROM_SIZE)
    {
      _gen_unit->extract_promoters_included_in(pos_C, pos_D, promoters_C);
      GeneticUnit::shift_promoters(promoters_C, len_D - len_B,  _gen_unit->get_dna()->get_length());
    }
    if (len_D >= PROM_SIZE)
    {
      _gen_unit->extract_promoters_included_in(pos_D, pos_E, promoters_D);
      GeneticUnit::shift_promoters(promoters_D, -len_B - len_C, _gen_unit->get_dna()->get_length());
    }

    // Reinsert the shifted promoters
    _gen_unit->insert_promoters(promoters_B);
    _gen_unit->insert_promoters(promoters_C);
    _gen_unit->insert_promoters(promoters_D);

    // 5) Look for new promoters including a breakpoint
    _gen_unit->look_for_new_promoters_around(len_A);
    _gen_unit->look_for_new_promoters_around(len_AD);
    _gen_unit->look_for_new_promoters_around(len_ADC);
    _gen_unit->look_for_new_promoters_around(len_ADCB);
  }
}

void ae_dna::ABCDE_to_ADBpCpE(int32_t pos_B, int32_t pos_C, int32_t pos_D, int32_t pos_E)
{
  // Rearrange the sequence from ABCDE to ADBpCpE (complex translocation with inversion
  // of segment defined between positions pos_B and pos_D)
  // Bp (resp Cp) stands for inverted B (resp C)
  //
  // Segments are identified by pos_x values as shown below.
  //
  // WARNING : Segment C includes nucleotide at pos_D // NOTE : WTF???
  //
  //         A      B        C       D        E                       A      D        Bp      Cp       E
  //      |----->=======[>=======>-------<]-------|        =>      |----->-------<]=======<=======<]-------|
  //          pos_B   pos_C    pos_D   pos_E
  //

  // Check points' order and range
  assert(pos_B >= 0 && pos_C >= pos_B && pos_D >= pos_C && pos_E >= pos_D && pos_E <= _length);

  // Compute segment lengths
  int32_t len_A = pos_B;
  int32_t len_B = pos_C - pos_B;
  int32_t len_C = pos_D - pos_C;
  int32_t len_D = pos_E - pos_D;
  int32_t len_E = _length - pos_E;
  int32_t len_AD   = len_A + len_D;
  int32_t len_ADB  = len_AD + len_B;
  int32_t len_ADBC = len_ADB + len_C;

  // Create new sequence
  char* new_genome = new char[_nb_blocks * BLOCK_SIZE];

  // Copy segments A and D
  memcpy(new_genome,         _data,            len_A * sizeof(char));
  memcpy(&new_genome[len_A], &_data[pos_D],  len_D * sizeof(char));


  // Build Bp and put it in the new genome
  char* inverted_segment = new char[len_B+1];

  for (int32_t i = 0, j = pos_C - 1 ; i < len_B ; i++, j--)
  {
    if (_data[j] == '0') inverted_segment[i] = '1';
    else                   inverted_segment[i] = '0';
  }
  inverted_segment[len_B] = '\0';

  memcpy(&new_genome[len_AD], inverted_segment, len_B * sizeof(char));

  delete [] inverted_segment;


  // Build Cp and put it in the new genome
  inverted_segment = new char[len_C+1];

  for (int32_t i = 0, j = pos_D - 1 ; i < len_C ; i++, j--)
  {
    if (_data[j] == '0') inverted_segment[i] = '1';
    else                   inverted_segment[i] = '0';
  }
  inverted_segment[len_C] = '\0';

  memcpy(&new_genome[len_ADB], inverted_segment, len_C * sizeof(char));

  delete [] inverted_segment;

  // Copy segment E into the new genome
  memcpy(&new_genome[len_ADBC], &_data[pos_E], len_E * sizeof(char));
  new_genome[_length] = '\0';


  // Replace sequence
  delete [] _data;
  _data = new_genome;


  // ========== Update promoter list ==========
  if (_length >= PROM_SIZE)
  {
    // Remove promoters that include a breakpoint
    _gen_unit->remove_promoters_around(pos_B);
    _gen_unit->remove_promoters_around(pos_C);
    _gen_unit->remove_promoters_around(pos_D);
    _gen_unit->remove_promoters_around(pos_E);

    // Create temporary lists for promoters to move and/or invert
    Promoters2Strands promoters_B = {{},{}};
    Promoters2Strands promoters_C = {{},{}};
    Promoters2Strands promoters_D = {{},{}};

    // 2) Extract promoters that are totally included in each segment to be moved (B, C and D)
    if (len_B >= PROM_SIZE)
    {
      _gen_unit->extract_promoters_included_in(pos_B, pos_C, promoters_B);
    }
    if (len_C >= PROM_SIZE)
    {
      _gen_unit->extract_promoters_included_in(pos_C, pos_D, promoters_C);
    }
    if (len_D >= PROM_SIZE)
    {
      _gen_unit->extract_promoters_included_in(pos_D, pos_E, promoters_D);
    }

    // 3a) Invert promoters of segments B and C
    GeneticUnit::invert_promoters(promoters_B, pos_B, pos_C);
    GeneticUnit::invert_promoters(promoters_C, pos_C, pos_D);

    // 3b) Shift these promoters positions
    GeneticUnit::shift_promoters(promoters_B, len_D,           _gen_unit->get_dna()->get_length());
    GeneticUnit::shift_promoters(promoters_C, len_D,           _gen_unit->get_dna()->get_length());
    GeneticUnit::shift_promoters(promoters_D, -len_B - len_C,  _gen_unit->get_dna()->get_length());

    // 4) Reinsert the shifted promoters
    _gen_unit->insert_promoters(promoters_C);
    _gen_unit->insert_promoters(promoters_B);
    _gen_unit->insert_promoters(promoters_D);

    // 5) Look for new promoters including a breakpoint
    _gen_unit->look_for_new_promoters_around(len_A);
    _gen_unit->look_for_new_promoters_around(len_AD);
    _gen_unit->look_for_new_promoters_around(len_ADB);
    _gen_unit->look_for_new_promoters_around(len_ADBC);
  }
}

void ae_dna::ABCDE_to_ACpDpBE(int32_t pos_B, int32_t pos_C, int32_t pos_D, int32_t pos_E)
{
  // Rearrange the sequence from ABCDE to ACpDpBE (complex translocation with inversion
  // of segment defined between positions pos_C and pos_E)
  // Cp (resp Dp) stands for inverted C (resp D)
  //
  // Segments are identified by pos_x values as shown below.
  //
  // WARNING : Segment D includes nucleotide at pos_E // NOTE : WTF???
  //
  //         A      B        C       D       E                        A       C'      D'       B       E
  //      |----<]-------->=======[>=======>-------|        =>      |-----<]=======>=======<]------->-------|
  //          pos_B    pos_C    pos_D   pos_E
  //

  // Check points' order and range
  assert(pos_B >= 0 && pos_C >= pos_B && pos_D >= pos_C && pos_E >= pos_D && pos_E <= _length);

  // Compute segment lengths
  int32_t len_A = pos_B;
  int32_t len_B = pos_C - pos_B;
  int32_t len_C = pos_D - pos_C;
  int32_t len_D = pos_E - pos_D;
  int32_t len_E = _length - pos_E;
  int32_t len_AC   = len_A + len_C;
  int32_t len_ACD  = len_AC + len_D;
  int32_t len_ACDB = len_ACD + len_B;

  // Create new sequence
  char* new_genome = new char[_nb_blocks * BLOCK_SIZE];

  // Copy segment A
  memcpy(new_genome, _data, len_A * sizeof(char));


  // Build Cp and put it in the new genome
  char* inverted_segment = new char[len_C+1];

  for (int32_t i = 0, j = pos_D - 1 ; i < len_C ; i++, j--)
  {
    if (_data[j] == '0') inverted_segment[i] = '1';
    else                   inverted_segment[i] = '0';
  }
  inverted_segment[len_C] = '\0';

  memcpy(&new_genome[len_A], inverted_segment, len_C * sizeof(char));

  delete [] inverted_segment;


  // Build Dp and put it in the new genome
  inverted_segment = new char[len_D+1];

  for (int32_t i = 0, j = pos_E - 1 ; i < len_D ; i++, j--)
  {
    if (_data[j] == '0') inverted_segment[i] = '1';
    else                   inverted_segment[i] = '0';
  }
  inverted_segment[len_D] = '\0';

  memcpy(&new_genome[len_AC], inverted_segment, len_D * sizeof(char));

  delete [] inverted_segment;

  // Copy segments B and E
  memcpy(&new_genome[len_ACD],   &_data[pos_B], len_B * sizeof(char));
  memcpy(&new_genome[len_ACDB],  &_data[pos_E], len_E * sizeof(char));
  new_genome[_length] = '\0';


  // Replace sequence
  delete [] _data;
  _data = new_genome;


  // ========== Update promoter list ==========
  // 1) Remove promoters that include a breakpoint
  // 2) Extract promoters that are totally included in each segment to be moved (B, C and D)
  // 3) Shift (and invert when needed) these promoters positions
  // 4) Reinsert the shifted promoters
  // 5) Look for new promoters including a breakpoint
  if (_length >= PROM_SIZE)
  {
    // 1) Remove promoters that include a breakpoint
    _gen_unit->remove_promoters_around(pos_B);
    _gen_unit->remove_promoters_around(pos_C);
    _gen_unit->remove_promoters_around(pos_D);
    _gen_unit->remove_promoters_around(pos_E);

    // Create temporary lists for promoters to move and/or invert
    Promoters2Strands promoters_B = {{},{}};
    Promoters2Strands promoters_C = {{},{}};
    Promoters2Strands promoters_D = {{},{}};

    // 2) Extract promoters that are totally included in each segment to be moved (B, C and D)
    if (len_B >= PROM_SIZE)
    {
      _gen_unit->extract_promoters_included_in(pos_B, pos_C, promoters_B);
    }
    if (len_C >= PROM_SIZE)
    {
      _gen_unit->extract_promoters_included_in(pos_C, pos_D, promoters_C);
    }
    if (len_D >= PROM_SIZE)
    {
      _gen_unit->extract_promoters_included_in(pos_D, pos_E, promoters_D);
    }

    // 3a) Invert promoters of segments C and D
    GeneticUnit::invert_promoters(promoters_C, pos_C, pos_D);
    GeneticUnit::invert_promoters(promoters_D, pos_D, pos_E);

    // 3b) Shift these promoters positions
    GeneticUnit::shift_promoters(promoters_B, len_C + len_D, _gen_unit->get_dna()->get_length());
    GeneticUnit::shift_promoters(promoters_C, -len_B,        _gen_unit->get_dna()->get_length());
    GeneticUnit::shift_promoters(promoters_D, -len_B,        _gen_unit->get_dna()->get_length());

    // 4) Reinsert the shifted promoters
    _gen_unit->insert_promoters(promoters_B);
    _gen_unit->insert_promoters(promoters_D);
    _gen_unit->insert_promoters(promoters_C);

    // 5) Look for new promoters including a breakpoint
    _gen_unit->look_for_new_promoters_around(len_A);
    _gen_unit->look_for_new_promoters_around(len_AC);
    _gen_unit->look_for_new_promoters_around(len_ACD);
    _gen_unit->look_for_new_promoters_around(len_ACDB);
  }
}

void ae_dna::inter_GU_ABCDE_to_ACDBE(int32_t pos_B, int32_t pos_C, int32_t pos_E)
{
  // Check points' order and range
  assert((pos_B >= 0 && pos_C >= pos_B) && (pos_E >= 0));

  if (pos_B != pos_C)
  {
    // Useful values
    Individual * indiv            = _indiv;
    GeneticUnit& chromosome     = indiv->get_genetic_unit_nonconst(0);
    GeneticUnit& plasmid        = indiv->get_genetic_unit_nonconst(1);
    GeneticUnit& destination_GU = (_gen_unit == &chromosome)? plasmid : chromosome;

    // Compute segment lengths
    int32_t len_A = pos_B;
    int32_t len_B = pos_C - pos_B;
    int32_t len_C = _length - pos_C;
    int32_t len_D = pos_E;
    int32_t len_E = destination_GU.get_dna()->get_length() - pos_E;
    int32_t len_AC = len_A + len_C;
    int32_t len_DB = len_D + len_B;
    int32_t len_DBE = len_DB + len_E;


    // Create the new sequence of this genetic unit
    int32_t tmp = ae_string::nb_blocks(len_AC);
    char* new_sequence_this = new char[tmp * BLOCK_SIZE];

    memcpy(new_sequence_this,           _data,          len_A * sizeof(char));
    memcpy(&new_sequence_this[len_A],   &_data[pos_C],  len_C * sizeof(char));
    new_sequence_this[len_AC]='\0';

    // Create the new sequence of the destination genetic unit
    tmp = ae_string::nb_blocks(len_DBE) * BLOCK_SIZE;
    char* new_sequence_dest = new char[tmp];
    char* dest_GU_former_seq = (char*) destination_GU.get_dna()->get_data();

    memcpy(new_sequence_dest,           dest_GU_former_seq,          len_D * sizeof(char));
    memcpy(&new_sequence_dest[len_D],   &_data[pos_B],               len_B * sizeof(char));
    memcpy(&new_sequence_dest[len_DB],  &dest_GU_former_seq[pos_E],  len_E * sizeof(char));
    new_sequence_dest[len_DBE]='\0';



    // ========== Update promoter list ==========
    // 1) Remove promoters that include a breakpoint
    // 2) Extract promoters that are totally included in each segment to be moved (B, C and E)
    // NB : Sequences have to be updated at this stage in order to have the correct lengths when managing new promoter positions
    // ........
    // 3) Shift these promoters positions
    // 4) Reinsert the shifted promoters
    // 5) Look for new promoters including a breakpoint


    // 1) Remove promoters that include a breakpoint
    _gen_unit->remove_promoters_around(pos_B);
    _gen_unit->remove_promoters_around(pos_C);
    destination_GU.remove_promoters_around(pos_E);

    // Create temporary lists for promoters to move and/or invert
    Promoters2Strands promoters_B = {{},{}};

    // 2) Extract promoters that are totally included in each segment to be moved (B, C and E)
    if (len_B >= PROM_SIZE)
    {
      _gen_unit->extract_promoters_included_in(pos_B, pos_C, promoters_B);
    }



    // ========== Replace sequences ==========
    set_data(new_sequence_this, len_AC);
    destination_GU.get_dna()->set_data(new_sequence_dest, len_DBE);


    // 3) Shift these promoters positions
    GeneticUnit::shift_promoters(promoters_B, len_D - len_A, destination_GU.get_dna()->get_length());

    // Reassign promoters to their new genetic unit
    for (auto& strand: {LEADING, LAGGING})
      for (auto& rna: promoters_B[strand])
        rna.set_genetic_unit(&destination_GU);

    // Shift the promoters of sequences C and E
    _gen_unit->move_all_promoters_after(pos_C, -len_B);
    destination_GU.move_all_promoters_after(pos_E, len_B);

    // 4) Reinsert the shifted promoters
    destination_GU.insert_promoters(promoters_B);

    // 5) Look for new promoters including a breakpoint
    _gen_unit->look_for_new_promoters_around(0);
    _gen_unit->look_for_new_promoters_around(len_A);
    destination_GU.look_for_new_promoters_around(0);
    destination_GU.look_for_new_promoters_around(len_D);
    destination_GU.look_for_new_promoters_around(len_DB);
  }
}

void ae_dna::inter_GU_ABCDE_to_BDCAE(int32_t pos_B, int32_t pos_C, int32_t pos_E)
{
  // Check points' order and range
  assert((pos_B >= 0 && pos_C >= pos_B) && (pos_E >= 0));

  // Compute segment lengths
  int32_t len_A = pos_B;
  int32_t len_B = pos_C - pos_B;
  int32_t len_C = _length - pos_C;
  //~ int32_t len_ABC = _length;
  int32_t len_DA = len_A + pos_E;

  inter_GU_ABCDE_to_ACDBE(0, pos_B, pos_E);
  inter_GU_ABCDE_to_ACDBE(len_B, (len_B+len_C), len_DA);
}
} // namespace aevol
