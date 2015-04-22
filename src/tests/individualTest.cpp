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
#include <gtest/gtest.h>
#include <string.h>



// =================================================================
//                            Project Files
// =================================================================
#include "ae_individual.h"
#include "macros.h"
#include "genetic_unit.h"
#include "ae_rna.h"
#include "ae_protein.h"
#include "ae_params_mut.h"



// ===========================================================================
//                             Declare Used Namespaces
// ===========================================================================
using namespace aevol;




//############################################################################
//                                                                           #
//                         Class individualTest                              #
//                                                                           #
//############################################################################
class individualTest : public testing::Test
{
 protected:
  virtual void SetUp(void);
  virtual void TearDown(void);

  ae_individual* indiv1;
  ae_individual* indiv2;
  ae_individual* indiv3;
  ae_individual* indiv4;
};

// ===========================================================================
//                                 Public Methods
// ===========================================================================
void individualTest::SetUp(void)
{
  // Build ad-hoc genomes
  // (and reverse to test the same things on the lagging strand.):
  //
  // indiv1: (AS + prom + AS + AG + AS + term + AS + prom + AS)
  // indiv2: reverse
  // indiv3: (AS + AG + AS + term + AS + prom + AS)
  // indiv4: reverse
  //
  // AS = Arbitrary Sequence
  // AG = Arbitrary Gene
  // Do not modify the sequences !
  char as[5][10] = {
    "0011",
    "11101",
    "110011",
    "11000",
    "000101"
  };
  char gene[255];
  sprintf(gene, "%s0011000100110110010001", SHINE_DAL_SEQ);
  char term[TERM_SIZE+1] = "01000001101";
  char prom[2][23] = {
    "0101010001110110010110", // dist from consensus: 2 => basal level: 0.6
    "0101011001110010010010"  // dist from consensus: 1 => basal level: 0.8
  };
  char* genome = new char[1024];
  sprintf( genome, "%s%s%s%s%s%s%s%s%s", as[0], prom[0], as[1], gene, as[2],
           term, as[3], prom[1], as[4]);

  // Build indiv1
  ae_params_mut params_mut;
  indiv1 = new ae_individual(NULL, NULL, NULL, &params_mut, 1.0, 10, 1000, false, 1, "anon-strain-1", 0);
  indiv1->add_GU(genome, strlen(genome));
  genome = NULL;

  // Do transcription and translation
  indiv1->do_transcription();
  indiv1->do_translation();



  // Build indiv2
  genome = indiv1->get_genetic_unit(0)->get_dna()->get_subsequence(0,0,LAGGING);
  indiv2 = new ae_individual(NULL, NULL, NULL, &params_mut, 1.0, 10, 1000, false, 1, "anon-strain-2", 0);
  indiv2->add_GU(genome, strlen(genome));
  genome = NULL;

  // Do transcription and translation
  indiv2->do_transcription();
  indiv2->do_translation();




  // Build indiv3
  genome = new char[1024];
  sprintf( genome, "%s%s%s%s%s%s%s", as[0], gene, as[1], term, as[2], prom[1], as[3]);
  indiv3 = new ae_individual(NULL, NULL, NULL, &params_mut, 1.0, 10, 1000, false, 1, "anon-strain-3", 0);
  indiv3->add_GU(genome, strlen(genome));
  genome = NULL;

  // Do transcription and translation
  indiv3->do_transcription();
  indiv3->do_translation();




  // Build indiv4
  genome = indiv3->get_genetic_unit(0)->get_dna()->get_subsequence(0,0,LAGGING);
  indiv4 = new ae_individual(NULL, NULL, NULL, &params_mut, 1.0, 10, 1000, false, 1, "anon-strain-4", 0);
  indiv4->add_GU(genome, strlen(genome));
  genome = NULL;

  // Do transcription and translation
  indiv4->do_transcription();
  indiv4->do_translation();


  // ***************************************************************************
  // The following commented code allows to print stuff about rnas and proteins

  // printf("%"PRId32" rnas and %"PRId32" prots\n", indiv4->get_rna_list()->get_nb_elts(), indiv4->get_protein_list()->get_nb_elts());
  // ae_list_node<ae_rna*>* rna_node = indiv4->get_rna_list()->get_first();
  // while (rna_node != NULL)
  // {
  //   printf("%s rna at pos %"PRId32" (%f, %"PRId32")\n",
  //           rna_node->get_obj()->get_strand() == LEADING ? "LEADING":"LAGGING",
  //           rna_node->get_obj()->get_promoter_pos(),
  //           rna_node->get_obj()->get_basal_level(),
  //           rna_node->get_obj()->get_transcript_length());

  //   rna_node = rna_node->get_next();
  // }

  // ae_list_node<ae_protein*>* protein_node = indiv4->get_protein_list()->get_first();
  // while (protein_node != NULL)
  // {
  //   printf("%s protein at pos %"PRId32" (length: %"PRId32", concentr: %f, nb_rnas: %"PRId32")\n",
  //           protein_node->get_obj()->get_strand() == LEADING ? "LEADING":"LAGGING",
  //           protein_node->get_obj()->get_shine_dal_pos(),
  //           protein_node->get_obj()->get_length(),
  //           protein_node->get_obj()->get_concentration(),
  //           protein_node->get_obj()->get_rna_list()->get_nb_elts());

  //   protein_node = protein_node->get_next();
  // }
}

void individualTest::TearDown(void)
{
  delete indiv1;
  delete indiv2;
  delete indiv3;
  delete indiv4;
}


TEST_F(individualTest, TestIndiv1)
{
  // Check genome size
  EXPECT_EQ(109, indiv1->get_amount_of_dna());
  EXPECT_EQ(109, indiv1->get_genetic_unit_seq_length(0) );

  // Check RNA list
  ae_list<ae_rna*>* rna_list = indiv1->get_rna_list();
  EXPECT_EQ(2, rna_list->get_nb_elts() );
  ae_rna* rna = rna_list->get_first()->get_obj();
  EXPECT_EQ(LEADING, rna->get_strand());
  EXPECT_EQ(4, rna->get_promoter_pos());
  EXPECT_FLOAT_EQ(0.6, rna->get_basal_level());
  EXPECT_EQ(50, rna->get_transcript_length());
  rna = rna_list->get_last()->get_obj();
  EXPECT_EQ(LEADING, rna->get_strand());
  EXPECT_EQ(81, rna->get_promoter_pos());
  EXPECT_FLOAT_EQ(0.8, rna->get_basal_level());
  EXPECT_EQ(82, rna->get_transcript_length());

  // Check protein list
  ae_list<ae_protein*>* prot_list = indiv1->get_protein_list();
  EXPECT_EQ(1, prot_list->get_nb_elts());
  ae_protein* prot = prot_list->get_first()->get_obj();
  EXPECT_EQ(LEADING, prot->get_strand());
  EXPECT_EQ(31, prot->get_shine_dal_pos());
  EXPECT_EQ(4, prot->get_length());
  EXPECT_FLOAT_EQ(1.4, prot->get_concentration());
  EXPECT_EQ(2, prot->get_rna_list()->get_nb_elts());
}

TEST_F(individualTest, TestIndiv2)
{
  // Check genome size
  EXPECT_EQ(109, indiv2->get_amount_of_dna());
  EXPECT_EQ(109, indiv2->get_genetic_unit_seq_length(0));

  // Check RNA list
  ae_list<ae_rna*>* rna_list = indiv2->get_rna_list();
  EXPECT_EQ(2, rna_list->get_nb_elts());
  ae_rna* rna = rna_list->get_first()->get_obj();
  EXPECT_EQ(LAGGING, rna->get_strand());
  EXPECT_EQ(104, rna->get_promoter_pos());
  EXPECT_FLOAT_EQ(0.6, rna->get_basal_level());
  EXPECT_EQ(50, rna->get_transcript_length());
  rna = rna_list->get_last()->get_obj();
  EXPECT_EQ(LAGGING, rna->get_strand());
  EXPECT_EQ(27, rna->get_promoter_pos());
  EXPECT_FLOAT_EQ(0.8, rna->get_basal_level());
  EXPECT_EQ(82, rna->get_transcript_length());

  // Check protein list
  ae_list<ae_protein*>* prot_list = indiv2->get_protein_list();
  EXPECT_EQ(1, prot_list->get_nb_elts());
  ae_protein* prot = prot_list->get_first()->get_obj();
  EXPECT_EQ(LAGGING, prot->get_strand());
  EXPECT_EQ(77, prot->get_shine_dal_pos());
  EXPECT_EQ(4, prot->get_length());
  EXPECT_FLOAT_EQ(1.4, prot->get_concentration());
  EXPECT_EQ(2, prot->get_rna_list()->get_nb_elts());
}

TEST_F(individualTest, TestIndiv3)
{
  // Check genome size
  EXPECT_EQ(81, indiv3->get_amount_of_dna());
  EXPECT_EQ(81, indiv3->get_genetic_unit_seq_length(0));

  // Check RNA list
  ae_list<ae_rna*>* rna_list = indiv3->get_rna_list();
  EXPECT_EQ(1, rna_list->get_nb_elts());
  ae_rna* rna = rna_list->get_first()->get_obj();
  EXPECT_EQ(LEADING, rna->get_strand());
  EXPECT_EQ(54, rna->get_promoter_pos());
  EXPECT_FLOAT_EQ(0.8, rna->get_basal_level());
  EXPECT_EQ(42, rna->get_transcript_length());

  // Check protein list
  ae_list<ae_protein*>* prot_list = indiv3->get_protein_list();
  EXPECT_EQ(1, prot_list->get_nb_elts());
  ae_protein* prot = prot_list->get_first()->get_obj();
  EXPECT_EQ(LEADING, prot->get_strand());
  EXPECT_EQ(4, prot->get_shine_dal_pos());
  EXPECT_EQ(4, prot->get_length());
  EXPECT_FLOAT_EQ(0.8, prot->get_concentration());
  EXPECT_EQ(1, prot->get_rna_list()->get_nb_elts());
}

TEST_F(individualTest, TestIndiv4)
{
  // Check genome size
  EXPECT_EQ(81, indiv4->get_amount_of_dna());
  EXPECT_EQ(81, indiv4->get_genetic_unit_seq_length(0));

  // Check RNA list
  ae_list<ae_rna*>* rna_list = indiv4->get_rna_list();
  EXPECT_EQ(1, rna_list->get_nb_elts());
  ae_rna* rna = rna_list->get_first()->get_obj();
  EXPECT_EQ(LAGGING, rna->get_strand());
  EXPECT_EQ(26, rna->get_promoter_pos());
  EXPECT_FLOAT_EQ(0.8, rna->get_basal_level());
  EXPECT_EQ(42, rna->get_transcript_length());

  // Check protein list
  ae_list<ae_protein*>* prot_list = indiv4->get_protein_list();
  EXPECT_EQ(1, prot_list->get_nb_elts());
  ae_protein* prot = prot_list->get_first()->get_obj();
  EXPECT_EQ(LAGGING, prot->get_strand());
  EXPECT_EQ(76, prot->get_shine_dal_pos());
  EXPECT_EQ(4, prot->get_length());
  EXPECT_FLOAT_EQ(0.8, prot->get_concentration());
  EXPECT_EQ(1, prot->get_rna_list()->get_nb_elts());
}

// ===========================================================================
//                                Protected Methods
// ===========================================================================

// ===========================================================================
//                              Non inline accessors
// ===========================================================================
