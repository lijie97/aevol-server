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
//                              Includes
// =================================================================
#include <inttypes.h>
#include <cstring>

#include <list>
#include <memory>

#include <gtest/gtest.h>

#include "Individual.h"
#include "macros.h"
#include "GeneticUnit.h"
#include "Rna.h"
#include "Protein.h"
#include "MutationParams.h"
#include "../libaevol/macros.h"
#include "ExpSetup.h"

using namespace aevol;

using std::list;




//############################################################################
//                                                                           #
//                         Class IndividualTest                              #
//                                                                           #
//############################################################################
class IndividualTest : public testing::Test
{
 protected:
  virtual void SetUp(void);
  virtual void TearDown(void);

  Individual* indiv1;
  Individual* indiv2;
  Individual* indiv3;
  Individual* indiv4;
};

// ===========================================================================
//                                 Public Methods
// ===========================================================================
void IndividualTest::SetUp(void)
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

  // Define a few arbitrary sequences
  char as[5][10] = {
    "0011",
    "11101",
    "110011",
    "11000",
    "000101"
  };

  // Define an arbitrary gene
  char gene[255];
  sprintf(gene, "%s0011000100110110010001", SHINE_DAL_SEQ);

  // Define an arbitrary terminator
  char term[TERM_SIZE+1] = "01000001101";

  // Define a few arbitrary promoters
  char prom[2][23] = {
    "0101010001110110010110", // dist from consensus: 2 => basal level: 0.6
    "0101011001110010010010"  // dist from consensus: 1 => basal level: 0.8
  };

  // Construct a gene with these arbitrary sequences
  char* genome = new char[1024];
  sprintf(genome, "%s%s%s%s%s%s%s%s%s", as[0], prom[0], as[1], gene, as[2],
           term, as[3], prom[1], as[4]);

  // Build indiv1
  MutationParams params_mut;
  indiv1 = new Individual(nullptr, nullptr, nullptr, std::make_shared<MutationParams>(params_mut), 1.0, 10, 1000, false, 1, "anon-strain-1", 0);
  indiv1->add_GU(genome, strlen(genome));
  genome = NULL;

  // Do transcription and translation
  indiv1->do_transcription();
  indiv1->do_translation();



  // Build indiv2
  // Reverse the whole genome
  genome = indiv1->genetic_unit(0).dna()->subsequence(0,0,LAGGING);
  indiv2 = new Individual(nullptr, nullptr, nullptr, std::make_shared<MutationParams>(params_mut), 1.0, 10, 1000, false, 1, "anon-strain-2", 0);
  indiv2->add_GU(genome, strlen(genome));
  genome = NULL;

  // Do transcription and translation
  indiv2->do_transcription();
  indiv2->do_translation();




  // Build indiv3
  genome = new char[1024];
  sprintf(genome, "%s%s%s%s%s%s%s", as[0], gene, as[1], term, as[2], prom[1], as[3]);
  indiv3 = new Individual(nullptr, nullptr, nullptr, std::make_shared<MutationParams>(params_mut), 1.0, 10, 1000, false, 1, "anon-strain-3", 0);
  indiv3->add_GU(genome, strlen(genome));
  genome = NULL;

  // Do transcription and translation
  indiv3->do_transcription();
  indiv3->do_translation();




  // Build indiv4
  genome = indiv3->genetic_unit(0).dna()->subsequence(0,0,LAGGING);
  indiv4 = new Individual(nullptr, nullptr, nullptr, std::make_shared<MutationParams>(params_mut), 1.0, 10, 1000, false, 1, "anon-strain-4", 0);
  indiv4->add_GU(genome, strlen(genome));
  genome = NULL;

  // Do transcription and translation
  indiv4->do_transcription();
  indiv4->do_translation();


  // ***************************************************************************
  // The following commented code allows to print stuff about rnas and proteins

  // printf("%"PRId32" rnas and %"PRId32" prots\n", indiv4->rna_list()->size(), indiv4->protein_list()->size());
  // list_node<Rna*>* rna_node = indiv4->rna_list()->get_first();
  // while (rna_node != NULL)
  // {
  //   printf("%s rna at pos %"PRId32" (%f, %"PRId32")\n",
  //           rna_node->get_obj()->strand() == LEADING ? "LEADING":"LAGGING",
  //           rna_node->get_obj()->promoter_pos(),
  //           rna_node->get_obj()->basal_level(),
  //           rna_node->get_obj()->transcript_length());

  //   rna_node = rna_node->get_next();
  // }

  // list_node<Protein*>* protein_node = indiv4->protein_list()->get_first();
  // while (protein_node != NULL)
  // {
  //   printf("%s protein at pos %"PRId32" (length: %"PRId32", concentr: %f, nb_rnas: %"PRId32")\n",
  //           protein_node->get_obj()->strand() == LEADING ? "LEADING":"LAGGING",
  //           protein_node->get_obj()->shine_dal_pos(),
  //           protein_node->get_obj()->length(),
  //           protein_node->get_obj()->concentration(),
  //           protein_node->get_obj()->rna_list()->size());

  ExpSetup* exp_s = new ExpSetup(nullptr);
  for (int i = 0; i <= 1; i++) {
    exp_s->set_fuzzy_flavor(i);

    if (FuzzyFactory::fuzzyFactory == NULL)
      FuzzyFactory::fuzzyFactory = new FuzzyFactory(exp_s);
    else {
      delete FuzzyFactory::fuzzyFactory;
      FuzzyFactory::fuzzyFactory = new FuzzyFactory(exp_s);
    }
    printf("F %d\n",i);
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

    // Define a few arbitrary sequences
    char as[5][10] = {
        "0011",
        "11101",
        "110011",
        "11000",
        "000101"
    };

    // Define an arbitrary gene
    char gene[255];
    sprintf(gene, "%s0011000100110110010001", SHINE_DAL_SEQ);

    // Define an arbitrary terminator
    char term[TERM_SIZE + 1] = "01000001101";

    // Define a few arbitrary promoters
    char prom[2][23] = {
        "0101010001110110010110", // dist from consensus: 2 => basal level: 0.6
        "0101011001110010010010"  // dist from consensus: 1 => basal level: 0.8
    };

    // Construct a gene with these arbitrary sequences
    char* genome = new char[1024];
    sprintf(genome, "%s%s%s%s%s%s%s%s%s", as[0], prom[0], as[1], gene, as[2],
            term, as[3], prom[1], as[4]);

    // Build indiv1
    MutationParams params_mut;
    indiv1 = new Individual(nullptr, nullptr, nullptr,
                            std::make_shared<MutationParams>(params_mut), 1.0,
                            10, 1000, false, 1, "anon-strain-1", 0);
    indiv1->add_GU(genome, strlen(genome));
    genome = NULL;

    // Do transcription and translation
    indiv1->do_transcription();
    indiv1->do_translation();



    // Build indiv2
    // Reverse the whole genome
    genome = indiv1->genetic_unit(0).dna()->subsequence(0, 0,
                                                                    LAGGING);
    indiv2 = new Individual(nullptr, nullptr, nullptr,
                            std::make_shared<MutationParams>(params_mut), 1.0,
                            10, 1000, false, 1, "anon-strain-2", 0);
    indiv2->add_GU(genome, strlen(genome));
    genome = NULL;

    // Do transcription and translation
    indiv2->do_transcription();
    indiv2->do_translation();




    // Build indiv3
    genome = new char[1024];
    sprintf(genome, "%s%s%s%s%s%s%s", as[0], gene, as[1], term, as[2], prom[1],
            as[3]);
    indiv3 = new Individual(nullptr, nullptr, nullptr,
                            std::make_shared<MutationParams>(params_mut), 1.0,
                            10, 1000, false, 1, "anon-strain-3", 0);
    indiv3->add_GU(genome, strlen(genome));
    genome = NULL;

    // Do transcription and translation
    indiv3->do_transcription();
    indiv3->do_translation();




    // Build indiv4
    genome = indiv3->genetic_unit(0).dna()->subsequence(0, 0,
                                                                    LAGGING);
    indiv4 = new Individual(nullptr, nullptr, nullptr,
                            std::make_shared<MutationParams>(params_mut), 1.0,
                            10, 1000, false, 1, "anon-strain-4", 0);
    indiv4->add_GU(genome, strlen(genome));
    genome = NULL;

    // Do transcription and translation
    indiv4->do_transcription();
    indiv4->do_translation();


    // ***************************************************************************
    // The following commented code allows to print stuff about rnas and proteins

    // printf("%"PRId32" rnas and %"PRId32" prots\n", indiv4->rna_list()->size(), indiv4->protein_list()->size());
    // list_node<Rna*>* rna_node = indiv4->rna_list()->get_first();
    // while (rna_node != NULL)
    // {
    //   printf("%s rna at pos %"PRId32" (%f, %"PRId32")\n",
    //           rna_node->get_obj()->strand() == LEADING ? "LEADING":"LAGGING",
    //           rna_node->get_obj()->promoter_pos(),
    //           rna_node->get_obj()->basal_level(),
    //           rna_node->get_obj()->transcript_length());

    //   rna_node = rna_node->get_next();
    // }

    // list_node<Protein*>* protein_node = indiv4->protein_list()->get_first();
    // while (protein_node != NULL)
    // {
    //   printf("%s protein at pos %"PRId32" (length: %"PRId32", concentr: %f, nb_rnas: %"PRId32")\n",
    //           protein_node->get_obj()->strand() == LEADING ? "LEADING":"LAGGING",
    //           protein_node->get_obj()->shine_dal_pos(),
    //           protein_node->get_obj()->length(),
    //           protein_node->get_obj()->concentration(),
    //           protein_node->get_obj()->rna_list()->size());

    //   protein_node = protein_node->get_next();
    // }
  }
}

void IndividualTest::TearDown(void)
{
  delete indiv1;
  delete indiv2;
  delete indiv3;
  delete indiv4;
}


TEST_F(IndividualTest, TestIndiv1)
{
  // Check that we have the right number of promoters, terminators etc
  // and at the right positions
  // "right" means those values we have computed by hand

  // Check genome size
  EXPECT_EQ(109, indiv1->amount_of_dna());
  EXPECT_EQ(109, indiv1->genetic_unit_seq_length(0));

  // Check RNA list
  list<const Rna*> rna_list = indiv1->rna_list();
  EXPECT_EQ(2, rna_list.size());
  const Rna* rna = rna_list.front();
  EXPECT_EQ(LEADING, rna->strand());
  EXPECT_EQ(4, rna->promoter_pos());
  EXPECT_FLOAT_EQ(0.6, rna->basal_level());
  EXPECT_EQ(50, rna->transcript_length());
  rna = rna_list.back();
  EXPECT_EQ(LEADING, rna->strand());
  EXPECT_EQ(81, rna->promoter_pos());
  EXPECT_FLOAT_EQ(0.8, rna->basal_level());
  EXPECT_EQ(82, rna->transcript_length());

  // Check protein list
  list<Protein*> prot_list = indiv1->protein_list();
  EXPECT_EQ(1, prot_list.size());
  Protein* prot = prot_list.front();
  EXPECT_EQ(LEADING, prot->strand());
  EXPECT_EQ(31, prot->shine_dal_pos());
  EXPECT_EQ(4, prot->length());
  EXPECT_FLOAT_EQ(1.4, prot->concentration());
  EXPECT_EQ(2, prot->rna_list().size());
}

TEST_F(IndividualTest, TestIndiv2)
{
  // Check genome size
  EXPECT_EQ(109, indiv2->amount_of_dna());
  EXPECT_EQ(109, indiv2->genetic_unit_seq_length(0));

  // Check RNA list
  list<const Rna*> rna_list = indiv2->rna_list();
  EXPECT_EQ(2, rna_list.size());
  const Rna* rna = rna_list.front();
  EXPECT_EQ(LAGGING, rna->strand());
  EXPECT_EQ(104, rna->promoter_pos());
  EXPECT_FLOAT_EQ(0.6, rna->basal_level());
  EXPECT_EQ(50, rna->transcript_length());
  rna = rna_list.back();
  EXPECT_EQ(LAGGING, rna->strand());
  EXPECT_EQ(27, rna->promoter_pos());
  EXPECT_FLOAT_EQ(0.8, rna->basal_level());
  EXPECT_EQ(82, rna->transcript_length());

  // Check protein list
  list<Protein*> prot_list = indiv2->protein_list();
  EXPECT_EQ(1, prot_list.size());
  Protein* prot = prot_list.front();
  EXPECT_EQ(LAGGING, prot->strand());
  EXPECT_EQ(77, prot->shine_dal_pos());
  EXPECT_EQ(4, prot->length());
  EXPECT_FLOAT_EQ(1.4, prot->concentration());
  EXPECT_EQ(2, prot->rna_list().size());
}

TEST_F(IndividualTest, TestIndiv3)
{
  // Check genome size
  EXPECT_EQ(81, indiv3->amount_of_dna());
  EXPECT_EQ(81, indiv3->genetic_unit_seq_length(0));

  // Check RNA list
  list<const Rna*> rna_list = indiv3->rna_list();
  EXPECT_EQ(1, rna_list.size());
  const Rna* rna = rna_list.front();
  EXPECT_EQ(LEADING, rna->strand());
  EXPECT_EQ(54, rna->promoter_pos());
  EXPECT_FLOAT_EQ(0.8, rna->basal_level());
  EXPECT_EQ(42, rna->transcript_length());

  // Check protein list
  list<Protein*> prot_list = indiv3->protein_list();
  EXPECT_EQ(1, prot_list.size());
  Protein* prot = prot_list.front();
  EXPECT_EQ(LEADING, prot->strand());
  EXPECT_EQ(4, prot->shine_dal_pos());
  EXPECT_EQ(4, prot->length());
  EXPECT_FLOAT_EQ(0.8, prot->concentration());
  EXPECT_EQ(1, prot->rna_list().size());
}

TEST_F(IndividualTest, TestIndiv4)
{
  // Check genome size
  EXPECT_EQ(81, indiv4->amount_of_dna());
  EXPECT_EQ(81, indiv4->genetic_unit_seq_length(0));

  // Check RNA list
  list<const Rna*> rna_list = indiv4->rna_list();
  EXPECT_EQ(1, rna_list.size());
  const Rna* rna = rna_list.front();
  EXPECT_EQ(LAGGING, rna->strand());
  EXPECT_EQ(26, rna->promoter_pos());
  EXPECT_FLOAT_EQ(0.8, rna->basal_level());
  EXPECT_EQ(42, rna->transcript_length());

  // Check protein list
  list<Protein*> prot_list = indiv4->protein_list();
  EXPECT_EQ(1, prot_list.size());
  Protein* prot = prot_list.front();
  EXPECT_EQ(LAGGING, prot->strand());
  EXPECT_EQ(76, prot->shine_dal_pos());
  EXPECT_EQ(4, prot->length());
  EXPECT_FLOAT_EQ(0.8, prot->concentration());
  EXPECT_EQ(1, prot->rna_list().size());
}

// ===========================================================================
//                                Protected Methods
// ===========================================================================

// ===========================================================================
//                              Non inline accessors
// ===========================================================================
