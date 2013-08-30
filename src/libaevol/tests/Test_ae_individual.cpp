//*****************************************************************************
//
//                         aevol - Artificial Evolution
//
// Copyright (C) 2004  LIRIS.
// Web: https://liris.cnrs.fr/
// E-mail: carole.knibbe@liris.cnrs.fr
// Original Authors : Guillaume Beslon, Carole Knibbe, Virginie Lefort
//                    David Parsons
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//*****************************************************************************
 
 
 
 
 
// ===========================================================================
//                               Include Libraries
// ===========================================================================
#include <string.h>



// ===========================================================================
//                             Include Project Files
// ===========================================================================
#include "Test_ae_individual.h"
#include <ae_macros.h>
#include <ae_genetic_unit.h>
#include <ae_rna.h>
#include <ae_protein.h>
#include <ae_params_mut.h>



// ===========================================================================
//                             Declare Used Namespaces
// ===========================================================================




//############################################################################
//                                                                           #
//                         Class Test_ae_individual                          #
//                                                                           #
//############################################################################
CPPUNIT_TEST_SUITE_REGISTRATION( Test_ae_individual );

// ===========================================================================
//                               Static attributes
// ===========================================================================

// ===========================================================================
//                                  Constructors
// ===========================================================================
Test_ae_individual::Test_ae_individual( void )
{
}

// ===========================================================================
//                                  Destructors
// ===========================================================================
Test_ae_individual::~Test_ae_individual( void )
{
}

// ===========================================================================
//                                   Operators
// ===========================================================================

// ===========================================================================
//                                 Public Methods
// ===========================================================================
void Test_ae_individual::setUp( void )
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
  indiv1 = new ae_individual(NULL, NULL, NULL, &params_mut, 1.0, 10, 1000, false, 0, 1, 0);
  indiv1->add_GU(genome, strlen(genome));
  genome = NULL;

  // Do transcription and translation
  indiv1->do_transcription();
  indiv1->do_translation();



  // Build indiv2
  genome = indiv1->get_genetic_unit(0)->get_dna()->get_subsequence(0,0,LAGGING);
  indiv2 = new ae_individual(NULL, NULL, NULL, &params_mut, 1.0, 10, 1000, false, 0, 1, 0);
  indiv2->add_GU(genome, strlen(genome));
  genome = NULL;

  // Do transcription and translation
  indiv2->do_transcription();
  indiv2->do_translation();

  


  // Build indiv3
  genome = new char[1024];
  sprintf( genome, "%s%s%s%s%s%s%s", as[0], gene, as[1], term, as[2], prom[1], as[3]);
  indiv3 = new ae_individual(NULL, NULL, NULL, &params_mut, 1.0, 10, 1000, false, 0, 1, 0);
  indiv3->add_GU(genome, strlen(genome));
  genome = NULL;

  // Do transcription and translation
  indiv3->do_transcription();
  indiv3->do_translation();




  // Build indiv4
  genome = indiv3->get_genetic_unit(0)->get_dna()->get_subsequence(0,0,LAGGING);
  indiv4 = new ae_individual(NULL, NULL, NULL, &params_mut, 1.0, 10, 1000, false, 0, 1, 0);
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

void Test_ae_individual::tearDown( void )
{
  delete indiv1;
  delete indiv2;
  delete indiv3;
  delete indiv4;
}

void Test_ae_individual::test1( void )
{
  // Check genome size
  CPPUNIT_ASSERT( indiv1->get_amount_of_dna() == 109 );
  CPPUNIT_ASSERT( indiv1->get_genetic_unit_seq_length(0) == 109 );

  // Check RNA list
  ae_list<ae_rna*>* rna_list = indiv1->get_rna_list();
  CPPUNIT_ASSERT( rna_list->get_nb_elts() == 2 );
  ae_rna* rna = rna_list->get_first()->get_obj();
  CPPUNIT_ASSERT( rna->get_strand() == LEADING );
  CPPUNIT_ASSERT( rna->get_promoter_pos() == 4 );
  CPPUNIT_ASSERT( rna->get_basal_level() == 0.6 );
  CPPUNIT_ASSERT( rna->get_transcript_length() == 50 );
  rna = rna_list->get_last()->get_obj();
  CPPUNIT_ASSERT( rna->get_strand() == LEADING );
  CPPUNIT_ASSERT( rna->get_promoter_pos() == 81 );
  CPPUNIT_ASSERT( rna->get_basal_level() == 0.8 );
  CPPUNIT_ASSERT( rna->get_transcript_length() == 82 );

  // Check protein list
  ae_list<ae_protein*>* prot_list = indiv1->get_protein_list();
  CPPUNIT_ASSERT( prot_list->get_nb_elts() == 1 );
  ae_protein* prot = prot_list->get_first()->get_obj();
  CPPUNIT_ASSERT( prot->get_strand() == LEADING );
  CPPUNIT_ASSERT( prot->get_shine_dal_pos() == 31 );
  CPPUNIT_ASSERT( prot->get_length() == 4 );
  CPPUNIT_ASSERT( prot->get_concentration() == 1.4 );
  CPPUNIT_ASSERT( prot->get_rna_list()->get_nb_elts() == 2 );
}

void Test_ae_individual::test2( void )
{
  // Check genome size
  CPPUNIT_ASSERT( indiv2->get_amount_of_dna() == 109 );
  CPPUNIT_ASSERT( indiv2->get_genetic_unit_seq_length(0) == 109 );

  // Check RNA list
  ae_list<ae_rna*>* rna_list = indiv2->get_rna_list();
  CPPUNIT_ASSERT( rna_list->get_nb_elts() == 2 );
  ae_rna* rna = rna_list->get_first()->get_obj();
  CPPUNIT_ASSERT( rna->get_strand() == LAGGING );
  CPPUNIT_ASSERT( rna->get_promoter_pos() == 104 );
  CPPUNIT_ASSERT( rna->get_basal_level() == 0.6 );
  CPPUNIT_ASSERT( rna->get_transcript_length() == 50 );
  rna = rna_list->get_last()->get_obj();
  CPPUNIT_ASSERT( rna->get_strand() == LAGGING );
  CPPUNIT_ASSERT( rna->get_promoter_pos() == 27 );
  CPPUNIT_ASSERT( rna->get_basal_level() == 0.8 );
  CPPUNIT_ASSERT( rna->get_transcript_length() == 82 );

  // Check protein list
  ae_list<ae_protein*>* prot_list = indiv2->get_protein_list();
  CPPUNIT_ASSERT( prot_list->get_nb_elts() == 1 );
  ae_protein* prot = prot_list->get_first()->get_obj();
  CPPUNIT_ASSERT( prot->get_strand() == LAGGING );
  CPPUNIT_ASSERT( prot->get_shine_dal_pos() == 77 );
  CPPUNIT_ASSERT( prot->get_length() == 4 );
  CPPUNIT_ASSERT( prot->get_concentration() == 1.4 );
  CPPUNIT_ASSERT( prot->get_rna_list()->get_nb_elts() == 2 );
}

void Test_ae_individual::test3( void )
{
  // Check genome size
  CPPUNIT_ASSERT( indiv3->get_amount_of_dna() == 81 );
  CPPUNIT_ASSERT( indiv3->get_genetic_unit_seq_length(0) == 81 );

  // Check RNA list
  ae_list<ae_rna*>* rna_list = indiv3->get_rna_list();
  CPPUNIT_ASSERT( rna_list->get_nb_elts() == 1 );
  ae_rna* rna = rna_list->get_first()->get_obj();
  CPPUNIT_ASSERT( rna->get_strand() == LEADING );
  CPPUNIT_ASSERT( rna->get_promoter_pos() == 54 );
  CPPUNIT_ASSERT( rna->get_basal_level() == 0.8 );
  CPPUNIT_ASSERT( rna->get_transcript_length() == 42 );

  // Check protein list
  ae_list<ae_protein*>* prot_list = indiv3->get_protein_list();
  CPPUNIT_ASSERT( prot_list->get_nb_elts() == 1 );
  ae_protein* prot = prot_list->get_first()->get_obj();
  CPPUNIT_ASSERT( prot->get_strand() == LEADING );
  CPPUNIT_ASSERT( prot->get_shine_dal_pos() == 4 );
  CPPUNIT_ASSERT( prot->get_length() == 4 );
  CPPUNIT_ASSERT( prot->get_concentration() == 0.8 );
  CPPUNIT_ASSERT( prot->get_rna_list()->get_nb_elts() == 1 );
}

void Test_ae_individual::test4( void )
{
  // Check genome size
  CPPUNIT_ASSERT( indiv4->get_amount_of_dna() == 81 );
  CPPUNIT_ASSERT( indiv4->get_genetic_unit_seq_length(0) == 81 );

  // Check RNA list
  ae_list<ae_rna*>* rna_list = indiv4->get_rna_list();
  CPPUNIT_ASSERT( rna_list->get_nb_elts() == 1 );
  ae_rna* rna = rna_list->get_first()->get_obj();
  CPPUNIT_ASSERT( rna->get_strand() == LAGGING );
  CPPUNIT_ASSERT( rna->get_promoter_pos() == 26 );
  CPPUNIT_ASSERT( rna->get_basal_level() == 0.8 );
  CPPUNIT_ASSERT( rna->get_transcript_length() == 42 );

  // Check protein list
  ae_list<ae_protein*>* prot_list = indiv4->get_protein_list();
  CPPUNIT_ASSERT( prot_list->get_nb_elts() == 1 );
  ae_protein* prot = prot_list->get_first()->get_obj();
  CPPUNIT_ASSERT( prot->get_strand() == LAGGING );
  CPPUNIT_ASSERT( prot->get_shine_dal_pos() == 76 );
  CPPUNIT_ASSERT( prot->get_length() == 4 );
  CPPUNIT_ASSERT( prot->get_concentration() == 0.8 );
  CPPUNIT_ASSERT( prot->get_rna_list()->get_nb_elts() == 1 );
}

// ===========================================================================
//                                Protected Methods
// ===========================================================================

// ===========================================================================
//                              Non inline accessors
// ===========================================================================
