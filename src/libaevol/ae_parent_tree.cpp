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


/** \class
 *  \brief
 */


// =================================================================
//                              Libraries
// =================================================================
#include <stdio.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_parent_tree.h>
#include <ae_experiment.h>
#include <ae_population.h>
#include <ae_individual.h>
#include <ae_replication_report.h>
#include <ae_genetic_unit.h>
#ifdef __REGUL
  #include <ae_influence_R.h>
  #include <ae_protein_R.h>
#endif





//##############################################################################
//                                                                             #
//                                Class ae_parent_tree                         #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================
int32_t ae_parent_tree::NO_PARENT = -1;

// =================================================================
//                             Constructors
// =================================================================

// =================================================================
//                             Destructors
// =================================================================

// =================================================================
//                            Public Methods
// =================================================================

/*
// ONLY USE FOR DEBUG PURPOSE
void ae_parent_tree::_dump( void ) {
  FILE *f = fopen("parent_tree.dat", "w+");
  int32_t *j = new int32_t;
  int32_t n = ae_common :: sim -> get_num_gener(), r;
  printf("++ DUMP(%d)\n", n);
  for(int32_t a = 0; a < ae_common :: init_pop_size; a++) {
    for(int32_t b = 0; b < ae_common :: init_pop_size; b++) {
      r = get_LCA(n, a, b, j);
      fprintf(f, "rel(%04d, %04d) = %03d (%04d)\n", a, b, r, *j);
    }
  }
  printf("-- DUMP(%d)\n", n);
  fclose(f);
}
*/

void ae_parent_tree::insert_current_generation( ) {
  int32_t n = ae_common::sim->get_num_gener();
  int32_t i = NO_PARENT; // indiv index
  int32_t j = NO_PARENT; // parent index
  ae_list_node* indiv_node = ae_common :: sim -> get_pop() -> get_indivs() -> get_first();
  ae_individual* indiv = NULL;
  ae_replication_report* report = NULL;
  printf("inserting (%d)...\n", n);
  while( indiv_node != NULL ) {

    indiv = (ae_individual*) indiv_node -> get_obj();
    report = indiv -> get_replic_report();
    i = report -> get_index() ;
    j = report -> get_parent_index() ;

    parent[n][i] = j;
    indiv_node = indiv_node -> get_next();
  }
  printf("inserted (%d).\n", n);
  //~ _dump(); // Commented because it would not compile (declaration commented and marked as "DEBUG ONLY")
  
}

// =================================================================
//                           Protected Methods
// =================================================================
