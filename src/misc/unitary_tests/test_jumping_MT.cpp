//*****************************************************************************
//
//                         aevol - Artificial Evolution
//
// Copyright (C) 2004 LIRIS.
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
#include <time.h>
#include <zlib.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_jumping_mt.h>





int main( int argc, char* argv[] )
{
  // Time coutning stuff
  timespec deb, fin, res;
  clock_t start;
  double exec;
  double elapsed;
  clock_getres(CLOCK_REALTIME, &res);
  
  // Create a reference generator
  sfmt_t an_sfmt;
  sfmt_t* ref = &an_sfmt;
  sfmt_init_gen_rand( ref, 25 );
  
  // Create a wrapped generator
  ae_jumping_mt* my_jumping_mt = new ae_jumping_mt( 25 );
  
  // Generate a few numbers and check
  printf( "******************** Reals in [0, 1)\n" );
  for ( int i = 0 ; i < 5 ; i++ )
  {
    printf( "%f\t\t%f\n", sfmt_genrand_real2( ref ), my_jumping_mt->random() );
  }
  printf( "******************** Integers in [0, 100)\n" );
  for ( int i = 0 ; i < 5 ; i++ )
  {
    printf( "%"PRId32"\t\t", (int32_t)( 100.0 * (((double)sfmt_genrand_uint32( ref )) / MT_RAND_MAX_PLUS_1) ) );
    printf( "%"PRId32"\n", my_jumping_mt->random( 100 ) );
  }
  
  // Jump ahead manually
	start = clock();
  clock_gettime(CLOCK_REALTIME, &deb);
  for ( uint64_t i = 0 ; i < 2000000 ; i++ )
  {
    sfmt_genrand_real2( ref );
  }
  clock_gettime(CLOCK_REALTIME, &fin);
	elapsed = clock() - start;
	elapsed = elapsed * 1000 / CLOCKS_PER_SEC;
  exec = (double) ((fin.tv_nsec - deb.tv_nsec) / 1000);
  printf( " %u:%u %lu:%lu %lu:%lu exec %f microsecond\n", res.tv_sec, res.tv_nsec, deb.tv_sec, deb.tv_nsec, fin.tv_nsec, fin.tv_sec, exec);
  printf( "%f ms\n", elapsed );
  
  // Jump ahead with polynomial method
	start = clock();
  clock_gettime(CLOCK_REALTIME, &deb);
  my_jumping_mt->jump();
  clock_gettime(CLOCK_REALTIME, &fin);
	elapsed = clock() - start;
	elapsed = elapsed * 1000 / CLOCKS_PER_SEC;
  exec = (double) ((fin.tv_nsec - deb.tv_nsec) / 1000);
  printf(" %u:%u %lu:%lu %lu:%lu exec %f microsecond\n", res.tv_sec, res.tv_nsec, deb.tv_sec, deb.tv_nsec, fin.tv_nsec, fin.tv_sec, exec);
  printf( "%f ms\n", elapsed );
  
  
  // Generate a few numbers and check
  printf( "******************** Reals in [0, 1)\n" );
  for ( int i = 0 ; i < 20 ; i++ )
  {
    printf( "%f\t\t%f\n", sfmt_genrand_real2( ref ), my_jumping_mt->random() );
  }
  printf( "******************** Integers in [0, 100)\n" );
  for ( int i = 0 ; i < 20 ; i++ )
  {
    printf( "%"PRId32"\t\t", (int32_t)( 100.0 * (((double)sfmt_genrand_uint32( ref )) / MT_RAND_MAX_PLUS_1) ) );
    printf( "%"PRId32"\n", my_jumping_mt->random( 100 ) );
  }
  
  //~ gzFile* myFile = (gzFile*) gzopen( "prng.ae", "w" );
  //~ my_jumping_mt->save( myFile );
  //~ gzclose( myFile );
  
  //~ // Generate a few numbers and check
  //~ printf( "******************** Reals in [0, 1)\n" );
  //~ for ( int i = 0 ; i < 20 ; i++ )
  //~ {
    //~ printf( "%f\n", my_jumping_mt->random() );
  //~ }
  //~ printf( "******************** Integers in [0, 100)\n" );
  //~ for ( int i = 0 ; i < 20 ; i++ )
  //~ {
    //~ printf( "%"PRId32"\n", my_jumping_mt->random( 100 ) );
  //~ }
  
  //~ gzFile* myOtherFile = (gzFile*) gzopen( "prng.ae", "r" );
  //~ ae_jumping_mt* my_other_jumping_mt = new ae_jumping_mt( myOtherFile );
  //~ gzclose( myOtherFile );
  
  //~ // Generate a few numbers and check
  //~ printf( "******************** Reals in [0, 1)\n" );
  //~ for ( int i = 0 ; i < 20 ; i++ )
  //~ {
    //~ printf( "%f\n", my_other_jumping_mt->random() );
  //~ }
  //~ printf( "******************** Integers in [0, 100)\n" );
  //~ for ( int i = 0 ; i < 20 ; i++ )
  //~ {
    //~ printf( "%"PRId32"\n", my_other_jumping_mt->random( 100 ) );
  //~ }
}








