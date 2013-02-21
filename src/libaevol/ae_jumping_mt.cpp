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
#include <string.h>
#include <time.h>



// =================================================================
//                            Project Files
// =================================================================
#include <ae_jumping_mt.h>
#include <ae_jump_poly.h>




//##############################################################################
//                                                                             #
//                             Class ae_jumping_mt                             #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================
int32_t ae_jumping_mt::nb_jumps   = 0;
double  ae_jumping_mt::jump_time  = 0;

// =================================================================
//                             Constructors
// =================================================================
/*!
  Create a generator initialized with a simple uint32_t
 */
ae_jumping_mt::ae_jumping_mt( const uint32_t& simple_seed )
{
  _sfmt = new sfmt_t();
  sfmt_init_gen_rand( _sfmt, simple_seed );
  
  // Jump to get rid of the initializatino skew
  jump();
}

/*!
  Create a copy of an existing generator
 */
ae_jumping_mt::ae_jumping_mt( const ae_jumping_mt& model )
{
  _sfmt = new sfmt_t();
  memcpy(_sfmt->state, model._sfmt->state, SFMT_N*sizeof(_sfmt->state[0]));
  _sfmt->idx = model._sfmt->idx;
}

/*!
  Load a generator from a gz backup file
 */
ae_jumping_mt::ae_jumping_mt( gzFile backup_file )
{
  _sfmt = new sfmt_t();
  gzread( backup_file, _sfmt->state, SFMT_N * sizeof( _sfmt->state[0] ) );
  gzread( backup_file, &(_sfmt->idx), sizeof( _sfmt->idx ) );
}

// =================================================================
//                             Destructors
// =================================================================
ae_jumping_mt::~ae_jumping_mt( void )
{
  delete _sfmt;
}

// =================================================================
//                            Public Methods
// =================================================================
/*!
  Jump Ahead by a predefined jump length
 */
void ae_jumping_mt::jump( void )
{
  //~ clock_t start, end;
  //~ start = clock();
  
  #ifdef TRIVIAL_METHOD_JUMP_SIZE
    for ( int i = 0 ; i < TRIVIAL_METHOD_JUMP_SIZE ; i++ )
    {
      sfmt_genrand_real2( _sfmt );
    }
  #else
    SFMT_jump( _sfmt, jump_poly );
  #endif
    
    
  //~ int nb_sup_steps = 100 * (double)rand() / ((double)RAND_MAX + 1);
  //~ for ( int i = 0 ; i < nb_sup_steps ; i++ )
  //~ {
    //~ sfmt_genrand_real2( _sfmt );
  //~ }
    
  
  
  //~ end = clock();
  //~ jump_time += (end - start) * 1000 / CLOCKS_PER_SEC;
  //~ nb_jumps++;
}

/*!
  Binomial drawing of parameter (nb_drawings, prob).

  Number of successes out of nb_drawings trials each of probability prob.
 */
int32_t ae_jumping_mt::binomial_random( int32_t nb_drawings, double prob )
{
  int32_t nb_success;

  // The binomial distribution is invariant under changing 
  // ProbSuccess to 1-ProbSuccess, if we also change the answer to 
  // NbTrials minus itself; we ll remember to do this below.
  double p;
  if ( prob <= 0.5 ) p = prob;
  else p = 1.0 - prob;

  // mean of the deviate to be produced   
  double mean = nb_drawings * p; 

  
  if ( nb_drawings < 25 )
  // Use the direct method while NbTrials is not too large. 
  // This can require up to 25 calls to the uniform random. 
  {
    nb_success = 0; 
    for ( int32_t j = 1 ; j <= nb_drawings ; j++ ) 
    {
      if ( random() < p ) nb_success++; 
    }
  }
  else if ( mean < 1.0 )
  // If fewer than one event is expected out of 25 or more trials, 
  // then the distribution is quite accurately Poisson. Use direct Poisson method. 
  {
    double g = exp( -mean ); 
    double t = 1.0; 
    int32_t j;
    for ( j = 0; j <= nb_drawings ; j++ ) 
    { 
      t = t * random(); 
      if ( t < g ) break; 
    }

    if ( j <= nb_drawings ) nb_success = j;
    else nb_success = nb_drawings;
  }

  else
  // Use the rejection method.
  {
    double en     = nb_drawings; 
    double oldg   = gammln( en + 1.0 ); 
    double pc     = 1.0 - p; 
    double plog   = log( p ); 
    double pclog  = log( pc ); 
    
    // rejection method with a Lorentzian comparison function.
    double sq = sqrt( 2.0 * mean * pc );  
    double angle, y, em, t;
    do 
    {
      do 
      {
        angle = M_PI * random(); 
        y = tan( angle ); 
        em = sq*y + mean; 
      } while ( em < 0.0 || em >= (en + 1.0) ); // Reject. 

      em = floor( em ); // Trick for integer-valued distribution. 
      t = 1.2 * sq * (1.0 + y*y)
              * exp( oldg - gammln(em + 1.0) - gammln(en - em + 1.0) + em * plog + (en - em) * pclog ); 

    } while ( random() > t ); // Reject. This happens about 1.5 times per deviate, on average. 
    
    nb_success = (int32_t) rint( em ); 
  }
  
  
  // Undo the symmetry transformation. 
  if ( p != prob ) nb_success = nb_drawings - nb_success;

  return nb_success;
}

double ae_jumping_mt::gaussian_random( void )
{
  double x1, x2;
  double r = 0;
  do 
  {
    x1 = 2.0 * random() - 1.0;
    x2 = 2.0 * random() - 1.0;
    
    r = x1*x1 + x2*x2; // (x1,x2) must be in the unit circle
  } 
  while ( (r >= 1.0) || (r == 0) );
  
  r = sqrt( (-2.0 * log(r)) / r ); // Box-muller transformation
  
  return x1 * r;
}

int8_t ae_jumping_mt::roulette_random( double* probs, int8_t nb_elts )
{
  double pick_one = 0.0;
  while ( pick_one == 0 )
  {
    pick_one = random();
  }
  
  int8_t found_org = 0;
  
  pick_one -= probs[0];
  while ( pick_one > 0 )
  {
    assert(found_org<nb_elts-1);
    pick_one -= probs[++found_org];
  }
  return found_org;
}

void ae_jumping_mt::multinomial_drawing( int32_t* destination, double* source, int32_t nb_drawings, int32_t nb_colors )
{
  //    This function generates a vector of random variates, each with the
  //    binomial distribution (source code from http://www.agner.org/random/).
  //    The multinomial distribution is the distribution you get when drawing
  //    balls from an urn with more than two colors, with replacement.

  //    Parameters:
  //    destination:    An output array to receive the number of balls of each 
  //                    color. Must have space for at least 'nb_colors' elements.
  //    source:         An input array containing the probability or fraction 
  //                    of each color in the urn. Must have 'nb_colors' elements.
  //                    All elements must be non-negative. The sum doesn't have
  //                    to be 1, but the sum must be positive.
  //    nb_drawings:    The number of balls drawn from the urn.                   
  //    nb_colors:      The number of possible colors. 
  
  
  if ( nb_drawings < 0 || nb_colors < 0 ) 
  {
    printf( "%s:%d: error: Negative parameter in multinomial function.\n", __FILE__, __LINE__ );
    assert( false );
    exit( EXIT_FAILURE );
  }
  if ( nb_colors == 0 ) return;
  
  // compute sum of probabilities
  double sum = 0.0;
  double p;
  for ( int32_t i = 0 ; i < nb_colors ; i++ ) 
  { 
    p = source[i];
    if ( p < 0 )
    {
      printf( "%s:%d: error: Negative parameter in multinomial function.\n", __FILE__, __LINE__ );
      assert( false );
      exit( EXIT_FAILURE );
    }
    sum += p;
  }
  if ( sum == 0 && nb_drawings > 0 ) 
  {
    printf( "Zero sum in multinomial function\n" );
    assert( false );
    exit( EXIT_FAILURE );
  }

  int32_t x;
  int32_t n = nb_drawings;
  for ( int32_t i = 0 ; i < nb_colors - 1 ; i++ ) 
  { 
    // generate output by calling binomial (nb_colors-1) times
    p = source[i];
    if ( sum <= p ) 
    {
      // this fixes two problems:
      // 1. prevent division by 0 when sum = 0
      // 2. prevent p/sum getting bigger than 1 in case of rounding errors
      x = n;
    }
    else 
    {    
      x = binomial_random( n, p/sum );
    }
    n = n - x; 
    sum = sum - p;
    destination[i] = x;
  }

  // get the last one
  destination[nb_colors-1] = n;
}

void ae_jumping_mt::save( gzFile backup_file ) const
{
  gzwrite( backup_file, _sfmt->state, SFMT_N * sizeof( _sfmt->state[0] ) );
  gzwrite( backup_file, &(_sfmt->idx), sizeof( _sfmt->idx ) );
}


// =================================================================
//                           Protected Methods
// =================================================================
double ae_jumping_mt::gammln( double X )
// Returns the value ln[gamma(X)] for X.
// The gamma function is defined by the integral  gamma(z) = int(0, +inf, t^(z-1).e^(-t)dt). 
// When the argument z is an integer, the gamma function is just the familiar factorial 
// function, but offset by one, n! = gamma(n + 1).
{
  double x, y, tmp, ser; 
  static double cof[6] = {  76.18009172947146,
                            -86.50532032941677,
                            24.01409824083091,
                            -1.231739572450155,
                            0.1208650973866179e-2,
                            -0.5395239384953e-5 }; 

  y = x = X;
  tmp = x + 5.5; 
  tmp -= (x+0.5) * log(tmp);
  ser = 1.000000000190015; 
  
  for ( int8_t j = 0 ; j <= 5 ; j++ )
  {
    ser += cof[j] / ++y; 
  }

  return -tmp + log(2.5066282746310005 * ser / x); 
}

// =================================================================
//                          Non inline accessors
// =================================================================
