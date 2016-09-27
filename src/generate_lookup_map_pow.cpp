//
// Created by arrouan on 13/09/16.
//
#include <cstdlib>
#include <cstdio>
#include <getopt.h>

#include "raevol/Rna_R.h"


using namespace aevol;

int main(int argc, char* argv[])
{
  double hill_shape_theta;
  double hill_shape_n;

  // 2) Define allowed options
  const char * options_list = "n:t:";
  static struct option long_options_list[] = {
      { "hill_shape_n",     required_argument,  NULL, 'n' },
      { "hill_shape_theta",      required_argument,  NULL, 't' },
      { 0, 0 }
  };


  // 3) Get actual values of the command-line options
  int option;
  while ((option = getopt_long(argc, argv, options_list, long_options_list, NULL)) != -1)
  {
    switch (option)
    {
      case 'n' :
      {
        hill_shape_n = atof(optarg);
        break;
      }
      case 't' :
      {
        hill_shape_theta = atof(optarg);
        break;
      }
      default :
      {
        // An error message is printed in getopt_long, we just need to exit
        exit(EXIT_FAILURE);
      }
    }
  }

  printf("Theta %f N %f\n",hill_shape_theta,hill_shape_n);
#ifdef __PROXY_POW_APPROX
  Rna_R::create_lookup_table(hill_shape_theta,hill_shape_n);
#endif
}
