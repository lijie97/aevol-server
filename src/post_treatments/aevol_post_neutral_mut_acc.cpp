//
// Created by duazel on 16/04/2020.
//
#include "IOJson.h"
#include "Individual.h"
#include "neutral_mutation_exp.h"
#include "neutral_mutation_output.h"

#include <iostream>
#include <cstdio>
#include <cstring>
#include <getopt.h>
#include <string>
#include <vector>

bool verbose = false;
uint32_t seed_prng = 0;
std::string inputFile;
std::string outputFile;
unsigned int number_generation;

void print_help(char *prog_path) {
  // Get the program file-name in prog_name (strip prog_path of the path)
  char *prog_name; // No new, it will point to somewhere inside prog_path
  if ((prog_name = strrchr(prog_path, '/'))) {
    prog_name++;
  } else {
    prog_name = prog_path;
  }

  printf("******************************************************************************\n");
  printf("*                                                                            *\n");
  printf("*                        aevol - Artificial Evolution                        *\n");
  printf("*                                                                            *\n");
  printf("* Aevol is a simulation platform that allows one to let populations of       *\n");
  printf("* digital organisms evolve in different conditions and study experimentally  *\n");
  printf("* the mechanisms responsible for the structuration of the genome and the     *\n");
  printf("* transcriptome.                                                             *\n");
  printf("*                                                                            *\n");
  printf("******************************************************************************\n");
  printf("\n");
  printf("%s: Accumulate neutral mutations on an individual.\n",
         prog_name);
  printf("\n");
  printf("Usage : %s -h or --help\n", prog_name);
  printf("   or : %s -V or --version\n", prog_name);
  printf("\nOptions\n");
  printf("  -n, --nb_generations\n\tnumber of generations\n");
  printf("  -h, --help\n\tprint this help, then exit\n");
  printf("  -V, --version\n\tprint version number, then exit\n");
}

void interpret_cmd_line_options(int argc, char **argv) {
  const char *short_options = "hVv";
  static struct option long_options[] = {{"help", no_argument, nullptr, 'h'},
                                         {"version", no_argument, nullptr, 'V'},
                                         {"verbose", no_argument, nullptr, 'v'},
                                         {"nb_generations", required_argument, nullptr, 'n'},
                                         {"seed", required_argument, nullptr, 's'},
                                         {"file", required_argument, nullptr, 'f'},
                                         {"output_file", required_argument, nullptr, 'o'},
                                         {0, 0, 0, 0}};

  int option;
  while ((option = getopt_long(argc, argv, short_options, long_options,
                               nullptr)) != -1) {
    switch (option) {
    case 'h' :
      print_help(argv[0]);
      exit(EXIT_SUCCESS);
    case 'V' :
      Utils::PrintAevolVersion();
      exit(EXIT_SUCCESS);
    case 'v':
      verbose = true;
      break;
    case 'n':
      number_generation = atol(optarg);
      break;
    case 's':
      seed_prng = atol(optarg);
      break;
    case 'f':
      inputFile.assign(optarg);
      break;
    case 'o':
      outputFile.assign(optarg);
      break;
    default:
      // An error message is printed in getopt_long, we just need to exit
      exit(EXIT_FAILURE);
    }
  }

  if (number_generation == 0) {
    Utils::ExitWithUsrMsg("The parameter number_generation is needed");
  }
  if (seed_prng == 0) {
    Utils::ExitWithUsrMsg("The parameter seed is needed");
  }

}

int main(int argc, char ** argv) {
  inputFile = "input.json";
  outputFile = "output.json";
  interpret_cmd_line_options(argc, argv);

  std::cout << argv[0] << " started with :" << std::endl;
  std::cout << " - number_generation = " << number_generation << std::endl;
  std::cout << " - seed = " << seed_prng << std::endl;
  std::cout << " - input file = " << inputFile << std::endl;
  std::cout << " - output file = " << outputFile << std::endl;
  IOJson inputJson(inputFile);

  out::init("result.txt", "mutation.txt");

  run_generations(number_generation, inputJson.individuals()[0]);
//  std::vector<Individual> indiv_vector;
//  indiv_vector.emplace_back(*indiv);

//  inputJson.individuals().clear();
//  inputJson.setIndividuals(indiv_vector);

//  inputJson.write(outputFile);
  return 0;
}
