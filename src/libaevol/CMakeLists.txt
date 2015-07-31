# ============================================================================
# Tell cmake about subdirectories to look into
# ============================================================================
add_subdirectory(SFMT-src-1.4)


# ============================================================================
# Find packages
# ============================================================================
find_package(ZLIB REQUIRED)


# ============================================================================
# Add PWD to include directories
# ============================================================================
include_directories(".")

# ============================================================================
# Add targets add their dependencies
# ============================================================================
add_library(aevol
  aevol.h
  Alignment.cpp
  Alignment.h
  Codon.cpp
  Codon.h
  Dna.cpp
  Dna.h
  Dump.cpp
  Dump.h
  ae_enums.h
  PhenotypicSegment.h
  ExpManager.cpp
  ExpManager.h
  ExpManager_X11.cpp
  ExpManager_X11.h
  ExpSetup.cpp
  ExpSetup.h
  Gaussian.h
  GeneMutation.cpp
  GeneMutation.h
  GeneTree.cpp
  GeneTree.h
  GeneTreeNode.cpp
  GeneTreeNode.h
  GeneticUnit.cpp
  GeneticUnit.h
  Habitat.cpp
  Habitat.h
  GridCell.cpp
  GridCell.h
  Individual.cpp
  Individual.h
  Individual_X11.cpp
  Individual_X11.h
  IndividualFactory.cpp
  IndividualFactory.h
  JumpingMT.cpp
  JumpingMT.h
  JumpPoly.h
  Logging.cpp
  Logging.h
  Mutation.cpp
  Mutation.h
  Observable.cpp
  Observable.h
  Observer.h
  OutputManager.cpp
  OutputManager.h
  MutationParams.cpp
  MutationParams.h
  Phenotype.h
  Protein.cpp
  Protein.h
  MersenneTwister.cpp
  MersenneTwister.h
  ReplicationReport.cpp
  ReplicationReport.h
  DnaReplicationReport.cpp
  DnaReplicationReport.h
  Rna.cpp
  Rna.h
  Selection.cpp
  Selection.h
  StatRecord.cpp
  StatRecord.h
  Stats.cpp
  Stats.h
  ae_string.cpp
  ae_string.h
  Tree.cpp
  Tree.h
  Utils.h
  VisAVis.cpp
  VisAVis.h
  X11Window.cpp
  X11Window.h
  ParameterLine.cpp
  ParameterLine.h
  Fuzzy.cpp
  Fuzzy.h
  Metrics.h
  macros.h
  NonCodingMetrics.h
  ParamLoader.cpp
  ParamLoader.h
  PhenotypicTargetHandler.cpp
  PhenotypicTargetHandler.h
  PhenotypicTarget.cpp
  PhenotypicTarget.h
  Point.cpp
  Point.h
  Time.cpp
  Time.h
  World.cpp
  World.h
  ObservableEvent.h FuzzyFactory.h FuzzyFactory.cpp AbstractFuzzy.h AbstractFuzzy.cpp HybridFuzzy.h HybridFuzzy.cpp)


# ============================================================================
# Use GSL and zlib
# ============================================================================
target_link_libraries(aevol PUBLIC ${ZLIB_LIBRARY})
target_link_libraries(aevol PUBLIC sfmt)


# ============================================================================
# Make STDC MACROS available (for fixed width integers)
# ============================================================================
target_compile_definitions(aevol
  PUBLIC __STDC_FORMAT_MACROS
  PUBLIC __STDC_CONSTANT_MACROS)


# ============================================================================
# Use C++11
# ============================================================================
target_compile_options(aevol PUBLIC "-std=c++14")
