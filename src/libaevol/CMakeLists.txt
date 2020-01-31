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
include_directories("./mutation")

# ============================================================================
# Add targets add their dependencies
# ============================================================================
set(aevol_sources
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
#  mutation/GeneMutation.cpp
#  mutation/GeneMutation.h
#  GeneTree.cpp
#  GeneTree.h
#  GeneTreeNode.cpp
#  GeneTreeNode.h
  GeneticUnit.cpp
  GeneticUnit.h
  Habitat.cpp
  Habitat.h
  HabitatFactory.cpp
  HabitatFactory.h
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
  mutation/Mutation.cpp
  mutation/Mutation.h
  Observable.cpp
  Observable.h
  Observer.h
  OutputManager.cpp
  OutputManager.h
  mutation/MutationParams.cpp
  mutation/MutationParams.h
  Phenotype.h
  Protein.cpp
  Protein.h
  #MersenneTwister.cpp
  #MersenneTwister.h
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
  Utils.cpp
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
  AeTime.cpp
  AeTime.h
  World.cpp
  World.h
  ObservableEvent.h
  FuzzyFactory.h
  FuzzyFactory.cpp
  AbstractFuzzy.h
  AbstractFuzzy.cpp
  HybridFuzzy.h
  HybridFuzzy.cpp
  mutation/PointMutation.cpp
  mutation/PointMutation.h
  mutation/SmallInsertion.cpp
  mutation/SmallInsertion.h
  mutation/SmallDeletion.cpp
  mutation/SmallDeletion.h
  mutation/LocalMutation.h
  mutation/Rearrangement.h
  mutation/Duplication.cpp
  mutation/Duplication.h
  mutation/Deletion.cpp
  mutation/Deletion.h
  mutation/Translocation.cpp
  mutation/Translocation.h
  mutation/Inversion.cpp
  mutation/Inversion.h
  mutation/HorizontalTransfer.h
  mutation/InsertionHT.cpp
  mutation/InsertionHT.h
  mutation/ReplacementHT.cpp
  mutation/ReplacementHT.h make_unique.h ae_logger.h
  SIMD_Individual.h
  SIMD_Individual.cpp
  SaveWorld.h
  SaveWorld.cpp
  LightTree.h
  LightTree.cpp
  Dna_SIMD.h
  Dna_SIMD.cpp
  mutation/DnaMutator.h
  mutation/DnaMutator.cpp
  Stats_SIMD.h
  Stats_SIMD.cpp
  mutation/MutationEvent.h
  mutation/MutationEvent.cpp
  AncestorStats.h
  AncestorStats.cpp)
  
if ( ${with-raevol} )
    set(aevol_sources ${aevol_sources} 
          raevol/Habitat_R.h
          raevol/Habitat_R.cpp
          raevol/Individual_R.h
          raevol/Individual_R.cpp
          raevol/Individual_R_X11.h
          raevol/Individual_R_X11.cpp
          raevol/PhenotypicTarget_R.h
          raevol/PhenotypicTarget_R.cpp
          raevol/PhenotypicTargetHandler_R.h
          raevol/PhenotypicTargetHandler_R.cpp
          raevol/Protein_R.h
          raevol/Protein_R.cpp
          raevol/Rna_R.h
          raevol/Rna_R.cpp 
	  #"SIMD_Abstract Metadata.cpp" "SIMD_Abstract Metadata.h" SIMD_Map_Metadata.cpp SIMD_Map_Metadata.h
	  )
endif ()

add_library(aevol ${aevol_sources})

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
