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
include_directories("./phenotype")
include_directories("./biochemistry")
include_directories("./rng")
include_directories("./world")
include_directories("./population")
include_directories("./stats")
include_directories("./cuda")
include_directories("./io")
include_directories("./gui")
include_directories("./utils")
# ============================================================================
# Add targets add their dependencies
# ============================================================================
set(aevol_sources
  aevol.h
  biochemistry/Alignment.cpp
  biochemistry/Alignment.h
  biochemistry/Codon.cpp
  biochemistry/Codon.h
  biochemistry/Dna.cpp
  biochemistry/Dna.h
  io/Dump.cpp
  io/Dump.h
  ae_enums.h
  phenotype/PhenotypicSegment.h
  world/ExpManager.cpp
  world/ExpManager.h
  world/ExpManager_X11.cpp
  world/ExpManager_X11.h
  world/ExpSetup.cpp
  world/ExpSetup.h
  phenotype/Gaussian.h
#  mutation/GeneMutation.cpp
#  mutation/GeneMutation.h
#  io/GeneTree.cpp
#  io/GeneTree.h
#  io/GeneTreeNode.cpp
#  io/GeneTreeNode.h
  biochemistry/GeneticUnit.cpp
  biochemistry/GeneticUnit.h
  world/Habitat.cpp
  world/Habitat.h
  world/HabitatFactory.cpp
  world/HabitatFactory.h
  world/GridCell.cpp
  world/GridCell.h
  rng/JumpingMT.cpp
  rng/JumpingMT.h
  rng/JumpPoly.h
  population/Individual.cpp
  population/Individual.h
  population/Individual_X11.cpp
  population/Individual_X11.h
  population/IndividualFactory.cpp
  population/IndividualFactory.h
  io/Logging.cpp
  io/Logging.h
  mutation/Mutation.cpp
  mutation/Mutation.h
  utils/Observable.cpp
  utils/Observable.h
  utils/Observer.h
  OutputManager.cpp
  OutputManager.h
  mutation/MutationParams.cpp
  mutation/MutationParams.h
  phenotype/Phenotype.h
  biochemistry/Protein.cpp
  biochemistry/Protein.h
  #rng/MersenneTwister.cpp
  #rng/MersenneTwister.h
  population/ReplicationReport.cpp
  population/ReplicationReport.h
  biochemistry/DnaReplicationReport.cpp
  biochemistry/DnaReplicationReport.h
  biochemistry/Rna.cpp
  biochemistry/Rna.h
  population/Selection.cpp
  population/Selection.h
  biochemistry/ae_string.cpp
  biochemistry/ae_string.h
  stats/StatRecord.cpp
  stats/StatRecord.h
  stats/Stats.cpp
  stats/Stats.h
  io/Tree.cpp
  io/Tree.h
  utils/Utils.cpp
  utils/Utils.h
  biochemistry/VisAVis.cpp
  biochemistry/VisAVis.h
  io/ParameterLine.cpp
  io/ParameterLine.h
  gui/X11Window.cpp
  gui/X11Window.h
  phenotype/Fuzzy.cpp
  phenotype/Fuzzy.h
  stats/Metrics.h
  macros.h
  stats/NonCodingMetrics.h
  io/ParamLoader.cpp
  io/ParamLoader.h
  phenotype/PhenotypicTargetHandler.cpp
  phenotype/PhenotypicTargetHandler.h
  phenotype/PhenotypicTarget.cpp
  phenotype/PhenotypicTarget.h
  phenotype/Point.cpp
  phenotype/Point.h
  AeTime.cpp
  AeTime.h
  world/World.cpp
  world/World.h
  utils//ObservableEvent.h
  phenotype/FuzzyFactory.h
  phenotype/FuzzyFactory.cpp
  phenotype/AbstractFuzzy.h
  phenotype/AbstractFuzzy.cpp
  phenotype/HybridFuzzy.h
  phenotype/HybridFuzzy.cpp
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
  mutation/ReplacementHT.h
  utils/make_unique.h
  ae_logger.h
  population/SIMD_Individual.h
  population/SIMD_Individual.cpp
  world/SaveWorld.h
  world/SaveWorld.cpp
  io/LightTree.h
  io/LightTree.cpp
  biochemistry/Dna_SIMD.h
  biochemistry/Dna_SIMD.cpp
  mutation/DnaMutator.h
  mutation/DnaMutator.cpp
  stats/Stats_SIMD.h
  stats/Stats_SIMD.cpp
  mutation/MutationEvent.h
  mutation/MutationEvent.cpp
  stats/AncestorStats.h
  stats/AncestorStats.cpp)
  
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
