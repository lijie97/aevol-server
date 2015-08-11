############################################
#                 Variables                #
############################################
AM_CPPFLAGS = $(AEVOLCPPFLAGS)
AM_CXXFLAGS = $(AEVOLCXXFLAGS)
AM_LDFLAGS = $(AEVOLLDFLAGS)

SUBDIRS = SFMT-src-1.4 .

########################
# Set library to build #
########################
noinst_LIBRARIES = libaevol.a

############################################
# C99 exact-width integer specific support #
############################################
# __STDC_FORMAT_MACROS allows us to use exact-width integer format specifiers e.g. PRId32 (for printf etc)
# __STDC_CONSTANT_MACROS allows us to define exact-width integer macros with e.g. INT32_C( <value> )
# __STDC_LIMIT_MACROS allows us to use exact-width integer limit macros with e.g. INT32_MAX
AM_CPPFLAGS += -D__STDC_FORMAT_MACROS -D__STDC_CONSTANT_MACROS -D__STDC_LIMIT_MACROS

common_headers = \
aevol.h \
Alignment.h \
Codon.h \
Dna.h \
Dump.h \
ae_enums.h \
PhenotypicSegment.h \
ExpManager.h \
ExpSetup.h \
Fuzzy.h \
Gaussian.h \
GeneticUnit.h \
GridCell.h \
Individual.h \
IndividualFactory.h \
Logging.h \
macros.h \
Mutation.h \
LocalMutation.h \
PointMutation.h \
SmallInsertion.h \
SmallDeletion.h \
Rearrangement.h \
Duplication.h \
Deletion.h \
Translocation.h \
Inversion.h \
HT.h \
InsertionHT.h \
ReplacementHT.h \
Observable.h \
ObservableEvent.h \
Observer.h \
OutputManager.h \
MutationParams.h \
Phenotype.h \
PhenotypicTarget.h \
PhenotypicTargetHandler.h \
Point.h \
Protein.h \
JumpingMT.h \
JumpPoly.h \
ReplicationReport.h \
DnaReplicationReport.h \
Rna.h \
Stats.h \
StatRecord.h \
NonCodingMetrics.h \
Metrics.h \
ae_string.h \
Tree.h \
Utils.h \
VisAVis.h \
Selection.h \
Habitat.h \
ParamLoader.h \
ParameterLine.h \
AeTime.h \
World.h \
make_unique.h

common_sources = \
Alignment.cpp \
Codon.cpp \
Dna.cpp \
Dump.cpp \
ExpManager.cpp \
ExpSetup.cpp \
GeneticUnit.cpp \
GridCell.cpp \
Individual.cpp \
IndividualFactory.cpp \
Logging.cpp \
Mutation.cpp \
PointMutation.cpp \
SmallInsertion.cpp \
SmallDeletion.cpp \
Duplication.cpp \
Deletion.cpp \
Translocation.cpp \
Inversion.cpp \
InsertionHT.cpp \
ReplacementHT.cpp \
Observable.cpp \
OutputManager.cpp \
MutationParams.cpp \
PhenotypicTarget.cpp \
PhenotypicTargetHandler.cpp \
Point.cpp \
Protein.cpp \
JumpingMT.cpp \
ReplicationReport.cpp \
DnaReplicationReport.cpp \
Rna.cpp \
Stats.cpp \
StatRecord.cpp \
ae_string.cpp \
Tree.cpp \
Utils.cpp \
VisAVis.cpp \
Selection.cpp \
Fuzzy.cpp \
Habitat.cpp \
ParamLoader.cpp \
ParameterLine.cpp \
AeTime.cpp \
World.cpp

X11_headers = \
Individual_X11.h \
ExpManager_X11.h \
X11Window.h

X11_sources = \
Individual_X11.cpp \
ExpManager_X11.cpp \
X11Window.cpp

# The target_CPPFLAGS for each target may look useless, it is here to force automake to create
# prefixed object files so that we don't use an object file built for another set of macrodefinitions
# e.g. We don't want to use an object file created with NO_X when building the X11 lib

if WITH_X
#libaevol_a_HEADERS = $(common_headers) $(X11_headers)
  libaevol_a_SOURCES = $(common_sources) $(X11_sources)
  libaevol_a_SOURCES += $(common_headers) $(X11_headers)
else # WITHOUT_X
#libaevol_a_HEADERS = $(common_headers)
  libaevol_a_SOURCES = $(common_sources)
  libaevol_a_SOURCES += $(common_headers)
endif # X
