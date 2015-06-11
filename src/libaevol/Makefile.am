############################################
#                 Variables                #
############################################
AM_CPPFLAGS = $(AEVOLCPPFLAGS)
AM_CXXFLAGS = $(AEVOLCXXFLAGS)
AM_LDFLAGS = $(AEVOLLDFLAGS)

if WITH_CPPUNIT
  SUBDIRS = SFMT-src-1.4 . tests
else # WITHOUT_CPPUNIT
  SUBDIRS = SFMT-src-1.4 .
endif

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
ae_logs.h \
macros.h \
ae_mutation.h \
ae_output_manager.h \
ae_params_mut.h \
phenotype.h \
phenotypic_target.h \
phenotypic_target_handler.h \
point.h \
Protein.h \
JumpingMT.h \
ae_jump_poly.h \
ae_replication_report.h \
dna_replic_report.h \
Rna.h \
ae_stats.h \
ae_stat_record.h \
ae_string.h \
ae_tree.h \
ae_utils.h \
ae_vis_a_vis.h \
ae_selection.h \
GeneTree.h \
GeneTreeNode.h \
GeneMutation.h \
habitat.h \
param_loader.h \
f_line.h \
Time.h \
world.h

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
ae_logs.cpp \
ae_mutation.cpp \
ae_output_manager.cpp \
ae_params_mut.cpp \
phenotypic_target.cpp \
phenotypic_target_handler.cpp \
point.cpp \
Protein.cpp \
JumpingMT.cpp \
ae_replication_report.cpp \
dna_replic_report.cpp \
Rna.cpp \
ae_stats.cpp \
ae_stat_record.cpp \
ae_string.cpp \
ae_tree.cpp \
ae_vis_a_vis.cpp \
ae_selection.cpp \
GeneTree.cpp \
GeneTreeNode.cpp \
GeneMutation.cpp \
Fuzzy.cpp \
habitat.cpp \
param_loader.cpp \
f_line.cpp \
Time.cpp \
world.cpp

X11_headers = \
Individual_X11.h \
ExpManager_X11.h \
ae_X11_window.h

X11_sources = \
Individual_X11.cpp \
ExpManager_X11.cpp \
ae_X11_window.cpp

# The target_CPPFLAGS for each target may look useless, it is here to force automake to create
# prefixed object files so that we don't use an object file built for another set of macrodefinitions
# e.g. We don't want to use an object file created with NO_X when building the X11 lib

libaevol_a_CPPFLAGS = $(AM_CPPFLAGS)

if WITH_X
#libaevol_a_HEADERS = $(common_headers) $(X11_headers)
  libaevol_a_SOURCES = $(common_sources) $(X11_sources)
  libaevol_a_SOURCES += $(common_headers) $(X11_headers)
else # WITHOUT_X
#libaevol_a_HEADERS = $(common_headers)
  libaevol_a_SOURCES = $(common_sources)
  libaevol_a_SOURCES += $(common_headers)
endif # X
