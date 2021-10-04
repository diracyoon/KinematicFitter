LIBRARY = KinematicFitter
OBJDIR = obj
DEPDIR = $(OBJDIR)/dep
SRCDIR = src
INCDIR = include
#ADDITIONAL_ROOTMAPLIBRARY= -rml libMatrix.so -rml libPhysics.so -rml libPhysicsToolsKinFitter.so
#ADDITIONAL_ROOTMAPLIBRARY= -rml libMatrix.so -rml libPhysics.so
ADDITIONAL_ROOTMAPLIBRARY= -rml libPhysicsToolsKinFitter.so -rml libDataFormats.so

include $(SKFlat_WD)/makefile.common

INCLUDES += -I$(SKFlat_WD)/DataFormats/include/
INCLUDES += -I$(SKFlat_WD)/external/KinematicFitter/include
INCLUDES += -I/cvmfs/cms.cern.ch/$(SCRAM_ARCH)/cms/$(cmsswrel)/src/
#INCLUDES += -I/cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/lhapdf/6.2.1-fmblme/include/
