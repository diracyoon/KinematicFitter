LIBRARY = KinematicFitter
OBJDIR = obj
DEPDIR = $(OBJDIR)/dep
SRCDIR = src
INCDIR = include
ADDITIONAL_ROOTMAPLIBRARY= -rml libPhysicsToolsKinFitter.so -rml libDataFormats.so -rml libAnalyzerTools.so -rml libTMVA.so -rml libonnxruntime.so

include $(SKFlat_WD)/makefile.common

#INCLUDES += -I$(SKFlat_WD)/external/onnxruntime-linux-x64-1.19.2/include
INCLUDES += -I$(SKFlat_WD)/external/onnxruntime/include
INCLUDES += -I$(SKFlat_WD)/external/onnxruntime/include/onnxruntime/core/session
INCLUDES += -I$(SKFlat_WD)/DataFormats/include/
INCLUDES += -I$(SKFlat_WD)/external/KinematicFitter/include
INCLUDES += -I$(SKFlat_WD)/Analyzers/include
INCLUDES += -I/cvmfs/cms.cern.ch/$(SCRAM_ARCH)/cms/$(cmsswrel)/src/

LDFLAGS += -L$(SKFlat_WD)/lib -lDataFormats
LDFLAGS += -L$(SKFlat_WD)/external/onnxruntime/build/Linux/Release/ -lonnxruntime
