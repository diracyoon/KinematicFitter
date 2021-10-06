<<<<<<< HEAD
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
=======
all: DataFormats AnalyzerTools GEScaleSyst KinematicFitter Analyzers Archive

DataFormats::
	(cd DataFormats; make)
	(mvexist.sh DataFormats/src/DataFormats_Dict_rdict.pcm lib/)
	(mvexist.sh DataFormats/libDataFormats.rootmap lib/)

AnalyzerTools::
	(cd AnalyzerTools; make)
	(mvexist.sh AnalyzerTools/src/AnalyzerTools_Dict_rdict.pcm lib/)
	(mvexist.sh AnalyzerTools/libAnalyzerTools.rootmap lib/)

GEScaleSyst::
	(cd external/GEScaleSyst; make)
	(mvexist.sh external/GEScaleSyst/GEScaleSyst_Dict_rdict.pcm lib/)
	(mvexist.sh external/GEScaleSyst/libGEScaleSyst.rootmap lib/)

Analyzers::
	(cd Analyzers; make)
	(mvexist.sh Analyzers/src/Analyzers_Dict_rdict.pcm lib/)
	(mvexist.sh Analyzers/libAnalyzers.rootmap lib/)

KinematicFitter::
	(cd external/KinematicFitter; make)
	(mvexist.sh external/KinematicFitter/src/KinematicFitter_Dict_rdict.pcm lib/)
	(mvexist.sh external/KinematicFitter/libKinematicFitter.rootmap lib/)

Archive::
	(tar -zcf lib/DataFormats.tar.gz DataFormats)
	(tar -zcf lib/AnalyzerTools.tar.gz AnalyzerTools)
	(tar -zcf lib/GEScaleSyst.tar.gz external/GEScaleSyst/GEScaleSyst.*)
	(tar -zcf lib/Analyzers.tar.gz Analyzers)

clean::
	(cd DataFormats; make clean)
	(cd AnalyzerTools; make clean)
	(cd external/GEScaleSyst; make clean)
	(cd external/KinematicFitter; make clean)	
	(cd Analyzers; make clean)

distclean::
	(cd DataFormats; make distclean)
	(cd AnalyzerTools; make distclean)
	(cd external/GEScaleSyst; make distclean)
	(cd external/KinematicFitter; make distclean)	
	(cd Analyzers; make distclean)

LibTarFile = tar/lib.tar.gz
$(LibTarFile): $(wildcard ./lib/*)
	tar -czf $@ ./lib/*

DataTarFile = tar/data.tar.gz
DataFiles = $(shell find data/$(SKFlatV)/ -type f -name '*')
$(DataTarFile): $(DataFiles)
	tar --exclude=data/$(SKFlatV)/Sample -czf $@ ./data/$(SKFlatV)/

CondorTar: $(LibTarFile) $(DataTarFile)
>>>>>>> c8f36c74f4592948a5bc4fad405ef028f76f13fa
