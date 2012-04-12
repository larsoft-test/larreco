PACKAGE     := TrackFinderTest
LIB_TYPE    := shared
LIB         := lib$(PACKAGE)
LIBCXXFILES := PropTest.cxx

LIBLINK     := -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -l$(PACKAGE)

TBINS := SurfaceTest.x TrackTest.x LATest.x
SIMPLEBINS := $(TBINS)
FLAGS := $(addprefix $(workdir),$(TBINS:.x=.success)) $(workdir)PropTest.success

test: tbin $(FLAGS)

$(workdir)PropTest.success: $(libdir)libTrackFinderTest.so

$(bindir)PropTest.x:
	touch $@

$(FLAGS): $(workdir)%.success: $(curdir)%.sh $(bindir)%.x
	echo "<**testing**> $*"
	if [ -f $(curdir)$*.fcl ]; then cp $(curdir)$*.fcl $(workdir); fi;
	cd $(workdir); if $< > $*.out 2> $*.err; then touch $@;else echo "Test failed."; fi


########################################################################
include SoftRelTools/arch_spec_root.mk
include SoftRelTools/arch_spec_art.mk
include SoftRelTools/arch_spec_pgsql.mk
include SoftRelTools/standard.mk

override BINLIBS += -lTrackFinder -lGeometry_service -lart_Framework_Services_Registry -lart_Framework_Core -lCintex

override LIBLIBS += $(LOADLIBES) -L$(ART_LIB) -lart_Framework_Services_Optional_TFileService_service -L$(ROOTSYS)/lib -lHist -lRIO -lGeom -lXMLIO -lGdml -lTree -lGraf -lGui -lGed -lGeomBuilder -lGraf3d -lGpad -lGeomPainter -L$(ROOTSYS)/lib -lGeom -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -lTrackFinder -lUtilities -lGeometry
