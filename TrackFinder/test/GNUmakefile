TBINS := SurfaceTest.x TrackTest.x PropTest.x LATest.x
SIMPLEBINS := $(TBINS)
FLAGS := $(addprefix $(workdir),$(TBINS:.x=.success))

test: tbin $(FLAGS)

$(FLAGS): $(workdir)%.success: $(curdir)%.sh $(bindir)%.x
	echo "<**testing**> $*"
	cd $(workdir); if $< > $*.out 2> $*.err; then touch $@;else echo "Test failed."; fi


########################################################################
include SoftRelTools/arch_spec_root.mk
include SoftRelTools/arch_spec_art.mk
include SoftRelTools/arch_spec_pgsql.mk
include SoftRelTools/standard.mk

override BINLIBS += -lTrackFinder -lGeometry_service -lart_Framework_Services_Registry -lart_Framework_Core -lCintex
