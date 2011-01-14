#
# $Id: GNUmakefile,v 1.1.1.1 2009/05/13 01:38:51 katori Exp $
#
# kinga,josh,echurch@fnal.gov
#
#
include SoftRelTools/arch_spec_root.mk

PACKAGE     := TrackFinder
LIB_TYPE    := shared
LIB         := lib$(PACKAGE)
LIBCXXFILES := $(wildcard *.cxx)
JOBFILES    := $(wildcard *.py)

LIBLINK    := $(LOADLIBES) -L/grid/fermiapp/uboone/code/externals-1.1.0/root/v5_18_00f/Linux64bit-2-6-2-5-GCC_4_1_2/lib/Linux2.6-GCC/ -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -lGeometry -lSimulation -lRawData -lFilters -lRecoBase -lGenfit -lVMC -l$(PACKAGE)

########################################################################
include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_art.mk
