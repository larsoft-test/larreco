#
# $Id: GNUmakefile,v 1.1.1.1 2009/05/13 01:38:51 katori Exp $
#
# brebel@fnal.gov
#
#
include SoftRelTools/arch_spec_root.mk
include SoftRelTools/arch_spec_art.mk

PACKAGE     := HitFinder
LIB_TYPE    := shared
LIB         := lib$(PACKAGE)
LIBCXXFILES := $(wildcard *.cxx)
JOBFILES    := $(wildcard *.py)

LIBLINK    := $(LOADLIBES) -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -lGeometry -l$(PACKAGE)

########################################################################
include SoftRelTools/standard.mk
