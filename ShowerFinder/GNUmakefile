#
# $Id: GNUmakefile,v 1.2 2008/07/14 19:18:07 brebel Exp $
#
# rahaman@fnal.gov
#
#
include SoftRelTools/arch_spec_root.mk

PACKAGE     := ShowerFinder
LIB_TYPE    := shared
LIB         := lib$(PACKAGE)
LIBCXXFILES := $(wildcard *.cxx)
JOBFILES    := $(wildcard *.py)

LIBLINK    := $(LOADLIBES) -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -lRawData -l$(PACKAGE)

########################################################################
include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_art.mk
