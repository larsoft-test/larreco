# Configuration file for CalData
#
# $Id: ex02.py,v 1.5 2010/05/18 21:24:21 kutschke Exp $
#
# Original author E. Church
#
# Spacing is not signficant in this file.

# Define the default configuration for the framework.
import FWCore.ParameterSet.python.Config as caldata

# Give this job a name.  
process = caldata.Process("CalData")

# Maximum number of events to do.
process.maxEvents = caldata.untracked.PSet(
    input = caldata.untracked.int32(9) # See if this works to run fewer than are in input file.
)

# Load the standard message logger configuration.
# Threshold=Info. Limit of 5 per category; then exponential backoff.
#process.load("Config/MessageLogger_cfi")

# Load the service that manages root files for histograms.
process.TFileService = caldata.Service(
    "TFileService",
    fileName = caldata.string("caldata_hist.root"),
    closeFileFast = caldata.untracked.bool(False)
)

process.Timing = caldata.Service("Timing");

# Define the geometry.
process.Geometry = caldata.Service(
    "Geometry",
    SurfaceY = caldata.double(130.0e2), # in cm
    Name     = caldata.string("argoneut"),
    GDML     = caldata.FileInPath("Geometry/gdml/argoneut.gdml")
)


# Define the FFT service
process.LArFFT = caldata.Service(
    "LArFFT",
    FFTSize   = caldata.int32(4096),#double the number of samples for ArgoNeuT to deal with long exponential tail
    FFTOption = caldata.string("ES")
)

# Service to get my MC events, which were run up through DetSim.
process.source = caldata.Source("PoolSource",
                                fileNames = caldata.untracked.vstring("out/genie_gen.root")
                                )

process.caldataCal = caldata.EDProducer(
    "CalWire",
    DetSimModuleLabel  = caldata.string("wiresim"),
    ResponseFile       = caldata.string("/grid/fermiapp/lbne/lar/aux/ArgoResponse1.2.root"),
    ExponentialEndBins = caldata.int32(300),
    UseRawData         = caldata.int32(0)
    )

process.ffthit = caldata.EDProducer(
    "FFTHitFinder",
    CalDataModuleLabel   = caldata.string("caldataCal"),
    MinSigInd       = caldata.double(8.0),
    MinSigCol       = caldata.double(11.0),
    IndWidth        = caldata.double(5.0),
    ColWidth        = caldata.double(7.5),
    Drift           = caldata.double(0.03069),
    # TPC's drift Velocity in units cm/(ADC Sample Time)
    POffset         = caldata.double(20.0),
    OOffset         = caldata.double(24.0),
    MaxMultiHit     = caldata.int32(3)
    )


# Write the events to the output file.
process.output = caldata.OutputModule(
    "PoolOutputModule",
    fileName = caldata.untracked.string('file:calhitf_gen.root'),
)

####### End of the section that defines and configures modules.#########

# Tell the system to execute all paths. Services, source, output are implied ....
process.doit = caldata.EndPath( process.caldataCal*process.ffthit*process.output )

