# Configuration file for Emptyfilter
#
# $Id: ex02.py,v 1.5 2010/05/18 21:24:21 kutschke Exp $
#
# Original author E. Church
#
# Spacing is not signficant in this file.

# Define the default configuration for the framework.
import FWCore.ParameterSet.python.Config as emptyfilter

# Give this job a name.  
process = emptyfilter.Process("Filters")

# Maximum number of events to do.
process.maxEvents = emptyfilter.untracked.PSet(
    input = emptyfilter.untracked.int32(9) # See if this works to run fewer than are in input file.
)

# Load the standard message logger configuration.
# Threshold=Info. Limit of 5 per category; then exponential backoff.
#process.load("Config/MessageLogger_cfi")

# Load the service that manages root files for histograms.
process.TFileService = emptyfilter.Service(
    "TFileService",
    fileName = emptyfilter.string("emptyfilter_hist.root"),
    closeFileFast = emptyfilter.untracked.bool(False)
)

process.Timing = emptyfilter.Service("Timing");

# Define the geometry.
process.Geometry = emptyfilter.Service(
    "Geometry",
    SurfaceY = emptyfilter.double(130.0e2), # in cm
    Name     = emptyfilter.string("argoneut"),
    GDML     = emptyfilter.FileInPath("Geometry/gdml/argoneut.gdml")
)


# Define the FFT service
process.LArFFT = emptyfilter.Service(
    "LArFFT",
    FFTSize   = emptyfilter.int32(4096),#double the number of samples for ArgoNeuT to deal with long exponential tail
    FFTOption = emptyfilter.string("ES")
)

# Service to get my MC events, which were run up through DetSim.
process.source = emptyfilter.Source("PoolSource",
                                fileNames = emptyfilter.untracked.vstring("single_gen.root")
                                )


process.caldataCal = emptyfilter.EDProducer(
    "CalWire",
    DigitModuleLabel  = emptyfilter.string("wiresim"),
    ResponseFile       = emptyfilter.string("/grid/fermiapp/lbne/lar/aux/ArgoResponse1.2.root"),
    ExponentialEndBins = emptyfilter.int32(300),
    UseRawData         = emptyfilter.int32(0)
    )

process.ffthit = emptyfilter.EDProducer(
    "FFTHitFinder",
    CalDataModuleLabel   = emptyfilter.string("caldataCal"),
    MinSigInd       = emptyfilter.double(6.0),
    MinSigCol       = emptyfilter.double(11.0),
    IndWidth        = emptyfilter.double(5.0),
    ColWidth        = emptyfilter.double(7.5),
    Drift           = emptyfilter.double(0.03069),
    # TPC's drift Velocity in units cm/(ADC Sample Time)
    POffset         = emptyfilter.double(20.0),
    OOffset         = emptyfilter.double(24.0),
    MaxMultiHit     = emptyfilter.int32(3)
    )

process.emptyfilter = emptyfilter.EDProducer(
    "EmptyFilter",
    HitsModuleLabel   = emptyfilter.string("ffthit"),
    MinHits          = emptyfilter.int32(30),
    MinIonization    = emptyfilter.double(750.0)
    )


# Write the events to the output file.
process.output = emptyfilter.OutputModule(
    "PoolOutputModule",
    fileName = emptyfilter.untracked.string('file:emptyfil_hists.root'),
)

####### End of the section that defines and configures modules.#########

# Tell the system to execute all paths. Services, source, output are implied ....
process.doit = emptyfilter.EndPath( process.caldataCal*process.ffthit*process.emptyfilter*process.output )

