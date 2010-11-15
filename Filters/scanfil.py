# Configuration file for Emptyfilter
#
# $Id: ex02.py,v 1.5 2010/05/18 21:24:21 kutschke Exp $
#
# Original author E. Church
#
# Spacing is not signficant in this file.

# Define the default configuration for the framework.
import FWCore.ParameterSet.python.Config as scanfilter

# Give this job a name.  
process = scanfilter.Process("Filters")

# Maximum number of events to do.
process.maxEvents = scanfilter.untracked.PSet(
    input = scanfilter.untracked.int32(30655) # See if this works to run fewer than are in input file.
)

# Load the standard message logger configuration.
# Threshold=Info. Limit of 5 per category; then exponential backoff.
#process.load("Config/MessageLogger_cfi")

# Load the service that manages root files for histograms.
process.TFileService = scanfilter.Service(
    "TFileService",
    fileName = scanfilter.string("scanfilter_hist.root"),
    closeFileFast = scanfilter.untracked.bool(False)
)

process.Timing = scanfilter.Service("Timing");

# Define the geometry.
process.Geometry = scanfilter.Service(
    "Geometry",
    SurfaceY = scanfilter.double(130.0e2), # in cm
    Name     = scanfilter.string("argoneut"),
    GDML     = scanfilter.FileInPath("Geometry/gdml/argoneut.gdml")
)




# Service to get my MC events, which were run up through DetSim.

# Define the FFT service
process.LArFFT = scanfilter.Service(
    "LArFFT",
    FFTSize   = scanfilter.int32(4096),#double the number of samples for ArgoNeuT to deal with long exponential tail
    FFTOption = scanfilter.string("ES")
)

process.source = scanfilter.Source("PoolSource",
                                fileNames = scanfilter.untracked.vstring("/argoneut/data/simplescan_data/fbtR650_handscanned_spitz7.root")
                                )

process.scanfilter = scanfilter.EDProducer(
    "ScanFilter",
    ScanModuleLabel   = scanfilter.string("merge"),
    DigitModuleLabel   = scanfilter.string("source"),
    Neutrino_req      = scanfilter.int32(1),  #0=no neutrinos, 1=maybe neutrino, 2=neutrino. 2 is most stringent, 0 is least stringent  
    NumShowers_req    = scanfilter.int32(0),  #Maximum number of showers required to pass
    NumTracks_req     = scanfilter.int32(2)   #Maximum number of tracks in any plane required to pass  
    )

# process.caldataCal = scanfilter.EDProducer(
#     "CalWire",
#     DigitModuleLabel  = scanfilter.string("scanfilter"),
#     ResponseFile       = scanfilter.string("/grid/fermiapp/lbne/lar/aux/ArgoResponse1.2.root"),
#     ExponentialEndBins = scanfilter.int32(300),
#     UseRawData         = scanfilter.int32(0)
#     )
# 
# process.ffthit = scanfilter.EDProducer(
#     "FFTHitFinder",
#     CalDataModuleLabel   = scanfilter.string("caldataCal"),
#     MinSigInd       = scanfilter.double(6.0),
#     MinSigCol       = scanfilter.double(11.0),
#     IndWidth        = scanfilter.double(5.0),
#     ColWidth        = scanfilter.double(7.5),
#     Drift           = scanfilter.double(0.03069),
#     # TPC's drift Velocity in units cm/(ADC Sample Time)
#     POffset         = scanfilter.double(20.0),
#     OOffset         = scanfilter.double(24.0),
#     MaxMultiHit     = scanfilter.int32(3)
#     )


# Write the events to the output file.
process.output = scanfilter.OutputModule(
    "PoolOutputModule",
    fileName = scanfilter.untracked.string('scanfilter_gen.root'),
)

####### End of the section that defines and configures modules.#########

# Tell the system to execute all paths. Services, source, output are implied ....
process.doit = scanfilter.EndPath(process.scanfilter*process.output )

