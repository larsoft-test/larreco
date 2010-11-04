# Configuration file for ClusterFinder
#
# $Id: ex02.py,v 1.5 2010/05/18 21:24:21 kutschke Exp $
#
# Original author J. Spitz   
#   joshua.spitz@yale.edu
#
# Spacing is not signficant in this file.

# Define the default configuration for the framework.
import FWCore.ParameterSet.python.Config as clusterfinder

# Give this job a name.  
process = clusterfinder.Process("ClusterFinder")

# Maximum number of events to do.
process.maxEvents = clusterfinder.untracked.PSet(
    input = clusterfinder.untracked.int32(1) # See if this works to run fewer than are in input file.
)

# Load the standard message logger configuration.
# Threshold=Info. Limit of 5 per category; then exponential backoff.
#process.load("Config/MessageLogger_cfi")

# Load the service that manages root files for histograms.
process.TFileService = clusterfinder.Service(
    "TFileService",
    fileName = clusterfinder.string("hough_hist.root"),
    closeFileFast = clusterfinder.untracked.bool(False)
)

process.Timing = clusterfinder.Service("Timing");

# Define the geometry.
process.Geometry = clusterfinder.Service(
    "Geometry",
    SurfaceY = clusterfinder.double(130.0e2), # in cm
    Name     = clusterfinder.string("argoneut"),
    GDML     = clusterfinder.FileInPath("Geometry/gdml/argoneut.gdml")
)


# Define the FFT service
process.LArFFT = clusterfinder.Service(
    "LArFFT",
    FFTSize   = clusterfinder.int32(4096),#double the number of samples for ArgoNeuT to deal with long exponential tail
    FFTOption = clusterfinder.string("ES")
)

# Service to get my MC events, which were run up through DetSim.
process.source = clusterfinder.Source("PoolSource",
                                fileNames = clusterfinder.untracked.vstring("/argoneut/app/users/spitz7/larsoft_new2/genie_ART_6GeV-c_100-muons_gen.root")
                                )

process.caldataCal = clusterfinder.EDProducer(
    "CalWire",
    DigitModuleLabel  = clusterfinder.string("wiresim"),
    ResponseFile       = clusterfinder.string("/grid/fermiapp/lbne/lar/aux/ArgoResponse1.2.root"),
    ExponentialEndBins = clusterfinder.int32(300),
    UseRawData         = clusterfinder.int32(0)
    )

process.ffthit = clusterfinder.EDProducer(
    "FFTHitFinder",
    CalDataModuleLabel   = clusterfinder.string("caldataCal"),
    MinSigInd       = clusterfinder.double(6.0),
    MinSigCol       = clusterfinder.double(11.0),
    IndWidth        = clusterfinder.double(5.0),
    ColWidth        = clusterfinder.double(7.5),
    Drift           = clusterfinder.double(0.03069),
    # TPC's drift Velocity in units cm/(ADC Sample Time)
    POffset         = clusterfinder.double(20.0),
    OOffset         = clusterfinder.double(24.0),
    MaxMultiHit     = clusterfinder.int32(3)
    )
    
 
process.dbscan = clusterfinder.EDProducer(
    "DBcluster",
    HitsModuleLabel   = clusterfinder.string("ffthit"),
    eps               = clusterfinder.double(1.0),
    eps2              = clusterfinder.double(0.75),
    minPts            = clusterfinder.int32(2)
    )  

process.hough = clusterfinder.EDProducer(
    "HoughLineFinder",
    DBScanModuleLabel         = clusterfinder.string("dbscan"),
    MaxLines                 = clusterfinder.int32(5),
    MinHits                  = clusterfinder.int32(3),
    SaveAccumulator          = clusterfinder.int32(0),
    NumAngleCells            = clusterfinder.int32(10000),
    RhoResolutionFactor      = clusterfinder.int32(10),
    SmootherSigma            = clusterfinder.double(0.),
    MaxDistance              = clusterfinder.double(5.),
    RhoZeroOutRange          = clusterfinder.int32(0),
    ThetaZeroOutRange        = clusterfinder.int32(0),
    HitsPerCluster           = clusterfinder.int32(1)
    )

# Write the events to the output file.
process.output = clusterfinder.OutputModule(
    "PoolOutputModule",
    fileName = clusterfinder.untracked.string('hough_gen.root'),
)

####### End of the section that defines and configures modules.#########

# Tell the system to execute all paths. Services, source, output are implied ....
process.doit = clusterfinder.EndPath(process.caldataCal*process.ffthit*process.dbscan*process.hough*process.output )

