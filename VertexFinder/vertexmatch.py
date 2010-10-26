# Configuration file for VertexMatch
#
# $Id: ex02.py,v 1.5 2010/05/18 21:24:21 kutschke Exp $
#
# Original author J. Spitz   
#   joshua.spitz@yale.edu
#
# Spacing is not signficant in this file.

# Define the default configuration for the framework.
import FWCore.ParameterSet.python.Config as vertexfinder

# Give this job a name.  
process = vertexfinder.Process("VertexMatch")

# Maximum number of events to do.
process.maxEvents = vertexfinder.untracked.PSet(
    input = vertexfinder.untracked.int32(3) # See if this works to run fewer than are in input file.
)

# Load the standard message logger configuration.
# Threshold=Info. Limit of 5 per category; then exponential backoff.
#process.load("Config/MessageLogger_cfi")

# Load the service that manages root files for histograms.
process.TFileService = vertexfinder.Service(
    "TFileService",
    fileName = vertexfinder.string("vertexfinder_hist.root"),
    closeFileFast = vertexfinder.untracked.bool(False)
)

process.Timing = vertexfinder.Service("Timing");

# Define the geometry.
process.Geometry = vertexfinder.Service(
    "Geometry",
    SurfaceY = vertexfinder.double(130.0e2), # in cm
    Name     = vertexfinder.string("argoneut"),
    GDML     = vertexfinder.FileInPath("Geometry/gdml/argoneut.gdml")
)


# Define the FFT service
process.LArFFT = vertexfinder.Service(
    "LArFFT",
    FFTSize   = vertexfinder.int32(4096),#double the number of samples for ArgoNeuT to deal with long exponential tail
    FFTOption = vertexfinder.string("ES")
)

# Service to get my MC events, which were run up through DetSim.
process.source = vertexfinder.Source("PoolSource",
                                fileNames = vertexfinder.untracked.vstring("/argoneut/app/users/echurch/larsoft/ART-SRT/out/genie_ART_6GeV-c_muons_gen.root")
                                )

process.caldataCal = vertexfinder.EDProducer(
    "CalWire",
    DigitModuleLabel  = vertexfinder.string("wiresim"),
    ResponseFile       = vertexfinder.string("/grid/fermiapp/lbne/lar/aux/ArgoResponse1.2.root"),
    ExponentialEndBins = vertexfinder.int32(300),
    UseRawData         = vertexfinder.int32(0)
    )

process.ffthit = vertexfinder.EDProducer(
    "FFTHitFinder",
    CalDataModuleLabel   = vertexfinder.string("caldataCal"),
    MinSigInd       = vertexfinder.double(8.0),
    MinSigCol       = vertexfinder.double(11.0),
    IndWidth        = vertexfinder.double(5.0),
    ColWidth        = vertexfinder.double(7.5),
    Drift           = vertexfinder.double(0.03069),
    # TPC's drift Velocity in units cm/(ADC Sample Time)
    POffset         = vertexfinder.double(20.0),
    OOffset         = vertexfinder.double(24.0),
    MaxMultiHit     = vertexfinder.int32(3)
    )

process.hough = vertexfinder.EDProducer(
    "HoughLineFinder",
    DBScanModuleLabel        = vertexfinder.string("dbscan"),
    MaxLines                 = vertexfinder.int32(5),
    MinHits                  = vertexfinder.int32(3),
    SaveAccumulator          = vertexfinder.int32(0),
    NumAngleCells            = vertexfinder.int32(10000),
    RhoResolutionFactor      = vertexfinder.int32(10),
    SmootherSigma            = vertexfinder.double(0.),
    MaxDistance              = vertexfinder.double(5.),
    RhoZeroOutRange          = vertexfinder.int32(0),
    ThetaZeroOutRange        = vertexfinder.int32(0),
    HitsPerCluster           = vertexfinder.int32(1)
    )
    
process.harris = vertexfinder.EDProducer(
    "HarrisVertexFinder",
    HitsModuleLabel      = vertexfinder.string("ffthit"),
    TimeBins             = vertexfinder.int32(256),
    MaxCorners           = vertexfinder.int32(20),
    Gsigma               = vertexfinder.double(1.),
    Window               = vertexfinder.int32(5),
    Threshold            = vertexfinder.double(0.1),
    SaveVertexMap        = vertexfinder.int32(-1)
    )  
    
process.vertexmatch = vertexfinder.EDProducer(
    "VertexMatch",    
    HoughModuleLabel       = vertexfinder.string("hough"),
    VertexModuleLabel      = vertexfinder.string("harris"),
    MaxDistance            = vertexfinder.double(30.)
    )     

# Write the events to the output file.
process.output = vertexfinder.OutputModule(
    "PoolOutputModule",
    fileName = vertexfinder.untracked.string('vertexmatch_gen.root'),
)

####### End of the section that defines and configures modules.#########

# Tell the system to execute all paths. Services, source, output are implied ....
process.doit = vertexfinder.EndPath(process.caldataCal*process.ffthit*process.hough*process.harris*process.vertexmatch*process.output )

