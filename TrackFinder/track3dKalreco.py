# Configuration file for Track3Dreco	
#
# Spacing is not signficant in this file.

# Define the default configuration for the framework.
import FWCore.ParameterSet.python.Config as trackfinder

# Give this job a name.  
process = trackfinder.Process("TrackFinder")

# Maximum number of events to do.
process.maxEvents = trackfinder.untracked.PSet(
    input = trackfinder.untracked.int32(10) 
)

# Load the standard message logger configuration.
# Threshold=Info. Limit of 5 per category; then exponential backoff.
#process.load("Config/MessageLogger_cfi")

# Load the service that manages root files for histograms.
process.TFileService = trackfinder.Service(
    "TFileService",
    fileName = trackfinder.string("out/genie_Argo_500-5000MeV-c_Kal_30degSpread_muons_tracked_hist.root"),
#    track3dreco_ArgoNeuT_CCQE_hist.root
    closeFileFast = trackfinder.untracked.bool(False)
)

process.Timing = trackfinder.Service("Timing");

# Define the geometry.
process.Geometry = trackfinder.Service(
    "Geometry",
    SurfaceY = trackfinder.double(130.0e2), # in cm
    Name     = trackfinder.string("argoneut"),
    GDML     = trackfinder.FileInPath("Geometry/gdml/argoneut.gdml")
)


# Define the FFT service
process.LArFFT = trackfinder.Service(
    "LArFFT",
    FFTSize   = trackfinder.int32(4096),#double the number of samples for ArgoNeuT to deal with long exponential tail
    FFTOption = trackfinder.string("ES")
)

# Define the LAr property service
process.LArProperties = trackfinder.Service(
    "LArProperties",
    Efield       = trackfinder.double(0.5),   #kV/cm
    Temperature  = trackfinder.double(87.),     #kelvin
    Electronlifetime   = trackfinder.double(1000000000.)     #microseconds
)


process.source = trackfinder.Source("PoolSource",
#                                fileNames = trackfinder.untracked.vstring("out/R650_D20090928_T112738.root"),
#                                    fileNames = trackfinder.untracked.vstring("out/genie_ART_400MeV-c_electrons_gen.root")
                                    fileNames = trackfinder.untracked.vstring("out/genie_ART_500-5500MeV-c_30degSpread_muons_gen.root")
#                                    genie_ART_400MeV-c_muons_gen.root"
                                    # 1502,4715,5552,5529 -- CCQE
                                    #skipEvents = trackfinder.untracked.uint32(1)
#                                eventsToProcess = trackfinder.untracked.vEventRange('650:1502-650:1503','650:4714-650:4716'),
                                )

process.caldataCal = trackfinder.EDProducer(
    "CalWire",
#    DigitModuleLabel  = trackfinder.string("source"),
    DigitModuleLabel  = trackfinder.string("wiresim"),
    ResponseFile       = trackfinder.string("/grid/fermiapp/lbne/lar/aux/ArgoResponse1.2.root"),
    ExponentialEndBins = trackfinder.int32(300),
    UseRawData         = trackfinder.int32(0)
    )

process.ffthit = trackfinder.EDProducer(
    "FFTHitFinder",
    CalDataModuleLabel   = trackfinder.string("caldataCal"),
    MinSigInd       = trackfinder.double(6.0),
    MinSigCol       = trackfinder.double(11.0),
    IndWidth        = trackfinder.double(5.0),
    ColWidth        = trackfinder.double(7.5),
    Drift           = trackfinder.double(0.03069),
    # TPC's drift Velocity in units cm/(ADC Sample Time)
    POffset         = trackfinder.double(20.0),
    OOffset         = trackfinder.double(24.0),
    MaxMultiHit     = trackfinder.int32(3)
    )
    
 
process.dbscan = trackfinder.EDProducer(
    "DBcluster",
    HitsModuleLabel   = trackfinder.string("ffthit"),
    eps               = trackfinder.double(1.0),
    eps2              = trackfinder.double(0.75),
    minPts            = trackfinder.int32(2)
    )

process.dbscanana = trackfinder.EDAnalyzer(
    "DBclusterAna",
    DigitModuleLabel =trackfinder.string("wiresim"),
     HitsModuleLabel=trackfinder.string("ffthit"),
     LArG4ModuleLabel=trackfinder.string("largeant"),
     CalDataModuleLabel=trackfinder.string("caldataCal"),
     GenieGenModuleLabel=trackfinder.string("singlegen"),
     ClusterFinderModuleLabel=trackfinder.string("dbscan")
    )    


process.hough = trackfinder.EDProducer(
    "HoughLineFinder",
    DBScanModuleLabel        = trackfinder.string("dbscan"),
    MaxLines                 = trackfinder.int32(5),
    MinHits                  = trackfinder.int32(3),
    SaveAccumulator          = trackfinder.int32(0),
    NumAngleCells            = trackfinder.int32(10000),
    RhoResolutionFactor      = trackfinder.int32(10),
    SmootherSigma            = trackfinder.double(0.),
    MaxDistance              = trackfinder.double(5.),
    RhoZeroOutRange          = trackfinder.int32(0),
    ThetaZeroOutRange        = trackfinder.int32(0),
    HitsPerCluster           = trackfinder.int32(1)
    )

process.linemerger = trackfinder.EDProducer(
    "LineMerger",
    ClusterModuleLabel       = trackfinder.string("hough"),
    Slope                    = trackfinder.double(1.5),
    Intercept                = trackfinder.double(102.0)
    )	

process.track3DKal = trackfinder.EDProducer(
    "Track3DKalman",
    HitModuleLabel       = trackfinder.string("ffthit"),
    GenieGenModuleLabel  = trackfinder.string("singlegen"), # comment out for real data.
    TMatch               = trackfinder.int32(10), # 11 is delicate. 10 cuts space hits by fac of 10,
    #                                               12 increases by factor two. EC, 5-Jan-2011.
    Chi2DOFmax           = trackfinder.double(10.0),
    GenfPRINT            = trackfinder.bool(True)
    )	


process.calitalians = trackfinder.EDAnalyzer(
    "CaloArgoItaliano",
    TrackModuleLabel       = trackfinder.string("track3DKal"),
    SpacePointDimension    = trackfinder.int32(2000),
    HitCollectionDimension    = trackfinder.int32(1000),
    HitInductionDimension    = trackfinder.int32(1000)
    )	


# Write the events to the output file.
process.output = trackfinder.OutputModule(
    "PoolOutputModule",
#    fileName = trackfinder.untracked.string('out/track3dreco_ArgoNeuT_CCQE.root'),
    fileName = trackfinder.untracked.string('out/genie_Argo_500-5500MeV-c_Kal_30degSpread_muons_tracked_gen.root')
)

####### End of the section that defines and configures modules.#########

# Tell the system to execute all paths. Services, source, output are implied ....
process.doit = trackfinder.EndPath(process.caldataCal*process.ffthit*process.track3DKal*process.calitalians*process.output )
#process.doit = trackfinder.EndPath(process.caldataCal*process.ffthit*process.dbscan*process.tracksoderberg*process.output )

