////////////////////////////////////////////////////////////////////////
//
// HoughLineFinder class
//
// joshua.spitz@yale.edu
//
//  This algorithm is designed to find lines (Houghclusters) from clusters found by DBSCAN after deconvolution and hit finding.
//  The algorithm is based on: 
//  Queisser, A. "Computing the Hough Transform", C/C++ Users Journal 21, 12 (Dec. 2003).
//  Niblack, W. and Petkovic, D. On Improving the Accuracy of the Hough Transform", Machine Vision and Applications 3, 87 (1990)  
////////////////////////////////////////////////////////////////////////

#include "ClusterFinder/HoughLineService.h"
#include "ClusterFinder/HoughLineModule.h"
extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}
// ROOT includes
#include <TCanvas.h>
#include "TDatabasePDG.h"
#include "TSystem.h"

#include <sstream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <vector>

#include "art/Framework/Core/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Persistency/Common/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Core/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
 
#include "RawData/RawDigit.h"
#include "Filters/ChannelFilter.h"
#include "SimulationBase/simbase.h"
#include "RecoBase/recobase.h"
#include "Geometry/geo.h"


cluster::HoughLineModule::HoughLineModule(fhicl::ParameterSet const& pset) : 
  fDBScanModuleLabel       (pset.get< std::string >("DBScanModuleLabel"))
  
{
  produces< std::vector<recob::Cluster> >();
}

cluster::HoughLineModule::~HoughLineModule()
{
}



void cluster::HoughLineModule::produce(art::Event& evt)
{

  //////////////////////////////////////////////////////
  // here is how to get a collection of objects out of the file
  // and connect it to a art::Handle
  //////////////////////////////////////////////////////
  // Read in the clusterList object(s).
  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fDBScanModuleLabel,clusterListHandle);

  art::PtrVector<recob::Cluster> clusIn;
  for(unsigned int ii = 0; ii < clusterListHandle->size(); ++ii)
    {
      art::Ptr<recob::Cluster> cluster(clusterListHandle, ii);
      clusIn.push_back(cluster);
    }
  
  art::ServiceHandle<cluster::HoughLineService> hls;
  
  // make a std::vector<recob::Cluster> for the output of the 
  // Hough Transform
  std::vector<recob::Cluster> clusOut;
  
  size_t numclus = hls->Transform(clusIn, clusOut);

  mf::LogDebug("HoughLineClusters") << "found " << numclus << "clusters with HoughLineService";

  //Point to a collection of clusters to output.
  std::auto_ptr<std::vector<recob::Cluster> > ccol(new std::vector<recob::Cluster>(clusOut));

  std::sort(ccol->begin(),ccol->end());//sort before Putting

  std::cout << std::setfill('-') << std::setw(175) << "-" << std::endl;
  std::cout << std::setfill(' ');
  std::cout << "HoughLineModule Summary:" << std::endl;
  for(int i = 0; i<ccol->size(); ++i) std::cout << ccol->at(i) << std::endl;
  std::cout << std::endl;

  evt.put(ccol);
 
    

 
}



