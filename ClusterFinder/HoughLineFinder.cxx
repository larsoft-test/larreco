////////////////////////////////////////////////////////////////////////
//
// \file HoughLineFinder.cxx
//
// \author joshua.spitz@yale.edu
//
//  This algorithm is designed to find lines (Houghclusters) from clusters found by DBSCAN 
//  after deconvolution and hit finding.
//  The algorithm is based on: 
//  Queisser, A. "Computing the Hough Transform", C/C++ Users Journal 21, 12 (Dec. 2003).
//  Niblack, W. and Petkovic, D. On Improving the Accuracy of the Hough Transform", Machine 
//  Vision and Applications 3, 87 (1990)  
////////////////////////////////////////////////////////////////////////

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}
#include <sstream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <vector>

// ROOT includes
#include <TCanvas.h>
#include "TDatabasePDG.h"
#include "TSystem.h"

// ART includes
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

// LArSoft includes 
#include "RawData/RawDigit.h"
#include "Filters/ChannelFilter.h"
#include "SimulationBase/simbase.h"
#include "RecoBase/recobase.h"
#include "Utilities/AssociationUtil.h"
#include "Geometry/geo.h"
#include "ClusterFinder/HoughBaseAlg.h"
#include "ClusterFinder/HoughLineFinder.h"

//------------------------------------------------------------------------------
cluster::HoughLineFinder::HoughLineFinder(fhicl::ParameterSet const& pset) 
  : fHLAlg(pset.get< fhicl::ParameterSet >("HoughBaseAlg"))
{
  this->reconfigure(pset);
  produces< std::vector<recob::Cluster> >();
  produces< art::Assns<recob::Cluster, recob::Hit> >();
}

//------------------------------------------------------------------------------
cluster::HoughLineFinder::~HoughLineFinder()
{
}

//------------------------------------------------------------------------------
void cluster::HoughLineFinder::reconfigure(fhicl::ParameterSet const& p)
{
  fDBScanModuleLabel = p.get< std::string >("DBScanModuleLabel");
  fHLAlg.reconfigure(p.get< fhicl::ParameterSet >("HoughBaseAlg"));
}

//------------------------------------------------------------------------------
void cluster::HoughLineFinder::produce(art::Event& evt)
{

  //////////////////////////////////////////////////////
  // here is how to get a collection of objects out of the file
  // and connect it to a art::Handle
  //////////////////////////////////////////////////////
  // Read in the clusterList object(s).
  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fDBScanModuleLabel,clusterListHandle);

  //art::PtrVector<recob::Cluster> clusIn;
  std::vector<art::Ptr<recob::Cluster> > clusIn;
  for(unsigned int ii = 0; ii < clusterListHandle->size(); ++ii){
    art::Ptr<recob::Cluster> cluster(clusterListHandle, ii);
    clusIn.push_back(cluster);
  }
  
  // make a std::vector<recob::Cluster> for the output of the 
  // Hough Transform and a std::vector< art::PtrVector<recob::Hit> >
  // to hold the associated hits
  std::vector<recob::Cluster>               clusOut;
  std::vector< art::PtrVector<recob::Hit> > clusHitsOut;
  
  size_t numclus = fHLAlg.Transform(clusIn, clusOut, clusHitsOut, evt, fDBScanModuleLabel);
    
  //size_t Transform(std::vector<art::Ptr<recob::Cluster> >           & clusIn,
                          //std::vector<recob::Cluster>                      & ccol,  
		     //std::vector< art::PtrVector<recob::Hit> >      & clusHitsOut,
		     //art::Event                                const& evt,
		     //std::string                               const& label);

  LOG_DEBUG("HoughLineClusters") << "found " << numclus << "clusters with HoughBaseAlg";

  //Point to a collection of clusters to output.
  std::auto_ptr<std::vector<recob::Cluster> > ccol(new std::vector<recob::Cluster>(clusOut));
  std::auto_ptr< art::Assns<recob::Cluster, recob::Hit> > assn(new art::Assns<recob::Cluster, recob::Hit>);

  mf::LogVerbatim("Summary") << std::setfill('-') << std::setw(175) << "-" << std::setfill(' ');
  mf::LogVerbatim("Summary") << "HoughLineFinder Summary:";
  for(size_t i = 0; i < ccol->size(); ++i){
    mf::LogVerbatim("Summary") << ccol->at(i);

    // associat the hits to this cluster
    util::CreateAssn(*this, evt, *(ccol.get()), clusHitsOut[i], *(assn.get()), i);
  }

  evt.put(ccol);
  evt.put(assn);
  return;
}



