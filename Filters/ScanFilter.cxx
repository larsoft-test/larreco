////////////////////////////////////////////////////////////////////////
//
// ScanFilter class
//
// joshua.spitz@yale.edu
//
////////////////////////////////////////////////////////////////////////

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}
#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>

//Framework Includes
#include "FWCore/Framework/interface/Event.h" 
#include "FWCore/ParameterSet/interface/ParameterSet.h" 
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h" 
#include "DataFormats/Common/interface/Handle.h" 
#include "DataFormats/Common/interface/Ptr.h" 
#include "DataFormats/Common/interface/PtrVector.h" 
#include "FWCore/Framework/interface/MakerMacros.h" 
#include "FWCore/ServiceRegistry/interface/Service.h" 
#include "FWCore/Services/interface/TFileService.h" 
#include "FWCore/Framework/interface/TFileDirectory.h" 
#include "FWCore/MessageLogger/interface/MessageLogger.h" 
#include "FWCore/ServiceRegistry/interface/ServiceMaker.h" 


//Larsoft Includes
#include "ScanFilter.h"
#include "RecoBase/recobase.h"
#include "T962_MergeData/ScanInfo.h"
#include "Geometry/geo.h"
#include "Geometry/WireGeo.h"

namespace filt{

  //-------------------------------------------------
  ScanFilter::ScanFilter(edm::ParameterSet const & pset) : 
    fScanModuleLabel      (pset.getParameter< std::string > ("ScanModuleLabel")),
    fNeutrino_req         (pset.getParameter< int >("Neutrino_req")),
    fNumShowers_req       (pset.getParameter< int >("NumShowers_req")),
    fNumTracks_req        (pset.getParameter< int >("NumTracks_req"))
  {   
  }

  //-------------------------------------------------
  ScanFilter::~ScanFilter()
  {
  }

  void ScanFilter::beginJob(const edm::EventSetup&)
  {
    return;
  }

  //-------------------------------------------------
  bool ScanFilter::filter(edm::Event &evt, edm::EventSetup const&)
  { 

    int failFlag = 1;
    int run = evt.id().run();
    int event = evt.id().event();
    
    edm::PtrVector<merge::ScanInfo> scanIn;
    scanIn.clear();

    edm::Service<geo::Geometry> geom;
    edm::Handle< std::vector<merge::ScanInfo> > scanHandle;
    evt.getByLabel(fScanModuleLabel,scanHandle);

    for(unsigned int i = 0; i < scanHandle->size(); ++i){     
     edm::Ptr<merge::ScanInfo> scaninfo(scanHandle, i);
      scanIn.push_back(scaninfo);     
    }

    for(unsigned int i = 0; i < scanIn.size(); ++i){
    if(scanIn[i]->Get_IsNeutrino()>=fNeutrino_req 
    && scanIn[i]->Get_NumShower()<=fNumShowers_req  
    && (scanIn[i]->Get_TrackInd()<=fNumTracks_req||scanIn[i]->Get_TrackCol()<=fNumTracks_req) 
    && scanIn[i]->Get_Run()==run 
    && scanIn[i]->Get_Event()==event)
    failFlag=0;
    }
 
    if(failFlag>0)
    return  false;

    return true;
  }
	      
} //end namespace
