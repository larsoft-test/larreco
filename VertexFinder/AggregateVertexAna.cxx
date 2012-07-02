////////////////////////////////////////////////////////////////////////
// $Id: AggVertex.cxx,v 1.1 2010/09/02 17:25:11 echurch Exp $
//
// AggregateVertexAna class
//
//  This class will produce the already-discovered vertices
//  containing a vector of track pointers  which we associate
//  to them. 
//
// echurch@fnal.gov
//
////////////////////////////////////////////////////////////////////////
extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

// Framework includes
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

#include "VertexFinder/AggregateVertexAna.h"
#include "Utilities/AssociationUtil.h"

#include "TH1.h"

namespace vertex{

  //-----------------------------------------------
  AggregateVertexAna::AggregateVertexAna(fhicl::ParameterSet const& pset)
    : fHitModuleLabel     (pset.get< std::string >("FFFTHitModuleLabel"))
    , fTrack3DModuleLabel (pset.get< std::string >("Track3DModuleLabel"))
    , fEndPointModuleLabel(pset.get< std::string >("EndPointModuleLabel"))
    , fVertexModuleLabel  (pset.get< std::string >("VertexModuleLabel"))
  {


  }

  //-----------------------------------------------
  AggregateVertexAna::~AggregateVertexAna()
  {
  }

  //-----------------------------------------------
  void AggregateVertexAna::beginJob()
  {
    art::ServiceHandle<art::TFileService> tfs;

    HnVtxes   = tfs->make<TH1F>("Num Vertices","Num Vertices",8,-0.5,7.5);
    HVtxSep   = tfs->make<TH1F>("Vertices spacing","Vertices spacing",20,0.001,5.0);
    HVtxRZ    = tfs->make<TH2F>("Vtx in RZ","Vtx in RZ",20,-50.0,+50.0,20,0.0,50.0);
    HnTrksVtx = tfs->make<TH1F>("Tracks per vtx","Tracks per vtx",8,-0.5,7.5);

    return;
  }

  //-----------------------------------------------
  void AggregateVertexAna::analyze(const art::Event& evt) 
  {
    art::Handle< std::vector<recob::Hit> > hitListHandle;
    evt.getByLabel(fHitModuleLabel,hitListHandle);
    for(unsigned int ii = 0; ii < hitListHandle->size(); ++ii){
      art::Ptr<recob::Hit> hit(hitListHandle, ii);
      fhitlist.push_back(hit); // class member
    }

    art::Handle< std::vector<recob::EndPoint2D> > epListHandle;
    evt.getByLabel(fEndPointModuleLabel,epListHandle);
    for(unsigned int ii = 0; ii < epListHandle->size(); ++ii){
      art::Ptr<recob::EndPoint2D> ep(epListHandle, ii);
      feplist.push_back(ep); // class member
    }

    art::Handle< std::vector<recob::Track> > trackListHandle;
    evt.getByLabel(fTrack3DModuleLabel,trackListHandle);
    for(unsigned int ii = 0; ii < trackListHandle->size(); ++ii){
      art::Ptr<recob::Track> track(trackListHandle, ii);
      ftracklist.push_back(track); // class member
    }

    art::Handle< std::vector<recob::Vertex> > vertexListHandle;
    evt.getByLabel(fVertexModuleLabel,vertexListHandle);
    for(unsigned int ii = 0; ii < vertexListHandle->size(); ++ii){
      art::Ptr<recob::Vertex> vertex(vertexListHandle, ii);
      fVertexlist.push_back(vertex); // class member
    }

    HnVtxes->Fill(feplist.size(),1);  

    art::FindManyP<recob::Track> fmt(vertexListHandle, evt, fVertexModuleLabel);
    art::FindManyP<recob::Hit>   fmh(vertexListHandle, evt, fVertexModuleLabel);

    for(size_t v1 = 0; v1 < fVertexlist.size(); ++v1)  {            
      
      std::vector< art::Ptr<recob::Track> > tvlist = fmt.at(v1);
      
      HnTrksVtx->Fill(tvlist.size(),1);
      
      if(tvlist.size() < 1) continue;
      
      std::vector< art::Ptr<recob::Hit> > hitvlist = fmh.at(v1);
      
      // Hits no longer has XYZ() method. To get 3d hit position info I'm going to have to 
      // loop on all the SpacePoints and loop on all Hits from there till it matches
      // one from this vertex. This affects the two Fill() efforts below. EC, 19-Nov-2010.      
      art::PtrVector<recob::Hit>::const_iterator hitv = hitvlist.begin();
      
      for(size_t v2 = v1+1; v2 < fVertexlist.size(); ++v2){            
	
	std::vector< art::Ptr<recob::Hit> > hitvlist2 = fmh.at(v2);
	
	std::vector< art::Ptr<recob::Hit> >::const_iterator hitv2 = hitvlist2.begin();
	
	// These two whiles should be each precisely one iteration long.
	while( hitv != hitvlist.end() ){
	  while( hitv2 != hitvlist2.end() ){
	    TVector3 dist;
	    mf::LogInfo("AggregateVertexAna") << "AggregateVertexAna: dist is " << dist.Mag() << ".";
	    HVtxSep->Fill(dist.Mag(),1);
	    hitv2++;
	  }
	  hitv++;
	}// end loop over hitv entries
	
      }// end loop over v2
    }// end loop over v1

    return;
  }// end analyze
}// end namespace
