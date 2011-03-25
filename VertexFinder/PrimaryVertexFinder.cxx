////////////////////////////////////////////////////////////////////////
//
// PrimaryVertexFinder class
//
// saima@ksu.edu
//
// This algorithm is designed to reconstruct the vertices by finding the
// 3d distance between the very first space points of each track 
// 
// This is Preliminary Work and needs modifications
// ////////////////////////////////////////////////////////////////////////


#include <iostream>

// Framework includes
#include "art/Framework/Core/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Persistency/Common/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "VertexFinder/PrimaryVertexFinder.h"

#include <iomanip>
#include <ios>
#include <sstream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include "TMath.h"


#include "RawData/RawDigit.h"
#include "Filters/ChannelFilter.h"
#include "SimulationBase/simbase.h"
#include "RecoBase/recobase.h"
#include "Geometry/geo.h"
#include "TH2.h"
#include "TVectorD.h"
#include "TGeoManager.h"


// //------------------------------------------------------------------------------
static bool sp_sort_z(const recob::SpacePoint& sp1, const recob::SpacePoint& sp2)
{
  const double* xyz1 = sp1.XYZ();
  const double* xyz2 = sp2.XYZ();
  return xyz1[2] < xyz2[2];
}
// //------------------------------------------------------------------------------
static bool sp_sort_zz(const std::vector<recob::SpacePoint>& sp_vec1, const std::vector<recob::SpacePoint>& sp_vec2)
{
  const double* xyz1 = sp_vec1[0].XYZ();
  const double* xyz2 = sp_vec2[0].XYZ();
  return xyz1[2] < xyz2[2];
}
// //------------------------------------------------------------------------------

namespace vertex{

//-----------------------------------------------------------------------------
  PrimaryVertexFinder::PrimaryVertexFinder(fhicl::ParameterSet const& pset) :
    fTrackModuleLabel  (pset.get< std::string >("TrackModuleLabel")),
    fVertexWindow      (pset.get<double     >  ("VertexWindow"))
  {
    produces< std::vector<recob::Vertex> >();
  }

//-----------------------------------------------------------------------------
  PrimaryVertexFinder::~PrimaryVertexFinder()
  {
  }

  //---------------------------------------------------------------------------
  void PrimaryVertexFinder::beginJob(){
    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;
    
    //    fNoVertices= tfs->make<TH2F>("fNoVertices", ";Event No; No of vertices", 100,0, 100, 30, 0, 30);
    
  }

// //-----------------------------------------------------------------------------
  void PrimaryVertexFinder::produce(art::Event& evt)
  {
    
    art::ServiceHandle<geo::Geometry> geom;
    
    //std::cout << "I am in Primary vertex finder " << std::endl;
        
    art::Handle< std::vector<recob::Track> > trackListHandle;
    evt.getByLabel(fTrackModuleLabel,trackListHandle);
    
    //Point to a collection of vertices to output.
    std::auto_ptr<std::vector<recob::Vertex> > vcol(new std::vector<recob::Vertex>);
    

    art::PtrVector<recob::Track> trkIn;
    for(unsigned int ii = 0; ii < trackListHandle->size(); ++ii)
      {
	art::Ptr<recob::Track> track(trackListHandle, ii);
	trkIn.push_back(track);
      }

    //std::cout << "number of tracks in this event = " << trkIn.size() << std::endl;

    std::vector<recob::SpacePoint> spacepoints;  // space points associated to each track 
    std::vector<recob::SpacePoint> startpoints_vec; // first space point of each track

    for(unsigned int j=0; j<trkIn.size();++j) { //loop over tracks
      spacepoints=trkIn[j]->SpacePoints();
      
      //std::cout<<"PrimaryVertexFinder got "<< spacepoints.size() <<" 3D spacepoint(s) from Track3Dreco.cxx"<<std::endl;
            
      for(unsigned int i=0; i<1; ++i) { // save the first SpacePoint of each Track... from now the SpacePoint ID represents the Track ID!!
	spacepoints[i].SetID(startpoints_vec.size());
	startpoints_vec.push_back(spacepoints[i]);
      }
    }// loop over tracks

    //get the collection of track IDs that have matched start points... (vertex_collection_int)
    std::vector<std::vector<int> > vertex_collection_int; 

    for (unsigned int i=0; i<startpoints_vec.size(); ++i){
      for (unsigned int j=i+1; j<startpoints_vec.size(); ++j){
	//std::cout << "distance between " << i << " and " << j << " = " << StartPointSeperation(startpoints_vec[i], startpoints_vec[j]) << std::endl;
	if(StartPointSeperation(startpoints_vec[i], startpoints_vec[j])<fVertexWindow){
	  if(!IsInVertexCollection(i, j, vertex_collection_int)){
	    std::vector<int> newvertex_int;
	    newvertex_int.push_back(i);
	    newvertex_int.push_back(j);
	    vertex_collection_int.push_back(newvertex_int);
	    //newvertex.clear();
	  }
	  else
	    {
	      int index = IndexInVertexCollection(i, j, vertex_collection_int);
	      //std::cout << "index where a new vertex will be added = " << index << std::endl;
	      if(!IsInNewVertex(i, vertex_collection_int[index])){
	      vertex_collection_int[index].push_back(i);
	      }
	     if(!IsInNewVertex(j, vertex_collection_int[index])){
	      vertex_collection_int[index].push_back(j);
	     }
	    }
	}
      }
    }

    //now add the unmatched track IDs to the collection
    for(int i=0; i<startpoints_vec.size(); i++){
      if(!IsInVertexCollection(i, i, vertex_collection_int)){
	std::vector<int> temp;
	temp.push_back(i);	
	vertex_collection_int.push_back(temp);
      }
    }

    //make spvertex_collection of 'SpacePoints' instead of spacepoint/Track IDs 
    std::vector<recob::SpacePoint> spvertex; 
    std::vector<std::vector<recob::SpacePoint> > spvertex_collection;

    for(unsigned int i=0; i<vertex_collection_int.size(); i++){
      for(std::vector<int>::iterator itr = vertex_collection_int[i].begin(); itr < vertex_collection_int[i].end(); ++itr){
	spvertex.push_back(startpoints_vec[*itr]);
      }
      std::sort(spvertex.begin(), spvertex.end(), sp_sort_z); // sort the start points in z
      spvertex_collection.push_back(spvertex);
      spvertex.clear();
    }

    std::sort(spvertex_collection.begin(), spvertex_collection.end(), sp_sort_zz); // sort the collection in z

    art::PtrVector<recob::Track> vTracks_vec;
    art::PtrVector<recob::Shower> vShowers_vec;
    
    for(unsigned int i=0; i < spvertex_collection.size() ; i++){ 
      double x = 0.;
      double y = 0.;
      double z = 0.;
      int elemsize = 0.;
      for(std::vector<recob::SpacePoint>::iterator itr = spvertex_collection[i].begin(); itr < spvertex_collection[i].end(); ++itr){
	vTracks_vec.push_back(trkIn[(*itr).ID()]);

	//calculate sum of x, y and z of a vertex
	x += ((*itr).XYZ())[0];
	y += ((*itr).XYZ())[1];
	z += ((*itr).XYZ())[2];
	elemsize = spvertex_collection[i].size();
      }

      double avgx = x/elemsize;
      double avgy = y/elemsize;
      double avgz = z/elemsize;
      
      Double_t vtxcoord[3];
	  vtxcoord[0] = avgx;
	  vtxcoord[1] = avgy;
	  vtxcoord[2] = avgz;    
            
	  recob::Vertex the3Dvertex(vTracks_vec, vShowers_vec, vtxcoord);
	  vcol->push_back(the3Dvertex);
	  vTracks_vec.clear();
    }
   
    mf::LogVerbatim("Summary") << std::setfill('-') << std::setw(175) << "-" << std::setfill(' ');
    mf::LogVerbatim("Summary") << "PrimaryVertexFinder Summary:";
    for(int i = 0; i<vcol->size(); ++i) mf::LogVerbatim("Summary") << vcol->at(i) ;

    evt.put(vcol);
    
  } // end of produce
} // end of vertex namespace

// //-----------------------------------------------------------------------------
double vertex::PrimaryVertexFinder::StartPointSeperation(recob::SpacePoint sp1, recob::SpacePoint sp2)
{
  double x= (sp2.XYZ()[0])-(sp1.XYZ()[0]);
  double y= (sp2.XYZ()[1])-(sp1.XYZ()[1]);
  double z= (sp2.XYZ()[2])-(sp1.XYZ()[2]);
  double distance = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
  return distance;
}
// //---------------------------------------------------------------------------------
bool vertex::PrimaryVertexFinder::IsInVertexCollection(int a, int b, std::vector<std::vector<int> > vertex_collection)
{
  int flag = 0;
  
  for(int i = 0; i < vertex_collection.size() ; i++){
    for(std::vector<int>::iterator itr = vertex_collection[i].begin(); itr < vertex_collection[i].end(); ++itr){
      if (a == *itr || b == *itr){
	flag = 1;
	break;
      }
    }
  }
  if(flag==1)
    return true;
  return false; 
}
// //------------------------------------------------------------------------------
int vertex::PrimaryVertexFinder::IndexInVertexCollection(int a, int b, std::vector<std::vector<int> > vertex_collection)
{
  int index;
  for(int i = 0; i < vertex_collection.size() ; i++){
    for(std::vector<int>::iterator itr = vertex_collection[i].begin(); itr < vertex_collection[i].end(); ++itr){
      if (a == *itr || b == *itr)
	index = i; 
    }
  }
  return index;
}
// //------------------------------------------------------------------------------
bool vertex::PrimaryVertexFinder::IsInNewVertex(int a, std::vector<int> newvertex)
{
  int flag = 0;
  for(int i = 0; i < newvertex.size() ; i++){
    if (a == newvertex[i]){
      flag = 1;
      break;
    }
  }
  
  if(flag==1)
    return true;
  return false; 
}
// //------------------------------------------------------------------------------
