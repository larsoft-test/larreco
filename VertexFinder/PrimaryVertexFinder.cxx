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
bool sort_pred2(const std::pair<art::Ptr<recob::Track>,double>& left, const std::pair<art::Ptr<recob::Track>,double>& right)
{
  return left.second < right.second;
}

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

    std::cout<<std::endl;
    std::cout<<"------------------------------------------------------------------------------"<<std::endl;
    //   std::cout << "run    : " << evt.Header().Run() << std::endl;
    //   std::cout << "subrun : " << evt.Header().Subrun() << std::endl;
    //std::cout << "event  : " << evt.Header().Event() << std::endl;
    
    std::cout << "event  : " << evt.id().event() << std::endl;
    
    
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

    std::cout << "number of tracks in this event = " << trkIn.size() << std::endl;

    std::vector<recob::SpacePoint> spacepoints;  // space points associated to each track 
    std::vector<recob::SpacePoint> startpoints_vec; // first space point of each track

    std::vector<double> start;
    std::vector<double> end;
    double startcos[3]={0,0,0};
    double endcos[3]={0,0,0};

    std::vector <TVector3> startvec;
    TVector3 startXYZ;
    
    std::vector <TVector3> endvec;
    TVector3 endXYZ;

    std::vector <TVector3> dircosvec;
    TVector3 dircosXYZ;

    for(unsigned int j=0; j<trkIn.size();++j) { //loop over tracks
      spacepoints=trkIn[j]->SpacePoints();

     
      trkIn[j]->Extent(start, end);
      trkIn[j]->Direction(startcos, endcos);

      startXYZ.SetXYZ(start[0],start[1],start[2]);
      endXYZ.SetXYZ(end[0],end[1],end[2]);
      dircosXYZ.SetXYZ(startcos[0],startcos[1],startcos[2]);
 
      startvec.push_back(startXYZ);
      endvec.push_back(endXYZ);
      dircosvec.push_back(dircosXYZ);
      
      start.clear();
      end.clear();
      startcos[3] = 0.;
      
      std::cout<<"PrimaryVertexFinder got "<< spacepoints.size() <<" 3D spacepoint(s) from Track3Dreco.cxx"<<std::endl;
            
      for(unsigned int i=0; i<1; ++i) { // save the first SpacePoint of each Track... from now the SpacePoint ID represents the Track ID!!
	spacepoints[i].SetID(startpoints_vec.size());
	startpoints_vec.push_back(spacepoints[i]);
      }
    }// loop over tracks

    for(unsigned int i=0; i<startvec.size(); i++){
      std::cout << "Tvector3 start point = " << std::endl; 
      startvec[i].Print();
    }
    for(unsigned int i=0; i<dircosvec.size(); i++){
      std::cout << "Tvector3 dir cos = " << std::endl; 
      dircosvec[i].Print();
    }

     //see the start points
     for(unsigned int i=0; i<startpoints_vec.size(); ++i){
 //       double posX = startpoints_vec[i].XYZ()[0];
 //       double posY = startpoints_vec[i].XYZ()[1];
 //       double posZ = startpoints_vec[i].XYZ()[2];
 //       std::cout << "start point = " << posX << ", " << posY << ", " << posZ << std::endl;
       std::cout << "Start Points from Track3Dreco.cxx = " << startpoints_vec[i] << std::endl;
     }


    //get the collection of track IDs that have matched start points... (vertex_collection_int)
//     std::vector<std::vector<int> > vertex_collection_int; 

//     for (unsigned int i=0; i<startpoints_vec.size(); ++i){
//       for (unsigned int j=i+1; j<startpoints_vec.size(); ++j){
// 	std::cout << "distance between " << i << " and " << j << " = " << StartPointSeperation(startpoints_vec[i], startpoints_vec[j]) << std::endl;
// 	if(StartPointSeperation(startpoints_vec[i], startpoints_vec[j])<fVertexWindow){
// 	  if((!IsInVertexCollection(i, vertex_collection_int)) && (!IsInVertexCollection(j, vertex_collection_int))){
// 	    std::vector<int> newvertex_int;
// 	    newvertex_int.push_back(i);
// 	    newvertex_int.push_back(j);
// 	    vertex_collection_int.push_back(newvertex_int);
// 	    //newvertex.clear();
// 	  }
// 	  else
// 	    {
// 	      int index = IndexInVertexCollection(i, j, vertex_collection_int);
// 	      //std::cout << "index where a new vertex will be added = " << index << std::endl;
// 	      if(!IsInNewVertex(i, vertex_collection_int[index])){
// 	      vertex_collection_int[index].push_back(i);
// 	      }
// 	     if(!IsInNewVertex(j, vertex_collection_int[index])){
// 	      vertex_collection_int[index].push_back(j);
// 	     }
// 	    }
// 	}
//       }
//     }


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




    std::vector<std::vector<int> > vertex_collection_int; 

    for (unsigned int i=0; i<trkIn.size(); ++i){
      for (unsigned int j=i+1; j<trkIn.size(); ++j){
	std::cout << "distance between " << i << " and " << j << " = " << StartPointSeperation(startpoints_vec[i], startpoints_vec[j]) << std::endl;
	double GAMMA = gammavalue(startvec[i], startvec[j], dircosvec[i], dircosvec[j]);
	std::cout << "gamma value = " << GAMMA << std::endl; 
	double ALPHA = alphavalue(GAMMA, startvec[i], startvec[j], dircosvec[i], dircosvec[j]);
	std::cout << "alpha value = " << ALPHA << std::endl;
	double MINDIST = MinDist(ALPHA, GAMMA, startvec[i], startvec[j], dircosvec[i], dircosvec[j]);
	std::cout << "MINIMUM DISTANCE = " << MINDIST << std::endl;
	TVector3 TRACK1POINT = PointOnExtendedTrack(ALPHA, startvec[i], dircosvec[i]);
	TVector3 TRACK2POINT = PointOnExtendedTrack(GAMMA, startvec[j], dircosvec[j]);
	std::cout << "POINTS ON THE TRACKS ARE:: " << std::endl;
	TRACK1POINT.Print();
	TRACK2POINT.Print();
	
	if(StartPointSeperation(startpoints_vec[i], startpoints_vec[j])<fVertexWindow){ ///// correct this



	  if((!IsInVertexCollection(i, vertex_collection_int)) && (!IsInVertexCollection(j, vertex_collection_int))){
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


    std::cout<< "*********************************************************" << std::endl;
    
    std::vector< std::pair<art::Ptr<recob::Track>, double> > trackpair;
    for(unsigned int i = 0; i<trkIn.size(); ++i){
      double length = (endvec[i]-startvec[i]).Mag();
      std::cout << "Track length calculated = " << length << std::endl;
      trackpair.push_back(std::pair<art::Ptr<recob::Track>,double>(trkIn[i],length));
    }
    for(unsigned int i = 0; i<trkIn.size(); ++i){
    std::cout << "track id is  = " << (trackpair[i].first)->ID() << std::endl;
    std::cout << "track length = " << (trackpair[i].second) << std::endl;
    }
    std::sort(trackpair.rbegin(), trackpair.rend(),sort_pred2);
    for(unsigned int i = 0; i<trkIn.size(); ++i){
      std::cout << "AFTER SORTING" << std::endl;
    std::cout << "track id is  = " << (trackpair[i].first)->ID() << std::endl;
    std::cout << "track length = " << (trackpair[i].second) << std::endl;
    }
    
    std::cout<< "*********************************************************" << std::endl;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //now add the unmatched track IDs to the collection
    for(int i=0; i<startpoints_vec.size(); i++){
      if(!IsInVertexCollection(i, vertex_collection_int)){
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
     
    // see the veretx collection
    for(unsigned int i=0; i<spvertex_collection.size(); ++i){
      for(std::vector<recob::SpacePoint>::iterator itr = spvertex_collection[i].begin(); itr < spvertex_collection[i].end(); ++itr){
	std::cout << "Start vertex vector = " << *itr << std::endl;
      }
      std::cout << "---------" << std::endl;
    }


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
bool vertex::PrimaryVertexFinder::IsInVertexCollection(int a, std::vector<std::vector<int> > vertex_collection)
{
  int flag = 0;
  
  for(int i = 0; i < vertex_collection.size() ; i++){
    for(std::vector<int>::iterator itr = vertex_collection[i].begin(); itr < vertex_collection[i].end(); ++itr){
      if (a == *itr){
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
double vertex::PrimaryVertexFinder::gammavalue(TVector3 startpoint1, TVector3 startpoint2, TVector3 dircos1, TVector3 dircos2)
{
  double gamma = ((startpoint1*dircos2)-(startpoint2*dircos2)+((dircos1*dircos2)*(startpoint2*dircos1))-((dircos1*dircos2)*(startpoint1*dircos1)))/(1-((dircos1*dircos2)*(dircos1*dircos2)));

  return gamma;
}
// //------------------------------------------------------------------------------
double vertex::PrimaryVertexFinder::alphavalue(double gamma, TVector3 startpoint1, TVector3 startpoint2, TVector3 dircos1, TVector3 dircos2)
{
  double alpha = (gamma*(dircos1*dircos2)) + (startpoint2*dircos1) - (startpoint1*dircos1);

  return alpha;
}
// //------------------------------------------------------------------------------
double vertex::PrimaryVertexFinder::MinDist(double alpha, double gamma, TVector3 startpoint1, TVector3 startpoint2, TVector3 dircos1, TVector3 dircos2)
{
  TVector3 mindis_vector = startpoint1 - startpoint2 + alpha*dircos1 - gamma*dircos2;
  double mindis = mindis_vector.Mag();
  return mindis;
}
// //------------------------------------------------------------------------------
TVector3 vertex::PrimaryVertexFinder::PointOnExtendedTrack(double alphagamma, TVector3 startpoint,  TVector3 dircos)
{
  TVector3 PointOnExtendedTrack = startpoint + (alphagamma * dircos);
  return PointOnExtendedTrack;
}
// //------------------------------------------------------------------------------

