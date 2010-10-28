////////////////////////////////////////////////////////////////////////
//
// ShowerFinder class
//
// roxanne.guenette@yale.edu
//
//  This algorithm is designed to find all the showers in an event.
////////////////////////////////////////////////////////////////////////

#include "ShowerFinder/ShowerFinder.h"

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}
#include <sstream>
#include <math.h>
#include <algorithm>

// Framework includes
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

// LArSoft Includes
#include "RawData/RawDigit.h"
//#include "HitFinder/FFTHitFinder.h"
#include "ClusterFinder/DBcluster.h"
#include "ClusterFinder/HoughLineFinder.h"
#include "VertexFinder/VertexMatch.h"
#include "Geometry/geo.h"
#include "RecoBase/recobase.h"
//#include "Utilities/LArFFT.h"

// ROOT 
#include "TMath.h"
#include "TF1.h"
#include "TH1D.h"


namespace shwf{

  //-------------------------------------------------
  ShowerFinder::ShowerFinder(edm::ParameterSet const& pset) :
    fvertices (pset.getParameter<std::string > ("./vertices")),
    fclusters (pset.getParameter<std::string > ("./clusters")),
    fhoughlines (pset.getParameter<std::string > ("./houghliness")),
    fhits (pset.getParameter<std::string > ("./hits")),
    fRcone (pset.getParameter<double > ("Rcone")), //Radius of the cone
    fLcone (pset.getParameter<double > ("Lcone"))  // Length (perpendicular to the base) of the cone
    
  {
    produces< std::vector<recob::Shower> >();
  }
  
  
  //-------------------------------------------------

  ShowerFinder::~ShowerFinder()
  {
  }
  
  //-------------------------------------------------
  void ShowerFinder::beginJob(const edm::EventSetup&)
  {
  }
  void ShowerFinder::endJob()
  {
  }

  //  
  //-------------------------------------------------
  void ShowerFinder::produce(edm::Event& evt, edm::EventSetup const&)
  { 
    
    
    //////////////////////////////////////////////////////
    // Make a std::auto_ptr<> for the thing you want to put into the event
    // because that handles the memory management for you
    //////////////////////////////////////////////////////
   
    //std::auto_ptr<std::vector<recob::Shower> > show_ID(new std::vector<recob::Shower>);
    //std::auto_ptr<std::vector<recob::Shower> > show_vcoord(new std::vector<recob::Shower>);
    std::auto_ptr<std::vector<recob::Shower> > showercol(new std::vector<recob::Shower>);


    // Read in the vertex List object(s).
    edm::Handle< std::vector<recob::Vertex> > vertexListHandle;
    evt.getByLabel(fvertices,vertexListHandle);
    // Read in the hough line List object(s).
    edm::Handle< std::vector<recob::Cluster> > houghListHandle;
    evt.getByLabel(fhoughlines,houghListHandle);
    // Read in the cluster List object(s).
    edm::Handle< std::vector<recob::Cluster> > clusterListHandle;
    evt.getByLabel(fclusters,clusterListHandle);
    // Read in the hit List object(s).
    edm::Handle< std::vector<recob::Hit> > hitcol;
    evt.getByLabel(fhits,hitcol);
    
    edm::PtrVector<recob::Cluster> protoShowers; //vector of clusters associated to a cone
        edm::PtrVector<recob::Hit> sHits;//hits associated to showers

    recob::Shower shower(protoShowers);


    edm::PtrVector<recob::Hit> vertexhits; //hits in the vertex
    edm::PtrVector<recob::Hit> clusterhits; //hits in the cluster
    edm::PtrVector<recob::Hit> hlhits; //hits in the hough Lines

 
    edm::Service<geo::Geometry> geom;
    
    unsigned int channel,plane,wire;
    
    //This vector will contain all strong and strongest vertices
    edm::PtrVector<recob::Vertex> vertSel;
    
    //This loop is going over all the vertices in the event 
    //and be interested in ONLY strong and strongest vertices.

    for(unsigned int iv = 0; iv < vertexListHandle->size(); ++iv)
      {
	edm::Ptr<recob::Vertex> vertex(vertexListHandle, iv);
	if(vertex->Strength() == 4 || vertex->Strength() == 3) vertSel.push_back(vertex);
	else continue;
      }

    //Definition of the geometry of the cone (which is basically a triangle)
    double scan_angle = 0; //angle of the scan steps
    double xa_cone = 0; // x coordinate of the cone's apex (wire number)
    double ya_cone = 0; // y coordinate of the cone's apex (drift time)
    double x1_cone = 0; // x coordinate of the cone's top left point (wire number)
    double y1_cone = 0; // y coordinate of the cone's top left point (drift time)
    double x2_cone = 0; // x coordinate of the cone's top right point (wire number)
    double y2_cone = 0; // y coordinate of the cone's top right point  (drift time)

    double fScone = sqrt( (fRcone*fRcone) + (fLcone*fLcone)); //The length of the side of the cone

    double cone_angle = TMath::ATan(fRcone / fLcone); // Opening angle of the cone (defined from input parameters)
    double compl_angle = 0;  

    int n_scan =(int)(TMath::Pi() / cone_angle); 

    double x_hit = 0; //x coordinate of hit
    double y_hit = 0; //y coordinate of hit

    int hits_cluster_counter = 0; //count the number of hits in a cluster that is inside a cone
    int hits_cluster_Total = 0; //The total number of hits in a cluster  
 
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //TO DO FOR EVERY PLANE!!!!!!!!!!!!!
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    //For EVERY vertex, the algorithm is going to scan the plane to find clusters contained in the scanning cones
    for(unsigned int ivert = 0; ivert < vertSel.size(); ++ivert)
      {
	//get the coordinates of the vertex for the summit of the cone
	xa_cone = vertSel[ivert]->Wire();
	ya_cone = vertSel[ivert]->DriftTime();
	  

	//Beginning of the scan!
	for(int iscan = 0; iscan < n_scan; iscan++){

	  //define the scan anlge
	  scan_angle = (TMath::Pi()/2.0) - (iscan*cone_angle);
	  
	  //get the complementary angle for geometry puurposes
	  compl_angle = (TMath::Pi()/2.0) - scan_angle - cone_angle;  
	  
	  //Calculate the coordinates of thre top left corner of the cone
	  x1_cone = xa_cone + fScone*(TMath::Sin(compl_angle));
	  y1_cone = ya_cone + fScone*(TMath::Cos(compl_angle));
	  
	  //Calculate the coordinates of thre top right corner of the cone
	  x2_cone = xa_cone + fScone*(TMath::Sin(scan_angle - cone_angle));
	  y2_cone = ya_cone + fScone*(TMath::Cos(scan_angle - cone_angle));



	  //Looking if a cluster is in this cone (loop over all hits of all clusters)
	  for(int iclust = 0; iclust < clusterListHandle->size(); iclust++){
	    edm::Ptr<recob::Cluster> clust(clusterListHandle, iclust);
	    
	    //GEt the hits vector from the cluster?????? To correct!!!
	    clusterhits = clust->Hits(1, -1);

	    //Loop over ALL hits in the cluster. Looking if the cluster's hit is comprised in the cone
	    for(int ihits = 0; ihits < clusterhits.size(); ihits++){

	      x_hit = clusterhits[ihits]->Wire()->RawDigit()->Channel();
	      y_hit = clusterhits[ihits]->CrossingTime();

	      
	      // Check in hits is INSIDE cone
	      if((x_hit >= xa_cone && x_hit <= x1_cone) && (x_hit >= xa_cone && x_hit <= x2_cone) && (y_hit >= ya_cone && y_hit <= y1_cone) && (y_hit >= ya_cone && y_hit <= y2_cone)) hits_cluster_counter++;


	    }//end hits loop

	    //If there is more than 50% if the cluster INSIDE the cone, this is a protoshower
	    if((hits_cluster_counter / clusterhits.size()) >= 0.5){

	      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	      //NEED TO TAKE OUT THE HOUGH LINES FROM THE PROTOSHOWERS!!!!!  
	      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	      //pshower_clust.push_back(clust); //this vector contains all the cluster in the cone proto-shower
	      protoShowers.push_back(clust); //this vector contains all the cluster in the cone proto-shower
	    }


	  } //end cluster loop
	  
	} //end scan loop
	
      } //end vertices loop

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //NEED TO SEPARATE THE SHOWERS FROM THE DIFFERENT VERTEX!!!!  
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    for(int is = protoShowers.size(); is > 0; is--){
      //check if protoshower from further vertex is also contained in vertex nearer... TODO!!!!
      //if the shower is stand alone ok, else, erase the next one
      //shower.SetID(is);
      //shower.SetVertexCoord(xa_cone, ya_cone);

      showercol->push_back(shower);
    }

    evt.put(showercol);

  } // end of produce
} // end of namespace
