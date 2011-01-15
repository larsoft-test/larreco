////////////////////////////////////////////////////////////////////////
//
// Track3Dreco class
//
// maddalena.antonello@lngs.infn.it
// ornella.palamara@lngs.infn.it
// ART port and edits: soderber,echurch@fnal.gov
//  This algorithm is designed to reconstruct 3D tracks through a simple 
//  2D-track matching algorithm
////////////////////////////////////////////////////////////////////////

// C++ includes
#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>

// Framework includes
#include "FWCore/Framework/interface/Event.h" 
#include "FWCore/ParameterSet/interface/ParameterSet.h" 
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h" 
#include "DataFormats/Common/interface/Handle.h" 
#include "DataFormats/Common/interface/View.h" 
#include "DataFormats/Common/interface/Ptr.h" 
#include "DataFormats/Common/interface/PtrVector.h" 
#include "FWCore/Framework/interface/MakerMacros.h" 
#include "FWCore/ServiceRegistry/interface/Service.h" 
#include "FWCore/Services/interface/TFileService.h" 
#include "FWCore/Framework/interface/TFileDirectory.h" 
#include "FWCore/MessageLogger/interface/MessageLogger.h" 
#include "FWCore/ServiceRegistry/interface/ServiceMaker.h" 

// LArSoft includes
#include "Track3Dreco.h"
#include "Geometry/geo.h"
#include "Geometry/WireGeo.h"
#include "Geometry/VolumeUtility.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Track.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/SpacePoint.h"
#include "Utilities/LArProperties.h"

// ROOT includes
#include "TVectorD.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMath.h"


//-------------------------------------------------
trkf::Track3Dreco::Track3Dreco(edm::ParameterSet const& pset) :
  fClusterModuleLabel     (pset.getParameter< std::string >("ClusterModuleLabel")),
  ftmatch                 (pset.getParameter< int    >("TMatch")),
  fchi2dof                (pset.getParameter< double >("Chi2DOFmax"))
{
  produces< std::vector<recob::Track> >();
}

//-------------------------------------------------
trkf::Track3Dreco::~Track3Dreco()
{
}

//-------------------------------------------------
void trkf::Track3Dreco::beginJob(const edm::EventSetup&)
{
  

}
void trkf::Track3Dreco::endJob()
{
}

//------------------------------------------------------------------------------------//
void trkf::Track3Dreco::produce(edm::Event& evt, edm::EventSetup const&)
{ 


  // get services
  edm::Service<geo::Geometry> geom;
  edm::Service<util::LArProperties> larprop;

  //////////////////////////////////////////////////////
  // Make a std::auto_ptr<> for the thing you want to put into the event
  // because that handles the memory management for you
  //////////////////////////////////////////////////////
  std::auto_ptr<std::vector<recob::Track> > tcol(new std::vector<recob::Track>);

  // define TPC parameters
  TString tpcName = geom->GetLArTPCVolumeName();

  //TPC dimensions
  //double m_TPCHalfZ = m_tpcVolumeUtility->GetHalfZ();
  double m_TPCHalfZ = geom->DetLength()-5.0;

  //  double YC =  (m_TPCHalfZ-5.)*2.; // TPC height in cm
  double YC =  (geom->DetHalfHeight()-0.5715)*2.; // TPC height in cm
  double Angle = geom->Plane(1).Wire(0).ThetaZ(false)-TMath::Pi()/2.; // wire angle with respect to the vertical direction
  // Parameters temporary defined here, but possibly to be retrieved somewhere in the code
  double timetick = 0.198;    //time sample in us
  double presamplings = 60.;
  const double wireShift=50.; // half the number of wires from the Induction(Collection) plane intersecting with a wire from the Collection(Induction) plane.
  double plane_pitch = geom->PlanePitch(0,1);   //wire plane pitch in cm 
  double wire_pitch = geom->WirePitch(0,1,0);    //wire pitch in cm
  double Efield_drift = 0.5;  // Electric Field in the drift region in kV/cm
  double Efield_SI = 0.7;     // Electric Field between Shield and Induction planes in kV/cm
  double Efield_IC = 0.9;     // Electric Field between Induction and Collection planes in kV/cm
  double Temperature = 87.6;  // LAr Temperature in K

  double driftvelocity = larprop->DriftVelocity(Efield_drift,Temperature);    //drift velocity in the drift region (cm/us)
  double driftvelocity_SI = larprop->DriftVelocity(Efield_SI,Temperature);    //drift velocity between shield and induction (cm/us)
  double driftvelocity_IC = larprop->DriftVelocity(Efield_IC,Temperature);    //drift velocity between induction and collection (cm/us)
  double timepitch = driftvelocity*timetick;                         //time sample (cm) 
  double tSI = plane_pitch/driftvelocity_SI/timetick;                   //drift time between Shield and Collection planes (time samples)
  double tIC = plane_pitch/driftvelocity_IC/timetick;                //drift time between Induction and Collection planes (time samples)


  // get input Cluster object(s).
  edm::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fClusterModuleLabel,clusterListHandle);

  
  // Declare some vectors..
  // Induction
  std::vector<double> Iwirefirsts;       // in cm
  std::vector<double> Iwirelasts;        // in cm
  std::vector<double> Itimefirsts;       // in cm
  std::vector<double> Itimelasts;        // in cm
  std::vector<double> Itimefirsts_line;  // in cm
  std::vector<double> Itimelasts_line;   // in cm
  std::vector < edm::PtrVector<recob::Hit> > IclusHitlists;
  std::vector<unsigned int> Icluster_count; 

  // Collection
  std::vector<double> Cwirefirsts;       // in cm
  std::vector<double> Cwirelasts;        // in cm
  std::vector<double> Ctimefirsts;       // in cm
  std::vector<double> Ctimelasts;        // in cm
  std::vector<double> Ctimefirsts_line;  // in cm
  std::vector<double> Ctimelasts_line;   // in cm
  std::vector< edm::PtrVector < recob::Hit> > CclusHitlists;
  std::vector<unsigned int> Ccluster_count; 

  for(unsigned int ii = 0; ii < clusterListHandle->size(); ++ii)
    {
      edm::Ptr<recob::Cluster> cl(clusterListHandle, ii);
      
      /////////////////////////
      //////////// 2D track FIT
      /////////////////////////
      
      // Figure out which View the cluster belongs to 
      int clPlane = cl->View()-1;

      // Some variables for the hit
      unsigned int channel;  //channel number
      float time;            //hit time at maximum
      unsigned int wire;              //hit wire number
      unsigned int plane;             //hit plane number
      
      
      edm::PtrVector<recob::Hit> hitlist;
      hitlist = cl->Hits( clPlane, -1);
      // std::sort(hitlist.begin(), hitlist.end(), hit_sort_2d); //sort hit by wire
      
      TGraph *the2Dtrack = new TGraph(hitlist.size());
      
      std::vector<double> wires;
      std::vector<double> times;
     
      
      int np=0;
      for(edm::PtrVectorItr<recob::Hit> theHit = hitlist.begin(); theHit != hitlist.end();  theHit++) //loop over cluster hits
	{
	  //recover the Hit
	  //      recob::Hit* theHit = (recob::Hit*)(*hitIter);
	  time = (*theHit)->CrossingTime() ;
	 
	  time -= presamplings;
	  
	  channel = (*theHit)->Channel();
	  geom->ChannelToWire(channel,plane,wire);
	  
	  //correct for the distance between wire planes
	  if(plane==0) time -= tSI;         // Induction
	  if(plane==1) time -= (tSI+tIC);   // Collection
	
	  //transform hit wire and time into cm
	  double wire_cm = (double)((wire+1) * wire_pitch); 
	  double time_cm = (double)(time * timepitch);
	  wires.push_back(wire_cm);
	  times.push_back(time_cm);
	
	  the2Dtrack->SetPoint(np,wire_cm,time_cm);
	  np++;
	}//end of loop over cluster hits
    
      // fit the 2Dtrack and get some info to store
      the2Dtrack->Fit("pol1","Q");
      TF1 *pol1=(TF1*) the2Dtrack->GetFunction("pol1");
      double Chi2 = pol1->GetChisquare();
      double NDF = pol1->GetNDF();
      double par[2];
      pol1->GetParameters(par);
      double intercept = par[0];
      double slope = par[1];
  
      double w0 = wires.front();      // first hit wire (cm)
      double w1 = wires.back();        // last hit wire (cm)
      double t0 = times.front();      // first hit time (cm)
      double t1 = times.back();        // last hit time (cm)
      double t0_line = intercept + (w0)*slope;       // time coordinate at wire w0 on the fit line (cm)  
      double t1_line = intercept + (w1)*slope;       // time coordinate at wire w1 on the fit line (cm)
    


      // actually store the 2Dtrack info
      switch(plane){
      case 0:
	Iwirefirsts.push_back(w0);
	Iwirelasts.push_back(w1);
	Itimefirsts.push_back(t0);
	Itimelasts.push_back(t1); 
	Itimefirsts_line.push_back(t0_line);
	Itimelasts_line.push_back(t1_line);    
	IclusHitlists.push_back(hitlist);
	Icluster_count.push_back(ii);
	break;
      case 1:
	Cwirefirsts.push_back(w0);
	Cwirelasts.push_back(w1);
	Ctimefirsts.push_back(t0);
	Ctimelasts.push_back(t1);
	Ctimefirsts_line.push_back(t0_line);
	Ctimelasts_line.push_back(t1_line);
	CclusHitlists.push_back(hitlist);
	Ccluster_count.push_back(ii);
	break;   
      }
      //delete the2Dtrack;
      delete pol1;
    }// end of loop over all input clusters
  
  /////////////////////////////////////////////////////
  /////// 2D Track Matching and 3D Track Reconstruction
  /////////////////////////////////////////////////////

  for(unsigned int collectionIter=0; collectionIter < CclusHitlists.size();collectionIter++){  //loop over Collection view 2D tracks
    // Recover previously stored info
    double Cw0 = Cwirefirsts[collectionIter];
    double Cw1 = Cwirelasts[collectionIter];
    double Ct0 = Ctimefirsts[collectionIter];
    double Ct1 = Ctimelasts[collectionIter];
    double Ct0_line = Ctimefirsts_line[collectionIter];
    double Ct1_line = Ctimelasts_line[collectionIter];
    edm::PtrVector<recob::Hit> hitsCtrk = CclusHitlists[collectionIter];

    for(unsigned int inductionIter=0;inductionIter<IclusHitlists.size();inductionIter++){   //loop over Induction view 2D tracks
      // Recover previously stored info
      double Iw0 = Iwirefirsts[inductionIter];
      double Iw1 = Iwirelasts[inductionIter];
      double It0 = Itimefirsts[inductionIter];
      double It1 = Itimelasts[inductionIter];
      double It0_line = Itimefirsts_line[inductionIter];
      double It1_line = Itimelasts_line[inductionIter];
      edm::PtrVector<recob::Hit> hitsItrk = IclusHitlists[inductionIter];

      // match 2D tracks
      if((fabs(Ct0_line-It0_line)<ftmatch*timepitch) && (fabs(Ct1_line-It1_line)<ftmatch*timepitch)){ 
	std::cout<<"-----> Track "<<collectionIter<< " Collection associated with track "<<inductionIter<< " Induction"<<std::endl;
	
       
        // Reconstruct the 3D track
	TVector3 XYZ0;  // track origin or interaction vertex
	XYZ0.SetXYZ(Ct0_line,(Cw0-Iw0)/(2.*TMath::Sin(Angle)),(Cw0+Iw0)/(2.*TMath::Cos(Angle))-YC/2.*TMath::Tan(Angle));

	//compute track startpoint and endpoint in Local co-ordinate system 
	TVector3 startpointVec(XYZ0.X(),XYZ0.Y(),XYZ0.Z());
 	TVector3 endpointVec(Ct1_line,(Cw1-Iw1)/(2.*TMath::Sin(Angle)),(Cw1+Iw1)/(2.*TMath::Cos(Angle))-YC/2.*TMath::Tan(Angle));
	//	TVector3 startpointVecLocal = m_tpcVolumeUtility->WorldToLocal(startpointVec);
	TVector3 startpointVecLocal = geom->Plane(1).WorldToLocal(startpointVec);
	//TVector3 endpointVecLocal = m_tpcVolumeUtility->WorldToLocal(endpointVec);
	TVector3 endpointVecLocal = geom->Plane(1).WorldToLocal(endpointVec);
	Double_t startpoint[3],endpoint[3];
	startpoint[0] = startpointVecLocal[0];
	startpoint[1] = startpointVecLocal[1];
	startpoint[2] = startpointVecLocal[2];
	endpoint[0]   = endpointVecLocal[0];
	endpoint[1]   = endpointVecLocal[1];
	endpoint[2]   = endpointVecLocal[2]; 

	//compute track (normalized) cosine directions in the World co-ordinate system
	TVector3 DirCos = endpointVecLocal - startpointVecLocal;
	DirCos.SetMag(1.0);//normalize vector

	//compute 3D track parameters (theta, phi, track pitch length) in the World co-ordinate system         
	double Theta; // Zenith angle with respect to y (vertical) axis in radians 
	double Phi;   // Azimuth angle with respect to x axis in radians
	Theta = TMath::ACos(DirCos.Y());
	Phi = TMath::ATan(DirCos.Z()/DirCos.X());
	Phi = Phi<0. ? Phi+TMath::Pi() : Phi ; // solve the ambiguities due to tangent periodicity	
	double cosgammaC; // angle between the track direction and the direction of the Collection wire pitch
	cosgammaC = TMath::Sin(Theta)*TMath::Sin(Phi)*TMath::Cos(Angle)+TMath::Cos(Theta)*TMath::Sin(Angle);
	double TrackPitchC  = wire_pitch/cosgammaC; // Collection Track Pitch length (i.e. the effective length of the track seen by the Collection wire) in cm  	
	double cosgammaI; // angle between the track direction and the direction of the Induction wire pitch
	cosgammaI = TMath::Sin(Theta)*TMath::Sin(Phi)*TMath::Cos(Angle)-TMath::Cos(Theta)*TMath::Sin(Angle);
	double TrackPitchI  = wire_pitch/cosgammaI; // Induction Track Pitch length (i.e. the effective length of the track seen by the Induction wire) in cm
	//double TrackLength_line = TMath::Sqrt(TMath::Power((startpoint[0]-endpoint[0]),2.)+TMath::Power((startpoint[1]-endpoint[1]),2.)+TMath::Power((startpoint[2]-endpoint[2]),2.));

	edm::Ptr <recob::Cluster> cl1(clusterListHandle,Icluster_count[inductionIter]);
	edm::Ptr <recob::Cluster> cl2(clusterListHandle,Ccluster_count[collectionIter]);
	edm::PtrVector<recob::Cluster> clustersPerTrack;
	clustersPerTrack.push_back(cl1);
	clustersPerTrack.push_back(cl2);

 
	/////////////////////////////
	// Match hits
	////////////////////////////

	//create collection of spacepoints that will be used when creating the Track object
	std::vector<recob::SpacePoint> spacepoints;
	

	edm::PtrVector<recob::Hit> minhits = hitsCtrk.size() <= hitsItrk.size() ? hitsCtrk : hitsItrk;
	edm::PtrVector<recob::Hit> maxhits = hitsItrk.size() <= hitsCtrk.size() ? hitsCtrk : hitsItrk;


	bool maxhitsMatch[maxhits.size()];
	for(unsigned int it=0;it<maxhits.size();it++) maxhitsMatch[it] = false;

	std::vector<recob::Hit*> hits3Dmatched;
	// For the matching start from the view where the track projection presents less hits
	unsigned int imaximum = 0;
	for(unsigned int imin=0;imin<minhits.size();imin++){ //loop over hits
	  //get wire - time coordinate of the hit
	  unsigned int channel,wire,plane1,plane2;
	  channel = minhits[imin]->Channel();
	  geom->ChannelToWire(channel,plane1,wire);
	  // get the wire-time co-ordinates of the hit to be matched
	  double w1 = (double)((wire+1)*wire_pitch);
 	  double t1 = plane1==1?(double)((minhits[imin]->CrossingTime()-presamplings-(tSI+tIC))*timepitch):(double)((minhits[imin]->CrossingTime()-presamplings-tSI)*timepitch); //in cm	  

	  //get the track origin co-ordinates in the two views
	  TVector2 minVtx2D;
	  plane1==1 ? minVtx2D.Set(Ct0,Cw0): minVtx2D.Set(It0,Iw0);
	  TVector2 maxVtx2D;
	  plane1==1 ? maxVtx2D.Set(It0,Iw0): maxVtx2D.Set(Ct0,Cw0);
	  
	  // get the track last hit co-ordinates in the two views
	  channel = minhits[minhits.size()-1]->Channel();
	  geom->ChannelToWire(channel,plane1,wire);
	  double w1_last = plane1==1?Cw1:Iw1;
	  double t1_last = plane1==1?Ct1:It1;  
	  channel = maxhits[maxhits.size()-1]->Channel();
	  geom->ChannelToWire(channel,plane2,wire);
	  double w2_last =plane2==1?Cw1:Iw1;
	  double t2_last =plane2==1?Ct1:It1;
	  
	  // compute the track length in the two views
	  double minLength = TMath::Sqrt(TMath::Power(t1_last-minVtx2D.X(),2)+TMath::Power(w1_last-minVtx2D.Y(),2));
	  double maxLength = TMath::Sqrt(TMath::Power(t2_last-maxVtx2D.X(),2)+TMath::Power(w2_last-maxVtx2D.Y(),2));
	  
	  //compute the distance of the hit (imin) from the relative track origin
	  double minDistance = maxLength* TMath::Sqrt(TMath::Power(t1-minVtx2D.X(),2)+TMath::Power(w1-minVtx2D.Y(),2))/minLength;
	  
	  //core matching algorithm
	  double difference = 9999999.;	  

	  for(unsigned int imax = imaximum; imax < maxhits.size(); imax++){ //loop over hits of the other view
	    if(!maxhitsMatch[imax]){
	      //get wire - time coordinate of the hit
	      channel = maxhits[imax]->Channel();
	      geom->ChannelToWire(channel,plane2,wire);
	      double w2 = (double)((wire+1)*wire_pitch);
	      double t2 = plane2==1?(double)((maxhits[imax]->CrossingTime()-presamplings-(tSI+tIC))*timepitch):(double)((maxhits[imax]->CrossingTime()-presamplings-tSI)*timepitch); //in cm
	      
	      bool timematch = (fabs(t1-t2)<ftmatch*timepitch);
	      bool wirematch = (fabs(w1-w2)<wireShift*wire_pitch);
	      
	      double maxDistance = TMath::Sqrt(TMath::Power(t2-maxVtx2D.X(),2)+TMath::Power(w2-maxVtx2D.Y(),2));
	      if (wirematch && timematch && fabs(maxDistance-minDistance)<difference) {
		difference = fabs(maxDistance-minDistance);
		imaximum = imax;
	      }
	    }
	  }
	  maxhitsMatch[imaximum]=true;
	  
	  edm::PtrVector<recob::Hit> sp_hits;
	  if(difference!= 9999999.){
	    sp_hits.push_back(minhits[imin]);
	    sp_hits.push_back(maxhits[imaximum]);
	  }
	  
	  // Get the time-wire co-ordinates of the matched hit
	  channel =  maxhits[imaximum]->Channel();
	  geom->ChannelToWire(channel,plane2,wire);
	  double w1_match = (double)((wire+1)*wire_pitch);
	  double t1_match = plane2==1?(double)((maxhits[imaximum]->CrossingTime()-presamplings-(tSI+tIC))*timepitch):(double)((maxhits[imaximum]->CrossingTime()-presamplings-tSI)*timepitch);

	  // create the 3D hit, compute its co-ordinates and add it to the 3D hits list	  
	  double Ct = plane1==1?t1:t1_match;
	  double Cw = plane1==1?w1:w1_match;
	  double Iw = plane1==1?w1_match:w1;

	  const TVector3 hit3d(Ct,(Cw-Iw)/(2.*TMath::Sin(Angle)),(Cw+Iw)/(2.*TMath::Cos(Angle))-YC/2.*TMath::Tan(Angle)); 
	  const TVector3 hit3dLocal = geom->Plane(plane1).WorldToLocal(hit3d);// m_tpcVolumeUtility->WorldToLocal(hit3d);
	  Double_t hitcoord[3];
	  hitcoord[0] = hit3dLocal.X();
	  hitcoord[1] = hit3dLocal.Y();
	  hitcoord[2] = hit3dLocal.Z();           

	  recob::SpacePoint mysp(sp_hits);//3d point at end of track
	  mysp.SetXYZ(hitcoord);
	  spacepoints.push_back(mysp);
	  
	}//loop over min-hits
      

	// Add the 3D track to the vector of the reconstructed tracks
        if(spacepoints.size()>0 || clustersPerTrack.size()>0){
	 
	  recob::Track  the3DTrack(clustersPerTrack,spacepoints);
	  double dircos[3];
	  DirCos.GetXYZ(dircos);
	  the3DTrack.SetDirection(dircos,dircos);
	  the3DTrack.SetTrackPitch(TrackPitchI,geo::kU);
	  the3DTrack.SetTrackPitch(TrackPitchC,geo::kV);

	  tcol->push_back(the3DTrack);
	}
      } //close match 2D tracks


    }//close loop over Induction view 2D tracks
    
  }//close loop over Collection xxview 2D tracks
  
  std::cout<<"Run "<<evt.run()<<" Event "<<evt.id().event()<<std::endl;
  std::cout<<"TRACK3DRECO found "<< tcol->size() <<" 3D track(s)"<<std::endl;
  
  evt.put(tcol);
  

}
