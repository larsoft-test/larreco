////////////////////////////////////////////////////////////////////////
//
// Track3Dreco class
//
// maddalena.antonello@lngs.infn.it
// ornella.palamara@lngs.infn.it
// ART port and edits: echurch@fnal.gov
//  This algorithm is designed to reconstruct 3D tracks through a simple 
//  2D-track matching algorithm
////////////////////////////////////////////////////////////////////////

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

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
#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include "Track3Dreco.h"
#include "Geometry/geo.h"
#include "Geometry/WireGeo.h"
#include "Geometry/VolumeUtility.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Track.h"
#include "RecoBase/Cluster.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompSVD.h"
#include "TH2F.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMath.h"


static bool hit_sort_2d(const recob::Hit* h1, const recob::Hit* h2)
{
  return h1->Wire()->RawDigit()->Channel() < h2->Wire()->RawDigit()->Channel();
}

//-------------------------------------------------
trkf::Track3Dreco::Track3Dreco(edm::ParameterSet const& pset) :
  fDBScanModuleLabel     (pset.getParameter< std::string >("DBScanModuleLabel")),
  ftmatch                (pset.getParameter< int >("TMatch"))
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


  // get the geometry
  edm::Service<geo::Geometry> geom;

  //Get Clusters
  // Read in the clusterList object(s).
  edm::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fDBScanModuleLabel,clusterListHandle);
  std::auto_ptr<std::vector<recob::Track> > tcol(new std::vector<recob::Track>);
  edm::PtrVector<recob::Cluster> clusterlist;
  for(unsigned int ii = 0; ii < clusterListHandle->size(); ++ii)
    {
      edm::Ptr<recob::Cluster> cluster(clusterListHandle, ii);
      clusterlist.push_back(cluster);
    }


  // TPC parameters
  TString tpcName = geom->GetLArTPCVolumeName();
  geo::VolumeUtility* m_tpcVolumeUtility = new geo::VolumeUtility( tpcName );
  //TPC dimensions
  double m_TPCHalfZ = m_tpcVolumeUtility->GetHalfZ();
  double YC =  (m_TPCHalfZ-5.)*2.; // TPC height in cm
  double Angle = geom->Plane(1).Wire(0).ThetaZ(false)-TMath::Pi()/2.; // wire angle with respect to the vertical direction
  // Parameters temporary defined here, but possibly to be retrieved somewhere in the code
  double timetick = 0.198;    //time sample in us
  double presamplings = 60.;
  const double wireShift=50.; // half the number of wires from the Induction(Collection) plane intersecting with a wire from the Collection(Induction) plane.
  double plane_pitch = 0.4;   //wire plane pitch in cm 
  double wire_pitch = 0.4;    //wire pitch in cm
  //double Efield_drift = 0.52;  // Electric Field in the drift region in kV/cm
  double Efield_SI = 0.7;     // Electric Field between Shield and Induction planes in kV/cm
  double Efield_IC = 0.9;     // Electric Field between Induction and Collection planes in kV/cm
  double Temperature = 87.6;  // LAr Temperature in K
  //  double driftvelocity = DriftVelocity(Efield_drift,Temperature);    //drift velocity in the drift region (cm/us)
  double driftvelocity = 0.153;
  double driftvelocity_SI = DriftVelocity(Efield_SI,Temperature);    //drift velocity between shield and induction (cm/us)
  double driftvelocity_IC = DriftVelocity(Efield_IC,Temperature);    //drift velocity between induction and collection (cm/us)
  double timepitch = driftvelocity*timetick;                         //time sample (cm) 
  double tSI = plane_pitch/driftvelocity_SI/timetick;                   //drift time between Shield and Collection planes (time samples)
  double tIC = plane_pitch/driftvelocity_IC/timetick;                //drift time between Induction and Collection planes (time samples)


  // Declare some vectors..
  // Induction
  std::vector<double> Islopes; 
  std::vector<double> Iintercepts;       // in cm
  std::vector<double> Ichi2s;
  std::vector<double> Indfs;
  std::vector<double> InormChi2s;  
  std::vector<double> Icrosstfirsts;     // in samples
  std::vector<double> Icrosstlasts;      // in samples
  std::vector<double> Iwirefirsts;       // in cm
  std::vector<double> Iwirelasts;        // in cm
  std::vector<double> Itimefirsts;       // in cm
  std::vector<double> Itimelasts;        // in cm
  std::vector<double> Itimefirsts_line;  // in cm
  std::vector<double> Itimelasts_line;   // in cm
  std::vector<int>    IclusterIDs;
  std::vector < edm::PtrVector<recob::Hit> > IclusHitlists;
  // Collection
  std::vector<double> Cslopes;
  std::vector<double> Cintercepts;       // in cm
  std::vector<double> Cchi2s;
  std::vector<double> Cndfs;
  std::vector<double> CnormChi2s;
  std::vector<double> Ccrosstfirsts;     // in samples
  std::vector<double> Ccrosstlasts;      // in samples
  std::vector<double> Cwirefirsts;       // in cm
  std::vector<double> Cwirelasts;        // in cm
  std::vector<double> Ctimefirsts;       // in cm
  std::vector<double> Ctimelasts;        // in cm
  std::vector<double> Ctimefirsts_line;  // in cm
  std::vector<double> Ctimelasts_line;   // in cm
  std::vector<int>    CclusterIDs;
  std::vector< edm::PtrVector < recob::Hit> > CclusHitlists;

  
  /////////////////////////
  //////////// 2D track FIT
  /////////////////////////

  TGraph* the2Dtrack=0;
  edm::PtrVectorItr<recob::Cluster> cl = clusterlist.begin();
  for(cl; cl != clusterlist.end(); cl++) 
    {
      
      // Figure out which View the cluster belongs to 
      int clPlane = (*cl)->View()-1;

      // Some variables for the hit
      unsigned int channel;  //channel number
      float time;            //hit time at maximum
      unsigned int wire;              //hit wire number
      unsigned int plane;             //hit plane number
   

      edm::PtrVector<recob::Hit> hitlist;
      hitlist = (*cl)->Hits( clPlane, -1);
      std::sort(hitlist.begin(), hitlist.end(), hit_sort_2d); //sort hit by wire
    
    the2Dtrack = new TGraph(hitlist.size());

    std::vector<double> wires;
    std::vector<double> times;
    std::vector<double> crossingtimes;

    int np=0;
    for(edm::PtrVectorItr<recob::Hit> theHit = hitlist.begin(); theHit != hitlist.end();  theHit++) //loop over cluster hits
      {
      //recover the Hit
	//      recob::Hit* theHit = (recob::Hit*)(*hitIter);
	time = (*theHit)->CrossingTime() ;
      crossingtimes.push_back(time);
      time -= presamplings;
      edm::Ptr <recob::Wire> theWire = (*theHit)->Wire();
      channel = theWire->RawDigit()->Channel();
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
    int nw = wires.size();
    double w0 = wires[0];      // first hit wire (cm)
    double w1 = wires[nw-1];   // last hit wire (cm)
    double t0 = times[0];      // first hit time (cm)
    double t1 = times[nw-1];   // last hit time (cm)
    double t0_line = intercept + (w0)*slope; // time coordinate at wire w0 on the fit line (cm)  
    double t1_line = intercept + (w1)*slope; // time coordinate at wire w1 on the fit line (cm)
    double crosst0 = crossingtimes[0];       // hit crossing time at wire w0 (samples)
    double crosst1 = crossingtimes[nw-1];    // hit crossing time at wire w1 (samples)


    // actually store the 2Dtrack info
    switch(plane){
    case 0:
      Iintercepts.push_back(intercept);
      Islopes.push_back(slope);
      Ichi2s.push_back(Chi2);
      Indfs.push_back(NDF);
      InormChi2s.push_back(Chi2/NDF);
      Icrosstfirsts.push_back(crosst0);
      Icrosstlasts.push_back(crosst1);
      Iwirefirsts.push_back(w0);
      Iwirelasts.push_back(w1);
      Itimefirsts.push_back(t0);
      Itimelasts.push_back(t1); 
      Itimefirsts_line.push_back(t0_line);
      Itimelasts_line.push_back(t1_line);    
      IclusterIDs.push_back((*cl)->ID());
      IclusHitlists.push_back(hitlist);
      break;
    case 1:
      Cintercepts.push_back(intercept);
      Cslopes.push_back(slope);
      Cchi2s.push_back(Chi2);
      Cndfs.push_back(NDF);
      CnormChi2s.push_back(Chi2/NDF);
      Ccrosstfirsts.push_back(crosst0);
      Ccrosstlasts.push_back(crosst1);
      Cwirefirsts.push_back(w0);
      Cwirelasts.push_back(w1);
      Ctimefirsts.push_back(t0);
      Ctimelasts.push_back(t1);
      Ctimefirsts_line.push_back(t0_line);
      Ctimelasts_line.push_back(t1_line);
      CclusterIDs.push_back((*cl)->ID());
      CclusHitlists.push_back(hitlist);
      break;   
    }
    
    delete pol1;
  }// end of loop over all cluster created by  HoughLineFinder

  

  /////////////////////////////////////////////////////
  /////// 2D Track Matching and 3D Track Reconstruction
  /////////////////////////////////////////////////////

  unsigned int Associations = 0;
  //  std::vector<recob::Track *> track3DVector;  //holds 3D tracks to be saved
  for(unsigned int collectionIter=0; collectionIter < CclusterIDs.size();collectionIter++){  //loop over Collection view 2D tracks
    // Recover previously stored info
    double Ccrosst0 = Ccrosstfirsts[collectionIter];
    double Ccrosst1 = Ccrosstlasts[collectionIter];
    double Cw0 = Cwirefirsts[collectionIter];
    double Cw1 = Cwirelasts[collectionIter];
    double Ct0 = Ctimefirsts[collectionIter];
    double Ct1 = Ctimelasts[collectionIter];
    double Ct0_line = Ctimefirsts_line[collectionIter];
    double Ct1_line = Ctimelasts_line[collectionIter];
    edm::PtrVector<recob::Hit> hitsCtrk = CclusHitlists[collectionIter];

    for(unsigned int inductionIter=0;inductionIter<IclusterIDs.size();inductionIter++){   //loop over Induction view 2D tracks
      // Recover previously stored info
      double Icrosst0 = Icrosstfirsts[inductionIter];
      double Icrosst1 = Icrosstlasts[inductionIter];
      double Iw0 = Iwirefirsts[inductionIter];
      double Iw1 = Iwirelasts[inductionIter];
      double It0 = Itimefirsts[inductionIter];
      double It1 = Itimelasts[inductionIter];
      double It0_line = Itimefirsts_line[inductionIter];
      double It1_line = Itimelasts_line[inductionIter];
      edm::PtrVector<recob::Hit> hitsItrk = IclusHitlists[inductionIter];

      // match 2D tracks
      if((fabs(Ct0_line-It0_line)<ftmatch*timepitch) && (fabs(Ct1_line-It1_line)<ftmatch*timepitch)){ 
	//std::cout<<"-----> Track "<<collectionIter<< " Collection associated with track "<<inductionIter<< " Induction"<<std::endl;
	Associations++;

        // Reconstruct the 3D track
	TVector3 XYZ0;  // track origin or interaction vertex
	XYZ0.SetXYZ(Ct0_line,(Cw0-Iw0)/(2.*TMath::Sin(Angle)),(Cw0+Iw0)/(2.*TMath::Cos(Angle))-YC/2.*TMath::Tan(Angle));


	//compute track startpoint and endpoint in Local co-ordinate system 
	const TVector3 startpointVec(XYZ0.X(),XYZ0.Y(),XYZ0.Z());
 	const TVector3 endpointVec(Ct1_line,(Cw1-Iw1)/(2.*TMath::Sin(Angle)),(Cw1+Iw1)/(2.*TMath::Cos(Angle))-YC/2.*TMath::Tan(Angle));
	const TVector3 startpointVecLocal = m_tpcVolumeUtility->WorldToLocal(startpointVec);
	const TVector3 endpointVecLocal = m_tpcVolumeUtility->WorldToLocal(endpointVec);
	Double_t startpoint[3],endpoint[3];
	startpoint[0] = startpointVecLocal[0];
	startpoint[1] = startpointVecLocal[1];
	startpoint[2] = startpointVecLocal[2];
	endpoint[0]   = endpointVecLocal[0];
	endpoint[1]   = endpointVecLocal[1];
	endpoint[2]   = endpointVecLocal[2]; 

	//compute track (normalized) cosine directions in the World co-ordinate system
        double DirCos[3];
        DirCos[0] = endpointVec.X()-startpointVec.X();
        DirCos[1] = endpointVec.Y()-startpointVec.Y();
        DirCos[2] = endpointVec.Z()-startpointVec.Z();
	double Norma = TMath::Sqrt(TMath::Power(DirCos[0],2.)+TMath::Power(DirCos[1],2.)+TMath::Power(DirCos[2],2.));
	DirCos[0] = DirCos[0]/Norma; 
	DirCos[1] = DirCos[1]/Norma;
	DirCos[2] = DirCos[2]/Norma;

	//compute 3D track parameters (theta, phi, track pitch length) in the World co-ordinate system         
	double Theta; // Zenith angle with respect to y (vertical) axis in radians 
	double Phi;   // Azimuth angle with respect to x axis in radians
	Theta = TMath::ACos(DirCos[1]);
	Phi = TMath::ATan(DirCos[2]/DirCos[0]);
	Phi = Phi<0. ? Phi+TMath::Pi() : Phi ; // solve the ambiguities due to tangent periodicity	
	double cosgammaC; // angle between the track direction and the direction of the Collection wire pitch
	cosgammaC = TMath::Sin(Theta)*TMath::Sin(Phi)*TMath::Cos(Angle)+TMath::Cos(Theta)*TMath::Sin(Angle);
	double TrackPitchC  = wire_pitch/cosgammaC; // Collection Track Pitch length (i.e. the effective length of the track seen by the Collection wire) in cm  	
	double cosgammaI; // angle between the track direction and the direction of the Induction wire pitch
	cosgammaI = TMath::Sin(Theta)*TMath::Sin(Phi)*TMath::Cos(Angle)-TMath::Cos(Theta)*TMath::Sin(Angle);
	double TrackPitchI  = wire_pitch/cosgammaI; // Induction Track Pitch length (i.e. the effective length of the track seen by the Induction wire) in cm
	//double TrackLength_line = TMath::Sqrt(TMath::Power((startpoint[0]-endpoint[0]),2.)+TMath::Power((startpoint[1]-endpoint[1]),2.)+TMath::Power((startpoint[2]-endpoint[2]),2.));

	/* This is how this works. Maddalena creates a 3D trk and then packs
	 2 k3D clusters into it. The first one contains the first and last hits
	 on the new trk. The second one contains all the matching ones in 
	 between.I will now add 2 more clusters. They will be kU, kT clusters. 
	 EC, 27-Aug-2010.
	*/ 
	//create the the 3D track
	std::vector<recob::Hit*> hits3D;
	edm::PtrVector <recob::Hit> hits3Dpv;
	edm::Ptr<recob::Hit> hit1();
 	hit1->SetView(geo::k3D);
	hit1->SetXYZ(startpoint);
 	hits3Dpv.push_back(hit1);    
	edm::Ptr<recob::Hit> hit2();
 	hit2->SetView(geo::k3D);
	hit2->SetXYZ(endpoint);
  	hits3Dpv.push_back(hit2);

 	recob::Cluster cl3D(hits3Dpv);
	cl3D.SetID(Associations);
	recob::Track*  the3DTrack = new recob::Track();
	// Eric's next chunk here to include original hits. EC, 27-Aug-2010.
	// Playing along with the detector-specific coding for now.
	recob::Cluster* clOrigC = new recob::Cluster(hitsCtrk);
	clOrigC->SetID(Associations);
 	the3DTrack->Add(clOrigC);
 	recob::Cluster* clOrigI = new recob::Cluster(hitsItrk);
	clOrigI->SetID(Associations);
 	the3DTrack->Add(clOrigI);

	the3DTrack->SetDirection(DirCos,DirCos);
	the3DTrack->SetPitch(TrackPitchI,geo::kU);
	the3DTrack->SetPitch(TrackPitchC,geo::kV);
 	the3DTrack->Add(cl3D);



	/////////////////////////////
	// Match hits
	////////////////////////////

	std::vector<const recob::Hit*> minhits = hitsCtrk.size() <= hitsItrk.size() ? hitsCtrk : hitsItrk;
 	std::vector<const recob::Hit*> maxhits = hitsItrk.size() <= hitsCtrk.size() ? hitsCtrk : hitsItrk;



	bool maxhitsMatch[maxhits.size()];
	for(unsigned int it=0;it<maxhits.size();it++) maxhitsMatch[it] = false;

	std::vector<const recob::Hit*> hits3Dmatched;
	// For the matching start from the view where the track projection presents less hits
	unsigned int imaximum = 0;
	for(unsigned int imin=0;imin<minhits.size();imin++){ //loop over hits
	  //get wire - time coordinate of the hit
	  recob::Wire* theWire = (recob::Wire*)minhits[imin]->Wire();
	  int channel,wire,plane1,plane2;
	  channel = theWire->RawDigit()->Channel();
	  geom->ChannelToWire(channel,plane1,wire);
	  // get the wire-time co-ordinates of the hit to be matched
	  double w1 = (double)((wire+1)*wire_pitch);
 	  double t1 = plane1==1?(double)((minhits[imin]->CrossingTime()-presamplings-(tSI+tIC))*timepitch):(double)((minhits[imin]->CrossingTime()-presamplings-tSI)*timepitch); //in cm	  
	  double crossTime1 = minhits[imin]->CrossingTime();
	  double MIPs1 = minhits[imin]->MIPs();

	  //get the track origin co-ordinates in the two views
	  TVector2 minVtx2D;
	  plane1==1 ? minVtx2D.Set(Ct0,Cw0): minVtx2D.Set(It0,Iw0);
	  TVector2 maxVtx2D;
	  plane1==1 ? maxVtx2D.Set(It0,Iw0): maxVtx2D.Set(Ct0,Cw0);
	  
	  // get the track last hit co-ordinates in the two views
	  theWire = (recob::Wire*)minhits[minhits.size()-1]->Wire();
	  channel = theWire->RawDigit()->Channel();
	  geom->ChannelToWire(channel,plane1,wire);
	  double w1_last = plane1==1?Cw1:Iw1;
	  double t1_last = plane1==1?Ct1:It1;  
	  theWire = (recob::Wire*)maxhits[maxhits.size()-1]->Wire();
	  channel = theWire->RawDigit()->Channel();
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
	      theWire = (recob::Wire*)maxhits[imax]->Wire();
	      channel = theWire->RawDigit()->Channel();
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
	  
	  
	  // Get the time-wire co-ordinates of the matched hit
	  theWire = (recob::Wire*)maxhits[imaximum]->Wire();
	  channel = theWire->RawDigit()->Channel();
	  geom->ChannelToWire(channel,plane2,wire);
	  double w1_match = (double)((wire+1)*wire_pitch);
	  double t1_match = plane2==1?(double)((maxhits[imaximum]->CrossingTime()-presamplings-(tSI+tIC))*timepitch):(double)((maxhits[imaximum]->CrossingTime()-presamplings-tSI)*timepitch);
	  double crossTime1_match = maxhits[imaximum]->CrossingTime();
	  double MIPs1_match = maxhits[imaximum]->MIPs();

	  // create the 3D hit, compute its co-ordinates and add it to the 3D hits list	  
	  double Ct = plane1==1?t1:t1_match;
	  double Cw = plane1==1?w1:w1_match;
	  double Iw = plane1==1?w1_match:w1;
	  double crossTime = plane1==1?crossTime1:crossTime1_match;
	  double MIPs = plane1==1?MIPs1:MIPs1_match;

	  const TVector3 hit3d(Ct,(Cw-Iw)/(2.*TMath::Sin(Angle)),(Cw+Iw)/(2.*TMath::Cos(Angle))-YC/2.*TMath::Tan(Angle)); 
	  recob::Hit* hit = new recob::Hit();
	  hit->SetView(geo::k3D);
	  hit->SetCrossingTime(crossTime);
	  hit->SetMIPs(MIPs);
	  const TVector3 hit3dLocal = m_tpcVolumeUtility->WorldToLocal(hit3d);
	  Double_t hitcoord[3];
	  hitcoord[0] = hit3dLocal.X();
	  hitcoord[1] = hit3dLocal.Y();
	  hitcoord[2] = hit3dLocal.Z();
	  hit->SetXYZ(hitcoord);           
	  hits3Dmatched.push_back(hit);
	}

	/// Output file for Direction Cosines 
// 	if(hits3Dmatched.size()>15){
// 	  std::ofstream outfile;
// 	  outfile.open("dircos_bis.txt",std::ios::out | std::ios::app);
// 	  outfile<<"\n";
// 	  outfile<<"run          : "<<evt.Header().Run()<<"\n";
// 	  outfile<<"event        : "<<evt.Header().Event()<<"\n";
// 	  outfile<<"trk3D ID     : "<<Associations<<"\n";
// 	  outfile<<"Direction Cosine : "<< DirCos[0]<<" "<<DirCos[1]<<" "<<DirCos[2]<<"\n";
// 	  outfile<<"Theta (deg)      : "<< Theta/TMath::Pi()*180. <<"\n";
// 	  outfile<<"Phi   (deg)      : "<< Phi/TMath::Pi()*180. <<"\n";
// 	  outfile<<"Start Point (cm) : "<< startpointVec.X()<<" "<<startpointVec.Y()<<" "<<startpointVec.Z()<<"\n";
// 	  outfile<<"End Point   (cm) : "<< endpointVec.X()<<" "<<endpointVec.Y()<<" "<<endpointVec.Z()<<"\n";
// 	  outfile.close();
// 	}
	//
	/// Output file for Drift Velocity 
	if(hits3Dmatched.size()>15){
	  std::ofstream outf;
	  outf.open("vdrift.txt",std::ios::out | std::ios::app);
	  outf<<"\n";
	  outf<<"run          : "<<evt.Header().Run()<<"\n";
	  outf<<"event        : "<<evt.Header().Event()<<"\n";
	  outf<<"trk3D ID     : "<<Associations<<"\n";
	  outf<<"Direction Cosine : "<< DirCos[0]<<" "<<DirCos[1]<<" "<<DirCos[2]<<"\n";
	  outf<<"Theta (deg)      : "<< Theta/TMath::Pi()*180. <<"\n";
	  outf<<"Phi   (deg)      : "<< Phi/TMath::Pi()*180. <<"\n";
	  outf<<"Start Point (cm) : "<< startpointVec.X()<<" "<<startpointVec.Y()<<" "<<startpointVec.Z()<<"\n";
	  outf<<"SP CrossT_C (sam): "<< Ccrosst0<<"\n";
	  outf<<"SP CrossT_I (sam): "<< Icrosst0<<"\n";
	  outf<<"End Point   (cm) : "<< endpointVec.X()<<" "<<endpointVec.Y()<<" "<<endpointVec.Z()<<"\n";
	  outf<<"EP CrossT_C (sam): "<< Ccrosst1<<"\n";
	  outf<<"EP CrossT_I (sam): "<< Icrosst1<<"\n";
	  outf.close();
	}
	//

	// Add the 3D track to the vector of the reconstructed tracks
        if(hits3Dmatched.size()>0){
	  recob::Cluster* cl3Dmatched= new recob::Cluster(hits3Dmatched);
	  cl3Dmatched->SetID(Associations);
	  the3DTrack->Add(cl3Dmatched);
	  tcol->push_back(the3DTrack);
	}
      } //close match 2D tracks
      
    }//close loop over Induction view 2D tracks
    
  }//close loop over Collection xxview 2D tracks
  
  std::cout<<"Run "<<evt.Header().Run()<<" Event "<<evt.Header().Event()<<std::endl;
  std::cout<<"TRACK3DRECO found "<<Associations<<" 3D track(s)"<<std::endl;
  
  evt.put(tcol);
  

}
//------------------------------------------------------------------------------------//
double trkf::Track3Dreco::DriftVelocity(double Efield, double Temperature){
  
  // Dirft Velocity as a function of Electric Field and LAr Temperature
  // from : W. Walkowiak, NIM A 449 (2000) 288-294
  
  double vd;
  double field = Efield; //Electric field kV/cm   default = 0.5
  double T = Temperature;
  
  
  double P1,P2,P3,P4,P5,P6,T0;
  P1=-0.01481; // K^-1
  P2=-0.0075; // K^-1
  P3=0.141;//(kV/cm)^-1
  P4=12.4;//kV/cm
  P5=1.627;//(kV/cm)^-P6
  P6=0.317;
  T0 = 90.371; // K

  vd=(P1*(T-T0)+1)*(P3*field*TMath::Log(1+P4/field) + P5*TMath::Power(field,P6))+P2*(T-T0);

  vd /= 10.;

  return vd;// in cm/us
}
