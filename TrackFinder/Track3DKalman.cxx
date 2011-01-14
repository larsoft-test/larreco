////////////////////////////////////////////////////////////////////////
//
// Track3DKalman class
//
// echurch@fnal.gov
//
//  This algorithm is designed to reconstruct 3D tracks through  
//  GENFIT's Kalman filter.
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
#include "Track3DKalman.h"
#include "Geometry/geo.h"
#include "Geometry/WireGeo.h"
#include "Geometry/VolumeUtility.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Track.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/SpacePoint.h"
#include "Utilities/LArProperties.h"
#include "SimulationBase/simbase.h"
#include "Simulation/sim.h"
#include "Simulation/SimListUtils.h"


// ROOT includes
#include "TVectorD.h"
#include "TFile.h"
#include "TGeoManager.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMath.h"

// GENFIT includes
#include "Genfit/GFException.h"
#include "Genfit/GFAbsTrackRep.h"
#include "Genfit/GeaneTrackRep2.h"
#include "Genfit/RKTrackRep.h"
#include "Genfit/GFConstField.h"
#include "Genfit/GFFieldManager.h"
#include "Genfit/PointHit.h"
#include "Genfit/GFTrack.h"
#include "Genfit/GFKalman.h"
#include "Genfit/GFDaf.h"

static bool sp_sort_3dz(const recob::SpacePoint& h1, const recob::SpacePoint& h2)
{
  const double* xyz1 = h1.XYZ();
  const double* xyz2 = h2.XYZ();
  return xyz1[2] < xyz2[2];
}
static bool sp_sort_3dx(const recob::SpacePoint& h1, const recob::SpacePoint& h2)
{
  const double* xyz1 = h1.XYZ();
  const double* xyz2 = h2.XYZ();
  return xyz1[0] > xyz2[0];
}

//-------------------------------------------------
trkf::Track3DKalman::Track3DKalman(edm::ParameterSet const& pset) :
  fHitModuleLabel     (pset.getParameter< std::string >("HitModuleLabel")),
  fGenieGenModuleLabel     (pset.getParameter< std::string >("GenieGenModuleLabel")),
  ftmatch                 (pset.getParameter< int    >("TMatch")),
  fchi2dof                (pset.getParameter< double >("Chi2DOFmax")),
  fGenfPRINT                (pset.getParameter< bool >("GenfPRINT"))
{
  produces< std::vector<recob::Track> >();

    // set the random number seed
    fRandom = dynamic_cast<TRandom3*>(gRandom);
    if ( fRandom == 0 ){
      fRandom = new TRandom3(time(0));
      gRandom = fRandom;
    }
}

//-------------------------------------------------
trkf::Track3DKalman::~Track3DKalman()
{
}

//-------------------------------------------------
void trkf::Track3DKalman::beginJob(const edm::EventSetup&)
{
  // Purposely, don't do this with auto-mojo of the TFileService.
  // But, I change my mind and do wanna do it this way now. EC, 7-Jan-2011.
  edm::Service<edm::TFileService> tfs;
  stMCT  = new TMatrixT<Double_t>(5,1);
  covMCT = new TMatrixT<Double_t>(5,5);
  stREC  = new TMatrixT<Double_t>(5,1);
  covREC = new TMatrixT<Double_t>(5,5);
  pMCT = new Double_t[4];
  pREC = new Double_t[4];

//   //TFile fileGENFIT("GENFITout.root","RECREATE");
  tree = tfs->make<TTree>("GENFITttree","GENFITttree");
//   //tree->Branch("stMCT",&stMCT,"stMCT[5]/F"); // "TMatrixT<Double_t>"
  tree->Branch("stMCT","TMatrixD",&stMCT,64000,0);
  //tree->Branch("covMCT",&covMCT,"covMCT[25]/F");
  tree->Branch("covMCT","TMatrixD",&covMCT,64000,0);
  //tree->Branch("stREC",&stREC,"stREC[5]/F");
  tree->Branch("stREC","TMatrixD",&stREC,64000,0);
  //tree->Branch("covREC",&covREC,"covREC[25]/F");
  tree->Branch("covREC","TMatrixD",&covREC,64000,0);
  tree->Branch("chi2",&chi2,"chi2/F");
  tree->Branch("nfail",&nfail,"nfail/I");
  tree->Branch("ndf",&ndf,"ndf/I");
  tree->Branch("chi2ndf",&chi2ndf,"chi2ndf/F");
  tree->Branch("pMCT",&pMCT,"pMCT[4]/F");
  tree->Branch("pREC",&pREC,"pREC[4]/F");

  //TGeoManager* geomGENFIT = new TGeoManager("Geometry", "Geane geometry");
  //TGeoManager::Import("config/genfitGeom.root");
  //  gROOT->Macro("config/Geane.C"); 
  
// Let's presume my already-defined geometry and physics list will do what's needed.
  // get services
  edm::Service<geo::Geometry> geom;
  edm::Service<util::LArProperties> larprop;

}

void trkf::Track3DKalman::endJob()
{
  //tree->Write();
  //fileGENFIT->Close();
}


//------------------------------------------------------------------------------------//
void trkf::Track3DKalman::produce(edm::Event& evt, edm::EventSetup const&)
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
  double tSI = plane_pitch/driftvelocity_SI/timetick;    //drift time between Shield and Induction planes (time samples)
  double tIC = plane_pitch/driftvelocity_IC/timetick;    //drift time between Induction and Collection planes (time samples)


  // get input Hit object(s).
  edm::Handle< std::vector<recob::Hit> > hitListHandle;
  evt.getByLabel(fHitModuleLabel,hitListHandle);

  edm::PtrVector<simb::MCTruth> mclist;
  if (!evt.isRealData())
    {
      std::cout << "Track3DKalman: This is MC." << std::endl;

      edm::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
      evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle);

      for (unsigned int ii = 0; ii <  mctruthListHandle->size(); ++ii)
	{
	  edm::Ptr<simb::MCTruth> mctparticle(mctruthListHandle,ii);
	  mclist.push_back(mctparticle);
	}
    }

  //create collection of spacepoints that will be used when creating the Track object
  std::vector<recob::SpacePoint> spacepoints;
  for(unsigned int ii = 0; ii < hitListHandle->size(); ++ii)
    {
      edm::Ptr<recob::Hit> hit(hitListHandle, ii);
      // Figure out which View the hit belongs to 
      geo::View_t view = hit->View();

      // Some variables for the hit
      unsigned int channel;  //channel number
      float time;            //hit time at maximum
      unsigned int wire;              //hit wire number
      unsigned int plane;             //hit plane number
      
      std::vector<double> wires;
      std::vector<double> times;
           
      time = hit->CrossingTime() ;
      time -= presamplings;
      channel = hit->Channel();
      geom->ChannelToWire(channel,plane,wire);
	  
      //correct for the distance between wire planes
      if(plane==0) time -= tSI;         // Induction
      if(plane==1) time -= (tSI+tIC);   // Collection
	
      //transform hit wire and time into cm
      double wire_cm = (double)((wire+1) * wire_pitch); 
      double time_cm = (double)(time * timepitch);
      wires.push_back(wire_cm);
      times.push_back(time_cm);

      // Loop over next set of hits and only consider if we're in the other plane.
      std::vector<unsigned int> hitIndex;
      for(unsigned int jj = ii+1; jj < hitListHandle->size(); ++jj)
	{
	  edm::Ptr<recob::Hit> hit2(hitListHandle, jj);
	  unsigned int channel2;  //channel number
	  float time2;            //hit time at maximum
	  unsigned int wire2;              //hit wire number
	  unsigned int plane2;             //hit plane number
      
	  std::vector<double> wires2;
	  std::vector<double> times2;
           
	  time2 = hit2->CrossingTime() ;
	  time2 -= presamplings;
	  channel2 = hit2->Channel();
	  geom->ChannelToWire(channel2,plane2,wire2);
	  if (plane2 == plane) continue;

	  //correct for the distance between wire planes
	  if(plane2==0) time2 -= tSI;         // Induction
	  if(plane2==1) time2 -= (tSI+tIC);   // Collection
	
	  //transform hit wire and time into cm
	  double wire2_cm = (double)((wire2+1) * wire_pitch); 
	  double time2_cm = (double)(time2 * timepitch);
	  wires2.push_back(wire2_cm);
	  times2.push_back(time2_cm);

	  
	  //compute 3D track parameters (theta, phi, track pitch length) in the World co-ordinate system
	  /* Keep this code chunk for copying from later when I need pitch for the Track properties.
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
	  */

	      
	  // bool timematch = (fabs(t1-t2)<ftmatch*timepitch); // Italians.
	  // EC, # of ticks drift between the 2 planes.
	  bool timematch = fabs(time-time2)/timetick < fabs(ftmatch*tIC); 
	  bool wirematch = fabs(wire-wire2) < fabs(wireShift*wire_pitch);
	  
	  //std::cout<<"Track3DKalman . timematch lhs " << fabs(time-time2)/timetick<< "wirematch lhs " << fabs(wire-wire2) << "." <<std::endl;
	  //	  std::cout<<"Track3DKalman . timematch rhs " << ftmatch*tIC << "wirematch rhs " << wireShift*wire_pitch<< "." <<std::endl;
	  if (wirematch && timematch)
	    {
	      //std::cout<<"Track3DKalman found a pair of matching 2D hit(s) in opposite planes. wirematch " << wirematch << " timematch " << timematch << "." <<std::endl;
	      hitIndex.push_back(jj);
	    }
	} // end loop on jj hits
      

	  // One loop just to fill sp_hits, cuz I want the full set of matched hits to
	  // go with pair of hits that form a SpacePoint. EC, 21-Dec-2010.
      edm::PtrVector<recob::Hit> sp_hits;
      for (unsigned int kk=0; kk<hitIndex.size() ; ++kk)
	{
	  edm::Ptr<recob::Hit> hit3(hitListHandle, hitIndex[kk]);	  
	  sp_hits.push_back(hit3);
	}
      // Now pair up hit with hit3 to make a SpacePoint.
      for (unsigned int kk=0; kk<hitIndex.size() ; ++kk)
	{
	  edm::Ptr<recob::Hit> hit3(hitListHandle, hitIndex[kk]);	  
	  unsigned int channel2;  //channel number
	  float time2;            //hit time at maximum
	  unsigned int wire2;              //hit wire number
	  unsigned int plane2;             //hit plane number
	  
	  std::vector<double> wires2;
	  std::vector<double> times2;
          
	  time2 = hit3->CrossingTime() ;
	  time2 -= presamplings;
	  channel2 = hit3->Channel();
	  geom->ChannelToWire(channel2,plane2,wire2);
	  
	  //correct for the distance between wire planes
	  if(plane2==0) time2 -= tSI;         // Induction
	  if(plane2==1) time2 -= (tSI+tIC);   // Collection
	  
	  //transform hit wire and time into cm
	  double wire2_cm = (double)((wire2+1) * wire_pitch); 
	  double time2_cm = (double)(time2 * timepitch);
	  
	  // create the 3D hit, compute its co-ordinates and add it to the 3D hits list. 1 is Colln.
	  double Ct = plane==1?time_cm:time2_cm;
	  double Cw = plane==1?wire_cm:wire2_cm;
	  double Iw = plane==1?wire2_cm:wire_cm;

	  const TVector3 hit3d(Ct,(Cw-Iw)/(2.*TMath::Sin(Angle)),(Cw+Iw)/(2.*TMath::Cos(Angle))-YC/2.*TMath::Tan(Angle)); 
	  const TVector3 hit3dLocal = geom->Plane(plane).WorldToLocal(hit3d);
	  Double_t hitcoord[3];
	  // The -25 gets us to volTPC coords, I claim.
	  hitcoord[2] = hit3d.X()-25.;  //  Temporarily, EC, 3-Jan-2011: hit3dLocal->hit3d
	  hitcoord[1] = hit3d.Y();  //   Local
	  hitcoord[0] = hit3d.Z()+20;  //   Local      
	  
	  recob::SpacePoint mysp(sp_hits);
	  mysp.SetXYZ(hitcoord);
	  spacepoints.push_back(mysp);
	  
	}
	  	  
    
    }//close loop over 1st plane hits.
  


  std::cout<<"Track3DKalman found "<< spacepoints.size() <<" 3D spacepoint(s)."<<std::endl;
// Add the 3D track to the vector of the reconstructed tracks
  if(spacepoints.size()>0)
    {
      // Insert the GENFIT/Kalman stuff here then make the tracks. Units are cm, GeV.
      const double posSmearXY = 0.1;
      const double posSmearZ = .2;
      const double momSmear = .2;
      const double resolution = 0.5; // dunno, 5 mm
      const int numIT = 3; // 3->1, EC, 6-Jan-2011. Back, 7-Jan-2011.

      TVector3 pos(39.0,0.0,0.0); // Doesn't get used. EC, 7-Jan-2011.
      TVector3 mom(-4,0.0,0.);
      TVector3 posErr(.5,.5,.5); // resolution. 5mm, which is too large.
      TVector3 momErr(3,.5,.5);

      mom.SetMag(1.);

      TVector3 posM(pos);
      posM.SetX(gRandom->Gaus(posM.X(),posSmearXY*posM.X()));
      posM.SetY(gRandom->Gaus(posM.Y(),posSmearXY*posM.Y()));
      posM.SetZ(gRandom->Gaus(posM.Z(),posSmearZ*posM.Z()));
      TVector3 momM(mom);
      momM.SetX(gRandom->Gaus(momM.X(),momSmear*momM.X()));
      momM.SetY(gRandom->Gaus(momM.Y(),momSmear*momM.Y()));
      momM.SetZ(gRandom->Gaus(momM.Z(),momSmear*momM.Z()));
      std::cout << "Track3DKalman: sort spacepoints by x (volTPC coords) for tracker's sake." << std::endl;
      //      std::sort(spacepoints.begin(), spacepoints.end(), sp_sort_3dz);
      std::sort(spacepoints.begin(), spacepoints.end(), sp_sort_3dx); // Reverse sort!

      genf::GFFieldManager::getInstance()->init(new genf::GFConstField(0.,0.,0.0));
      genf::GFDetPlane planeG((TVector3)(spacepoints.at(0).XYZ()),momM);
      genf::GFAbsTrackRep *repMC;
      TVector3 MCOrigin;
      TVector3 MCMomentum;

      // This is strictly for MC
      if (!evt.isRealData())
	{
	  // Below breaks are stupid, I realize.
	  for( unsigned int ii = 0; ii < mclist.size(); ++ii )
	    {
	    //edm::Ptr<const simb::MCTruth> mc(mctruthListHandle,i);
	    edm::Ptr<simb::MCTruth> mc(mclist[ii]);
	    for(int jj = 0; jj < mc->NParticles(); ++jj)
	      {
		simb::MCParticle part(mc->GetParticle(jj));
		std::cout<<"FROM MC TRUTH, the particle's pdg code is: "<<part.PdgCode();
		std::cout<<" with energy = "<<part.E()<<std::endl;
		std::cout<<" and vtx and momentum in Global (not volTPC) coords are "<<std::endl;
		MCOrigin.SetXYZ(part.Vx(),part.Vy(),part.Vz()); // V for Vertex
		MCMomentum.SetXYZ(part.Px(),part.Py(),part.Pz());
		MCOrigin.Print();
		MCMomentum.Print();
		repMC = new genf::RKTrackRep(MCOrigin,
					     MCMomentum,
					     posErr,
					     momErr,
					     part.PdgCode());  
		break;
	      }
	    break;
	  }
	  //for saving of MC truth    
	  stMCT->ResizeTo(repMC->getState());
	  *stMCT = repMC->getState();
	  covMCT-> ResizeTo(repMC->getCov());
	  *covMCT = repMC->getCov();
	  std::cout<<"Track3DKalman: repMC, covMC are ... " << std::endl;;
	  repMC->getState().Print();
	  repMC->getCov().Print();

	} // !isRealData
      

      std::cout<<"Track3DKalman about to do GAbsTrackRep."<<std::endl;
      // Initialize with 1st spacepoint location and a guess at the momentum.
      genf::GFAbsTrackRep *rep = new genf::RKTrackRep(//posM-.5/momM.Mag()*momM,
						      (TVector3)(spacepoints.at(0).XYZ()),
						      momM,
						      posErr,
						      momErr,
						      13);  // mu- hypothesis
      std::cout<<"Track3DKalman: about to do GFTrack. repDim is " << rep->getDim() <<std::endl;
      std::cout<<"Track3DKalman: number of spacepoints is " << spacepoints.size() <<std::endl;
      genf::GFTrack fitTrack(rep);//initialized with smeared rep
      // Gonna sort in x cuz I want to essentially transform here to volTPC coords.
      // volTPC coords, cuz that's what the Geant3/Geane stepper wants, as that's its understanding
      // from the Geant4 geometry, which it'll use. EC, 7-Jan-2011.
      int ihit = 0;
      for (unsigned int point=0;point<spacepoints.size();++point)
	{
	  TVector3 spt3 = (TVector3)(spacepoints.at(point).XYZ());
	  if (point>0)
	    {
	      // Icarus paper suggests we may wanna decimate our data in order to give
	      // trackfitter a better idea of multiple-scattering. EC, 7-Jan-2011.
	      //if (fabs(spt3[0]-spacepoints.at(point-1).XYZ()[0]) < 2) continue;
	    }
	  std::cout<<"Track3DKalman ihit xyz..." << spt3[0]<<","<< spt3[1]<<","<< spt3[2] <<std::endl;
	  fitTrack.addHit(new genf::PointHit(spt3,resolution),
			  1,//dummy detector id
			  ihit++
			  );
	}

      std::cout<<"Track3DKalman about to do GFKalman."<<std::endl;
      genf::GFKalman k;
      //k.setBlowUpFactor(500); // Instead of 500 out of box. EC, 6-Jan-2011.
      //k.setInitialDirection(+1); // Instead of 1 out of box. EC, 6-Jan-2011.
      k.setNumIterations(numIT);
      std::cout<<"Track3DKalman back from setNumIterations."<<std::endl;
      try{
	//	std::cout<<"Track3DKalman about to processTrack."<<std::endl;
	k.processTrack(&fitTrack);
	//std::cout<<"Track3DKalman back from processTrack."<<std::endl;
      }
      catch(GFException& e){
	std::cout<<"Track3DKalman just caught a GFException."<<std::endl;
	e.what();
	std::cerr<<"Exceptions won't be further handled ->exit(1) "<<__LINE__<<std::endl;
	 //	exit(1);
      }

      if(rep->getStatusFlag()==0) // 0 is successful completion
	{
	  if(fGenfPRINT)std::cout << __FILE__ << " " << __LINE__ << std::endl;
	  if(fGenfPRINT)std::cout << "Track3DKalman.cxx: Original plane:" << std::endl;
	  planeG.Print();
	  if(fGenfPRINT)std::cout << "Track3DKalman.cxx: Current (fit) reference Plane:" << std::endl;
	  rep->getReferencePlane().Print();
	  if(planeG!=rep->getReferencePlane()) 
	    if(fGenfPRINT)std::cout<<"Track3DKalman: Original hit plane (not surprisingly) not current reference Plane!"<<std::endl;
	}

      stREC->ResizeTo(rep->getState());
      *stREC = rep->getState();
      covREC->ResizeTo(rep->getCov());
      *covREC = rep->getCov();
      if(fGenfPRINT)
	{
	  std::cout << "Track3DKalman.cxx: Final State and Cov:" << std::endl;
	  stREC->Print();
	  covREC->Print();
	}
      chi2 = rep->getChiSqu();
      ndf = rep->getNDF();
      nfail = fitTrack.getFailedHits();
      chi2ndf = chi2/ndf;
      
      for (int ii=0;ii<3;++ii)
	{
	  pMCT[ii] = MCMomentum[ii]/MCMomentum.Mag();
	  pREC[ii] = rep->getReferencePlane().getNormal()[ii];
	}
      pMCT[3] = MCMomentum.Mag();
      pREC[3] = -1.0/(*stREC)[0][0];
      
      std::cout<<"Track3DKalman about to do tree->Fill(). Chi2/ndf is " << chi2/ndf << ". All in volTPC coords .... pMCT[0-3] is " << -pMCT[2] << ", " << pMCT[1] << ", " << -pMCT[0] << ", " << pMCT[3] << ". pREC[0-3] is " << pREC[0] << ", "<< pREC[1] << ", " << pREC[2] << ", " << pREC[3] << "." <<std::endl;
      tree->Fill();


      // May need a new Track constructor... EC, 23-Dec-2010.
      /*recob::Track  the3DTrack(clustersPerTrack,spacepoints);
      double dircos[3];
      DirCos.GetXYZ(dircos);
      the3DTrack.SetDirection(dircos,dircos);
      the3DTrack.SetTrackPitch(TrackPitchI,geo::kU);
      the3DTrack.SetTrackPitch(TrackPitchC,geo::kV);
      
      tcol->push_back(the3DTrack);
      */
    } // spacepoints.size()>0


  std::cout<<"Run "<<evt.run()<<" Event "<<evt.id().event()<<std::endl;
  std::cout<<"Track3DKalman found "<< tcol->size() <<" 3D track(s)"<<std::endl;
  
  //  evt.put(tcol);
  

}
