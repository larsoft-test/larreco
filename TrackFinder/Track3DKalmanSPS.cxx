////////////////////////////////////////////////////////////////////////
//
// \file Track3DKalmanSPS.cxx
//
// \author echurch@fnal.gov
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
#include "art/Framework/Core/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Persistency/Common/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Core/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

// LArSoft includes
#include "TrackFinder/SpacePointService.h"
#include "TrackFinder/Track3DKalmanSPS.h"
#include "Geometry/geo.h"
#include "RecoBase/recobase.h"
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

//-------------------------------------------------
trkf::Track3DKalmanSPS::Track3DKalmanSPS(fhicl::ParameterSet const& pset) 
{

    this->reconfigure(pset);

    produces< std::vector<recob::Track> >();

    // set the random number seed
    fRandom = dynamic_cast<TRandom3*>(gRandom);
    if ( fRandom == 0 ){
      fRandom = new TRandom3(time(0));
      gRandom = fRandom;
    }
}

void trkf::Track3DKalmanSPS::reconfigure(fhicl::ParameterSet pset) 
  {

    fClusterModuleLabel   = pset.get< std::string >("ClusterModuleLabel");
    fGenieGenModuleLabel   = pset.get< std::string >("GenieGenModuleLabel");
    fPosErr                = pset.get< std::vector < double >  >("PosErr3");   // resolution. cm
    fMomErr                = pset.get< std::vector < double >  >("MomErr3");   // GeV
    fMomStart              = pset.get< std::vector < double >  >("MomStart3"); // Will be unit norm'd.
    fGenfPRINT             = pset.get< bool >("GenfPRINT");

  }

//-------------------------------------------------
trkf::Track3DKalmanSPS::~Track3DKalmanSPS()
{

  /*
    delete stMCT ;
    delete covMCT;
    delete stREC;
    delete covREC;
    
    delete fpMCT;
    delete fpREC;
    delete fpRECL;
    delete fpRECt3D;
  */
}

//-------------------------------------------------
void trkf::Track3DKalmanSPS::beginJob()
{


  art::ServiceHandle<art::TFileService> tfs;
  
  stMCT  = new TMatrixT<Double_t>(5,1);
  covMCT = new TMatrixT<Double_t>(5,5);
  stREC  = new TMatrixT<Double_t>(5,1);
  covREC = new TMatrixT<Double_t>(5,5);
  
  fpMCT = new Float_t[4];
  fpREC = new Float_t[4];
  fpRECL = new Float_t[4];
  fpRECt3D = new Float_t[4];
  fDimSize = 400; // if necessary will get this from pset in constructor.

  fshx = new Float_t[fDimSize];
  fshy = new Float_t[fDimSize];
  fshz = new Float_t[fDimSize];

//   //TFile fileGENFIT("GENFITout.root","RECREATE");

  
  tree = tfs->make<TTree>("GENFITttree","GENFITttree");
  //tree->Branch("stMCT",&stMCT,"stMCT[5]/F"); // "TMatrixT<Double_t>"

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
  tree->Branch("evtNo",&evtt,"evtNo/I");
  tree->Branch("chi2ndf",&chi2ndf,"chi2ndf/F");

  tree->Branch("trkNo",&nTrks,"trkNo/I");
  tree->Branch("ptsNo",&fptsNo,"ptsNo/I");
  tree->Branch("shx",fshx,"shx[ptsNo]/F");
  tree->Branch("shy",fshy,"shy[ptsNo]/F");
  tree->Branch("shz",fshz,"shz[ptsNo]/F");

  tree->Branch("pMCT",fpMCT,"pMCT[4]/F");
  tree->Branch("pRECKalF",fpREC,"pRECKalF[4]/F");
  tree->Branch("pRECKalL",fpRECL,"pRECKalL[4]/F");
  tree->Branch("pRECt3D",fpRECt3D,"pRECt3D[4]/F");
  

  //TGeoManager* geomGENFIT = new TGeoManager("Geometry", "Geane geometry");
  //TGeoManager::Import("config/genfitGeom.root");
  //  gROOT->Macro("config/Geane.C"); 
 
}

void trkf::Track3DKalmanSPS::endJob()
{
  if (!rep) delete rep;
  if (!repMC) delete repMC;
}


//------------------------------------------------------------------------------------//
void trkf::Track3DKalmanSPS::produce(art::Event& evt)
{ 

  rep=0;
  repMC=0;

  // get services
  art::ServiceHandle<geo::Geometry> geom;
  art::ServiceHandle<util::LArProperties> larprop;

  //////////////////////////////////////////////////////
  // Make a std::auto_ptr<> for the thing you want to put into the event
  // because that handles the memory management for you
  //////////////////////////////////////////////////////
  std::auto_ptr<std::vector<recob::Track> > tcol(new std::vector<recob::Track>);

  // define TPC parameters
  TString tpcName = geom->GetLArTPCVolumeName();


  // get input Hit object(s).
  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fClusterModuleLabel,clusterListHandle);
  art::ServiceHandle<trkf::SpacePointService> sps;

  art::PtrVector<simb::MCTruth> mclist;
  if (!evt.isRealData())
    {

      //      std::cout << "Track3DKalmanSPS: This is MC." << std::endl;
      // std::cout<<"Run "<<evt.run()<<" Event "<<evt.id().event()<<std::endl;

      art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
      evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle);

      for (unsigned int ii = 0; ii <  mctruthListHandle->size(); ++ii)
	{
	  art::Ptr<simb::MCTruth> mctparticle(mctruthListHandle,ii);
	  mclist.push_back(mctparticle);
	}
    }

  //create collection of spacepoints that will be used when creating the Track object
  std::vector<recob::SpacePoint> spacepoints;
  art::PtrVector<recob::Cluster> clusterIn;
  // std::cout<<"Run "<<evt.run()<<" Event "<<evt.id().event()<<std::endl;
  mf::LogInfo("Track3DKalmanSPS: ") << "There are " <<  clusterListHandle->size() << " Clusters in this event (over all the planes).";


  for(unsigned int ii = 0; ii < clusterListHandle->size(); ++ii)
    {
      art::Ptr<recob::Cluster> cluster(clusterListHandle, ii);
      clusterIn.push_back(cluster);
    }

      TVector3 MCOrigin;
      TVector3 MCMomentum;
      // TVector3 posErr(.05,.05,.05); // resolution. 0.5mm
      // TVector3 momErr(.1,.1,0.2);   // GeV
      TVector3 posErr(fPosErr[0],fPosErr[1],fPosErr[2]); // resolution. 0.5mm
      TVector3 momErr(fMomErr[0],fMomErr[1],fMomErr[2]);   // GeV

      // This is strictly for MC
      if (!evt.isRealData())
	{
	  // Below breaks are stupid, I realize. But rather than keep all the MC
	  // particles I just take the first primary, e.g., muon and only keep its
	  // info in the branches of the Ttree. I could generalize later, ...
	  for( unsigned int ii = 0; ii < mclist.size(); ++ii )
	    {
	    //art::Ptr<const simb::MCTruth> mc(mctruthListHandle,i);
	    art::Ptr<simb::MCTruth> mc(mclist[ii]);
	    for(int jj = 0; jj < mc->NParticles(); ++jj)
	      {
		simb::MCParticle part(mc->GetParticle(jj));
		mf::LogInfo("Track3DKalmanSPS: ") << "FROM MC TRUTH, the particle's pdg code is: "<<part.PdgCode()<< " with energy = "<<part.E() <<", with energy = "<<part.E()<< " and vtx and momentum in Global (not volTPC) coords are " ;
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
	  mf::LogInfo("Track3DKalmanSPS: ") <<" repMC, covMC are ... " ;
	  repMC->getState().Print();
	  repMC->getCov().Print();

	} // !isRealData
      
      //      art::PtrVectorItr<recob::Track> trackIter = trackIn.begin();
      
      // loop over tracks is obsolesced here in SPS. Instead, loop over clusters.
      // Take its hits, call makeSpacePoints(), out of which we'll get our spacepoints.

      for (art::PtrVectorItr<recob::Cluster> iclus = clusterIn.begin();
	   iclus != clusterIn.end(); ++iclus)
	{
	  for (art::PtrVectorItr<recob::Cluster> jclus = iclus+1;
	       jclus != clusterIn.end(); ++jclus)
	    {
	      for (art::PtrVectorItr<recob::Cluster> kclus = jclus+1;
		   kclus != clusterIn.end(); ++kclus)
		{
		  spacepoints.clear();
		  // Just concatenate all hits into one hits vector. makeSpacePoints()
		  // itself will enforce separate planes, etc., requirements. 

		  art::PtrVector<recob::Hit> hits;
		  art::PtrVector<recob::Hit> hitsi;
		  art::PtrVector<recob::Hit> hitsj;
		  art::PtrVector<recob::Hit> hitsk;

		  hitsi = (*iclus)->Hits();
		  hitsj = (*jclus)->Hits();
		  hitsk = (*kclus)->Hits();
		  
		  for (art::PtrVectorItr<recob::Hit> ihit = hitsi.begin(); ihit != hitsi.end(); ++ihit)
		    {
		      hits.push_back(*ihit);
		    }
		  for (art::PtrVectorItr<recob::Hit> ihit = hitsj.begin(); ihit != hitsj.end(); ++ihit)
		    {
		      hits.push_back(*ihit);
		    }
		  for (art::PtrVectorItr<recob::Hit> ihit = hitsk.begin(); ihit != hitsk.end(); ++ihit)
		    {
		      hits.push_back(*ihit);
		    }
		  
		  if (hits.size()==0) break;
		  // Add the 3D track to the vector of the reconstructed tracks
		  sps->makeSpacePoints(hits,spacepoints,&evt);
		  // Insert the GENFIT/Kalman stuff here then make the tracks. Units are cm, GeV.
		  mf::LogInfo("Track3DKalmanSPS: ") << "found "<< spacepoints.size() <<" 3D spacepoint(s) for Cluster combo ", iclus, jclus, kclus, " .";
		  
		  if(spacepoints.size()>0)
		    {

		      const double resolution = 0.5; // dunno, 5 mm
		      const int numIT = 9; // 3->1, EC, 6-Jan-2011. Back, 7-Jan-2011.
		      
	    //TVector3 mom(0.0,0.0,2.0);
		      TVector3 mom(fMomStart[0],fMomStart[1],fMomStart[2]);
	    //mom.SetMag(1.);
		      TVector3 momM(mom);
		      momM.SetX(gRandom->Gaus(momM.X(),momErr.X()/* *momM.X() */));
		      momM.SetY(gRandom->Gaus(momM.Y(),momErr.Y()/* *momM.Y() */));
		      momM.SetZ(gRandom->Gaus(momM.Z(),momErr.Z()/* *momM.Z() */));



		      std::sort(spacepoints.begin(), spacepoints.end(), sp_sort_3dz);
	   
	    
		      genf::GFFieldManager::getInstance()->init(new genf::GFConstField(0.,0.,0.0));
		      genf::GFDetPlane planeG((TVector3)(spacepoints[0].XYZ()),momM);
	    

      //      std::cout<<"Track3DKalmanSPS about to do GAbsTrackRep."<<std::endl;
      // Initialize with 1st spacepoint location and a guess at the momentum.
		      rep = new genf::RKTrackRep(//posM-.5/momM.Mag()*momM,
						 (TVector3)(spacepoints[0].XYZ()),
						 momM,
						 posErr,
						 momErr,
						 13);  // mu- hypothesis
      //      std::cout<<"Track3DKalmanSPS: about to do GFTrack. repDim is " << rep->getDim() <<std::endl;


		      genf::GFTrack fitTrack(rep);//initialized with smeared rep
      // Gonna sort in x cuz I want to essentially transform here to volTPC coords.
      // volTPC coords, cuz that's what the Geant3/Geane stepper wants, as that's its understanding
      // from the Geant4 geometry, which it'll use. EC, 7-Jan-2011.
		      int ihit = 0;
		      fptsNo = 0;
		      for (unsigned int point=0;point<spacepoints.size();++point)
			{
			  
			  TVector3 spt3 = (TVector3)(spacepoints[point].XYZ());
			  if (point%20) // Jump out of loop except on every 20th pt.
			    {
			      //continue;
			      // Icarus paper suggests we may wanna decimate our data in order to give
			      // trackfitter a better idea of multiple-scattering. EC, 7-Jan-2011.
			      //if (fabs(spt3[0]-spacepoints.at(point-1).XYZ()[0]) < 2) continue;
			    }
			  if (fptsNo<fDimSize)
			    {
			      fshx[fptsNo] = spt3[0];
			      fshy[fptsNo] = spt3[1];
			      fshz[fptsNo] = spt3[2];
			    }
			  fptsNo++;


			  mf::LogDebug("Track3DKalmanSPS: ") << "ihit xyz..." << spt3[0]<<","<< spt3[1]<<","<< spt3[2];
			  fitTrack.addHit(new genf::PointHit(spt3,resolution),
					  1,//dummy detector id
					  ihit++
					  );
			} // end loop over spacepoints.

	    //      std::cout<<"Track3DKalmanSPS about to do GFKalman."<<std::endl;
		      genf::GFKalman k;
      //k.setBlowUpFactor(500); // Instead of 500 out of box. EC, 6-Jan-2011.
      //k.setInitialDirection(+1); // Instead of 1 out of box. EC, 6-Jan-2011.
		      k.setNumIterations(numIT);
      //      std::cout<<"Track3DKalmanSPS back from setNumIterations."<<std::endl;
		      try{
	      //	std::cout<<"Track3DKalmanSPS about to processTrack."<<std::endl;
			k.processTrack(&fitTrack);
			//std::cout<<"Track3DKalmanSPS back from processTrack."<<std::endl;
		      }
		      catch(GFException& e){
			mf::LogError("Track3DKalmanSPS: ") << "just caught a GFException."<<std::endl;
			e.what();
			mf::LogError("Track3DKalmanSPS: ") << "Exceptions won't be further handled ->exit(1) "<<__LINE__;

			//	exit(1);
		      }

		      if(rep->getStatusFlag()==0) // 0 is successful completion
			{
			  mf::LogDebug("Track3DKalmanSPS: ") << __FILE__ << " " << __LINE__ ;
			  mf::LogDebug("Track3DKalmanSPS: ") << "Track3DKalmanSPS.cxx: Original plane:";
			  
			  if(fGenfPRINT) planeG.Print();
			  mf::LogDebug("Track3DKalmanSPS: ") << "Current (fit) reference Plane:";
			  if(fGenfPRINT) rep->getReferencePlane().Print();
		
			  mf::LogDebug("Track3DKalmanSPS: ") << "Track3DKalmanSPS.cxx: Last reference Plane:";
			  if(fGenfPRINT) rep->getLastPlane().Print();
		
			  if(fGenfPRINT) 
			    {
			      if(planeG!=rep->getReferencePlane()) 
				mf::LogDebug("Track3DKalmanSPS: ")	<<"Track3DKalmanSPS: Original hit plane (not surprisingly) not current reference Plane!"<<std::endl;
			    }

			  stREC->ResizeTo(rep->getState());
			  *stREC = rep->getState();
			  covREC->ResizeTo(rep->getCov());
			  *covREC = rep->getCov();
			  if(fGenfPRINT)
			    {
			      mf::LogDebug("Track3DKalmanSPS: ") << " Final State and Cov:";
			      stREC->Print();
			      covREC->Print();
			    }
			  chi2 = rep->getChiSqu();
			  ndf = rep->getNDF();
			  nfail = fitTrack.getFailedHits();
			  chi2ndf = chi2/ndf;
			  double dircoss[3],dircose[3];
			  //		(*trackIter)->Direction(dircoss,dircose);      
		
			  for (int ii=0;ii<3;++ii)
			    {
			      fpMCT[ii] = MCMomentum[ii]/MCMomentum.Mag();
			      fpREC[ii] = rep->getReferencePlane().getNormal()[ii];
			      fpRECL[ii] = rep->getLastPlane().getNormal()[ii];
			      fpRECt3D[ii] = dircoss[ii];
			    }
			  fpMCT[3] = MCMomentum.Mag();
			  fpREC[3] = -1.0/(*stREC)[0][0];
      
			  evtt = (unsigned int) evt.id().event();
			  mf::LogInfo("Track3DKalmanSPS: ") << "Track3DKalmanSPS about to do tree->Fill(). Chi2/ndf is " << chi2/ndf << ". All in volTPC coords .... pMCT[0-3] is " << fpMCT[0] << ", " << fpMCT[1] << ", " << fpMCT[2] << ", " << fpMCT[3] << ". pREC[0-3] is " << fpREC[0] << ", "<< fpREC[1] << ", " << fpREC[2] << ", " << fpREC[3] << ".";
	  
			  tree->Fill();


			  art::PtrVector<recob::Cluster> clusters;
			  clusters.push_back(*iclus);
			  clusters.push_back(*jclus);
			  clusters.push_back(*kclus);
		

			  recob::Track  the3DTrack(clusters,spacepoints);
			  double dircosF[3];
			  double dircosL[3];
			  for (int ii=0;ii<3;++ii)
			    {
			      dircosF[ii] = fpREC[ii];
			      dircosL[ii] = fpRECL[ii];
			    }
			  the3DTrack.SetDirection(dircosF,dircosL);
			  the3DTrack.SetID(tcol->size());
			  tcol->push_back(the3DTrack);
			} // getStatusFlag 
		    } // spacepoints.size()>0
	
		  
		  evt.put(tcol);

		}  // loop on kclus
	    } // jclus
	} // iclus
}
