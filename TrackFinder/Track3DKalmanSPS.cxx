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
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

// LArSoft includes
#include "TrackFinder/SpacePointService.h"
#include "TrackFinder/Track3DKalmanSPS.h"
#include "Geometry/geo.h"
#include "RecoBase/recobase.h"
#include "Utilities/LArProperties.h"
#include "SimulationBase/simbase.h"
#include "Simulation/sim.h"
#include "Utilities/AssociationUtil.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"


// ROOT includes
#include "TVectorD.h"
#include "TFile.h"
#include "TGeoManager.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMath.h"
#include "TPrincipal.h"

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
#include "Utilities/AssociationUtil.h"

static bool sp_sort_3dz(const recob::SpacePoint& h1, const recob::SpacePoint& h2)
{
  const double* xyz1 = h1.XYZ();
  const double* xyz2 = h2.XYZ();
  return xyz1[2] < xyz2[2];
}

//-------------------------------------------------
trkf::Track3DKalmanSPS::Track3DKalmanSPS(fhicl::ParameterSet const& pset) :
  fDoFit(true), fDecimate(1)
{

    this->reconfigure(pset);

    produces< std::vector<recob::Track> >();
    produces<art::Assns<recob::Track, recob::Cluster> >();
    // get the random number seed, use a random default if not specified    
    // in the configuration file.  
    unsigned int seed = pset.get< unsigned int >("Seed", sim::GetRandomNumberSeed());

    createEngine( seed );

}

void trkf::Track3DKalmanSPS::reconfigure(fhicl::ParameterSet const& pset) 
{
  
  fClusterModuleLabel    = pset.get< std::string >("ClusterModuleLabel");
  fProngModuleLabel      = pset.get< std::string >("ProngModuleLabel");
  fGenieGenModuleLabel   = pset.get< std::string >("GenieGenModuleLabel");
  fG4ModuleLabel         = pset.get< std::string >("G4ModuleLabel");
  fPosErr                = pset.get< std::vector < double >  >("PosErr3");   // resolution. cm
  fMomErr                = pset.get< std::vector < double >  >("MomErr3");   // GeV
  fMomStart              = pset.get< std::vector < double >  >("MomStart3"); // 
  fPerpLim               = pset.get< double  >("PerpLimit"); // 
  fDoFit                 = pset.get< bool  >("DoFit", true); // 
  fDecimate              = pset.get< int  >("Decimate", 1); // 
  fGenfPRINT             = pset.get< bool >("GenfPRINT");
  
}

//-------------------------------------------------

trkf::Track3DKalmanSPS::~Track3DKalmanSPS()
{
}


//-------------------------------------------------
void trkf::Track3DKalmanSPS::beginJob()
{


  art::ServiceHandle<art::TFileService> tfs;
  
  stMCT  = new TMatrixT<Double_t>(5,1);
  covMCT = new TMatrixT<Double_t>(5,5);
  stREC  = new TMatrixT<Double_t>(5,1);
  covREC = new TMatrixT<Double_t>(5,5);
  
  fpMCMom = new Float_t[4];
  fpMCPos = new Float_t[4];
  fpREC = new Float_t[4];
  fpRECL = new Float_t[4];
  fpRECLE = new Float_t[4];
  fpRECt3D = new Float_t[4];
  fState0 = new Float_t[5];
  fCov0 = new Float_t[25];
  fDimSize = 20000; // if necessary will get this from pset in constructor.

  fshx = new Float_t[fDimSize];
  fshy = new Float_t[fDimSize];
  fshz = new Float_t[fDimSize];
  feshx = new Float_t[fDimSize];
  feshy = new Float_t[fDimSize];
  feshz = new Float_t[fDimSize];
  feshyz = new Float_t[fDimSize];
  fupdate = new Float_t[fDimSize];
  fth  = new Float_t[fDimSize];
  feth = new Float_t[fDimSize];
  fedudw = new Float_t[fDimSize];
  fedvdw = new Float_t[fDimSize];
  feu = new Float_t[fDimSize];
  fev = new Float_t[fDimSize];
  fsep = new Float_t[fDimSize];
  
  fPC1 = new Float_t[3];
  fPC2 = new Float_t[3];
  fPC3 = new Float_t[3];
  fPCmeans = new Float_t[3];
  fPCsigmas = new Float_t[3];
  fPCevals = new Float_t[3];

//   //TFile fileGENFIT("GENFITout.root","RECREATE");

  
  tree = tfs->make<TTree>("GENFITttree","GENFITttree");
  //tree->Branch("stMCT",&stMCT,"stMCT[5]/F"); // "TMatrixT<Double_t>"

  //tree->Branch("stMCT","TMatrixD",&stMCT,64000,0);
  //tree->Branch("covMCT",covMCT,"covMCT[25]/F");
  tree->Branch("covMCT","TMatrixD",&covMCT,64000,0);
  tree->Branch("stREC",fState0,"stREC[5]/F");
  //tree->Branch("stREC","TMatrixD",&stREC,64000,0);
  tree->Branch("covREC",fCov0,"covREC[25]/F");
  //tree->Branch("covREC","TMatrixD",&covREC,64000,0);
  
  
  tree->Branch("chi2",&chi2,"chi2/F");
  tree->Branch("nfail",&nfail,"nfail/I");
  tree->Branch("ndf",&ndf,"ndf/I");
  tree->Branch("evtNo",&evtt,"evtNo/I");
  tree->Branch("chi2ndf",&chi2ndf,"chi2ndf/F");

  tree->Branch("trkNo",&nTrks,"trkNo/I");
  tree->Branch("ptsNo",&fptsNo,"ptsNo/I");
  tree->Branch("cont",&fcont,"cont/I"); //O? Yes, O. Not 0, not L, ...
  tree->Branch("shx",fshx,"shx[ptsNo]/F");
  tree->Branch("shy",fshy,"shy[ptsNo]/F");
  tree->Branch("shz",fshz,"shz[ptsNo]/F");
  tree->Branch("sep",fsep,"sep[ptsNo]/F");
  tree->Branch("eshx",feshx,"eshx[ptsNo]/F");
  tree->Branch("eshy",feshy,"eshy[ptsNo]/F");
  tree->Branch("eshz",feshz,"eshz[ptsNo]/F");
  tree->Branch("eshyz",feshyz,"eshyz[ptsNo]/F");  
  tree->Branch("update",fupdate,"update[ptsNo]/F");
  tree->Branch("th",fth,"th[ptsNo]/F");  
  tree->Branch("eth",feth,"eth[ptsNo]/F");
  tree->Branch("edudw",fedudw,"edudw[ptsNo]/F");
  tree->Branch("edvdw",fedvdw,"edvdw[ptsNo]/F");
  tree->Branch("eu",feu,"eu[ptsNo]/F");
  tree->Branch("ev",fev,"ev[ptsNo]/F");


  tree->Branch("pcMeans", fPCmeans,"pcMeans[3]/F");
  tree->Branch("pcSigmas",fPCsigmas,"pcSigmas[3]/F");
  tree->Branch("pcEvals", fPCevals,"pcEvals[3]/F");
  tree->Branch("pcEvec1",fPC1,"pcEvec1[3]/F");
  tree->Branch("pcEvec2",fPC2,"pcEvec2[3]/F");
  tree->Branch("pcEvec3",fPC3,"pcEvec3[3]/F");


  tree->Branch("pMCMom",fpMCMom,"pMCMom[4]/F");
  tree->Branch("pMCPos",fpMCPos,"pMCPos[4]/F");
  tree->Branch("pRECKalF",fpREC,"pRECKalF[4]/F");
  tree->Branch("pRECKalL",fpRECL,"pRECKalL[4]/F");
  tree->Branch("pRECKalLE",fpRECLE,"pRECKalLE[4]/F");
  tree->Branch("pRECt3D",fpRECt3D,"pRECt3D[4]/F");
  

  //TGeoManager* geomGENFIT = new TGeoManager("Geometry", "Geane geometry");
  //TGeoManager::Import("config/genfitGeom.root");
  //  gROOT->Macro("config/Geane.C"); 
 
}

void trkf::Track3DKalmanSPS::endJob()
{
  if (!rep) delete rep;
  if (!repMC) delete repMC;

  /*
  //  not sure why I can't do these, but at least some cause seg faults.
  delete[] stMCT;
  delete[] covMCT;
  delete[] stREC;
  delete[] covREC;
  */

  delete[] fpREC;
  delete[] fpRECL;
  delete[] fpRECLE;
  delete[] fpRECt3D;
  delete[] fState0;
  delete[] fCov0;

  delete[] fshx;
  delete[] fshy;
  delete[] fshz;
  delete[] feshx;
  delete[] feshy;
  delete[] feshyz;
  delete[] feshz;
  delete[] fupdate;
  delete[] fth;
  delete[] feth;
  delete[] fedudw;
  delete[] fedvdw;
  delete[] feu;
  delete[] fev;
  delete[] fsep;

  delete[] fPCmeans;
  delete[] fPCsigmas;
  delete[] fPCevals;
  delete[] fPC1;
  delete[] fPC2;
  delete[] fPC3;

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
  std::auto_ptr< art::Assns<recob::Track, recob::Cluster> > assn(new art::Assns<recob::Track, recob::Cluster>); 
  unsigned int tcnt = 0;

  // define TPC parameters
  TString tpcName = geom->GetLArTPCVolumeName();


  // get input Hit object(s).
  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fClusterModuleLabel,clusterListHandle);
  art::Handle< std::vector<recob::Prong> > prongListHandle;
  evt.getByLabel(fProngModuleLabel,prongListHandle);


  art::PtrVector<simb::MCTruth> mclist;
  std::vector<const sim::SimChannel*> simChannelHandle;
  
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

      evt.getView(fG4ModuleLabel, simChannelHandle);      
    }

  //create collection of spacepoints that will be used when creating the Track object

  art::PtrVector<recob::Cluster> clusterIn;
  art::PtrVector<recob::Prong> prongIn;
  // std::cout<<"Run "<<evt.run()<<" Event "<<evt.id().event()<<std::endl;
  mf::LogInfo("Track3DKalmanSPS: ") << "There are " <<  prongListHandle->size() << " Prongs (spacepoint clumps) in this event.";

  art::ServiceHandle<art::RandomNumberGenerator> rng;
  CLHEP::HepRandomEngine &engine = rng->getEngine();
  CLHEP::RandGaussQ gauss(engine);
  
  for(unsigned int ii = 0; ii < clusterListHandle->size(); ++ii)
    {
      art::Ptr<recob::Cluster> cluster(clusterListHandle, ii);
      clusterIn.push_back(cluster);
    }
  for(unsigned int ii = 0; ii < prongListHandle->size(); ++ii)
    {
      art::Ptr<recob::Prong> prong(prongListHandle, ii);
      prongIn.push_back(prong);
    }

      TVector3 MCOrigin;
      TVector3 MCMomentum;
      // TVector3 posErr(.05,.05,.05); // resolution. 0.5mm
      // TVector3 momErr(.1,.1,0.2);   // GeV
      TVector3 posErr(fPosErr[0],fPosErr[1],fPosErr[2]); // resolution. 0.5mm
      TVector3 momErr(fMomErr[0],fMomErr[1],fMomErr[2]);   // GeV
      TVector3 momErrFit(fMomErr[0],fMomErr[1],fMomErr[2]);   // GeV

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
      nTrks = 0;

      // loop on Prongs
      size_t cntp(0);
      art::PtrVector<recob::Prong>::const_iterator pprong = prongIn.begin();
      while (pprong!=prongIn.end()) 
	{
	  const std::vector<recob::SpacePoint>& spacepoints = (*pprong)->SpacePoints();

	  if (spacepoints.size()<10) { pprong++; continue;} // for now...
		  
		  mf::LogInfo("Track3DKalmanSPS: ")<<"\n\t found "<<spacepoints.size()<<" 3D spacepoint(s) for this prong \n";
				  
		  const double resolution = posErr.Mag(); 
		  const int numIT = 5; 
		  // Let's find track's principle components.
		  // We will sort along that direction, rather than z.
		  // Further, we will skip outliers away from main axis.
		  TVectorD* evals; 
		  TMatrixD* evecs; 
		  TVectorD* means; 
		  TVectorD* sigmas; 

		  TPrincipal* principal = new TPrincipal(3,"ND");
		  Double_t* data = new Double_t[3];
		  evals = new TVectorD(3);
		  evecs = new TMatrixD(3,3);
		  means = new TVectorD(3);
		  sigmas= new TVectorD(3);
		  // I need to shuffle these around, so use copy constructor
		  // to make non-const version spacepointss.
		  std::vector<recob::SpacePoint> spacepointss(spacepoints);
		  std::sort(spacepointss.begin(), spacepointss.end(), sp_sort_3dz);
		  int nTailPoints = 0; // 100;
		  for (unsigned int point=0;point<(spacepointss.size()-nTailPoints);++point)
		    {
		      data[0] = spacepointss[point].XYZ()[0];
		      data[1] = spacepointss[point].XYZ()[1];
		      data[2] = spacepointss[point].XYZ()[2];
		      //		      std::cout << "Spacepoint " << point << " added:" << spacepointss[point].XYZ()[0]<< ", " << spacepointss[point].XYZ()[1]<< ", " << spacepointss[point].XYZ()[2]<< ". " << std::endl;
		      principal->AddRow(data);
		    }
		  delete [] data;
		  principal->MakePrincipals();
		  /*
		  principal->Test();
		  principal->MakeHistograms();
		  principal->Print("MSEV");
		  */
		  evals = (TVectorD*)(principal->GetEigenValues());
		  evecs = (TMatrixD*)(principal->GetEigenVectors());
		  means = (TVectorD*)(principal->GetMeanValues());
		  sigmas= (TVectorD*)(principal->GetSigmas());

		  std::vector<TVector3*> pcs;
		  Double_t* pc = new Double_t[3];
		  for (unsigned int point=0;point<spacepointss.size();++point)
		    {
		      principal->X2P((Double_t *)(spacepointss[point].XYZ()),pc);
		      pcs.push_back((TVector3 *)pc);
		    }
		  delete [] pc;
		  Double_t tmp[3], tmp2[3];
		  principal->X2P((Double_t *)(means->GetMatrixArray()),tmp);
		  principal->X2P((Double_t *)(sigmas->GetMatrixArray()),tmp2);
		  for (unsigned int ii=0;ii<3;++ii)
		    {
		      fPCmeans[ii] = (Float_t )(tmp[ii]);
		      fPCsigmas[ii] = (Float_t )(tmp2[ii]);
		      fPCevals[ii] = (Float_t )(evals->GetMatrixArray())[ii];
		      // This method requires apparently pulling all 9
		      // elements. Maybe 3 works. 
		      // Certainly, w can't be a scalar, I discovered.
		      double w[9];
		      evecs->ExtractRow(ii,0,w);
		      fPC1[ii] = w[0];
		      fPC2[ii] = w[1];
		      fPC3[ii] = w[2];
		    }
		  Double_t tmp3[3], tmp4[3], tmp5[3];
		  principal->X2P((Double_t *)fPC1,tmp3);
		  principal->X2P((Double_t *)fPC2,tmp4);
		  principal->X2P((Double_t *)fPC3,tmp5);


		  // 21-Sep. Use a mip approximation assuming muons, assuming straight lines
		  // and a small angle wrt beam. 2.2 MeV/cm. And set momErrFit to momM/5.
		  fMomStart[0] = spacepointss[spacepointss.size()-1].XYZ()[0] - spacepointss[0].XYZ()[0];
		  fMomStart[1] = spacepointss[spacepointss.size()-1].XYZ()[1] - spacepointss[0].XYZ()[1];
		  fMomStart[2] = spacepointss[spacepointss.size()-1].XYZ()[2] - spacepointss[0].XYZ()[2];
		  // This presumes a muon. 
		  TVector3 mom(0.0022*fMomStart[0],0.0022*fMomStart[1],0.0022*fMomStart[2]);
		  // Over-estimate by just enough (10%).
		  mom.SetMag(1.1 * mom.Mag()); 
		  // My true 0.5 GeV/c muons need a yet bigger over-estimate.
		  //if (mom.Mag()<0.7) mom.SetMag(1.2*mom.Mag());  
		  //		  if (mom.Mag()>2.0) mom.SetMag(10.0*mom.Mag());  
		  //		  mom.SetMag(3*mom.Mag()); // EC, 15-Feb-2012. TEMPORARY!!!
		  // If any point is outside TPC, this track is uncontained.
		  // Try a higher momentum starting value in that case.
		  bool uncontained(false);
		  double close(10.); // cm
		  if (
		      spacepointss[spacepointss.size()-1].XYZ()[0] > (2.*geom->DetHalfWidth(0,0)-close) || spacepointss[spacepointss.size()-1].XYZ()[0] < close ||
		      spacepointss[0].XYZ()[0] > (2.*geom->DetHalfWidth(0,0)-close) || spacepointss[0].XYZ()[0] < close ||
		      spacepointss[spacepointss.size()-1].XYZ()[1] > (2.*geom->DetHalfHeight(0,0)-close) || (spacepointss[spacepointss.size()-1].XYZ()[1] < -2.*geom->DetHalfHeight(0,0)+close) ||
		      spacepointss[0].XYZ()[1] > (2.*geom->DetHalfHeight(0,0)-close) || spacepointss[0].XYZ()[1] < (-2.*geom->DetHalfHeight(0,0)+close) ||
		      spacepointss[spacepointss.size()-1].XYZ()[2] > (geom->DetLength(0,0)-close) || spacepointss[spacepointss.size()-1].XYZ()[2] < close ||
		      spacepointss[0].XYZ()[2] > (geom->DetLength(0,0)-close) || spacepointss[0].XYZ()[2] < close
		      )
		    uncontained = true; 

		  if (uncontained) 
		    {		      
		      // Big enough to not run out of gas right at end of
		      // track and give large angular deviations which
		      // will kill the fit.
		      mom.SetMag(3.0 * mom.Mag()); 
		      std::cout<<"Track3DKalmanSPS: Uncontained track ... Run "<<evt.run()<<" Event "<<evt.id().event()<<std::endl;
		    }
		  fcont = (int) (!uncontained);
		  TVector3 momM(mom);
		  TVector3 momErrFit(momM[0]/100.0,
		  		     momM[1]/100.0,
		  		     momM[2]/100.0);   // GeV
	    
		  genf::GFFieldManager::getInstance()->init(new genf::GFConstField(0.,0.,0.0));
		  genf::GFDetPlane planeG((TVector3)(spacepointss[0].XYZ()),momM);
	    

		  //      std::cout<<"Track3DKalmanSPS about to do GAbsTrackRep."<<std::endl;
		  // Initialize with 1st spacepoint location and ...
		  rep = new genf::RKTrackRep(//posM-.5/momM.Mag()*momM,
					     (TVector3)(spacepointss[0].XYZ()),
					     momM,
					     posErr,
					     momErrFit,
					     -13);  // mu+ hypothesis
		  //      std::cout<<"Track3DKalmanSPS: about to do GFTrack. repDim is " << rep->getDim() <<std::endl;


		  genf::GFTrack fitTrack(rep);//initialized with smeared rep
		  // Gonna sort in z cuz I want to essentially transform here to volTPC coords.
		  // volTPC coords, cuz that's what the Geant3/Geane stepper wants, as that's its understanding
		  // from the Geant4 geometry, which it'll use. EC, 7-Jan-2011.
		  int ihit = 0;
		  fptsNo = 0;
		  for (unsigned int point=0;point<spacepointss.size();++point)
		    {
		      double sep;
		      // Calculate the distance in 2nd and 3rd PCs and
		      // reject spt if it's too far out. Remember, the 
		      // sigmas are sqrt(eigenvals).
		      double tmp[3];
		      principal->X2P((Double_t *)(spacepointss[point].XYZ()),tmp);
		      sep = sqrt(tmp[1]*tmp[1]/fPCevals[1]+tmp[2]*tmp[2]/fPCevals[2]);
		      if ((fabs(sep) > fPerpLim) && (point<(spacepointss.size()-nTailPoints)))
		      {
			//			std::cout << "Spacepoint " << point << " DROPPED!!!:" << spacepointss[point].XYZ()[0]<< ", " << spacepointss[point].XYZ()[1]<< ", " << spacepointss[point].XYZ()[2]<< ". " << std::endl;
			continue;
		      }
		      if (point%fDecimate) // Jump out of loop except on every fDecimate^th pt. fDecimate==1 never sees continue.
			{
			  continue;
			}
		      TVector3 spt3 = (TVector3)(spacepointss[point].XYZ());
		      std::vector <double> err3;
		      err3.push_back(spacepointss[point].ErrXYZ()[0]);
		      err3.push_back(spacepointss[point].ErrXYZ()[2]);
		      err3.push_back(spacepointss[point].ErrXYZ()[4]);
		      err3.push_back(spacepointss[point].ErrXYZ()[5]); // lower triangle diags.
		      if (fptsNo<fDimSize)
			{
			  fshx[fptsNo] = spt3[0];
			  fshy[fptsNo] = spt3[1];
			  fshz[fptsNo] = spt3[2];
			  feshx[fptsNo] = err3[0];
			  feshy[fptsNo] = err3[1];
			  feshz[fptsNo] = err3[3];
			  feshyz[fptsNo] = err3[2];
			  fsep[fptsNo] = sep;
			  if (fptsNo>1)
			    {
			      TVector3 pointer(fshx[fptsNo]-fshx[fptsNo-1],fshy[fptsNo]-fshy[fptsNo-1],fshz[fptsNo]-fshz[fptsNo-1]);
			      TVector3 pointerPrev(fshx[fptsNo-1]-fshx[fptsNo-2],fshy[fptsNo-1]-fshy[fptsNo-2],fshz[fptsNo-1]-fshz[fptsNo-2]);
			      fth[fptsNo] = (pointer.Unit()).Angle(pointerPrev.Unit());
			    }
			  feth[fptsNo] = 0.0;
			  fedudw[fptsNo] = 0.0;
			  fedvdw[fptsNo] = 0.0;
			  feu[fptsNo] = 0.0;
			  fev[fptsNo] = 0.0;
			  fupdate[fptsNo] = 0.0;
			}


		      mf::LogDebug("Track3DKalmanSPS: ") << "ihit xyz..." << spt3[0]<<","<< spt3[1]<<","<< spt3[2];

		      fitTrack.addHit(new genf::PointHit(spt3,err3),
				      1,//dummy detector id
				      ihit++
				      );
		      fptsNo++;
		    } // end loop over spacepoints.

		  mf::LogInfo("Track3DKalmanSPS: ") << "Fitting on " << fptsNo << " points.";
		  //      std::cout<<"Track3DKalmanSPS about to do GFKalman."<<std::endl;
		  genf::GFKalman k;
		  k.setBlowUpFactor(500); // 500 out of box. EC, 6-Jan-2011.
		  k.setMomHigh(20.0); // Don't fit above this many GeV.
		  k.setMomLow(0.1);   // Don't fit below this many GeV.
		  k.setMaxUpdate(0.1); // 0 out abs(update) bigger than this.
		  k.setInitialDirection(+1); // Instead of 1 out of box. EC, 6-Jan-2011.
		  k.setNumIterations(numIT);
		  bool skipFill = false;
		  //      std::cout<<"Track3DKalmanSPS back from setNumIterations."<<std::endl;
		  try{
		    //	std::cout<<"Track3DKalmanSPS about to processTrack."<<std::endl;
		    if (fDoFit) k.processTrack(&fitTrack);
		    std::vector < TMatrixT<Double_t> > measCov(fitTrack.getMeasuredCov());
		    std::vector < TMatrixT<Double_t> > measUpdate(fitTrack.getMeasuredUpdate());
		    for (unsigned int ihit=0; ihit<fptsNo; ihit++)
		      {
			feth[ihit] = (Float_t ) (measCov.at(ihit)[0][0]); // eth
			fedudw[ihit] = (Float_t ) (measCov.at(ihit)[1][1]); 
			fedvdw[ihit] = (Float_t ) (measCov.at(ihit)[2][2]); 
			feu[ihit] = (Float_t ) (measCov.at(ihit)[3][3]); 
			fev[ihit] = (Float_t ) (measCov.at(ihit)[4][4]);
			fupdate[ihit] = (Float_t ) (measUpdate.at(ihit)[0][0]);
		      }
		    //std::cout<<"Track3DKalmanSPS back from processTrack."<<std::endl;
		  }
		  //catch(GFException& e){
		  catch(cet::exception &e){
		    mf::LogError("Track3DKalmanSPS: ") << "just caught a cet::exception."<<std::endl;
		    e.what();
		    mf::LogError("Track3DKalmanSPS: ") << "Exceptions won't be further handled, line: "<<__LINE__;
		    mf::LogError("Track3DKalmanSPS: ") << "Skip filling big chunks of the TTree, line: "<<__LINE__;
		    skipFill = true;
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
		      
		      if (!skipFill)
			{
			  stREC->ResizeTo(rep->getState());
			  *stREC = rep->getState();
			  covREC->ResizeTo(rep->getCov());
			  *covREC = rep->getCov();
			  double dum[5];
			  double dum2[5];
			  for (unsigned int ii=0;ii<5;ii++)
			    {
			      stREC->ExtractRow(ii,0,dum);
			      fState0[ii] = dum[0];
			      covREC->ExtractRow(ii,0,dum2);
			      for (unsigned int jj=0;jj<5;jj++)
				{
				  fCov0[ii*5+jj] = dum2[jj];
				}
			    }
			  if(fGenfPRINT)
			    {
			      mf::LogInfo("Track3DKalmanSPS: ") << " First State and Cov:";
			      stREC->Print();
			      covREC->Print();
			    }
			  chi2 = (Float_t)(rep->getChiSqu());
			  ndf = rep->getNDF();
			  nfail = fitTrack.getFailedHits();
			  chi2ndf = (Float_t)(chi2/ndf);
			  double dircoss[3],dircose[3];
			  //		(*trackIter)->Direction(dircoss,dircose);      
			  for (int ii=0;ii<3;++ii)
			    {
			      fpMCMom[ii] = MCMomentum[ii];
			      fpMCPos[ii] = MCOrigin[ii];
			      fpREC[ii] = rep->getMom(rep->getReferencePlane())[ii];
			      // fpRECLE is an extrap forward from first hit. Captures MS'ing 
			      // only macroscopically. Not really what's desired.
			      // fpRECL is actually the last plane's momentum estimate.
			      //fpRECLE[ii] = rep->getMom(rep->getLastPlane())[ii]; 
			      fpRECL[ii] = (dynamic_cast<genf::RKTrackRep *> (rep))->getMomLast(rep->getLastPlane())[ii]; 
			      fpRECt3D[ii] = dircoss[ii];
			    }
			  fpMCMom[3] = MCMomentum.Mag();
			  
			  fpREC[3] = rep->getMom(rep->getReferencePlane()).Mag();
			  fpRECL[3] = (dynamic_cast<genf::RKTrackRep *> (rep))->getMomLast(rep->getLastPlane()).Mag();
			  //fpRECLE[3] = rep->getMom(rep->getLastPlane()).Mag();

			  nTrks++;

			  mf::LogInfo("Track3DKalmanSPS: ") << "Track3DKalmanSPS about to do tree->Fill(). Chi2/ndf is " << chi2/ndf << ". All in volTPC coords .... pMCT[0-3] is " << fpMCMom[0] << ", " << fpMCMom[1] << ", " << fpMCMom[2] << ", " << fpMCMom[3] << ". pREC[0-3] is " << fpREC[0] << ", "<< fpREC[1] << ", " << fpREC[2] << ", " << fpREC[3] <<  ". pRECL[0-3] is " << fpRECL[0] << ", "<< fpRECL[1] << ", " << fpRECL[2] << ", " << fpRECL[3] << ".";
		      
			} // end !skipFill

		      evtt = (unsigned int) evt.id().event();
		      tree->Fill();
		      
		      art::PtrVector<recob::Cluster> clusters;
		      // Use Assns to get the clusters for the spacepoints.
		      clusters = util::FindManyP<recob::Cluster>(prongIn, evt, fProngModuleLabel, cntp);
		      cntp++;

		      recob::Track  the3DTrack(clusters,spacepointss);
		      double dircosF[3];
		      double dircosL[3];
		      
		      for (int ii=0;ii<3;++ii)
			{
			  dircosF[ii] = fpREC[ii]/fpREC[3];
			  dircosL[ii] = fpRECL[ii]/fpRECL[3];
			}
		      // Add the 3D track to the vector of the reconstructed tracks		      
		      the3DTrack.SetDirection(dircosF,dircosL);
		      the3DTrack.SetID(tcnt++);
		      tcol->push_back(the3DTrack);
		      util::CreateAssn(*this, evt, *(tcol.get()), clusters,*(assn.get()));
		      
		    } // getStatusFlag 

		  if (!rep) delete rep;

	  pprong++;
	  
	} // end while on prongs.

      if (!repMC) delete repMC;

      if (tcol->size()) 
	{ 
	  evt.put(tcol); 
	  evt.put(assn);
	}
}
