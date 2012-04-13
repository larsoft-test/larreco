//
// Name: SeedFinder.cxx
//
// Purpose: Implementation file for module SeedFinder.
//
// Ben Jones, MIT
//

#include <vector>
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "TrackFinder/SeedFinder.h"
#include "TrackFinder/SpacePointService.h"
#include "Geometry/geo.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "RecoBase/Hit.h"
#include "RecoBase/Seed.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Track.h"
#include "RecoBase/Prong.h"
#include "Utilities/AssociationUtil.h"

namespace trkf {

  SeedFinder::SeedFinder(const fhicl::ParameterSet& pset) :
    fFilter(true),
    fMerge(false)
  {
    reconfigure(pset);
    produces<std::vector<recob::Seed> >();
  
  
    mf::LogInfo("SeedFinder") 
      << "SeedFinder configured with the following parameters:\n"
      << "  ClusterModuleLabel = " << fClusterModuleLabel << "\n"
      << "  Filter = " << fFilter << "\n"
      << "  Merge = " << fMerge;
  }

  SeedFinder::~SeedFinder()
  {
  }

  void SeedFinder::reconfigure(fhicl::ParameterSet const& pset)
  {
    fClusterModuleLabel    = pset.get<std::string>("ClusterModuleLabel");
    fFilter                = pset.get<bool>("Filter");
    fMerge                 = pset.get<bool>("Merge");
    fSeedMode              = pset.get<int>("SeedMode");
    fSeedLength            = pset.get<double>("SeedLength");
    fMinPointsInCluster    = pset.get<unsigned int>("MinPointsInCluster");
    fMinPointsInSeed       = pset.get<unsigned int>("MinPointsInSeed");
    fAngularDev            = pset.get<double>("AngularDev");


  }

  void SeedFinder::beginJob()
  {}

  void SeedFinder::produce(art::Event& evt)
  {
    
    std::auto_ptr<std::vector<recob::Seed> > seeds(new std::vector<recob::Seed>);

    std::vector<std::vector<recob::SpacePoint> > SpacePointVectors = GetSpacePointsFromClusters(fClusterModuleLabel, evt);
    if(SpacePointVectors.size() > 0)
      {
	for(std::vector<std::vector<recob::SpacePoint> >::const_iterator it=SpacePointVectors.begin(); 
	    it!=SpacePointVectors.end();
	    it++)
	  {
	    std::vector<recob::Seed*> SeedFinderOutput;
	    if(fSeedMode==0)
	      SeedFinderOutput = FindSeedExhaustively(*it);
	    else if(fSeedMode==1)
	      SeedFinderOutput = FindAsManySeedsAsPossible(*it);
	    else 
	      throw cet::exception("SeedFinder") << 
		"Unkown seed mode " << fSeedMode<<"\n";

	    if(SeedFinderOutput.size()>0)
	      for(unsigned int i=0; i!=SeedFinderOutput.size(); i++)
		if(SeedFinderOutput.at(i)->IsValid())
		  seeds->push_back(*(SeedFinderOutput.at(i)));	    
	  }			  
      }
    else
      std::cout<<"Seed finder made no seeds : no space points in event"<<std::endl;
    
    evt.put(seeds);
    
  }



  // Given a cluster label, produce spacepoints from the event
  //  very possible this spacepoint supply will be swapped out later.

  std::vector<std::vector<recob::SpacePoint> > SeedFinder::GetSpacePointsFromClusters(std::string ClusterModuleLabel, art::Event& evt)
  {
    // Get Services.

    art::ServiceHandle<trkf::SpacePointService> sptsvc;
    art::ServiceHandle<geo::Geometry> geom;

    art::Handle< std::vector<recob::Cluster> > clusterh;
    evt.getByLabel(ClusterModuleLabel, clusterh);

    // Get clusters.

    std::vector<std::vector<recob::SpacePoint> > SpacePointVectors;

    // Make a double or triple loop over clusters in distinct views
    // (depending on minimum number of views configured in SpacePointService).

    if(clusterh.isValid()) {



      art::PtrVector<recob::Hit> hits;      

      // Loop over first cluster.

      int nclus = clusterh->size();
      for(int iclus = 0; iclus < nclus; ++iclus) {
	art::Ptr<recob::Cluster> piclus(clusterh, iclus);
	geo::View_t iview = piclus->View();

	// Test first view.

	if((iview == geo::kU && sptsvc->enableU()) ||
	   (iview == geo::kV && sptsvc->enableV()) ||
	   (iview == geo::kW && sptsvc->enableW())) {

	  // Store hits from first view into hit vector.

	  art::PtrVector<recob::Hit> ihits = piclus->Hits();
	  unsigned int nihits = ihits.size();
	  hits.clear();
	  hits.reserve(nihits);
	  for(art::PtrVector<recob::Hit>::const_iterator i = ihits.begin();
	      i != ihits.end(); ++i)
	    hits.push_back(*i);
	  
	  // Loop over second cluster.

	  for(int jclus = 0; jclus < iclus; ++jclus) {
	    art::Ptr<recob::Cluster> pjclus(clusterh, jclus);
	    geo::View_t jview = pjclus->View();

	    // Test second view.

	    if(((jview == geo::kU && sptsvc->enableU()) ||
		(jview == geo::kV && sptsvc->enableV()) ||
		(jview == geo::kW && sptsvc->enableW()))
	       && jview != iview) {

	      // Store hits from second view into hit vector.

	      art::PtrVector<recob::Hit> jhits = pjclus->Hits();
	      unsigned int njhits = jhits.size();
	      assert(hits.size() >= nihits);
	      //hits.resize(nihits);
	      while(hits.size() > nihits)
		hits.pop_back();
	      assert(hits.size() == nihits);
	      hits.reserve(nihits + njhits);
	      for(art::PtrVector<recob::Hit>::const_iterator j = jhits.begin();
		  j != jhits.end(); ++j)
		hits.push_back(*j);
	  

	      // Loop over third cluster.

	      for(int kclus = 0; kclus < jclus; ++kclus) {
		art::Ptr<recob::Cluster> pkclus(clusterh, kclus);
		geo::View_t kview = pkclus->View();

		// Test third view.

		if(((kview == geo::kU && sptsvc->enableU()) ||
		    (kview == geo::kV && sptsvc->enableV()) ||
		    (kview == geo::kW && sptsvc->enableW()))
		   && kview != iview && kview != jview) {

		  // Store hits from third view into hit vector.

		  art::PtrVector<recob::Hit> khits = pkclus->Hits();
		  unsigned int nkhits = khits.size();
		  assert(hits.size() >= nihits + njhits);
		  //hits.resize(nihits + njhits);
		  while(hits.size() > nihits + njhits)
		    hits.pop_back();
		  assert(hits.size() == nihits + njhits);
		  hits.reserve(nihits + njhits + nkhits);
		  for(art::PtrVector<recob::Hit>::const_iterator k = khits.begin();
		      k != khits.end(); ++k)
		    hits.push_back(*k);

		  // Make three-view space points.

		  std::vector<recob::SpacePoint> spts;
		  sptsvc->makeSpacePoints(hits, spts,
					  fFilter, fMerge, 0., 0.);

		  if(spts.size() > 0) {
		    SpacePointVectors.push_back(spts);
		  }
		}
	      }
	    }
	  }
	}
      }

    }
    return SpacePointVectors;
  }





  //
  // Take a collection of spacepoints and return a
  // vector of as many straight line seeds as possible
  //
  std::vector<recob::Seed *> SeedFinder::FindAsManySeedsAsPossible(std::vector<recob::SpacePoint> Points)
  {
    std::vector<recob::Seed*> ReturnVector;

    recob::Seed* TrackSeedThisCombo;
    bool KeepChopping=true;
    if(Points.size()<fMinPointsInCluster) KeepChopping=false;
    while(KeepChopping)
      {
	TrackSeedThisCombo = FindSeedAtEnd(Points);
	if(TrackSeedThisCombo->IsValid())
	  {
	    ReturnVector.push_back(TrackSeedThisCombo);
	  }
	// if enough left, chop off some points and try again
	if(Points.size() > fMinPointsInSeed)
	  Points.resize(Points.size()-fMinPointsInSeed);
	else  // or if there is nothing left to chop, leave the loop
	  KeepChopping=false;

      }

    return ReturnVector;
  }






  //
  // Take a collection of spacepoints sorted in Z and find
  // exactly one straight seed, as high in Z as possible.
  //
  std::vector<recob::Seed*> SeedFinder::FindSeedExhaustively(std::vector<recob::SpacePoint> Points)
  {
    recob::Seed* TrackSeedThisCombo;
    bool KeepChopping=true;
    bool FoundValidSeed=false;
    if(Points.size()<=fMinPointsInCluster) KeepChopping=false;
    while(KeepChopping)
      {
	TrackSeedThisCombo = FindSeedAtEnd(Points);
	if(TrackSeedThisCombo->IsValid())
	  {
	    KeepChopping=false;
	    FoundValidSeed=true;
	  }
	else
	  {
	    // if we found an invalid or weak seed, chop off some points and try again
  if(Points.size() > fMinPointsInSeed)
    Points.resize(Points.size()-fMinPointsInSeed);
  else  // or if there is nothing left to chop, leave the loop
    KeepChopping=false;
	  }
      }
    if(!FoundValidSeed)
      {
	TrackSeedThisCombo = new recob::Seed();
      }
    else
      {
	double SeedPoint[3], SeedErr[3];
	TrackSeedThisCombo->GetPoint(SeedPoint,SeedErr);
	std::cout<<"FindSeedEx : " << SeedPoint[0]<< " " <<
	  SeedPoint[1]<< " " <<
	  SeedPoint[2]<< " " << std::endl;
      }
    std::vector<recob::Seed*> ReturnVector;
    ReturnVector.push_back(TrackSeedThisCombo);
    return ReturnVector;
  }


  //
  // Start with a vector of hits ordered in Z.
  // Put together a seed using the last N hits, find its
  // centre, direction and strength. Return this seed for
  // further scrutiny.
  //
  recob::Seed * SeedFinder::FindSeedAtEnd(std::vector<recob::SpacePoint> Points)
  {

    recob::Seed * ReturnSeed;

    std::vector<int> PointsInRange;
    PointsInRange.clear();

    TVector3 HighestZPoint( Points.at(Points.size()-1).XYZ()[0],
			    Points.at(Points.size()-1).XYZ()[1],
			    Points.at(Points.size()-1).XYZ()[2] );

    TVector3 CentreOfPoints(0.,0.,0.);

    for(int index=Points.size()-1; index!=-1; --index)
      TVector3 CentreOfPoints(0.,0.,0.);

    for(int index=Points.size()-1; index!=-1; --index)
      {
	// first check z, then check total distance
	//  (much faster, since most will be out of range in z anyway)


	if( ( HighestZPoint[2] - Points.at(index).XYZ()[2] ) < fSeedLength)
	  {
	    if( ( ( pow(HighestZPoint[0]-Points.at(index).XYZ()[0],2) +
		    pow(HighestZPoint[1]-Points.at(index).XYZ()[1],2) +
		    pow(HighestZPoint[2]-Points.at(index).XYZ()[2],2) ) < pow(fSeedLength,2)))
	      {
		PointsInRange.push_back(index);
		CentreOfPoints[0]+=Points.at(index).XYZ()[0];
		CentreOfPoints[1]+=Points.at(index).XYZ()[1];
		CentreOfPoints[2]+=Points.at(index).XYZ()[2];
	      }
	  }
      }

    unsigned int NPoints = PointsInRange.size();
    if(NPoints==0) return new recob::Seed();
    CentreOfPoints = (1./float(NPoints)) * CentreOfPoints;

    if(NPoints>fMinPointsInSeed)
      {
	float costheta=0, phi=0, costheta2=0, phi2=0;

	for(unsigned int i=0; i!=PointsInRange.size(); i++)
	  {
	    TVector3 ThisPoint(Points.at(PointsInRange.at(i)).XYZ()[0],
			       Points.at(PointsInRange.at(i)).XYZ()[1],
			       Points.at(PointsInRange.at(i)).XYZ()[2]);
	    TVector3 Step = CentreOfPoints-ThisPoint;

	    if(Step.Z()<0) Step=-Step;


	    phi       += Step.Phi();
	    phi2      += pow(Step.Phi(),2);

	    // Need to check range of costheta to avoid getting a nan.
	    //  (if CosTheta=1, floating point errors mess us up)
	    if(Step.CosTheta()<0.9999)
	      {
		costheta     += Step.CosTheta();
		costheta2    += pow(Step.CosTheta(),2);
	      }
	    else
	      {
		costheta+=1;
		costheta2+=1;
	      }
	    //      std::cout<<phi<< " " <<phi2<< " " <<costheta << " " <<costheta2<\
	      <std::endl;
	  }


	float sigtheta    = pow((costheta2 / NPoints - pow( costheta/NPoints, 2 )), 0.5);
	float sigphi      = pow((phi2      / NPoints - pow( phi/NPoints,   2 )),    0.5);
	float meantheta   = acos(costheta      / NPoints);
	float meanphi     = phi                / NPoints;

	// Total angular deviation (spherical polars)
	float AngularDev =  pow(pow(sigtheta,2) + pow(sin(sigtheta) * sigphi, 2), 0.5);

	//	std::cout<<"SeedFinder debug: " << AngularDev<<std::endl;

	if(AngularDev<fAngularDev)
	  {

	    double PtArray[3], DirArray[3];

	    PtArray[0]=CentreOfPoints.X();
	    PtArray[1]=CentreOfPoints.Y();
	    PtArray[2]=CentreOfPoints.Z();

	    // calculate direction from average over hits

	    TVector3 SeedDirection;

	    SeedDirection.SetMagThetaPhi(fSeedLength, meantheta, meanphi);

	    DirArray[0]=SeedDirection.X();
	    DirArray[1]=SeedDirection.Y();
	    DirArray[2]=SeedDirection.Z();

	    ReturnSeed = new recob::Seed(PtArray,DirArray);
	  }
        else ReturnSeed = new recob::Seed();
      }
    else ReturnSeed = new recob::Seed();

    return ReturnSeed;
  }











  void SeedFinder::endJob()
  {

  }
}
