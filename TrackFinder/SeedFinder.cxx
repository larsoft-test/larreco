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
#include "Geometry/geo.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Seed.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Track.h"
#include "RecoBase/Prong.h"
#include "Utilities/AssociationUtil.h"
#include "Utilities/DetectorProperties.h"
#include "TMatrixD.h"
#include "TVectorD.h"

namespace trkf {

  SeedFinder::SeedFinder(const fhicl::ParameterSet& pset) :
    fSptalg(pset.get<fhicl::ParameterSet>("SpacePointAlg")),
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
    fhicl::ParameterSet seedConfig = pset.get<fhicl::ParameterSet>("SeedConfig");

    fSptalg                = seedConfig.get<fhicl::ParameterSet>("SpacePointAlg");
    fClusterModuleLabel    = seedConfig.get<std::string>("ClusterModuleLabel");
    fHitModuleLabel        = seedConfig.get<std::string>("HitModuleLabel");
    fFilter                = seedConfig.get<bool>("Filter");
    fMerge                 = seedConfig.get<bool>("Merge");
    fSeedMode              = seedConfig.get<int>("SeedMode");
    fSource                = seedConfig.get<int>("Source");
    fSeedLength            = seedConfig.get<double>("SeedLength");
    fMinPointsInCluster    = seedConfig.get<unsigned int>("MinPointsInCluster");
    fMinPointsInSeed       = seedConfig.get<unsigned int>("MinPointsInSeed");
    fAngularDev            = seedConfig.get<double>("AngularDev");
    fRefits                = seedConfig.get<double>("Refits");
  }

  void SeedFinder::beginJob()
  {}

  void SeedFinder::produce(art::Event& evt)
  {
    
    std::auto_ptr<std::vector<recob::Seed> > seeds(new std::vector<recob::Seed>);

    std::vector<std::vector<recob::SpacePoint> > SpacePointVectors;
    if(fSource==1)
      SpacePointVectors = GetSpacePointsFromClusters(fClusterModuleLabel, evt);
    else if(fSource==2)
      {
	SpacePointVectors = GetSpacePointsFromHits(fHitModuleLabel, evt);
	std::cout<<"Getting space points from hits"<<std::endl;
      }
    else	  
      {
	throw cet::exception("SeedFinder") << 
	  "Unkown source mode " << fSource<<"\n";
      }
    if(SpacePointVectors.size() > 0)
      {
	for(std::vector<std::vector<recob::SpacePoint> >::const_iterator it=SpacePointVectors.begin(); 
	    it!=SpacePointVectors.end();
	    it++)
	  {
	    std::vector<std::vector<recob::SpacePoint> > PointsUsed;
	    std::vector<recob::Seed*>                    SeedFinderOutput;
	    if(fSeedMode==0)
	      SeedFinderOutput = FindSeedExhaustively(*it, PointsUsed);
	    else if(fSeedMode==1)
	      SeedFinderOutput = FindAsManySeedsAsPossible(*it, PointsUsed);
	    else 
	      throw cet::exception("SeedFinder") << 
		"Unkown seed mode " << fSeedMode<<"\n";

	   
	    if(SeedFinderOutput.size()>0)
	      {
		PointsUsed.resize(SeedFinderOutput.size());
		for(unsigned int i=0; i!=SeedFinderOutput.size(); i++)
		if(SeedFinderOutput.at(i)->IsValid())
		  {
		    if(fRefits>0)
		      RefitSeed(SeedFinderOutput.at(i), PointsUsed.at(i));		    
		    seeds->push_back(*SeedFinderOutput.at(i));
		  }
	      }
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

    art::ServiceHandle<geo::Geometry> geom;

    art::Handle< std::vector<recob::Cluster> > clusterh;
    evt.getByLabel(ClusterModuleLabel, clusterh);

    // Get clusters.

    std::vector<std::vector<recob::SpacePoint> > SpacePointVectors;

    // Make a double or triple loop over clusters in distinct views
    // (depending on minimum number of views configured in SpacePointAlg).

    if(clusterh.isValid()) {



      art::PtrVector<recob::Hit> hits;      

      // Loop over first cluster.

      int nclus = clusterh->size();
      for(int iclus = 0; iclus < nclus; ++iclus) {
	art::Ptr<recob::Cluster> piclus(clusterh, iclus);
	geo::View_t iview = piclus->View();

	// Test first view.

	if((iview == geo::kU && fSptalg.enableU()) ||
	   (iview == geo::kV && fSptalg.enableV()) ||
	   (iview == geo::kW && fSptalg.enableW())) {

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

	    if(((jview == geo::kU && fSptalg.enableU()) ||
		(jview == geo::kV && fSptalg.enableV()) ||
		(jview == geo::kW && fSptalg.enableW()))
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

		if(((kview == geo::kU && fSptalg.enableU()) ||
		    (kview == geo::kV && fSptalg.enableV()) ||
		    (kview == geo::kW && fSptalg.enableW()))
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
		  fSptalg.makeSpacePoints(hits, spts,
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


  // Extract vector of hits from event

  art::PtrVector<recob::Hit> SeedFinder::GetHitsFromEvent(std::string HitModuleLabel, art::Event & evt)
  {
    art::PtrVector<recob::Hit> TheHits;
    art::Handle< std::vector<recob::Hit> > hith;
    evt.getByLabel(HitModuleLabel, hith);
    for(unsigned int i=0; i<hith->size(); ++i)
      {
	art::Ptr<recob::Hit> hit(hith, i);
	TheHits.push_back(hit);
      }
    
    return TheHits;
  }


  // Generate spacepoints from hits in event to feed SeedFinder
  std::vector<std::vector<recob::SpacePoint> > SeedFinder::GetSpacePointsFromHits(std::string HitModuleLabel, art::Event& evt)
  {
    std::vector<std::vector<recob::SpacePoint> > ReturnVec;
    art::PtrVector<recob::Hit> TheHits = GetHitsFromEvent(HitModuleLabel, evt);
    std::vector<recob::SpacePoint> TheSPs = GetSpacePointsFromHitVector(TheHits);
    if(TheSPs.size()>0)
      {
	ReturnVec.push_back(TheSPs);
      }
    return ReturnVec;
  }
  


  //
  // Given a collection of hits, make space points.
  // The sp vector fills the first entry of a nested
  // vector to keep compatability with the cluster
  // based approach
  //
  
  std::vector<recob::SpacePoint>  SeedFinder::GetSpacePointsFromHitVector(art::PtrVector<recob::Hit> Hits)
  {
    std::vector<recob::SpacePoint> ReturnVec;
    
    fSptalg.makeSpacePoints(Hits, ReturnVec,
    			    fFilter, fMerge, 0., 0.);

 
    return ReturnVec;
      
  }



  //
  // Take a collection of spacepoints and return a
  // vector of as many straight line seeds as possible
  //
  std::vector<recob::Seed *> SeedFinder::FindAsManySeedsAsPossible(std::vector<recob::SpacePoint> Points, std::vector<std::vector<recob::SpacePoint> > & PointsUsed)
  {
    std::vector<recob::Seed*> ReturnVector;
    PointsUsed.clear();

    recob::Seed* TrackSeedThisCombo;
    bool KeepChopping=true;
    if(Points.size()<fMinPointsInCluster) KeepChopping=false;
    while(KeepChopping)
      {
	std::vector<int> ToChop;
	std::vector<recob::SpacePoint> SPs;

	TrackSeedThisCombo = FindSeedAtEnd(Points,SPs,ToChop);
	if(TrackSeedThisCombo->IsValid())
	  {
	    ReturnVector.push_back(TrackSeedThisCombo);
	    PointsUsed.push_back(SPs);
	  }
	// if enough left, chop off some points and try again
	if((Points.size()-ToChop.size()) > fMinPointsInSeed)
	  {
	    sort(ToChop.begin(), ToChop.end());
	    for(int i=(int)ToChop.size()-1; i!=-1; i--) 
	      {
		Points.erase(Points.begin()+ToChop.at(i));
	      }
	  }
	else  // or if there is nothing left to chop, leave the loop
	  KeepChopping=false;

      }

    return ReturnVector;
  }






  //
  // Take a collection of spacepoints sorted in Z and find
  // exactly one straight seed, as high in Z as possible.
  //
  std::vector<recob::Seed*> SeedFinder::FindSeedExhaustively(std::vector<recob::SpacePoint> Points, std::vector<std::vector<recob::SpacePoint> >& PointsUsed)
  {
    recob::Seed* TrackSeedThisCombo;
    bool KeepChopping=true;
    bool FoundValidSeed=false;
    if(Points.size()<=fMinPointsInCluster) KeepChopping=false;
    while(KeepChopping)
      {
	std::vector<recob::SpacePoint> SPs;
	TrackSeedThisCombo = FindSeedAtEnd(Points, SPs);
	if(TrackSeedThisCombo->IsValid())
	  {
	    KeepChopping=false;
	    PointsUsed.push_back(SPs);
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
  recob::Seed * SeedFinder::FindSeedAtEnd(std::vector<recob::SpacePoint> Points, std::vector<recob::SpacePoint>& PointsUsed)
  {
    PointsUsed.clear();
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
	    if( ( pow(HighestZPoint[0]-Points.at(index).XYZ()[0],2) + 
		  pow(HighestZPoint[1]-Points.at(index).XYZ()[1],2) + 
		  pow(HighestZPoint[2]-Points.at(index).XYZ()[2],2) ) < pow(fSeedLength,2))
	      {
		PointsInRange.push_back(index);
		CentreOfPoints[0]+=Points.at(index).XYZ()[0];
		CentreOfPoints[1]+=Points.at(index).XYZ()[1];
		CentreOfPoints[2]+=Points.at(index).XYZ()[2];
		PointsUsed.push_back(Points.at(index));
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



  //
  // Start with a vector of hits ordered in Z.
  // Put together a seed using the last N hits, find its
  // centre, direction and strength. Return this seed for
  // further scrutiny.
  //
recob::Seed * SeedFinder::FindSeedAtEnd(std::vector<recob::SpacePoint> Points,std::vector<recob::SpacePoint>& PointsUsed, std::vector<int>& ToThrow)
  {

    recob::Seed * ReturnSeed;
    PointsUsed.clear();
    std::vector<int> PointsInRange;
    PointsInRange.clear();
    ToThrow.clear();
    
    TVector3 HighestZPoint( Points.at(Points.size()-2).XYZ()[0],
			    Points.at(Points.size()-2).XYZ()[1],
			    Points.at(Points.size()-2).XYZ()[2] );

  
    TVector3 CentreOfPoints(0.,0.,0.);

    for(int index=Points.size()-1; index!=-1; --index)
      {
	// first check z, then check total distance
	//  (much faster, since most will be out of range in z anyway)

	if( ( HighestZPoint[2] - Points.at(index).XYZ()[2] ) < fSeedLength)
	  {
	    double DistanceToHighZ =pow( pow(HighestZPoint[0]-Points.at(index).XYZ()[0],2) +
				      pow(HighestZPoint[1]-Points.at(index).XYZ()[1],2) +
					 pow(HighestZPoint[2]-Points.at(index).XYZ()[2],2),0.5 ); 
	      if( DistanceToHighZ < fSeedLength)
		{
		  PointsInRange.push_back(index);
		  CentreOfPoints[0]+=Points.at(index).XYZ()[0];
		  CentreOfPoints[1]+=Points.at(index).XYZ()[1];
		  CentreOfPoints[2]+=Points.at(index).XYZ()[2];
		  PointsUsed.push_back(Points.at(index));
		}
	  }
      }
    
    
    unsigned int NPoints = PointsInRange.size();
    if(NPoints==0) return new recob::Seed();
    CentreOfPoints = (1./float(NPoints)) * CentreOfPoints;

    std::cout<<"Trying seed at " << CentreOfPoints[0]<<" " <<
      CentreOfPoints[1]<<" " << CentreOfPoints[2]<<" " <<NPoints<<std::endl;

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
    
    ToThrow=PointsInRange;
    
    return ReturnSeed;
  }

  

  void SeedFinder::RefitSeed(recob::Seed * TheSeed, std::vector<recob::SpacePoint> SpacePoints)
  {

    std::cout<<"Refit module called on vector of " << SpacePoints.size() << " space points " << std::endl;

    // Get the services we need
    art::ServiceHandle<geo::Geometry>            geom;
    art::ServiceHandle<util::DetectorProperties> det;

    if((geom->NTPC()!=1)||geom->NTPC()!=1)
      {
	throw cet::exception("SeedFinder : Refit only works for 1 tpc, 1 cryostat detector")<<
	  "seed refitting feature not yet developped for multi cryostat or "<<
	  "multi TPC detector - see TrackFinder/SeedFinder.cxx"<<std::endl;	
      }
    size_t Planes = geom->TPC(0).Nplanes();
   
   
    // Get this hits in each view that made this seed

    std::map<int,art::PtrVector<recob::Hit> > HitMap;
    std::map<unsigned short, bool>            HitsClaimed;


    // for each spacepoint
    for(std::vector<recob::SpacePoint>::const_iterator itSP=SpacePoints.begin();
	itSP!=SpacePoints.end(); itSP++)
      {
	// get hits from each plane (ensuring each used once only)
	for(size_t plane=0; plane!=Planes; plane++)
	  {
	    art::PtrVector<recob::Hit> HitsThisSP = itSP->Hits(plane);
	    for(art::PtrVector<recob::Hit>::const_iterator itHit=HitsThisSP.begin();
		itHit!=HitsThisSP.end(); itHit++)
	      {
		if(!HitsClaimed[(*itHit)->Channel()])
		  {
		    HitMap[plane].push_back(*itHit);
		    HitsClaimed[(*itHit)->Channel()]=true;
		  }
		
	      }
	  }
      }


    // Begin iterative refit procedure

    // Loop the prescribed number of times
    for(int loop=0; loop!=fRefits; loop++)
      {
	std::cout<<"SeedFinder running refit " << loop<<std::endl;
	for(size_t plane=0; plane!=Planes; plane++)
	  {
	    art::PtrVector<recob::Hit> HitsThisPlane = HitMap[plane];
	    

	    std::cout<<"Making centre of mass refit" << std::endl;
	    // Adjust centre of mass
	    //-----------------------

	    // Get seed central point in wire, time coordinates

	    double SeedCentralWire=0, SeedCentralTime=0;
	    
	    double SeedPt[3], Err[3], SeedDir[3];
	    TheSeed->GetPoint(    SeedPt,   Err);
	    TheSeed->GetDirection(SeedDir,  Err);

	    TVector3 SeedPoint(SeedPt[0],      SeedPt[1],  SeedPt[2]);
	    TVector3 SeedDirection(SeedDir[0], SeedDir[1], SeedDir[2]);

	    int centralchannel = geom->NearestChannel(SeedPt, plane, 0, 0);
	   
	    unsigned int c,t,p,w ;
	    geom->ChannelToWire(centralchannel, c,t,p,w);
	    SeedCentralWire = w;
	    SeedCentralTime = det->ConvertXToTicks(SeedPt[0],plane,t,c);

	    
	    //Get the centre of mass of hits in this view
	    

	    double HitCentralTime=0, HitCentralWire=0, TotalCharge=0;
	    for(art::PtrVector<recob::Hit>::const_iterator itHit=HitsThisPlane.begin(); itHit!=HitsThisPlane.end(); itHit++)
	      {
		unsigned int Wire;
		geom->ChannelToWire((*itHit)->Channel(), c, t, p, Wire);		
		HitCentralTime  += (*itHit)->PeakTime() * (*itHit)->Charge();
		HitCentralWire  += Wire                 * (*itHit)->Charge();
		TotalCharge     += (*itHit)->Charge();
	      }
	    
	    HitCentralTime/=TotalCharge;
	    HitCentralWire/=TotalCharge;


	    // Move the new centre of mass half way between these two
	    double CentralWireShift = 0.5*(HitCentralWire - SeedCentralWire);
	    double CentralTimeShift = 0.5*(HitCentralTime - SeedCentralTime);
	    
	    std::cout<<"SeedFinder suggests moving COM from " << SeedCentralWire <<", " << SeedCentralTime;
	    std::cout<<" to " << HitCentralWire<<" " << HitCentralTime<<std::endl;		
	    
	    // Find the direction in which we are allowed to shift the central point (perp to wires)
	    double Wire1End1[3],Wire1End2[3];

	    geom->WireEndPoints(0,0,plane,1, Wire1End1, Wire1End2);
	    
	    TVector3 WireVec(Wire1End2[0]-Wire1End1[0],
			     Wire1End2[1]-Wire1End1[1],
			     Wire1End2[2]-Wire1End1[2]);
	    
	    WireVec=WireVec.Unit();

	    TVector3 XVec(1,0,0);
	    
	    TVector3 PlaneNormDirection=(WireVec.Cross(XVec)).Unit();
	    
	    double WirePitch = geom->WirePitch();
	   

	    // Move seed point

	    SeedPt[0] = det->ConvertTicksToX(SeedCentralTime + CentralTimeShift,plane,0,0);
	    SeedPt[1] = SeedPoint[1] + CentralWireShift * WirePitch * PlaneNormDirection[1] ;
	    SeedPt[2] = SeedPoint[2] + CentralWireShift * WirePitch * PlaneNormDirection[2] ;

	    TheSeed->SetPoint(SeedPt);

	    for(int i=0; i!=3; i++)
	      SeedPoint[i]=SeedPt[i];
	    
	    SeedCentralWire += CentralWireShift;
	    SeedCentralTime += CentralTimeShift;

	    

	    // Adjust direction
	    //-----------------------
	    

	    std::cout<<"Making direction refit" << std::endl;

	    double OneTickInX = det->GetXTicksCoefficient();
	    
	    double d2=0, d=0, N=0, dx=0, x=0;


	    double Theta=0;
	    TotalCharge=0;
	    for(art::PtrVector<recob::Hit>::const_iterator itHit=HitsThisPlane.begin(); itHit!=HitsThisPlane.end(); itHit++)
	      {
		unsigned int Wire;
		double Time;
		double XDisp;
		double DDisp;

		geom->ChannelToWire((*itHit)->Channel(), c, t, p, Wire);
		DDisp = (double(Wire)-SeedCentralWire)*WirePitch;
		
		Time=(*itHit)->PeakTime();
		XDisp = Time*OneTickInX;
		
		x  += XDisp;
		d  += DDisp;
		dx += DDisp*XDisp;
		N  += 1;
		d2 += DDisp*DDisp;

		TotalCharge+=(*itHit)->Charge();
	
 		Theta += atan( (Time-SeedCentralTime)*OneTickInX
			       /  ((double(Wire)-SeedCentralWire) * WirePitch)) * (*itHit)->Charge();
	      }
	    Theta/=TotalCharge;

	    TMatrixD LeastSquaresLHS(2,2);
	    LeastSquaresLHS[0][0]=d2;
	    LeastSquaresLHS[0][1]=d;
	    LeastSquaresLHS[1][0]=d;
	    LeastSquaresLHS[1][1]=N;

	    TVectorD LeastSquaresRHS(2);
	    LeastSquaresRHS[0]=dx;
	    LeastSquaresRHS[1]=x;

	    LeastSquaresLHS.Invert();

	    TVectorD LeastSquaresResult(2);
	    LeastSquaresResult = LeastSquaresLHS*LeastSquaresRHS;

	    // parameters for t = ac + b
	    double a = LeastSquaresResult[0];
	    double b = LeastSquaresResult[1];


	    // Project seed direction vector into parts parallel and perp to pitch

	    double SeedPlaneComp     = SeedDirection.Dot(PlaneNormDirection.Unit());
	    double SeedTimeComp      = SeedDirection.Dot(XVec.Unit());
	    double SeedOutOfPlaneComp= SeedDirection.Dot(WireVec.Unit());

	    double SeedLengthInPlane = pow( pow(SeedPlaneComp,2)+pow(SeedTimeComp,2), 0.5);
	    
	    std::cout<<"Previous  a in plane : " << SeedTimeComp / SeedPlaneComp<<std::endl;
	    std::cout<<"Suggested a in plane : " << a <<std::endl;


	    double Newa=0.5*(a+SeedTimeComp/SeedPlaneComp);
	    
	    double OriginalSeedLength = SeedDirection.Mag();
	    
	    // Set seed direction to this theta without changing length
	    SeedTimeComp  = Newa*SeedPlaneComp;
	    
	    // build the new direction from these 3 orthogonal components
	    SeedDirection = 
	      SeedTimeComp       * XVec +
	      SeedPlaneComp      * PlaneNormDirection + 
	      SeedOutOfPlaneComp * WireVec.Unit(); 
	    
	    SeedDirection.SetMag(OriginalSeedLength);

	    // Set the seed direction accordingly
	    for(int i=0; i!=3; i++)
	      SeedDir[i]=SeedDirection[i];
	    
	    TheSeed->SetDirection(SeedDir);
	    
	  } // next plane
      } // next iteration
    
	
  }





  void SeedFinder::endJob()
  {

  }
}
