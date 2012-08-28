//
// Name: SeedFinderAlgorithm.cxx
//
// Purpose: Implementation file for module SeedFinderAlgorithm.
//
// Ben Jones, MIT
//

#include <vector>
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "TrackFinder/SeedFinderAlgorithm.h"
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

  //----------------------------------------------------------------------------
  SeedFinderAlgorithm::SeedFinderAlgorithm(const fhicl::ParameterSet& pset) :
    fSptalg(pset.get<fhicl::ParameterSet>("SpacePointAlg"))
  {
    reconfigure(pset);
  }

  //----------------------------------------------------------------------------
  SeedFinderAlgorithm::~SeedFinderAlgorithm()
  {
  }

  //----------------------------------------------------------------------------
  void SeedFinderAlgorithm::reconfigure(fhicl::ParameterSet const& pset)
  {
 
    fSptalg                = pset.get<fhicl::ParameterSet>("SpacePointAlg");

    fInitSeedLength        = pset.get<double>("InitSeedLength");
    fMinPointsInSeed       = pset.get<unsigned int>("MinPointsInSeed");
    fAngularDev            = pset.get<double>("AngularDev");

    fRefits                = pset.get<double>("Refits");

    fExtendThresh          = pset.get<double>("ExtendThresh");
    fExtendStep            = pset.get<double>("ExtendStep");
    fExtendResolution      = pset.get<double>("ExtendResolution");

    fMaxViewRMS.resize(3);
    fMaxViewRMS            = pset.get<std::vector<double> >("MaxViewRMS"); 
  }




  //------------------------------------------------------------
  // Find seeds in a collection of collections of SPs
  //   (for cluster method)

  std::vector<recob::Seed> SeedFinderAlgorithm::FindSeeds(std::vector<std::vector<recob::SpacePoint> > const& InputPoints, std::vector<std::vector<recob::SpacePoint> >& CataloguedSPs)
  {
    CataloguedSPs.clear();
    std::vector<recob::Seed> ReturnVector;
    for(size_t i=0; i!=InputPoints.size(); ++i)
      {
	std::vector<std::vector<recob::SpacePoint> > SPsThisCombo;
	std::vector<recob::Seed> SeedsThisCombo = FindSeeds(InputPoints.at(i), SPsThisCombo);
	for(size_t j=0; j!=SeedsThisCombo.size(); ++j)
	  {
	    ReturnVector.push_back(  SeedsThisCombo[j] ) ;
	    CataloguedSPs.push_back( SPsThisCombo[j]   ) ;
	  }
      }

    return ReturnVector;
  }

  

  //------------------------------------------------------------
  // Use the spacepoint service to turn a vector of hits 
  //  into a vector of spacepoints
  //
  
  std::vector<recob::SpacePoint>  SeedFinderAlgorithm::GetSpacePointsFromHitVector(art::PtrVector<recob::Hit>  Hits)
  {
    std::vector<recob::SpacePoint> ReturnVec;
    
    fSptalg.makeSpacePoints(Hits, ReturnVec);
    //   std::cout<<"Making SPs : " <<Hits.size()<< " " << ReturnVec.size()<< std::endl;
 
    return ReturnVec;
      
  }


  

  //------------------------------------------------------------
  // Given a set of spacepoints, find seeds, and catalogue
  //  spacepoints by the seeds they formed
  //
  std::vector<recob::Seed> SeedFinderAlgorithm::FindSeeds(std::vector<recob::SpacePoint> const& AllSpacePoints, std::vector<std::vector<recob::SpacePoint> > & PointsInSeeds)
  {
    // Vector of seeds found to return
    std::vector<recob::Seed>       ReturnVector;

    // This vector keeps track of the status of each point.  
    // The key is the position in the AllSpacePoints vector.
    // The value is 0: point unused, 1: point used in seed, 2: point thrown but unused 
    std::map<int,int>               PointStatus;
    
    
    // Keep track of how many SPs we used already 
    int TotalSPsUsed=0;
    int TotalNoOfSPs=int(AllSpacePoints.size());
    
    // Empty the relevant vectors
    ReturnVector.clear();
    PointStatus.clear();
    PointsInSeeds.clear();

    // Follow this loop until all SPs used up
    bool KeepChopping=true;

    while(KeepChopping)
      {
	// This vector keeps a list of the points used in this seed
	std::vector<int> PointsUsed;

	// Find exactly one seed, starting at high Z
	recob::Seed* TheSeed = FindSeedAtEnd(AllSpacePoints, PointStatus, PointsUsed);

	// If it was a good seed, collect up the relevant spacepoints
	// and add the seed to the return vector 
	if(TheSeed->IsValid())
	  {
	    ReturnVector.push_back(*TheSeed);
	    
	    std::vector<recob::SpacePoint> SPs;
	    for(size_t i=0; i!=PointsUsed.size(); ++i) 
	      SPs.push_back(AllSpacePoints.at(PointsUsed.at(i)));
	    PointsInSeeds.push_back(SPs);
	  }
	// Update the status of the spacepoints we used in this attempt
	if(TheSeed->IsValid())
	  for(size_t i=0; i!=PointsUsed.size(); PointStatus[PointsUsed.at(i++)]=1);
	else
	  for(size_t i=0; i!=PointsUsed.size(); PointStatus[PointsUsed.at(i++)]=2);
	
	TotalSPsUsed=0;
	for(size_t i=0; i!=PointStatus.size(); ++i)
	  {
	    if(PointStatus[i]!=0) TotalSPsUsed++;
	  }

	if((TotalNoOfSPs-TotalSPsUsed)<fMinPointsInSeed)
	  KeepChopping=false;
      }
    return ReturnVector;
  }


  //------------------------------------------------------------
  // Try to find one seed at the high Z end of a set of spacepoints 
  //

  recob::Seed * SeedFinderAlgorithm::FindSeedAtEnd(std::vector<recob::SpacePoint> const& Points, std::map<int,int>& PointStatus, std::vector<int>& PointsInRange)
  {
    // This pointer will be returned later
    recob::Seed * ReturnSeed;
    
    // Keep track of spacepoints we used, not just their IDs
    std::vector<recob::SpacePoint> PointsUsed;
   
    // Clear output vector
    PointsInRange.clear();
   
    // Loop through hits looking for highest Z seedable point
    TVector3 HighestZPoint;
    bool NoPointFound=true;
    int counter = Points.size()-1; 
    while(NoPointFound==true)
      {
	if(PointStatus[counter]==0)
	  {
	    HighestZPoint = TVector3(Points.at(counter).XYZ()[0],
				     Points.at(counter).XYZ()[1],
				     Points.at(counter).XYZ()[2]);
	    NoPointFound=false;
	  }
	else
	  counter--;
      }


    // Now we have the high Z point, loop through collecting
    // near enough hits.  We look 2 seed lengths away, since 
    // the seed is bidirectional from the central point 
    TVector3 CentreOfPoints(0.,0.,0.);
    double TwiceLength = 2.0*fInitSeedLength;
    
    for(int index=Points.size()-1; index!=-1; --index)
      {
	if(PointStatus[index]==0)
	  {
	    // first check z, then check total distance
	    //  (much faster, since most will be out of range in z anyway)
	    if( ( HighestZPoint[2] - Points.at(index).XYZ()[2] ) < TwiceLength)
	      {
		double DistanceToHighZ =pow( pow(HighestZPoint[0]-Points.at(index).XYZ()[0],2) +
					     pow(HighestZPoint[1]-Points.at(index).XYZ()[1],2) +
					     pow(HighestZPoint[2]-Points.at(index).XYZ()[2],2),0.5 ); 
		if( DistanceToHighZ < TwiceLength)
		  {
		    PointsInRange.push_back(index);
		    CentreOfPoints[0]+=Points.at(index).XYZ()[0];
		    CentreOfPoints[1]+=Points.at(index).XYZ()[1];
		    CentreOfPoints[2]+=Points.at(index).XYZ()[2];
		    PointsUsed.push_back(Points.at(index));
		    //  std::cout<<"In range : " << index<<std::endl;
		  }
	      }
	    else break;
	  }
      }
    
    
    // Check we have enough points in here to form a seed,
    // otherwise return a dud
    unsigned int NPoints = PointsInRange.size();
    if(NPoints<fMinPointsInSeed) return new recob::Seed();
    CentreOfPoints = (1./float(NPoints)) * CentreOfPoints;
    
    std::cout<<"Trying seed at " << CentreOfPoints[0]<<" " <<
      CentreOfPoints[1]<<" " << CentreOfPoints[2]<<" " <<NPoints<<std::endl;
    
    // See if seed points have some linearity
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
    
    // If seed is tight enough angularly    
    if(AngularDev<fAngularDev)
      {
	double PtArray[3], DirArray[3];
	
	PtArray[0]=CentreOfPoints.X();
	PtArray[1]=CentreOfPoints.Y();
	PtArray[2]=CentreOfPoints.Z();
	

	// calculate direction from average over hits
	
	TVector3 SeedDirection;
	
	SeedDirection.SetMagThetaPhi(fInitSeedLength, meantheta, meanphi);
	

	DirArray[0]=SeedDirection.X();
	DirArray[1]=SeedDirection.Y();
	DirArray[2]=SeedDirection.Z();
	

	ReturnSeed = new recob::Seed(PtArray,DirArray);
    }
    else
      {
	return new recob::Seed();
      }
	
    bool ThrowOutSeed = false;
   
      // If we use extendable seeds, go through extend and refit procedure
    if(fExtendThresh>0)
      {
	ThrowOutSeed = ExtendSeed(ReturnSeed, Points, PointStatus, PointsInRange);
      }

    // Otherwise, make optional refit of fixed length seed
    else
      {	
	if(fRefits>0)
	  RefitSeed(ReturnSeed, PointsUsed);
      }
	
    if(fMaxViewRMS.at(0)>0)
      {
	std::vector<double> RMS = GetHitRMS(ReturnSeed, PointsUsed);
	std::cout<<"RMS vector size : " << RMS.size() <<", contents ";
	for(size_t j=0; j!=fMaxViewRMS.size(); j++)
	  {
	    std::cout<<RMS.at(j)<<" ";
	    if(fMaxViewRMS.at(j)<RMS.at(j)) ThrowOutSeed=true;
	  }
	std::cout<<std::endl; 
      }


    // If the seed is marked as bad, return a dud, otherwise
    //  return the ReturnSeed pointer
    if(!ThrowOutSeed)
      return ReturnSeed;
    else
      return new recob::Seed;
  }






  //-----------------------------------------------------
  // Try walking a seed out past its end to find more hits.
  // After every step, refit to find optimal centre and direction. 
  //

  bool SeedFinderAlgorithm::ExtendSeed(recob::Seed * TheSeed, std::vector<recob::SpacePoint> const& AllSpacePoints, std::map<int,int>& PointStatus, std::vector<int>& PointsUsed)
  {
    std::cout<<"Begin ExtendSeed method"<<std::endl;
   
    // Get the actual spacepoints for the IDs provided
    std::vector<recob::SpacePoint>   SPsUsed        = ExtractSpacePoints(AllSpacePoints, PointsUsed);

    // This is the seed we will return - initially make a fresh seed with the same coordinates as the input seed
    recob::Seed * BestSeed = new recob::Seed(*TheSeed);

    // Get direction and centre information for this seed
    double ThisDir[3], ThisPt[3], ThisErr[3];
    BestSeed->GetDirection( ThisDir, ThisErr );
    BestSeed->GetPoint(     ThisPt,  ThisErr );

    TVector3 VecDir( ThisDir[0], ThisDir[1], ThisDir[2] );
    TVector3 VecPt(  ThisPt[0],  ThisPt[1],  ThisPt[2]  );

    // Keep a record of some quantities for the best seed found
    //  so far.  We are only allowed to extend if we meet a set
    //  of seed quality criteria relative to the best seed so far found.
    std::vector<double>              BestRMS        = GetHitRMS(TheSeed, SPsUsed);
    std::vector<double>              ThisRMS        = BestRMS;
    double BestAveRMS =  pow(pow(ThisRMS.at(0),2)+pow(ThisRMS.at(1),2)+pow(ThisRMS.at(2),2),0.5);
    double ThisAveRMS =  BestAveRMS;

    int NoOfHits = CountHits(SPsUsed);

    int ThisN         = NoOfHits;
    int BestN         = ThisN;

    double ThisdNdx   = double(NoOfHits) / VecDir.Mag();
    double BestdNdx   = ThisdNdx;


    // We extend the seed in both directions.  Backward first:

    bool KeepExtending=true;
    while(KeepExtending!=false)
      {
        // Get data from existing seed
        BestSeed->GetDirection( ThisDir, ThisErr );
	BestSeed->GetPoint(     ThisPt,  ThisErr );

        VecDir = TVector3( ThisDir[0], ThisDir[1], ThisDir[2] );
	VecPt  = TVector3(  ThisPt[0],  ThisPt[1],  ThisPt[2]  );

        // Make seed extension
        VecDir.SetMag(VecDir.Mag() + fExtendStep);
	VecPt = VecPt + fExtendStep * VecDir.Unit();
        for(int i=0; i!=3; ++i)
          {
            ThisDir[i] = VecDir[i];
            ThisPt[i]  = VecPt[i];
          }

	recob::Seed* TheNewSeed = new recob::Seed(ThisPt, ThisDir, ThisErr, ThisErr);


        // Find nearby spacepoints and refit
	std::vector<int> NearbySPs               = DetermineNearbySPs(TheNewSeed, AllSpacePoints, PointStatus, fExtendResolution);
	if(NearbySPs.size()<3) return true;
	std::cout<<"Size of SP vec in ext "<< NearbySPs.size()<<std::endl;
	std::vector<recob::SpacePoint> ThePoints = ExtractSpacePoints(AllSpacePoints, NearbySPs);

	RefitSeed(TheNewSeed,ThePoints);

        // With refitted seed, count number of hits and work out dNdx and RMS
        NearbySPs =  DetermineNearbySPs(TheNewSeed, AllSpacePoints, PointStatus, fExtendResolution);
	if(NearbySPs.size()<2) return true;
        ThePoints = ExtractSpacePoints(AllSpacePoints, NearbySPs);

        NoOfHits = CountHits(ThePoints);
	std::cout<<"SPs: " <<NearbySPs.size()<<" " << ThePoints.size()<< " Hits: " << NoOfHits<<std::endl;
        ThisN        = NoOfHits;
        ThisdNdx     = double(NoOfHits) / VecDir.Mag();
        ThisRMS      = GetHitRMS(TheNewSeed,ThePoints);
        ThisAveRMS   = pow(pow(ThisRMS.at(0),2)+pow(ThisRMS.at(1),2)+pow(ThisRMS.at(2),2),0.5);

	std::cout<<" Deciding whether to increase: " << ThisdNdx<<" " <<BestdNdx<<", " << ThisAveRMS<<" " <<BestAveRMS<<std::endl;

        // Decide whether to keep the extended seed
        if((ThisdNdx > BestdNdx*fExtendThresh)&&(ThisAveRMS < BestAveRMS/fExtendThresh)&&(ThisN>BestN))
          {
            BestAveRMS = ThisAveRMS;
	    BestdNdx   = ThisdNdx;
            BestN      = ThisN;

            delete BestSeed;
            BestSeed   = TheNewSeed;
          }
	else
          KeepExtending=false;

      }

    // Then extend forward
    KeepExtending=true;
    while(KeepExtending!=false)
      {
	// Get data from existing seed
        BestSeed->GetDirection( ThisDir, ThisErr );
        BestSeed->GetPoint(     ThisPt,  ThisErr );

        VecDir = TVector3( ThisDir[0], ThisDir[1], ThisDir[2] );
        VecPt  = TVector3(  ThisPt[0],  ThisPt[1],  ThisPt[2]  );


        // Make seed extension
	VecDir.SetMag(VecDir.Mag() + fExtendStep);
        VecPt = VecPt - fExtendStep * VecDir.Unit();
	for(int i=0; i!=3; ++i)
          {
            ThisDir[i] = VecDir[i];
            ThisPt[i]  = VecPt[i];
          }
	recob::Seed* TheNewSeed = new recob::Seed(ThisPt, ThisDir, ThisErr, ThisErr);

	// Find nearby spacepoints and for refit
	std::vector<int> NearbySPs =  DetermineNearbySPs(TheNewSeed, AllSpacePoints, PointStatus, fExtendResolution);
	std::vector<recob::SpacePoint> ThePoints = ExtractSpacePoints(AllSpacePoints, NearbySPs);
        RefitSeed(TheNewSeed,ThePoints);

        // With refitted seed, count number of hits and work out dNdx and RMS
	NearbySPs    = DetermineNearbySPs(TheNewSeed, AllSpacePoints, PointStatus, fExtendResolution);
        ThePoints    = ExtractSpacePoints(AllSpacePoints, NearbySPs);

        NoOfHits = CountHits(ThePoints);
        ThisN        = NoOfHits;
        ThisdNdx     = double(NoOfHits) / VecDir.Mag();
        ThisRMS      = GetHitRMS(TheNewSeed,ThePoints);
        ThisAveRMS   = pow(pow(ThisRMS.at(0),2)+pow(ThisRMS.at(1),2)+pow(ThisRMS.at(2),2),0.5);

	std::cout<<" Deciding whether to increase: " << ThisdNdx<<" " <<BestdNdx<<std::endl;

        // Decide whether to keep the extended seed
        if((ThisdNdx > BestdNdx*fExtendThresh)&&(ThisAveRMS < BestAveRMS/fExtendThresh)&&(ThisN>BestN))
          {
            BestRMS    = ThisRMS;
            BestN      = ThisN;
            BestAveRMS = ThisAveRMS;
            BestdNdx   = ThisdNdx;

            delete BestSeed;
            BestSeed   = TheNewSeed;
          }
        else
          {
	    PointsUsed =  DetermineNearbySPs(TheNewSeed, AllSpacePoints, PointStatus, fExtendResolution);
            for(size_t i=0; i!=NearbySPs.size(); ++i) PointStatus[NearbySPs.at(i)]=1;

            KeepExtending=false;
          }
      }
    std::cout<<"Done extending.  Final length : " << BestSeed->GetLength()<<std::endl;

    BestSeed->GetDirection( ThisDir, ThisErr);
    BestSeed->GetPoint(     ThisPt,  ThisErr);

    TheSeed->SetDirection(  ThisDir, ThisErr);
    TheSeed->SetPoint(      ThisPt,  ThisErr);

    bool ReturnVal=false;
    /*
    std::cout<<"RMS Vals after extend: ";
    for(size_t i=0; i!=BestRMS.size(); ++i)
      {
	std::cout << BestRMS.at(i) << " " << fMaxViewRMS.at(i)<<", " ;
        //      if(BestRMS.at(i)>fMaxViewRMS.at(i)) ReturnVal=true;
      }
    std::cout<<std::endl;
    */
    return ReturnVal;
  }


  //-----------------------------------------------------
  // Given a seed, find which spacepoints from a collection
  //  are within some distance
  //

  std::vector<int> SeedFinderAlgorithm::DetermineNearbySPs(recob::Seed* TheSeed, std::vector<recob::SpacePoint> AllSpacePoints, std::map<int, int> PointStatus, double MaxDistance)
  {
    std::vector<int> ReturnVector;

    double SeedPt[3], Err[3];
    TheSeed->GetPoint( SeedPt,  Err );
    
   
    
    double SeedLength = TheSeed->GetLength();

    for(int i=AllSpacePoints.size()-1; i!=-1; --i)
      {
	if(PointStatus[i]!=1)
	  {
	    // check Z separation first - helps prevent doing unnecessary calculations
	    //  std::cout<<"NearbySP Dist " << TheSeed->GetDistanceFrom(AllSpacePoints.at(i))<<std::endl;
	    if(( SeedPt[2]-AllSpacePoints.at(i).XYZ()[2] ) > SeedLength) 
	      break;
	    
 	    else if( TheSeed->GetDistanceFrom(AllSpacePoints.at(i))
		     < MaxDistance )
	      ReturnVector.push_back(i);
	    
	  }
	
      }

    return ReturnVector;
  }



  //-----------------------------------------------------
  // Find the RMS distance of the hits in a set of spacepoints
  //  from the spine of a seed in each view.
  //

  std::vector<double> SeedFinderAlgorithm::GetHitRMS(recob::Seed * TheSeed, std::vector<recob::SpacePoint> SpacePoints)
  {
    // Vector to return later, of one RMS value for each view
    std::vector<double>    ReturnVec;

    // Services we need
    art::ServiceHandle<geo::Geometry>                geom;
    art::ServiceHandle<util::DetectorProperties>     det;

    // At present, algorithm works only for 1 TPC, 1 cryostat
    unsigned int c=0, t=0;

    std::map<int,art::PtrVector<recob::Hit> > HitMap;
    std::map<unsigned short, bool>            HitsClaimed;
    size_t Planes = geom->TPC(0).Nplanes();

    // for each spacepoint
    for(std::vector<recob::SpacePoint>::const_iterator itSP=SpacePoints.begin();
        itSP!=SpacePoints.end(); itSP++)
      {
	art::PtrVector<recob::Hit> HitsThisSP = fSptalg.getAssociatedHits(*itSP);
	for(art::PtrVector<recob::Hit>::const_iterator itHit=HitsThisSP.begin();
	    itHit!=HitsThisSP.end(); itHit++)
	  {
	    if(!HitsClaimed[(*itHit)->Channel()])
	      {
		//		std::cout<<"RMS : hit in view " << (*itHit)->View()<<std::endl;
		HitMap[(*itHit)->View()].push_back(*itHit);
		HitsClaimed[(*itHit)->Channel()]=true;
	      }
	    
	  }
      }
   
    // Find the RMS in each plane
    for(size_t planeno=0; planeno!=Planes; planeno++)
      {
	int view = geom->TPC(0).Plane(planeno).View();
	
	art::PtrVector<recob::Hit> HitsThisPlane = HitMap[view];
	double SeedCentralWire=0, SeedCentralTime=0;

	double SeedPt[3], Err[3], SeedDir[3];
	TheSeed->GetPoint(      SeedPt,   Err);
	TheSeed->GetDirection(  SeedDir,  Err);

	TVector3 SeedPoint(     SeedPt[0],     SeedPt[1],   SeedPt[2]);
        TVector3 SeedDirection( SeedDir[0],    SeedDir[1],  SeedDir[2]);

	// Determine useful directions for this plane
	double Wire0End1[3],Wire0End2[3];
	double Wire1End1[3],Wire1End2[3];


	geom->WireEndPoints(c,t,planeno,0, Wire0End1, Wire0End2);
	geom->WireEndPoints(c,t,planeno,1, Wire1End1, Wire1End2);

	TVector3 Wire0End1Vec(Wire0End1[0],Wire0End1[1],Wire0End1[2]);
	TVector3 Wire0End2Vec(Wire0End2[0],Wire0End2[1],Wire0End2[2]);
	TVector3 Wire1End1Vec(Wire1End1[0],Wire1End1[1],Wire1End1[2]);
	TVector3 Wire1End2Vec(Wire1End2[0],Wire1End2[1],Wire1End2[2]);

	TVector3 WireVec    = (Wire0End1Vec-Wire0End2Vec).Unit();
	TVector3 XVec(1,0,0);
	TVector3 PlaneNormDirection=(WireVec.Cross(XVec)).Unit();
	if( (Wire1End1Vec-Wire0End1Vec).Dot(PlaneNormDirection)<0)
	  PlaneNormDirection = -PlaneNormDirection;

	double   WirePitch  = geom->WirePitch();
	double   WireOffset = PlaneNormDirection.Dot(Wire1End1);

	SeedCentralWire = (SeedPoint.Dot(PlaneNormDirection) - WireOffset) / WirePitch;

        double SeedDirInPlane = SeedDirection.Dot(PlaneNormDirection);
        double SeedDirInTime  = SeedDirection.Dot(XVec);

        SeedCentralTime = det->ConvertXToTicks(SeedPt[0],planeno,t,c);

        double OneTickInX = det->GetXTicksCoefficient();

        double LSTot=0;
	int    Count=0;

	double a = SeedDirInTime / SeedDirInPlane;


	unsigned int p=0;
	for(art::PtrVector<recob::Hit>::const_iterator itHit=HitsThisPlane.begin(); itHit!=HitsThisPlane.end(); itHit++)
          {
            unsigned int Wire;
            double Time;
            double XDisp;
            double DDisp;

            geom->ChannelToWire((*itHit)->Channel(), c, t, p, Wire);
            DDisp = (double(Wire)-SeedCentralWire)*WirePitch;

            Time=(*itHit)->PeakTime();
            XDisp = (Time-SeedCentralTime)*OneTickInX;

            LSTot += pow(XDisp - a * DDisp, 2);
            Count++;
          }
	ReturnVec.push_back(pow(LSTot/Count,0.5));
      }
    return ReturnVec;
  }


  //-----------------------------------------------------
  // Given a seed and some spacepoints, refit the seed
  //  to the hits in the SPs
  //

  void SeedFinderAlgorithm::RefitSeed(recob::Seed * TheSeed, std::vector<recob::SpacePoint> SpacePoints)
  {

    std::cout<<"Refit module called on vector of " << SpacePoints.size() << " space points " << std::endl;
    //  std::cout<<"Beginning of refit: "<<std::endl;
    //  TheSeed->Print();
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
	art::PtrVector<recob::Hit> HitsThisSP = fSptalg.getAssociatedHits(*itSP);
	for(art::PtrVector<recob::Hit>::const_iterator itHit=HitsThisSP.begin();
	    itHit!=HitsThisSP.end(); itHit++)
	  {
	    //  std::cout<<"Hit view " << (*itHit)->View()<<std::endl;
	    if(!HitsClaimed[(*itHit)->Channel()])
	      {
		HitMap[(*itHit)->View()].push_back(*itHit);
		HitsClaimed[(*itHit)->Channel()]=true;
	      }
	  }
      }

    double                  WirePitch  = geom->WirePitch();
    TVector3                XVec       = TVector3(1,0,0);

    std::vector<double>     WireOffsets;            WireOffsets.resize(Planes);
    std::vector<TVector3>   PlaneNormDirections;    PlaneNormDirections.resize(Planes);
    std::vector<TVector3>   WireVecs;               WireVecs.resize(Planes);
    std::vector<double>     HitCentralTimes;        HitCentralTimes.resize(Planes);
    std::vector<double>     HitCentralWires;        HitCentralWires.resize(Planes);
    std::vector<int>        NHits;                  NHits.resize(Planes);

    double AverageXWeighted;
    double AverageXDemocratic;
    int    TotalHits;

    // Calculate geometrical stuff for each plane
    for(size_t planeno=0; planeno!=Planes; planeno++)
      {
	// Determine useful directions for this plane
	int view = geom->TPC(0).Plane(planeno).View();

	double Wire0End1[3],Wire0End2[3];
	double Wire1End1[3],Wire1End2[3];
	geom->WireEndPoints(0,0,planeno,0, Wire0End1, Wire0End2);
	geom->WireEndPoints(0,0,planeno,1, Wire1End1, Wire1End2);
	
	TVector3 Wire0End1Vec(Wire0End1[0],Wire0End1[1],Wire0End1[2]);
	TVector3 Wire0End2Vec(Wire0End2[0],Wire0End2[1],Wire0End2[2]);
	TVector3 Wire1End1Vec(Wire1End1[0],Wire1End1[1],Wire1End1[2]);
	TVector3 Wire1End2Vec(Wire1End2[0],Wire1End2[1],Wire1End2[2]);
	
	WireVecs[view]            = (Wire0End1Vec-Wire0End2Vec).Unit();
        PlaneNormDirections[view] = (WireVecs[view].Cross(XVec)).Unit();
	if( (Wire1End1Vec-Wire0End1Vec).Dot(PlaneNormDirections[view])<0)
	  PlaneNormDirections[view] = -PlaneNormDirections[view];
	
	WireOffsets[view]         = PlaneNormDirections.at(view).Dot(Wire0End1Vec);
	
	art::PtrVector<recob::Hit> HitsThisPlane = HitMap[view]; 
	unsigned int c=0, t=0, p=0;
	for(art::PtrVector<recob::Hit>::const_iterator itHit=HitsThisPlane.begin(); itHit!=HitsThisPlane.end(); itHit++)
	  {
	    unsigned int Wire;
	    geom->ChannelToWire((*itHit)->Channel(), c, t, p, Wire);
	    HitCentralTimes[view]  += (*itHit)->PeakTime();
	    HitCentralWires[view]  += Wire;
	  }
	
	HitCentralTimes[view]/=HitsThisPlane.size();
	HitCentralWires[view]/=HitsThisPlane.size();
	NHits[view]=HitsThisPlane.size();
	
	AverageXWeighted      += det->ConvertTicksToX(HitCentralTimes[view],planeno,0,0) * NHits[view];
	AverageXDemocratic    += det->ConvertTicksToX(HitCentralTimes[view],planeno,0,0);
	TotalHits             += NHits[view];
      }
    
    // Times are easy - no coordinate degeneracy
    AverageXDemocratic /= Planes;
    AverageXWeighted   /= TotalHits;
    

    // Begin iterative refit procedure
    // Loop the prescribed number of times
    for(int loop=0; loop!=fRefits; loop++)
      {
        // std::cout<<"SeedFinder running refit " << loop<<std::endl;
        for(size_t view=0; view!=Planes; view++)
          {
	    // Adjust central seed point
	    // --------------------------

            // Get seed central point in this view



            double SeedPt[3], Err[3], SeedDir[3];
            TheSeed->GetPoint(    SeedPt,   Err);
            TheSeed->GetDirection(SeedDir,  Err);

            TVector3 SeedPoint(SeedPt[0],      SeedPt[1],  SeedPt[2]);
            TVector3 SeedDirection(SeedDir[0], SeedDir[1], SeedDir[2]);

	    //	    std::cout<<"Moving seed point from "<<SeedPt[0]<<" " << SeedPt[1]<<" " <<SeedPt[2]<<std::endl;

	    // Find what the seed central wire was before
            double SeedCentralWire = (PlaneNormDirections.at(view).Dot(SeedPoint) - WireOffsets.at(view)) / WirePitch;
            
	    //    std::cout<<"Central wires: " << HitCentralWires[view]<< " " << SeedCentralWire<<std::endl;

	    // Move it half way to that predicted from hits
	  
	    double InPlaneShift = 0.5*(HitCentralWires[view] - SeedCentralWire) * WirePitch;

	    SeedCentralWire = 0.5*(HitCentralWires[view] + SeedCentralWire);          



            SeedPt[0] = AverageXDemocratic;
            SeedPt[1] = SeedPt[1] + InPlaneShift * PlaneNormDirections.at(view)[1] ;
            SeedPt[2] = SeedPt[2] + InPlaneShift * PlaneNormDirections.at(view)[2] ;

	    // Update the seed
            TheSeed->SetPoint(SeedPt);

	    //  std::cout<<"Moving seed point to "<<SeedPt[0]<<" " << SeedPt[1]<<" " <<SeedPt[2]<<std::endl<<std::endl;

            for(int i=0; i!=3; i++)
              SeedPoint[i]=SeedPt[i];



            // Adjust direction
            //-----------------------

            double d2=0, d=0, N=0, dx=0, x=0;

	    art::PtrVector<recob::Hit> HitsThisPlane = HitMap[view]; 
	    unsigned int c=0, t=0, p=0;
	    // std::cout<<"No of hits : " << HitsThisPlane.size()<<std::endl;
            for(art::PtrVector<recob::Hit>::const_iterator itHit=HitsThisPlane.begin(); itHit!=HitsThisPlane.end(); itHit++)
              {
               unsigned int Wire;
	       double Time;
	       double XDisp;
	       
	       double DDisp;

	       geom->ChannelToWire((*itHit)->Channel(), c, t, p, Wire);
	       DDisp = (double(Wire)-SeedCentralWire)*WirePitch;

	       Time=(*itHit)->PeakTime();
	       XDisp = det->ConvertTicksToX(Time,p,t,c) - AverageXDemocratic;

	       x  += XDisp;
	       d  += DDisp;
	       dx += DDisp*XDisp;
	       N  += 1;
	       d2 += DDisp*DDisp;

              }
	    TMatrixD LeastSquaresLHS(2,2);
            LeastSquaresLHS[0][0]=d2;
	    LeastSquaresLHS[0][1]=d;
            LeastSquaresLHS[1][0]=d;
            LeastSquaresLHS[1][1]=N;

	    TVectorD LeastSquaresRHS(2);
	    LeastSquaresRHS[0]=dx;
            LeastSquaresRHS[1]=x;
	    try{
	      LeastSquaresLHS.Invert();
	    }
	    catch(...)
	      {
		return;
	      }
            TVectorD LeastSquaresResult(2);
	    LeastSquaresResult = LeastSquaresLHS*LeastSquaresRHS;
	    // parameters for t = ac + b
            double a = LeastSquaresResult[0];
            //   double b = LeastSquaresResult[1];


	    // Project seed direction vector into parts parallel and perp to pitch

            double SeedPlaneComp     = SeedDirection.Dot(PlaneNormDirections.at(view).Unit());
            double SeedTimeComp      = SeedDirection.Dot(XVec.Unit());
	    double SeedOutOfPlaneComp= SeedDirection.Dot(WireVecs[view].Unit());


            double Newa=0.5*(a+SeedTimeComp/SeedPlaneComp);

            double OriginalSeedLength = SeedDirection.Mag();

            // Set seed direction to this theta without changing length
            SeedPlaneComp  = SeedTimeComp / Newa;

            // build the new direction from these 3 orthogonal components
            SeedDirection =
              SeedTimeComp       * XVec +
              SeedPlaneComp      * PlaneNormDirections.at(view) +
              SeedOutOfPlaneComp * WireVecs.at(view);

            SeedDirection.SetMag(OriginalSeedLength);
            //  std::cout<<"Moving seed dir from "<<SeedDir[0]<< " " <<SeedDir[1]<< " " <<SeedDir[2]<<std::endl;
            // Set the seed direction accordingly
            for(int i=0; i!=3; i++)
              SeedDir[i]=SeedDirection[i];
            //      std::cout<<"to "<<SeedDir[0]<< " " <<SeedDir[1]<< " " <<SeedDir[2]<<std::endl;
            TheSeed->SetDirection(SeedDir);

          } // next plane
      } // next iteration

    //    std::cout<<"End of refit: "<<std::endl;
    //  TheSeed->Print();
  
  }






  //------------------------------------------------------------
  // Given a set of spacepoints, count how many unique hits
  //  are present

  double SeedFinderAlgorithm::CountHits(std::vector<recob::SpacePoint>  SpacePoints)
  {
    art::ServiceHandle<geo::Geometry>            geom;

    // A map of HitID to true/false whether it has been counted already
    std::map<unsigned short, bool>            HitsClaimed;
    HitsClaimed.clear();
    
    // For each spacepoint, check hit contents
    for(std::vector<recob::SpacePoint>::const_iterator itSP=SpacePoints.begin();
	itSP!=SpacePoints.end(); itSP++)
      {
	unsigned int p = 0, w=0, t=0, c=0;
	
	art::PtrVector<recob::Hit> HitsThisSP = fSptalg.getAssociatedHits((*itSP));
	for(art::PtrVector<recob::Hit>::const_iterator itHit=HitsThisSP.begin();
	    itHit!=HitsThisSP.end(); itHit++)
	  {
	    geom->ChannelToWire((*itHit)->Channel(), c, t, p, w);
	    HitsClaimed[(*itHit)->Channel()]=true;
	  }
      }
    return HitsClaimed.size();
  }

  std::vector<recob::SpacePoint> SeedFinderAlgorithm::ExtractSpacePoints(std::vector<recob::SpacePoint> AllPoints, std::vector<int> IDsToExtract)
  {
    std::vector<recob::SpacePoint> ReturnVec;
    for(size_t i=0; i!=IDsToExtract.size(); ++i)
      ReturnVec.push_back(AllPoints.at(IDsToExtract.at(i)));

    return ReturnVec;
  }




}
