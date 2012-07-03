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
    fhicl::ParameterSet seedConfig = pset.get<fhicl::ParameterSet>("SeedConfig");

    fSptalg                = seedConfig.get<fhicl::ParameterSet>("SpacePointAlg");

    fInitSeedLength        = seedConfig.get<double>("InitSeedLength");
    fMinPointsInSeed       = seedConfig.get<unsigned int>("MinPointsInSeed");
    fAngularDev            = seedConfig.get<double>("AngularDev");

    fRefits                = seedConfig.get<double>("Refits");

    fExtendThresh          = seedConfig.get<double>("ExtendThresh");
    fExtendStep            = seedConfig.get<double>("ExtendStep");
    fExtendResolution      = seedConfig.get<double>("ExtendResolution");

    fMaxViewRMS.resize(3);
    fMaxViewRMS            = seedConfig.get<std::vector<double> >("MaxViewRMS"); 
  }

 


  

  //------------------------------------------------------------
  // Use the spacepoint service to turn a vector of hits 
  //  into a vector of spacepoints
  //
  
  std::vector<recob::SpacePoint>  SeedFinderAlgorithm::GetSpacePointsFromHitVector(art::PtrVector<recob::Hit>  Hits)
  {
    std::vector<recob::SpacePoint> ReturnVec;
    
    fSptalg.makeSpacePoints(Hits, ReturnVec,
    			    fSptalg.filter(), fSptalg.merge(), 0., 0.);

 
    return ReturnVec;
      
  }


  //------------------------------------------------------------
  // Given a set of spacepoints, find seeds, and catalogue
  //  spacepoints by the seeds they formed
  //
  std::vector<recob::Seed *> SeedFinderAlgorithm::FindSeeds(std::vector<recob::SpacePoint> const& AllSpacePoints, std::vector<std::vector<recob::SpacePoint> > & PointsInSeeds)
  {
    // Vector of seeds found to return
    std::vector<recob::Seed*>       ReturnVector;

    // This vector keeps track of the status of each point.  
    // The key is the position in the AllSpacePoints vector.
    // The value is 0: point unused, 1: point used in seed, 2: point thrown but unused 
    std::map<int,int>               PointStatus;

    // Keep track of how many SPs we used already 
    int TotalSPsUsed=0;
    
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
	    ReturnVector.push_back(TheSeed);
	    
	    std::vector<recob::SpacePoint> SPs;
	    for(size_t i=0; i!=PointsUsed.size(); ++i) 
	      SPs.push_back(AllSpacePoints.at(PointsUsed.at(i)));
	    PointsInSeeds.push_back(SPs);
	  }
	
	// Update the status of the spacepoints we used in this attempt
	if(TheSeed->IsValid())
	  for(size_t i=0; i!=PointsUsed.size(); PointsUsed.at(i++)=1);
	else
	  for(size_t i=0; i!=PointsUsed.size(); PointsUsed.at(i++)=2);

	// Keep track of how many we used
	TotalSPsUsed+=PointsUsed.size();

	// We run out of points when there are less than the minimum seed
	//  count unused. Then we break the loop.
	if( (AllSpacePoints.size() - TotalSPsUsed) < fMinPointsInSeed)
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
    else return new recob::Seed();
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









  bool SeedFinderAlgorithm::ExtendSeed(recob::Seed * TheSeed, std::vector<recob::SpacePoint> const& AllSpacePoints, std::map<int,int>& PointStatus, std::vector<int>& PointsUsed)
  {
    std::cout<<"Begin ExtendSeed method"<<std::endl;
    std::vector<int>                 PointsIncluded = PointsUsed;
    std::vector<recob::SpacePoint>   SPsUsed        = ExtractSpacePoints(AllSpacePoints, PointsUsed);

    recob::Seed * BestSeed = new recob::Seed(*TheSeed);

    double ThisDir[3], ThisPt[3], ThisErr[3];
    BestSeed->GetDirection( ThisDir, ThisErr );
    BestSeed->GetPoint(     ThisPt,  ThisErr );

    TVector3 VecDir( ThisDir[0], ThisDir[1], ThisDir[2] );
    TVector3 VecPt(  ThisPt[0],  ThisPt[1],  ThisPt[2]  );

    std::vector<double>              BestRMS        = GetHitRMS(TheSeed, SPsUsed);
    std::vector<double>              ThisRMS        = BestRMS;

    double BestAveRMS =  pow(pow(ThisRMS.at(0),2)+pow(ThisRMS.at(1),2)+pow(ThisRMS.at(2),2),0.5);
    double ThisAveRMS =  BestAveRMS;

    double BestdNdx   = double(SPsUsed.size()) / VecDir.Mag();
    double ThisdNdx   = BestdNdx;
    int BestN=SPsUsed.size();
    int ThisN=BestN;

    bool KeepExtending=true;
    // First extend backward


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


        // Find nearby spacepoints and for refit
	std::vector<int> NearbySPs =  DetermineNearbySPs(TheNewSeed, AllSpacePoints, PointStatus, fExtendResolution);
	std::vector<recob::SpacePoint> ThePoints = ExtractSpacePoints(AllSpacePoints, NearbySPs);

	RefitSeed(TheNewSeed,ThePoints);

        // With refitted seed, count number of hits and work out dNdx and RMS
        NearbySPs =  DetermineNearbySPs(TheNewSeed, AllSpacePoints, PointStatus, fExtendResolution);
        ThePoints = ExtractSpacePoints(AllSpacePoints, NearbySPs);

        int NoOfHits = CountHits(ThePoints);
        ThisN      = NoOfHits;
        ThisdNdx   = double(NoOfHits) / VecDir.Mag();
        ThisRMS    = GetHitRMS(TheNewSeed,ThePoints);
        ThisAveRMS = pow(pow(ThisRMS.at(0),2)+pow(ThisRMS.at(1),2)+pow(ThisRMS.at(2),2),0.5);

	std::cout<<" Deciding whether to increase: " << ThisdNdx<<" " <<BestdNdx<<", " << ThisAveRMS<<" " <<BestAveRMS<<std::endl;

        // Decide whether to keep the extended seed
        if((ThisdNdx > BestdNdx*fExtendThresh)&&(ThisAveRMS < BestAveRMS/fExtendThresh)&&(ThisN>BestN))
          {
            // if both then we keep the extended seed
            double SeedPos[3], Err[3];
            BestSeed->GetPoint(SeedPos,Err);

            BestAveRMS = ThisAveRMS;
	    BestdNdx   = ThisdNdx;
            BestN      = ThisN;

            delete BestSeed;
	    std::cout<<"Moving seed from " <<SeedPos[0]<<" " <<SeedPos[1]<< " " <<SeedPos[2]<<std::endl;
            BestSeed   = TheNewSeed;
            BestSeed->GetPoint(SeedPos,Err);
	    std::cout<<"to " <<SeedPos[0]<<" " <<SeedPos[1]<< " " <<SeedPos[2]<<std::endl;
          }
	else
          KeepExtending=false;

      }

    std::vector<int> FinalSPs;

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
	std::cout<<"   Before Refit: " << std::endl;
        TheNewSeed->Print();
        RefitSeed(TheNewSeed,ThePoints);
	std::cout<<"   After Refit: " << std::endl;
	TheNewSeed->Print();
        // With refitted seed, count number of hits and work out dNdx and RMS
	NearbySPs =  DetermineNearbySPs(TheNewSeed, AllSpacePoints, PointStatus, fExtendResolution);
        ThePoints = ExtractSpacePoints(AllSpacePoints, NearbySPs);

        int NoOfHits = CountHits(ThePoints);
        ThisN      = NoOfHits;
        ThisdNdx   = double(NoOfHits) / VecDir.Mag();
        ThisRMS    = GetHitRMS(TheNewSeed,ThePoints);
        ThisAveRMS = pow(pow(ThisRMS.at(0),2)+pow(ThisRMS.at(1),2)+pow(ThisRMS.at(2),2),0.5);

	std::cout<<" Deciding whether to increase: " << ThisdNdx<<" " <<BestdNdx<<std::endl;

        // Decide whether to keep the extended seed
        if((ThisdNdx > BestdNdx*fExtendThresh)&&(ThisAveRMS < BestAveRMS/fExtendThresh)&&(ThisN>BestN))
          {
            // if both then we keep the extended seed
            BestRMS=ThisRMS;
            BestN = ThisN;
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


    BestSeed->GetDirection(ThisDir, ThisErr);
    BestSeed->GetPoint(ThisPt,ThisErr);

    TheSeed->SetDirection(ThisDir,ThisErr);
    TheSeed->SetPoint(ThisPt,ThisErr);

    std::cout<<"Best seed at end of extend:" <<std::endl;
    TheSeed->Print();

    bool ReturnVal=false;

    std::cout<<"RMS Vals : ";
    for(size_t i=0; i!=BestRMS.size(); ++i)
      {
	std::cout << BestRMS.at(i) << " " << fMaxViewRMS.at(i)<<", " ;
        //      if(BestRMS.at(i)>fMaxViewRMS.at(i)) ReturnVal=true;
      }
    std::cout<<std::endl;
    return ReturnVal;
  }



  std::vector<int> SeedFinderAlgorithm::DetermineNearbySPs(recob::Seed* TheSeed, std::vector<recob::SpacePoint> AllSpacePoints, std::map<int, int> PointStatus, double ExtendResolution)
  {
    std::vector<int> ReturnVector;

    double SeedPt[3], SeedDir[3], Err[3];
    TheSeed->GetPoint( SeedPt,  Err );
    TheSeed->GetPoint( SeedDir, Err );

    double SeedLength = TheSeed->GetLength();

    for(int i=AllSpacePoints.size()-1; i!=-1; --i)
      {
	if(PointStatus[i]!=1)
	  {
	    // std::cout<<"NearbySPs Routine " << TheSeed->GetDistanceFrom(AllSpacePoints.at(i))<<std::endl;
	    if(( SeedPt[2]-AllSpacePoints.at(i).XYZ()[2] ) > SeedLength) break;
	    if( TheSeed->GetDistanceFrom(AllSpacePoints.at(i))
		< ExtendResolution )
	      ReturnVector.push_back(i);

	  }
      }

    return ReturnVector;
  }





  std::vector<double> SeedFinderAlgorithm::GetHitRMS(recob::Seed * TheSeed, std::vector<recob::SpacePoint> SpacePoints)
  {
    std::vector<double> ReturnVec;
    art::ServiceHandle<geo::Geometry>            geom;
    art::ServiceHandle<util::DetectorProperties>            det;


    std::map<int,art::PtrVector<recob::Hit> > HitMap;
    std::map<unsigned short, bool>            HitsClaimed;
    size_t Planes = geom->TPC(0).Nplanes();

    // for each spacepoint
    for(std::vector<recob::SpacePoint>::const_iterator itSP=SpacePoints.begin();
        itSP!=SpacePoints.end(); itSP++)
      {
        // get hits from each plane (ensuring each used once only)
	for(size_t plane=0; plane!=Planes; plane++)
	  {
	    art::PtrVector<recob::Hit> HitsThisSP = fSptalg.getAssociatedHits(*itSP);
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
    for(size_t plane=0; plane!=Planes; plane++)
      {
	art::PtrVector<recob::Hit> HitsThisPlane = HitMap[plane];
	double SeedCentralWire=0, SeedCentralTime=0;

	double SeedPt[3], Err[3], SeedDir[3];
	TheSeed->GetPoint(    SeedPt,   Err);
	TheSeed->GetDirection(SeedDir,  Err);

	TVector3 SeedPoint(SeedPt[0],      SeedPt[1],  SeedPt[2]);
        TVector3 SeedDirection(SeedDir[0], SeedDir[1], SeedDir[2]);

	double Wire1End1[3],Wire1End2[3];
        geom->WireEndPoints(0,0,plane,1, Wire1End1, Wire1End2);
	TVector3 WireVec(Wire1End2[0]-Wire1End1[0],
			 Wire1End2[1]-Wire1End1[1],
                         Wire1End2[2]-Wire1End1[2]);

        WireVec=WireVec.Unit();

        TVector3 XVec(1,0,0);

        TVector3 PlaneNormDirection=(WireVec.Cross(XVec)).Unit();

        double SeedDirInPlane = SeedDirection.Dot(PlaneNormDirection);
        double SeedDirInTime  = SeedDirection.Dot(XVec);


        int centralchannel = geom->NearestChannel(SeedPt, plane, 0, 0);

        unsigned int c,t,p,w ;
        geom->ChannelToWire(centralchannel, c,t,p,w);
        SeedCentralWire = w;
        SeedCentralTime = det->ConvertXToTicks(SeedPt[0],plane,t,c);

        double WirePitch = geom->WirePitch();
        double OneTickInX = det->GetXTicksCoefficient();

        double LSTot=0;
	int    Count=0;

	double a = SeedDirInTime / SeedDirInPlane;



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


  



  //----------------------------------------------------------------------------
  void SeedFinderAlgorithm::RefitSeed(recob::Seed * TheSeed, std::vector<recob::SpacePoint> SpacePoints)
  {

    std::cout<<"Refit module called on vector of " << SpacePoints.size() << " space points " << std::endl;

    // Get the services we need
    art::ServiceHandle<geo::Geometry>            geom;
    art::ServiceHandle<util::DetectorProperties> det;

    if((geom->NTPC()!=1)||geom->NTPC()!=1)
      {
	throw cet::exception("SeedFinderAlgorithm : Refit only works for 1 tpc, 1 cryostat detector")<<
	  "seed refitting feature not yet developped for multi cryostat or "<<
	  "multi TPC detector - see TrackFinder/SeedFinderAlgorithm.cxx"<<std::endl;	
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
	    unsigned int p = 0;
	    unsigned int w = 0;
	    unsigned int t = 0;
	    unsigned int c = 0;
	    
	    art::PtrVector<recob::Hit> HitsThisSP = fSptalg.getAssociatedHits(*itSP);
	    for(art::PtrVector<recob::Hit>::const_iterator itHit=HitsThisSP.begin();
		itHit!=HitsThisSP.end(); itHit++)
	      {
		geom->ChannelToWire((*itHit)->Channel(), c, t, p, w);
		if(p != plane) continue;

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
	// std::cout<<"SeedFinderAlgorithm running refit " << loop<<std::endl;
	for(size_t plane=0; plane!=Planes; plane++)
	  {
	    art::PtrVector<recob::Hit> HitsThisPlane = HitMap[plane];
	    

	    //   std::cout<<"Making centre of mass refit" << std::endl;
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
	    
	    //    std::cout<<"SeedFinderAlgorithm suggests moving COM from " << SeedCentralWire <<", " << SeedCentralTime;
	    //    std::cout<<" to " << HitCentralWire<<" " << HitCentralTime<<std::endl;		
	    
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
	    

	    //	    std::cout<<"Making direction refit" << std::endl;

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
	    //double b = LeastSquaresResult[1];


	    // Project seed direction vector into parts parallel and perp to pitch

	    double SeedPlaneComp     = SeedDirection.Dot(PlaneNormDirection.Unit());
	    double SeedTimeComp      = SeedDirection.Dot(XVec.Unit());
	    double SeedOutOfPlaneComp= SeedDirection.Dot(WireVec.Unit());

	    //double SeedLengthInPlane = pow( pow(SeedPlaneComp,2)+pow(SeedTimeComp,2), 0.5);
	    
	    //	    std::cout<<"Previous  a in plane : " << SeedTimeComp / SeedPlaneComp<<std::endl;
	    //   std::cout<<"Suggested a in plane : " << a <<std::endl;


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
