#include "art/Persistency/Common/PtrVector.h"

#ifndef BEZIERTRACKEVBUILDER_H
#define BEZIERTRACKEVBUILDER_H

//
// Name: BezierTrackJoiner.h
//
// Purpose: Header file for module BezierTrackJoiner.  This modules makes
//          bezier tracks out of seed collections
//
// Configuration parameters.
//
// SeedModuleLabel;     // Cluster module label (e.g. "dbcluster").
// HitModuleLabel;      // Hit module label (e.g. "FFTHitFinder")
//
// Ben Jones, MIT
//

#include "art/Framework/Core/EDProducer.h"

namespace recob
{
  class Track;
  class Hit;
}

class TVector3;

namespace trkf {

  class BezierTrack;

  class BezierTrackJoiner : public art::EDProducer
  {
  public:
 
    // Constructors, destructor

    explicit BezierTrackJoiner(fhicl::ParameterSet const& pset);
    virtual ~BezierTrackJoiner();

    
    // Overrides.

    void reconfigure(fhicl::ParameterSet const& pset);
    void beginJob();
    void produce(art::Event& evt);
    void endJob();
    
    void GetImpact(TVector3 t1pt, TVector3 t1dir, TVector3 t2pt, TVector3 t2dir, double& ImpactParam, double& Dist1, double& Dist2);
    

  private:

    // Fcl Attributes.

    std::string fTrackModuleLabel;
     
    double fJoinThreshold;
    double fVertexAngle;
    double fExtrapDistance;
  };
}

#endif // SPACEPOINTFINDER_H



#include "art/Framework/Core/ModuleMacros.h" 


namespace trkf {
  DEFINE_ART_MODULE(BezierTrackJoiner);
}



//
// Name: BezierTrackJoiner.cxx
//
// Purpose: Implementation file for module BezierTrackJoiner.
//
// Ben Jones, MIT
//

#include <vector>
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "RecoBase/Hit.h"
#include "RecoBase/Seed.h"
#include "RecoBase/Track.h"
#include "RecoBase/Vertex.h"
#include "RecoObjects/BezierTrack.h"
#include "Utilities/AssociationUtil.h"



namespace trkf {

  //-------------------------------------------------------------------------
  BezierTrackJoiner::BezierTrackJoiner(const fhicl::ParameterSet& pset)
  {
    reconfigure(pset);
    produces< std::vector<recob::Track> >();
    produces< std::vector<recob::Vertex> >();
    produces< art::Assns<recob::Vertex, recob::Track> >();
  }

  //-------------------------------------------------------------------------
  BezierTrackJoiner::~BezierTrackJoiner()
  {
  }

  //-------------------------------------------------------------------------
  void BezierTrackJoiner::reconfigure(fhicl::ParameterSet const& pset)
  {
    fTrackModuleLabel   = pset.get<std::string>("TrackModuleLabel");
    fJoinThreshold      = pset.get<double>("JoinThreshold");
    fVertexAngle        = pset.get<double>("VertexAngle");
    fExtrapDistance     = pset.get<double>("ExtrapDistance");
  }

  //-------------------------------------------------------------------------
  void BezierTrackJoiner::beginJob()
  {}


  //-------------------------------------------------------------------------
  void BezierTrackJoiner::produce(art::Event& evt)
  {



    // Declare products to store
    
    std::unique_ptr< std::vector<recob::Track > > tracksout ( new std::vector<recob::Track>);
    std::unique_ptr< std::vector<recob::Vertex > > verticesout ( new std::vector<recob::Vertex>);
    std::unique_ptr< art::Assns<recob::Vertex, recob::Track > > tvassn( new art::Assns<recob::Vertex, recob::Track>);
    
        
    mf::LogVerbatim("BezierTrackJoiner") << "Getting tracks ";

    // Extract track objects from event

    art::Handle< std::vector<recob::Track> > trackh;
    evt.getByLabel(fTrackModuleLabel, trackh);

    art::PtrVector<recob::Track>   Tracks;
    for(unsigned int i=0; i < trackh->size(); ++i)
      {
	art::Ptr<recob::Track> track(trackh,i);
	Tracks.push_back(track);
      }

    


    mf::LogVerbatim("BezierTrackJoiner") << "Make bez tracks ";
    // Make bezier tracks
    std::vector<trkf::BezierTrack> BTracks;
    BTracks.clear();
    for(size_t i=0; i!=Tracks.size(); i++)
      BTracks.push_back(BezierTrack(*Tracks.at(i)));


    std::vector<trkf::BezierTrack> JoinedTracks = BTracks;
    
    std::vector<TVector3> TrackEnd1s;
    std::vector<TVector3> TrackEnd2s;
    std::vector<TVector3> TrackDir1s;
    std::vector<TVector3> TrackDir2s;


    for(size_t i=0; i!=BTracks.size(); ++i)
      {
	TrackEnd1s.push_back( BTracks.at(i).GetTrackPointV(0));
	TrackEnd2s.push_back( BTracks.at(i).GetTrackPointV(1));
	TrackDir1s.push_back( BTracks.at(i).GetTrackDirectionV(0));
	TrackDir2s.push_back( BTracks.at(i).GetTrackDirectionV(1));	
      }
 
    std::vector<std::map<int,bool> >     TracksMeetAtPoints;    
    std::vector<std::vector<TVector3> >  MeetPoints;
     
    // The connection[i][j] gives impact parameters for track ends
    //  The sign tells us which end of track j 
    for(size_t i=0; i!=BTracks.size(); ++i)
      {
	for(size_t j=0; j!=BTracks.size(); ++j)
	  {
	    if(i!=j)
	      {
		double impact, disti, distj;
		GetImpact(TrackEnd1s[i],TrackDir1s[i],TrackEnd1s[j],TrackDir1s[j], impact, disti, distj);

		if((impact < fJoinThreshold)
		   && (fabs(disti) < fExtrapDistance )
		   && (fabs(distj) < fExtrapDistance )
		   && (disti < 0 )
		   && (distj < 0 ))
		  {
		    bool ThisTrackAlreadyMet =false;
		    TVector3 ExtrapPoint1 = TrackEnd1s[i] + TrackDir1s[j].Unit()*disti;
		    TVector3 ExtrapPoint2 = TrackEnd1s[j] + TrackDir1s[j].Unit()*distj;
		    
		    for(size_t pt=0; pt!=TracksMeetAtPoints.size(); ++pt)
		      {
			if(TracksMeetAtPoints.at(pt)[-i]==true)
			  {
			    TracksMeetAtPoints.at(pt)[-j]=true;
			    MeetPoints.at(pt).push_back(ExtrapPoint1);
			    MeetPoints.at(pt).push_back(ExtrapPoint2);
			    ThisTrackAlreadyMet=true;
			    
			  }
			else if(TracksMeetAtPoints.at(pt)[-j]==true)
			  {
			    TracksMeetAtPoints.at(pt)[-i]=true;
			    MeetPoints.at(pt).push_back(ExtrapPoint1);
			    MeetPoints.at(pt).push_back(ExtrapPoint2);
			    ThisTrackAlreadyMet=true;
			  }
		      }
		    if(!ThisTrackAlreadyMet)
		      {
			size_t pt = TracksMeetAtPoints.size();
			TracksMeetAtPoints.push_back(std::map<int,bool>());
			TracksMeetAtPoints.at(pt)[-i]=true;
			TracksMeetAtPoints.at(pt)[-j]=true;
			MeetPoints.push_back(std::vector<TVector3>());
			MeetPoints.at(pt).push_back(ExtrapPoint1);
			MeetPoints.at(pt).push_back(ExtrapPoint2);
		      }

		  }

		GetImpact(TrackEnd2s[i],TrackDir2s[i],TrackEnd2s[j],TrackDir2s[j], impact, disti, distj);
		if((impact < fJoinThreshold)
		   && (fabs(disti) < fExtrapDistance )
		   && (fabs(distj) < fExtrapDistance )
		   && (disti > 0 )
		   && (distj > 0 ))
		  {
		    bool ThisTrackAlreadyMet =false;
		    TVector3 ExtrapPoint1 = TrackEnd1s[i] + TrackDir2s[j].Unit()*disti;
		    TVector3 ExtrapPoint2 = TrackEnd2s[j] + TrackDir2s[j].Unit()*distj;
		    
		    for(size_t pt=0; pt!=TracksMeetAtPoints.size(); ++pt)
		      {
			if(TracksMeetAtPoints.at(pt)[i]==true)
			  {
			    TracksMeetAtPoints.at(pt)[j]=true;
			    MeetPoints.at(pt).push_back(ExtrapPoint1);
			    MeetPoints.at(pt).push_back(ExtrapPoint2);
			    ThisTrackAlreadyMet=true;
			    
			  }
			else if(TracksMeetAtPoints.at(pt)[j]==true)
			  {
			    TracksMeetAtPoints.at(pt)[i]=true;
			    MeetPoints.at(pt).push_back(ExtrapPoint1);
			    MeetPoints.at(pt).push_back(ExtrapPoint2);
			    ThisTrackAlreadyMet=true;
			  }
		      }
		    if(!ThisTrackAlreadyMet)
		      {
			size_t pt = TracksMeetAtPoints.size();
			TracksMeetAtPoints.push_back(std::map<int,bool>());
			TracksMeetAtPoints.at(pt)[i]=true;
			TracksMeetAtPoints.at(pt)[j]=true;
			MeetPoints.push_back(std::vector<TVector3>());
			MeetPoints.at(pt).push_back(ExtrapPoint1);
			MeetPoints.at(pt).push_back(ExtrapPoint2);
		      }
		    
		    
		  }
		
		GetImpact(TrackEnd1s[i],TrackDir1s[i],TrackEnd2s[j],TrackDir2s[j], impact, disti, distj);
		if((impact < fJoinThreshold)
		   && (fabs(disti) < fExtrapDistance )
		   && (fabs(distj) < fExtrapDistance )
		   && (disti < 0 )
		   && (distj > 0 ))
		  {
		    bool ThisTrackAlreadyMet =false;
		    TVector3 ExtrapPoint1 = TrackEnd1s[i] + TrackDir2s[j].Unit()*disti;
		    TVector3 ExtrapPoint2 = TrackEnd2s[j] + TrackDir2s[j].Unit()*distj;
		    
		    for(size_t pt=0; pt!=TracksMeetAtPoints.size(); ++pt)
		      {
			if(TracksMeetAtPoints.at(pt)[-i]==true)
			  {
			    TracksMeetAtPoints.at(pt)[j]=true;
			    MeetPoints.at(pt).push_back(ExtrapPoint1);
			    MeetPoints.at(pt).push_back(ExtrapPoint2);
			    ThisTrackAlreadyMet=true;
			    
			  }
			else if(TracksMeetAtPoints.at(pt)[j]==true)
			  {
			    TracksMeetAtPoints.at(pt)[-i]=true;
			    MeetPoints.at(pt).push_back(ExtrapPoint1);
			    MeetPoints.at(pt).push_back(ExtrapPoint2);
			    ThisTrackAlreadyMet=true;
			  }
		      }
		    if(!ThisTrackAlreadyMet)
		      {
			size_t pt = TracksMeetAtPoints.size();
			TracksMeetAtPoints.push_back(std::map<int,bool>());
			TracksMeetAtPoints.at(pt)[-i]=true;
			TracksMeetAtPoints.at(pt)[j]=true;
			MeetPoints.push_back(std::vector<TVector3>());
			MeetPoints.at(pt).push_back(ExtrapPoint1);
			MeetPoints.at(pt).push_back(ExtrapPoint2);
		      }
	  
				    
		  }
		
	      }
	  }
      }



    for(size_t pt=0; pt!=TracksMeetAtPoints.size(); ++pt)
      {
	std::cout<<"Making adjustments for meet pt " <<pt<<std::endl;
	TVector3 FinalMeetPt;
	for(size_t i=0; i!=MeetPoints.at(pt).size(); ++i)
	  {
	    FinalMeetPt += MeetPoints.at(pt).at(i);
	  } 
	FinalMeetPt *= (1./float(MeetPoints.at(pt).size()));
	
	for(size_t i=0; i!=BTracks.size(); ++i)
	  {
	    if(TracksMeetAtPoints.at(pt)[i]==true)
	      {
		TVector3 NewDirection = (FinalMeetPt - TrackEnd2s[i])*0.5;
		double NewPos[3];
		double NewDir[3];
		double Err[3];
		
		for(size_t n=0; n!=3; ++n)
		  {
		    NewPos[n]=FinalMeetPt[n];
		    NewDir[n]=NewDirection[n];
		  }
		
		recob::Seed NewSeed(NewPos, NewDir, Err, Err);
		BTracks.at(i).GetSeedVector().push_back(NewSeed);
	      }
	    else if(TracksMeetAtPoints.at(pt)[-i]==true)
	      {
		TVector3 NewDirection = -(FinalMeetPt - TrackEnd2s[i])*0.5;
		double NewPos[3];
		double NewDir[3];
		double Err[3];
		
		for(size_t n=0; n!=3; ++n)
		  {
		    NewPos[n]=FinalMeetPt[n];
		    NewDir[n]=NewDirection[n];
		  }
		
		recob::Seed NewSeed(NewPos, NewDir, Err, Err);
		BTracks.at(i).GetSeedVector().insert(BTracks.at(i).GetSeedVector().begin(), NewSeed);

	      }
	  }
      }
    
    for(size_t i=0; i!=BTracks.size(); ++i)
      {
	tracksout->push_back(*(BTracks.at(i).GetBaseTrack()));
      }
    

    evt.put(std::move(tracksout));
    evt.put(std::move(verticesout));
    evt.put(std::move(tvassn));
    

  }

  void BezierTrackJoiner::GetImpact(TVector3 t1pt, TVector3 t1dir, TVector3 t2pt, TVector3 t2dir, double& ImpactParam, double& Dist1, double& Dist2)
  {
    double lambda1 = ((t2pt-t1pt).Cross(t2dir)).Dot(t1dir.Cross(t2dir)) / (t1dir.Cross(t2dir)).Mag2();
    double lambda2 = ((t1pt-t2pt).Cross(t1dir)).Dot(t2dir.Cross(t1dir)) / (t2dir.Cross(t1dir)).Mag2();
    TVector3 c = t1pt + t1dir*lambda1 - t2pt - t2dir*lambda2;
    
    ImpactParam = c.Mag();
    Dist1 = lambda1 * t1dir.Mag();
    Dist2 = lambda2 * t2dir.Mag();
  }



   
  //-------------------------------------------------------------------------
  void BezierTrackJoiner::endJob()
  {

  }
}
