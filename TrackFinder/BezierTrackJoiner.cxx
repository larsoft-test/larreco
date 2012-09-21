//
// Name: BezierTrackJoiner.cxx
//
// Purpose: Implementation file for module BezierTrackJoiner.
//
// Ben Jones, MIT
//

#include <vector>
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "TrackFinder/BezierTrackJoiner.h"
#include "Geometry/geo.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "RecoBase/Hit.h"
#include "RecoBase/Seed.h"
#include "RecoBase/Track.h"
#include "RecoBase/Prong.h"
#include "RecoBase/Vertex.h"
#include "TrackFinder/BezierTrack.h"
#include "Utilities/AssociationUtil.h"



namespace trkf {

  //-------------------------------------------------------------------------
  BezierTrackJoiner::BezierTrackJoiner(const fhicl::ParameterSet& pset)
  {
    reconfigure(pset);
    produces< std::vector<recob::Track> >();
    produces< std::vector<recob::Vertex> >();
    produces< art::Assns<recob::Track, recob::Hit> >();
    produces< art::Assns<recob::Vertex, recob::Track> >();
    produces< art::Assns<recob::Vertex, recob::Hit> >();
  }

  //-------------------------------------------------------------------------
  BezierTrackJoiner::~BezierTrackJoiner()
  {
  }

  //-------------------------------------------------------------------------
  void BezierTrackJoiner::reconfigure(fhicl::ParameterSet const& pset)
  {
    fTrackModuleLabel   = pset.get<std::string>("TrackModuleLabel");
    fHitModuleLabel     = pset.get<std::string>("HitModuleLabel");
    fJoinThreshold      = pset.get<double>("JoinThreshold");
    fVertexAngle        = pset.get<double>("VertexAngle");
    fDistanceForAngle   = pset.get<double>("DistanceForAngle");   
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
    std::unique_ptr< art::Assns<recob::Track, recob::Hit > > assn( new art::Assns<recob::Track, recob::Hit>);
    std::unique_ptr< art::Assns<recob::Vertex, recob::Track > > tvassn( new art::Assns<recob::Vertex, recob::Track>);
    std::unique_ptr< art::Assns<recob::Vertex, recob::Hit > > vhassn( new art::Assns<recob::Vertex, recob::Hit>);
    
        
    mf::LogVerbatim("BezierTrackJoiner") << "Getting tracks ";

    // Extract track objects from event

    art::Handle< std::vector<recob::Track> > trackh;
    evt.getByLabel(fTrackModuleLabel, trackh);

    art::FindManyP<recob::Hit> fmh(trackh, evt, fTrackModuleLabel);

    std::vector< std::vector<art::Ptr<recob::Hit> > > trackHits;
    art::PtrVector<recob::Track>   Tracks;
    Tracks.clear();
    for(unsigned int i=0; i < trackh->size(); ++i)
      {
	art::Ptr<recob::Track> track(trackh,i);
	Tracks.push_back(track);

	// get the hits associated with this track
	trackHits.push_back(fmh.at(i));
      }

    


    mf::LogVerbatim("BezierTrackJoiner") << "Make bez tracks ";
    // Make bezier tracks
    std::vector<trkf::BezierTrack*> BTracks;
    BTracks.clear();
    for(size_t i=0; i!=Tracks.size(); i++)
      BTracks.push_back(new BezierTrack(*Tracks.at(i)));


    std::vector<trkf::BezierTrack*> JoinedTracks = BTracks;

    size_t PrevSize=BTracks.size()+1;
     
    // If we joined some tracks last time, try again
    //  until no more joins
    while(JoinedTracks.size()<PrevSize)
      {
	PrevSize = JoinedTracks.size();
	JoinedTracks = MakeTouchJoins( JoinedTracks, fJoinThreshold);
	mf::LogVerbatim("BezierTrackJoiner") << "Making touch joins  : " << PrevSize << " " << JoinedTracks.size();

      }

    PrevSize=BTracks.size()+1;


    while(JoinedTracks.size()<PrevSize)
      {
	PrevSize = JoinedTracks.size();
	JoinedTracks = MakeGlancingJoins( JoinedTracks, fJoinThreshold);
	mf::LogVerbatim("BezierTrackJoiner") << "Making glancing joins  : " << PrevSize << " " << JoinedTracks.size();
	
      }
    

    std::vector<std::vector<int> > TracksToVertices;
    std::vector<TVector3> Vertices = MakeVertexCandidates(JoinedTracks, TracksToVertices);

    
    

    for(size_t i=0; i!=JoinedTracks.size(); ++i)
      {
	tracksout->push_back(JoinedTracks.at(i)->GetBaseTrack());

	///\todo need to associate hits to the tracks
      }

    for(size_t i=0; i!=Vertices.size(); ++i)
      {
	double xyz[3];
	for(size_t j=0; j!=3; j++)
	  xyz[j]=Vertices.at(i)[j];
	
	verticesout->push_back(recob::Vertex(xyz, i)); 

	// associate tracks to this vertex
	for(size_t t = 0; t < TracksToVertices[i].size(); ++t)
	  util::CreateAssn(*this, evt, *verticesout, *tracksout, *tvassn, 
			   (size_t)TracksToVertices[i][t], (size_t)TracksToVertices[i][t]+1);

	///\todo need to associate hits to the vertices

      }
    
    mf::LogVerbatim("BezierTrackJoiner") << "storing products ";

    // Store the fruits of our labour in the event
    mf::LogVerbatim("BezierTrackJoiner") << "Storing in evt";
    evt.put(std::move(tracksout));
    evt.put(std::move(verticesout));
    evt.put(std::move(tvassn));
    

  }


  //-------------------------------------------------------------------------
  std::vector<BezierTrack*> BezierTrackJoiner::MakeTouchJoins(std::vector<BezierTrack*> BTracks, double fJoinThreshold)
  {

    std::cout<<"looping tracks " <<std::endl;
    std::map<int,bool> UsedTracks;
    std::vector<BezierTrack*> ReturnVector;
    
    ReturnVector.clear();
    UsedTracks.clear();

    // Loop through track pairs looking for joins
    for(size_t i=1; i<BTracks.size(); ++i)
      for(size_t j=0; j<BTracks.size(); ++j)
	{
	  if((!UsedTracks[i])&&(!UsedTracks[j]))
	    {
	      TVector3 Endi1 = BTracks.at(i)->GetTrackPointV(1);
	      TVector3 Endj0 = BTracks.at(j)->GetTrackPointV(0);
	   
	      
	      std::cout<<"EndSeparation " <<i<<" " << j <<" " << (Endi1-Endj0).Mag()<<std::endl;
		  
	      if((Endi1-Endj0).Mag()<fJoinThreshold)
		{  
		  TVector3 Diri = BTracks.at(i)->GetTrackPointV(1) - BTracks.at(i)->GetTrackPointV(1-fDistanceForAngle);
		  TVector3 Dirj = (BTracks.at(j)->GetTrackPointV(1) - BTracks.at(j)->GetTrackPointV(1-fDistanceForAngle));

		  std::cout<<"checking angle  "<<Diri.Angle(Dirj) <<std::endl;

		  if(Diri.Angle(Dirj)<fVertexAngle)
		    {
		      std::cout<<"Making a join"<<std::endl;
		      UsedTracks[i]=UsedTracks[j]=true;
		      recob::Track Combined = BTracks.at(i)->GetJoinedBaseTrack(BTracks.at(j));		      		      
		      ReturnVector.push_back(new BezierTrack(Combined));
		    } // Angle small enough
		} // within join distance threshold
	    } // unused tracks
	} // loop i,j
    
    
    // Add into the return vector any that were not joined
    for(size_t i=0; i!=BTracks.size(); ++i)
      if(!UsedTracks[i]){
	ReturnVector.push_back(BTracks.at(i));
      }

    return ReturnVector;
  }

  /*
  std::vector<BezierTrack*> BezierTrackJoiner::MakeExtrapJoins(std::vector<BezierTrack*> BTracks, double fMissThreshold, double fAngleThreshold, double fJumpDistance, double fLengthFraction)
  {
    std::map<int,bool> UsedTracks;
    std::vector<BezierTrack*> ReturnVector;
    
    ReturnVector.clear();
    UsedTracks.clear();

    // Loop through track pairs looking for joins
    for(size_t i=1; i<BTracks.size(); ++i)
      for(size_t j=0; j<BTracks.size(); ++j)
	{
	  if((!UsedTracks[i])&&(!UsedTracks[j])&&(i!=j))
	    {
	      TVector3 Endi1 = BTracks.at(i)->GetTrackPointV(1);
	      TVector3 Diri1 = (Endi1 - BTracks.at(i)->GetTrackPointV(1-fLengthFraction)).Unit();
	      TVector3 Endj0 = BTracks.at(j)->GetTrackPointV(0);
	      TVector3 Dirj0 = (BTracks.at(i)->GetTrackPointV(fLengthFraction) - Endj0).Unit();
	      
	      if((Diri1.Angle(Dirj0)) < fAngleThreshold)
		{
		  if(evd

		}
	    }
	}
  
  }

  
  */

  //-------------------------------------------------------------------------
  std::vector<BezierTrack*> BezierTrackJoiner::MakeGlancingJoins(std::vector<BezierTrack*> BTracks, double fJoinThreshold)
  {

    std::cout<<"looping tracks " <<std::endl;
    std::map<int,bool> UsedTracks;
    std::vector<BezierTrack*> ReturnVector;
    
    ReturnVector.clear();
    UsedTracks.clear();

    // Loop through track pairs looking for joins
    for(size_t i=1; i<BTracks.size(); ++i)
      for(size_t j=0; j<BTracks.size(); ++j)
	{
	  if((!UsedTracks[i])&&(!UsedTracks[j])&&(i!=j))
	    {
	      TVector3 Endi1 = BTracks.at(i)->GetTrackPointV(1);
	      TVector3 Endj0 = BTracks.at(j)->GetTrackPointV(0);
	    
	      double s1,d1, s2,d2;
	      BTracks.at(i)->GetClosestApproach(Endj0,s1,d1);
	      BTracks.at(j)->GetClosestApproach(Endi1,s2,d2);

	      std::cout<<"Glancing join distances  : "<<d1<<" " << d2<<std::endl;
	      std::cout<<"Glancing join s values   : "<<s1<<" " << s2<<std::endl;
	      
	      if((d1<d2)&&(d1<fJoinThreshold))
		{
		  std::cout<<"Segs from track end " << BTracks.at(i)->WhichSegment(s1)<<" " <<BTracks.at(i)->NSegments()<<std::endl;
		  if(BTracks.at(i)->WhichSegment(s1)==(BTracks.at(i)->NSegments()-1))
		    {
		      recob::Track TheTrack 
			= BTracks.at(i)->GetJoinedPartBaseTrack(BTracks.at(j),
								0, BTracks.at(i)->WhichSegment(s1),
								0, BTracks.at(j)->NSegments()-1);
		      UsedTracks[i]=UsedTracks[j]=true;		      
		      ReturnVector.push_back(new BezierTrack(TheTrack));
		    }
		}
	      if((d2<d1)&&(d2<fJoinThreshold))
		{
		  std::cout<<"Segs from track end " << BTracks.at(j)->WhichSegment(s2)<<" " << BTracks.at(j)->NSegments()<<std::endl;
		  if(BTracks.at(j)->WhichSegment(s2)==0)
		    {
		      recob::Track TheTrack 
			= BTracks.at(i)->GetJoinedPartBaseTrack(BTracks.at(j),
								0, BTracks.at(i)->NSegments()-1,
								1, BTracks.at(j)->NSegments()-1);
		      UsedTracks[i]=UsedTracks[j]=true;
		      ReturnVector.push_back(new BezierTrack(TheTrack));
		    }		     
		}    
	      
	    }
      
	}    
    
    // Add into the return vector any that were not joined
    for(size_t i=0; i!=BTracks.size(); ++i)
      if(!UsedTracks[i])
	ReturnVector.push_back(BTracks.at(i));
    
    return ReturnVector;
  }


  //-------------------------------------------------------------------------
  std::vector<TVector3> BezierTrackJoiner::MakeVertexCandidates(std::vector<BezierTrack*> BTracks,std::vector<std::vector<int> >& TracksToVertices)
  {
    std::cout<<"looping tracks " <<std::endl;
    std::map<int,bool> UsedTracks;
    std::vector<TVector3> ReturnVector;
    
    ReturnVector.clear();
    UsedTracks.clear();

    // Loop through track pairs looking for joins
    for(size_t i=1; i<BTracks.size(); ++i)
      {
	TVector3 Endi0 = BTracks.at(i)->GetTrackPointV(0);

	// first check if track meets any existing vertices
	for(size_t jvtx=0; jvtx!=ReturnVector.size(); jvtx++)
	  {
	    if((Endi0-ReturnVector.at(jvtx)).Mag()<fJoinThreshold)
	      {
		TracksToVertices.at(jvtx).push_back(i);
		UsedTracks[i]=true;		
	      }	    
	  }

	// then check if track meets any other tracks to form
	//  a new vertex
	for(size_t j=0; j<BTracks.size(); ++j)
	  {
	    if((!UsedTracks[i])&&(!UsedTracks[j])&&(i!=j))
	      {

		TVector3 Endj0 = BTracks.at(j)->GetTrackPointV(0);
		
		if((Endi0-Endj0).Mag()<fJoinThreshold)
		  {
		    TVector3 VertexPos=0.5*(Endi0+Endj0);
		    ReturnVector.push_back(VertexPos);
		    std::cout<<"Found Vertex at " <<VertexPos[0]<<
		      " " <<  VertexPos[1]<< " " << VertexPos[2]<<std::endl;
		    std::vector<int> TracksThisVtx;
		    TracksThisVtx.push_back(i);
		    TracksThisVtx.push_back(j);

		    UsedTracks[i]=UsedTracks[j]=true;
		    
		    TracksToVertices.push_back(TracksThisVtx);
		    
		    
		  }
	    }
	}
      }
    return ReturnVector;
    
  }

   
  //-------------------------------------------------------------------------
  void BezierTrackJoiner::endJob()
  {

  }
}
