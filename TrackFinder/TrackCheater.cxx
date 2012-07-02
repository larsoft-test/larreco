////////////////////////////////////////////////////////////////////////
/// \file    TrackCheater.cxx
/// \brief   create perfectly reconstructed prongs
/// \version $Id: GeometryTest.cxx,v 1.1 2011/02/17 01:45:48 brebel Exp $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#include <vector>

// ROOT includes
#include "TVector3.h"

// LArSoft includes
#include "MCCheater/BackTracker.h"
#include "TrackFinder/TrackCheater.h"
#include "RecoBase/recobase.h"
#include "Utilities/AssociationUtil.h"
#include "Simulation/sim.h"
#include "Simulation/SimListUtils.h"
#include "Utilities/PhysicalConstants.h"
#include "Utilities/DetectorProperties.h"

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace trkf{

  //--------------------------------------------------------------------
  TrackCheater::TrackCheater(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);

    produces< std::vector<recob::Track>                   >();
    produces< std::vector<recob::SpacePoint>              >();
    produces< art::Assns<recob::Track, recob::Cluster>    >();
    produces< art::Assns<recob::Track, recob::SpacePoint> >();
    produces< art::Assns<recob::Track, recob::Hit>        >();
  }

  //--------------------------------------------------------------------
  TrackCheater::~TrackCheater()
  {
  }

  //--------------------------------------------------------------------
  void TrackCheater::reconfigure(fhicl::ParameterSet const& pset)
  {
    fCheatedClusterLabel = pset.get< std::string >("CheatedClusterLabel", "cluster" );
    fG4ModuleLabel       = pset.get< std::string >("G4ModuleLabel",       "largeant");

    return;
  }

  //--------------------------------------------------------------------
  void TrackCheater::produce(art::Event& evt)
  {
    art::ServiceHandle<cheat::BackTracker>       bt;
    art::ServiceHandle<geo::Geometry>            geo;
    art::ServiceHandle<util::DetectorProperties> detp;

    // grab the clusters that have been reconstructed
    art::Handle< std::vector<recob::Cluster> > clustercol;
    evt.getByLabel(fCheatedClusterLabel, clustercol);

    // make a vector of them - we aren't writing anything out to a file
    // so no need for a art::PtrVector here
    std::vector< art::Ptr<recob::Cluster> > clusters;
    art::fill_ptr_vector(clusters, clustercol);
    
    // loop over the clusters and figure out which particle contributed to each one
    std::vector< art::Ptr<recob::Cluster> >::iterator itr = clusters.begin();

    // make a map of vectors of art::Ptrs keyed by eveID values
    std::map< int, std::vector< art::Ptr<recob::Cluster> > > eveClusterMap;
    std::map< int, std::vector< art::Ptr<recob::Cluster> > >::iterator clusterMapItr = eveClusterMap.begin();

    // loop over all clusters and fill in the map
    while( itr != clusters.end() ){

      // in the ClusterCheater module we set the cluster ID to be 
      // the eve particle track ID*1000 + plane*100 + tpc number.  The
      // floor function on the cluster ID / 1000 will give us
      // the eve track ID
      int eveID = floor((*itr)->ID()/1000.);

      clusterMapItr = eveClusterMap.find(eveID);
	
      // is this id already in the map, if so extend the collection 
      // by one cluster, otherwise make a new collection and put it in
      // the map
      if( clusterMapItr != eveClusterMap.end() ){
	  ((*clusterMapItr).second).push_back((*itr));
      }
      else{
	std::vector< art::Ptr<recob::Cluster> > clustervec;
	clustervec.push_back(*itr);
	eveClusterMap[eveID] = clustervec;
      }

      itr++;
    }// end loop over clusters

    // loop over the map and make prongs
    std::auto_ptr< std::vector<recob::Track> >                   trackcol(new std::vector<recob::Track>);
    std::auto_ptr< std::vector<recob::SpacePoint> >              spcol  (new std::vector<recob::SpacePoint>);
    std::auto_ptr< art::Assns<recob::Track, recob::SpacePoint> > tspassn(new art::Assns<recob::Track, recob::SpacePoint>);
    std::auto_ptr< art::Assns<recob::Track, recob::Cluster> >    tcassn (new art::Assns<recob::Track, recob::Cluster>);
    std::auto_ptr< art::Assns<recob::Track, recob::Hit> >        thassn (new art::Assns<recob::Track, recob::Hit>);

    for(clusterMapItr = eveClusterMap.begin(); clusterMapItr != eveClusterMap.end(); clusterMapItr++){

      // separate out the hits for each particle into the different views
      std::vector< art::Ptr<recob::Cluster> > eveClusters( (*clusterMapItr).second );

      art::PtrVector<recob::Cluster> ptrvs;

      sim::Particle *part = bt->ParticleList()[(*clusterMapItr).first];

      // is this prong electro-magnetic in nature or 
      // hadronic/muonic?  EM --> shower, everything else is a track
      if( abs(part->PdgCode()) != 11  &&
	  abs(part->PdgCode()) != 22  &&
	  abs(part->PdgCode()) != 111 ){

	mf::LogInfo("TrackCheater") << "prong of " << (*clusterMapItr).first 
				    << " is a track with pdg code "
				    << part->PdgCode();

	for(size_t c = 0; c < eveClusters.size(); ++c) ptrvs.push_back(eveClusters[c]);

	// get the particle for this track ID
	size_t numtraj = part->NumberTrajectoryPoints();
	std::vector<TVector3> points(numtraj);
	std::vector<TVector3> dirs(numtraj);

	size_t nviews = geo->Nviews();
	std::vector< std::vector<double> > dQdx(nviews);

	// loop over the particle trajectory
	size_t spStart = spcol->size();
	for(size_t t = 0; t < numtraj; ++t){
	  points[t].SetXYZ(part->Vx(t), part->Vy(t), part->Vz(t));
	  dirs[t].SetXYZ(part->Px(t)/part->P(), part->Py(t)/part->P(), part->Pz(t)/part->P());

	  // use the same dQdx in each view and just make it direct energy -> charge conversion
	  double eLoss = 0.;
	  double dx    = 0.;
	  if(t > 0){
	    eLoss = fabs(part->E(t)  - part->E(t-1))*util::kGeVToElectrons*detp->ElectronsToADC();
	    dx    = fabs(part->Vx(t) - part->Vx(t-1));
	  }
	  dQdx[0].push_back(eLoss/dx);
	  dQdx[1].push_back(eLoss/dx);
	  dQdx[2].push_back(eLoss/dx);

	  double xyz[3]    = {part->Vx(t), part->Vy(t), part->Vz(t)};
	  double xyzerr[6] = {1.e-3};
	  double chisqr    = 0.9;

	  // make the space point and set its ID and XYZ
	  recob::SpacePoint sp(&xyz[0], xyzerr, chisqr, eveClusters[0]->ID()*10000 + t);
	  spcol->push_back(sp);
	}
	
	size_t spEnd = spcol->size();
	
	// add a track to the collection.  Make the track
	// ID be the same as the track ID for the eve particle
	std::vector<double> momentum(2);
	momentum[0] = part->P();
	momentum[1] = part->P(numtraj-1);
	trackcol->push_back(recob::Track(points, dirs, dQdx, momentum, (*clusterMapItr).first));

	// associate the track with its clusters
	util::CreateAssn(*this, evt, *trackcol, ptrvs, *tcassn);

	art::FindManyP<recob::Hit> fmh(ptrvs, evt, fCheatedClusterLabel);

	// assume the input tracks were previously associated with hits
	for(size_t p = 0; p < ptrvs.size(); ++p){
	  std::vector< art::Ptr<recob::Hit> > hits = fmh.at(p);
	  util::CreateAssn(*this, evt, *trackcol, hits, *thassn);
	}

	// associate the track to the space points
	util::CreateAssn(*this, evt, *trackcol, *spcol, *tspassn, spStart, spEnd);

	mf::LogInfo("TrackCheater") << "adding track: \n" 
				     << trackcol->back()
				     << "\nto collection.";

      }//end if this is a track

    } // end loop over the map

    evt.put(trackcol);
    evt.put(spcol);
    evt.put(tcassn);
    evt.put(thassn);
    evt.put(tspassn);

    return;

  } // end produce

} // end namespace
