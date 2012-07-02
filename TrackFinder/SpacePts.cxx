////////////////////////////////////////////////////////////////////////
//
// \file SpacePts.cxx
//
////////////////////////////////////////////////////////////////////////

// C++ includes
#include <math.h>
#include <algorithm>

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
#include "TrackFinder/SpacePts.h"
#include "Geometry/geo.h"
#include "RecoBase/recobase.h"
#include "Utilities/LArProperties.h"
#include "Utilities/AssociationUtil.h"

// ROOT includes
#include "TVectorD.h"
#include "TMath.h"
#include "TGraph.h"
#include "TF1.h"

namespace trkf {
struct SortByWire 
{
  bool operator() (art::Ptr<recob::Hit> const& h1, art::Ptr<recob::Hit> const& h2) const 
  { return 
      h1->Wire()->RawDigit()->Channel() < 
      h2->Wire()->RawDigit()->Channel() ;
  }
};
} // end namespace

///\todo: SpacePts looks very detector specific, can it be cleaned up?

//-------------------------------------------------
trkf::SpacePts::SpacePts(fhicl::ParameterSet const& pset)
{
  this->reconfigure(pset);
  
  produces< std::vector<recob::Track>                   >();
  produces< std::vector<recob::SpacePoint>              >();
  produces< art::Assns<recob::Track, recob::SpacePoint> >();
  produces< art::Assns<recob::Track, recob::Cluster>    >();
  produces< art::Assns<recob::Track, recob::Hit>        >();
  
}

//-------------------------------------------------
trkf::SpacePts::~SpacePts()
{
}

void trkf::SpacePts::reconfigure(fhicl::ParameterSet const& pset) 
{
  fPreSamplings           = pset.get< double >("TicksOffset");
  ftmatch                 = pset.get< int    >("TMatch");
  fClusterModuleLabel     = pset.get< std::string >("ClusterModuleLabel");
  fEndPoint2DModuleLabel  = pset.get< std::string >("EndPoint2DModuleLabel");
  fvertexclusterWindow    = pset.get< double >("vertexclusterWindow");
}

//-------------------------------------------------
void trkf::SpacePts::beginJob()
{
}

void trkf::SpacePts::endJob()
{
}

//------------------------------------------------------------------------------------//
void trkf::SpacePts::produce(art::Event& evt)
{ 

  
  // get services
  art::ServiceHandle<geo::Geometry> geom;
  art::ServiceHandle<util::LArProperties> larprop;
  
  //////////////////////////////////////////////////////
  // Make a std::auto_ptr<> for the thing you want to put into the event
  // because that handles the memory management for you
  //////////////////////////////////////////////////////
  std::auto_ptr<std::vector<recob::Track>      >              tcol (new std::vector<recob::Track>);	   
  std::auto_ptr<std::vector<recob::SpacePoint> > 	      spcol(new std::vector<recob::SpacePoint>);
  std::auto_ptr<art::Assns<recob::Track, recob::SpacePoint> > tspassn(new art::Assns<recob::Track, recob::SpacePoint>);
  std::auto_ptr<art::Assns<recob::Track, recob::Cluster> >    tcassn(new art::Assns<recob::Track, recob::Cluster>);
  std::auto_ptr<art::Assns<recob::Track, recob::Hit> >        thassn(new art::Assns<recob::Track, recob::Hit>);
  
  // define TPC parameters
  TString tpcName = geom->GetLArTPCVolumeName();
  
  //TPC dimensions
  double YC =  (geom->DetHalfHeight())*2.; // TPC height in cm
  double Angle = geom->Plane(1).Wire(0).ThetaZ(false)-TMath::Pi()/2.; // wire angle with respect to the vertical direction
  // Parameters temporary defined here, but possibly to be retrieved somewhere in the code
  double timetick = 0.198;    //time sample in us
  double presamplings = fPreSamplings; // 60.;
  const double wireShift=50.; // half the number of wires from the Induction(Collection) plane intersecting with a wire from the Collection(Induction) plane.
  double plane_pitch = geom->PlanePitch(0,1);   //wire plane pitch in cm 
  double wire_pitch = geom->WirePitch(0,1,0);    //wire pitch in cm
  double Efield_drift = 0.5;  // Electric Field in the drift region in kV/cm
  double Efield_SI = 0.7;     // Electric Field between Shield and Induction planes in kV/cm
  double Efield_IC = 0.9;     // Electric Field between Induction and Collection planes in kV/cm
  double Temperature = 90.;  // LAr Temperature in K
  
  double driftvelocity = larprop->DriftVelocity(Efield_drift,Temperature);    //drift velocity in the drift region (cm/us)
  double driftvelocity_SI = larprop->DriftVelocity(Efield_SI,Temperature);    //drift velocity between shield and induction (cm/us)
  double driftvelocity_IC = larprop->DriftVelocity(Efield_IC,Temperature);    //drift velocity between induction and collection (cm/us)
  double timepitch = driftvelocity*timetick;                         //time sample (cm) 
  double tSI = plane_pitch/driftvelocity_SI/timetick;                   //drift time between Shield and Collection planes (time samples)
  double tIC = plane_pitch/driftvelocity_IC/timetick;                //drift time between Induction and Collection planes (time samples)
  
  
  // get input Cluster object(s).
  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fClusterModuleLabel,clusterListHandle);
  
  // get input EndPoint2D object(s).
  art::Handle< std::vector<recob::EndPoint2D> > endpointListHandle;
  evt.getByLabel(fEndPoint2DModuleLabel,endpointListHandle); 
  
  art::PtrVector<recob::EndPoint2D> endpointlist;
  if(evt.getByLabel(fEndPoint2DModuleLabel,endpointListHandle))
    for (unsigned int i = 0; i < endpointListHandle->size(); ++i){
      art::Ptr<recob::EndPoint2D> endpointHolder(endpointListHandle,i);
      endpointlist.push_back(endpointHolder);
    }
   
   // Declare some vectors..
   // Induction
   std::vector<double> Iwirefirsts;       // in cm
   std::vector<double> Iwirelasts;        // in cm
   std::vector<double> Itimefirsts;       // in cm
   std::vector<double> Itimelasts;        // in cm
   std::vector<double> Itimefirsts_line;  // in cm
   std::vector<double> Itimelasts_line;   // in cm
   std::vector < std::vector< art::Ptr<recob::Hit> > > IclusHitlists;
   std::vector<unsigned int> Icluster_count; 

   // Collection
   std::vector<double> Cwirefirsts;       // in cm
   std::vector<double> Cwirelasts;        // in cm
   std::vector<double> Ctimefirsts;       // in cm
   std::vector<double> Ctimelasts;        // in cm
   std::vector<double> Ctimefirsts_line;  // in cm
   std::vector<double> Ctimelasts_line;   // in cm
   std::vector< std::vector< art::Ptr<recob::Hit> > > CclusHitlists;
   std::vector<unsigned int> Ccluster_count; 

   art::FindManyP<recob::Hit> fm(clusterListHandle, evt, fClusterModuleLabel);

   for(unsigned int ii = 0; ii < clusterListHandle->size(); ++ii){

      art::Ptr<recob::Cluster> cl(clusterListHandle, ii);
      
      // Figure out which View the cluster belongs to 
      /*
      int clPlane = cl->View()-1;
      std::cout << "SpacePts: Plane/View/SignalType ..." << clPlane <<" "<<cl->View()<<" "<<geom->Plane(clPlane).SignalType() << std::endl; 
      */
      // Gaaaaaah! Change me soon!!! But, for now, 
      // let's just chuck one plane's worth of info. EC, 30-Mar-2011.
      ///\todo This is really horrendous code.  Never, Never, Never test on Detector name!!!!
      if (cl->View() == geo::kW && !geom->GetDetectorName().Contains(geo::kArgoNeuT)) continue; 
      mf::LogWarning("SpacePts:") << "!!! 3 Plane detectors only using kV and kW Views for now!!!!";
       //only consider merged-lines that are associated with the vertex.
       //this helps get rid of through-going muon background -spitz                  
       int vtx2d_w = -99999;
       double vtx2d_t = -99999;
       bool found2dvtx = false;
       
      for (unsigned int j = 0; j<endpointlist.size();j++){
	if (endpointlist[j]->View() == cl->View()){
	  vtx2d_w = endpointlist[j]->WireNum();
	  vtx2d_t = endpointlist[j]->DriftTime();
	  found2dvtx = true;
	  break;
	}
      }
      if (found2dvtx){
	double w = cl->StartPos()[0];
	double t = cl->StartPos()[1];
	double dtdw = cl->dTdW();
	double t_vtx = t+dtdw*(vtx2d_w-w);
	double dis = TMath::Abs(vtx2d_t-t_vtx);
	if (dis>fvertexclusterWindow)	  continue;	
      }
      //else continue; //what to do if a 2D vertex is not found? perhaps vertex finder was not even run.
    
     // Some variables for the hit
      unsigned int channel;  //channel number
      float time;            //hit time at maximum
      unsigned int wire;     //hit wire number
      unsigned int plane;    //hit plane number
      unsigned int tpc;      //hit tpc number
      unsigned int cstat;    //hit cryostat number
      
      
      std::vector< art::Ptr<recob::Hit> > hitlist = fm.at(ii);
      std::sort(hitlist.begin(), hitlist.end(), trkf::SortByWire());
      
      TGraph *the2Dtrack = new TGraph(hitlist.size());
      
      std::vector<double> wires;
      std::vector<double> times;
     
      
      int np=0;
      //loop over cluster hits
      for(std::vector< art::Ptr<recob::Hit> >::const_iterator theHit = hitlist.begin(); theHit != hitlist.end();  theHit++){
	//recover the Hit
	//      recob::Hit* theHit = (recob::Hit*)(*hitIter);
	time = (*theHit)->PeakTime() ;
	
	time -= presamplings;
	
	channel = (*theHit)->Channel();
	geom->ChannelToWire(channel,cstat, tpc,plane,wire);
	
	//correct for the distance between wire planes
	//          if(plane==0) time -= tSI;         // Induction
	//          if(plane==1) time -= (tSI+tIC);   // Collection
	
	
	if(geom->Cryostat(cstat).TPC(tpc).Plane(plane).SignalType() == geo::kCollection) 
	  time -= tIC;   // Collection
	//transform hit wire and time into cm
	double wire_cm = 0.; 
	if(geom->Cryostat(cstat).TPC(tpc).Plane(plane).SignalType() == geo::kInduction)
	  wire_cm = (double)((wire+3.95) * wire_pitch);          
	else
	  wire_cm = (double)((wire+1.84) * wire_pitch);
	
	//double time_cm = (double)(time * timepitch);
	double time_cm;
	if(time>tSI) time_cm = (double)( (time-tSI)*timepitch + tSI*driftvelocity_SI*timetick);
	else time_cm = time*driftvelocity_SI*timetick;
	
	wires.push_back(wire_cm);
	times.push_back(time_cm);
	
	the2Dtrack->SetPoint(np,wire_cm,time_cm);
	np++;
      }//end of loop over cluster hits
      
      // fit the 2Dtrack and get some info to store
      try{
	the2Dtrack->Fit("pol1","Q");
      }
      catch(...){
	std::cout<<"The 2D track fit failed"<<std::endl;
	continue;
      }
      
      TF1 *pol1=(TF1*) the2Dtrack->GetFunction("pol1");
      double par[2];
      pol1->GetParameters(par);
      double intercept = par[0];
      double slope = par[1];
      
  
      double w0 = wires.front();      // first hit wire (cm)
      double w1 = wires.back();        // last hit wire (cm)
      double t0 = times.front();      // first hit time (cm)
      double t1 = times.back();        // last hit time (cm)
      double t0_line = intercept + (w0)*slope;// time coordinate at wire w0 on the fit line (cm)  
      double t1_line = intercept + (w1)*slope;// time coordinate at wire w1 on the fit line (cm)
    


      // actually store the 2Dtrack info
      switch(geom->Cryostat(cstat).TPC(tpc).Plane(plane).SignalType()){
         case geo::kInduction:
            Iwirefirsts.push_back(w0);
            Iwirelasts.push_back(w1);
            Itimefirsts.push_back(t0);
            Itimelasts.push_back(t1); 
            Itimefirsts_line.push_back(t0_line);
            Itimelasts_line.push_back(t1_line);    
            IclusHitlists.push_back(hitlist);
            Icluster_count.push_back(ii);
            break;
         case geo::kCollection:
            Cwirefirsts.push_back(w0);
            Cwirelasts.push_back(w1);
            Ctimefirsts.push_back(t0);
            Ctimelasts.push_back(t1);
            Ctimefirsts_line.push_back(t0_line);
            Ctimelasts_line.push_back(t1_line);
            CclusHitlists.push_back(hitlist);
            Ccluster_count.push_back(ii);
            break;   
         case geo::kMysteryType:
	    break;
      }
      delete pol1;
   }// end of loop over all input clusters
  
   /////////////////////////////////////////////////////
   /////// 2D Track Matching and 3D Track Reconstruction
   /////////////////////////////////////////////////////

   for(unsigned int collectionIter=0; collectionIter < CclusHitlists.size();collectionIter++){  //loop over Collection view 2D tracks
      // Recover previously stored info
      double Cw0 = Cwirefirsts[collectionIter];
      double Cw1 = Cwirelasts[collectionIter];
      //double Ct0 = Ctimefirsts[collectionIter];
      //double Ct1 = Ctimelasts[collectionIter];
      double Ct0_line = Ctimefirsts_line[collectionIter];
      double Ct1_line = Ctimelasts_line[collectionIter];
      std::vector< art::Ptr<recob::Hit> > hitsCtrk = CclusHitlists[collectionIter];

      double collLength = TMath::Sqrt( TMath::Power(Ct1_line - Ct0_line,2) + TMath::Power(Cw1 - Cw0,2));

      for(unsigned int inductionIter=0;inductionIter<IclusHitlists.size();inductionIter++){   //loop over Induction view 2D tracks
         // Recover previously stored info
         double Iw0 = Iwirefirsts[inductionIter];
         double Iw1 = Iwirelasts[inductionIter];
         //double It0 = Itimefirsts[inductionIter];
         //double It1 = Itimelasts[inductionIter];
         double It0_line = Itimefirsts_line[inductionIter];
         double It1_line = Itimelasts_line[inductionIter];
         std::vector< art::Ptr<recob::Hit> > hitsItrk = IclusHitlists[inductionIter];

         double indLength = TMath::Sqrt( TMath::Power(It1_line - It0_line,2) + TMath::Power(Iw1 - Iw0,2)); 

         bool forward_match = ((fabs(Ct0_line-It0_line)<ftmatch*timepitch) && (fabs(Ct1_line-It1_line)<ftmatch*timepitch));
         bool backward_match = ((fabs(Ct0_line-It1_line)<ftmatch*timepitch) && (fabs(Ct1_line-It0_line)<ftmatch*timepitch));
	 

         if(forward_match || backward_match ){ 	

            // Reconstruct the 3D track
            TVector3 XYZ0, XYZ1;  // track endpoints
            if(forward_match){
               XYZ0.SetXYZ(Ct0_line,(Cw0-Iw0)/(2.*TMath::Sin(Angle)),(Cw0+Iw0)/(2.*TMath::Cos(Angle))-YC/2.*TMath::Tan(Angle));
               XYZ1.SetXYZ(Ct1_line,(Cw1-Iw1)/(2.*TMath::Sin(Angle)),(Cw1+Iw1)/(2.*TMath::Cos(Angle))-YC/2.*TMath::Tan(Angle));
            }
            else{
               XYZ0.SetXYZ(Ct0_line,(Cw0-Iw1)/(2.*TMath::Sin(Angle)),(Cw0+Iw1)/(2.*TMath::Cos(Angle))-YC/2.*TMath::Tan(Angle));
               XYZ1.SetXYZ(Ct1_line,(Cw1-Iw0)/(2.*TMath::Sin(Angle)),(Cw1+Iw0)/(2.*TMath::Cos(Angle))-YC/2.*TMath::Tan(Angle));
            }

            //compute track direction in Local co-ordinate system
            //WARNING:  There is an ambiguity introduced here for the case of backwards-going tracks.  
            //If available, vertex info. could sort this out.
            TVector3 startpointVec,endpointVec;
            TVector2 collVtx, indVtx;
            if(XYZ0.Z() <= XYZ1.Z()){
                startpointVec.SetXYZ(XYZ0.X(),XYZ0.Y(),XYZ0.Z());
                endpointVec.SetXYZ(XYZ1.X(),XYZ1.Y(),XYZ1.Z());
                if(forward_match){
                    collVtx.Set(Ct0_line,Cw0);
                    indVtx.Set(It0_line,Iw0);
                }else{
                    collVtx.Set(Ct0_line,Cw0);
                    indVtx.Set(It1_line,Iw1);
                }
            }
            else{
                startpointVec.SetXYZ(XYZ1.X(),XYZ1.Y(),XYZ1.Z());
                endpointVec.SetXYZ(XYZ0.X(),XYZ0.Y(),XYZ0.Z());
                if(forward_match){
                    collVtx.Set(Ct1_line,Cw1);
                    indVtx.Set(It1_line,Iw1);
                }else{
                    collVtx.Set(Ct1_line,Cw1);
                    indVtx.Set(It0_line,Iw0);
                }
            }

            //compute track (normalized) cosine directions in the TPC co-ordinate system
            TVector3 DirCos = endpointVec - startpointVec;
            
            //SetMag casues a crash if the magnitude of the vector is zero
            try
            {
            DirCos.SetMag(1.0);//normalize vector
            }
            catch(...){std::cout<<"The Spacepoint is infinitely small"<<std::endl;
            continue;
            }

            art::Ptr <recob::Cluster> cl1(clusterListHandle,Icluster_count[inductionIter]);
            art::Ptr <recob::Cluster> cl2(clusterListHandle,Ccluster_count[collectionIter]);
            art::PtrVector<recob::Cluster> clustersPerTrack;
            clustersPerTrack.push_back(cl1);
            clustersPerTrack.push_back(cl2);

 
            /////////////////////////////
            // Match hits
            ////////////////////////////

            //create collection of spacepoints that will be used when creating the Track object
            std::vector<recob::SpacePoint> spacepoints;
	

            std::vector< art::Ptr<recob::Hit> > minhits = hitsCtrk.size() <= hitsItrk.size() ? hitsCtrk : hitsItrk;
	    std::vector< art::Ptr<recob::Hit> > maxhits = hitsItrk.size() <= hitsCtrk.size() ? hitsCtrk : hitsItrk;


            bool maxhitsMatch[maxhits.size()];
            for(unsigned int it=0;it<maxhits.size();it++) maxhitsMatch[it] = false;

            std::vector<recob::Hit*> hits3Dmatched;
            // For the matching start from the view where the track projection presents less hits
            unsigned int imaximum = 0;
            for(unsigned int imin=0;imin<minhits.size();imin++){ //loop over hits
               //get wire - time coordinate of the hit
	      unsigned int channel,wire,plane1,plane2,tpc,cstat;
               channel = minhits[imin]->Channel();
               geom->ChannelToWire(channel,cstat,tpc,plane1,wire);
               // get the wire-time co-ordinates of the hit to be matched
               //double w1 = (double)((wire+1)*wire_pitch);
               double w1=0;
               
               //the 3.95 and 1.84 below are the ArgoNeuT TPC offsets for the induction and collection plane, respectively and are in units of wire pitch.
               if(geom->Cryostat(cstat).TPC(tpc).Plane(plane1).SignalType() == geo::kInduction)
               w1 = (double)((wire+3.95) * wire_pitch);          
               else
               w1 = (double)((wire+1.84) * wire_pitch);
               
               double temptime1 = minhits[imin]->PeakTime()-presamplings;
               if(geom->Cryostat(cstat).TPC(tpc).Plane(plane1).SignalType() == geo::kCollection) temptime1 -= tIC;
               double t1;// = plane1==1?(double)((minhits[imin]->PeakTime()-presamplings-tIC)*timepitch):(double)((minhits[imin]->PeakTime()-presamplings)*timepitch); //in cm
               if(temptime1>tSI) t1 = (double)( (temptime1-tSI)*timepitch + tSI*driftvelocity_SI*timetick);
               else t1 = temptime1*driftvelocity_SI*timetick;

               //get the track origin co-ordinates in the two views
               TVector2 minVtx2D;
               (geom->Plane(plane1,tpc).SignalType() == geo::kCollection) ? minVtx2D.Set(collVtx.X(),collVtx.Y()): minVtx2D.Set(indVtx.X(),indVtx.Y());
               TVector2 maxVtx2D;
               (geom->Plane(plane1,tpc).SignalType() == geo::kCollection) ? maxVtx2D.Set(indVtx.X(),indVtx.Y()): maxVtx2D.Set(collVtx.X(),collVtx.Y());
               
               double ratio = (collLength>indLength) ? collLength/indLength : indLength/collLength;	  

               //compute the distance of the hit (imin) from the relative track origin
               double minDistance = ratio*TMath::Sqrt(TMath::Power(t1-minVtx2D.X(),2) + TMath::Power(w1-minVtx2D.Y(),2));
	  
	  
               //core matching algorithm
               double difference = 9999999.;	  

               for(unsigned int imax = 0; imax < maxhits.size(); imax++){ //loop over hits of the other view
                  if(!maxhitsMatch[imax]){
                     //get wire - time coordinate of the hit
                     channel = maxhits[imax]->Channel();
                     geom->ChannelToWire(channel,cstat,tpc,plane2,wire);
                     //double w2 = (double)((wire+1)*wire_pitch);
                     double w2=0.;
                     if(geom->Cryostat(cstat).TPC(tpc).Plane(plane2).SignalType() == geo::kInduction)
                     w2 = (double)((wire+3.95) * wire_pitch);          
                     else
                     w2 = (double)((wire+1.84) * wire_pitch);
                     
                     double temptime2 = maxhits[imax]->PeakTime()-presamplings;
                     if(geom->Cryostat(cstat).TPC(tpc).Plane(plane2).SignalType() == geo::kCollection) temptime2 -= tIC;
                     double t2;
                     if(temptime2>tSI) t2 = (double)( (temptime2-tSI)*timepitch + tSI*driftvelocity_SI*timetick);
                     else t2 = temptime2*driftvelocity_SI*timetick;
                   
                     bool timematch = (fabs(t1-t2)<ftmatch*timepitch);
                     bool wirematch = (fabs(w1-w2)<wireShift*wire_pitch);
	      
                     double maxDistance = TMath::Sqrt(TMath::Power(t2-maxVtx2D.X(),2)+TMath::Power(w2-maxVtx2D.Y(),2));
                     if (wirematch && timematch && fabs(maxDistance-minDistance)<difference) {
                        difference = fabs(maxDistance-minDistance);
                        imaximum = imax;
                     }
                  }
               }
               maxhitsMatch[imaximum]=true;
	  
               art::PtrVector<recob::Hit> sp_hits;
               if(difference!= 9999999.){
                  sp_hits.push_back(minhits[imin]);
                  sp_hits.push_back(maxhits[imaximum]);
               }
	  
               // Get the time-wire co-ordinates of the matched hit
               channel =  maxhits[imaximum]->Channel();
               geom->ChannelToWire(channel,cstat,tpc,plane2,wire);
             
               //double w1_match = (double)((wire+1)*wire_pitch);  
               double w1_match=0.;
               if(geom->Cryostat(cstat).TPC(tpc).Plane(plane2).SignalType() == geo::kInduction)
               w1_match = (double)((wire+3.95) * wire_pitch);          
               else
               w1_match = (double)((wire+1.84) * wire_pitch);
               
               double temptime3 = maxhits[imaximum]->PeakTime()-presamplings;
               if(geom->Cryostat(cstat).TPC(tpc).Plane(plane2).SignalType() == geo::kCollection) temptime3 -= tIC;
               double t1_match;
               if(temptime3>tSI) t1_match = (double)( (temptime3-tSI)*timepitch + tSI*driftvelocity_SI*timetick);
               else t1_match = temptime3*driftvelocity_SI*timetick;
             
               // create the 3D hit, compute its co-ordinates and add it to the 3D hits list	  
               double Ct = geom->Cryostat(cstat).TPC(tpc).Plane(plane1).SignalType()==geo::kCollection?t1:t1_match;
               double Cw = geom->Cryostat(cstat).TPC(tpc).Plane(plane1).SignalType()==geo::kCollection?w1:w1_match;
               double Iw = geom->Cryostat(cstat).TPC(tpc).Plane(plane1).SignalType()==geo::kCollection?w1_match:w1;

               const TVector3 hit3d(Ct,(Cw-Iw)/(2.*TMath::Sin(Angle)),(Cw+Iw)/(2.*TMath::Cos(Angle))-YC/2.*TMath::Tan(Angle)); 
               
               
               Double_t hitcoord[3];       
               hitcoord[0] = hit3d.X();
               hitcoord[1] = hit3d.Y();
               hitcoord[2] = hit3d.Z();

               /*
               double yy,zz;
               if(geom->ChannelsIntersect(geom->PlaneWireToChannel(0,(int)((Iw/wire_pitch)-3.95)),      geom->PlaneWireToChannel(1,(int)((Cw/wire_pitch)-1.84)),yy,zz))
               {
               //channelsintersect provides a slightly more accurate set of y and z coordinates. use channelsintersect in case the wires in question do cross.
               hitcoord[1] = yy;
               hitcoord[2] = zz;               
               mf::LogInfo("SpacePts: ") << "SpacePoint adding xyz ..." << hitcoord[0] <<","<< hitcoord[1] <<","<< hitcoord[2];	       
// 	           std::cout<<"wire 1: "<<(Iw/wire_pitch)-3.95<<" "<<(Cw/wire_pitch)-1.84<<std::endl;
//                std::cout<<"Intersect: "<<yy<<" "<<zz<<std::endl;
	           }
	           else
	           continue;
               */
                   
	       double err[6] = {util::kBogusD};
               recob::SpacePoint mysp(hitcoord, err, spacepoints.size());//3d point at end of track
	       // Don't add a spacepoint right on top of the last one.
	       const double eps(0.1); // 1mm
	       if (spacepoints.size()>=1){
		 TVector3 magNew(mysp.XYZ()[0],mysp.XYZ()[1],mysp.XYZ()[2]);
		 TVector3 magLast(spacepoints[spacepoints.size()-1].XYZ()[0],
				  spacepoints[spacepoints.size()-1].XYZ()[1],
				  spacepoints[spacepoints.size()-1].XYZ()[2]);
		 if (!(magNew.Mag()>=magLast.Mag()+eps || 
		       magNew.Mag()<=magLast.Mag()-eps) )
		   continue;
	       }
	       spacepoints.push_back(mysp);	  

            }//loop over min-hits
      

            // Add the 3D track to the vector of the reconstructed tracks
            if(spacepoints.size()>0){

	      // make a vector of the trajectory points along the track
	      std::vector<TVector3> xyz(spacepoints.size());
	      size_t spStart = spcol->size();
	      for(size_t s = 0; s < spacepoints.size(); ++s){
		xyz[s] = TVector3(spacepoints[s].XYZ());
		spcol->push_back(spacepoints[s]);		
	      }
	      size_t spEnd = spcol->size();
		
	      ///\todo really should fill the direction cosines with unique values 
	      std::vector<TVector3> dircos(spacepoints.size(), DirCos);

	      std::vector< std::vector<double> > dQdx;
	      std::vector<double> mom(2, util::kBogusD);
	      recob::Track  the3DTrack(xyz, dircos, dQdx, mom, tcol->size());
	      tcol->push_back(the3DTrack);

	      // make associations between the track and space points
	      util::CreateAssn(*this, evt, *tcol, *spcol, *tspassn, spStart, spEnd);

	      // now the track and clusters
	      util::CreateAssn(*this, evt, *tcol, clustersPerTrack, *tcassn);

	      // and the hits and track
	      art::FindManyP<recob::Hit> fmh(clustersPerTrack, evt, fClusterModuleLabel);
	      for(size_t cpt = 0; cpt < clustersPerTrack.size(); ++cpt)
		util::CreateAssn(*this, evt, *tcol, fmh.at(cpt), *thassn);

            }
         } //close match 2D tracks

      }//close loop over Induction view 2D tracks
    
   }//close loop over Collection xxview 2D tracks


   ///\todo the tracks made in this module should really be associated to Hits, Clusters and SpacePoints
   // but don't do it now as it isn't clear anyone is actually using the module

   mf::LogVerbatim("Summary") << std::setfill('-') << std::setw(175) << "-" << std::setfill(' ');
   mf::LogVerbatim("Summary") << "SpacePts Summary:";
   for(unsigned int i = 0; i<tcol->size(); ++i) mf::LogVerbatim("Summary") << tcol->at(i) ;
 
   evt.put(tcol);
   evt.put(spcol);
   evt.put(tspassn);
   evt.put(tcassn);
   evt.put(thassn);

}
