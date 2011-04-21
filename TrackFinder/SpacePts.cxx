////////////////////////////////////////////////////////////////////////
//
// \file SpacePts.cxx
//
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
#include "SpacePts.h"
#include "Geometry/geo.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/recobase.h"
#include "Utilities/LArProperties.h"

// ROOT includes
#include "TVectorD.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMath.h"

namespace trkf {
struct SortByWire 
{
  bool operator() (recob::Hit const& h1, recob::Hit const& h2) const 
  { return 
      h1.Wire()->RawDigit()->Channel() < 
      h2.Wire()->RawDigit()->Channel() ;
  }
};
} // end namespace

//-------------------------------------------------
trkf::SpacePts::SpacePts(fhicl::ParameterSet const& pset)
{
  this->reconfigure(pset);

   produces< std::vector<recob::Track> >();
}

//-------------------------------------------------
trkf::SpacePts::~SpacePts()
{
}

void trkf::SpacePts::reconfigure(fhicl::ParameterSet pset) 
{
  fPreSamplings           = pset.get< double >("TicksOffset");
  ftmatch                 = pset.get< int    >("TMatch");
  fClusterModuleLabel     = pset.get< std::string >("ClusterModuleLabel");

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
   std::auto_ptr<std::vector<recob::Track> > tcol(new std::vector<recob::Track>);

   // define TPC parameters
   TString tpcName = geom->GetLArTPCVolumeName();

   //TPC dimensions
   //double m_TPCHalfZ = m_tpcVolumeUtility->GetHalfZ();
   double m_TPCHalfZ = geom->DetLength()-5.0;

   //  double YC =  (m_TPCHalfZ-5.)*2.; // TPC height in cm
   double YC =  (geom->DetHalfHeight()-0.5715)*2.; // TPC height in cm
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
   double Temperature = 87.6;  // LAr Temperature in K

   double driftvelocity = larprop->DriftVelocity(Efield_drift,Temperature);    //drift velocity in the drift region (cm/us)
   double driftvelocity_SI = larprop->DriftVelocity(Efield_SI,Temperature);    //drift velocity between shield and induction (cm/us)
   double driftvelocity_IC = larprop->DriftVelocity(Efield_IC,Temperature);    //drift velocity between induction and collection (cm/us)
   double timepitch = driftvelocity*timetick;                         //time sample (cm) 
   double tSI = plane_pitch/driftvelocity_SI/timetick;                   //drift time between Shield and Collection planes (time samples)
   double tIC = plane_pitch/driftvelocity_IC/timetick;                //drift time between Induction and Collection planes (time samples)


   // get input Cluster object(s).
   art::Handle< std::vector<recob::Cluster> > clusterListHandle;
   evt.getByLabel(fClusterModuleLabel,clusterListHandle);

  
   // Declare some vectors..
   // Induction
   std::vector<double> Iwirefirsts;       // in cm
   std::vector<double> Iwirelasts;        // in cm
   std::vector<double> Itimefirsts;       // in cm
   std::vector<double> Itimelasts;        // in cm
   std::vector<double> Itimefirsts_line;  // in cm
   std::vector<double> Itimelasts_line;   // in cm
   std::vector < art::PtrVector<recob::Hit> > IclusHitlists;
   std::vector<unsigned int> Icluster_count; 

   // Collection
   std::vector<double> Cwirefirsts;       // in cm
   std::vector<double> Cwirelasts;        // in cm
   std::vector<double> Ctimefirsts;       // in cm
   std::vector<double> Ctimelasts;        // in cm
   std::vector<double> Ctimefirsts_line;  // in cm
   std::vector<double> Ctimelasts_line;   // in cm
   std::vector< art::PtrVector < recob::Hit> > CclusHitlists;
   std::vector<unsigned int> Ccluster_count; 

   for(unsigned int ii = 0; ii < clusterListHandle->size(); ++ii)
   {

      art::Ptr<recob::Cluster> cl(clusterListHandle, ii);
      
      // Figure out which View the cluster belongs to 
      int clPlane = cl->View()-1;
      // Gaaaaaah! Change me soon!!! But, for now, 
      // let's just chuck one plane's worth of info. EC, 30-Mar-2011.
      if (cl->View() == geo::kW) continue; //kW, you'd think.

      //      std::cout << "SpacePts: Plane/View/SignalType ..." << clPlane <<" "<<cl->View()<<" "<<geom->Plane(clPlane).SignalType() << std::endl;
      // Some variables for the hit
      unsigned int channel;  //channel number
      float time;            //hit time at maximum
      unsigned int wire;              //hit wire number
      unsigned int plane;             //hit plane number
      
      
      art::PtrVector<recob::Hit> hitlist;
      hitlist = cl->Hits();
      hitlist.sort(trkf::SortByWire());
      

      
      std::vector<double> wires;
      std::vector<double> times;
     
      
      int np=0;
      for(art::PtrVectorItr<recob::Hit> theHit = hitlist.begin(); theHit != hitlist.end();  theHit++) //loop over cluster hits
      {
         //recover the Hit
         //      recob::Hit* theHit = (recob::Hit*)(*hitIter);
         time = (*theHit)->PeakTime() ;
	 
         time -= presamplings;
	  
         channel = (*theHit)->Channel();
         geom->ChannelToWire(channel,plane,wire);
	  

         //correct for the distance between wire planes
         if(plane==0) time -= tSI;         // Induction
         if(plane==1) time -= (tSI+tIC);   // Collection
	
         //transform hit wire and time into cm
         double wire_cm = (double)((wire+1) * wire_pitch); 
         double time_cm = (double)(time * timepitch);
         wires.push_back(wire_cm);
         times.push_back(time_cm);
	
         np++;
      }//end of loop over cluster hits
    
      
  
      double w0 = wires.front();      // first hit wire (cm)
      double w1 = wires.back();        // last hit wire (cm)
      double t0 = times.front();      // first hit time (cm)
      double t1 = times.back();        // last hit time (cm)
    


      // actually store the 2Dtrack info
      switch(plane){
         case 0:
            Iwirefirsts.push_back(w0);
            Iwirelasts.push_back(w1);
            Itimefirsts.push_back(t0);
            Itimelasts.push_back(t1); 
            IclusHitlists.push_back(hitlist);
            Icluster_count.push_back(ii);
            break;
         case 1:
            Cwirefirsts.push_back(w0);
            Cwirelasts.push_back(w1);
            Ctimefirsts.push_back(t0);
            Ctimelasts.push_back(t1);
            CclusHitlists.push_back(hitlist);
            Ccluster_count.push_back(ii);
            break;   
      }


   }// end of loop over all input clusters
  
   /////////////////////////////////////////////////////
   /////// 2D Track Matching and 3D Track Reconstruction
   /////////////////////////////////////////////////////

   for(unsigned int collectionIter=0; collectionIter < CclusHitlists.size();collectionIter++){  //loop over Collection view 2D tracks
      // Recover previously stored info
      double Cw0 = Cwirefirsts[collectionIter];
      double Cw1 = Cwirelasts[collectionIter];
      double Ct0 = Ctimefirsts[collectionIter];
      double Ct1 = Ctimelasts[collectionIter];
      art::PtrVector<recob::Hit> hitsCtrk = CclusHitlists[collectionIter];

      for(unsigned int inductionIter=0;inductionIter<IclusHitlists.size();inductionIter++){   //loop over Induction view 2D tracks
         // Recover previously stored info
         double Iw0 = Iwirefirsts[inductionIter];
         double Iw1 = Iwirelasts[inductionIter];
         double It0 = Itimefirsts[inductionIter];
         double It1 = Itimelasts[inductionIter];
         art::PtrVector<recob::Hit> hitsItrk = IclusHitlists[inductionIter];

	 // My 1000s below. EC, 11-Apr-2011
         if((fabs(Ct0-It0)<ftmatch*timepitch) && (fabs(Ct1-It1)<ftmatch*timepitch)){ 
            //std::cout<<"-----> Track "<<collectionIter<< " Collection associated with track "<<inductionIter<< " Induction"<<std::endl;
	
       
            // Reconstruct the 3D track
            TVector3 XYZ0;  // track origin or interaction vertex
            XYZ0.SetXYZ(Ct0,(Cw0-Iw0)/(2.*TMath::Sin(Angle)),(Cw0+Iw0)/(2.*TMath::Cos(Angle))-YC/2.*TMath::Tan(Angle));

            //compute track startpoint and endpoint in Local co-ordinate system 
            TVector3 startpointVec(XYZ0.X(),XYZ0.Y(),XYZ0.Z());
            TVector3 endpointVec(Ct1,(Cw1-Iw1)/(2.*TMath::Sin(Angle)),(Cw1+Iw1)/(2.*TMath::Cos(Angle))-YC/2.*TMath::Tan(Angle));

            //compute track (normalized) cosine directions in the TPC co-ordinate system
            TVector3 DirCos = endpointVec - startpointVec;
            DirCos.SetMag(1.0);//normalize vector

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
	

            art::PtrVector<recob::Hit> minhits = hitsCtrk.size() <= hitsItrk.size() ? hitsCtrk : hitsItrk;
            art::PtrVector<recob::Hit> maxhits = hitsItrk.size() <= hitsCtrk.size() ? hitsCtrk : hitsItrk;


            bool maxhitsMatch[maxhits.size()];
            for(unsigned int it=0;it<maxhits.size();it++) maxhitsMatch[it] = false;

            std::vector<recob::Hit*> hits3Dmatched;
            // For the matching start from the view where the track projection presents less hits
            unsigned int imaximum = 0;
            for(unsigned int imin=0;imin<minhits.size();imin++){ //loop over hits
               //get wire - time coordinate of the hit
               unsigned int channel,wire,plane1,plane2;
               channel = minhits[imin]->Channel();
               geom->ChannelToWire(channel,plane1,wire);
               // get the wire-time co-ordinates of the hit to be matched
               double w1 = (double)((wire+1)*wire_pitch);
               double t1 = plane1==1?(double)((minhits[imin]->PeakTime()-presamplings-(tSI+tIC))*timepitch):(double)((minhits[imin]->PeakTime()-presamplings-tSI)*timepitch); //in cm	  

               //get the track origin co-ordinates in the two views
               TVector2 minVtx2D;
               plane1==1 ? minVtx2D.Set(Ct0,Cw0): minVtx2D.Set(It0,Iw0);
               TVector2 maxVtx2D;
               plane1==1 ? maxVtx2D.Set(It0,Iw0): maxVtx2D.Set(Ct0,Cw0);
	  
               // get the track last hit co-ordinates in the two views
               channel = minhits[minhits.size()-1]->Channel();
               geom->ChannelToWire(channel,plane1,wire);
               double w1_last = plane1==1?Cw1:Iw1;
               double t1_last = plane1==1?Ct1:It1;  
               channel = maxhits[maxhits.size()-1]->Channel();
               geom->ChannelToWire(channel,plane2,wire);
               double w2_last =plane2==1?Cw1:Iw1;
               double t2_last =plane2==1?Ct1:It1;
	  
               // compute the track length in the two views
               double minLength = TMath::Sqrt(TMath::Power(t1_last-minVtx2D.X(),2)+TMath::Power(w1_last-minVtx2D.Y(),2));
               double maxLength = TMath::Sqrt(TMath::Power(t2_last-maxVtx2D.X(),2)+TMath::Power(w2_last-maxVtx2D.Y(),2));
	  
               //compute the distance of the hit (imin) from the relative track origin
               double minDistance = maxLength* TMath::Sqrt(TMath::Power(t1-minVtx2D.X(),2)+TMath::Power(w1-minVtx2D.Y(),2))/minLength;
	  
               //core matching algorithm
               double difference = 9999999.;	  

               for(unsigned int imax = imaximum; imax < maxhits.size(); imax++){ //loop over hits of the other view
                  if(!maxhitsMatch[imax]){
                     //get wire - time coordinate of the hit
                     channel = maxhits[imax]->Channel();
                     geom->ChannelToWire(channel,plane2,wire);
                     double w2 = (double)((wire+1)*wire_pitch);
                     double t2 = plane2==1?(double)((maxhits[imax]->PeakTime()-presamplings-(tSI+tIC))*timepitch):(double)((maxhits[imax]->PeakTime()-presamplings-tSI)*timepitch); //in cm
	      
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
               geom->ChannelToWire(channel,plane2,wire);
               double w1_match = (double)((wire+1)*wire_pitch);
               double t1_match = plane2==1?(double)((maxhits[imaximum]->PeakTime()-presamplings-(tSI+tIC))*timepitch):(double)((maxhits[imaximum]->PeakTime()-presamplings-tSI)*timepitch);

               // create the 3D hit, compute its co-ordinates and add it to the 3D hits list	  
               double Ct = plane1==1?t1:t1_match;
               double Cw = plane1==1?w1:w1_match;
               double Iw = plane1==1?w1_match:w1;

               const TVector3 hit3d(Ct,(Cw-Iw)/(2.*TMath::Sin(Angle)),(Cw+Iw)/(2.*TMath::Cos(Angle))-YC/2.*TMath::Tan(Angle)); 
               Double_t hitcoord[3];       
               hitcoord[0] = hit3d.X();
               hitcoord[1] = hit3d.Y();
               hitcoord[2] = hit3d.Z();         
	       mf::LogInfo("SpacePts: ") << "SpacePoint adding xyz ..." << hitcoord[0] <<","<< hitcoord[1] <<","<< hitcoord[2];
               recob::SpacePoint mysp(sp_hits);//3d point at end of track
               mysp.SetXYZ(hitcoord);
               mysp.SetID(spacepoints.size());
               spacepoints.push_back(mysp);
	  
            }//loop over min-hits
      

            // Add the 3D track to the vector of the reconstructed tracks
            if(spacepoints.size()>0 || clustersPerTrack.size()>0){
	 
               recob::Track  the3DTrack(clustersPerTrack,spacepoints);
               double dircos[3];
               DirCos.GetXYZ(dircos);
               the3DTrack.SetDirection(dircos,dircos);
               the3DTrack.SetID(tcol->size());
               tcol->push_back(the3DTrack);
            }
         } //close match 2D tracks


      }//close loop over Induction view 2D tracks
    
   }//close loop over Collection xxview 2D tracks

   mf::LogVerbatim("Summary") << std::setfill('-') << std::setw(175) << "-" << std::setfill(' ');
   mf::LogVerbatim("Summary") << "SpacePts Summary:";
   for(int i = 0; i<tcol->size(); ++i) mf::LogVerbatim("Summary") << tcol->at(i) ;
 
   evt.put(tcol);
  

}
