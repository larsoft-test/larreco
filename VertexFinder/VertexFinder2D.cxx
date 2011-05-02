////////////////////////////////////////////////////////////////////////
//
// VertexFinder2D class
//
// tjyang@fnal.gov
//
// This algorithm is designed to reconstruct the vertices using the
// 2D cluster information
// 
// This is Preliminary Work and needs modifications
// ////////////////////////////////////////////////////////////////////////


#include <iostream>

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
#include "VertexFinder/VertexFinder2D.h"
#include <iomanip>
#include <ios>
#include <sstream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <TVector3.h>
#include <vector>



#include "RawData/RawDigit.h"
#include "Filters/ChannelFilter.h"
#include "SimulationBase/simbase.h"
#include "RecoBase/recobase.h"
#include "Geometry/geo.h"
#include "Utilities/LArProperties.h"
#include "TH1.h"
#include "TH2.h"
#include "TVectorD.h"
#include "TGeoManager.h"
#include "TMath.h"
#include "TGraph.h"
#include "TF1.h"

namespace vertex{

//-----------------------------------------------------------------------------
  VertexFinder2D::VertexFinder2D(fhicl::ParameterSet const& pset)
  {  
    this->reconfigure(pset);    
    produces< std::vector<recob::Vertex> >();
  }
//-----------------------------------------------------------------------------
  VertexFinder2D::~VertexFinder2D()
  {
  }

  //---------------------------------------------------------------------------
  void VertexFinder2D::reconfigure(fhicl::ParameterSet p) 
  {
    fClusterModuleLabel  = p.get< std::string >("ClusterModuleLabel");
    return;
  }
  //-------------------------------------------------------------------------
  void VertexFinder2D::beginJob(){
    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;
  }

// //-----------------------------------------------------------------------------
  void VertexFinder2D::produce(art::Event& evt)
  {

    
    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<util::LArProperties> larprop;

    // define TPC parameters
    TString tpcName = geom->GetLArTPCVolumeName();
    
    //  double YC =  (m_TPCHalfZ-5.)*2.; // TPC height in cm
    double YC =  (geom->DetHalfHeight()-5.0)*2.; // *ArgoNeuT* TPC active-volume height in cm
    double Angle = geom->Plane(1).Wire(0).ThetaZ(false)-TMath::Pi()/2.; // wire angle with respect to the vertical direction
    // Parameters temporary defined here, but possibly to be retrieved somewhere in the code
    double timetick = 0.198;    //time sample in us
    double presamplings = 60.;
    //const double wireShift=50.; // half the number of wires from the Induction(Collection) plane intersecting with a wire from the Collection(Induction) plane.
    double plane_pitch = geom->PlanePitch(0,1);   //wire plane pitch in cm 
    double wire_pitch = geom->WirePitch(0,1,0);    //wire pitch in cm
    double Efield_drift = 0.5;  // Electric Field in the drift region in kV/cm
    double Efield_SI = 0.7;     // Electric Field between Shield and Induction planes in kV/cm
    double Efield_IC = 0.9;     // Electric Field between Induction and Collection planes in kV/cm
    double Temperature = 87.6;  // LAr Temperature in K
    
    double driftvelocity = larprop->DriftVelocity(Efield_drift,Temperature); //drift velocity in the drift 
    //region (cm/us)
    double driftvelocity_SI = larprop->DriftVelocity(Efield_SI,Temperature); //drift velocity between shield 
    //and induction (cm/us)
    double driftvelocity_IC = larprop->DriftVelocity(Efield_IC,Temperature); //drift velocity between induction 
    //and collection (cm/us)
    double timepitch = driftvelocity*timetick;                               //time sample (cm) 
    double tSI = plane_pitch/driftvelocity_SI/timetick;                      //drift time between Shield and 
    //Collection planes (time samples)
    double tIC = plane_pitch/driftvelocity_IC/timetick;                      //drift time between Induction and 
    //Collection planes (time samples)
    
    
    art::Handle< std::vector<recob::Cluster> > clusterListHandle;
    evt.getByLabel(fClusterModuleLabel,clusterListHandle);

    art::PtrVector<recob::Cluster> clusters;
    for (unsigned int ii = 0; ii <  clusterListHandle->size(); ++ii){
      art::Ptr<recob::Cluster> clusterHolder(clusterListHandle,ii);
      clusters.push_back(clusterHolder);
    }


    //Point to a collection of vertices to output.
    std::auto_ptr<std::vector<recob::Vertex> > vcol(new std::vector<recob::Vertex>);

    int nplanes = geom->Nplanes();
    
    std::vector<int> Cls[nplanes]; //index to clusters in each view
    std::vector<double> dtdwstart;

    //loop over clusters
    for(unsigned int iclu=0; iclu<clusters.size();++iclu){

      switch(clusters[iclu]->View()){

      case geo::kU :
	Cls[0].push_back(iclu);
	break;
      case geo::kV :
	Cls[1].push_back(iclu);
	break;
      case geo::kW :
	Cls[2].push_back(iclu);
	break;
      default :
	break;
      }

      unsigned int channel;
      unsigned int wire;
      unsigned int plane;
      double time;
      
      std::vector<double> wires;
      std::vector<double> times;
      
      art::PtrVector<recob::Hit> hit;
      hit = clusters[iclu]->Hits();
      int n = 0;
      for (unsigned int i = 0; i<hit.size(); i++){
	time = hit[i]->PeakTime();
	channel = hit[i]->Channel();
	geom->ChannelToWire(channel,plane,wire);
	wires.push_back(wire);
	times.push_back(time);
	n++;
	//if (n==10) break;
      }
      if (n>=2){
	TGraph *the2Dtrack = new TGraph(n,&wires[0],&times[0]);
	the2Dtrack->Fit("pol1","Q");
	TF1 *pol1=(TF1*) the2Dtrack->GetFunction("pol1");
	double par[2];
	pol1->GetParameters(par);
	//std::cout<<iclu<<" "<<par[1]<<" "<<clusters[iclu]->dTdW()<<std::endl;
	dtdwstart.push_back(par[1]);
	delete the2Dtrack;
      }
      else dtdwstart.push_back(clusters[iclu]->dTdW());
    }
    
    std::vector<int> cluvtx[nplanes];
    std::vector<double> vtx_w;
    std::vector<double> vtx_t;
    
    for (int i = 0; i<nplanes; i++){
      if (Cls[i].size()>=1){
	//find the longest two clusters
	int c1 = -1;
	int c2 = -1;
	double ww0 = -999;
	double wb1 = -999;
	double we1 = -999;
	double wb2 = -999;
	double we2 = -999;
	double tt0 = -999;
	double tt1 = -999;
	double tt2 = -999;
	double dtdw1 = -999;
	double dtdw2 = -999;
	double lclu1 = -999;
	double lclu2 = -999;
	for (unsigned j = 0; j<Cls[i].size(); j++){
	  double lclu = sqrt(pow((clusters[Cls[i][j]]->StartPos()[0]-clusters[Cls[i][j]]->EndPos()[0])*13.5,2)+pow(clusters[Cls[i][j]]->StartPos()[1]-clusters[Cls[i][j]]->EndPos()[1],2));
	  bool rev = false;
	  if (c1!=-1){
	    double wb = clusters[Cls[i][j]]->StartPos()[0];
	    //double we = clusters[Cls[i][j]]->EndPos()[0];
	    double tt = clusters[Cls[i][j]]->StartPos()[1];
	    double dtdw = dtdwstart[Cls[i][j]];
	    tt0 = (dtdw1*dtdw*(wb1-wb)+dtdw1*tt-dtdw*tt1)/(dtdw1-dtdw);
	    ww0 = (tt-tt1+dtdw1*wb1-dtdw*wb)/(dtdw1-dtdw);
	    if (fabs(wb1-ww0)>fabs(we1-ww0)) rev = true;//reverse cluster dir
	  }
	  if (lclu1<lclu){
	    //this is to remove delta rays from the muon cluster
	    if (c1!=-1&&((!rev&&ww0<wb1+15)||(rev&&ww0>we1-15))){
	      lclu2 = lclu1;
	      c2 = c1;
	      wb2 = wb1;
	      we2 = we1;
	      tt2 = tt1;
	      dtdw2 = dtdw1;
	    }
	    lclu1 = lclu;
	    c1 = Cls[i][j];
	    wb1 = clusters[Cls[i][j]]->StartPos()[0];
	    we1 = clusters[Cls[i][j]]->EndPos()[0];
	    tt1 = clusters[Cls[i][j]]->StartPos()[1];
	    dtdw1 = dtdwstart[Cls[i][j]];
	  }
	  else if (lclu2<lclu){
	    if (ww0<wb1+15){
	      lclu2 = lclu;
	      c2 = Cls[i][j];
	    }
	  }
	}
	if (c1!=-1&&c2!=-1){
	  cluvtx[i].push_back(c1);
	  cluvtx[i].push_back(c2);
	  //std::cout<<i<<" "<<c1<<" "<<c2<<std::endl;
	  double w1 = clusters[c1]->StartPos()[0];
	  double t1 = clusters[c1]->StartPos()[1];
	  double k1 = dtdwstart[c1];
	  double w2 = clusters[c2]->StartPos()[0];
	  double t2 = clusters[c2]->StartPos()[1];
	  double k2 = dtdwstart[c2];
	  //std::cout<<k1<<" "<<k2<<" "<<k1-k2<<std::endl;
	  //calculate the vertex
	  if (fabs(k1-k2)<1e-10) continue;
	  double t0 = (k1*k2*(w1-w2)+k1*t2-k2*t1)/(k1-k2);
	  double w0 = (t2-t1+k1*w1-k2*w2)/(k1-k2);
	  vtx_w.push_back(w0);
	  vtx_t.push_back(t0);
	}
	else if (Cls[i].size()>=1){
	  cluvtx[i].push_back(Cls[i][0]);
	  vtx_w.push_back(clusters[Cls[i][0]]->StartPos()[0]);
	  vtx_t.push_back(clusters[Cls[i][0]]->StartPos()[1]);
	}
      }
      else {//no cluster found
	vtx_w.push_back(-1);
	vtx_t.push_back(-1);
      }
    }

    Double_t vtxcoord[3];
    if (Cls[0].size()>0&&Cls[1].size()>0){//ignore w view
      double Iw0 = (vtx_w[0]+1)*wire_pitch;
      double Cw0 = (vtx_w[1]+1)*wire_pitch;
      double It0 = vtx_t[0] - presamplings;
      It0 -= tSI;
      It0 *= timepitch;
      double Ct0 = vtx_t[1] - presamplings;
      Ct0 -= tSI+tIC;
      Ct0 *= timepitch;
      vtxcoord[0] = Ct0;
      vtxcoord[1] = (Cw0-Iw0)/(2.*TMath::Sin(Angle));
      vtxcoord[2] = (Cw0+Iw0)/(2.*TMath::Cos(Angle))-YC/2.*TMath::Tan(Angle);
    }
    else{
      vtxcoord[0] = -99999;
      vtxcoord[1] = -99999;
      vtxcoord[2] = -99999;
    }
    
    //need to implement this
    art::PtrVector<recob::Track> vTracks_vec;
    art::PtrVector<recob::Shower> vShowers_vec;
    
    recob::Vertex the3Dvertex(vTracks_vec, vShowers_vec, vtxcoord);
    vcol->push_back(the3Dvertex);

    evt.put(vcol);
    
  } // end of produce
} // end of vertex namespace

// //-----------------------------------------------------------------------------
