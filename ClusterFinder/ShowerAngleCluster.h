////////////////////////////////////////////////////////////////////////
/// \file  ShowerAngleCluster.h
/// \brief
///
///
/// \version $Id: SingleGen.cxx,v 1.4 2010/03/29 09:54:01 brebel Exp $
/// \author:  biagio, ART port: echurch, detector agnostic + branch into cluster module: andrzejs
////////////////////////////////////////////////////////////////////////
#ifndef SHOWERANGLECLUSTER_H
#define SHOWERANGLECLUSTER_H

#include "art/Framework/Core/EDProducer.h" // include the proper bit of the framework

#include <vector>
#include <string>

#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>

#include "TH1F.h"
#include "TTree.h"
#include "TMath.h"

#include "RecoBase/Cluster.h"

namespace cluster {

  class ShowerAngleCluster : public art::EDProducer {
    
  public:

    /**METHODS global*/
    explicit ShowerAngleCluster(fhicl::ParameterSet const& pset);/**Constructor*/
    virtual ~ShowerAngleCluster();                               /**Destructor*/
    void beginJob();                                     
    void reconfigure(fhicl::ParameterSet const& pset);
    void produce(art::Event& evt);                       /**Routine that finds the cluster and sets the dTdW of the 2D shower*/
   
    void   AngularDistribution(art::PtrVector < recob::Hit> hitlist); /**Calculate 2D angle histograms, provided vertex is know */ 
    void   Get2DVariables(art::PtrVector < recob::Hit> hitlist);   /** Calculate 2D variables to be saved into 	  */ 
    void   FitAngularDistributions(int plane);  /** Get actual 2D omega angle from the histograms  */
  

    void   GetVertexN(art::Event& evt); /**Get vertex from MCtruth information. Temporary. */


    
//    float LastWire[2];  // last wire of the shower
//    float LastTime[2];  // last t_hit of the shower
    
   


      
  
    const static float pi            = 3.1415;
    
    const static float ftimetick      =  0.198; // time sample in us
    const static float presamplings  = 60.;
    const static float fdriftvelocity =  0.157;
    const static int    alpha           = 5;     // parameter (how many RMs (of the anglular distribution) is large the cone of the shower)

 private:

  
   int fRun,fEvent,fSubRun;
   

 // 2D slope and intercept of the shower axis
    double slope[3];       // in cm, cm
 //   float intercept[3];   // in cm, cm
 //   float slope_wt[3];       // in wire, tick
 //   float intercept_wt[3];   // in wire,tick

   
//input labels:
 
  std::string fClusterModuleLabel;
  std::string fVertexCLusterModuleLabel;
    std::string fMCGeneratorLabel;
    std::string fLarGeantlabel;
  int fUseMCVertex;

 

  double xyz_vertex[3];


 
  double fMean_wire_pitch ;   // wire pitch in cm

    std::vector<double> fOmega_Mean;    // Mean value of the 2D angular distribution (0=Ind - 1=Coll) cm,cm
    std::vector<double> fOmega_RMS;     // RMS of the 2D angular distribution  (0=Ind - 1=Coll) cm, cm
    
    std::vector<double> fOmega_Mean_reb;    // Mean value of the 2D angular Rebinned by 4
    std::vector<double> fOmega_RMS_reb;     // RMS of the 2D angular distribution  Rebinned by 4
    std::vector<double> fOmega_Mean_Mean;    // Mean value of the 2D angular use mean instead of maximum
 //   std::vector<double> fOmega_Mean_RMS;     // RMS of the 2D angular distribution use mean instead of maximum

   
    std::vector<double> fOmega_wt_Mean; // Mean value of the angular distribution (0=Ind - 1=Coll) wire,time
    std::vector<double> fOmega_wt_RMS;  // RMS of the angular distribution  (0=Ind - 1=Coll) wire,time

  
std::vector<std::vector<double> > fSingleEvtAngle;  // vector to show the plane omega distributions 
std::vector<std::vector<double> > fSingleEvtAngleVal;  //vector to show the plane omega distributions


    std::vector<unsigned int> fWire_vertex;  // wire coordinate of vertex for each plane
    std::vector<double> fTime_vertex;  // time coordinate of vertex for each plane

    std::vector<unsigned int> fWire_last;  // wire coordinate of last point for each plane
    std::vector<double> fTime_last;  // time coordinate of last point for each plane

 std::vector<unsigned int> fChannel_vertex;  // channel coordinate of vertex for each plane
 std::vector<unsigned int> fChannel_last;  // channel coordinate of vertex for each plane

 

  unsigned int tpc;    //tpc type
  unsigned int fNPlanes; // number of planes  

  TH1F* fh_theta[3]; 
  std::vector< TH1F * > fh_omega_evt_reb;
  std::vector< TH1F * > fh_omega_evt;
  TH1F* fh_theta_wt[3]; 
  TTree* ftree_cluster;



  }; // class ShowerAngleCluster

}

#endif // SHOWERANGLECLUSTER_H

