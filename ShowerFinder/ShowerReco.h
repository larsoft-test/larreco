////////////////////////////////////////////////////////////////////////
/// \file  ShowerReco.h
/// \brief
///
///
/// \version $Id: SingleGen.cxx,v 1.4 2010/03/29 09:54:01 brebel Exp $
/// \author:  biagio, ART port: echurch
////////////////////////////////////////////////////////////////////////
#ifndef SHOWERRECO_H
#define SHOWERRECO_H

#include "art/Framework/Core/EDProducer.h" // include the proper bit of the framework

#include <vector>
#include <string>

#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>

#include "TH1F.h"
#include "TTree.h"

#include "RecoBase/Cluster.h"

namespace shwf {

  class ShowerReco : public art::EDProducer {
    
  public:

    /**METHODS global*/
    explicit ShowerReco(fhicl::ParameterSet const& pset);/**Constructor*/
    virtual ~ShowerReco();                               /**Destructor*/
    void beginJob();                                     
    void reconfigure(fhicl::ParameterSet const& pset);
    void produce(art::Event& evt);                       /**Actual routine that reconstruct the shower*/
   
    int    Get3Daxis(float thetaI, float thetaC, float Wire_vertexI, float Wire_vertexC, float Time_vertex); // in rad
    int    Get2Dvariables(float Wire_vertexI_wt, float Wire_vertexC_wt, float Time_I_wt, float Time_C_wt); // in rad
    void   GetVertex(); // Get shower vertex
    void   GetPitchLength(float theta, float phi); //Get pitch length of both planes
    void   AngularDistributionI(art::PtrVector < recob::Hit> hitlistInd); //Get angular distribution of the shower (Induction plane)
    void   AngularDistributionC(art::PtrVector < recob::Hit> hitlistCol);  //Get angular distribution of the shower (Collection plane) 
    void   FitAngularDistributions(); 
    void   LongTransEnergyI(art::PtrVector < recob::Hit> hitlistInd); //Longtudinal and transverse enegry of the shower (Induction plane)
    void   LongTransEnergyC(art::PtrVector < recob::Hit> hitlistCol); //Longtudinal and transverse enegry of the shower (Collection plane)


    float theta_Mean[2];    // Mean value of the angular distribution (0=Ind - 1=Coll) cm,cm
    float theta_RMS[2];     // RMS of the angular distribution  (0=Ind - 1=Coll) cm, cm
    float theta_wt_Mean[2]; // Mean value of the angular distribution (0=Ind - 1=Coll) wire,time
    float theta_wt_RMS[2];  // RMS of the angular distribution  (0=Ind - 1=Coll) wire,time
   
    float LastWire[2];  // last wire of the shower
    float LastTime[2];  // last t_hit of the shower
    
    // 2D slope and intercept of the shower axis
    float slope[3];       // in cm, cm
    float intercept[3];   // in cm, cm
    float slope_wt[3];       // in wire, tick
    float intercept_wt[3];   // in wire,tick

    float Pitch[3];
      
    float theta, phi; // Director angles
    float wireI1; // Wire Vertex position (Induction)
    float timeI1; // Time Vertex position (Induction)
    float wire1C; // Wire Vertex poistion (Collection)
    float time1C; // Time Vertex position (Collection)
    const static float pi            = 3.141519;
    float fmean_wire_pitch ;   // wire pitch in cm
    const static float ftimetick      =  0.198; // time sample in us
    const static float presamplings  = 60.;
    const static float fdriftvelocity =  0.157;
    const static int    alpha           = 5;     // parameter (how many RMs (of the anglular distribution) is large the cone of the shower)

 private:

  std::string fShwrOutput;
  std::string fClusterModuleLabel;
  float IdEdx4cm; // dedx of the first 4cm of the shower
  float CdEdx4cm; //dedx of the first 4cm of the shower 
  int IdedxCounter; // dedx of the first 4cm of the shower
  int CdedxCounter; //dedx of the first 4cm of the shower 
  
  // Some variables for the hit
  unsigned int channel_I, channel_C;  // channel number
  float time_I,time_C; // hit time at maximum
  unsigned int wire_I, wire_C;    //hit wire number
  unsigned int plane;    //hit plane number
  unsigned int tpc;    //hit plane number
  
  TH1F* fh_theta[3]; 
  TH1F* fh_theta_wt[3]; 
  TH1F* fh_dedx[3]; 
  TH1F* fsh_nrg[3];
  TH1F* fsh_Tnrg[3];
  TH1F* fsh_long_hit[3];
  TTree* ftree_shwf;

  // Save enegry in a file
  std::ofstream myfile;
  
  // Variables to get the angular distribution of the hits
  float AI, BI, thetaI; // Induction  plane
  float AC, BC, thetaC; // Collection plane


  }; // class ShowerReco

}

#endif // SHOWERRECO_H

