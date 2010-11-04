#ifndef SHOWERRECO_H
#define SHOWERRECO_H

#include "FWCore/Framework/interface/EDProducer.h" // include the proper bit of the framework

#include <vector>
#include <string>

#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>

#include "TH1F.h"

#include "RecoBase/Cluster.h"



namespace shwf {

  class ShowerReco : public edm::EDProducer {
    
  public:
    
    explicit ShowerReco(edm::ParameterSet const& pset); 
    virtual ~ShowerReco();
    
    void produce(edm::Event& evt, edm::EventSetup const&);
    void beginJob(edm::EventSetup const&);
   

 private:

    double DriftVelocity(double Efield, double Temperature); // in cm/us
    int    Get3Daxis(std::vector<double> theta_polar, std::vector<double> Wire_vertex, std::vector<double> Time_vertex); // in rad
    int    Get2Dvariables(std::vector<double> Wire_vertex_wt, std::vector<double> Time_wt); // in rad
    void   GetVertex(); // Get shower vertex
    void   GetPitchLength(double theta, double phi); //Get pitch length of both planes

    void   AngularDistributionI(edm::PtrVector<recob::Hit> hitlistInd); //Get angular distribution of the shower (Induction plane)
    void   AngularDistributionC(edm::PtrVector<recob::Hit> hitlistCol);  //Get angular distribution of the shower (Collection plane) 
    void   FitAngularDistributions(); 
    void   LongTransEnergyI(edm::PtrVector<recob::Hit> hitlistInd); //Longtudinal and transverse enegry of the shower (Induction plane)
    void   LongTransEnergyC(edm::PtrVector<recob::Hit> hitlistCol); //Longtudinal and transverse enegry of the shower (Collection plane)

    void   SortHitList(edm::PtrVector<recob::Hit> hitlist);




    std::string     fclusterModuleLabel;

    //std::string     fhitsModuleLabel;

    std::vector<double> theta_Mean;    // Mean value of the angular distribution cm,cm
    std::vector<double> theta_RMS;     // RMS of the angular distribution cm, cm
    std::vector<double> theta_wt_Mean; // Mean value of the angular distribution wire,time
    std::vector<double> theta_wt_RMS;  // RMS of the angular distribution wire,time

    std::vector<int>    LastWire;  // last wire of the shower
    std::vector<int>    LastTime;  // last t_hit of the shower

    // 2D slope and intercept of the shower axis per plane
    std::vector<double> slope;       // in cm, cm
    std::vector<double> intercept;   // in cm, cm
    std::vector<double> slope_wt;       // in wire, tick
    std::vector<double> intercept_wt;   // in wire,tick

    std::vector<double> TrackPitch;

    double theta, phi; // Polar Coordinates of the shower axis
    std::vector<float> wire_vertex;  // Wire Vertex position
    std::vector<float> time_vertex; // Time Vertex position

    const static double pi = M_PI;

//TODO SET CORRRECT CALL ->DISCUSS WITH BIAGIO WHAT IS NEEDEDFOR WHAT

    const static double timetick      =  0.198; // time sample in us
    const static double presamplings  = 60.;
    const static double driftvelocity =  0.157;
    const static int    alpha           = 8;     // parameter (how many RMs (of the anglular distribution) is large the cone of the shower) 

    // Some variables for the hit
    std::vector<unsigned int> channel;  // channel number
    std::vector<float> time_hit; // hit time at maximum
    std::vector<int> wire_hit;    //hit wire number
    int plane_hit;             //hit plane number

    TH1F* h_theta[3]; 
    TH1F* h_theta_wt[3]; 
    TH1F* sh_nrg[3];
    TH1F* sh_Tnrg[3];
    TH1F* sh_long_hit[3];

    
    // Variables to get the angular distribution of the hits
    std::vector<double> a_polar, b_polar, theta_polar;

  protected: 
    

  }; // class ShowerReco

}

#endif // SHOWERRECO_H

