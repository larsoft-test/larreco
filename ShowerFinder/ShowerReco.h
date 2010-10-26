#ifndef SHOWERRECO_H
#define SHOWERRECO_H

#include "FWCore/Framework/interface/EDProducer.h" // include the proper bit of the framework

#include <vector>
#include <string>

//still needed, has to be put elsewhere?

#include "TH1.h"

namespace shwf {
   
  class Shower : public edm::EDProducer {
    
  public:
    
    explicit Shower(edm::ParameterSet const& pset); 
    virtual ~Shower();
    
    void produce(edm::Event& evt, edm::EventSetup const&);
    void beginJob(edm::EventSetup const&);

    
    double DriftVelocity(double Efield, double Temperature); // in cm/us
    int    Get3Daxis(double thetaI, double thetaC, double Wire_vertexI, double Wire_vertexC, double Time_vertex); // in rad
    int    Get2Dvariables(double Wire_vertexI_wt, double Wire_vertexC_wt, double Time_I_wt, double Time_C_wt); // in rad
    void   GetVertex(); // Get shower vertex
    void   GetPitchLength(double theta, double phi); //Get pitch length of both planes

//TO BE DONE :: WHAT IS TO BE CHANGED? 
    void   AngularDistributionI(std::vector<const recob::Hit*> hitlistInd, geo::Geometry* geom); //Get angular distribution of the shower (Induction plane)
    void   AngularDistributionC(std::vector<const recob::Hit*> hitlistCol, geo::Geometry* geom);  //Get angular distribution of the shower (Collection plane) 
    void   FitAngularDistributions(); 
    void   LongTransEnergyI(std::vector<const recob::Hit*> hitlistInd,  geo::Geometry* geom); //Longtudinal and transverse enegry of the shower (Induction plane)
    void   LongTransEnergyC(std::vector<const recob::Hit*> hitlistCol,  geo::Geometry* geom); //Longtudinal and transverse enegry of the shower (Collection plane)
    void   WriteHistos(edm::EventHandle& evt);
    float  McReleasedEnergy(const edm::EventHandle& evt);

    std::string     clusters;
    std::string     showers;  

    double theta_Mean[2];    // Mean value of the angular distribution (0=Ind - 1=Coll) cm,cm
    double theta_RMS[2];     // RMS of the angular distribution  (0=Ind - 1=Coll) cm, cm
    double theta_wt_Mean[2]; // Mean value of the angular distribution (0=Ind - 1=Coll) wire,time
    double theta_wt_RMS[2];  // RMS of the angular distribution  (0=Ind - 1=Coll) wire,time

    int    LastWire[2];  // last wire of the shower
    int    LastTime[2];  // last t_hit of the shower

    // 2D slope and intercept of the shower axis
    double slope[2];       // in cm, cm
    double intercept[2];   // in cm, cm
    double slope_wt[2];       // in wire, tick
    double intercept_wt[2];   // in wire,tick

    double TrackPitch[2];

    double theta, phi; // Director angles
    float wireI1; // Wire Vertex position (Induction)
    float timeI1; // Time Vertex position (Induction)
    float wire1C; // Wire Vertex poistion (Collection)
    float time1C; // Time Vertex position (Collection)
    const static double pi            = 3.14;
    const static double wire_pitch    =  0.4;   // wire pitch in cm
    const static double timetick      =  0.198; // time sample in us
    const static double presamplings  = 60.;
    const static double driftvelocity =  0.157;
    const static int    alpha           = 8;     // parameter (how many RMs (of the anglular distribution) is large the cone of the shower) 
    

 private:


    // Some variables for the hit
    unsigned int channel_I, channel_C;  // channel number
    float time_I,time_C; // hit time at maximum
    int wire_I, wire_C;    //hit wire number
    int plane;             //hit plane number

    TH1F* h_theta[3]; 
    TH1F* h_theta_wt[3]; 
    TH1F* sh_nrg[3];
    TH1F* sh_Tnrg[3];
    TH1F* sh_long_hit[3];

    // Save enegry in a file

    std::ofstream myfile;
    
    // Variables to get the angular distribution of the hits
    double AI, BI, thetaI; // Induction  plane
    double AC, BC, thetaC; // Collection plane
   
    //THOMAS MC
    float m_totenergy;

  protected: 
    

  }; // class ShowerReco

}

#endif // SHOWERRECO_H

