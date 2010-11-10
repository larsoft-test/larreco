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

    //METHODS THAT ALREADY WORK



    //METOHDS THAT NEEDS TO BE DONE
    edm::PtrVector<recob::Hit> SortHitList(edm::PtrVector<recob::Hit> hitlist);  // Sort a list of hit, maybe can be put in recob::hits when finished



    double DriftVelocity(double Efield, double Temperature); // Return the drift velocity in in cm/us
    int    Get3Daxis(std::vector<double> theta_polar, std::vector<double> Wire_vertex, std::vector<double> Time_vertex); // Get the 3d Shower axis in rad
    int    Get2Dvariables(std::vector<double> Wire_vertex_wt, std::vector<double> Time_wt); // Get the 2d view per plane in rad
    void   GetVertex(); // Get shower vertex
    void   GetPitchLength(double theta, double phi); //Get pitch length of both planes

    void   AngularDistribution(edm::PtrVector<recob::Hit> hitlist, int plane); //Get angular distribution of the shower 
    void   FitAngularDistributions(); 
    void   LongTransEnergy(edm::PtrVector<recob::Hit> hitlist, int plane); //Longtudinal and transverse enegry of the shower 

    std::string     fclusterModuleLabel;
    //std::string     fhitsModuleLabel;

    //truly needed global
    unsigned int planes;                //Number of planes, will be set by edm::Geometry in beginJob
    std::vector<unsigned int> wires;    //Number of wires per plane, will be set by edm::Geometry in beginJob

    std::vector<double> plane_pitch;    //pitch between the wire planes, will be set by edm::Geometry in beginJob
    std::vector<double> wire_pitch;     //pitch between each 2 wires in one plane, will be set by edm::Geometry in beginJob
    std::vector<double> mean_wire_pitch;//mean pitch between the wires in one plane, will be set by edm::Geometry in beginJob



    const static double pi = M_PI;


    double timetick; // time sample in us
    double driftvelocity; // drift velocity of electrons a liquids temperature in cm/us


    //to be checked if these can be defined also local 

    std::vector<double> theta_Mean;    // Mean value of the angular distribution cm,cm
    std::vector<double> theta_RMS;     // RMS of the angular distribution cm, cm
    std::vector<double> theta_wt_Mean; // Mean value of the angular distribution wire,time
    std::vector<double> theta_wt_RMS;  // RMS of the angular distribution wire,time

    std::vector<int>    LastWire;      // last wire of the shower
    std::vector<int>    LastTime;      // last t_hit of the shower

    // 2D slope and intercept of the shower axis per plane
    std::vector<double> slope;         // in cm, cm
    std::vector<double> intercept;     // in cm, cm
    std::vector<double> slope_wt;      // in wire, tick
    std::vector<double> intercept_wt;  // in wire,tick

    double theta, phi; // Polar Coordinates of the shower axis
    std::vector<float> wire_vertex;    // Wire Vertex position
    std::vector<float> time_vertex;    // Time Vertex position




//TODO SET CORRRECT innititialization ->DISCUSS WITH BIAGIO WHAT IS NEEDED FOR WHAT & with Brian how to set the parameterfile



    // Some variables for the hit
    std::vector<unsigned int> channel;  // channel number
    std::vector<float> time_hit; // hit time at maximum
    std::vector<int> wire_hit;    //hit wire number
    int plane_hit;             //hit plane number

    std::vector<TH1F*> h_theta; 
    std::vector<TH1F*> h_theta_wt; 
    std::vector<TH1F*> sh_nrg;
    std::vector<TH1F*> sh_Tnrg;
    std::vector<TH1F*> sh_long_hit;

    
    // Variables to get the angular distribution of the hits
    std::vector<double> a_polar, b_polar, theta_polar;



  protected: 
    

  }; // class ShowerReco

}

#endif // SHOWERRECO_H

