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
#include "TTree.h"

#include "RecoBase/Cluster.h"



namespace shwf {

  class ShowerReco : public edm::EDProducer {
    
  public:

    /**METHODS global*/
    explicit ShowerReco(edm::ParameterSet const& pset);    /**Constructor*/
    virtual ~ShowerReco();                                 /**Deconstructor*/
    void beginJob(edm::EventSetup const&);                 /**Needed by Art as initialzer*/ 
    void produce(edm::Event& evt, edm::EventSetup const&); /**Actual routine that reconstruct the shower*/
   

 private:

    /**METHODS private*/
    void  AngularDistribution(std::vector<edm::PtrVector<recob::Hit> > hitlist); /**Get angular distribution of the shower in each plane */
    void  FitAngularDistributions();                                             /**Fit the angular distibution*/ 
    void  Get3Daxis();                                                           /** Get the 3d Shower axis*/
    void  LongTransEnergy(std::vector<edm::PtrVector<recob::Hit> > hitlist);     /**Longitudinal and transverse enegry of the shower*/ 
    void  Get2Dvariables();                                                      /** Get the 2d view per plane in rad*/


    /**Variables needed as global*/
    std::vector<unsigned int> wires;      /**Number of wires per plane, will be set by edm::Geometry in beginJob*/
    std::vector<double> plane_pitch;      /**Pitch between the wire planes, will be set by edm::Geometry in beginJob*/
    std::vector<double> wire_pitch;       /**Pitch between each 2 wires in one plane, will be set by edm::Geometry in beginJob*/
    std::vector<double> mean_wire_pitch;  /**Mean pitch between the wires in one plane, will be set by edm::Geometry in beginJob*/
    std::vector<unsigned int> wire_vertex;/**Vertex position in #wire */
    std::vector<float> time_vertex;       /**Vertex position in time*/
    std::vector<double> theta_Mean;       /**Mean value of the angular shower distribution cm,cm*/
    std::vector<double> theta_RMS;        /**RMS of the angular shower distribution cm, cm*/
    double theta, phi;                    /**Polar Coordinates of the shower axis*/
    const static double pi = M_PI;        /**short handle to access Pi*/
    const static int alpha = 8;           /**Parameter: #RMS of the anglular distribution which are part of the shower cone*/
    std::vector<int> LastWire;            /**Last wire_hit of the shower*/
    std::vector<float> LastTime;          /**Last time_hit of the shower*/
    std::vector<double> totCharge;        /**Total Charge per plane in MIPs*/
    std::vector<double> slope_2d;         /** 2d information from reconstructed shower for line*/
    std::vector<double> intercept_2d;     /** with form time = slope * wire + intercept in cm, cm per plane */


 
    TTree* tree;
    std::string     fclusterModuleLabel;

    /**Histograms to return the shower reconstruction*/
    std::vector<TH1F*> h_theta;     /**Histo for the angular distribution theta of the shower*/
    std::vector<TH1F*> h_theta_wt;  /**Histo for the angular distribution theta wire of the shower*/
    std::vector<TH1F*> sh_nrg;      /**Histo for the longitudinal energy distribution of the shower*/
    std::vector<TH1F*> sh_Tnrg;     /**Histo for the longitudinal energy distribution of the shower*/
    std::vector<TH1F*> sh_long_hit; /**Histo for the Transverse HIT distribution of the shower*/

    /**To be checked how to get by external handle...*/

    double timetick; // time sample in us
    double driftvelocity; // drift velocity of electrons a liquids temperature in cm/us

  protected: 
    

  }; // class ShowerReco

}

#endif // SHOWERRECO_H

