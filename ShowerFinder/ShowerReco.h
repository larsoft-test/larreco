////////////////////////////////////////////////////////////////////////
/// \file  ShowerReco.h
/// \brief
///
///
/// \version $Id: SingleGen.cxx,v 1.4 2010/03/29 09:54:01 brebel Exp $
/// \author  biagio
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
    virtual ~ShowerReco();                               /**Deconstructor*/
    void beginJob();                                     /**Needed by Art as initialzer*/ 
    void produce(art::Event& evt);                       /**Actual routine that reconstruct the shower*/
   

 private:

    /**METHODS private*/
    void  AngularDistribution(std::vector<art::PtrVector<recob::Hit> > hitlist); /**Get angular distribution of the shower in each plane */
    void  FitAngularDistributions();                                             /**Fit the angular distibution*/ 
    void  Get3Daxis();                                                           /** Get the 3d Shower axis*/
    void  LongTransEnergy(std::vector<art::PtrVector<recob::Hit> > hitlist);     /**Longitudinal and transverse enegry of the shower*/ 
    void  Get2Dvariables();                                                      /** Get the 2d view per plane in rad*/


    /**Variables needed as global*/
    std::vector<unsigned int> fwires;      /**Number of wires per plane, will be set by edm::Geometry in beginJob*/
    std::vector<double> fplane_pitch;      /**Pitch between the wire planes, will be set by edm::Geometry in beginJob*/
    std::vector<double> fwire_pitch;       /**Pitch between each 2 wires in one plane, will be set by edm::Geometry in beginJob*/
    std::vector<double> fmean_wire_pitch;  /**Mean pitch between the wires in one plane, will be set by edm::Geometry in beginJob*/
    std::vector<unsigned int> fwire_vertex;/**Vertex position in #wire */
    std::vector<float> ftime_vertex;       /**Vertex position in time*/
    std::vector<double> ftheta_Mean;       /**Mean value of the angular shower distribution cm,cm*/
    std::vector<double> ftheta_RMS;        /**RMS of the angular shower distribution cm, cm*/
    double ftheta, fphi;                    /**Polar Coordinates of the shower axis*/
    const static double fpi = M_PI;        /**short handle to access Pi*/
    const static int falpha = 8;           /**Parameter: #RMS of the anglular distribution which are part of the shower cone*/
    std::vector<int> fLastWire;            /**Last wire_hit of the shower*/
    std::vector<float> fLastTime;          /**Last time_hit of the shower*/
    std::vector<double> ftotCharge;        /**Total Charge per plane in MIPs*/
    std::vector<double> fslope_2d;         /** 2d information from reconstructed shower for line*/
    std::vector<double> fintercept_2d;     /** with form time = slope * wire + intercept in cm, cm per plane */

 
    TTree* ftree_shwf;
    std::string     fclusterModuleLabel;

    /**Histograms to return the shower reconstruction*/
    std::vector<TH1F*> fh_theta;     /**Histo for the angular distribution theta of the shower*/
    std::vector<TH1F*> fh_theta_wt;  /**Histo for the angular distribution theta wire of the shower*/
    std::vector<TH1F*> fsh_nrg;      /**Histo for the longitudinal energy distribution of the shower*/
    std::vector<TH1F*> fsh_Tnrg;     /**Histo for the longitudinal energy distribution of the shower*/
    std::vector<TH1F*> fsh_long_hit; /**Histo for the Transverse HIT distribution of the shower*/

    /**To be checked how to get by external handle...*/

    double ftimetick; // time sample in us
    double fdriftvelocity; // drift velocity of electrons a liquids temperature in cm/us

  protected: 
    

  }; // class ShowerReco

}

#endif // SHOWERRECO_H

