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
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include <vector>
#include <string>

#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>

#include "TH1F.h"
#include "TTree.h"
#include "TMath.h"
#include "TF1.h"
#include "TH2F.h"

#include "Utilities/LArProperties.h"
#include "Utilities/GeometryUtilities.h"
#include "Utilities/DetectorProperties.h"
#include "RecoBase/Cluster.h"

namespace cluster {

  class ShowerAngleCluster : public art::EDProducer {
    
  public:

    /**METHODS global*/
    explicit ShowerAngleCluster(fhicl::ParameterSet const& pset);/**Constructor*/
    virtual ~ShowerAngleCluster();                               /**Destructor*/
    void beginJob();                                     
    void beginRun(art::Run& run);
    void reconfigure(fhicl::ParameterSet const& pset);
    void produce(art::Event& evt);                       /**Routine that finds the cluster and sets the dTdW of the 2D shower*/
   
    void   Find2DAxisRough(unsigned int NClust,std::vector < art::Ptr < recob::Hit> > hitlist); /**Calculate 2D angle histograms, provided vertex is know */ 
    void   Get2DVariables(unsigned int nClust,std::vector < art::Ptr < recob::Hit> > hitlist);   /** Calculate 2D variables to be saved into 	  */ 
    void   FitAngularDistributions(std::vector < art::Ptr < recob::Hit> > hitlist);  /** Get actual 2D omega angle from the histograms  */
  

    void   GetVertexN(art::Event& evt); /**Get vertex from MCtruth information. Temporary. */


    
    //    float LastWire[2];  // last wire of the shower
    //    float LastTime[2];  // last t_hit of the shower
    
     
    

  private:

    
    double fWiretoCm,fTimetoCm,fWireTimetoCmCm;
    
    std::vector< unsigned int >fNWires;
    double fNTimes;
    
    
    
    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<util::DetectorProperties> detp;
    art::ServiceHandle<util::LArProperties> larp;
    art::ServiceHandle<geo::Geometry> geo;
    util::GeometryUtilities gser;
   
    float fTimeTick; // time sample in us
    float fPresamplings;
    float fDriftVelocity;
    const static int    alpha           = 5;     // parameter (how many RMs (of the anglular distribution) is large the cone of the shower)
   
   
    int fRun,fEvent,fSubRun;
   

    
    std::vector< unsigned int > fMinWire,fMaxWire;
    std::vector < double > fMinTime,fMaxTime;
    
    //input parameter labels:
 
    std::string fClusterModuleLabel;
    std::string fVertexCLusterModuleLabel;
    std::string fMCGeneratorLabel;
    std::string fLarGeantlabel;
    int fUseMCVertex;
  
  
    // 2D slope and intercept of the shower axis
    std::vector <double> slope;       // in wire, time
    std::vector <double> slope_cm;       // in cm, cm
    std::vector <double> xangle;       // in cm, cm
   
    std::vector <double> lineslope;   // in wire, time

    std::vector <double>  lineinterc;   

   
    std::vector<double> xyz_vertex;
    std::vector<double> xyz_vertex_fit;

 
    double fWirePitch ;   // wire pitch in cm

    std::vector<double> fOmega_Mean;    // Mean value of the 2D angular distribution (0=Ind - 1=Coll) cm,cm
    std::vector<double> fOmega_RMS;     // RMS of the 2D angular distribution  (0=Ind - 1=Coll) cm, cm
    
    //std::vector<double> fOmega_Mean_reb;    // Mean value of the 2D angular Rebinned by 4
    //std::vector<double> fOmega_RMS_reb;     // RMS of the 2D angular distribution  Rebinned by 4
    //std::vector<double> fOmega_Mean_Mean;    // Mean value of the 2D angular use mean instead of maximum
 
   // std::vector<double> fOmega_Mean_line;    // Mean value of the 2D angular distribution obtained from linear fit
   // std::vector<double> fOmega_RMS_line;     // RMS of the 2D angular distribution obtained from linear fit

   
   // std::vector<double> fOmega_wt_Mean; // Mean value of the angular distribution (0=Ind - 1=Coll) wire,time
   // std::vector<double> fOmega_wt_RMS;  // RMS of the angular distribution  (0=Ind - 1=Coll) wire,time

  
    std::vector<std::vector<double> > fSingleEvtAngle;  // vector to show the plane omega distributions 
    std::vector<std::vector<double> > fSingleEvtAngleVal;  //vector to show the plane omega distributions


    std::vector<std::vector<double> > fShowerWidthProfile2D;  // vector to show the plane shower Width distribution 
    std::vector<std::vector<double> > fShowerChargeProfile2D;  //vector to show the plane shower Charge distribution
    std::vector<std::vector<double> > fShowerPosition2D;  //vector to store the positions of hit values stored in the previous two vectors.

    std::vector<unsigned int> fWireVertex;  // wire coordinate of vertex for each plane
    std::vector<double> fTimeVertex;  // time coordinate of vertex for each plane

    std::vector<unsigned int> fWireLast;  // wire coordinate of last point for each plane
    std::vector<double> fTimeLast;  // time coordinate of last point for each plane

    
    std::vector<double>  test_wire_start,test_time_start;

    std::vector<double> fRMS_wire;
    std::vector<double> fRMS_time;
    std::vector<double> fChisq;
    std::vector<double> fcorrelation;
    std::vector<double> fcovariance;



    // unsigned int tpc;    //tpc type
    unsigned int fNPlanes; // number of planes  

    TH1F* fh_theta[3]; 
   // std::vector< TH1F * > fh_omega_evt_reb;
    std::vector< TH1F * > fh_omega_evt;
    TH1F *  fh_omega_single;
    TH1F* fh_theta_wt[3]; 
    TTree* ftree_cluster;
    std::vector < TF1 * > linefit_cm;
    std::vector < TF1 * > linefit2_cm;	
    std::vector < TH2F * > tgx;
    std::vector < TH2F * > tgx2;

    // TH1F * hithist[3];
    void CalculateAxisParameters(unsigned nClust, std::vector < art::Ptr < recob::Hit> >  hitlist,double wstart,double tstart,double wend,double tend);
    
    void Find2DStartPoints(int nClust,std::vector< art::Ptr < recob::Hit> > hitlist);

    void Find2DBestPlanes(std::vector<int> &best_planes);
  
    void Find_Extreme_Intercepts(std::vector < art::Ptr<recob::Hit> > hitlist,
			     double perpslope,
			     double &inter_high,
			     double &inter_low);  
    
    art::Ptr<recob::Hit> FindClosestHit(std::vector < art::Ptr < recob::Hit> > hitlist,
					unsigned int wire,
					double time);  
    
    void RefineStartPoints(unsigned int nClust, std::vector< art::Ptr < recob::Hit> >  hitlist, double  wire_start,double time_start);
    
    double Get2DAngleForHit( art::Ptr<recob::Hit> starthit,std::vector < art::Ptr < recob::Hit> > hitlist);
    
    double Get2DAngleForHit( unsigned int wire, double time,std::vector < art::Ptr < recob::Hit> > hitlist);
    
    
    int GetPlaneAndTPC(art::Ptr<recob::Hit> a,unsigned int &p,unsigned int &cs,unsigned int &t,unsigned int &w);
    
    void FindMinMaxWiresTimes(unsigned int nClust,std::vector < art::Ptr < recob::Hit> >  hitlist);
    
    void ClearandResizeVectors(unsigned int nClusters);
    
 //   double Get2DDistance(double wire1,double time1,double wire2,double time2);
  
 //   double Get2DPitchDistance(double angle,double inwire,double wire);
    
 //   int GetPointOnLine(double slope,double intercept,double wire1,double time1,double &wireout,double &timeout);
    
 //   int GetPointOnLine(double slope,double wirestart,double timestart,double wire1,double time1,double &wireout,double &timeout);
    
 //   int GetPointOnLineWSlopes(double slope,double intercept,double ort_intercept,double &wireout,double &timeout);
    
  
    
    //temporary
    int mcpdg;
    double mcenergy;
    double mcphi;
    double mctheta;
    std::vector< unsigned int> mcwirevertex;  // wire coordinate of vertex for each plane
    std::vector< double> mctimevertex;  // time coordinate of vertex for each plane

    bool startflag;
    bool endflag;
    bool matchflag;
    
    
  }; // class ShowerAngleCluster

}

#endif // SHOWERANGLECLUSTER_H

