////////////////////////////////////////////////////////////////////////
//
// ShowerReco class - V.1.0 - 11/16/2010 
//
// biagio.rossi@lhep.unibe.ch   (FWMK : argoneut specific)
// thomas.strauss@lhep.unibe.ch (ART  : general detector)
//
// This algorithm is designed to reconstruct showers
// 
///////////////////////////////////////////////////////////////////////

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>

// Framework includes
#include "FWCore/Framework/interface/Event.h" 
#include "FWCore/ParameterSet/interface/ParameterSet.h" 
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h" 
#include "DataFormats/Common/interface/Handle.h" 
#include "DataFormats/Common/interface/Ptr.h" 
#include "DataFormats/Common/interface/PtrVector.h" 
#include "FWCore/Framework/interface/MakerMacros.h" 
#include "FWCore/ServiceRegistry/interface/Service.h" 
#include "FWCore/Services/interface/TFileService.h" 
#include "FWCore/Framework/interface/TFileDirectory.h" 
#include "FWCore/MessageLogger/interface/MessageLogger.h" 
#include "FWCore/ServiceRegistry/interface/ServiceMaker.h" 

#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompSVD.h"
#include "TH2F.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"

// LArSoft includes
#include "Simulation/sim.h"
#include "ShowerFinder/ShowerReco.h"
#include "Geometry/geo.h"
#include "RecoBase/recobase.h"


// ***************** //
shwf::ShowerReco::ShowerReco(edm::ParameterSet const& pset) : 

  fclusterModuleLabel (pset.getParameter< std::string >("clusters"))
  
{
  produces< std::vector<recob::Shower> >();
}

// ***************** //
shwf::ShowerReco::~ShowerReco()
{
}

// ***************** //
void shwf::ShowerReco::beginJob(edm::EventSetup const&)
{
  /** Get Geometry*/
  edm::Service<geo::Geometry> geo;
  unsigned int planes = geo->Nplanes();
  
  double mean_wirepitch;

  for (unsigned int i=0;i<planes;++i){
    mean_wirepitch = 0.;
    plane_pitch.push_back(geo->PlanePitch(i,i+1));
    wires.push_back(geo->Nwires(i));
    for (unsigned int j=0;j<wires[i]-1;++j){
      wire_pitch.push_back(geo->WirePitch(j,j+1,i));
      mean_wirepitch+=geo->WirePitch(j,j+1,i);
    }
    mean_wire_pitch.push_back(mean_wirepitch/geo->Nwires(i));
  }
  /**Get TFileService and define output Histograms*/
  edm::Service<edm::TFileService> tfs;

  /** Create Histos names*/
  char tit_h_theta[128] = {0};
  char tit_h_theta_wt[128] = {0};
  char sh_long_tit[128] = {0};
  char sh_tit[128] = {0};
  char shT_tit[128] = {0};
  
  int nbins;
  for(unsigned int i=0;i<planes;++i){

    /**Histos for the angular distribution theta of the shower*/
    sprintf(&tit_h_theta[0],"h_theta_%i",i);
    h_theta.push_back(tfs->make<TH1F>(tit_h_theta,"Theta distribution",180,-180., 180.));

    /**Histos for the angular distribution theta wire of the shower*/
    sprintf(&tit_h_theta[0],"theta_wire_%i",i);
    h_theta_wt.push_back(tfs->make<TH1F>(tit_h_theta_wt,"Theta wire distribution",45,-180., 180.));
 
    /**Histos for the longitudinal energy distribution of the shower */
    sprintf(&sh_tit[0],"sh_nrg1_%i",i);                   /**number of wires used,min wire,max_wire you need for the anlysis*/
    sh_nrg.push_back(tfs->make<TH1F>(sh_tit,"energy reco",wires[i],               0.,      wires[i]*mean_wire_pitch[i]));

    /**Histos for the transverse energy distribution of the shower*/
    sprintf(&shT_tit[0],"shT_nrg1_%i",i);                   /**units are ticks most lickely, but how did Biagio get size of it???*/
    sh_Tnrg.push_back(tfs->make<TH1F>(shT_tit,"energy reco",80,-40.,40.));
  
    /**Histos for the Transverse HIT distribution of the shower*/
    nbins = (int)(wires[i]*mean_wire_pitch[i]);
    sprintf(&sh_long_tit[0],"sh_long_hit_%i",i);                           /**nbins,min wire,max_wire you need for the analysis*/
    sh_long_hit.push_back(tfs->make<TH1F>(sh_long_tit,"longitudinal hit reco",nbins, 0.,     wires[i]*mean_wire_pitch[i]));
  }

  tree =tfs->make<TTree>("ShowerReco","Results");/**All-knowing tree with reconstruction information*/
  tree->Branch("theta_Mean","std::vector<double>", &theta_Mean);
  tree->Branch("theta_RMS","std::vector<double>", &theta_RMS);
  tree->Branch("theta","double",&theta);
  tree->Branch("phi","double",&phi);
  tree->Branch("totCharge","std::vector<double>",&totCharge);
  tree->Branch("mean_wire_pitch","std::vector<double>", &mean_wire_pitch);
  tree->Branch("wire_vertex","std::vector<unsigned int>", &wire_vertex);
  tree->Branch("time_vertex","std::vector<float>", &time_vertex);
  tree->Branch("slope_2d"," std::vector<double>", &slope_2d);
  tree->Branch("intercept_2d","std::vector<double>", &intercept_2d);


}

// ***************** //
void shwf::ShowerReco::produce(edm::Event& evt, edm::EventSetup const&)
{ 

  /**TODO THIS VALUES SHOULD BE PARAMETERS OF THE MODULE, evtl from database
  *double Efield_SI        =  0.7;     // Electric Field between Shield and Induction planes in kV/cm
  *double Efield_IC        =  0.9;     // Electric Field between Induction and Collection planes in kV/cm
  *double Temperature      = 87.6;  // LAr Temperature in K
  *check if there can be a replacement later for the product, not needed now
  *double timepitch        = driftvelocity*timetick;                  //time sample (cm) */

  /** Get Geometry */
  edm::Service<geo::Geometry> geo;
  //edm::Service<util::LArProperties> larprop;
  unsigned int planes = geo->Nplanes();
  /**TODO GET THESE FROM A PARAMETERSET/DATABASE*/
  timetick      =  0.198; //get from parameterset
  driftvelocity =  0.157;  //get from paramtereset 9either k and V) 
  //driftvelocity = larprob->DriftVelocity(Efield_SI,Temperature);
  
  /**Get Clusters*/
  edm::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fclusterModuleLabel,clusterListHandle);


  edm::PtrVector<recob::Hit> hitlist_helpwire; /** Define hitlist of all hits per wire*/
  edm::PtrVector<recob::Hit> hitlist_helpplane; /** Define hitlist to sort by #wire in hitlist_helpwire*/
  std::vector<edm::PtrVector<recob::Hit> > hitlist_plane; // Define vector of hitslists for each Induction and Collection planes

  /**Loop to fill the std::vector<edm:PtrVector<recob::Hit>> with the hits of each plane*/
  for(unsigned int iplane = 0; iplane < planes; ++iplane){//Loop over planes
    geo::View_t view = geo->Plane(iplane).View();

    hitlist_helpplane.clear();
    for(unsigned int iwire=0; iwire<wires[iplane];++iwire){//Loop over wires
      hitlist_helpwire.clear();
      for(int iclust = 0; iclust < clusterListHandle->size(); iclust++){//Loop over cluster
        edm::Ptr<recob::Cluster> clust(clusterListHandle, iclust);
        if(view == clust->View()){//If correct plane
          hitlist_helpplane = clust->Hits(iplane,iwire);//Fill with Hits from the wire (hitlist_help is automatically sorted)
        }
      for (int ihits =0; ihits<hitlist_helpplane.size(); ++ihits){//work around since clust->Hits return a std::vector,recob::Hit*> object
        edm::Ptr<recob::Hit> hit_help ;
        hit_help = hitlist_helpplane[ihits];
        hitlist_helpwire.push_back(hit_help);
        }
      }
    }
    hitlist_plane.push_back(hitlist_helpplane);
  }

  /**Get the angular distribution of the shower*/
  AngularDistribution(hitlist_plane);
  /**Fit the angular distributions*/
  FitAngularDistributions();
  /**Get the 3d axis of the shower cone */
  Get3Daxis();
  /**Determine the enrgy release in each plane*/
  LongTransEnergy(hitlist_plane);
  /**Get back the 2d variables*/
  Get2Dvariables();
  /**Fill the output tree with all information */
  tree->Fill();
}

// ***************** //
void shwf::ShowerReco::AngularDistribution(std::vector<edm::PtrVector<recob::Hit> > hitlist_all){ 

  /** Determine the angular distribution of the shower energy in each plane */
  edm::Service<geo::Geometry> geo;
  unsigned int planes = geo->Nplanes();
  unsigned int iplane_hit,iwire_hit,ichannel_hit;
  float ftime_hit;
  double a_polar,b_polar,theta_polar;

  for(unsigned int iplane = 0; iplane < planes; ++iplane){/**Loop over planes*/
    edm::PtrVector<recob::Hit> hitlist = hitlist_all[iplane];
    for (int ihits=0;ihits<hitlist.size();++ihits){
      edm::Ptr<recob::Hit> hit_help ;
      hit_help = hitlist[ihits]; 
      const recob::Hit* hit = hit_help.get(); /** Retrieve info of the hits*/
      ftime_hit = hit->CrossingTime(); /** Hit crossing time*/
      ichannel_hit = hit->Wire()->RawDigit()->Channel();
      geo->ChannelToWire(ichannel_hit, iplane_hit, iwire_hit); /** Get Information of wire/plane of the hit */

      /** Determine the vertex*/
      if (ihits==0) {    
        wire_vertex.push_back(iwire_hit); 
        time_vertex.push_back(ftime_hit);
        }

      /** moving to polar coordinates (cm,cm coordinate) */
      b_polar = (iwire_hit - wire_vertex[iplane])* mean_wire_pitch[iplane]+ mean_wire_pitch[iplane]; /**in cm*/
      a_polar = (ftime_hit - time_vertex[iplane])* timetick *driftvelocity; /** in cm*/ 
      theta_polar = asin(a_polar/sqrt(pow(a_polar,2)+pow(b_polar,2))); /**in rad*/
      theta_polar = 180*theta_polar/pi; /** in deg*/

      /** Filling the histo (angle, energy of the hit)*/
      h_theta[iplane]->Fill(theta_polar, hit_help->MIPs()); /** angle in deg (cm,cm) coordinate*/
      }

    std::cout << "Plane " << iplane << "VertexWire= " << wire_vertex[iplane] << "   Time= " << time_vertex[iplane] << std::endl;

    }
  return;
}

// ***************** //
void shwf::ShowerReco::FitAngularDistributions(){
  /** Fit function of the angular distribution (cm,cm)*/
  edm::Service<geo::Geometry> geo;
  unsigned int planes = geo->Nplanes();
  TF1 *gau = new TF1("gaus","gaus",-60, 60);

  for(unsigned int iplane = 0; iplane < planes; ++iplane){
    h_theta[iplane]->Fit("gaus","QR");        /** Fit of the angular distribution*/
    theta_Mean[iplane] = gau->GetParameter(1);/** Mean value of the fit */
    theta_RMS[iplane] = gau->GetParameter(2); /** RMS of the fit of the angular distribution in deg*/
  }
}

// ***************** //
void shwf::ShowerReco::Get3Daxis(){
  edm::Service<geo::Geometry> geo;
  unsigned int planes = geo->Nplanes();
  double a_polar,b_polar,theta_polar;
  float wire_help,time_help;
  double slope[2];         /** Need only 2 * 2d information to reconstruct 3d for line in 3d*/
  double intercept[2];     /** with form time = slope * wire + interceptin cm, cm */

  for(unsigned int iplane = 0; iplane < planes; ++iplane){/**Loop over planes*/
    /** moving to polar coordinates (cm,cm coordinate) */
    b_polar = wire_vertex[iplane]* mean_wire_pitch[iplane]+ mean_wire_pitch[iplane]; /**in cm*/
    a_polar = time_vertex[iplane]* timetick *driftvelocity; /** in cm*/ 
    theta_polar = asin(a_polar/sqrt(pow(a_polar,2)+pow(b_polar,2))); /**in rad*/
    theta_polar = 180*theta_polar/pi; /** in deg*/

    /**Calculate the shower cone direction in each plane time=slope[plane]*wire + intercept[plane]*/
    if(geo->Plane(iplane).SignalType() == geo::kCollection){
      slope[1] = tan((theta_polar*pi/180));
      wire_help = wire_vertex[iplane] *mean_wire_pitch[iplane]; // in cm 
      time_help = time_vertex[iplane] *timetick *driftvelocity; // in cm
      intercept[1] = time_help-(slope[1]*wire_vertex[iplane]);
      }
    else if (geo->Plane(iplane).SignalType() == geo::kInduction) {
      slope[0] = tan((theta_polar*pi/180));
      wire_help = wire_vertex[iplane] *mean_wire_pitch[iplane]; // in cm 
      time_help = time_vertex[iplane] *timetick *driftvelocity; // in cm
      intercept[0] = time_help-(slope[0]*wire_vertex[iplane]);
      }
    }

  /**Calculate 3d axis*/
  double l(0),m(0),n(0);
  double angle_rad = 30* pi /180;
  l = 1;
  m = (1/(2*sin(angle_rad)))*((1/slope[1])-(1/slope[0]));
  n = (1/(2*cos(angle_rad)))*((1/slope[1])+(1/slope[0]));
 
  phi   = atan(n/l);
  theta = acos(m/(sqrt(pow(l,2)+pow(m,2)+pow(n,2))));

  phi = phi<0. ? phi+TMath::Pi() : phi ; /** solve the ambiguities due to tangent periodicity*/

  phi = phi*180/pi;
  theta = theta*180/pi;
  std::cout << "Shower axis: Phi=" <<phi << "   theta=" << theta <<std::endl;
  return;
}

// ***************** //
void shwf::ShowerReco::LongTransEnergy(std::vector<edm::PtrVector<recob::Hit> > hitlist_all){

  edm::Service<geo::Geometry> geo;
  unsigned int planes = geo->Nplanes();
  unsigned int iplane_hit,iwire_hit,ichannel_hit;
  float ftime_hit;
  double a_polar,b_polar,theta_polar, wire_rot, time_rot;
  double low_th, high_th, theta_sh, totcharge;   /** thresholds for adding up only the energy of the shower in a cone of alpha*RMS the angular distribution of the shower itself*/

  for(unsigned int iplane = 0; iplane < planes; ++iplane){/**Loop over planes*/
    totcharge =0.;
    edm::PtrVector<recob::Hit> hitlist = hitlist_all[iplane];
    for (int ihits=0;ihits<hitlist.size();++ihits){
      edm::Ptr<recob::Hit> hit_help ;
      hit_help = hitlist[ihits]; 
      const recob::Hit* hit = hit_help.get(); /** Retrieve info of the hits*/
      ftime_hit = hit->CrossingTime(); /** Hit crossing time*/
      ichannel_hit = hit->Wire()->RawDigit()->Channel();
      geo->ChannelToWire(ichannel_hit, iplane_hit, iwire_hit); /** Get Information of wire/plane of the hit */

      b_polar = (iwire_hit - wire_vertex[iplane])* mean_wire_pitch[iplane]+ mean_wire_pitch[iplane]; /**in cm*/
      a_polar = (ftime_hit - time_vertex[iplane])* timetick *driftvelocity; /** in cm*/ 
      theta_polar = asin(a_polar/sqrt(pow(a_polar,2)+pow(b_polar,2))); /**in rad*/
      theta_polar = 180*theta_polar/pi; /** in deg*/

      low_th  = theta_Mean[iplane]-(alpha*theta_RMS[iplane]);
      high_th = theta_Mean[iplane]+(alpha*theta_RMS[iplane]);

      if((theta_polar<low_th) || (theta_polar>high_th)) continue;

      if( (theta_polar>(theta_Mean[iplane]-0.1))&&(theta_polar<(theta_Mean[iplane]+0.1)) ){
        LastWire[iplane] = iwire_hit;
        LastTime[iplane] = ftime_hit; 
      }

      theta_sh  = theta_Mean[iplane]*pi/180; // in radians
      wire_rot = (iwire_hit - wire_vertex[iplane])* mean_wire_pitch[iplane] *cos(theta_sh) + (ftime_hit - time_vertex[iplane])* timetick *driftvelocity *sin(theta_sh);
      time_rot = (ftime_hit - time_vertex[iplane])* timetick *driftvelocity /cos(theta_sh) - (iwire_hit - wire_vertex[iplane])* mean_wire_pitch[iplane] *sin(theta_sh) 
                 - (ftime_hit - time_vertex[iplane])* timetick *driftvelocity *tan(theta_sh) *sin(theta_sh);   
    
      totcharge +=hit->MIPs(); // Sum the energy of all the hits
    
      sh_nrg[iplane]->Fill(wire_rot, hit->MIPs()); // Fill histo longitudinal distr of the energy (IND)
      sh_Tnrg[iplane]->Fill(time_rot, hit->MIPs());// Fill histo transverse   distr of the energy (IND)
      sh_long_hit[iplane]->Fill(wire_rot, 1);
    
    }
    totCharge.push_back(totcharge);
    std::cout << "In plane" << iplane << "totCharge (Mips)= " << totcharge<< std::endl;
  }
  return;
}

// ***************** //
void shwf::ShowerReco::Get2Dvariables(){
  edm::Service<geo::Geometry> geo;
  unsigned int planes = geo->Nplanes();
  for(unsigned int iplane = 0; iplane < planes; ++iplane){/**Loop over planes*/  
    slope_2d.push_back((LastTime[iplane]-time_vertex[iplane]) /(LastWire[iplane]-wire_vertex[iplane])); 
    intercept_2d.push_back(time_vertex[iplane] - wire_vertex[iplane]*slope_2d[iplane]);
    std::cout << "Cone in plane " << iplane << "slope " << slope_2d[iplane] << "intercept " << intercept_2d[iplane] << std::endl;
    }
  return;
}
