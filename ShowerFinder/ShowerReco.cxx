////////////////////////////////////////////////////////////////////////
//
// \file ShowerReco.cxx
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
#include "art/Framework/Core/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Persistency/Common/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Core/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

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
shwf::ShowerReco::ShowerReco(fhicl::ParameterSet const& pset) : 

  fclusterModuleLabel (pset.get< std::string >("clusters"))
  
{
  produces< std::vector<recob::Shower> >();
}

// ***************** //
shwf::ShowerReco::~ShowerReco()
{
}

// ***************** //
void shwf::ShowerReco::beginJob()
{
  /** Get Geometry*/
  art::ServiceHandle<geo::Geometry> geo;
  unsigned int planes = geo->Nplanes();
  
  double mean_wirepitch;

  for (unsigned int i=0;i<planes;++i){
    mean_wirepitch = 0.;
    fplane_pitch.push_back(geo->PlanePitch(i,i+1));
    fwires.push_back(geo->Nwires(i));
    for (unsigned int j=0;j<fwires[i]-1;++j){
      fwire_pitch.push_back(geo->WirePitch(j,j+1,i));
      mean_wirepitch+=geo->WirePitch(j,j+1,i);
    }
    fmean_wire_pitch.push_back(mean_wirepitch/geo->Nwires(i));
  }
  /**Get TFileService and define output Histograms*/
  art::ServiceHandle<art::TFileService> tfs;

  /** Create Histos names*/
  char tit_h_theta[128] = {0};
  char tit_h_theta_wt[128] = {0};
  char sh_long_tit[128] = {0};
  char sh_tit[128] = {0};
  char shT_tit[128] = {0};
  
  int nbins;
  for(unsigned int i=0;i<planes;++i){

    /**Histos for the angular distribution theta of the shower*/
    sprintf(&tit_h_theta[0],"fh_theta_%i",i);
    fh_theta.push_back(tfs->make<TH1F>(tit_h_theta,"Theta distribution",180,-180., 180.));

    /**Histos for the angular distribution theta wire of the shower*/
    sprintf(&tit_h_theta[0],"ftheta_wire_%i",i);
    fh_theta_wt.push_back(tfs->make<TH1F>(tit_h_theta_wt,"Theta wire distribution",45,-180., 180.));
 
    /**Histos for the longitudinal energy distribution of the shower */
    sprintf(&sh_tit[0],"fsh_nrg1_%i",i);                   /**number of wires used,min wire,max_wire you need for the anlysis*/
    fsh_nrg.push_back(tfs->make<TH1F>(sh_tit,"energy reco",fwires[i],               0.,      fwires[i]*fmean_wire_pitch[i]));

    /**Histos for the transverse energy distribution of the shower*/
    sprintf(&shT_tit[0],"fshT_nrg1_%i",i);                   /**units are ticks most lickely, but how did Biagio get size of it???*/
    fsh_Tnrg.push_back(tfs->make<TH1F>(shT_tit,"energy reco",80,-40.,40.));
  
    /**Histos for the Transverse HIT distribution of the shower*/
    nbins = (int)(fwires[i]*fmean_wire_pitch[i]);
    sprintf(&sh_long_tit[0],"fsh_long_hit_%i",i);                           /**nbins,min wire,max_wire you need for the analysis*/
    fsh_long_hit.push_back(tfs->make<TH1F>(sh_long_tit,"longitudinal hit reco",nbins, 0.,     fwires[i]*fmean_wire_pitch[i]));
  }

  ftree_shwf =tfs->make<TTree>("ShowerReco","Results");/**All-knowing tree with reconstruction information*/
  ftree_shwf->Branch("ftheta_Mean","std::vector<double>", &ftheta_Mean);
  ftree_shwf->Branch("ftheta_RMS","std::vector<double>", &ftheta_RMS);
  ftree_shwf->Branch("ftheta","double",&ftheta);
  ftree_shwf->Branch("fphi","double",&fphi);
  ftree_shwf->Branch("ftotCharge","std::vector<double>",&ftotCharge);
  ftree_shwf->Branch("fmean_wire_pitch","std::vector<double>", &fmean_wire_pitch);
  ftree_shwf->Branch("fwire_vertex","std::vector<unsigned int>", &fwire_vertex);
  ftree_shwf->Branch("ftime_vertex","std::vector<float>", &ftime_vertex);
  ftree_shwf->Branch("fslope_2d"," std::vector<double>", &fslope_2d);
  ftree_shwf->Branch("fintercept_2d","std::vector<double>", &fintercept_2d);


}

// ***************** //
void shwf::ShowerReco::produce(art::Event& evt)
{ 

  /**TODO THIS VALUES SHOULD BE PARAMETERS OF THE MODULE, evtl from database
  *double Efield_SI        =  0.7;     // Electric Field between Shield and Induction planes in kV/cm
  *double Efield_IC        =  0.9;     // Electric Field between Induction and Collection planes in kV/cm
  *double Temperature      = 87.6;  // LAr Temperature in K
  *check if there can be a replacement later for the product, not needed now
  *double timepitch        = fdriftvelocity*ftimetick;                  //time sample (cm) */

  /** Get Geometry */
  art::ServiceHandle<geo::Geometry> geo;
  //art::ServiceHandle<util::LArProperties> larprop;
  unsigned int planes = geo->Nplanes();
  /**TODO GET THESE FROM A PARAMETERSET/DATABASE*/
  ftimetick      =  0.198; //get from parameterset
  fdriftvelocity =  0.157;  //get from paramtereset 9either k and V) 
  //fdriftvelocity = larprob->DriftVelocity(Efield_SI,Temperature);
  
  /**Get Clusters*/
  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fclusterModuleLabel,clusterListHandle);


  art::PtrVector<recob::Hit> hitlist_helpwire; /** Define hitlist of all hits per wire*/
  art::PtrVector<recob::Hit> hitlist_helpplane; /** Define hitlist to sort by #wire in hitlist_helpwire*/
  std::vector<art::PtrVector<recob::Hit> > hitlist_plane; // Define vector of hitslists for each Induction and Collection planes

  /**Loop to fill the std::vector<edm:PtrVector<recob::Hit>> with the hits of each plane*/
  for(unsigned int iplane = 0; iplane < planes; ++iplane){//Loop over planes
    geo::View_t view = geo->Plane(iplane).View();

    hitlist_helpplane.clear();
    for(unsigned int iwire=0; iwire<fwires[iplane];++iwire){//Loop over wires
      hitlist_helpwire.clear();
      for(int iclust = 0; iclust < clusterListHandle->size(); iclust++){//Loop over cluster
        art::Ptr<recob::Cluster> clust(clusterListHandle, iclust);
        if(view == clust->View()){//If correct plane
          hitlist_helpplane = clust->Hits(iplane,iwire);//Fill with Hits from the wire (hitlist_help is automatically sorted)
        }
      for (int ihits =0; ihits<hitlist_helpplane.size(); ++ihits){//work around since clust->Hits return a std::vector,recob::Hit*> object
        art::Ptr<recob::Hit> hit_help ;
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
  ftree_shwf->Fill();
}

// ***************** //
void shwf::ShowerReco::AngularDistribution(std::vector<art::PtrVector<recob::Hit> > hitlist_all){ 

  /** Determine the angular distribution of the shower energy in each plane */
  art::ServiceHandle<geo::Geometry> geo;
  unsigned int planes = geo->Nplanes();
  unsigned int iplane_hit,iwire_hit,ichannel_hit;
  float ftime_hit;
  double a_polar,b_polar,theta_polar;

  for(unsigned int iplane = 0; iplane < planes; ++iplane){/**Loop over planes*/
    art::PtrVector<recob::Hit> hitlist = hitlist_all[iplane];
    for (int ihits=0;ihits<hitlist.size();++ihits){
      art::Ptr<recob::Hit> hit_help ;
      hit_help = hitlist[ihits]; 
      const recob::Hit* hit = hit_help.get(); /** Retrieve info of the hits*/
      ftime_hit = hit->PeakTime(); /** Hit crossing time*/
      ichannel_hit = hit->Wire()->RawDigit()->Channel();
      geo->ChannelToWire(ichannel_hit, iplane_hit, iwire_hit); /** Get Information of wire/plane of the hit */

      /** Determine the vertex*/
      if (ihits==0) {    
        fwire_vertex.push_back(iwire_hit); 
        ftime_vertex.push_back(ftime_hit);
        }

      /** moving to polar coordinates (cm,cm coordinate) */
      b_polar = (iwire_hit - fwire_vertex[iplane])* fmean_wire_pitch[iplane]+ fmean_wire_pitch[iplane]; /**in cm*/
      a_polar = (ftime_hit - ftime_vertex[iplane])* ftimetick *fdriftvelocity; /** in cm*/ 
      theta_polar = asin(a_polar/sqrt(pow(a_polar,2)+pow(b_polar,2))); /**in rad*/
      theta_polar = 180*theta_polar/fpi; /** in deg*/

      /** Filling the histo (angle, energy of the hit)*/
      fh_theta[iplane]->Fill(theta_polar, hit_help->Charge()); /** angle in deg (cm,cm) coordinate*/
      }

    std::cout << "Plane " << iplane << "VertexWire= " << fwire_vertex[iplane] << "   Time= " << ftime_vertex[iplane] << std::endl;

    }
  return;
}

// ***************** //
void shwf::ShowerReco::FitAngularDistributions(){
  /** Fit function of the angular distribution (cm,cm)*/
  art::ServiceHandle<geo::Geometry> geo;
  unsigned int planes = geo->Nplanes();
  TF1 *gau = new TF1("gaus","gaus",-60, 60);

  for(unsigned int iplane = 0; iplane < planes; ++iplane){
    fh_theta[iplane]->Fit("gaus","QR");        /** Fit of the angular distribution*/
    ftheta_Mean[iplane] = gau->GetParameter(1);/** Mean value of the fit */
    ftheta_RMS[iplane] = gau->GetParameter(2); /** RMS of the fit of the angular distribution in deg*/
  }
}

// ***************** //
void shwf::ShowerReco::Get3Daxis(){
  art::ServiceHandle<geo::Geometry> geo;
  unsigned int planes = geo->Nplanes();
  double a_polar,b_polar,theta_polar;
  float wire_help,time_help;
  double slope[2];         /** Need only 2 * 2d information to reconstruct 3d for line in 3d*/
  double intercept[2];     /** with form time = slope * wire + interceptin cm, cm */

  for(unsigned int iplane = 0; iplane < planes; ++iplane){/**Loop over planes*/
    /** moving to polar coordinates (cm,cm coordinate) */
    b_polar = fwire_vertex[iplane]* fmean_wire_pitch[iplane]+ fmean_wire_pitch[iplane]; /**in cm*/
    a_polar = ftime_vertex[iplane]* ftimetick *fdriftvelocity; /** in cm*/ 
    theta_polar = asin(a_polar/sqrt(pow(a_polar,2)+pow(b_polar,2))); /**in rad*/
    theta_polar = 180*theta_polar/fpi; /** in deg*/

    /**Calculate the shower cone direction in each plane time=slope[plane]*wire + intercept[plane]*/
    if(geo->Plane(iplane).SignalType() == geo::kCollection){
      slope[1] = tan((theta_polar*fpi/180));
      wire_help = fwire_vertex[iplane] *fmean_wire_pitch[iplane]; // in cm 
      time_help = ftime_vertex[iplane] *ftimetick *fdriftvelocity; // in cm
      intercept[1] = time_help-(slope[1]*fwire_vertex[iplane]);
      }
    else if (geo->Plane(iplane).SignalType() == geo::kInduction) {
      slope[0] = tan((theta_polar*fpi/180));
      wire_help = fwire_vertex[iplane] *fmean_wire_pitch[iplane]; // in cm 
      time_help = ftime_vertex[iplane] *ftimetick *fdriftvelocity; // in cm
      intercept[0] = time_help-(slope[0]*fwire_vertex[iplane]);
      }
    }

  /**Calculate 3d axis*/
  double l(0),m(0),n(0);
  double angle_rad = 30* fpi /180;
  l = 1;
  m = (1/(2*sin(angle_rad)))*((1/slope[1])-(1/slope[0]));
  n = (1/(2*cos(angle_rad)))*((1/slope[1])+(1/slope[0]));
 
  fphi   = atan(n/l);
  ftheta = acos(m/(sqrt(pow(l,2)+pow(m,2)+pow(n,2))));

  fphi = fphi<0. ? fphi+TMath::Pi() : fphi ; /** solve the ambiguities due to tangent periodicity*/

  fphi = fphi*180/fpi;
  ftheta = ftheta*180/fpi;
  std::cout << "Shower axis: fphi=" <<fphi << "   ftheta=" << ftheta <<std::endl;
  return;
}

// ***************** //
void shwf::ShowerReco::LongTransEnergy(std::vector<art::PtrVector<recob::Hit> > hitlist_all){

  art::ServiceHandle<geo::Geometry> geo;
  unsigned int planes = geo->Nplanes();
  unsigned int iplane_hit,iwire_hit,ichannel_hit;
  float ftime_hit;
  double a_polar,b_polar,theta_polar, wire_rot, time_rot;
  double low_th, high_th, theta_sh, totcharge;   /** thresholds for adding up only the energy of the shower in a cone of falpha*RMS the angular distribution of the shower itself*/

  for(unsigned int iplane = 0; iplane < planes; ++iplane){/**Loop over planes*/
    totcharge =0.;
    art::PtrVector<recob::Hit> hitlist = hitlist_all[iplane];
    for (int ihits=0;ihits<hitlist.size();++ihits){
      art::Ptr<recob::Hit> hit_help ;
      hit_help = hitlist[ihits]; 
      const recob::Hit* hit = hit_help.get(); /** Retrieve info of the hits*/
      ftime_hit = hit->PeakTime(); /** Hit crossing time*/
      ichannel_hit = hit->Wire()->RawDigit()->Channel();
      geo->ChannelToWire(ichannel_hit, iplane_hit, iwire_hit); /** Get Information of wire/plane of the hit */

      b_polar = (iwire_hit - fwire_vertex[iplane])* fmean_wire_pitch[iplane]+ fmean_wire_pitch[iplane]; /**in cm*/
      a_polar = (ftime_hit - ftime_vertex[iplane])* ftimetick *fdriftvelocity; /** in cm*/ 
      theta_polar = asin(a_polar/sqrt(pow(a_polar,2)+pow(b_polar,2))); /**in rad*/
      theta_polar = 180*theta_polar/fpi; /** in deg*/

      low_th  = ftheta_Mean[iplane]-(falpha*ftheta_RMS[iplane]);
      high_th = ftheta_Mean[iplane]+(falpha*ftheta_RMS[iplane]);

      if((theta_polar<low_th) || (theta_polar>high_th)) continue;

      if( (theta_polar>(ftheta_Mean[iplane]-0.1))&&(theta_polar<(ftheta_Mean[iplane]+0.1)) ){
        fLastWire[iplane] = iwire_hit;
        fLastTime[iplane] = ftime_hit; 
      }

      theta_sh  = ftheta_Mean[iplane]*fpi/180; // in radians
      wire_rot = (iwire_hit - fwire_vertex[iplane])* fmean_wire_pitch[iplane] *cos(theta_sh) + (ftime_hit - ftime_vertex[iplane])* ftimetick *fdriftvelocity *sin(theta_sh);
      time_rot = (ftime_hit - ftime_vertex[iplane])* ftimetick *fdriftvelocity /cos(theta_sh) - (iwire_hit - fwire_vertex[iplane])* fmean_wire_pitch[iplane] *sin(theta_sh) 
                 - (ftime_hit - ftime_vertex[iplane])* ftimetick *fdriftvelocity *tan(theta_sh) *sin(theta_sh);   
    
      totcharge +=hit->Charge(); // Sum the energy of all the hits
    
      fsh_nrg[iplane]->Fill(wire_rot, hit->Charge()); // Fill histo longitudinal distr of the energy (IND)
      fsh_Tnrg[iplane]->Fill(time_rot, hit->Charge());// Fill histo transverse   distr of the energy (IND)
      fsh_long_hit[iplane]->Fill(wire_rot, 1);
    
    }
    ftotCharge.push_back(totcharge);
    std::cout << "In plane" << iplane << "ftotCharge (Mips)= " << totcharge<< std::endl;
  }
  return;
}

// ***************** //
void shwf::ShowerReco::Get2Dvariables(){
  art::ServiceHandle<geo::Geometry> geo;
  unsigned int planes = geo->Nplanes();
  for(unsigned int iplane = 0; iplane < planes; ++iplane){/**Loop over planes*/  
    fslope_2d.push_back((fLastTime[iplane]-ftime_vertex[iplane]) /(fLastWire[iplane]-fwire_vertex[iplane])); 
    fintercept_2d.push_back(ftime_vertex[iplane] - fwire_vertex[iplane]*fslope_2d[iplane]);
    std::cout << "Cone in plane " << iplane << "slope " << fslope_2d[iplane] << "intercept " << fintercept_2d[iplane] << std::endl;
    }
  return;
}
