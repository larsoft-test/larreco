////////////////////////////////////////////////////////////////////////
//
// ShowerReco class - V.0.2 - 11/04/2010 
//
// biagio.rossi@lhep.unibe.ch   (FWMK)
// thomas.strauss@lhep.unibe.ch (ART)
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

//-------------------------------------------------
shwf::ShowerReco::ShowerReco(edm::ParameterSet const& pset) : 

  fclusterModuleLabel (pset.getParameter< std::string >("clusters"))
  
{
  produces< std::vector<recob::Shower> >();
}

//-------------------------------------------------
shwf::ShowerReco::~ShowerReco()
{
}

//-------------------------------------------------
void shwf::ShowerReco::beginJob(edm::EventSetup const&)
{
  // Get Geometry
  edm::Service<geo::Geometry> geo;
  planes = geo->Nplanes();
  
  double mean_wirepitch;

  for (unsigned int i=0;i<planes-1;++i){
    mean_wirepitch = 0.;
    plane_pitch.push_back(geo->PlanePitch(i,i+1));
    wires.push_back(geo->Nwires(i));
    for (unsigned int j=0;j<wires[i]-1;++j){
      wire_pitch.push_back(geo->WirePitch(j,j+1,i));
      mean_wirepitch+=geo->WirePitch(j,j+1,i);
    }
    mean_wire_pitch.push_back(mean_wirepitch/geo->Nwires(i));
  }
  //Get TFileService and define output Histograms
  edm::Service<edm::TFileService> tfs;

  // Create Histos names
  char tit_h_theta[128] = {0};
  char tit_h_theta_wt[128] = {0};
  char sh_long_tit[128] = {0};
  char sh_tit[128] = {0};
  char shT_tit[128] = {0};
  
//p equals number of planes, make vector of histograms, so p=0 to max_planes below
  int nbins;
  for(unsigned int i=0;i<planes-1;++i){

    //Histos for the angular distribution of the shower
    sprintf(&tit_h_theta[0],"h_theta_%i",i);
    h_theta.push_back(tfs->make<TH1F>(tit_h_theta,"Theta distribution",180,-180., 180.));

    sprintf(&tit_h_theta[0],"theta_wire_%i",i);
    h_theta_wt.push_back(tfs->make<TH1F>(tit_h_theta_wt,"Theta wire distribution",45,-180., 180.));
 
    // Histos for the longitudinal energy distribution of the shower 
    sprintf(&sh_tit[0],"sh_nrg1_%i",i);                   //number of wires used  //min wire    //max_wire you need for the anlysis
    sh_nrg.push_back(tfs->make<TH1F>(sh_tit,"energy reco",wires[i],               0.,           wires[i]*mean_wire_pitch[i]));

    //Histo for the transverse energy distribution of the shower
    sprintf(&shT_tit[0],"shT_nrg1_%i",i);                   //units are ticks most lickely, but how did biagio get it???
    sh_Tnrg.push_back(tfs->make<TH1F>(shT_tit,"energy reco",80,-40.,40.));
  
    //Histo for the Transverse HIT distribution of the shower
    nbins = (int)(wires[i]*mean_wire_pitch[i]);
    sprintf(&sh_long_tit[0],"sh_long_hit_%i",i);                                   //min wire    //max_wire you need for the anlysis
    sh_long_hit.push_back(tfs->make<TH1F>(sh_long_tit,"longitudinal hit reco",nbins, 0.,         wires[i]*mean_wire_pitch[i]));
  }
}

//------------------------------------------------------------------------------------//
void shwf::ShowerReco::produce(edm::Event& evt, edm::EventSetup const&)
{ 

  //TODO THIS VALUES SHOULD BE PARAMETERS OF THE MODULE, evtl from database
  //double Efield_SI        =  0.7;     // Electric Field between Shield and Induction planes in kV/cm
  //double Efield_IC        =  0.9;     // Electric Field between Induction and Collection planes in kV/cm
  //double Temperature      = 87.6;  // LAr Temperature in K
  //check if there can be a replacement later for the product, not needed now
  //double timepitch        = driftvelocity*timetick;                  //time sample (cm) 

  // Get Geometry
  edm::Service<geo::Geometry> geo;

  timetick      =  0.198; //get from parameterset
  driftvelocity =  0.157;  //get from paramtereset 9either k and V) 
  //driftvelocity = DriftVelocity(Efield_SI,Temperature);
  
  //Get Clusters
  edm::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fclusterModuleLabel,clusterListHandle);


  edm::PtrVector<recob::Hit> hitlist_helpwire; // Define hitlist of all hits per wire
  edm::PtrVector<recob::Hit> hitlist_helpplane; // Define hitlist of all hits per plane
  std::vector<edm::PtrVector<recob::Hit> > hitlist_plane; // Define vector of hitls for each Induction and Collection planes
  int i,ii;

  for(unsigned int iplane = 0; iplane < geo->Nplanes(); ++iplane){//Loop over planes
    geo::View_t view = geo->Plane(i).View();
    i = iplane;
    hitlist_helpplane.clear();
    for(unsigned int iwire=0; iwire<wires[i];++iwire){//Loop over wires
      ii = iwire;
      hitlist_helpwire.clear();
      for(int iclust = 0; iclust < clusterListHandle->size(); iclust++){//Loop over cluster
        edm::Ptr<recob::Cluster> clust(clusterListHandle, iclust);
        if(view == clust->View()){//If correct plane
          hitlist_helpplane = clust->Hits(i,ii);//Fill with Hits from the wire (hitlist_help is automatically sorted)
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

/*
 
  
  // Save enegry in a file
  myfile.open ("/data/larsoft/releases/thomas/shower_energy.txt", std::ios::app);

  AngularDistributionI(hitlistInd);
  AngularDistributionC(hitlistCol);
  FitAngularDistributions();
  Get3Daxis(theta_Mean[0], theta_Mean[1], wire_vertex, wire_vertex, time_vertex);
  LongTransEnergyI(hitlistInd);
  LongTransEnergyC(hitlistCol);
  Get2Dvariables(wire_vertex, wire_vertex, time_vertex, time_vertex);

//line 275 - Use the Geometry to define the number of planes, ie geo::Geometry::Nplanes(), rather than hardcoding it to 2 14

  double dir2D[4] = {0};
  dir2D[0] = slope_wt[0];
  dir2D[1] = intercept_wt[0]; //in ticks
  dir2D[2] = slope_wt[1];
  dir2D[3] = intercept_wt[1]; //in ticks
  
  // 3D directions
  double dir3D[4] = {0};
  dir3D[0] = slope[0];
  dir3D[1] = intercept[0]; //in cm
  dir3D[0] = slope[1]; 
  dir3D[1] = intercept[1]; //in cm
 
  // Shower vertex
  double vertex[4] = {0};
  vertex[0] = wire_vertex/wirepitch;
  vertex[1] = wire_vertex/wirepitch;
  vertex[2] = time_vertex/(timetick*driftvelocity);
  vertex[3] = time_vertex/(timetick*driftvelocity);

  // Shower END 
  double sh_end[4] = {0};
  sh_end[0] = LastWire[0];
  sh_end[1] = LastWire[1];
  sh_end[2] = LastTime[0];
  sh_end[3] = LastTime[1];

  //  std::cout << "LastWire  " << LastWire[0]<<"   time="<<LastTime[1]<<std::endl;
  
  //std::vector<recob::Shower *> Shower3DVector(2);  //holds the 3D shower to be saved 
  for(int p=0;p<2;p++)Shower3DVector[p] = new recob::Shower();
  geo::View_t PlaneI = geo::kU;
  Shower3DVector[0]->SetPitch( Pitch[0], PlaneI);
  geo::View_t PlaneC = geo::kV;
  Shower3DVector[1]->SetPitch( Pitch[1], PlaneC);
  
  // Assign 2D slopes and intercepts in wire, ticks
  Shower3DVector[0]->fDir2D[0] = dir2D[0];
  Shower3DVector[0]->fDir2D[1] = dir2D[1];
  Shower3DVector[0]->fDir2D[2] = dir2D[2];
  Shower3DVector[0]->fDir2D[3] = dir2D[3];

  // Assign 3D slopes and intercepts in cm, cm
  Shower3DVector[0]->fDir3D[0] = dir3D[0];
  Shower3DVector[0]->fDir3D[1] = dir3D[1];
  Shower3DVector[0]->fDir3D[2] = dir3D[2];
  Shower3DVector[0]->fDir3D[3] = dir3D[3];

  // Assign vertwx coordinate
  Shower3DVector[0]->fVertex[0] = vertex[0];
  Shower3DVector[0]->fVertex[1] = vertex[1];
  Shower3DVector[0]->fVertex[2] = vertex[2];
  Shower3DVector[0]->fVertex[3] = vertex[3];

  // Assign lasthit of the shower coordinate
  Shower3DVector[0]->fLastWire[0] = LastWire[0];
  Shower3DVector[0]->fLastWire[1] = LastWire[1];
  Shower3DVector[0]->fLastTime[0] = LastTime[0];
  Shower3DVector[0]->fLastTime[1] = LastTime[1];


  */
}

/*



//BRIAN PASS THE CORRECT HISTGRAMS AND COMBINE INTO ONLY ONE FUNCTION< INSTEAD OF COLL+ind.

// ***************** //
void shwf::ShowerReco::AngularDistributionI(edm::PtrVector<recob::Hit> hitlistInd){ 

  // Angular distribution of the energy of the shower - Induction plane
  edm::Service<geo::Geometry> geo;
  int   loopI = 0; // Flag
  for(std::vector<const recob::Hit>::iterator hitIterInd = hitlistInd.begin(); hitIterInd != hitlistInd.end();  hitIterInd++){
    recob::Hit* theHit_I = (recob::Hit*)(*hitIterInd); // Retrieve info of the hits

//BRIAN edm::Ptr<recob::Hit> = *hitIterInd    

    time_hit = theHit_I->CrossingTime(); // Hit crossing time

//BRIAN geo->ChannelToWire(theHit_I->WIre()->RawDigit()->Channel(), plane, geom->wire_I)


    recob::Wire* thewire_hit = (recob::Wire*)theHit_I->Wire(); // Retrive info from the Wire
    channel_hit = thewire_hit->RawDigit()->Channel();
    geo->ChannelToWire(channel_hit, plane_hit, wire_hit);
    
    //if(wire_hit<61)continue;
    
    // Determine the vertex
    if(loopI==0){
      wire_vertex = wire_hit;
      time_vertex = time_hit;
      //wire_vertex = 62;
      //time_vertex = 1045;
    } 
    
    // moving to polar coordinates (cm,cm coordinate)
    b_polar = (wire_hit - wire_vertex)*wirepitch+wirepitch; //in cm
    AI = (time_hit - time_vertex)*timetick*driftvelocity; // in cm 
    theta_polar = asin(AI/sqrt(pow(AI,2)+pow(b_polar,2))); //in rad
    theta_polar = 180*theta_polar/3.14; // in deg
    
    // Filling the histo (angle, energy of the hit)
    h_theta[0]->Fill(theta_polar, theHit_I->MIPs()); // angle in cm,cm coordinate
    
    // moving to polar coordinates (wire, tick coordinate)
    // int b_polar_wt = (wire_hit - wire_vertex)+1; // in wire
    //int AI_wt = (time_hit - time_vertex);   // in ticks 
    //float thetaI_wt = asin(AI_wt/sqrt(pow(AI_wt,2)+pow(b_polar_wt,2))); //in rad
    //thetaI_wt = 180*thetaI_wt/3.14; // in deg
    // Filling the histo (angle, energy of the hit)
    //h_theta_wt[0]->Fill(thetaI_wt, theHit_I->MIPs());//angle in wire,tick coordinate
    
    loopI++; // flag increment
  }
  
   std::cout << "VertexWireI= " << wire_vertex << "   VerTimeC= " << time_vertex << std::endl;
}
  

// ******************************* //

// Angular distribution of the energy of the shower - Collection view
void shwf::ShowerReco::AngularDistributionC(edm::PtrVector<recob::Hit> hitlistCol){
  edm::Service<geo::Geometry> geo;
  int    loopC = 0; // flag
  
  for(std::vector<const recob::Hit*>::iterator hitIterCol = hitlistCol.begin(); hitIterCol != hitlistCol.end();  hitIterCol++){
    recob::Hit* theHit_C = (recob::Hit*)(*hitIterCol);
    time_hit = theHit_C->CrossingTime();  

    recob::Wire* thewire_hit = (recob::Wire*)theHit_C->Wire();
    channel_hit = thewire_hit->RawDigit()->Channel();
    geo->ChannelToWire(channel_hit, plane_hit, wire_hit);
   
    if(loopC==0){
      wire_vertex = wire_hit;
      time_vertex = time_hit;
      //wire_vertex = 95;
      //time_vertex = 1060; 
    }

    //    std::cout << "VertexWireC= " << wire_vertex << "   VerTimeC= " << time_vertex << std::endl;
    
    // moving to polar coordinate
    b_polar = (wire_hit - wire_vertex)*wirepitch + wirepitch; // in cm
    AC = (time_hit - time_vertex)*timetick*driftvelocity; //in cm 
    thetaC = asin(AC/sqrt(pow(AC,2)+pow(b_polar,2)));
    thetaC = 180*thetaC/3.14;
   
    h_theta[1]->Fill(thetaC, theHit_C->MIPs()); // Filling the histo (angle, energy of the hit)
    
    // moving to polar coordinates (wire, tick coordinate)
    //b_polar = (wire_hit - wire_vertex); //in cm
    //AC = (time_hit - time_vertex); // in cm 
    //thetaC = asin(AC/sqrt(pow(AC,2)+pow(b_polar,2))); //in rad
    //thetaC = 180*thetaC/3.14; // in deg
    // Filling the histo (angle, energy of the hit)
    //h_theta_wt[1]->Fill(thetaC, theHit_C->MIPs());//angle in wire,tick coordinate
    
    loopC++; // flag counter
  }
  std::cout << "VertexWireC= " << wire_vertex << "   VerTimeC= " << time_vertex << std::endl;

  return (void)0;
}
 

// ***************** //

void shwf::ShowerReco::FitAngularDistributions(){
  // Fit function of the angular distribution (cm,cm)
  TF1 *gau = new TF1("gaus","gaus",-60, 60);
  h_theta[0]->Fit("gaus","QR"); // Fit of the angular distribution
  theta_Mean[0] = gau->GetParameter(1);// Mean value of the fit (Induction)
  theta_RMS[0] = gau->GetParameter(2); // RMS of the fit of the angular distribution (IND) in deg
 
  h_theta[1]->Fit("gaus","QR");
  theta_Mean[1] = gau->GetParameter(1);// Mean value of the fit (Collection)
  theta_RMS[1]  = gau->GetParameter(2);// RMS of the fit of the angular distribution (COL) in deg
 
  myfile << theta_RMS[0] << "    " <<theta_RMS[1];
  //std::cout << "Ind Theta(cm,cm)=" << theta_Mean[0] << " RMS=" << theta_RMS[0] <<std::endl;
  //std::cout << "Col theta(cm,cm)=" << theta_Mean[1] << " RMS=" << theta_RMS[1] <<std::endl;

  // Fit function of the angular distribution (wire,tick)
  //TF1 *gau_wt = new TF1("gaus_wt","gaus",-120, 120);
  //double MAX = h_theta_wt[0]->GetMaximumBin();
  //std::cout<< " MAX " << MAX << std::endl;
  //h_theta_wt[0]->Fit("gaus_wt","QR"); // Fit of the angular distribution
  //theta_wt_Mean[0] = gau_wt->GetParameter(1);// Mean value of the fit (Induction)
  //theta_wt_RMS[0]  = gau_wt->GetParameter(2); // RMS of the fit of the angular distribution (IND) in deg
 
  //h_theta_wt[1]->Fit("gaus_wt","QR");
  //theta_wt_Mean[1] = gau_wt->GetParameter(1);// Mean value of the fit (Collection)
  //theta_wt_RMS[1]  = gau_wt->GetParameter(2);// RMS of the fit of the angular distribution (COL) in deg

  //std::cout << "Ind Theta(w,t)=" << theta_wt_Mean[0] << " RMS=" << theta_wt_RMS[0] <<std::endl;
  //std::cout << "Col theta(w,t)=" << theta_wt_Mean[1] << " RMS=" << theta_wt_RMS[1] <<std::endl;
  
}
// ***************** //

void shwf::ShowerReco::LongTransEnergyI(edm::PtrVector<recob::Hit> hitlistInd){
  
  // Longitudinal energy of the shower (roto-translation) - Induction plane
  edm::Service<geo::Geometry> geo;
  double thetaI_sh;
  double wireI_rot, timeI_rot; // New coordinate after roto-transl
  double wireI_cm, timeI_cm;   // Wire and time info in cm
  double totInrg   = 0;        // tot enegry of the shower in induction
  int    loop_nrgI = 0;        // flag
  int    alpha     = 8;        // parameter (how many RMs (of the anglular distribution) is large the cone of the shower)
  double Low_thI, High_thI; // thresholds for adding up only the energy of the shower in a cone of alpha*RMS the angular distribution of the shower itself

  for(std::vector<const recob::Hit*>::iterator hitIterInd = hitlistInd.begin(); hitIterInd != hitlistInd.end();  hitIterInd++){
    recob::Hit* theHit_I = (recob::Hit*)(*hitIterInd);
    time_hit = theHit_I->CrossingTime() ;  

    recob::Wire* thewire_hit = (recob::Wire*)theHit_I->Wire();
    channel_hit = thewire_hit->RawDigit()->Channel();
    geo->ChannelToWire(channel_hit, plane_hit, wire_hit);
    
    //    if(wire_hit>205)continue;
    //wire_vertex = 68;
    //time_vertex = 480;
    wireI_cm = wire_hit * wirepitch; //in cm
    if(loop_nrgI==0)wire_vertex = wire_vertex * wirepitch; //in cm
    timeI_cm = time_hit *timetick*driftvelocity; //im cm
    if(loop_nrgI==0)time_vertex = time_vertex *timetick*driftvelocity; //in cm
 
    // moving to polar coordinates
    BI = (wireI_cm - wire_vertex)+wirepitch; //in cm
    AI = (timeI_cm - time_vertex); // in cm 
    theta_polar = asin(AI/sqrt(pow(AI,2)+pow(BI,2)));
    theta_polar = 180*theta_polar/3.14; // in deg
    //std::cout << " wire_vertex=" << wire_vertex << " BI= " << BI << "    theta_polar = " << theta_polar <<std::endl;
    
    Low_thI  = theta_Mean[0]-(alpha*theta_RMS[0]);
    High_thI = theta_Mean[0]+(alpha*theta_RMS[0]);
       
    //std::cout << "thresholds= " << Low_thI << "   " << High_thI << std::endl;
    //    if((theta_polar<Low_thI) || (theta_polar>High_thI)){
      //std::cout << " skipping hits outside 5*RMS out of the cone" << loop_nrgI++ << std::endl; 
      //continue;
      //}
    
    if( (theta_polar>(theta_Mean[0]-0.1))&&(theta_polar<(theta_Mean[0]+0.1)) ){
      LastWire[0] = (int)wire_hit;
      LastTime[0] = (int)time_hit; 
    }
    thetaI_sh  = theta_Mean[1]*3.14/180; // in radians
    wireI_rot = (wireI_cm - wire_vertex)*cos(thetaI_sh) + (timeI_cm-time_vertex)*sin(thetaI_sh);
    
    timeI_rot = timeI_cm*(1/cos(thetaI_sh)) - time_vertex*(1/cos(thetaI_sh)) +(wire_vertex -wireI_cm)*sin(thetaI_sh) + (time_vertex-timeI_cm)*tan(thetaI_sh)*sin(thetaI_sh);   
    
    totInrg +=theHit_I->MIPs(); // Sum the energy of all the hits
    
    sh_nrg[0]->Fill(wireI_rot, theHit_I->MIPs()); // Fill histo longitudinal distr of the energy (IND)
    sh_Tnrg[0]->Fill(timeI_rot, theHit_I->MIPs());// Fill histo transverse   distr of the energy (IND)
    
    sh_long_hit[0]->Fill(wireI_rot, 1);
    
    loop_nrgI++;
    
  }

  myfile << "   " << totInrg << "   ";

  std::cout << "TotInrg(GeV)= " << totInrg<< std::endl;
  
  //std::cout << " Last_wire Ind= " << LastWire[0] << " Last_timeI  " << LastTime[0]<< std::endl;
}

// -------------------------- //
void shwf::ShowerReco::LongTransEnergyC(edm::PtrVector<recob::Hit> hitlistCol){
  // alogorithm for energy vs dx of the shower (roto-translation) COLLECTION VIEW

  int    alpha     = 8;        // parameter (how many RMs (of the anglular distribution) is large the cone of the shower)
 
 edm::Service<geo::Geometry> geo;
  double thetaC_sh, wireC_rot, timeC_rot, wireC_cm, timeC_cm;
  int loop_nrgC = 0;
  double Low_thC, High_thC;
  
  double totCnrg   = 0; // tot enegry of the shower in collection
  for(std::vector<const recob::Hit*>::iterator hitIterCol = hitlistCol.begin(); hitIterCol != hitlistCol.end();  hitIterCol++){
    recob::Hit* theHit_C = (recob::Hit*)(*hitIterCol);
    time_hit = theHit_C->CrossingTime() ;  

    recob::Wire* thewire_hit = (recob::Wire*)theHit_C->Wire();
    channel_hit = thewire_hit->RawDigit()->Channel();
    geo->ChannelToWire(channel_hit, plane_hit, wire_hit);

    wireC_cm = wire_hit * wirepitch; //in cm
    if(loop_nrgC==0)wire_vertex = wire_vertex * wirepitch; //in cm
    timeC_cm = time_hit *timetick*driftvelocity; //in cm
    if(loop_nrgC==0)time_vertex = time_vertex *timetick*driftvelocity; //in cmm
    //std::cout << "Vertex wire_vertex= " << wire_vertex << " VerTimeC=" << time_vertex << std::endl;
    //std::cout << "wire_hit= " << wireC_cm << " TimeC=" << timeC_cm << std::endl;
    
    // moving to polar coordinates
    b_polar = (wireC_cm - wire_vertex)+wirepitch; //in cm
    AC = (timeC_cm - time_vertex); // in cm 
    thetaC = asin(AC/sqrt(pow(AC,2)+pow(b_polar,2)));
    thetaC = 180*thetaC/3.14; // in deg
    //std::cout << " wire_vertex=" << wire_vertex << " b_polar= " << b_polar << "    theta_polar = " << theta_polar <<std::endl;
    
    Low_thC  = theta_Mean[1]-(alpha*theta_RMS[1]);
    High_thC = theta_Mean[1]+(alpha*theta_RMS[1]);
       
    //std::cout << "thresholds= " << Low_thI << "   " << High_thI << std::endl;
    if((thetaC<Low_thC) || (thetaC>High_thC)) continue;
     
    if( (thetaC>(theta_Mean[1]-0.1))&&(thetaC<(theta_Mean[1]+0.1)) ){
      LastWire[1] = wire_hit;
      LastTime[1] = time_hit;
    }
    thetaC_sh = theta_Mean[1]*3.14/180;    
    wireC_rot = (wireC_cm - wire_vertex)*cos(thetaC_sh) + (timeC_cm-time_vertex)*sin(thetaC_sh);  
    
    timeC_rot = timeC_cm*(1/cos(thetaC_sh)) - time_vertex*(1/cos(thetaC_sh)) +(wire_vertex -wireC_cm)*sin(thetaC_sh) + (time_vertex-timeC_cm)*tan(thetaC_sh)*sin(thetaC_sh);   
    
    //std::cout << "sin thetaC= " << sin(thetaC_sh) << " costC=" << cos(thetaC_sh) << std::endl;
    //std::cout << "wireC_rot  " << wireC_rot << "  " << thetaC_sh << std::endl;
    
    //    if(wire_hit>155&&wire_hit<170)std::cout << "Mips[" << wire_hit << "]=" << theHit_C->MIPs() <<std::endl;
    totCnrg +=theHit_C->MIPs(); // Sum the energy of all the hits
    sh_Tnrg[1]->Fill(timeC_rot,theHit_C->MIPs());
    sh_nrg[1]->Fill(wireC_rot, theHit_C->MIPs());
    //std::cout << " wireCrot= " << wireC_rot << "  hitC= " << theHit_C->MIPs()<<std::endl;
    sh_long_hit[1]->Fill(wireC_rot, 1);
    loop_nrgC++;
  }

  std::cout << "TotCnrg= " << totCnrg << "  in GeV="  <<(((((totCnrg/7.3)*6300)*100)/60)*23.6)*pow(10,-9)<< std::endl;
  myfile << "   " << totCnrg << std::endl;
  myfile.close(); //std::cout << " Last_wire Coll= " << LastWire[1] << " Last_timeC  " << LastTime[1] << std::endl; 

}

//--------------------------------------------------//


//------------------------------------------------------------------------------------//  
int shwf::ShowerReco::Get3Daxis(double theta_polar, double thetaC, double Wire_vertexI, double Wire_vertexC, double Time_vertex){
  
  double timetick = 0.198;    //time sample in microsec
 
  // Making the two lines in the two views
  // time=a*wire + b
 
  slope[0] = tan((theta_polar*pi/180));
  slope[1] = tan((thetaC*pi/180));

  Wire_vertexI = Wire_vertexI *wirepitch; // in cm 
  Wire_vertexC = Wire_vertexC *wirepitch; // in cm
  Time_vertex  = Time_vertex  *timetick*driftvelocity; // in cm

  intercept[0] = Time_vertex - (slope[0]*Wire_vertexI);
  intercept[1] = Time_vertex - (slope[1]*Wire_vertexC);

  // std::cout << " Slope[0]=" <<slope[0] << "  Intercept[0]=" << intercept[0] <<std::endl;
  ///std::cout << " Slope[1]=" <<slope[1] << "  Intercept[1]=" << intercept[1] <<std::endl;  

  double l(0),m(0),n(0);
 
  double angle_rad = 30* pi /180;

  l = 1;
  m = (1/(2*sin(angle_rad)))*((1/slope[1])-(1/slope[0]));
  n = (1/(2*cos(angle_rad)))*((1/slope[1])+(1/slope[0]));
 
  // Director angles
  //double theta(0), phi(0); // Director angles
  phi   = atan(n/l);
  theta = acos(m/(sqrt(pow(l,2)+pow(m,2)+pow(n,2))));
  
  std::cout << "befor Phi=" <<phi*180/pi  << "   theta=" << theta*180/pi <<std::endl;
  phi = phi<0. ? phi+TMath::Pi() : phi ; // solve the ambiguities due to tangent periodicity
  //if(phi<0)phi=-1*phi;
  //if(phi>0){
    //phi = phi*180/pi;
    //phi = 180 - phi;
    //phi = phi*pi/180;
    //}
  std::cout << " Phi=" <<phi*180/pi  << "   theta=" << theta*180/pi <<std::endl;
  myfile << "   " <<phi*180/pi << "   ";
  myfile << "   " << theta*180/pi << "   ";
 
  GetPitchLength(theta, phi);

  return 0;
}

//------------------------------------------------------------------------------------//  
int shwf::ShowerReco::Get2Dvariables(double Wire_vertexI_wt, double Wire_vertexC_wt, double time_hit_wt, double time_hit_wt){
  
  Wire_vertexI_wt = Wire_vertexI_wt /wirepitch;
  Wire_vertexC_wt = Wire_vertexC_wt /wirepitch;
  time_hit_wt = time_hit_wt /(driftvelocity*timetick);
  time_hit_wt = time_hit_wt /(driftvelocity*timetick);
 
  //std::cout << " WvertexI=" <<Wire_vertexI_wt << "  WvertexC=" << Wire_vertexC_wt <<std::endl;
  //std::cout << " LastTIme[0]=" <<LastTime[0] << "  TimeVertexI=" << time_hit_wt <<std::endl;
  //std::cout << " LastTIme[1]=" <<LastTime[1] << "  TimeVertexC=" << time_hit_wt <<std::endl;

  // Making the two lines in the two views
  // time=a*wire + b
  slope_wt[0] = (LastTime[0]-time_hit_wt)/(LastWire[0]-Wire_vertexI_wt);
  slope_wt[1] = (LastTime[1]-time_hit_wt)/(LastWire[1]-Wire_vertexC_wt);
 
  intercept_wt[0] = time_hit_wt - Wire_vertexI_wt*slope_wt[0];
  intercept_wt[1] = time_hit_wt - Wire_vertexC_wt*slope_wt[1];

  //  slope_wt[0] = tan((theta_polar_wt*pi/180));
  //slope_wt[1] = tan((thetaC_wt*pi/180));
  
  //intercept_wt[0] = time_hit_wt - (slope_wt[0]*Wire_vertexI_wt);
  //intercept_wt[1] = time_hit_wt - (slope_wt[1]*Wire_vertexC_wt);
  
  //  std::cout << " Slope_wt[0]=" <<slope_wt[0] << "  Intercept_wt[0]=" << intercept_wt[0] <<std::endl;
  //std::cout << " Slope_wt[1]=" <<slope_wt[1] << "  Intercept_wt[1]=" << intercept_wt[1] <<std::endl;  

  return 0;
}


// automatically called by Get3Daxis
void shwf::ShowerReco::GetPitchLength(double theta, double phi){
  
//TODO replace with correctly call


  //std::cout << "theta_input=" <<theta << " Phi=" << phi <<std::endl;  

  double angle_rad = 30* pi /180;

  Pitch[1] = wire_pitch/((cos(angle_rad)*sin(theta)*sin(phi))+(sin(angle_rad)*cos(theta)));
  Pitch[0] = wire_pitch/((cos(angle_rad)*sin(theta)*sin(phi))-(sin(angle_rad)*cos(theta)));

  std::cout << "IND-Pitch[0]=" <<Pitch[0] << " COL-Pitch[1]=" << Pitch[1] <<std::endl;   

  myfile << Pitch[0] << "   ";
  myfile << Pitch[1] << "   ";


  return (void)0;
}


//-----------------------------------------------//
void GetVertex(){



  return (void)0;
}
 

//-----------------------------------------------//
double shwf::ShowerReco::DriftVelocity(double Efield, double Temperature){
  
  // Dirft Velocity as a function of Electric Field and LAr Temperature
  // from : W. Walkowiak, NIM A 449 (2000) 288-294
  
  double vd;
  double field = Efield; //Electric field kV/cm   default = 0.5
  double T = Temperature;
  
  
  double P1,P2,P3,P4,P5,P6,T0;
  P1=-0.01481; // K^-1
  P2=-0.0075; // K^-1
  P3=0.141;//(kV/cm)^-1
  P4=12.4;//kV/cm
  P5=1.627;//(kV/cm)^-P6
  P6=0.317;
  T0 = 90.371; // K

  vd=(P1*(T-T0)+1)*(P3*field*TMath::Log(1+P4/field) + P5*TMath::Power(field,P6))+P2*(T-T0);

  vd /= 10.;

  return vd;// in cm/us
}

 

*/
