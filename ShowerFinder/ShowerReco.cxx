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


shwf::ShowerReco::ShowerReco(fhicl::ParameterSet const& pset)
{
  this->reconfigure(pset);
  produces< std::vector<recob::Shower> >();
}

void shwf::ShowerReco::reconfigure(fhicl::ParameterSet pset) 
{
  fClusterModuleLabel = pset.get< std::string >("ClusterModuleLabel");
  fShwrOutput = pset.get< std::string >("ShowerOutputTxtFile");
}

// ***************** //
shwf::ShowerReco::~ShowerReco()
{
}

namespace shwf {
struct SortByWire 
{
  bool operator() (recob::Hit const& h1, recob::Hit const& h2) const 
  { return 
      h1.Wire()->RawDigit()->Channel() < 
      h2.Wire()->RawDigit()->Channel() ;
  }
};
}

// ***************** //
void shwf::ShowerReco::beginJob()
{
  /** Get Geometry*/
  art::ServiceHandle<geo::Geometry> geo;
  unsigned int planes = geo->Nplanes();
  fmean_wire_pitch = geo->WirePitch(0,1,0);    //wire pitch in cm

  /**Get TFileService and define output Histograms*/
  art::ServiceHandle<art::TFileService> tfs;

  /** Create Histos names*/
  char tit_dedx[128] = {0};
  char tit_h_theta[128] = {0};
  char tit_h_theta_wt[128] = {0};
  char sh_long_tit[128] = {0};
  char sh_tit[128] = {0};
  char shT_tit[128] = {0};
  
  int nbins;
  for(unsigned int i=0;i<planes;++i){


    //    sprintf(&tit_dedx[0],"fh_dedx_%.4i_%.4i_%i",i);
    sprintf(&tit_dedx[0],"fh_dedx_._%i",i);
    fh_dedx[i] = tfs->make<TH1F>(tit_dedx,"dEdx vs distance from vertex",120,0, 40);

    /**Histos for the angular distribution theta of the shower*/
    sprintf(&tit_h_theta[0],"fh_theta_%i",i);
    fh_theta[i] = tfs->make<TH1F>(tit_h_theta,"Theta distribution",720,-180., 180.);

    /**Histos for the angular distribution theta wire of the shower*/
    sprintf(&tit_h_theta[0],"ftheta_wire_%i",i);
    fh_theta_wt[i] = tfs->make<TH1F>(tit_h_theta_wt,"Theta wire distribution",720,-180., 180.);
 
    /**Histos for the longitudinal energy distribution of the shower */
    sprintf(&sh_tit[0],"fsh_nrg1_%i",i);                   /**number of wires used,min wire,max_wire you need for the anlysis*/
    fsh_nrg[i] = tfs->make<TH1F>(sh_tit,"energy reco",240,0.,240*fmean_wire_pitch);

    /**Histos for the transverse energy distribution of the shower*/
    sprintf(&shT_tit[0],"fshT_nrg1_%i",i);                   /**units are ticks most lickely, but how did Biagio get size of it???*/
    fsh_Tnrg[i] = tfs->make<TH1F>(shT_tit,"energy reco",80,-40.,40.);
  
    /**Histos for the Transverse HIT distribution of the shower*/
    nbins = (int)(240*fmean_wire_pitch);
    sprintf(&sh_long_tit[0],"fsh_long_hit_%i",i);                           /**nbins,min wire,max_wire you need for the analysis*/
    fsh_long_hit[i] = tfs->make<TH1F>(sh_long_tit,"longitudinal hit reco",nbins, 0.,     240*fmean_wire_pitch);
  }

  ftree_shwf =tfs->make<TTree>("ShowerReco","Results");/**All-knowing tree with reconstruction information*/
  /*
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
  */

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
  //fdriftvelocity = larprob->DriftVelocity(Efield_SI,Temperature);
  
  /**Get Clusters*/
  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fClusterModuleLabel,clusterListHandle);

   art::PtrVector < recob::Hit> hitlistCol;
   art::PtrVector < recob::Hit> hitlistInd;

    std::auto_ptr<std::vector<recob::Shower> > Shower3DVector(new std::vector<recob::Shower>);

  for(unsigned int ii = 0; ii < clusterListHandle->size(); ++ii)
    {

      art::Ptr<recob::Cluster> cl(clusterListHandle, ii);
      
      // Figure out which View the cluster belongs to 
      
      ///      int clPlane = cl->View()-1;

      art::PtrVector<recob::Hit> hitlist;
      hitlist = cl->Hits();
      hitlist.sort(shwf::SortByWire());
      unsigned int p(0),w(0), c(0), t(0); //c=channel, p=plane, w=wire

      for(art::PtrVectorItr<recob::Hit> a = hitlist.begin(); a != hitlist.end();  a++) //loop over cluster hits
      {
	c=(*a)->Wire()->RawDigit()->Channel(); 
	geo->ChannelToWire(c,t,p,w);

	if(geo->Plane(p,t).SignalType() == geo::kCollection)
	  {	  
	    hitlistCol.push_back(*a);

	  }
	else if (geo->Plane(p,t).SignalType() == geo::kInduction)
	  {
	    hitlistInd.push_back(*a);
	  } 
      }

      hitlistInd.sort(shwf::SortByWire());
      hitlistCol.sort(shwf::SortByWire());


    } // End loop on clusters.


     
  // Save energy in a file
  myfile.open (fShwrOutput.c_str(), std::ios::app);
  
  AngularDistributionI(hitlistInd); // 2D Direction of the shower Induction
  AngularDistributionC(hitlistCol); // 2D Direction of the shower Collection
  FitAngularDistributions();              // Fit of 2d distributions (for both planes)
  Get3Daxis(theta_Mean[0], theta_Mean[1], wireI1, wire1C, time1C); 
  LongTransEnergyI(hitlistInd); //Longitudinal and Transverse energy profile of the Shower induction
  LongTransEnergyC(hitlistCol); //Longitudinal and Transverse energy profile of the Shower induction
  Get2Dvariables(wireI1, wire1C, timeI1, time1C);

  /**Fill the output tree with all information */
  //  ftree_shwf->Fill();

  // This needs work, clearly.  
  //for(int p=0;p<2;p++)Shower3DVector->push_back(shower);
  //evt.put(Shower3DVector)

}


// ***************** //
void shwf::ShowerReco::FitAngularDistributions(){
  /** Fit function of the angular distribution (cm,cm)*/
  art::ServiceHandle<geo::Geometry> geo;
  unsigned int planes = geo->Nplanes();
  TF1 *gau = new TF1("gaus","gaus",-60, 60);

  /*
  for(unsigned int iplane = 0; iplane < planes; ++iplane){
    fh_theta[iplane]->Fit("gaus","QR");        // Fit of the angular distribution
    ftheta_Mean[iplane] = gau->GetParameter(1);// Mean value of the fit 
    ftheta_RMS[iplane] = gau->GetParameter(2); // RMS of the fit of the angular distribution in deg
  }
*/
}

// ***************** //


// ***************** //
int shwf::ShowerReco::Get2Dvariables(float Wire_vertexI_wt, float Wire_vertexC_wt, float Time_I_wt, float Time_C_wt){
  
  art::ServiceHandle<geo::Geometry> geo;
  // only needed for drawing the axis of the shower in the event display

  Wire_vertexI_wt = Wire_vertexI_wt /0.4;
  Wire_vertexC_wt = Wire_vertexC_wt /0.4;
  Time_I_wt = Time_I_wt /(fdriftvelocity*ftimetick);
  Time_C_wt = Time_C_wt /(fdriftvelocity*ftimetick);
 
  // Making the two lines in the two views
  slope_wt[0] = (LastTime[0]-Time_I_wt)/(LastWire[0]-Wire_vertexI_wt);
  slope_wt[1] = (LastTime[1]-Time_C_wt)/(LastWire[1]-Wire_vertexC_wt);
 
  intercept_wt[0] = Time_I_wt - Wire_vertexI_wt*slope_wt[0];
  intercept_wt[1] = Time_C_wt - Wire_vertexC_wt*slope_wt[1];

  myfile << std::endl;
  myfile.close(); 
  return 0;
}
 
int shwf::ShowerReco::Get3Daxis(float thetaI, float thetaC, float Wire_vertexI, float Wire_vertexC, float Time_vertex){
  
  //Get theta and phi (polar angles "direction of the shower")

  //float ftimetick = 0.198;    //time sample in microsec
 
  slope[0] = tan((thetaI*pi/180));
  slope[1] = tan((thetaC*pi/180));

  Wire_vertexI = Wire_vertexI *0.4; // in cm 
  Wire_vertexC = Wire_vertexC *0.4; // in cm
  Time_vertex  = Time_vertex  *ftimetick*fdriftvelocity; // in cm

  intercept[0] = Time_vertex - (slope[0]*Wire_vertexI);
  intercept[1] = Time_vertex - (slope[1]*Wire_vertexC);

 
  float l(0),m(0),n(0);
  // Get Geometry
  art::ServiceHandle<geo::Geometry> geom;
 
  float Angle = geom->Plane(1).Wire(0).ThetaZ(false)-TMath::Pi()/2.; // wire angle with respect to the vertical direction
  std::cout << "Angle= " <<  Angle<<std::endl;
  
  float angle_rad = 30* pi /180;
  angle_rad = Angle;
  
  l = 1;
  m = (1/(2*sin(angle_rad)))*((1/slope[1])-(1/slope[0]));
  n = (1/(2*cos(angle_rad)))*((1/slope[1])+(1/slope[0]));
 
  std::cout << "________________m= " << m << std::endl;

  // Director angles
  phi   = atan(n/l);
  theta = acos(m/(sqrt(pow(l,2)+pow(m,2)+pow(n,2))));

  float Phi = phi>0. ? (TMath::Pi()/2)-phi : fabs(phi)-(TMath::Pi()/2) ; // solve the ambiguities due to tangent periodicity
  float Theta=0;
  if(Phi>0)Theta = theta-(TMath::Pi()/2);
  if(Phi<0)Theta = (TMath::Pi()/2)-theta;

  if(slope[0]==0 || slope[1]==0){
    Phi   = 0;
    Theta = 0;
    theta = 0;
    phi = 0;
  }
  
  std::cout << " Phi=" <<Phi*180/pi << "   theta=" << Theta*180/pi <<std::endl;
  myfile << "   " <<Phi*180/pi << "   ";
  myfile << "   " <<Theta*180/pi << "   ";
  
  GetPitchLength(theta, phi); //Get pitch of (two) wire planes (only for Argoneut)
  return 0;
}

int counterI =0;
void shwf::ShowerReco::LongTransEnergyI(art::PtrVector < recob::Hit> hitlistInd){
  

  // Longitudinal energy of the shower (roto-translation) - Induction plane
  float thetaI_sh;
  float wireI_rot, timeI_rot; // New coordinate after roto-transl
  float wireI_cm, timeI_cm;   // Wire and time info in cm
  float totInrg   = 0;        // tot enegry of the shower in induction
  int    loop_nrgI = 0;        // flag
  float Low_thI, High_thI; // thresholds for adding up only the energy of the shower in a cone of alpha*RMS the angular distribution of the shower itself
  IdEdx4cm = 0; // dedx|cm of the shower in Induction
  IdedxCounter = 0;
  art::ServiceHandle<geo::Geometry> geom;
 
  for(art::PtrVectorItr<recob::Hit> hitIterInd = hitlistInd.begin(); hitIterInd != hitlistInd.end();  hitIterInd++){
    art::Ptr<recob::Hit> theHit_I = (*hitIterInd);
    time_I = theHit_I->PeakTime() ;  
    //time_I -= presamplings;
    art::Ptr<recob::Wire> theWire_I = theHit_I->Wire();
    channel_I = theWire_I->RawDigit()->Channel();
    geom->ChannelToWire(channel_I, tpc, plane, wire_I);

    //   if(wire_I>218)continue;
   
    wireI_cm = wire_I * 0.4; //in cm
    if(loop_nrgI==0)wireI1 = wireI1 * 0.4; //in cm
    timeI_cm = time_I *ftimetick*fdriftvelocity; //im cm
    if(loop_nrgI==0)timeI1 = timeI1 *ftimetick*fdriftvelocity; //in cm
 
    // moving to polar coordinates
    BI = (wireI_cm - wireI1)+0.4; //in cm
    AI = (timeI_cm - timeI1); // in cm 
    thetaI = asin(AI/sqrt(pow(AI,2)+pow(BI,2)));
    thetaI = 180*thetaI/3.14; // in deg
    //std::cout << " WireI1=" << wireI1 << " BI= " << BI << "    ThetaI = " << thetaI <<std::endl;
    
    Low_thI  = theta_Mean[0]-(alpha*theta_RMS[0]);
    High_thI = theta_Mean[0]+(alpha*theta_RMS[0]);
       
    counterI++;
    if(counterI==1)std::cout << "IND-thresholds= " << Low_thI << "   " << High_thI << std::endl;
    /*  if((thetaI<Low_thI) || (thetaI>High_thI)){
      //std::cout << " skipping hits outside 5*RMS out of the cone" << loop_nrgI++ << std::endl; 
      continue;
      }*/
    
    if( (thetaI>(theta_Mean[0]-1.0))&&(thetaI<(theta_Mean[0]+1.0)) ){
      LastWire[0] = (int)wire_I;
      LastTime[0] = (int)time_I; 
    }
    thetaI_sh  = theta_Mean[0]*3.14/180; // in radians
    wireI_rot = (wireI_cm - wireI1)*cos(thetaI_sh) + (timeI_cm-timeI1)*sin(thetaI_sh);
    
    timeI_rot = timeI_cm*(1/cos(thetaI_sh)) - timeI1*(1/cos(thetaI_sh)) +(wireI1 -wireI_cm)*sin(thetaI_sh) + (timeI1-timeI_cm)*tan(thetaI_sh)*sin(thetaI_sh);   
    
    totInrg +=theHit_I->Charge(); // Sum the energy of all the hits
    
    fsh_nrg[0]->Fill(wireI_rot, theHit_I->Charge()); // Fill histo longitudinal distr of the energy (IND)
    fsh_Tnrg[0]->Fill(timeI_rot, theHit_I->Charge());// Fill histo transverse   distr of the energy (IND)
    loop_nrgI++;

    if( ((((wireI_cm-wireI1)/0.4)*Pitch[0])<10.8)&&((((wireI_cm-wireI1)/0.4)*Pitch[0])>0.8)){ 
      //      if(IdedxCounter<=10)
      IdEdx4cm+=theHit_I->Charge();//Sum the energy of all the hits
      IdedxCounter++;
      std::cout << "dedxI= " << IdEdx4cm<<std::endl;
      //if(theHit_I->Charge()<0)std::cout << "Inrg= " << theHit_I->Charge()<<std::endl;
    }
  }
  
  myfile << "   " << totInrg << "   ";
  
  std::cout<< "    TotInrg(GeV)= " << totInrg<< std::endl;
  
  //std::cout << " Last_wire Ind= " << LastWire[0] << " Last_timeI  " << LastTime[0]<< std::endl;
}

int counter=0;   
 // -------------------------- //
void shwf::ShowerReco::LongTransEnergyC(art::PtrVector < recob::Hit> hitlistCol){
  // alogorithm for energy vs dx of the shower (roto-translation) COLLECTION VIEW
  float thetaC_sh, wireC_rot, timeC_rot, wireC_cm, timeC_cm;
  int loop_nrgC = 0;
  float Low_thC, High_thC;
  
  float totCnrg = 0; // tot enegry of the shower in collection
  CdEdx4cm = 0; // tot enegry of the shower in collection
  CdedxCounter = 0;
  art::ServiceHandle<geo::Geometry> geom;

  for(art::PtrVectorItr<recob::Hit> hitIterCol = hitlistCol.begin(); hitIterCol != hitlistCol.end();  hitIterCol++){
    art::Ptr<recob::Hit> theHit_C = (*hitIterCol);
    time_C = theHit_C->PeakTime() ;  
    //time_C -= (presamplings+10.1);
    art::Ptr<recob::Wire> theWire_C = theHit_C->Wire();
    channel_C = theWire_C->RawDigit()->Channel();
    geom->ChannelToWire(channel_C, tpc, plane, wire_C);

    //    if(time_C<1020)continue;

    wireC_cm = wire_C * 0.4; //in cm
    if(loop_nrgC==0)wire1C = wire1C * 0.4; //in cm
    timeC_cm = time_C *ftimetick*fdriftvelocity; //in cm
    if(loop_nrgC==0)time1C = time1C *ftimetick*fdriftvelocity; //in cmm
    //std::cout << "Vertex wire1C= " << wire1C << " VerTimeC=" << time1C << std::endl;
    //std::cout << "wire_C= " << wireC_cm << " TimeC=" << timeC_cm << std::endl;
    
    // moving to polar coordinates
    BC = (wireC_cm - wire1C)+0.4; //in cm
    AC = (timeC_cm - time1C); // in cm 
    thetaC = asin(AC/sqrt(pow(AC,2)+pow(BC,2)));
    thetaC = 180*thetaC/3.14; // in deg
    //std::cout << " WireI1=" << wireI1 << " BI= " << BI << "    ThetaI = " << thetaI <<std::endl;
    
    Low_thC  = theta_Mean[1]-(alpha*theta_RMS[1]);
    High_thC = theta_Mean[1]+(alpha*theta_RMS[1]);
 
    counter++;
    if(counter==1)std::cout << "CoLL-thresholds= " << Low_thC << "   " << High_thC << std::endl;
    //if((thetaC<Low_thC) || (thetaC>High_thC)) continue;
     
    if( (thetaC>(theta_Mean[1]-1.0))&&(thetaC<(theta_Mean[1]+1.0)) ){
      LastWire[1] = wire_C;
      LastTime[1] = time_C;
    }
    thetaC_sh = theta_Mean[1]*3.14/180;    
    wireC_rot = (wireC_cm - wire1C)*cos(thetaC_sh) + (timeC_cm-time1C)*sin(thetaC_sh);  
    
    timeC_rot = timeC_cm*(1/cos(thetaC_sh)) - time1C*(1/cos(thetaC_sh)) +(wire1C -wireC_cm)*sin(thetaC_sh) + (time1C-timeC_cm)*tan(thetaC_sh)*sin(thetaC_sh);   
    
    double lifetime = exp(-time_C*ftimetick/735);

    totCnrg +=(theHit_C->Charge()/lifetime); // Sum the energy of all the hits
    fsh_Tnrg[1]->Fill(timeC_rot,theHit_C->Charge());
    fsh_nrg[1]->Fill(wireC_rot, theHit_C->Charge());
    loop_nrgC++;
    
    if( ((((wireC_cm-wire1C)/0.4)*Pitch[1])<40.4)&&((((wireC_cm-wire1C)/0.4)*Pitch[1])>0.8)){ 
      //if(CdedxCounter<=10)
      CdEdx4cm+=(theHit_C->Charge()/lifetime);//Sum the energy of all the hits
      CdedxCounter++;
      std::cout << " COLL- hit-Vertex=" << ((wireC_cm-wire1C)/0.4)*Pitch[1] << " dedx|(MIPs/cm)=" << theHit_C->Charge()/Pitch[1]<< std::endl;
      fh_dedx[1]->Fill( ((wireC_cm-wire1C)/0.4)*Pitch[1], ((((theHit_C->Charge()/lifetime)/Pitch[1])*10/7)/7.6)*6250*23.6*pow(10,-6) )  ;


    }
    
  }
  //  std::cout << "COLECTION - dEdx|4cm= " << CdEdx4cm<< "   TotCnrg= " << totCnrg << "  in GeV="  <<(((((totCnrg/7.3)*6300)*100)/60)*23.6)*pow(10,-9)<< std::endl;
  myfile << "   " << totCnrg;
}

//------------------------------------------------------------------------------------//  


void shwf::ShowerReco::AngularDistributionI(art::PtrVector < recob::Hit>  hitlistInd){ 
  art::ServiceHandle<geo::Geometry> geom;
  // Angular distribution of the energy of the shower - Induction plane
  int   loopI = 0; // Flag
  // this should changed on the loop on the cluster of the shower
  for(art::PtrVectorItr<recob::Hit> hitIterInd = hitlistInd.begin(); hitIterInd != hitlistInd.end();  hitIterInd++){
    art::Ptr<recob::Hit> theHit_I = (*hitIterInd); // Retrieve info of the hits
    
    time_I = theHit_I->PeakTime(); // Hit crossing time
    //time_I -= presamplings;
    art::Ptr<recob::Wire> theWire_I = theHit_I->Wire(); // Retrive info from the Wire
    channel_I = theWire_I->RawDigit()->Channel();
    geom->ChannelToWire(channel_I, tpc, plane, wire_I);
    
    //Here we should take GetVertex function to retrieve the vertex of the shower
    // Determine the vertex
    if(loopI==0){
      wireI1 = wire_I;
      timeI1 = time_I; 
    } 
    
    // moving to polar coordinates (cm,cm coordinate)
    BI = (wire_I - wireI1)*0.4+0.4; //in cm
    AI = (time_I - timeI1)*ftimetick*0.158; // in cm 
    thetaI = asin( (AI/sqrt(pow(AI,2)+pow(BI,2)))); //in rad
    thetaI = 180*thetaI/3.14; // in deg
    // Filling the histo (angle, energy of the hit)
    fh_theta[0]->Fill(thetaI, theHit_I->Charge()); // angle in cm,cm coordinate
    
  
    loopI++; // flag increment
  }
  
}
  

// ******************************* //

// Angular distribution of the energy of the shower - Collection view
void shwf::ShowerReco::AngularDistributionC(art::PtrVector < recob::Hit>  hitlistCol){
  
  int    loopC = 0; // flag
  art::ServiceHandle<geo::Geometry> geom;
  // this should changed on the loop on the cluster of the shower
  for(art::PtrVectorItr<recob::Hit> hitIterCol = hitlistCol.begin(); hitIterCol != hitlistCol.end();  hitIterCol++){
    art::Ptr<recob::Hit> theHit_C = (*hitIterCol);
    time_C = theHit_C->PeakTime();  
    //time_C -= (presamplings+10.1);
    art::Ptr<recob::Wire> theWire_C = theHit_C->Wire();
    channel_C = theWire_C->RawDigit()->Channel();
    geom->ChannelToWire(channel_C, tpc, plane, wire_C);
    
    //Here we should take GetVertex function to retrieve the vertex of the shower
    if(loopC==0){
      wire1C = wire_C;
      time1C = time_C;
    }

    //    std::cout << "VertexWireC= " << wire1C << "   VerTimeC= " << time1C << std::endl;
    
    // moving to polar coordinate
    BC = (wire_C - wire1C)*0.4 + 0.4; // in cm
    AC = (time_C - time1C)*ftimetick*0.158; //in cm 
    thetaC = asin(  AC/sqrt(pow(AC,2)+pow(BC,2)) );
    thetaC = 180*thetaC/3.14;
    fh_theta[1]->Fill(thetaC, theHit_C->Charge()); // Filling the histo (angle, energy of the hit)
    
    loopC++; // flag counter
  }
  std::cout << "VertexWireC= " << wire1C << "   VerTimeC= " << time1C << std::endl;

  return (void)0;
}

// automatically called by Get3Daxis
void shwf::ShowerReco::GetPitchLength(float theta, float phi){
  
  //Get pitch of 2 wire planes 
  //for generalization to n planes, different formulas for the geometrical transofrmation are needed (this is ArgoNeuT specific)

  // Get Geometry
  art::ServiceHandle<geo::Geometry> geom;

  // TPC parameters
  TString tpcName = geom->GetLArTPCVolumeName();

  float Angle = geom->Plane(1).Wire(0).ThetaZ(false)-TMath::Pi()/2.; // wire angle with respect to the vertical direction

  std::cout << "theta_input=" <<theta << " Phi=" << phi <<std::endl;  
  
  float angle_rad = (180+30)* pi /180;
  Pitch[1] = fabs(fmean_wire_pitch/((cos(angle_rad)*sin(theta)*sin(phi))+(sin(angle_rad)*cos(theta))));
  Pitch[0] = fabs(fmean_wire_pitch/((cos(angle_rad)*sin(theta)*sin(phi))-(sin(angle_rad)*cos(theta))));

  std::cout << "IND-Pitch[0]=" <<Pitch[0] << " COL-Pitch[1]=" << Pitch[1] <<std::endl;   

  myfile << Pitch[0] << "   ";
  myfile << Pitch[1] << "   ";
  myfile << "dEdxI=" <<IdEdx4cm<< "   ";
  myfile << CdEdx4cm<< "   "; 

  std::cout << "INDUCTION  - dEdx|4cm= " << IdEdx4cm<< "   " <<IdedxCounter<<std::endl;
  std::cout << "COLLECTION - dEdx|4cm= " << CdEdx4cm<< "   " <<CdedxCounter<<std::endl;

  return (void)0;
}


/*****************************************************/
void GetVertex(){

  // here we should GetThe verteces list
  // then understand which one belongs to a shower
  // getting the cluster that contain that vertex and pass everything to the main routine for the pghysical variables determination


  return (void)0;
}
