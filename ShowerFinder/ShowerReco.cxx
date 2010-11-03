////////////////////////////////////////////////////////////////////////
//
// ShowerReco class - V.0.1 - 13.10.2010
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
#include "Simulation/LArVoxelList.h"
#include "Simulation/LArVoxelID.h"
#include "Simulation/LArVoxelData.h"
#include "ShowerFinder/ShowerReco.h"
#include "Geometry/geo.h"
#include "Geometry/WireGeo.h"
#include "Geometry/VolumeUtility.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Track.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Shower.h"
#include "RecoBase/Cluster.h"


static bool hit_sort_2d(const recob::Hit *h1, const recob::Hit *h2){
  return h1->Wire()->RawDigit()->Channel() < h2->Wire()->RawDigit()->Channel();
  }




//-------------------------------------------------
shwf::Shower::Shower(edm::ParameterSet const& pset) : 


  fclusters (pset.getParameter< std::string >("./clusters")),
  fshowers  (pset.getParameter< std::string >("./shower"))
  
{
  produces< std::vector<recob::Shower> >();
}


//-------------------------------------------------
shwf::Shower::~Shower()
{
}

void shwf::Shower::endJob()
{
}


//-------------------------------------------------
void shwf::Shower::beginJob(edm::EventSetup const&)
{
   edm::Service<edm::TFileService> tfs;
   edm::Service<geo::Geometry> geo;


//TODO change name of histograms so it works well...

  // Histos for the angular distribution of the shower
  char tit_h_theta[128] = {0};
  for(int p=0;p<2;p++){
    sprintf(&tit_h_theta[0],"h_theta_%i",p);
    h_theta[p] = new TH1F(tit_h_theta,"Theta distribution",180,-180., 180.);
  }
  char tit_h_theta_wt[128] = {0};
  for(int p=0;p<2;p++){
    sprintf(&tit_h_theta[0],"theta_wire_%i",p);
    h_theta_wt[p] = new TH1F(tit_h_theta_wt,"Theta wire distribution",45,-180., 180.);
  }
  
  // Histos for the longitudinal energy distribution of the shower 
  char sh_tit[128] = {0};
  for(int p=0;p<2;p++){
    sprintf(&sh_tit[0],"sh_nrg1_%i",p);
    sh_nrg[p] = new TH1F(sh_tit,"energy reco",240,0.,240*0.4);
  }
  
  //Histo for the Transverse energy distribution of the shower
  char shT_tit[128] = {0};
  for(int p=0;p<2;p++){
    sprintf(&shT_tit[0],"shT_nrg1_%i",p);
    sh_Tnrg[p] = new TH1F(shT_tit,"energy reco",80,-40.,40.);
  }
  
  //Histo for the Transverse HIT distribution of the shower
  char sh_long_tit[128] = {0};
  for(int p=0;p<2;p++){
    sprintf(&sh_long_tit[0],"sh_long_hit_%i",p);                   //96 = 240*0.4 to avoid error in make command
    sh_long_hit[p] = new TH1F(sh_long_tit,"longitudinal hit reco",96, 0., 240*0.4);
  }




}



//------------------------------------------------------------------------------------//
void shwf::Shower::produce(edm::Event& evt, edm::EventSetup const&)
{ 

  //ART CONVERSION
  //TODO HOW TO GET THE reco::Objects

  // Get Geometry
  edm::Service<geo::Geometry> geom;
  // TPC parameters
  TString tpcName = geom->GetLArTPCVolumeName();
  //geo::VolumeUtility* m_tpcVolumeUtility = new geo::VolumeUtility( tpcName );
  
  //TPC dimensions
  //  double m_TPCHalfZ = m_tpcVolumeUtility->GetHalfZ();

  // Parameters temporary defined here, but possibly to be retrieved somewhere in the code
  double plane_pitch      =  0.4;   //wire plane pitch in cm 
  double Efield_SI        =  0.7;     // Electric Field between Shield and Induction planes in kV/cm
  double Efield_IC        =  0.9;     // Electric Field between Induction and Collection planes in kV/cm
  double Temperature      = 87.6;  // LAr Temperature in K
  double driftvelocity_SI = DriftVelocity(Efield_SI,Temperature);    //drift velocity between shield and induction (cm/us)
  double driftvelocity_IC = DriftVelocity(Efield_IC,Temperature);    //drift velocity between induction and collection (cm/us)
  double timepitch        = driftvelocity*timetick;                  //time sample (cm) 
  double tSI              = plane_pitch/driftvelocity_SI/timetick;   //drift time between Shield and Collection planes (time samples)
  double tIC              = (plane_pitch/driftvelocity_IC/timetick); //drift time between Induction and Collection planes (time samples)

  
  
  //Get Clusters
  edm::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fclusters,clusterListHandle);

  edm::Handle< std::vector<recob::Hit> > hitcol;
  evt.getByLabel(fhits,hitcol);


  edm::PtrVector<recob::Hit> hitlist; // Define hitlist
  edm::PtrVector<recob::Hit>  hitlistInd, hitlistCol; // Define hitlist of Induction and Collection plane


  
 
  for(int iclust = 0; iclust < clusterListHandle->size(); iclust++)
    {
    edm::Ptr<recob::Cluster> clust(clusterListHandle, iclust);//Get cluster information


    hitlist = clust->Hits( 0, -1); // Assign the hitlist of the cluster to "hitlist variable" for INDUCTION plane

  /*  std::sort(hitlist.begin(), hitlist.end(), hit_sort_2d); //sort hit by wire    


   /* 
    //loop to fill "hitlistInd" that contain info only from INDUCTION plane hits
    int p(0),w(0), c(0); //c=channel, p=plane, w=wire

    for(int ihits = 0; ihits < hitlist.size(); ihits++)
      {
      
      c = hitlist[ihits]->Wire()->RawDigit()->Channel();
      geom->ChannelToWire(c,p,w);
      hitlistInd.push_back(ihits); 
    }
    
    hitlist = clust->Hits( 1, -1);  // Assign the hitlist of the cluster to "hitlist variable" for COLLECTION plane
    std::sort(hitlist.begin(), hitlist.end(), hit_sort_2d); //sort hit by wire  
    
    //loop to fill "hitlistCol" that contain info only from COLLECTION plane hits

    for(int ihits = 0; ihits < hitlist.size(); ihits++)
      {
      
      c = hitlist[ihits]->Wire()->RawDigit()->Channel();
      geom->ChannelToWire(c,p,w);
      hitlistCol.push_back(ihits); 
    }
    
    std::sort(hitlistInd.begin(), hitlistInd.end(), hit_sort_2d); //sort induction  hit by wire
    std::sort(hitlistCol.begin(), hitlistCol.end(), hit_sort_2d); //sort collection hit by wire 
*/

  } // end loop on DBscan cluster 
  
 
 /*
 
  
  // Save enegry in a file
  myfile.open ("/data/larsoft/releases/thomas/shower_energy.txt", std::ios::app);

  AngularDistributionI(hitlistInd, geom);
  AngularDistributionC(hitlistCol, geom);
  FitAngularDistributions();
  Get3Daxis(theta_Mean[0], theta_Mean[1], wireI1, wire1C, time1C);

  float McReleased =0.;
  McReleased = McReleasedEnergy(evt);
  if (McReleased>0) {  myfile << "MC  " << McReleased << "   "; }
  else {  myfile << "MC  " << 0.00 << "   "; }
  LongTransEnergyI(hitlistInd, geom);
  LongTransEnergyC(hitlistCol, geom);
  Get2Dvariables(wireI1, wire1C, timeI1, time1C);

  WriteHistos(evt);

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
  vertex[0] = wireI1/0.4;
  vertex[1] = wire1C/0.4;
  vertex[2] = timeI1/(timetick*driftvelocity);
  vertex[3] = time1C/(timetick*driftvelocity);

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

// ***************** //
float shwf::Shower::McReleasedEnergy(const edm::EventHandle& evt){
  //Thomas
  //Retrieve the total energy released by the MC simulation 

  float m_totenergy=0.;

  typedef std::vector<const sim::LArVoxelList* > larVoxelLists_t;
  larVoxelLists_t larVoxelLists;

  try {
    evt.DetSim().Get("",larVoxelLists);
  }
  catch ( edm::Exception& exception )
    {
      std::cout<< "No LarVoxelList objects found in DetSim() folder!" 
                << std::endl;
      std::cerr << __FILE__ << ", line " << __LINE__
		<< ": error in getting LArVoxelList objects DetSim() folder:"
		<< std::endl
		<< "    Error in " << exception.fFile << ", line " << exception.fLine 
		<< ": " << exception.fId
		<< std::endl;
      // Do no more processing for this event.
      return -1;
    }


  typedef std::vector<const sim::ParticleList*> particleLists_t;
  particleLists_t particleLists;
  
  // Read in a vector of BlorgList pointers. Note that your objects 
  // may not be in the DetSim folder; look up the location!
  try {
    evt.DetSim().Get("",particleLists);
    }
  catch(edm::Exception& exception)
    {
      std::cout<< "No particleList objects found in DetSim() folder!" 
               << std::endl;
      std::cerr << __FILE__ << ", line " << __LINE__
		<< ": error in getting particleList objects DetSim() folder:"
		<< std::endl
		<< "    Error in " << exception.fFile << ", line " << exception.fLine 
		<< ": " << exception.fId
		<< std::endl;
      // Do no more processing for this event.
      return -1;
    }

  for( particleLists_t::const_iterator k = particleLists.begin();  k != particleLists.end(); ++k )
    {
      const sim::ParticleList* particleList = (*k);

     for( larVoxelLists_t::const_iterator l = larVoxelLists.begin();  l != larVoxelLists.end(); ++l )
       {

         const sim::LArVoxelList* larVoxelList = (*l);
   
      for ( sim::LArVoxelList::const_iterator j = larVoxelList->begin(), end = larVoxelList->end(); j != end; ++j )
	{
	  // A LArVoxelList is a map, so each entry is a pair, with a
	  // first and second member.
	  const sim::LArVoxelID& voxelID = (*j).first;
	  const sim::LArVoxelData& data  = (*j).second;

          const sim::LArVoxelData& voxelData = larVoxelList->at( voxelID );

          int numberParticles = voxelData.NumberParticles();

          for ( int i = 0; i != numberParticles; ++i )
            {
              // Perhaps this is all you want: the deposited in the voxel by that particle.
              //double energy = voxelData.Energy(i);
  
              // If you want to know more about the particle:
              // Get the simulation track ID number of the particle.
              int trackID = voxelData.TrackID(i);
 
              // Fetch the particle object.
              const sim::Particle* particle = particleList->at( trackID );

              // Now you can access detailed particle information. For example:
              int pdg = particle->PdgCode();
             
              if (abs(pdg)==11)  m_totenergy += voxelData.Energy();
              //std::cout<< "voxelData.Energy()" << voxelData.Energy() << "data.Energy()" << data.Energy() << std::endl;
 
            }
          // But we're not yet done! There's typically "unassigned" energy in a voxel,
          // the sum of the particles whose energies were too small to include in the 
          // simulation output.
          double unassignedEnergy = voxelData.UnassignedEnergy();
          // ... process the left-over energy

          //std::cout<< "unassignedEnergy" << unassignedEnergy << std::endl;
	} // For each voxel

    } // For each LArVoxelList
}
  return m_totenergy;



}





// ***************** //
void shwf::Shower::AngularDistributionI(edm::PtrVector<recob::Hit> *hitlistInd,  geo::Geometry *geom){ 

  // Angular distribution of the energy of the shower - Induction plane
  int   loopI = 0; // Flag
  for(std::vector<const recob::Hit*>::iterator hitIterInd = hitlistInd.begin(); hitIterInd != hitlistInd.end();  hitIterInd++){
    recob::Hit* theHit_I = (recob::Hit*)(*hitIterInd); // Retrieve info of the hits
    
    time_I = theHit_I->CrossingTime(); // Hit crossing time
    //time_I -= presamplings;
    recob::Wire* theWire_I = (recob::Wire*)theHit_I->Wire(); // Retrive info from the Wire
    channel_I = theWire_I->RawDigit()->Channel();
    geom->ChannelToWire(channel_I, plane, wire_I);
    
    //if(wire_I<61)continue;
    
    // Determine the vertex
    if(loopI==0){
      wireI1 = wire_I;
      timeI1 = time_I;
      //wireI1 = 62;
      //timeI1 = 1045;
    } 
    
    // moving to polar coordinates (cm,cm coordinate)
    BI = (wire_I - wireI1)*0.4+0.4; //in cm
    AI = (time_I - timeI1)*timetick*0.158; // in cm 
    thetaI = asin(AI/sqrt(pow(AI,2)+pow(BI,2))); //in rad
    thetaI = 180*thetaI/3.14; // in deg
    
    // Filling the histo (angle, energy of the hit)
    h_theta[0]->Fill(thetaI, theHit_I->MIPs()); // angle in cm,cm coordinate
    
    // moving to polar coordinates (wire, tick coordinate)
    // int BI_wt = (wire_I - wireI1)+1; // in wire
    //int AI_wt = (time_I - timeI1);   // in ticks 
    //float thetaI_wt = asin(AI_wt/sqrt(pow(AI_wt,2)+pow(BI_wt,2))); //in rad
    //thetaI_wt = 180*thetaI_wt/3.14; // in deg
    // Filling the histo (angle, energy of the hit)
    //h_theta_wt[0]->Fill(thetaI_wt, theHit_I->MIPs());//angle in wire,tick coordinate
    
    loopI++; // flag increment
  }
  
   std::cout << "VertexWireI= " << wireI1 << "   VerTimeC= " << timeI1 << std::endl;
}
  

// ******************************* //

// Angular distribution of the energy of the shower - Collection view
void shwf::Shower::AngularDistributionC(edm::PtrVector<recob::Hit> *hitlistCol,  geo::Geometry *geom){
  
  int    loopC = 0; // flag
  
  for(std::vector<const recob::Hit*>::iterator hitIterCol = hitlistCol.begin(); hitIterCol != hitlistCol.end();  hitIterCol++){
    recob::Hit* theHit_C = (recob::Hit*)(*hitIterCol);
    time_C = theHit_C->CrossingTime();  
    //time_C -= (presamplings+10.1);
    recob::Wire* theWire_C = (recob::Wire*)theHit_C->Wire();
    channel_C = theWire_C->RawDigit()->Channel();
    geom->ChannelToWire(channel_C, plane, wire_C);
   
    if(loopC==0){
      wire1C = wire_C;
      time1C = time_C;
      //wire1C = 95;
      //time1C = 1060; 
    }

    //    std::cout << "VertexWireC= " << wire1C << "   VerTimeC= " << time1C << std::endl;
    
    // moving to polar coordinate
    BC = (wire_C - wire1C)*0.4 + 0.4; // in cm
    AC = (time_C - time1C)*timetick*0.158; //in cm 
    thetaC = asin(AC/sqrt(pow(AC,2)+pow(BC,2)));
    thetaC = 180*thetaC/3.14;
   
    h_theta[1]->Fill(thetaC, theHit_C->MIPs()); // Filling the histo (angle, energy of the hit)
    
    // moving to polar coordinates (wire, tick coordinate)
    //BC = (wire_C - wire1C); //in cm
    //AC = (time_C - time1C); // in cm 
    //thetaC = asin(AC/sqrt(pow(AC,2)+pow(BC,2))); //in rad
    //thetaC = 180*thetaC/3.14; // in deg
    // Filling the histo (angle, energy of the hit)
    //h_theta_wt[1]->Fill(thetaC, theHit_C->MIPs());//angle in wire,tick coordinate
    
    loopC++; // flag counter
  }
  std::cout << "VertexWireC= " << wire1C << "   VerTimeC= " << time1C << std::endl;

  return (void)0;
}
 

// ***************** //

void shwf::Shower::FitAngularDistributions(){
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

void shwf::Shower::LongTransEnergyI(edm::PtrVector<recob::Hit> *hitlistInd,  geo::Geometry *geom){
  
  // Longitudinal energy of the shower (roto-translation) - Induction plane
  double thetaI_sh;
  double wireI_rot, timeI_rot; // New coordinate after roto-transl
  double wireI_cm, timeI_cm;   // Wire and time info in cm
  double totInrg   = 0;        // tot enegry of the shower in induction
  int    loop_nrgI = 0;        // flag
  int    alpha     = 8;        // parameter (how many RMs (of the anglular distribution) is large the cone of the shower)
  double Low_thI, High_thI; // thresholds for adding up only the energy of the shower in a cone of alpha*RMS the angular distribution of the shower itself

  for(std::vector<const recob::Hit*>::iterator hitIterInd = hitlistInd.begin(); hitIterInd != hitlistInd.end();  hitIterInd++){
    recob::Hit* theHit_I = (recob::Hit*)(*hitIterInd);
    time_I = theHit_I->CrossingTime() ;  
    //time_I -= presamplings;
    recob::Wire* theWire_I = (recob::Wire*)theHit_I->Wire();
    channel_I = theWire_I->RawDigit()->Channel();
    geom->ChannelToWire(channel_I, plane, wire_I);
    
    //    if(wire_I>205)continue;
    //wireI1 = 68;
    //timeI1 = 480;
    wireI_cm = wire_I * 0.4; //in cm
    if(loop_nrgI==0)wireI1 = wireI1 * 0.4; //in cm
    timeI_cm = time_I *timetick*driftvelocity; //im cm
    if(loop_nrgI==0)timeI1 = timeI1 *timetick*driftvelocity; //in cm
 
    // moving to polar coordinates
    BI = (wireI_cm - wireI1)+0.4; //in cm
    AI = (timeI_cm - timeI1); // in cm 
    thetaI = asin(AI/sqrt(pow(AI,2)+pow(BI,2)));
    thetaI = 180*thetaI/3.14; // in deg
    //std::cout << " WireI1=" << wireI1 << " BI= " << BI << "    ThetaI = " << thetaI <<std::endl;
    
    Low_thI  = theta_Mean[0]-(alpha*theta_RMS[0]);
    High_thI = theta_Mean[0]+(alpha*theta_RMS[0]);
       
    //std::cout << "thresholds= " << Low_thI << "   " << High_thI << std::endl;
    //    if((thetaI<Low_thI) || (thetaI>High_thI)){
      //std::cout << " skipping hits outside 5*RMS out of the cone" << loop_nrgI++ << std::endl; 
      //continue;
      //}
    
    if( (thetaI>(theta_Mean[0]-0.1))&&(thetaI<(theta_Mean[0]+0.1)) ){
      LastWire[0] = (int)wire_I;
      LastTime[0] = (int)time_I; 
    }
    thetaI_sh  = theta_Mean[1]*3.14/180; // in radians
    wireI_rot = (wireI_cm - wireI1)*cos(thetaI_sh) + (timeI_cm-timeI1)*sin(thetaI_sh);
    
    timeI_rot = timeI_cm*(1/cos(thetaI_sh)) - timeI1*(1/cos(thetaI_sh)) +(wireI1 -wireI_cm)*sin(thetaI_sh) + (timeI1-timeI_cm)*tan(thetaI_sh)*sin(thetaI_sh);   
    
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
void shwf::Shower::LongTransEnergyC(edm::PtrVector<recob::Hit> *hitlistCol,  geo::Geometry *geom){
  // alogorithm for energy vs dx of the shower (roto-translation) COLLECTION VIEW
  double thetaC_sh, wireC_rot, timeC_rot, wireC_cm, timeC_cm;
  int loop_nrgC = 0;
  double Low_thC, High_thC;
  
  double totCnrg   = 0; // tot enegry of the shower in collection
  for(std::vector<const recob::Hit*>::iterator hitIterCol = hitlistCol.begin(); hitIterCol != hitlistCol.end();  hitIterCol++){
    recob::Hit* theHit_C = (recob::Hit*)(*hitIterCol);
    time_C = theHit_C->CrossingTime() ;  
    //time_C -= (presamplings+10.1);
    recob::Wire* theWire_C = (recob::Wire*)theHit_C->Wire();
    channel_C = theWire_C->RawDigit()->Channel();
    geom->ChannelToWire(channel_C, plane, wire_C);

    wireC_cm = wire_C * 0.4; //in cm
    if(loop_nrgC==0)wire1C = wire1C * 0.4; //in cm
    timeC_cm = time_C *timetick*driftvelocity; //in cm
    if(loop_nrgC==0)time1C = time1C *timetick*driftvelocity; //in cmm
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
       
    //std::cout << "thresholds= " << Low_thI << "   " << High_thI << std::endl;
    if((thetaC<Low_thC) || (thetaC>High_thC)) continue;
     
    if( (thetaC>(theta_Mean[1]-0.1))&&(thetaC<(theta_Mean[1]+0.1)) ){
      LastWire[1] = wire_C;
      LastTime[1] = time_C;
    }
    thetaC_sh = theta_Mean[1]*3.14/180;    
    wireC_rot = (wireC_cm - wire1C)*cos(thetaC_sh) + (timeC_cm-time1C)*sin(thetaC_sh);  
    
    timeC_rot = timeC_cm*(1/cos(thetaC_sh)) - time1C*(1/cos(thetaC_sh)) +(wire1C -wireC_cm)*sin(thetaC_sh) + (time1C-timeC_cm)*tan(thetaC_sh)*sin(thetaC_sh);   
    
    //std::cout << "sin thetaC= " << sin(thetaC_sh) << " costC=" << cos(thetaC_sh) << std::endl;
    //std::cout << "wireC_rot  " << wireC_rot << "  " << thetaC_sh << std::endl;
    
    //    if(wire_C>155&&wire_C<170)std::cout << "Mips[" << wire_C << "]=" << theHit_C->MIPs() <<std::endl;
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
int shwf::Shower::Get3Daxis(double thetaI, double thetaC, double Wire_vertexI, double Wire_vertexC, double Time_vertex){
  
  double timetick = 0.198;    //time sample in microsec
 
  // Making the two lines in the two views
  // time=a*wire + b
 
  slope[0] = tan((thetaI*pi/180));
  slope[1] = tan((thetaC*pi/180));

  Wire_vertexI = Wire_vertexI *0.4; // in cm 
  Wire_vertexC = Wire_vertexC *0.4; // in cm
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
int shwf::Shower::Get2Dvariables(double Wire_vertexI_wt, double Wire_vertexC_wt, double Time_I_wt, double Time_C_wt){
  
  Wire_vertexI_wt = Wire_vertexI_wt /0.4;
  Wire_vertexC_wt = Wire_vertexC_wt /0.4;
  Time_I_wt = Time_I_wt /(driftvelocity*timetick);
  Time_C_wt = Time_C_wt /(driftvelocity*timetick);
 
  //std::cout << " WvertexI=" <<Wire_vertexI_wt << "  WvertexC=" << Wire_vertexC_wt <<std::endl;
  //std::cout << " LastTIme[0]=" <<LastTime[0] << "  TimeVertexI=" << Time_I_wt <<std::endl;
  //std::cout << " LastTIme[1]=" <<LastTime[1] << "  TimeVertexC=" << Time_C_wt <<std::endl;

  // Making the two lines in the two views
  // time=a*wire + b
  slope_wt[0] = (LastTime[0]-Time_I_wt)/(LastWire[0]-Wire_vertexI_wt);
  slope_wt[1] = (LastTime[1]-Time_C_wt)/(LastWire[1]-Wire_vertexC_wt);
 
  intercept_wt[0] = Time_I_wt - Wire_vertexI_wt*slope_wt[0];
  intercept_wt[1] = Time_C_wt - Wire_vertexC_wt*slope_wt[1];

  //  slope_wt[0] = tan((thetaI_wt*pi/180));
  //slope_wt[1] = tan((thetaC_wt*pi/180));
  
  //intercept_wt[0] = Time_I_wt - (slope_wt[0]*Wire_vertexI_wt);
  //intercept_wt[1] = Time_C_wt - (slope_wt[1]*Wire_vertexC_wt);
  
  //  std::cout << " Slope_wt[0]=" <<slope_wt[0] << "  Intercept_wt[0]=" << intercept_wt[0] <<std::endl;
  //std::cout << " Slope_wt[1]=" <<slope_wt[1] << "  Intercept_wt[1]=" << intercept_wt[1] <<std::endl;  

  return 0;
}


// automatically called by Get3Daxis
void shwf::Shower::GetPitchLength(double theta, double phi){
  
//TODO replace with correctly call

  const static double wire_pitch    =  0.4;   // wire pitch in cm

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
double shwf::Shower::DriftVelocity(double Efield, double Temperature){
  
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

void shwf::Shower::WriteHistos(edm::EventHandle& evt){
  
  char Out_file[128] = "Shower_0000";
  sprintf(&Out_file[7],"%.4i.root", evt.Header().Run());
  sprintf(&Out_file[7+4],".root");
  
  TFile O_file(Out_file, "update" );
  std::cout << Out_file << std::endl;
  if (O_file.IsZombie()) {
    std::cout << "Error opening file" << std::endl;
    exit(-1);
  }  
  //h_theta[0]->Write(); 
  //h_theta[1]->Write(); 
  ////  h_theta_wt[0]->Write(); 
  ////h_theta_wt[1]->Write(); 
  //sh_nrg[0]->Write();
  //sh_nrg[1]->Write();
  //sh_Tnrg[0]->Write();
  //sh_Tnrg[1]->Write();
  //sh_long_hit[0]->Write();
  //sh_long_hit[1]->Write();
}
 

*/
