////////////////////////////////////////////////////////////////////////
//
// KingaCluster class
//
// kinga.partyka@yale.edu
//

////////////////////////////////////////////////////////////////////////

#include "KingaCluster.h"
extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}
// ROOT includes
#include <TCanvas.h>
#include "TDatabasePDG.h"
#include "TSystem.h"

#include <sstream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <vector>
#include <dirent.h>

#include "art/Framework/Core/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Persistency/Common/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

 
#include "RawData/RawDigit.h"
#include "Filters/ChannelFilter.h"
#include "SimulationBase/simbase.h"
#include "Simulation/sim.h"
#include "Simulation/SimListUtils.h"
#include "RecoBase/recobase.h"
#include "Geometry/geo.h"
#include "Utilities/LArProperties.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"

cluster::KingaCluster::KingaCluster(fhicl::ParameterSet const& pset) :
  fDBScanModuleLabel       (pset.get< std::string >("DBScanModuleLabel"))
  //fGenieGenModuleLabel          (pset.get< std::string >("GenieGenModuleLabel"))
  //fthreshold                (pset.get< int >("Threshold")),
 
{
  produces< std::vector<recob::Cluster> >();
}

cluster::KingaCluster::~KingaCluster()
{
}

//-------------------------------------------------
void cluster::KingaCluster::beginJob(){
  // get access to the TFile service
  art::ServiceHandle<art::TFileService> tfs;
  art::ServiceHandle<geo::Geometry> geo;
  //unsigned int planes = geo->Nplanes();
 
 
  fh_theta_ind= tfs->make<TH1F>("fh_theta_ind","theta angle in degrees, Induction Plane", 180,-180 ,180  );
  fh_theta_coll= tfs->make<TH1F>("fh_theta_coll","theta angle in degrees, Collection Plane", 180,-180 ,180  );
  fh_theta_ind_2D= tfs->make<TH1F>("fh_theta_ind_2D","theta angle in degrees, Induction Plane", 180,-180 ,180  );
  fh_theta_coll_2D= tfs->make<TH1F>("fh_theta_coll_2D","theta angle in degrees, Collection Plane", 180,-180 ,180  );
 fh_theta_ind_Area= tfs->make<TH1F>("fh_theta_ind_Area","Hit Area vs theta angle in degrees, Induction Plane", 180,-180 ,180  );
  fh_theta_coll_Area= tfs->make<TH1F>("fh_theta_coll_Area","Hit Area vs theta angle in degrees, Collection Plane", 180,-180 ,180  );



 
  
  // Create Histos names
 //  char tit_h_theta[128] = {0};
//   for(unsigned int i=0;i<planes;++i){
// 
//     /**Histos for the angular distribution theta of the shower*/
//     sprintf(&tit_h_theta[0],"fh_theta_%i",i);
//     fh_theta.push_back(tfs->make<TH1F>(tit_h_theta,"Theta distribution",180,-180., 180.));
//     }

}

//-----------------------------------------------------------------

void cluster::KingaCluster::produce(art::Event& evt)
{
std::cout<<"In KingaCluster::produce(art::Event& evt)"<<std::endl;
std::cout<<" Working on ";
std::cout << "Run: " << evt.run();
std::cout << " Event: " << evt.id().event() << std::endl;
  fpeaks_found=1;
  fMC=0; 
 int RunNo= evt.run();
 int EventNo=evt.id().event();
 art::ServiceHandle<util::LArProperties> larp;
// double electronlifetime=larp->ElectronLifetime();
 




if (evt.isRealData()) 
    {
      std::cout<<" YOU ARE WORKING WITH DATA !!!!!!!!!!!!!!!!!!!!! "<<std::endl;
    //let's see if this Run has been scanned and thus we can find vertex info for it:
    DIR *pDIR=0;
    struct dirent *entry=0;
    const char path[60]="/argoneut/data/simplescan_data/simplescan_text/";
  char path_[60];
  strcpy(path_,path);
    if(pDIR==opendir(path))
    {
       while(entry==readdir(pDIR))
       {
         if(strcmp(entry->d_name,".")!=0 && strcmp(entry->d_name, "..")!=0)
         {
          std::string file=entry->d_name;
          std::cout<<"*******file is "<<file<<std::endl;
         std::string::size_type pos_end=file.rfind(".");
         std::string no_string;
         no_string=file.substr(pos_end-3,3);
         std::cout<<"no_string= "<<no_string<<std::endl;
         std::istringstream stream(no_string);
         int No;
         stream>>No;
           if(RunNo==No)
           {
            std::cout<<"RUN NO "<<RunNo<<" was scanned, we can find a vertex info :) "<<std::endl;
            //Now, open the file and find the vertex info in it:
            //.........................................
            strcat(path_,entry->d_name);
            std::string k;
            int run,event,blah1,blah2,blah3,blah4,time_ind,time_coll,w_ind,w_coll;
            std::ifstream ScannedFile;
            ScannedFile.open(path_,std::ios::in);
            if(!ScannedFile.is_open()){std::cout<<" Couldn't open file named: "<<entry->d_name<<std::endl;}
            if(ScannedFile.is_open()){std::cout<<" openED file named: "<<entry->d_name<<std::endl;} 
             
             
            while(getline(ScannedFile,k))
            {
            std::istringstream ins;
            ins.clear();
            ins.str(k);
            ins>>run>>event>>blah1>>blah2>>blah3>>blah4>>time_ind>>time_coll>>w_ind>>w_coll;
            //std::cout<<run<<" "<<event<<std::endl;
            if(EventNo==event && RunNo==run)
            {
            ftime_vertex.push_back(time_ind);
            ftime_vertex.push_back(time_coll);
            fwire_vertex.push_back(w_ind);
            fwire_vertex.push_back(w_coll);
            
            std::cout<<"GOT VERTEX INFO FROM THE SCANNED FILE:"<<std::endl;
            std::cout<<"("<<time_ind<<" , "<<w_ind<<" )"<<std::endl;
            std::cout<<"("<<time_coll<<" , "<<w_coll<<" )"<<std::endl;
            
            
            break;
            }//event=event
            
            
            
            
            
            }//while getline
            
            
            
            
           //.........................................
            ScannedFile.close();
            break; 
         
           } //if the run was scanned
           
           
           
         }
       
       
       } //while it loops through all the entries
        closedir(pDIR);
    
    
    }//if can open directory
    else{ std::cout<<" THIS RUN HASN'T BEEN SCANNED! SORRY! "<<std::endl;}
    
    } //if realData




else {

std::cout<<" YOU ARE WORKING WITH MC !!!!!!!!!!!!!!!!!!!!! "<<std::endl;
 fMC=1;
 fGenieGenModuleLabel="generator";

art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle);
art::PtrVector<simb::MCTruth> mclist;
  for (unsigned int ii = 0; ii <  mctruthListHandle->size(); ++ii)
    {
      art::Ptr<simb::MCTruth> mctparticle(mctruthListHandle,ii);
      mclist.push_back(mctparticle);
    } 
    
    
    
    
for( unsigned int i = 0; i < mclist.size(); ++i ){

    art::Ptr<simb::MCTruth> mc(mclist[i]);

    simb::MCParticle neut(mc->GetParticle(i));

    // std::cout<<"vertex: "<<neut.Nu().Vx()<<" "<<neut.Nu().Vy()<<" "<<neut.Nu().Vz()<<std::endl;
    MCvertex[0] =neut.Vx();
    MCvertex[1] =neut.Vy();
    MCvertex[2] =neut.Vz();
    std::cout<<"MCvertex[0]= "<<MCvertex[0]<<std::endl;
    std::cout<<"driftvelocity= "<<larp->DriftVelocity(larp->Efield(),larp->Temperature())<<std::endl;
    double drifttick=(MCvertex[0]/larp->DriftVelocity(larp->Efield(),larp->Temperature()))*(1./.198);
    
    std::cout<<"%%%%%%%%%%%%%%%%%%   drifttick= "<<drifttick<<std::endl;
    ftime_vertex.push_back(drifttick);
    ftime_vertex.push_back(drifttick);
    

  }






}



  //////////////////////////////////////////////////////
  // here is how to get a collection of objects out of the file
  // and connect it to a art::Handle
  //////////////////////////////////////////////////////
  // Read in the clusterList object(s).
  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fDBScanModuleLabel,clusterListHandle);
  //Point to a collection of clusters to output.
  std::auto_ptr<std::vector<recob::Cluster> > ccol(new std::vector<recob::Cluster>);


  art::ServiceHandle<geo::Geometry> geom;
  art::PtrVector<recob::Cluster> clusIn;
 
  art::PtrVector<recob::Hit> hits;
  art::PtrVector<recob::Hit> clusterHits;
 

 
  for(unsigned int ii = 0; ii < clusterListHandle->size(); ++ii)
    {
      art::Ptr<recob::Cluster> cluster(clusterListHandle, ii);
      clusIn.push_back(cluster);
    }
std::cout<<"No of DBSCAN clusters= "<<clusIn.size()<<std::endl;



unsigned int p(0),w(0), channel(0);

std::cout<<"No of planes = "<<geom->Nplanes()<<std::endl;
  for(unsigned int plane = 0; plane < geom->Nplanes(); plane++) {

    for(unsigned int j=0; j<clusIn.size();++j) {
   
    hits=clusIn[j]->Hits();
    for(unsigned int i = 0; i< hits.size(); ++i){
    channel=hits[i]->Wire()->RawDigit()->Channel();
    geom->ChannelToWire(channel,p,w);

    if(p == plane){
    allhits.push_back(hits[i]);
    //std::cout<<"plane= "<<plane<<" wire= "<<w<<" time= "<<hits[i]->PeakTime()<<std::endl; 
    }
   
    }
   // std::cout<<"hits.size()= "<<hits.size()<<std::endl;
    }
   
    std::cout<<"allhits.size()="<<allhits.size()<<std::endl;
   
   
    //Now we have hits for the plane that we are on right now, so let's do some work:
   //maxBin.clear();
   std::cout<<"ATTENTION, STARTING WORK ON PLANE# "<<plane<<std::endl;
    AngularDistribution(plane);
    FindMax(plane);
    if(fpeaks_found==0){
    std::cout<<"KingaClusters FAILED on this event because no peaks were found. Perhaps your threshold for peak's height is too big. Goodbye! "<<std::endl;
    allhits.clear();
     maxBin.clear();
     maxBinValues.clear();
     SortedMaxBin.clear();
     MaxStartPoint.clear();
     MaxEndPoint.clear();
     MaxStartPointTheta.clear();
     MaxEndPointTheta.clear();
     HitsWithClusterID.clear();
     FinalPeaks.clear();
     OriginalmaxBinValues.clear();
     for(int bin=0; bin< fh_theta_ind_2D->GetNbinsX(); bin++){
   
   fh_theta_ind_2D->SetBinContent(bin,0);
   fh_theta_coll_2D->SetBinContent(bin,0);
   fh_theta_ind->SetBinContent(bin,0);
   fh_theta_coll->SetBinContent(bin,0);
   fh_theta_coll_Area->SetBinContent(bin,0);
   fh_theta_ind_Area->SetBinContent(bin,0);
   }
     
    return;
    }
    //FinalPeaks();
    FindClusters(plane);
    std::cout<<"HitsWithClusterID.size()= "<<HitsWithClusterID.size()<< "compare with allhits.size()= "<<allhits.size()<<std::endl;
    for(unsigned int ClusterNo=0; ClusterNo<MaxStartPoint.size();ClusterNo++) {
    
       for(unsigned int j=0; j<HitsWithClusterID.size();j++){
     
       if(HitsWithClusterID[j]==(ClusterNo+1)){
       
       clusterHits.push_back(allhits[j]);
       } //if
    
    
    
    } //loop over HitsWithClusterID

// let's look at the clusters produced:
std::cout<<"For Cluster # "<<ClusterNo<<" we have "<<clusterHits.size()<<" hits :"<<std::endl;
if(ClusterNo==4){
for(unsigned int i=0; i<clusterHits.size();i++){
channel=clusterHits[i]->Wire()->RawDigit()->Channel();
    geom->ChannelToWire(channel,p,w);
std::cout<<"wire ="<<w<<"  "<<clusterHits[i]->PeakTime();

}
}
std::cout<<std::endl;


//     //.................................
    if (clusterHits.size()>0)
	    {
	      ///!todo: need to define start and end positions for this cluster and slopes for dTdW, dQdW
	      unsigned int p = 0; 
	      unsigned int sw = 0;
	      unsigned int ew = 0;
	      geom->ChannelToWire(clusterHits[0]->Wire()->RawDigit()->Channel(), p, sw);
	      geom->ChannelToWire(clusterHits[clusterHits.size()-1]->Wire()->RawDigit()->Channel(), p, ew);

	      
	      
	      recob::Cluster cluster(clusterHits, 
				     sw*1., 0.,
				     clusterHits[0]->PeakTime(), clusterHits[0]->SigmaPeakTime(),
				     ew*1., 0.,
				     clusterHits[clusterHits.size()-1]->PeakTime(), clusterHits[clusterHits.size()-1]->SigmaPeakTime(),
				     -999., 0., 
				     -999., 0.,
				    ClusterNo);
std::cout<<"Produced Cluster #"<<ClusterNo<<std::endl;
	      ccol->push_back(cluster);
	      //std::cout<<"no of hits for this cluster is "<<clusterHits.size()<<std::endl;
	      clusterHits.clear();
	      //////
	    }
    
    
    
    
    
   } //clusters
    
     allhits.clear();
     maxBin.clear();
     maxBinValues.clear();
     SortedMaxBin.clear();
     MaxStartPoint.clear();
     MaxEndPoint.clear();
     MaxStartPointTheta.clear();
     MaxEndPointTheta.clear();
     HitsWithClusterID.clear();
     FinalPeaks.clear();
     OriginalmaxBinValues.clear();
     std::cout<<"Should be starting to work on the other plane now"<<std::endl;
    }//Planes
   
   evt.put(ccol);
   
  for(int bin=0; bin< fh_theta_ind_2D->GetNbinsX(); bin++){
   
   fh_theta_ind_2D->SetBinContent(bin,0);
   fh_theta_coll_2D->SetBinContent(bin,0);
   fh_theta_ind->SetBinContent(bin,0);
   fh_theta_coll->SetBinContent(bin,0);
   fh_theta_coll_Area->SetBinContent(bin,0);
   fh_theta_ind_Area->SetBinContent(bin,0);
   }
   
  return;
    }
    
//..............................................................  

void cluster::KingaCluster::AngularDistribution(int plane){   
 
 art::ServiceHandle<geo::Geometry> geom;
 if(fMC==1){
unsigned int channel2,plane2,wire2;  
plane2=plane;

  if(plane==0){
	MCvertex[0]=.3;//force time coordinate to be closer to induction plane 
	}
      else{
	MCvertex[0]=-.3;//force time coordinate to be closer to collection plane
     }
      channel2 = geom->NearestChannel(MCvertex);
      geom->ChannelToWire(channel2,plane2,wire2);   
   std::cout<<"%%%%%%%%%%%%%%%%%%   WIRE VERTEX IS: "<<wire2<<std::endl;
   fwire_vertex.push_back(wire2);
   
  } 
   
//std::vector<unsigned int> fwire_vertex,ftime_vertex;
double a_polar, b_polar,theta_polar;
 ftimetick      =  0.198; //get from parameterset
 fdriftvelocity =  0.157;  //get from paramtereset 9either k and V)
 fpi=3.141592653;

//Define vertex if it is not given by the first hit:

//622_2738:
 // fwire_vertex.push_back(62);
//  fwire_vertex.push_back(93);
//  ftime_vertex.push_back(1058);
//  ftime_vertex.push_back(1060);
//......................................

// fwire_vertex.push_back(62);
// fwire_vertex.push_back(72);
// ftime_vertex.push_back(460);
// ftime_vertex.push_back(460);

//pizero:
// fwire_vertex.push_back(83);
// fwire_vertex.push_back(87);
// ftime_vertex.push_back(420);
// ftime_vertex.push_back(430);
//.........................
//634_2660:
// fwire_vertex.push_back(2);
// fwire_vertex.push_back(2);
// ftime_vertex.push_back(1200);
// ftime_vertex.push_back(1200);

// 650_1728:
// fwire_vertex.push_back(137);
// fwire_vertex.push_back(180);
// ftime_vertex.push_back(710);
// ftime_vertex.push_back(710);

//650_2319:
// fwire_vertex.push_back(66);
// fwire_vertex.push_back(28);
// ftime_vertex.push_back(1225);
// ftime_vertex.push_back(1225);
//..............................
//628_6207:
// fwire_vertex.push_back(48);
// fwire_vertex.push_back(52);
// ftime_vertex.push_back(1350);
// ftime_vertex.push_back(1350);

//628_7100:
// fwire_vertex.push_back(92);
// fwire_vertex.push_back(110);
// ftime_vertex.push_back(850);
// ftime_vertex.push_back(850);




unsigned int channel=0, w=0;
unsigned int p=plane;

//std::cout<<"for PLANE "<<plane<<" fwire_vertex= "<<fwire_vertex[plane]<<" ftime_vertex= "<<ftime_vertex[plane]<<std::endl;


//std::cout<<"No of HITS for plane "<<plane<<" is: "<<allhits.size()<<std::endl;
 
  for(unsigned int i = 0; i< allhits.size(); ++i){
  
  // if(i==0){
//         fwire_vertex=allhits[i]->Wire()->RawDigit()->Channel();
//         ftime_vertex=allhits[i]->PeakTime();
//         std::cout<<"for PLANE "<<plane<<" fwire_vertex= "<<fwire_vertex<<" ftime_vertex= "<<ftime_vertex<<std::endl;
//           }   
 
 channel=allhits[i]->Wire()->RawDigit()->Channel();
 geom->ChannelToWire(channel,p,w);
 
  
      b_polar = (w - fwire_vertex[plane])* 0.4; /**in cm*/
      a_polar = (allhits[i]->PeakTime() - ftime_vertex[plane])* ftimetick *fdriftvelocity; /** in cm*/
      theta_polar = asin(a_polar/sqrt(pow(a_polar,2)+pow(b_polar,2))); /**in rad*/
      theta_polar = 180*theta_polar/fpi; /** in deg*/
     //fh_theta[plane]->Fill(theta_polar,allhits[i]->Charge());
     
if (plane==0 ) {

//std::cout<<"plane= "<<plane<<" theta= "<<theta_polar<<"  channel= "<<allhits[i]->Wire()->RawDigit()->Channel()<<" time= "<<allhits[i]->PeakTime()<<" fwire_vertex[plane]= "<<fwire_vertex[plane]<<" ftime_vertex[plane]= "<<ftime_vertex[plane]<<std::endl;


fh_theta_ind->Fill(theta_polar);
fh_theta_ind_2D->Fill(theta_polar,allhits[i]->Charge());
 fh_theta_ind_Area->Fill(theta_polar,(allhits[i]->EndTime()-allhits[i]->StartTime())* ftimetick *fdriftvelocity*0.4);
}
if (plane==1 ) {

//std::cout<<"plane= "<<plane<<"  w= "<<w<<" time= "<<allhits[i]->PeakTime()<<" theta= "<<theta_polar<<" fwire_vertex[plane]= "<<fwire_vertex[plane]<<" ftime_vertex[plane]= "<<ftime_vertex[plane];


// if(w>95 && w<118 && allhits[i]->PeakTime()>500 && allhits[i]->PeakTime()<800){
// std::cout<<"  ***********"<<std::endl;
// }
// else{ std::cout<<std::endl;}




fh_theta_coll->Fill(theta_polar);
fh_theta_coll_2D->Fill(theta_polar,allhits[i]->Charge());
fh_theta_coll_Area->Fill(theta_polar,(allhits[i]->EndTime()-allhits[i]->StartTime())* ftimetick *fdriftvelocity*0.4);

}
  }
 
  // for(int i=0;i<startTimes.size();i++){
//   std::cout<<"startTime is at bin = "<<startTimes[i]<<" and its value is "<<fh_theta_ind->GetBinContent(startTimes[i])<<std::endl;
//   }
//   for(int i=0;i<endTimes.size();i++){
//   std::cout<<"endTime is at bin = "<<endTimes[i]<<" and its value is "<<fh_theta_ind->GetBinContent(endTimes[i])<<std::endl;
//   }
 

    }
    
    
    
//..............................................................   
void cluster::KingaCluster::FindMax(int plane){  

 std::cout<<"No of bins= "<<fh_theta_ind->GetNbinsX()<<std::endl;
  std::cout<<" Bincontent= "<<fh_theta_ind->GetBinContent(48)<<std::endl;
  std::cout<<" Bincontent= "<<fh_theta_ind->GetBinContent(49)<<std::endl;
  std::cout<<" Bincontent= "<<fh_theta_ind->GetBinContent(50)<<std::endl;
// std::vector<int> PossibleFinalMax, FinalMax;
   std::vector<int> startTimes;  //stores time of 1st local minimum
    // std::vector<int> maxBin;    //stores time of local maximum
    std::vector<int> endTimes;    //stores time of 2nd local minimum
   bool maxFound; //Flag for whether a value>threshold has been found
  
   int minTimeHolder;
 
  startTimes.clear();
  maxBin.clear();
  endTimes.clear();
 std::cout<<"We have "<<allhits.size()<<" hits for plane "<<plane<<std::endl;

 //  double threshold=76*allhits.size()/2;
//   double MinThreshold=1000;
//  double threshold=6;
//   double MinThreshold=4;

// worked for pizero and nu_e :
 // double threshold=80000;
//   double MinThreshold=60000;
  //............................
  
  
  
  // double threshold=30000;
//   double MinThreshold=20000;

  //for lines:
  // double threshold=20000;
//   double MinThreshold=10000;

  double threshold=10000;
  double MinThreshold=5000;
  
  
  int time=1;
  int ValidPeak=0;
  std::cout<<"Threshold that a peak must have in order to be considered a peak = "<<threshold<<". For inflection points we must have the value to drop to "<<MinThreshold<<std::endl;
  
 //......................................................... 
 //collection plane:
 if (plane==1){
 
  for(int bin=1; bin<fh_theta_coll_2D->GetNbinsX()+1;bin++){
 
  if(fh_theta_coll_2D->GetBinContent(bin)>fh_theta_coll_2D->GetBinContent(bin+1) && fh_theta_coll_2D->GetBinContent(bin+1)<fh_theta_coll_2D->GetBinContent(bin+2)) {
      //only add points if we've already found a local max above threshold.
      if(maxFound) {
        endTimes.push_back(time+1);
        maxFound = false;
        //keep these in case new hit starts right away
        minTimeHolder = time+2;
       }
      else {minTimeHolder = time+1; }
    }
  //if not a minimum,-> test if we are at a local maximum
    //if so and the max value is above threshold add it and proceed.
    else if(fh_theta_coll_2D->GetBinContent(bin)<fh_theta_coll_2D->GetBinContent(bin+1) &&
        fh_theta_coll_2D->GetBinContent(bin+1)>fh_theta_coll_2D->GetBinContent(bin+2) &&
        fh_theta_coll_2D->GetBinContent(bin+1) > threshold) {
      maxFound = true;
      maxBin.push_back(time+1);
      startTimes.push_back(minTimeHolder);         
    }

    time++; 
     
  }
  
  if(maxBin.size()==0){
  
  std::cout<<"COLLECTION PLANE: "<<std::endl;
  std::cout<<" COULDN'T FIND ANY MAXIMA IN YOUR THETA DISTRIBUTION!!!! PROGRAM WILL NOT PROCEED!!!"<<std::endl;
    fpeaks_found=0;}
 if(maxBin.size()>0){     
     
 for(unsigned int i=0;i<maxBin.size();i++){
  std::cout<<"maxTime is at bin = "<<maxBin[i]<<" and its value is "<<fh_theta_coll_2D->GetBinContent(maxBin[i])<<std::endl;
  }
  
  
  // Lets make sure that the first bin in the maxBin corresponds to the highest peak:

//std::vector<double> maxBinValues;
  for(unsigned int i=0;i<maxBin.size();i++){
  maxBinValues.push_back(fh_theta_coll_2D->GetBinContent(maxBin[i]));
  OriginalmaxBinValues.push_back(fh_theta_coll_2D->GetBinContent(maxBin[i]));
  }
  std::cout<<"The largest is at position:  "<<std::distance(maxBinValues.begin(),std::max_element(maxBinValues.begin(),maxBinValues.end()))<<" which corresponds to bin #   "<<maxBin[std::distance(maxBinValues.begin(),std::max_element(maxBinValues.begin(),maxBinValues.end()))]<<" and its value= "<<*std::max_element(maxBinValues.begin(),maxBinValues.end())<<std::endl;
  
  //sort values from the largest to the smallest in maxBinValues, then find the corresponding bin numbers to create SortedMaxBin:
  
  sort(maxBinValues.begin(),maxBinValues.end());
 std::cout<<"maxBinValues after sort:"<<std::endl;
 for(unsigned int i=0;i<maxBinValues.size();i++){

   std::cout<<maxBinValues[i]<<std::endl;
 }

reverse (maxBinValues.begin(),maxBinValues.end());
 std::cout<<"maxBinValues in the correct order are now:"<<std::endl;
 for(unsigned int i=0;i<maxBinValues.size();i++){

   std::cout<<maxBinValues[i]<<std::endl;
 }

 for(unsigned int i=0; i<maxBinValues.size();i++){

   std::vector<double>::iterator pos=std::find( OriginalmaxBinValues.begin(), OriginalmaxBinValues.end(),maxBinValues[i]);
   SortedMaxBin.push_back(maxBin[pos-OriginalmaxBinValues.begin()]);

 }


  //create SortedMaxBin vector whose first element is the global max:
  
 //  SortedMaxBin.push_back(maxBin[std::distance(maxBinValues.begin(),std::max_element(maxBinValues.begin(),maxBinValues.end()))]);
//   for(int i=0; i<maxBin.size();i++){
//    if(i!=std::distance(maxBinValues.begin(),std::max_element(maxBinValues.begin(),maxBinValues.end())))
//   SortedMaxBin.push_back(maxBin[i]);
  
//   }
  
  std::cout<<"SortexMaxBin elements are: "<<std::endl;
for(unsigned int i=0; i<SortedMaxBin.size(); i++)
{

std::cout<<SortedMaxBin[i]<<std::endl;

}
// int ValidPeak=0;
  // loop over maxima and find where they start on the left and right side, form cluster for each:
  
  for(unsigned int maxNo=0; maxNo<SortedMaxBin.size();maxNo++){
  
  //loop over the ranges and make sure that your peaks don't fall into already formed clusters
  //std::cout<<"Right now we have "<<MaxStartPoint.size()<<" ranges"<<std::endl;
  if(MaxStartPoint.size()==0){ValidPeak=1;}
  for(unsigned int NoRange=0; NoRange<MaxStartPoint.size();NoRange++){
  //std::cout<<"Checking peak "<<SortedMaxBin[maxNo]<<std::endl;
  if(SortedMaxBin[maxNo]>MaxStartPoint[NoRange] && SortedMaxBin[maxNo]<MaxEndPoint[NoRange]){
  //maxNo++;
  //std::cout<<"this peak is out of the picture! --> "<<SortedMaxBin[maxNo]<<std::endl;
  ValidPeak=0;
  break;}
  else {ValidPeak=1;
  //std::cout<<"passed"<<std::endl;
  }
  }
   if(ValidPeak==1){
   
   std::cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<std::endl;
  FinalPeaks.push_back(SortedMaxBin[maxNo]);
   std::cout<<"We are working on peak at bin #"<<SortedMaxBin[maxNo]<<std::endl;
  //start at the peak and go left
      for(int LeftBin=SortedMaxBin[maxNo]-1;LeftBin>SortedMaxBin[maxNo]-30; LeftBin--)
      {
        if(fh_theta_coll_2D->GetBinContent(LeftBin)<MinThreshold && fh_theta_coll_2D->GetBinContent(LeftBin-1)>fh_theta_coll_2D->GetBinContent(LeftBin)) {
           MaxStartPoint.push_back(LeftBin);
          // std::cout<<"picked option 1, startin point @bin "<<LeftBin<<std::endl;
           std::cout<<"For peak at bin= "<<SortedMaxBin[maxNo]<<"("<<-180+2*SortedMaxBin[maxNo]<<" degrees)"<<" LeftBin="<<LeftBin<<"("<<-180+2*LeftBin<<" degrees)"<<" RightBin= ";
           break;
           }
           else if(fh_theta_coll_2D->GetBinContent(LeftBin)<MinThreshold && fh_theta_coll_2D->GetBinContent(LeftBin-1)==0){
           MaxStartPoint.push_back(LeftBin-1);
           //std::cout<<"picked option 2, startin point @bin "<<LeftBin-1<<std::endl;
           std::cout<<"For peak at bin= "<<SortedMaxBin[maxNo]<<"("<<-180+2*SortedMaxBin[maxNo]<<" degrees)"<<" LeftBin="<<LeftBin-1<<"("<<-180+2*(LeftBin-1)<<" degrees)"<<" RightBin= ";
           break;
           
           }
           else if (LeftBin==SortedMaxBin[maxNo]-29){std::cout<<" cannot find starting point of the peak!!!!!"<<std::endl;}
         
      }
     
   
   
   
  
  for(int RightBin=SortedMaxBin[maxNo]+1;RightBin<SortedMaxBin[maxNo]+30; RightBin++)
      {
        if(fh_theta_coll_2D->GetBinContent(RightBin)<MinThreshold && fh_theta_coll_2D->GetBinContent(RightBin+1)>fh_theta_coll_2D->GetBinContent(RightBin)) {
           MaxEndPoint.push_back(RightBin);
           std::cout<<RightBin<<"("<<-180+2*RightBin<<" degrees)"<<std::endl;
           break;
           }
         else if(fh_theta_coll_2D->GetBinContent(RightBin)<MinThreshold && fh_theta_coll_2D->GetBinContent(RightBin+1)==0){
           MaxEndPoint.push_back(RightBin+1);
           std::cout<<RightBin+1<<"("<<-180+2*(RightBin+1)<<" degrees)"<<std::endl;
           break;
           }
         else if(RightBin==SortedMaxBin[maxNo]+29){std::cout<<" cannot find end point of the peak!!!!!"<<std::endl;}
         
      }
      
      
      
      
  } //valid peak
  
  ValidPeak=0;
  
}//peaks
  
  

  } //if maxBin.size()>0
  
  
  
 } //plane 1
 startTimes.clear();
  //maxBin.clear();
  endTimes.clear();
  time=1;
  ValidPeak=0;
  
 //......................................................... 
 
 
 //induction plane
 if(plane==0){
 
 std::cout<<"No of bins= "<<fh_theta_ind_2D->GetNbinsX()<<std::endl;
for(int bin=1; bin<fh_theta_ind_2D->GetNbinsX()+1;bin++){
 
  if(fh_theta_ind_2D->GetBinContent(bin)>fh_theta_ind_2D->GetBinContent(bin+1) && fh_theta_ind_2D->GetBinContent(bin+1)<fh_theta_ind_2D->GetBinContent(bin+2)) {
      //only add points if we've already found a local max above threshold.
      if(maxFound) {
        endTimes.push_back(time+1);
        maxFound = false;
        //keep these in case new hit starts right away
        minTimeHolder = time+2;
       }
      else {minTimeHolder = time+1; }
    }
  //if not a minimum test if we are at a local maximum
    //if so and the max value is above threshold add it and proceed.
    else if(fh_theta_ind_2D->GetBinContent(bin)<fh_theta_ind_2D->GetBinContent(bin+1) &&
        fh_theta_ind_2D->GetBinContent(bin+1)>fh_theta_ind_2D->GetBinContent(bin+2) &&
        fh_theta_ind_2D->GetBinContent(bin+1) > threshold) {
      maxFound = true;
      maxBin.push_back(time+1);
      startTimes.push_back(minTimeHolder);         
    }
    time++;
 
  }
 std::cout<<"INDUCTION PLANE: "<<std::endl;

  if(maxBin.size()==0){
  std::cout<<" COULDN'T FIND ANY MAXIMA IN YOUR THETA DISTRIBUTION!!!! PROGRAM WILL NOT PROCEED!!!"<<std::endl;
    fpeaks_found=0;}
    
 if(maxBin.size()>0){   
 std::cout<<"maxBin.size()= "<<maxBin.size()<<std::endl;
  for(unsigned int i=0;i<maxBin.size();i++){
  std::cout<<"maxTime is at bin = "<<maxBin[i]<<" ("<<-180+2*maxBin[i]<<" degrees)"<<" and its value is "<<fh_theta_ind_2D->GetBinContent(maxBin[i])<<std::endl;
std::cout<<"...................................."<<std::endl;
 }//maxBin
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////





// Lets make sure that the first bin in the maxBin corresponds to the highest peak:

//std::vector<double> maxBinValues;
  for(unsigned int i=0;i<maxBin.size();i++){
  maxBinValues.push_back(fh_theta_ind_2D->GetBinContent(maxBin[i]));
  OriginalmaxBinValues.push_back(fh_theta_ind_2D->GetBinContent(maxBin[i]));

  }
  std::cout<<"The largest is at position:  "<<std::distance(maxBinValues.begin(),std::max_element(maxBinValues.begin(),maxBinValues.end()))<<" which corresponds to bin #   "<<maxBin[std::distance(maxBinValues.begin(),std::max_element(maxBinValues.begin(),maxBinValues.end()))]<<" and its value= "<<*std::max_element(maxBinValues.begin(),maxBinValues.end())<<std::endl;
  
  //sort values from the largest to the smallest in maxBinValues, then find the corresponding bin numbers to create SortedMaxBin:
  
  sort(maxBinValues.begin(),maxBinValues.end());
 std::cout<<"maxBinValues after sort:"<<std::endl;
 for(unsigned int i=0;i<maxBinValues.size();i++){

   std::cout<<maxBinValues[i]<<std::endl;
 }

reverse (maxBinValues.begin(),maxBinValues.end());
 std::cout<<"maxBinValues in the correct order are now:"<<std::endl;
 for(unsigned int i=0;i<maxBinValues.size();i++){

   std::cout<<maxBinValues[i]<<std::endl;
 }

 for(unsigned int i=0; i<maxBinValues.size();i++){

   std::vector<double>::iterator pos=std::find( OriginalmaxBinValues.begin(), OriginalmaxBinValues.end(),maxBinValues[i]);
   SortedMaxBin.push_back(maxBin[pos-OriginalmaxBinValues.begin()]);

 }
  
  
  
  
  
  
  
  
  //create SortedMaxBin vector whose first element is the global max:
  
  
  // std::vector<int> SortedMaxBin;
  // SortedMaxBin.push_back(maxBin[std::distance(maxBinValues.begin(),std::max_element(maxBinValues.begin(),maxBinValues.end()))]);
//   for(int i=0; i<maxBin.size();i++){
//    if(i!=std::distance(maxBinValues.begin(),std::max_element(maxBinValues.begin(),maxBinValues.end())))
//   SortedMaxBin.push_back(maxBin[i]);
//   
//   }
  
  std::cout<<"SortexMaxBin elements are: "<<std::endl;
for(unsigned int i=0; i<SortedMaxBin.size(); i++)
{

std::cout<<SortedMaxBin[i]<<std::endl;

}

  // loop over maxima and find where they start on the left and right side, form cluster for each:
  
  for(unsigned int maxNo=0; maxNo<SortedMaxBin.size();maxNo++){
  
  //loop over the ranges and make sure that your peaks don't fall into already formed clusters
  if(MaxStartPoint.size()==0){ValidPeak=1;}
  for(unsigned int NoRange=0; NoRange<MaxStartPoint.size();NoRange++){
  if(SortedMaxBin[maxNo]>MaxStartPoint[NoRange] && SortedMaxBin[maxNo]<MaxEndPoint[NoRange]){
  ValidPeak=0;
  break;}
  else {ValidPeak=1;}
  }
  
   if(ValidPeak==1){
   
  std::cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<std::endl;
  FinalPeaks.push_back(SortedMaxBin[maxNo]);
  //start at the peak and go left
      for(int LeftBin=SortedMaxBin[maxNo]-1;LeftBin>SortedMaxBin[maxNo]-30; LeftBin--)
      {
        if(fh_theta_ind_2D->GetBinContent(LeftBin)<MinThreshold && fh_theta_ind_2D->GetBinContent(LeftBin-1)>fh_theta_ind_2D->GetBinContent(LeftBin)) {
           MaxStartPoint.push_back(LeftBin);
            
           std::cout<<"For peak at bin= "<<SortedMaxBin[maxNo]<<"("<<-180+2*SortedMaxBin[maxNo]<<" degrees)"<<" LeftBin="<<LeftBin<<"("<<-180+2*LeftBin<<" degrees)"<<" RightBin= ";
        
           break;
           }
           else if(fh_theta_ind_2D->GetBinContent(LeftBin)<MinThreshold && fh_theta_ind_2D->GetBinContent(LeftBin-1)==0){
           MaxStartPoint.push_back(LeftBin-1);
           //std::cout<<"picked option 2, startin point @bin "<<LeftBin-1<<std::endl;
           std::cout<<"For peak at bin= "<<SortedMaxBin[maxNo]<<"("<<-180+2*SortedMaxBin[maxNo]<<" degrees)"<<" LeftBin="<<LeftBin-1<<"("<<-180+2*(LeftBin-1)<<" degrees)"<<" RightBin= ";
           break;
           
           }
           else if (LeftBin==SortedMaxBin[maxNo]-30){std::cout<<"cannot find starting point of the peak!!!!!"<<std::endl;}
         
      }
  
  
  for(int RightBin=SortedMaxBin[maxNo]+1;RightBin<SortedMaxBin[maxNo]+30; RightBin++)
      {
        if(fh_theta_ind_2D->GetBinContent(RightBin)<MinThreshold && fh_theta_ind_2D->GetBinContent(RightBin+1)>fh_theta_ind_2D->GetBinContent(RightBin)) {
           MaxEndPoint.push_back(RightBin);
           std::cout<<RightBin<<"("<<-180+2*RightBin<<" degrees)"<<std::endl;
           break;
           }
           else if(fh_theta_ind_2D->GetBinContent(RightBin)<MinThreshold && fh_theta_ind_2D->GetBinContent(RightBin+1)==0){
           MaxEndPoint.push_back(RightBin+1);
           std::cout<<RightBin+1<<"("<<-180+2*(RightBin+1)<<" degrees)"<<std::endl;
           break;
           }
           else if(RightBin==SortedMaxBin[maxNo]+20){std::cout<<"cannot find end point of the peak!!!!!"<<std::endl;}
         
      }
  
  }
  ValidPeak=0;
  
}//peaks
  
  
}// if maxBin.size()>0 this means that we can find peaks in the theta distribution 
  
  
 } //plane 0
 //......................................................... 
}   

//..............................................................   
// void cluster::KingaCluster::FinalPeaks(){  
// std::cout<<"In FinalPeaks()"<<std::endl;
// // for(int i=0;i<maxBin.size();i++){
// //   std::cout<<"maxTime is at bin = "<<maxBin[i]<<" and its value is "<<fh_theta_ind_2D->GetBinContent(maxBin[i])<<std::endl;
// //   }
// 
// 
// 
// }
void cluster::KingaCluster::FitAngularDistributions(){  

//Now Fit gaussian:
 
  // TF1 *gau = new TF1("gaus","gaus",-6, 6);
//
//     fh_theta_ind->Fit("gaus","QR");        /** Fit of the angular distribution*/
//     std::cout<<"Mean of the fit = "<< gau->GetParameter(1)<<std::endl;/** Mean value of the fit */
//     std::cout<<"RMS= "<< gau->GetParameter(2)<<std::endl; /** RMS of the fit of the angular distribution in deg*/

}

//..............................................................   


void cluster::KingaCluster::FindClusters(int plane){ 

std::cout<<"FORMING CLUSTERS NOW :) "<<std::endl;
std::cout<<"FinalPeaks are at bin(s):  ";
for(unsigned int i=0; i<FinalPeaks.size();i++)
{
std::cout<<FinalPeaks[i]<<" which corresponds to angle=  "<<-180+2*FinalPeaks[i]<<std::endl;

}


art::ServiceHandle<geo::Geometry> geom;
unsigned int channel=0, w=0;
unsigned int p=plane;


//int no_noise_hits=0;
std::cout<<"In FindClusters(int plane), we should be producing "<<MaxStartPoint.size()<<" clusters"<<std::endl;
double a_polar, b_polar,theta_polar;

//HitsWithClusterID.clear();

std::vector<double> DiffAngles;
DiffAngles.clear();



for(unsigned int i = 0; i< allhits.size(); ++i){

 channel=allhits[i]->Wire()->RawDigit()->Channel();
 geom->ChannelToWire(channel,p,w);
 
b_polar = (w - fwire_vertex[plane])* 0.4; /**in cm*/
      a_polar = (allhits[i]->PeakTime() - ftime_vertex[plane])* ftimetick *fdriftvelocity; /** in cm*/
      theta_polar = asin(a_polar/sqrt(pow(a_polar,2)+pow(b_polar,2))); /**in rad*/
      theta_polar = 180*theta_polar/fpi; /** in deg*/
for(unsigned int ClusterNo=0; ClusterNo<MaxStartPoint.size();ClusterNo++){

  if(theta_polar>=(-180+2*MaxStartPoint[ClusterNo]) && theta_polar<=(-180+2*MaxEndPoint[ClusterNo])){
  //want to start counting from 1, O is reserved for hits that will be marked as noise
   HitsWithClusterID.push_back(ClusterNo+1);
   break;}
   else if(ClusterNo==MaxStartPoint.size()-1){
   //decide where noise hits go
   
  //  if(w>95 && w<128 && allhits[i]->PeakTime()> 880 && allhits[i]->PeakTime()<1048){
//    std::cout<<"Noise hit at w= "<<w<<" t= "<<allhits[i]->PeakTime()<<" with theta_polar= "<<theta_polar;
//   // std::cout<<"FinalPeaks.size()= "<<FinalPeaks.size()<<std::endl;
//   }
   for(unsigned int peakNo=0;peakNo<FinalPeaks.size();peakNo++){
   DiffAngles.push_back(fabs(-180+2*FinalPeaks[peakNo]-theta_polar));
   // if(w>95 && w<128 && allhits[i]->PeakTime()> 880 && allhits[i]->PeakTime()<1048){
//    std::cout<<"diff for peak "<<peakNo<<" is "<<fabs(-180+2*FinalPeaks[peakNo]-theta_polar)<<std::endl;
//    }
   }
   //now take minimum of DiffAngles and find at which position it is at, this position corresponds to clusterNo +1 , because we don't want to mark hits with zero 
   
   int position=std::distance(DiffAngles.begin(),std::min_element(DiffAngles.begin(),DiffAngles.end()));
   
   HitsWithClusterID.push_back(position+1);
   // if(w>95 && w<128 && allhits[i]->PeakTime()> 880 && allhits[i]->PeakTime()<1048){
//    std::cout<<"  This hit is closest to cluster # "<<position+1<<std::endl;
//    }
   //no_noise_hits++;
   DiffAngles.clear();
   }


} //loop over all ranges


} //allhits


//std::cout<<"In FindClusters(int plane), marked "<<no_noise_hits<<" hits as NOISE"<<std::endl;
//std::cout<<"HitsWithClusterID contains the following clusterIDs:"<<std::endl;

//for(int i=0; i<HitsWithClusterID.size();i++){

//std::cout<<HitsWithClusterID[i]<<"  ";


//}

//std::cout<<std::endl;






}
