////////////////////////////////////////////////////////////////////////
//
// DBSCANfinder.cxx
//
// kinga.partyka@yale.edu
//
//  This algorithm finds clusters of hits, they can be of arbitrary shape.You need to specify 2(3) parameters: epsilon, epsilon2 and MinPoints as explained in the corresponding xml file.In my comments a 'point' reference appears quite often. A 'point' is basically a simple hit which only contains wire and time information. This algorithm is based on DBSCAN(Density Based Spatial Clustering of Applications with Noise): M. Ester, H.-P. Kriegel, J. Sander, and X. Xu, A density-based algorithm for discovering clusters in large spatial databases with noise, Second International Conference on Knowledge Discovery and Data Mining, pp. 226-231, AAAI Press. 1996. 
// ( Some of this code is from "Antonio Gulli's coding playground")  
// (Once you run it, you can then display the discovered clusters on event display, each cluser is going to appear in a different color.)
////////////////////////////////////////////////////////////////////////


//Framework includes:
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

#include "ClusterFinder/DBcluster.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "Geometry/geo.h"
#include "SimulationBase/simbase.h"
#include "Simulation/sim.h"
#include "RecoBase/recobase.h"
#include "TGeoManager.h"
#include "TH1.h"

//-------------------------------------------------
cluster::DBcluster::DBcluster(edm::ParameterSet const& pset) : 
  
  fhitsModuleLabel(pset.getParameter< std::string >("HitsModuleLabel")),
  fEps            (pset.getParameter< double      >("eps")            ),
  fEps2           (pset.getParameter< double      >("eps2")           ),
  fMinPts         (pset.getParameter< int         >("minPts")         )
  
  // std::cerr << "Set epsilon distance = " << fEps << std::endl;
  // std::cerr << "Set epsilon distance = " << fEps2 << std::endl;
  // std::cerr << "Set Minimum number of Points that have to be in a neighborhood of a point = " << fMinPts << std::endl;
  
{  
  produces<std::vector<recob::Cluster> >();
  
  
}

//-------------------------------------------------
cluster::DBcluster::~DBcluster()
{
}

//-------------------------------------------------
void cluster::DBcluster::beginJob(edm::EventSetup const&){
  // get access to the TFile service
  edm::Service<edm::TFileService> tfs;

  fhitwidth= tfs->make<TH1F>(" fhitwidth","width of hits in cm", 50000,0 ,5  );
  fhitwidth_ind_test= tfs->make<TH1F>("fhitwidth_ind_test","width of hits in cm", 50000,0 ,5  );
  fhitwidth_coll_test= tfs->make<TH1F>("fhitwidth_coll_test","width of hits in cm", 50000,0 ,5  );
    
}

double cluster::DBScan::getSimilarity(const std::vector<double> v1, const std::vector<double> v2){
  
   
  //for Euclidean distance comment everything out except this-->>>
  // return sqrt((v2[1]-v1[1])*(v2[1]-v1[1])+(v2[0]-v1[0])*(v2[0]-v1[0]));
  //------------------------------------------------------------------------
  // return fabs( v2[0]-v1[0]); //for rectangle
  //---------------------------------------------------------------------- 
  //Manhattan distance:
  //return fabs(v1[0]-v2[0])+fabs(v1[1]-v2[1]);
  //---------------------------------------------------------------------
  return (( v2[0]-v1[0])*( v2[0]-v1[0])); //for ellipse

}

//----------------------------------------------------------------
double cluster::DBScan::getSimilarity2(const std::vector<double> v1, const std::vector<double> v2){

  //-------------------------------------------
  //return fabs( v2[1]-v1[1]);//for rectangle
  //------------------------------------------
  return (( v2[1]-v1[1])*( v2[1]-v1[1]));//for ellipse
  
  
}

//----------------------------------------------------------------
double cluster::DBScan::getWidthFactor(const std::vector<double> v1, const std::vector<double> v2){
 
  //double k=0.13; //this number was determined by looking at flat muon hits' widths. The average width of these hits in cm is 0.505, so 4*2*(w1^2)=2.04 where w1=w2=0.505, e^2.044= 7.69. In order not to change the distance in time direction of the ellipse we want to make it equal to 1 for these hits. Thus the k factor is k=1/7.69=0.13//for coeff=4

  //std::cout<<"v1[2]= "<<v1[2]<<" v2[2]= "<<v2[2]<<" so the factor is : "<<(( v1[2]*v1[2])+( v2[2]*v2[2]))*k<<std::endl;

  //double k=0.78;
  //..................................................
  double k=0.1;//for 4.5 coeff
  double WFactor=(exp(4.6*(( v1[2]*v1[2])+( v2[2]*v2[2]))))*k;
  //........................................................
  //Let's try something different:
  // double k=1.96;
  //    double WFactor=(( v1[2]*v1[2])+( v2[2]*v2[2]))*k;
  if(WFactor>1){
    if(WFactor<6.25) {return WFactor;}//remember that we are increasing the distance in eps2 as sqrt of this number (i.e sqrt(6.25))
    else {return 6.25;}
   
  }
  else {return 1.0;}
  //.............................................
  
  
  
}

//----------------------------------------------------------------

	
std::vector<unsigned int> cluster::DBScan::findNeighbors( unsigned int pid, double threshold,double threshold2)
{
  // std::cout<<"EPSILON IS: "<<threshold<<std::endl;
  //std::cout<<"EPSILON2 IS: "<<threshold2<<std::endl;
  std::vector<unsigned int> ne;
  
  for ( int unsigned j=0; j < _sim.size(); j++)
    {
      
      //  if 	((pid != j ) && ((_sim[pid][j]) < threshold) )//for circle
      //----------------------------------------------------------------------------------------------    
      // if 	((pid != j ) && ((_sim[pid][j]) < threshold) && ((_sim2[pid][j]) < threshold2))//for rectangle
      //----------------------------------------------------------------------------------------------
      if 	((pid != j ) && (((_sim[pid][j])/ (threshold*threshold))+ ((_sim2[pid][j])/ (threshold2*threshold2*(_sim3[pid][j]))))<1) //ellipse
	//*(_sim3[pid][j])))) < 1 ) //for ellipse
	//-------------------------------------------------------------------------------------

	{
	
	  ne.push_back(j);
	  // std::cout<<"Neighbor of pid "<<pid<<" is"<<j<<" their WFactor= "<<_sim3[pid][j]<<std::endl;
	  for(int unsigned i=0; i<ne.size();i++){
	    // std::cout<<ne[i]<<" ";
	  }
	  // std::cout<<std::endl;
	  //std::cout << "sim(" << pid  << "," << j << ")" << _sim(pid, j) << ">" << threshold << std::endl;
	}

      //$$$
      else{
 	if(pid != j ){
	  // std::cout<<"point "<<pid<<" and "<<j<<" NOT neighbors b/c not in ellipse (WFactor= "<<_sim3[pid][j]<<" ), not less than 1 but "<<((_sim[pid][j])/ (threshold*threshold))+ ((_sim2[pid][j])/ (threshold2*threshold2*(_sim3[pid][j])))<<std::endl;
	    
	}
      }
      //$$$
    }

  //$$$
  //  std::cout<<"Neighbors of pid "<<pid<<" : ";
  // 	  for(int unsigned i=0; i<ne.size();i++){
  // 	     std::cout<<ne[i]<<" ";
  // 	  }
  // 	  std::cout<<std::endl;

  //$$$




  return ne;
};
//-----------------------------------------------------------------

void cluster::DBScan::computeSimilarity()
{
  int size = _ps.size();
  _sim.resize(size, std::vector<double>(size));
  //	std::cout<<"running computeSimilarity"<<std::endl;
  for ( int i=0; i < size; i++)
    {
      for ( int j=i+1; j < size; j++)
	{
	  _sim[j] [i] = _sim[i][ j] = getSimilarity(_ps[i], _ps[j]);
		
		
	  //	std::cout << "similarity for (" << i << ", " << j << ")=" << _sim[i] [j] << " ";
	}
      // std::cout << std::endl;
    }
};



//------------------------------------------------------------------
void cluster::DBScan::computeSimilarity2()
{
  int size = _ps.size();
  _sim2.resize(size, std::vector<double>(size));
  //	std::cout<<"running computeSimilarity"<<std::endl;
  for ( int i=0; i < size; i++)
    {
      for ( int j=i+1; j < size; j++)
	{

	  _sim2[j] [i] = _sim2[i][ j] = getSimilarity2(_ps[i], _ps[j]);
	  //	_sim[j] [i] = _sim[i][ j] = 0.9;
		
	  //	std::cout << "similarity for (" << i << ", " << j << ")=" << _sim[i] [j] << " ";
	}
      // std::cout << std::endl;
    }
};



//------------------------------------------------------------------

void cluster::DBScan::computeWidthFactor()
{
  int size = _ps.size();
  _sim3.resize(size, std::vector<double>(size));
       
  for ( int i=0; i < size; i++)
    {
      for ( int j=i+1; j < size; j++)
	{

	  _sim3[j] [i] = _sim3[i][ j] = getWidthFactor(_ps[i], _ps[j]);
		
	  //std::cout<<"width factor is "<<_sim3[i][j]<<" ";
	}
      // std::cout << std::endl;
    }
};



//------------------------------------------------------------------

//single point output
std::ostream& cluster::operator<<(std::ostream& o,const std::vector<double>& p)
{
  o << "{ ";
  
  for(unsigned int x=0;x<p.size();x++)
    {
      
      o<<" "<<p[x];
      // o<<"SIZE OF POINT IS: "<<p.size()<<" and the point is: "<<p[x];
    }
  o << " }, ";
  
  return o;
}

//--------------------------------------------------------------------

// clusters output
std::ostream& cluster::operator<<(std::ostream& o, const cluster::DBScan& cs)
{
   
  for(unsigned int i=0;i<cs._clusters.size();i++)
    {
      
      o<<"c("<<i+1<<")=";
      
      
      for(unsigned int j=0;j<cs._pointId_to_clusterId.size();j++)
	{
	  if (cs._pointId_to_clusterId[j]==(i+1)){
	      
	    o<<cs._ps[j];
	  }
	}//for
      o << std::endl;
    }//for
  return o;
}
//----------------------------------------------------------------
/////////////////////////////////////////////////////////////////
// This is the algorithm that finds clusters:
void cluster::DBScan::run_cluster() 

{
  
  
  unsigned int cid = 1;
  // foreach pid
  for ( unsigned int pid = 0; pid < _ps.size(); pid++)
    {
      // not already visited
      if (!_visited[pid]){  
	
	_visited[pid] = true;
	//me:
	//	std::cout<<"_visited vector is: ";
	for(unsigned int i=0;i<_visited.size();i++){
	  // std::cout<<_visited[i]<<" ";
	}
	//	std::cout<<std::endl;
	// get the neighbors
	std::vector<unsigned int> ne = findNeighbors(pid, _eps,_eps2);
	
	// not enough support -> mark as noise
	if (ne.size() < _minPts)
	  {
	    // std::cout<<"Not enough neighbors for point: "<<pid<<"(It has: "<<ne.size()<<" ) ! "<<std::endl;
	    _noise[pid] = true;
	    // std::cout<<"_noise vector is: ";
	    //for(unsigned int i=0;i<_noise.size();i++){
	    //  std::cout<<_noise[i]<<" "; }
	    // std::cout<<std::endl;
	  } 
	    
	else 
	  {
	    
	    // std::cout<<" Yes! point "<<pid<<" has enough neighbors ("<<ne.size()<<")"<<std::endl;
	    
	    // Add p to current cluster
	    
	    std::vector<unsigned int> c;              // a new cluster
	   
	    c.push_back(pid);   	// assign pid to cluster
	    //  std::cout<<"Cluster c now contains the following pids: ";
	    //  for(unsigned int i=0;i<c.size();i++){
	    // 	      // std::cout<<c[i]<<" ";
	    // 	    }
	    // std::cout<<std::endl;
	    _pointId_to_clusterId[pid]=cid;
	    // std::cout<<"_pointId_to_clusterId is now: ";
	    //  for(unsigned int i=0;i<_pointId_to_clusterId.size();i++){
	    // 	      std::cout<<_pointId_to_clusterId[i]<<" ";
	    // 	    }
	    //  std::cout<<std::endl;
	    // go to neighbors
	    for (unsigned int i = 0; i < ne.size(); i++)
	      {
		unsigned int nPid = ne[i];
		
		// not already visited
		if (!_visited[nPid])
		  {
		    _visited[nPid] = true;
		    // std::cout<<"_visited vector is now(2): ";
		    for(unsigned int i=0;i<_visited.size();i++){
		      // std::cout<<_visited[i]<<" ";
		    }
		    // std::cout<<std::endl;
		    // go to neighbors
		    std::vector<unsigned int> ne1 = findNeighbors(nPid, _eps, _eps2);
		    //   std::cout<<"found this many: "<<ne1.size()<<" neighbors for pointid "<<nPid<<". THey are: ";
		    //  for(unsigned int i=0;i<ne1.size();i++){
		    // 		      //  std::cout<<ne1[i]<<" ";
		    // 		    }
		    
		    //  std::cout<<std::endl;
		    // enough support
		    if (ne1.size() >= _minPts)
		      {
		       	
			//	std::cout<<"ok, so let's include pid: "<<nPid<<std::endl;
			// join
			
			for(unsigned int i=0;i<ne1.size();i++)
			  {
			    // join neighbord
			    //ne.push_back(n1); 
			    ne.push_back(ne1[i]); 
			    
			  }
			//	std::cout<<" neighbors for "<<nPid<<" are now: ";
			for(unsigned int i=0;i<ne.size();i++){
			  // std::cout<<ne[i]<<" ";
			}
			//	std::cout << std::endl;
		      }
		  }
		
		// not already assigned to a cluster
		if (!_pointId_to_clusterId[nPid])
		  {
		    //  std::cout<<"pointid: "<<nPid<<" hasn't been added to a cluster, so adding now!"<<std::endl;
		    
		    c.push_back(nPid);
		    _pointId_to_clusterId[nPid]=cid;
		  }
	      }
	    
	    //  std::cout<<"Our cluster #"<<cid<<" has the following pointids: ";
	    // 	    for(unsigned int i=0; i<c.size();i++){
	    // 	      std::cout<<c[i]<<" ";
	      
	    // 	    }
	    // 	     std::cout<<std::endl;
	    _clusters.push_back(c);
	    
	    
	    cid++;
	  }
      } // if (!visited
    } // for
  
  //  std::cout<<"final _noise vector is: ";
  //    for(unsigned int i=0;i<_noise.size();i++){
  //      std::cout<<_noise[i]<<" "; }
  //    std::cout<<std::endl;



  int noise=0;
  //no_hits=_noise.size();

  std::cout<<"NO OF HITS IS "<<_noise.size();
   
  
  for(unsigned int y=0;y< _pointId_to_clusterId.size();++y){

    if  (_pointId_to_clusterId[y]==0) noise++;


  }
  
  std::cout<<" , "<<noise<<" is noise"<<std::endl;
   
  std::cout<<"THE CURRENT NOISE LEVEL IS: "<<(double(noise)/double(_noise.size()))*100<<" %"<<std::endl;
  

}
//-----------------------------------------------------------------


void cluster::DBcluster::produce(edm::Event& evt, edm::EventSetup const&){
   
  //get a collection of clusters
   
  std::auto_ptr<std::vector<recob::Cluster> > ccol(new std::vector<recob::Cluster>);
    
    
    
  //std::cout << "event  : " << evt.Header().Event() << std::endl;
  edm::Service<geo::Geometry> geom;

  
  
  edm::Handle< std::vector<recob::Hit> > hitcol;
  evt.getByLabel(fhitsModuleLabel,hitcol);
  
  
  ///loop over all hits in the event and look for clusters (for each plane)
  
  
    edm::PtrVector<recob::Hit> allhits;
    edm::PtrVector<recob::Hit> clusterHits;

  unsigned int p(0),w(0), channel(0);
  for(int plane=0; plane<geom->Nplanes(); plane++)
    {
      for(unsigned int i = 0; i< hitcol->size(); ++i){
  
	edm::Ptr<recob::Hit> hit(hitcol, i);
  
  
	channel=hit->Wire()->RawDigit()->Channel();
	geom->ChannelToWire(channel,p,w);
    
	if(p == plane) allhits.push_back(hit);
   
      }  
  


      // std::cout<<"number of hits is: "<<hit.size()<<std::endl;
      //  double wire_n;
      //   double hit_time;

      // wire_n=hit[0]->Wire()->RawDigit()->Channel();
      //  hit_time=(hit[0]->StartTime()+hit[0]->EndTime())/2.;
      // std::cout<<"wire# : "<<wire_n<<std::endl;
      //  std::cout<<"time: "<<hit_time<<std::endl;



      // for(unsigned int i=0;i < hit.size(); i++){
      // std::cout<<"hit "<<i<<" : ";
      // wire_n=hit[i]->Wire()->RawDigit()->Channel();
      //  hit_time=(hit[i]->StartTime()+hit[i]->EndTime())/2.;
      //  std::cout<<"wire# : "<<wire_n<<" ";
      //  std::cout<<"time: "<<hit_time<<" ";
      //  std::cout<<std::endl;

      //}
      // std::cout<<"-----------------------------------"<<std::endl;
 
      //------------------------------------------------------------------
      // Determine spacing between wires (different for each detector)
      ///get 2 first wires and find their spacing (wire_dist)
 

      const geo::WireGeo& wire = geom->Plane(0).Wire(0);
      const double pos[3] = {0., 0.0, 0.};
      double posWorld0[3] = {0.};
      double posWorld1[3] = {0.};
      wire.LocalToWorld(pos, posWorld0);

      //std::cerr << "wire 0, plane 0 test pos in world coordinates is (" << posWorld0[0] << "," << posWorld0[1] << "," << posWorld0[2]<< ")" << std::endl;

      const geo::WireGeo& wire1 = geom->Plane(0).Wire(1);
      wire1.LocalToWorld(pos, posWorld1);

      //std::cerr << "wire 1, plane 0 test pos in world coordinates is (" << posWorld1[0] << "," << posWorld1[1] << "," << posWorld1[2]<< ")" << std::endl;
      double wire_dist =posWorld0[1]- posWorld1[1];
      //std::cout<<"wire spacing is: "<< wire_dist<<std::endl;
      //----------------------------------------------------------------
   
      std::vector<std::vector<double> > ps;

      // Filling "Points" with hit info
      //*******************************************************************
      //std::cout<<"in DBSCANfinder.cxx file number of hits="<<hit.size()<<std::endl;
      for (unsigned int j = 0; j < allhits.size(); j++)
	{
	  int dims=3;//our point is defined by 3 elements:wire#,center of the hit, and the hit width
	  std::vector<double> p(dims);
      
     
	  p[0] = (allhits[j]->Wire()->RawDigit()->Channel())*wire_dist;
	  p[1]=((allhits[j]->StartTime()+allhits[j]->EndTime())/2.)*0.03069;
	  p[2]=(allhits[j]->EndTime()-allhits[j]->StartTime())*0.03069;   //width of a hit in cm
	  // p[2]=(allhits[j]->EndTime()-allhits[j]->StartTime());   //width of a hit in cm
	  //std::cout<<"wire= "<<p[0]<<" width= "<< p[2]<<" (this # squared is: "<<p[2]*p[2]<<" )"<<std::endl;
	  fhitwidth->Fill(p[2]);
	  if(allhits[j]->Wire()->RawDigit()->Channel()<240){ fhitwidth_ind_test->Fill(p[2]);}
	  if(allhits[j]->Wire()->RawDigit()->Channel()>240){ fhitwidth_coll_test->Fill(p[2]);}
	  // std::cout<<"Point "<<j<<"= ("<<allhits[j]->Wire()->RawDigit()->Channel()<<", "<<(allhits[j]->StartTime()+allhits[j]->EndTime())/2.<<", "<<p[2]<<")"<<std::endl;
	  //  std::cout << p[i] << ' ';
	  
	  ps.push_back(p);
	  //  std::cout << std::endl;

      
	}
      //*******************************************************************

      cluster::DBScan clusters(ps, fEps,fEps2,fMinPts );

  
      clusters.computeSimilarity();
      clusters.computeSimilarity2();
      clusters.computeWidthFactor();

      clusters.run_cluster();



      //std::cout<<clusters;
      std::cout<<"DBSCAN found "<<clusters._clusters.size()<<" cluster(s)."<<std::endl;


      for(unsigned int i=0; i<clusters._clusters.size();i++)

	{
	  //recob::Cluster* reco_cl= new recob::Cluster();
	  for(unsigned int j=0;j<clusters._pointId_to_clusterId.size();j++){

	    if(clusters._pointId_to_clusterId[j]==(i+1)){

	      // reco_cl->Add(allhits[j]);
	      clusterHits.push_back(allhits[j]);
	    }
       
	  }
     
	  //ccol->push_back(recob::Cluster(clusterHits));
	  //recob::Cluster(clusterHits).SetID(i);
    
    
    
	  ////////
	  if (clusterHits.size()>0)
	    {
	      recob::Cluster cluster(clusterHits);
	      cluster.SetID(i);
	      ccol->push_back(cluster);
	      std::cout<<"no of hits for this cluster is "<<clusterHits.size()<<std::endl;
	      clusterHits.clear();
	      //////
	    }
   
   
   
   
	}

 
      allhits.clear();
    }


  if(ccol->size() == 0){
    std::cerr << "no clusters made for this event" << std::endl;
  }
  evt.put(ccol);
  return;
}






