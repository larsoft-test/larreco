////////////////////////////////////////////////////////////////////////
//
// FeatureVertexFinder class
//
// jasaadi@fnal.gov
//
// This algorithm is designed to reconstruct the vertices using the
// 2D cluster information and CornerFinderAlg to find 2-d and 3-d verticies
//
// We will utilize many of the methods found in VeretxFinder2d but modify them 
// to return vertex candidates and only a point that has matching in at least 
// 2 views will be used
//
//
//
// This is Preliminary Work and needs modifications
//
// ////////////////////////////////////////////////////////////////////////

// ##########################
// ### Basic C++ Includes ###
// ##########################
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <ios>
#include <sstream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <vector>

// ##########################
// ### Framework Includes ###
// ##########################
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ########################
// ### LArSoft Includes ###
// ########################
#include "RecoBase/EndPoint2D.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Track.h"
#include "RecoBase/Shower.h"
#include "RecoBase/Vertex.h"
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/AssociationUtil.h"
#include "RecoAlg/CornerFinderAlg.h"

// #####################
// ### ROOT Includes ###
// #####################
#include "TMath.h"
#include "TH1D.h"
#include "TVectorD.h"
#include "TGeoManager.h"
#include "TMath.h"
#include "TGraph.h"
#include "TF1.h"
#include "TVector3.h"

// =====================================================================================================
// =====================================================================================================
#ifndef FeatureVertexFinder_H
#define FeatureVertexFinder_H

class TH1D;

///vertex reconstruction
namespace vertex {
   
 class FeatureVertexFinder :  public art::EDProducer {
    
  public:
    
    explicit FeatureVertexFinder(fhicl::ParameterSet const& pset); 
    virtual ~FeatureVertexFinder();        
    void beginJob();
    void reconfigure(fhicl::ParameterSet const& p);

    
    void produce(art::Event& evt);

  private:
    
    TH1D *dtIC;
  
    std::string fClusterModuleLabel;
    // JA:Adding Corner Finder Alg
    cluster::CornerFinderAlg  fCorner;

  };
    
}//<---End namespace vertex

#endif // FeatureVertexFinder_H
// =====================================================================================================
// =====================================================================================================








// =====================================================================================================
// =====================================================================================================
namespace vertex{

//-----------------------------------------------------------------------------
  FeatureVertexFinder::FeatureVertexFinder(fhicl::ParameterSet const& pset):
  fCorner(pset.get<fhicl::ParameterSet>("CornerPset"))
  {  
    /*this->*/reconfigure(pset);    
    produces< std::vector<recob::Vertex> >();
    produces< std::vector<recob::EndPoint2D> >();
    produces< art::Assns<recob::EndPoint2D, recob::Hit> >();
    produces< art::Assns<recob::Vertex, recob::Hit> >();
    produces< art::Assns<recob::Vertex, recob::Shower> >();
    produces< art::Assns<recob::Vertex, recob::Track> >();
  }
//-----------------------------------------------------------------------------
  FeatureVertexFinder::~FeatureVertexFinder()
  {
  }

  //---------------------------------------------------------------------------
  void FeatureVertexFinder::reconfigure(fhicl::ParameterSet const& p) 
  {
    fClusterModuleLabel  = p.get< std::string >("ClusterModuleLabel");
    return;
  }
  //-------------------------------------------------------------------------
  void FeatureVertexFinder::beginJob(){
    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;
    dtIC = tfs->make<TH1D>("dtIC","It0-Ct0",100,-5,5);
    dtIC->Sumw2();
  }

// //-----------------------------------------------------------------------------
  void FeatureVertexFinder::produce(art::Event& evt)
  {
  
  
    std::cout<<"Top of FeatureVertexFinder Produce"<<std::endl;  
    // #########################
    // ### Geometry Services ###
    // #########################
    art::ServiceHandle<geo::Geometry> geom;
    
    // ###############################
    // ### LAr Properties Services ###
    // ###############################
    art::ServiceHandle<util::LArProperties> larprop;
    
    // ####################################
    // ### Detector Properties Services ###
    // ####################################
    art::ServiceHandle<util::DetectorProperties> detprop;
    
    
    // ###############################
    // ### Defining TPC Parameters ###
    // ###############################
    // === TPC Name === (JA: Not sure if I need this)
    TString tpcName = geom->GetLArTPCVolumeName();
    // === Length of the detector ===
    double DetectorLength =  (geom->DetHalfHeight())*2.;
    // === wire angle with respect to the vertical direction ===
    double Angle = geom->Plane(1).Wire(0).ThetaZ(false)-TMath::Pi()/2.; 
    
    // ==================================================
    // === Calculating the Timetick to CM conversion  ===
    // ==================================================
    double TimeTick = detprop->SamplingRate()/1000.; //<---To get units of microsecond...not nanosec
    double Efield_drift = larprop->Efield();      // Electric Field in the drift region in kV/cm
    double Temperature  = larprop->Temperature(); // LAr Temperature in K
    // === Drift velocity in the drift region (cm/us) ===
    double DriftVelocity = larprop->DriftVelocity(Efield_drift,Temperature); 
    // === Time sample (cm) === 
    double TimetoCm = DriftVelocity*TimeTick; 
    
    
    double presamplings = detprop->TriggerOffset(); //trigger offset
    
    // =============================================
    // === Wire Pitch variable (to be set later) ===
    // =============================================
    double wire_pitch   = 0;
    // ===========================================================
    // === Number of planes in this detector (to be set later) ===
    // ===========================================================
    int nplanes = 0;
    
    std::vector<double> vtx_wire = {0.};
    std::vector<double> vtx_time = {0.};
    std::vector<double> vtx_plane = {0.};
    
    
    
    
    
//-------------------------------------------------------------------------------------------------------------------------------------------------   
    // ####################################################################
    // ### Utilizing the CornerFinderAlg to get a list of EndPoint2d's  ###
    // ####################################################################
    
    std::vector<double>   feature_wire = {0.};
    std::vector<double>   feature_time = {0.};
    std::vector<double>   feature_plane = {0.};
    std::vector<uint32_t> feature_channel =  {0.};
    int nFeatures = 0;
    
    double x_feature[1000] = {0.}, y_feature[1000] = {0.}, z_feature[1000] = {0.};
    double yy = 0., zz = 0.;
    int n3dFeatures = 0;
    
    // ###################################################
    // ### Take in the raw information about the event ###
    // ###################################################
    fCorner.TakeInRaw(evt);
    
    // #######################################################
    // ### Push the features into a vector of Endpoint2d's ###
    // #######################################################
    std::vector<recob::EndPoint2D> EndPoints;
    fCorner.get_feature_points(EndPoints);
    
    // ########################################################
    // ### Loop over all the features and record their info ###
    // ########################################################
    for(size_t i=0; i!=EndPoints.size(); ++i)
    	{
	
	feature_wire.push_back(EndPoints.at(i).WireID().Wire);
	feature_time.push_back(EndPoints.at(i).DriftTime());
	feature_plane.push_back(EndPoints.at(i).WireID().Plane);
	feature_channel.push_back(geom->PlaneWireToChannel(EndPoints.at(i).WireID().Plane, EndPoints.at(i).WireID().Wire, 0, 0));
	nFeatures++;
	
	//std::cout<<"EndPoints View() = "<<EndPoints.at(i).WireID().Plane<<std::endl;
	//std::cout<<"EndPoints DriftTime() = "<<EndPoints.at(i).DriftTime()<<std::endl;
	//std::cout<<"EndPoints Wire = "<<EndPoints.at(i).WireID().Wire<<std::endl;
	//std::cout<<"EndPoints Channel = "<<geom->PlaneWireToChannel(EndPoints.at(i).WireID().Plane, EndPoints.at(i).WireID().Wire, 0, 0)<<std::endl;
        //EndPointTimesInPlane[EndPoints.at(i).View()].push_back(EndPoints.at(i).DriftTime());
	
	
       	}//<---End i loop finding 2d Features
     
    // ##############################################################################
    // ### Now loop over those features and see if they match up in between views ###
    // ##############################################################################
    
    for(int feature1 = nFeatures; feature1 > 0; feature1 --)
    	{
	for(int feature2 = 0; feature2 < feature1; feature2++)
		{
		// ##########################################################################
		// ### Check to make sure we are comparing features from different planes ###
		// ##########################################################################
		if(feature_plane[feature1] != feature_plane[feature2])
			{
			// ####################################################################
			// ### Checking to see if these features have intersecting channels ###
			// ####################################################################
			if( geom->ChannelsIntersect( feature_channel[feature2], feature_channel[feature1], yy, zz) &&
			    std::abs(feature_time[feature2] - feature_time[feature1]) < detprop->TimeOffsetV())
			    	{
				x_feature[n3dFeatures] = detprop->ConvertTicksToX(feature_time[feature1], feature_plane[feature1], 0, 0); //<---Hardcoding tpc and cryostat for now.
				y_feature[n3dFeatures] = yy;
				z_feature[n3dFeatures] = zz;
				n3dFeatures++;
				
				//std::cout<<std::endl;
				//std::cout<<"### Found a 3d Feature ###"<<std::endl;
				//std::cout<<"X = "<<x_feature<<" , Y = "<<y_feature<<" , Z = "<<z_feature<<std::endl;
				//std::cout<<std::endl;
				
	
				}//<---End Checking if the feature matches in 3d
			}//<---End Checking features are in different planes
		}//<---End feature2 loop
	}//<---End feature1 loop

 
//-------------------------------------------------------------------------------------------------------------------------------------------------    
    
    //std::cout<<"Making new collections for output"<<std::endl;
    //Point to a collection of vertices to output.
    std::unique_ptr<std::vector<recob::Vertex> >                 vcol(new std::vector<recob::Vertex>);          //3D vertex
    std::unique_ptr<std::vector<recob::EndPoint2D> >             epcol(new std::vector<recob::EndPoint2D>);  //2D vertex
    std::unique_ptr< art::Assns<recob::EndPoint2D, recob::Hit> > assnep(new art::Assns<recob::EndPoint2D, recob::Hit>);
    std::unique_ptr< art::Assns<recob::Vertex, recob::Shower> >  assnsh(new art::Assns<recob::Vertex, recob::Shower>);
    std::unique_ptr< art::Assns<recob::Vertex, recob::Track> >   assntr(new art::Assns<recob::Vertex, recob::Track>);
    std::unique_ptr< art::Assns<recob::Vertex, recob::Hit> >     assnh(new art::Assns<recob::Vertex, recob::Hit>);
    
    
    
    // ###################################################
    // ### Retreiving the Cluster Module for the event ###
    // ###################################################
    art::Handle< std::vector<recob::Cluster> > clusterListHandle;
    evt.getByLabel(fClusterModuleLabel,clusterListHandle);
    
    // ##################################
    // ### Filling the Cluster Vector ###
    // ##################################
    art::PtrVector<recob::Cluster> clusters;
    for (unsigned int ii = 0; ii <  clusterListHandle->size(); ++ii){
      art::Ptr<recob::Cluster> clusterHolder(clusterListHandle,ii);
      clusters.push_back(clusterHolder);
    }//<---End ii loop
    
    // #################################################
    // ### Finding hits associated with the clusters ###
    // #################################################
    art::FindManyP<recob::Hit> fmh(clusterListHandle, evt, fClusterModuleLabel);
    
    // ==============================================================================
    // ===  Loop over all the clusters that are found in the cluster list and fit ===
    // ==  the hits in those clusters to establish a vector of slopes (dtdwstart) ===
    // === that has the same number and order as our cluster PtrVector (clusters) ===
    // ==============================================================================
    std::vector<int> Cls[3]; //<---- Index to clusters in each view
    std::vector<double> dtdwstart; //<----Slope (delta Time Tick vs delta Wire) 
    std::vector<double> dtdwstartError; //<---Error on the slope    
    
    //std::cout<<"Getting ready to loop over clusters"<<std::endl;
    for(size_t iclu = 0; iclu < clusters.size(); ++iclu)
    	{
	//std::cout<<"Looping over the first cluster"<<std::endl;
	// ######################################################
	// ### Determine which view the current cluster is in ###
	// ######################################################
	switch(clusters[iclu]->View())
		{
	    	//std::cout<<"Inside the switch loop"<<std::endl;
	  	case geo::kU :
	    	  Cls[0].push_back(iclu);
	          break;
	  	case geo::kV :
	  	  Cls[1].push_back(iclu);
	  	  break;
	  	case geo::kZ :
	  	  Cls[2].push_back(iclu);
	  	  break;
	  	default :
	  	  break;
	  	}
	  
	// #############################################################
	// ### Filling wires and times into a TGraph for the cluster ###
	// #############################################################	  
	std::vector<double> wires;
	std::vector<double> times;
	  
	// ### Gathering the hits associated with the current cluster ###
	std::vector< art::Ptr<recob::Hit> > hit = fmh.at(iclu);
	
	// ### Counting the number of hits in the current cluster (n) ###
	int n = 0;
	// ############################################
	// ### Looping over the hits in the cluster ###
	// ############################################
	//std::cout<<"Getting ready to loop over this clusters hits"<<std::endl;
	for(size_t i = 0; i < hit.size(); ++i)
		{
		// ### JA: Should I convert wire and time to cm? Will that help the fitting? ###
		
	    	wires.push_back(hit[i]->WireID().Wire);
	    	times.push_back(hit[i]->PeakTime());
	    	++n;
	  	}//<---End loop over hits (i)
	// ################################################################
	// ### If there are 2 or more hits in the cluster fill a TGraph ###
	// ###         and fit a from a polynomial or order 1           ###
	// ################################################################
	if(n>=2)
		{
		// ###################################
		// ### Push the hits into a TGraph ###
		// ###################################
	    	TGraph *the2Dtrack = new TGraph(n,&wires[0],&times[0]);  
		// === Try to fit the TGraph with a 1st order polynomial ===
	    	try
			{
	      		the2Dtrack->Fit("pol1","Q");
	      		TF1 *pol1=(TF1*) the2Dtrack->GetFunction("pol1");
	      		double par[2];
			double parerror[2];
	      		pol1->GetParameters(par);
			parerror[1] = pol1->GetParError(1);
	      		//std::cout<<iclu<<" "<<par[1]<<" "<<clusters[iclu]->dTdW()<<std::endl;
	      
	      		// #######################################################################
	      		// ### Take change in time tick vs change in wire (dT/dW) from the fit ###
	      		// #######################################################################
	      		dtdwstart.push_back(par[1]);
			dtdwstartError.push_back(parerror[1]);
	    		}//<---End try to fit with a polynomial order 1
			
	    	// ############################################################
	    	// ### If the fitter fails just take dT/dW from the cluster ###
	    	// ############################################################
	    	catch(...)
			{
	      		mf::LogWarning("FeatureVertexFinder") << "Fitter failed, using the clusters default dTdW()";
	      		delete the2Dtrack;
	      		dtdwstart.push_back(clusters[iclu]->dTdW());
	      		continue;
	    		}
	    	delete the2Dtrack;
	  	}//<---End if the cluster has 2 or more hits
	// #################################################
	// ### If the cluster has fewer than 2 hits just ### 
	// ###      take the dT/dW from the cluster      ###
	// #################################################
	else {dtdwstart.push_back(clusters[iclu]->dTdW());}
	}//<---End loop over clusters iclu
	
	
    
    
    
    
    // ##########################################
    // ### Variables to save for each cluster ###
    // ##########################################
    int   nClustersFound = 0;				//<---Number of clusters to be evaluated
    int   n2dVertexCandidates = 0;			//<---Number of candidate 2d Vertices found
    float Clu_Plane[1000] = {0};         	     	//<---Plane of the current cluster
    float Clu_StartPos_Wire[1000]= {0};     		//<---Starting wire number of cluster
    float Clu_StartPosUncer_Wire[1000]= {0};		//<---Starting wire number uncertainty
    float Clu_StartPos_TimeTick[1000]= {0};      	//<---Starting TDC value of the cluster
    float Clu_StartPosUncer_TimeTick[1000]= {0}; 	//<---Starting TDC value uncertainty
  
    float Clu_EndPos_Wire[1000]= {0};       		//<---Ending wire number of cluster
    float Clu_EndPosUncer_Wire[1000]= {0};  		//<---Ending wire number uncertainty
    float Clu_EndPos_TimeTick[1000]= {0};        	//<---Ending TDC value of the cluster
    float Clu_EndPosUncer_TimeTick[1000]= {0};   	//<---Ending TDC value uncertainty
  
    float Clu_Slope[1000]= {0};	  	   		//<---Calculated Slope of the cluster (TDC/Wire)
    float Clu_SlopeUncer[1000]= {0};        		//<---Slope Error
    float Clu_Yintercept[1000]= {0};			//<---Clusters Y Intercept
    float Clu_Length[1000]= {0};			//<---Calculated Length of the cluster
    
    
    // ##############################
    // ### Looping over cryostats ###
    // ##############################
    for(size_t cstat = 0; cstat < geom->Ncryostats(); ++cstat)
    	{
    	// ##########################
      	// ### Looping over TPC's ###
      	// ##########################
      	for(size_t tpc = 0; tpc < geom->Cryostat(cstat).NTPC(); ++tpc)
		{
		nplanes = geom->Cryostat(cstat).TPC(tpc).Nplanes();
		
		// === Defining some variables to be used later JA: Update this ===
		std::vector<int> cluvtx[nplanes];
    		
    		// #################################
		// ### Loop over the wire planes ###
		// #################################
		for (int i = 0; i < geom->Cryostat(cstat).TPC(tpc).Nplanes(); ++i)
			{
    			//           geom->WirePitch(Wire1, Wire2, Plane#, TPC#, Cyro#);
    			wire_pitch = geom->WirePitch(0,1,i,tpc,cstat);
    			
			// ##############################################
	  		// ### If there is at least one cluster found ###
	  		// ##############################################
	  		if (Cls[i].size() >= 1)
				{
				// ##############################
	    			// ### Loop over each cluster ###
	    			// ##############################
	    			for (unsigned j = 0; j<Cls[i].size(); ++j)
					{
					
					// === Current Clusters Plane ===
					Clu_Plane[nClustersFound]  		  = clusters[Cls[i][j]]->View();
					//std::cout<<"Clu_Plane[nClustersFound] = "<<Clu_Plane[nClustersFound]<<std::endl;
					
					// === Current Clusters StartPos ===
					Clu_StartPos_Wire[nClustersFound]	  = clusters[Cls[i][j]]->StartPos()[0];
					Clu_StartPosUncer_Wire[nClustersFound]	  = clusters[Cls[i][j]]->SigmaStartPos()[0];
					Clu_StartPos_TimeTick[nClustersFound]	  = clusters[Cls[i][j]]->StartPos()[1];
					Clu_StartPosUncer_TimeTick[nClustersFound]= clusters[Cls[i][j]]->SigmaStartPos()[1];
					
					// === Current Clusters EndPos ===
					Clu_EndPos_Wire[nClustersFound]		  = clusters[Cls[i][j]]->EndPos()[0];
					Clu_EndPosUncer_Wire[nClustersFound]	  = clusters[Cls[i][j]]->SigmaEndPos()[0];
					Clu_EndPos_TimeTick[nClustersFound]	  = clusters[Cls[i][j]]->EndPos()[1];
					Clu_EndPosUncer_TimeTick[nClustersFound]  = clusters[Cls[i][j]]->SigmaEndPos()[1];
					
					// === Current Clusters Slope (In Wire and Time Tick)
					Clu_Slope[nClustersFound] 		  = dtdwstart[Cls[i][j]];
					Clu_SlopeUncer[nClustersFound]		  = dtdwstartError[Cls[i][j]];
					
					
					Clu_Length[nClustersFound] 		  = std::sqrt(pow((clusters[Cls[i][j]]->StartPos()[0]-clusters[Cls[i][j]]->EndPos()[0])*13.5,2)
					 						    +pow(clusters[Cls[i][j]]->StartPos()[1]-clusters[Cls[i][j]]->EndPos()[1],2));

					// ######################################################
					// ### Given a slope and a point find the y-intercept ###
					// ###                   c = y-mx                     ###
					// ######################################################
					Clu_Yintercept[nClustersFound] = Clu_StartPos_TimeTick[nClustersFound] - (Clu_Slope[nClustersFound] * Clu_StartPos_Wire[nClustersFound]);
					
					
					nClustersFound++;
    					}//<---End looping over clusters
					
					
    				}//<---End checking that we have at least one cluster found
			// ################################################################
			// ## If no clusters were found then put in dummy vertex values ###
			// ################################################################
			else
				{
	    			vtx_wire.push_back(-1);
	    			vtx_time.push_back(-1);
				vtx_plane.push_back(-1);
	  			}//<---End no clusters found else statement
			
			
			
			
			// ################################################################################
			// ### Now we try to find a 2-d vertex in the plane we are currently looking in ###
			// ################################################################################
			for (unsigned int n = nClustersFound; n > 0; n--)
				{
				// #######################################################
				// ###   Looping over the clusters starting from the   ###
				// ### first cluster and checking against the nCluster ###
				// #######################################################
				for (unsigned int m = 0; m < n; m++)
					{
					// ###########################################################
					// ### Checking to make sure clusters are in the same view ###
					// ###########################################################
					if(Clu_Plane[n] == Clu_Plane[m])
						{
						// --- Skip the vertex if the lines slope don't intercept ---
						if(Clu_Slope[m] - Clu_Slope[n] == 0){break;}
						
						
						// ============================================================
						// === X intersection = (yInt2 - yInt1) / (slope1 - slope2) ===
						float intersection_X = (Clu_Yintercept[n] - Clu_Yintercept[m]) / (Clu_Slope[m] - Clu_Slope[n]);
						// ================================================
						// === Y intersection = (slope1 * XInt) + yInt1 ===
						float intersection_Y = (Clu_Slope[m] * intersection_X) + Clu_Yintercept[m];
						
						
						//if(intersection_X < 1 || intersection_Y < 1) {continue;}
						//if( intersection_X > geom->Nwires(Clu_Plane[n],0,0) || 
						//intersection_X > geom->Nwires(Clu_Plane[m],0,0) || intersection_Y > detprop->NumberTimeSamples() ) {break;}
						
						// ##########################################################
						// ### Filling the vector of Vertex Wire, Time, and Plane ###
						// ##########################################################
						
						// -----------------------------------------------------------------------------
						// --- Skip this vertex if the X  and Y intersection is outside the detector ---
						// --- using geom->Nwires(plane,tpc,cyrostat) & detprop->NumberTimeSamples() ---
						// -----------------------------------------------------------------------------
						if( intersection_X > 1 && intersection_Y > 0 && 
						    (intersection_X < geom->Nwires(Clu_Plane[n],0,0) || intersection_X < geom->Nwires(Clu_Plane[m],0,0) ) &&
						    intersection_Y < detprop->NumberTimeSamples() )
							{
							/*std::cout<<std::endl;
							std::cout<<"Vertex Wire = "<<intersection_X<<" , Vertex Time Tick = "<<intersection_Y<<"  , Vertex Plane # = "<<Clu_Plane[m]<<std::endl;
							std::cout<<"geom->Nwires(Clu_Plane[n],0,0) = "<<geom->Nwires(Clu_Plane[n],0,0)<<std::endl;
							std::cout<<"geom->Nwires(Clu_Plane[m],0,0) = "<<geom->Nwires(Clu_Plane[m],0,0)<<std::endl;
							std::cout<<std::endl;*/
							
							vtx_wire.push_back(intersection_X);
							vtx_time.push_back(intersection_Y);
							vtx_plane.push_back(Clu_Plane[m]);
							n2dVertexCandidates++;
							}//<---End saving a "good 2d vertex" candidate
						
						}//<---End making sure we are in the same plane
					}//<---End m ++ loop
				}//<--- End n-- loop
			
			// === Zeroing the clusters for the current view
			nClustersFound = 0;
    			}//<---End loop over wireplanes
    		}//<---End looping over TPC's
    	}//<---End looping over cryostats
    
    double y_coord = 0, z_coord = 0;
    
    double x_3dVertex[1000] = {0.}, y_3dVertex[1000] = {0.}, z_3dVertex[1000] = {0.};
    int n3dVertex = 0;
    // --------------------------------------------------------------------------
    // ---   Having now found a very long list of potential 2-d end points    ---
    // --- we need to check if any of them match between planes and only keep ---
    // ---                       those that have matches                      ---
    // --------------------------------------------------------------------------
    for(unsigned int vtx = n2dVertexCandidates; vtx > 0; vtx--)
    	{
	for (unsigned int vtx1 = 0; vtx1 < vtx; vtx1++)
		{
		// ###########################################################################
		// ### Check to make sure we are comparing verticies from different planes ###
		// ###########################################################################
		if(vtx_plane[vtx1] != vtx_plane[vtx])
			{
			/*std::cout<<std::endl;
			std::cout<<"vtx_wire[vtx1] = "<<vtx_wire[vtx1]<<" , vtx_wire[vtx] = "<<vtx_wire[vtx]<<std::endl;
			std::cout<<"vtx_time[vtx1] = "<<vtx_time[vtx1]<<" , vtx_time[vtx] = "<<vtx_time[vtx]<<std::endl;
			std::cout<<"vtx_plane[vtx1] = "<<vtx_plane[vtx1]<<" , vtx_plane[vtx] = "<<vtx_plane[vtx]<<std::endl;
			std::cout<<std::endl;*/
			
			// === To figure out if these two verticies are from a common point 
			// === we need to check if the channels intersect and if they are 
			// === close in time ticks as well...to do this we have to do some 
			// === converting to use geom->PlaneWireToChannel(PlaneNo, Wire, tpc, cstat)
			// === JA: Need to include vtx tpc, and cstat to make detector agnositc
			
			unsigned int vtx1_plane   = vtx_plane[vtx1];
			unsigned int vtx1_wire    = vtx_wire[vtx1];
			uint32_t     vtx1_channel = geom->PlaneWireToChannel(vtx1_plane, vtx1_wire, 0, 0);
			
			unsigned int vtx2_plane   = vtx_plane[vtx];
			unsigned int vtx2_wire    = vtx_wire[vtx];
			uint32_t     vtx2_channel = geom->PlaneWireToChannel(vtx2_plane, vtx2_wire, 0, 0);
			
			// ##############################################################################
			// ### Check to see if the channels intersect and save the y and z coordinate ###
			// ##############################################################################
			bool match = false;
			
			try
				{
				match = geom->ChannelsIntersect( vtx1_channel, vtx2_channel, y_coord, z_coord);
				
				
				}
			catch(...)
				{
	      			mf::LogWarning("FeatureVertexFinder") << "match failed for some reason";
	      
	      			match = false;
	      			continue;
	    			}
			if( match )
				{
				//std::cout<<std::endl;
				//std::cout<<"Vertex Intersection Found!!!"<<std::endl;
				//std::cout<<"Offset between planes = "<<detprop->TimeOffsetV()<<std::endl;
				//std::cout<<std::endl;
				
				// #############################################
				// ### Now check if the matched channels are ###
				// ###   within 2*TimeOffset between planes  ###
				// #############################################
				if(std::abs(vtx_time[vtx1] - vtx_time[vtx]) < 2*detprop->TimeOffsetV())
					{
					//        detprop->ConvertTicksToX(ticks, plane, tpc, cryostat)
					x_3dVertex[n3dVertex] = detprop->ConvertTicksToX(vtx_time[vtx], vtx_plane[vtx], 0, 0); //<---Hardcoding tpc and cryostat for now.
					y_3dVertex[n3dVertex] = y_coord;
					z_3dVertex[n3dVertex] = z_coord;
					n3dVertex++;
					
					
					//std::cout<<std::endl;
					//std::cout<<" ============== REAL 3-D VERTEX FOUND ============"<<std::endl;
					//std::cout<<" std::abs(vtx_time[vtx1] - vtx_time[vtx] = "<<std::abs(vtx_time[vtx1] - vtx_time[vtx])<<std::endl;
					//std::cout<<"x_coord = "<<x_coord<<" , y_coord = "<<y_coord<<" , z_coord = "<<z_coord<<std::endl;
					//std::cout<<std::endl;
					
					
					
					
					}//<---End Checking if the vertices agree "well enough" in time tick
				
				}//<---End Checking if verticies intersect

			}//<--- End checking we are in different planes
		}//<---end vtx1 for loop
	}//<---End vtx for loop



       	// ######################################################################################################
       	// ##################### Now I am going to deal with the 3 cases for finding 3d verticies ###############
       	// ### Case1) If n3dVertex == 1: Then I only found one intersection point using the clusters found,
       	// ###                           in this case I will simply save to the event a EndPoint2d and a 3d
       	// ###                           vertex (JA: Need to add relevant associations later)
       	// ### Case2) If n3dVertex > 1: In this case I will attempt to find a 3d feature that is within 2cm
       	// ###                          in x, y, and z. If there is more than one of these points that match
       	// ###                          I will, for now (JA: Need to improve), take the one with the smaller 
       	// ###                          z coordinate (more upstream)
       	// ###
       	// ### Case3) If n3dVertex = 0: Still trying to figure out what to do here....JA:(will come back to)
       	// #######################################################################################################
       

//-----------------------------------------------------------------------------------------------------------------------------       
       	// ### Case 1...only one 3d vertex found ###
       	if (n3dVertex == 1)
		{
		// ##############################
    		// ### Looping over cryostats ###
    		// ##############################
    		for(size_t cstat = 0; cstat < geom->Ncryostats(); ++cstat)
    			{
    			// ##########################
      			// ### Looping over TPC's ###
      			// ##########################
      			for(size_t tpc = 0; tpc < geom->Cryostat(cstat).NTPC(); ++tpc)
				{
    				// #################################
				// ### Loop over the wire planes ###
				// #################################
				for (int i = 0; i < geom->Cryostat(cstat).TPC(tpc).Nplanes(); ++i)
					{
					double xyz[3] = {x_3dVertex[0], y_3dVertex[0], z_3dVertex[0]};
					double EndPoint2d_TimeTick = detprop->ConvertXToTicks(x_3dVertex[0],i, tpc, cstat);
					
					int EndPoint2d_Wire 	   = geom->NearestWire(xyz , i, tpc, cstat);
					int EndPoint2d_Channel     = geom->NearestChannel(xyz, i, tpc, cstat);
					geo::View_t View	   = geom->View(EndPoint2d_Channel);
					geo::WireID wireID(cstat,tpc,i,EndPoint2d_Wire);
					
					
					// ### Saving the 2d Vertex found ###
					recob::EndPoint2D vertex( EndPoint2d_TimeTick , //<---TimeTick
								  wireID ,		//<---geo::WireID
								  1 ,			//<---Vtx strength (JA: ?)
								  epcol->size() ,	//<---Vtx ID (JA: ?)
								  View ,		//<---Vtx View 	
								  1 );			//<---Vtx Strength (JA: Need to figure this one?)
					epcol->push_back(vertex);	
					
					
							   
					}//<---End loop over Planes
				}//<---End loop over tpc's
			}//<---End loop over cryostats
			
		// ############################
		// ### Saving the 3d vertex ###
		// ############################
		double xyz2[3] = {x_3dVertex[0], y_3dVertex[0], z_3dVertex[0]};
		recob::Vertex the3Dvertex(xyz2, vcol->size());
		vcol->push_back(the3Dvertex);
		}//<---End Case 1, only one 3d Vertex found

//----------------------------------------------------------------------------------------------------------------------------- 
	// ###      Case 2...multiple 3d verticies found       ###
	// ### To break the degeneracy start by looping over   ###
	// ###    the 3d verticies found, then loop over the   ###
	// ### 3d feature points found only keeping the points ###
	// ###       that are within 2 cm of the vertex        ###
	
	int vertexNumber = -1;
	double x_matched = 999, y_matched = 999, z_matched = 999;
	if (n3dVertex > 1)
		{
		// ##############################
		// ### Looping over verticies ###
		// ##############################
		for(int l = 0; l < n3dVertex; l++)
			{
			// #############################
			// ### Looping over features ###
			// #############################
			for (int f = 0; f < n3dFeatures; f++)
				{
				// ###########################################################
				// ### Looking for features and verticies within 3cm in 3d ###
				// ###########################################################
				if(std::abs(x_3dVertex[l] - x_feature[f]) <= 3 &&  
				   std::abs(y_3dVertex[l] - y_feature[f]) <= 3 &&
				   std::abs(z_3dVertex[l] - z_feature[f]) <= 3)
				   	{
					// ### Put in a catch against multiple matches ###
					if( vertexNumber != l && z_3dVertex[l] < z_matched) 
						{
						vertexNumber = l;
						x_matched = x_3dVertex[l];
						y_matched = y_3dVertex[l];
						z_matched = z_3dVertex[l];
						std::cout<<std::endl;
						std::cout<<" ###########################################"<<std::endl;
						std::cout<<" ### 3d match between feature and vertex ###"<<std::endl;
						std::cout<<std::endl;
						std::cout<<"x_3dVertex[l] = "<<x_3dVertex[l]<<" , y_3dVertex[l] = "<<y_3dVertex[l]<<" , z_3dVertex[l] = "<<z_3dVertex[l]<<std::endl;
						std::cout<<std::endl;
						
						}//<---End catch against multiple matches
					
		
					}//<---End finding a match between features and verticies 
				}//<--- End f for loop
			}//<--End l for loop
			
		// ##############################
    		// ### Looping over cryostats ###
    		// ##############################
    		for(size_t cstat = 0; cstat < geom->Ncryostats(); ++cstat)
    			{
    			// ##########################
      			// ### Looping over TPC's ###
      			// ##########################
      			for(size_t tpc = 0; tpc < geom->Cryostat(cstat).NTPC(); ++tpc)
				{
    				// #################################
				// ### Loop over the wire planes ###
				// #################################
				for (int i = 0; i < geom->Cryostat(cstat).TPC(tpc).Nplanes(); ++i)
					{
					double xyz[3] = {x_matched, y_matched, z_matched};
					double EndPoint2d_TimeTick = detprop->ConvertXToTicks(x_matched,i, tpc, cstat);
					
					int EndPoint2d_Wire 	   = geom->NearestWire(xyz , i, tpc, cstat);
					int EndPoint2d_Channel     = geom->NearestChannel(xyz, i, tpc, cstat);
					geo::View_t View	   = geom->View(EndPoint2d_Channel);
					geo::WireID wireID(cstat,tpc,i,EndPoint2d_Wire);
					
					
					// ### Saving the 2d Vertex found ###
					recob::EndPoint2D vertex( EndPoint2d_TimeTick , //<---TimeTick
								  wireID ,		//<---geo::WireID
								  1 ,			//<---Vtx strength (JA: ?)
								  epcol->size() ,	//<---Vtx ID (JA: ?)
								  View ,		//<---Vtx View 	
								  1 );			//<---Vtx Strength (JA: Need to figure this one?)
					epcol->push_back(vertex);	
					
					
							   
					}//<---End loop over Planes
				}//<---End loop over tpc's
			}//<---End loop over cryostats
			
		// ############################
		// ### Saving the 3d vertex ###
		// ############################
		double xyz2[3] = {x_matched, y_matched, z_matched};
		recob::Vertex the3Dvertex(xyz2, vcol->size());
		vcol->push_back(the3Dvertex);	
		
		}//<---End Case 2, many 3d Verticies found
		
		
		
//----------------------------------------------------------------------------------------------------------------------------- 
	
	// ### This is for the case where we didn't find a good match
	// ### in 3d for a vertex...still need to figure out what to do here
	if (n3dVertex == 0)
		{
		
		
		
		
		
		
		}	
       	
    
//----------------------------------------------------------------------------------------------------------------------------- 
    

    mf::LogVerbatim("Summary") << std::setfill('-') << std::setw(175) << "-" << std::setfill(' ');
    mf::LogVerbatim("Summary") << "FeatureVertexFinder Summary:";
    for(size_t i = 0; i<epcol->size(); ++i) mf::LogVerbatim("Summary") << epcol->at(i) ;
    for(size_t i = 0; i<vcol->size(); ++i) mf::LogVerbatim("Summary") << vcol->at(i) ;
    
    evt.put(std::move(epcol));
    evt.put(std::move(vcol));
    evt.put(std::move(assnep));
    evt.put(std::move(assntr));
    evt.put(std::move(assnsh));
    evt.put(std::move(assnh));

  } // end of produce
}


// =====================================================================================================
// =====================================================================================================






// =====================================================================================================
// =====================================================================================================
namespace vertex{

  DEFINE_ART_MODULE(FeatureVertexFinder);

} // end of vertex namespace
// =====================================================================================================
// =====================================================================================================

