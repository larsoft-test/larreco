////////////////////////////////////////////////////////////////////////
//
// CornerFinderAlg class
//
// yo-mama@state.edu
//
//  blah blah blah...put description here 
//  
//  
//  
////////////////////////////////////////////////////////////////////////


// ### Lots of includes...not sure if we need them all...just copy and pasting
#include <iostream>
#include <vector>
extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}
#include <sstream>
#include <fstream>
#include <math.h>
#include <algorithm>

#include "TMath.h"

// Framework includes
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "RawData/raw.h"
#include "RawData/RawDigit.h"
#include "Filters/ChannelFilter.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/Geometry.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Wire.h"
#include "ClusterFinder/CornerFinderAlg.h"



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


// NOTE: In the .h file I assumed this would belong in the cluster class....if 
// we decide otherwise we will need to search and replace for this


//-----------------------------------------------------------------------------
cluster::CornerFinderAlg::CornerFinderAlg(fhicl::ParameterSet const& pset) 
// ### This is where we take in the RawData from the .fcl file
  //:fRawDataLabel (pset.get<std::string>("RawDataModuleLabel"))
// ### Reading in CaldataModule ###
  //:fCalDataModuleLabel (pset.get<std::string>("CalDataModuleLabel")) 
{
  this->reconfigure(pset);
}

//-----------------------------------------------------------------------------
cluster::CornerFinderAlg::~CornerFinderAlg()
{
  delete WireData_histos;
  delete WireData_IDs;
}

//-----------------------------------------------------------------------------
void cluster::CornerFinderAlg::reconfigure(fhicl::ParameterSet const& p)
{
  // ### These are all the tuneable .fcl file parameters from the event ###
  fRawDataModuleLabel			 = p.get< std::string	 >("RawDataModuleLabel");
  fCalDataModuleLabel  			 = p.get< std::string 	 >("CalDataModuleLabel");
  fConversion_threshold     		 = p.get< float    	 >("Conversion_threshold");
  fConversion_bins_per_input_x  	 = p.get< int      	 >("Conversion_bins_per_input_x");
  fConversion_bins_per_input_y       	 = p.get< int      	 >("Conversion_bins_per_input_y");
  fDerivative_method        		 = p.get< std::string    >("Derivative_method");
  fCornerScore_neighborhood     	 = p.get< int		 >("CornerScore_neighborhood");
  fCornerScore_algorithm		 = p.get< std::string    >("CornerScore_algorithm");
  fMaxSuppress_neighborhood		 = p.get< int		 >("MaxSuppress_neighborhood");
  fMaxSuppress_threshold		 = p.get< int		 >("MaxSuppress_threshold");
  fIntegral_bin_threshold                = p.get< float          >("Integral_bin_threshold");
  fIntegral_fraction_threshold           = p.get< float          >("Integral_fraction_threshold");
}


//-----------------------------------------------------------------------------
void cluster::CornerFinderAlg::TakeInRaw( art::Event const&evt)

{
  // Use the TFile service in art
  art::ServiceHandle<art::TFileService> tfs;
  
  /* Get the list of wires.*/
  art::PtrVector<recob::Wire> WireObj;
  
  art::Handle< std::vector<recob::Wire> > wireVecHandle;
  try{

    evt.getByLabel(fCalDataModuleLabel,wireVecHandle);
    
    for(size_t wireIter = 0; wireIter < wireVecHandle->size(); wireIter++){
      art::Ptr<recob::Wire> wire(wireVecHandle, wireIter);
      WireObj.push_back(wire);
      
    }//<---End Looping over wires
    
  }//<---End try

  // Catch if things go badly
  catch(cet::exception& excep){
    std::cout<<"Bail out!"<<std::endl;;
  }
  
  // Getting the bins for the histograms in terms of number of wires and number of time ticks
  art::ServiceHandle<geo::Geometry> geom;
  const unsigned int nPlanes = geom->Nplanes();
  const unsigned int nTimeTicks = (WireObj.at(0))->NSignal();

  
  // Creating the histograms
  WireData_histos = new TH2F*[nPlanes];
  for (uint i_plane=0; i_plane < nPlanes; i_plane++){
    
    std::stringstream ss_tmp_name,ss_tmp_title;
    ss_tmp_name << "h_WireData_" << i_plane;
    ss_tmp_title << fCalDataModuleLabel << " wire data for plane " << i_plane << ";Wire Number;Time Tick";
    WireData_histos[i_plane] = tfs->make<TH2F>(ss_tmp_name.str().c_str(),ss_tmp_title.str().c_str(),
					       geom->Nwires(i_plane),0,geom->Nwires(i_plane),nTimeTicks,0,nTimeTicks);

  }
  
  
  /* For now, we need something to associate each wire in the histogram with a wire_id.
     This is not a beautiful way of handling this, but for now it should work. */
  WireData_IDs = new geo::WireID*[nPlanes];
  for(uint i_plane=0; i_plane < nPlanes; i_plane++)
    WireData_IDs[i_plane] = new geo::WireID[geom->Nwires(i_plane)];
  
  /* Now do the loop over the wires. */
  for (auto const& wire : WireObj) {
    
    
    std::vector<geo::WireID> possible_wireIDs = geom->ChannelToWire(wire->Channel());
    geo::WireID this_wireID;
    try { this_wireID = possible_wireIDs.at(0);}
    catch(cet::exception& excep) { std::cout << "Bail out! No Possible Wires!"<< std::endl; }
    
    uint i_plane = this_wireID.Plane;
    uint i_wire = this_wireID.Wire;

    WireData_IDs[i_plane][i_wire] = this_wireID;
    
    std::vector<float> signal(wire->Signal());
    for(uint i_time = 0; i_time < nTimeTicks; i_time++){
      WireData_histos[i_plane]->SetBinContent(i_wire,i_time,signal[i_time]);  
    }//<---End time loop
    
  }//<-- End loop over wires
  
  
}//<---End TakeInRaw



/// All other methods go below here.....................

//-----------------------------------------------------------------------------------
// This gives us a vecotr of EndPoint2D objects that correspond to possible corners
std::vector<recob::EndPoint2D> cluster::CornerFinderAlg::get_feature_points(){

  std::vector<recob::EndPoint2D> corner_vector;

  art::ServiceHandle<geo::Geometry> geom;
  const unsigned int nPlanes = geom->Nplanes();

  for(uint i=0; i < nPlanes; i++){
    geo::PlaneGeo pg = geom->Plane(nPlanes);
    attach_feature_points(WireData_histos[i],WireData_IDs[i],pg.View(),corner_vector);
  }

  return corner_vector;
}

//-----------------------------------------------------------------------------------
// This gives us a vecotr of EndPoint2D objects that correspond to possible corners
// Uses line integral score as corner strength
std::vector<recob::EndPoint2D> cluster::CornerFinderAlg::get_feature_points_LineIntegralScore(){

  std::vector<recob::EndPoint2D> corner_vector;

  art::ServiceHandle<geo::Geometry> geom;
  const unsigned int nPlanes = geom->Nplanes();

  for(unsigned int i=0; i < nPlanes; ++i){
    geo::PlaneGeo pg = geom->Plane(i);
    attach_feature_points_LineIntegralScore(WireData_histos[i],WireData_IDs[i],pg.View(),corner_vector);
  }

  return corner_vector;
}

//-----------------------------------------------------------------------------
// This makes our data histogram from a collection of wires
//TH2F* cluster::CornerFinderAlg::make_data_histogram(std::vector<recob::Wire> wires, geo::View_t view){
//}


//-----------------------------------------------------------------------------
// This puts on all the feature points in a given view, using a given data histogram
void cluster::CornerFinderAlg::attach_feature_points(TH2F *h_wire_data, geo::WireID *wireIDs, geo::View_t view, 
						     std::vector<recob::EndPoint2D> & corner_vector){

  // Use the TFile service in art
  art::ServiceHandle<art::TFileService> tfs;

  const int x_bins = h_wire_data->GetNbinsX();
  const float x_min = h_wire_data->GetXaxis()->GetBinLowEdge(1);
  const float x_max = h_wire_data->GetXaxis()->GetBinUpEdge(x_bins);

  const int y_bins = h_wire_data->GetNbinsY();
  const float y_min = h_wire_data->GetYaxis()->GetBinLowEdge(1);
  const float y_max = h_wire_data->GetYaxis()->GetBinUpEdge(y_bins);

  const int converted_y_bins = y_bins/fConversion_bins_per_input_y;
  const int converted_x_bins = x_bins/fConversion_bins_per_input_x;

  std::stringstream conversion_name; conversion_name << "h_conversion_" << view;
  std::stringstream dx_name; dx_name << "h_derivative_x_" << view;
  std::stringstream dy_name; dy_name << "h_derivative_x_" << view;
  std::stringstream cornerScore_name; cornerScore_name << "h_cornerScore_" << view;
  std::stringstream maxSuppress_name; maxSuppress_name << "h_maxSuppress_" << view;

  TH2F *h_conversion = tfs->make<TH2F>((conversion_name.str()).c_str(),"Image Conversion Histogram",
				       converted_x_bins,x_min,x_max,
				       converted_y_bins,y_min,y_max);
  TH2F *h_derivative_x = tfs->make<TH2F>((dx_name.str()).c_str(),"Partial Derivatives (x)",
					 converted_x_bins,x_min,x_max,
					 converted_y_bins,y_min,y_max);
  TH2F *h_derivative_y = tfs->make<TH2F>((dy_name.str()).c_str(),"Partial Derivatives (y)",
					 converted_x_bins,x_min,x_max,
					 converted_y_bins,y_min,y_max);
  TH2D *h_cornerScore = tfs->make<TH2D>((cornerScore_name.str()).c_str(),"Feature Point Corner Score",
					converted_x_bins,x_min,x_max,
					converted_y_bins,y_min,y_max);
  TH2D *h_maxSuppress = tfs->make<TH2D>((maxSuppress_name.str()).c_str(),"Corner Points (Maximum Suppressed)",
					converted_x_bins,x_min,x_max,
					converted_y_bins,y_min,y_max);

  create_image_histo(h_wire_data,h_conversion);  
  create_derivative_histograms(h_conversion,h_derivative_x,h_derivative_y);
  create_cornerScore_histogram(h_derivative_x,h_derivative_y,h_cornerScore);
  perform_maximum_suppression(h_cornerScore,corner_vector,wireIDs,view,h_maxSuppress);

}

//-----------------------------------------------------------------------------
// This puts on all the feature points in a given view, using a given data histogram
void cluster::CornerFinderAlg::attach_feature_points_LineIntegralScore(TH2F *h_wire_data, geo::WireID *wireIDs, geo::View_t view, 
								       std::vector<recob::EndPoint2D> & corner_vector){

  // Use the TFile service in art
  art::ServiceHandle<art::TFileService> tfs;

  const int x_bins = h_wire_data->GetNbinsX();
  const float x_min = h_wire_data->GetXaxis()->GetBinLowEdge(1);
  const float x_max = h_wire_data->GetXaxis()->GetBinUpEdge(x_bins);

  const int y_bins = h_wire_data->GetNbinsY();
  const float y_min = h_wire_data->GetYaxis()->GetBinLowEdge(1);
  const float y_max = h_wire_data->GetYaxis()->GetBinUpEdge(y_bins);

  const int converted_y_bins = y_bins/fConversion_bins_per_input_y;
  const int converted_x_bins = x_bins/fConversion_bins_per_input_x;

  std::stringstream conversion_name; conversion_name << "h_conversion_" << view;
  std::stringstream dx_name; dx_name << "h_derivative_x_" << view;
  std::stringstream dy_name; dy_name << "h_derivative_x_" << view;
  std::stringstream cornerScore_name; cornerScore_name << "h_cornerScore_" << view;
  std::stringstream maxSuppress_name; maxSuppress_name << "h_maxSuppress_" << view;

  TH2F *h_conversion = tfs->make<TH2F>((conversion_name.str()).c_str(),"Image Conversion Histogram",
				       converted_x_bins,x_min,x_max,
				       converted_y_bins,y_min,y_max);
  TH2F *h_derivative_x = tfs->make<TH2F>((dx_name.str()).c_str(),"Partial Derivatives (x)",
					 converted_x_bins,x_min,x_max,
					 converted_y_bins,y_min,y_max);
  TH2F *h_derivative_y = tfs->make<TH2F>((dy_name.str()).c_str(),"Partial Derivatives (y)",
					 converted_x_bins,x_min,x_max,
					 converted_y_bins,y_min,y_max);
  TH2D *h_cornerScore = tfs->make<TH2D>((cornerScore_name.str()).c_str(),"Feature Point Corner Score",
					converted_x_bins,x_min,x_max,
					converted_y_bins,y_min,y_max);
  TH2D *h_maxSuppress = tfs->make<TH2D>((maxSuppress_name.str()).c_str(),"Corner Points (Maximum Suppressed)",
					converted_x_bins,x_min,x_max,
					converted_y_bins,y_min,y_max);

  create_image_histo(h_wire_data,h_conversion);  
  create_derivative_histograms(h_conversion,h_derivative_x,h_derivative_y);
  create_cornerScore_histogram(h_derivative_x,h_derivative_y,h_cornerScore);

  std::vector<recob::EndPoint2D> corner_vector_tmp;
  perform_maximum_suppression(h_cornerScore,corner_vector_tmp,wireIDs,view,h_maxSuppress);

  std::stringstream LI_name; conversion_name << "h_lineIntegralScore_" << view;
  TH2F *h_lineIntegralScore =  tfs->make<TH2F>((LI_name.str()).c_str(),"Line Integral Score",
					       x_bins,x_min,x_max,
					       y_bins,y_min,y_max);
  calculate_line_integral_score(h_wire_data,corner_vector_tmp,corner_vector,h_lineIntegralScore);
    
}


//-----------------------------------------------------------------------------
// Convert to pixel
void cluster::CornerFinderAlg::create_image_histo(TH2F *h_wire_data, TH2F *h_conversion) {
  
  int counter_ybin=1;
  int counter_xbin=1;
  float temp_sum=0;

  for(int iy=1; iy<=h_conversion->GetNbinsY(); iy++){
    for(int ix=1; ix<=h_conversion->GetNbinsX(); ix++){

      temp_sum=0;
      for(int jy=counter_ybin; jy<counter_ybin+fConversion_bins_per_input_y; jy++){
	for(int jx=counter_xbin; jx<counter_xbin+fConversion_bins_per_input_x; jx++){
	  temp_sum += h_wire_data->GetBinContent(jx,jy);
	}
      }
      
      counter_ybin += fConversion_bins_per_input_y;
      counter_xbin += fConversion_bins_per_input_x;


      if( temp_sum > fConversion_threshold)
	h_conversion->SetBinContent(ix,iy,temp_sum);
    }
  }


  return;
}

//-----------------------------------------------------------------------------
// Derivative

void cluster::CornerFinderAlg::create_derivative_histograms(TH2F *h_conversion, TH2F *h_derivative_x, TH2F *h_derivative_y){

  const int x_bins = h_conversion->GetNbinsX();
  const int y_bins = h_conversion->GetNbinsY();

  for(int iy=1; iy<=y_bins; iy++){
    for(int ix=1; ix<=x_bins; ix++){
      
      //set to zero the derivatives near the corners...
      if ( iy <= fDerivative_neighborhood 
	   || ix <= fDerivative_neighborhood 
	   || iy >= (y_bins-fDerivative_neighborhood) 
	   || ix == (x_bins-fDerivative_neighborhood) )
	{
	  h_derivative_x->SetBinContent(ix,iy,0);
	  h_derivative_y->SetBinContent(ix,iy,0);
	}
      
      else
	{
	  h_derivative_x->SetBinContent(ix,iy,
					2*(h_conversion->GetBinContent(ix+1,iy)-h_conversion->GetBinContent(ix-1,iy))
					+ 1*(h_conversion->GetBinContent(ix+1,iy+1)-h_conversion->GetBinContent(ix-1,iy+1))
					+ 1*(h_conversion->GetBinContent(ix+1,iy-1)-h_conversion->GetBinContent(ix-1,iy-1)));
	  h_derivative_y->SetBinContent(ix,iy,
					2*(h_conversion->GetBinContent(ix,iy+1)-h_conversion->GetBinContent(ix,iy-1))
					+ 1*(h_conversion->GetBinContent(ix-1,iy+1)-h_conversion->GetBinContent(ix-1,iy-1))
					+ 1*(h_conversion->GetBinContent(ix+1,iy+1)-h_conversion->GetBinContent(ix+1,iy-1)));
	}
    }
  }

  return;

}


//-----------------------------------------------------------------------------
// Corner Score

void cluster::CornerFinderAlg::create_cornerScore_histogram(TH2F *h_derivative_x, TH2F *h_derivative_y, TH2D *h_cornerScore){

  const int x_bins = h_derivative_x->GetNbinsX();
  const int y_bins = h_derivative_y->GetNbinsY();

  const double noble_epsilon = 1e-5;
  //const double harris_kappa = 0.05;
  
  //the structure tensor elements
  float st_xx, st_xy, st_yy;

  for(int iy=1; iy<=y_bins; iy++){
    for(int ix=1; ix<=x_bins; ix++){
      
      //ignore things that will push us over the histogram edge
      if( ix-fCornerScore_neighborhood <= 0 
	  || ix+fCornerScore_neighborhood >= x_bins
	  || iy-fCornerScore_neighborhood <=0
	  || iy+fCornerScore_neighborhood>=y_bins)
	{
	  continue;
	}
      

      else{

	st_xx=0.; st_xy=0.; st_yy=0.;

	for(int jx=ix-fCornerScore_neighborhood; jx<=ix+fCornerScore_neighborhood; jx++){
	  for(int jy=iy-fCornerScore_neighborhood; jy<=iy+fCornerScore_neighborhood; jy++){

	    st_xx += h_derivative_x->GetBinContent(jx,jy)*h_derivative_x->GetBinContent(jx,jy);
	    st_yy += h_derivative_y->GetBinContent(jx,jy)*h_derivative_y->GetBinContent(jx,jy);
	    st_xy += h_derivative_x->GetBinContent(jx,jy)*h_derivative_y->GetBinContent(jx,jy);

	  }
	}

	
	h_cornerScore->SetBinContent(ix,iy,
				     (st_xx*st_yy-st_xy*st_xy) / (st_xx+st_yy + noble_epsilon));
      }

    }
  }

  return;
}


//-----------------------------------------------------------------------------
// Max Supress
size_t cluster::CornerFinderAlg::perform_maximum_suppression(TH2D *h_cornerScore, 
							     std::vector<recob::EndPoint2D> & corner_vector,
							     geo::WireID *wireIDs, 
							     geo::View_t view, 
							     TH2D *h_maxSuppress=NULL){

  const int x_bins = h_cornerScore->GetNbinsX();
  const int y_bins = h_cornerScore->GetNbinsY();

  double temp_max;
  bool temp_center_bin;

  for(int iy=1; iy<=y_bins; iy++){
    for(int ix=1; ix<=x_bins; ix++){

      if(h_cornerScore->GetBinContent(ix,iy) < fMaxSuppress_threshold)
	continue;

      temp_max = -1000;
      temp_center_bin = false;

      for(int jx=ix-fMaxSuppress_neighborhood; jx<=ix+fMaxSuppress_neighborhood; jx++){
	for(int jy=iy-fMaxSuppress_neighborhood; jy<=iy+fMaxSuppress_neighborhood; jy++){

	  if(h_cornerScore->GetBinContent(jx,jy) > temp_max){
	    temp_max = h_cornerScore->GetBinContent(jx,jy);
	    if(jx==ix && jy==iy) temp_center_bin=true;
	    else{ temp_center_bin=false; }
	  }

	}
      }

      if(temp_center_bin){
	
	float time_tick = fConversion_bins_per_input_y*(0.5 * (float)iy);
	int wire_number = fConversion_bins_per_input_x*(0.5 * ix);
	double totalQ = 0;
	int id = 0;
	recob::EndPoint2D corner(time_tick,
				 wireIDs[wire_number],
				 h_cornerScore->GetBinContent(ix,iy),
				 id,
				 view,
				 totalQ);
	corner_vector.push_back(corner);

	if(h_maxSuppress)
	  h_maxSuppress->SetBinContent(ix,iy,h_cornerScore->GetBinContent(ix,iy));	
      }
      
    }
  }
  
  return corner_vector.size();

}

/* Silly little function for doing a line integral type thing. Needs improvement. */
float cluster::CornerFinderAlg::line_integral(TH2F *hist, int begin_x, float begin_y, int end_x, float end_y, float threshold){

  int x1 = hist->GetXaxis()->FindBin( begin_x );
  int y1 = hist->GetYaxis()->FindBin( begin_y );
  int x2 = hist->GetXaxis()->FindBin( end_x );
  int y2 = hist->GetYaxis()->FindBin( end_y );

  if(x1==x2 && abs(y1-y2)<1e-5)
    return 0;

  if(x2<x1){
    int tmp = x2;
    x2 = x1;
    x1 = tmp;

    int tmp_y = y2;
    y2 = y1;
    y1 = tmp_y;
  }

  float fraction = 0;
  int bin_counter = 0;

  if(x2!=x1){

    float slope = (y2-y1)/((float)(x2-x1));

    for(int ix=x1; ix<=x2; ix++){

      int y_min,y_max;
      
      if(slope>=0){
	y_min = y1 + slope*(ix-x1);
	y_max = y1 + slope*(ix+1-x1);
      }
      else {
	y_max = (y1+1) + slope*(ix-x1);
	y_min = (y1+1) + slope*(ix+1-x1);
      }

      for(int iy=y_min; iy<=y_max; iy++){
	bin_counter++;

	if( hist->GetBinContent(ix,iy) > threshold )
	  fraction += 1.;
      }

    }
  }
  else{
    
    int y_min,y_max;
    if(y1<y2){
      y_min=y1; y_max=y2;
    }
    else{
      y_min=y2; y_max=y1;
    }
    for(int iy=y_min; iy<=y_max; iy++){
	bin_counter++;
	if( hist->GetBinContent(x1,iy) > threshold)
	  fraction += 1.;
      }

  }

  return fraction/bin_counter;
}

//-----------------------------------------------------------------------------
// Do the silly little line integral score thing
size_t cluster::CornerFinderAlg::calculate_line_integral_score( TH2F* h_wire_data, 
								std::vector<recob::EndPoint2D> const & corner_vector, 
								std::vector<recob::EndPoint2D> & corner_lineIntegralScore_vector,
								TH2F* h_lineIntegralScore){

  float score;

  for(auto const i_corner : corner_vector){

    score=0;
    
    for(auto const j_corner : corner_vector){


      if( line_integral(h_wire_data,
			i_corner.WireID().Wire,i_corner.DriftTime(),
			j_corner.WireID().Wire,j_corner.DriftTime(),
			fIntegral_bin_threshold) > fIntegral_fraction_threshold)
	{
	  score+=1.;
	}

    }

    recob::EndPoint2D corner(i_corner.DriftTime(),
			     i_corner.WireID(),
			     score,
			     i_corner.ID(),
			     i_corner.View(),
			     i_corner.Charge());

    corner_lineIntegralScore_vector.push_back(corner);
    

    h_lineIntegralScore->SetBinContent(h_wire_data->GetXaxis()->FindBin(i_corner.WireID().Wire),
				       h_wire_data->GetYaxis()->FindBin(i_corner.DriftTime()),
				       score);
    
  }
  
  return corner_lineIntegralScore_vector.size();
}



