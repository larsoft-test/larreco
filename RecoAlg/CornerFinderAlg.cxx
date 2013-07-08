////////////////////////////////////////////////////////////////////////
//
// CornerFinderAlg class
//
// wketchum@fnal.gov
//
// CornerFinder is meant to use image-processing techniques (mainly Harris-Stephens
// corner-finding) to find "corners" using the information from calibrated wires.
//  
//  Conversion_algorithm options:  
//     standard --- basically a copy of the calibrated wires
//     skeleton --- a thinned copy of the calibrated wires
//     binary   --- ticks above threshold get assigned a value 10*threshold, everything else = threshold
//     function --- apply a function (like a double-Gaussian) to a neighborhood around each tick
//       
//  Derivative options:
//     Sobel --- apply a Sobel mask (neighborhood of 1 or 2 supported)
//     local --- take slope from immediately neighboring bins (neighborhood of 1 supported)
//
//  CornerScore_algorithm options:
//     Noble  --- determinant / (trace + Noble_epsilon)
//     Harris --- determinant - (trace)^2 * Harris_kappa
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

#include "RawData/raw.h"
#include "RawData/RawDigit.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Wire.h"
#include "RecoAlg/CornerFinderAlg.h"

#include "RecoObjects/BezierTrack.h"

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
{
  this->reconfigure(pset);

  // set the sizes of the WireData_histos and WireData_IDs
  unsigned int nPlanes = fGeom->Nplanes();
  WireData_histos.resize(nPlanes);
  fConversion_histos.resize(nPlanes);
  fDerivativeX_histos.resize(nPlanes);
  fDerivativeY_histos.resize(nPlanes);
  fCornerScore_histos.resize(nPlanes);
  fMaxSuppress_histos.resize(nPlanes);

  /* For now, we need something to associate each wire in the histogram with a wire_id.
     This is not a beautiful way of handling this, but for now it should work. */
  WireData_IDs.resize(nPlanes);
  for(uint i_plane=0; i_plane < nPlanes; ++i_plane)
    WireData_IDs[i_plane].resize(fGeom->Nwires(i_plane));
  

}

//-----------------------------------------------------------------------------
cluster::CornerFinderAlg::~CornerFinderAlg()
{
  for (auto wd_histo : WireData_histos)
    delete wd_histo;

  for (auto histo : fConversion_histos) delete histo;
  for (auto histo : fDerivativeX_histos) delete histo;
  for (auto histo : fDerivativeY_histos) delete histo;
  for (auto histo : fCornerScore_histos) delete histo;
  for (auto histo : fMaxSuppress_histos) delete histo;

  WireData_histos.clear();
  WireData_IDs.clear();
}

//-----------------------------------------------------------------------------
void cluster::CornerFinderAlg::reconfigure(fhicl::ParameterSet const& p)
{
  // ### These are all the tuneable .fcl file parameters from the event ###
  fCalDataModuleLabel  			 = p.get< std::string 	 >("CalDataModuleLabel");
  fConversion_threshold     		 = p.get< float    	 >("Conversion_threshold");
  fConversion_bins_per_input_x  	 = p.get< int      	 >("Conversion_bins_per_input_x");
  fConversion_bins_per_input_y       	 = p.get< int      	 >("Conversion_bins_per_input_y");
  fConversion_algorithm                  = p.get< std::string    >("Conversion_algorithm");
  fConversion_func                       = p.get< std::string    >("Conversion_function");
  fConversion_func_neighborhood     	 = p.get< int		 >("Conversion_func_neighborhood");
  fDerivative_method        		 = p.get< std::string    >("Derivative_method");
  fDerivative_neighborhood     	         = p.get< int		 >("Derivative_neighborhood");
  fCornerScore_neighborhood     	 = p.get< int		 >("CornerScore_neighborhood");
  fCornerScore_algorithm		 = p.get< std::string    >("CornerScore_algorithm");
  fCornerScore_Noble_epsilon		 = p.get< float          >("CornerScore_Noble_epsilon");
  fCornerScore_Harris_kappa		 = p.get< float          >("CornerScore_Harris_kappa");
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
  const unsigned int nTimeTicks = (WireObj.at(0))->NSignal();

  event_number = evt.event();
  run_number = evt.run();

  // Creating the histograms
  for (uint i_plane=0; i_plane < fGeom->Nplanes(); i_plane++){
    
    std::stringstream ss_tmp_name,ss_tmp_title;
    ss_tmp_name << "h_WireData_" << i_plane << "_" << run_number << "_" << event_number;
    ss_tmp_title << fCalDataModuleLabel << " wire data for plane " << i_plane << ", Run " << run_number << ", Event " << event_number << ";Wire Number;Time Tick";

    if(WireData_histos[i_plane]) {
      WireData_histos[i_plane]->Reset();
      WireData_histos[i_plane]->SetName(ss_tmp_name.str().c_str());
      WireData_histos[i_plane]->SetTitle(ss_tmp_title.str().c_str());
    }
    else
      WireData_histos[i_plane] = new TH2F(ss_tmp_name.str().c_str(),
					  ss_tmp_title.str().c_str(),
					  fGeom->Nwires(i_plane),
					  0,
					  fGeom->Nwires(i_plane),
					  nTimeTicks,
					  0,
					  nTimeTicks);
    /*
      WireData_histos[i_plane] = tfs->make<TH2F>(ss_tmp_name.str().c_str(),
                                                 ss_tmp_title.str().c_str(),
                                                 fGeom->Nwires(i_plane),
                                                 0,
                                                 fGeom->Nwires(i_plane),
                                                 nTimeTicks,
                                                 0,
                                                 nTimeTicks);
    */

  }
  
  /* Now do the loop over the wires. */
  for (auto const& wire : WireObj) {
    
    
    std::vector<geo::WireID> possible_wireIDs = fGeom->ChannelToWire(wire->Channel());
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
void cluster::CornerFinderAlg::get_feature_points(std::vector<recob::EndPoint2D> & corner_vector){

  for(unsigned int cstat = 0; cstat < fGeom->Ncryostats(); ++cstat){
    for(unsigned int tpc = 0; tpc < fGeom->Cryostat(cstat).NTPC(); ++tpc){
      for(unsigned int plane = 0; plane < fGeom->Cryostat(cstat).TPC(tpc).Nplanes(); ++plane){
	attach_feature_points(WireData_histos[plane],WireData_IDs[plane],fGeom->Cryostat(cstat).TPC(tpc).Plane(plane).View(),corner_vector);
      }
    }
  }

}

//-----------------------------------------------------------------------------------
// This gives us a vecotr of EndPoint2D objects that correspond to possible corners
// Uses line integral score as corner strength
void cluster::CornerFinderAlg::get_feature_points_LineIntegralScore(std::vector<recob::EndPoint2D> & corner_vector){

  for(unsigned int cstat = 0; cstat < fGeom->Ncryostats(); ++cstat){
    for(unsigned int tpc = 0; tpc < fGeom->Cryostat(cstat).NTPC(); ++tpc){
      for(unsigned int plane = 0; plane < fGeom->Cryostat(cstat).TPC(tpc).Nplanes(); ++plane){
	attach_feature_points_LineIntegralScore(WireData_histos[plane],WireData_IDs[plane],fGeom->Cryostat(cstat).TPC(tpc).Plane(plane).View(),corner_vector);
      }
    }
  }

}

//-----------------------------------------------------------------------------
// This puts on all the feature points in a given view, using a given data histogram
void cluster::CornerFinderAlg::attach_feature_points(TH2F *h_wire_data, 
						     std::vector<geo::WireID> wireIDs, 
						     geo::View_t view, 
						     std::vector<recob::EndPoint2D> & corner_vector){


  const int x_bins = h_wire_data->GetNbinsX();
  const float x_min = h_wire_data->GetXaxis()->GetBinLowEdge(1);
  const float x_max = h_wire_data->GetXaxis()->GetBinUpEdge(x_bins);

  const int y_bins = h_wire_data->GetNbinsY();
  const float y_min = h_wire_data->GetYaxis()->GetBinLowEdge(1);
  const float y_max = h_wire_data->GetYaxis()->GetBinUpEdge(y_bins);

  const int converted_y_bins = y_bins/fConversion_bins_per_input_y;
  const int converted_x_bins = x_bins/fConversion_bins_per_input_x;

  std::stringstream conversion_name;  conversion_name  << "h_conversion_"   << view << "_" << run_number << "_" << event_number;
  std::stringstream dx_name;          dx_name          << "h_derivative_x_" << view << "_" << run_number << "_" << event_number;
  std::stringstream dy_name;          dy_name          << "h_derivative_y_" << view << "_" << run_number << "_" << event_number;
  std::stringstream cornerScore_name; cornerScore_name << "h_cornerScore_"  << view << "_" << run_number << "_" << event_number;
  std::stringstream maxSuppress_name; maxSuppress_name << "h_maxSuppress_"  << view << "_" << run_number << "_" << event_number;

  if(fConversion_histos[view]){
    fConversion_histos[view]->Reset();
    fConversion_histos[view]->SetName(conversion_name.str().c_str());
  }
  else
    fConversion_histos[view] = new TH2F(conversion_name.str().c_str(),"Image Conversion Histogram",
					   converted_x_bins,x_min,x_max,
					   converted_y_bins,y_min,y_max);
  
  if(fDerivativeX_histos[view]){
    fDerivativeX_histos[view]->Reset();
    fDerivativeX_histos[view]->SetName(dx_name.str().c_str());
  }
  else
    fDerivativeX_histos[view] = new TH2F(dx_name.str().c_str(),"Partial Derivatives (x)",
					    converted_x_bins,x_min,x_max,
					    converted_y_bins,y_min,y_max);
  
  if(fDerivativeY_histos[view]){
    fDerivativeY_histos[view]->Reset();
    fDerivativeY_histos[view]->SetName(dy_name.str().c_str());
  }
  else
    fDerivativeY_histos[view] = new TH2F(dy_name.str().c_str(),"Partial Derivatives (y)",
					    converted_x_bins,x_min,x_max,
					    converted_y_bins,y_min,y_max);
  
  if(fCornerScore_histos[view]){
    fCornerScore_histos[view]->Reset();
    fCornerScore_histos[view]->SetName(cornerScore_name.str().c_str());
  }
  else
    fCornerScore_histos[view] = new TH2D(cornerScore_name.str().c_str(),"Corner Score",
					    converted_x_bins,x_min,x_max,
					    converted_y_bins,y_min,y_max);

  if(fMaxSuppress_histos[view]){
    fMaxSuppress_histos[view]->Reset();
    fMaxSuppress_histos[view]->SetName(maxSuppress_name.str().c_str());
  }
  else
    fMaxSuppress_histos[view] = new TH2D(maxSuppress_name.str().c_str(),"Corner Points (Maximum Suppressed)",
					    converted_x_bins,x_min,x_max,
					    converted_y_bins,y_min,y_max);
  
  create_image_histo(h_wire_data,fConversion_histos[view]);  
  create_derivative_histograms(fConversion_histos[view],fDerivativeX_histos[view],fDerivativeY_histos[view]);
  create_cornerScore_histogram(fDerivativeX_histos[view],fDerivativeY_histos[view],fCornerScore_histos[view]);
  perform_maximum_suppression(fCornerScore_histos[view],corner_vector,wireIDs,view,fMaxSuppress_histos[view]);
}

//-----------------------------------------------------------------------------
// This puts on all the feature points in a given view, using a given data histogram
void cluster::CornerFinderAlg::attach_feature_points_LineIntegralScore(TH2F *h_wire_data, 
								       std::vector<geo::WireID> wireIDs, 
								       geo::View_t view, 
								       std::vector<recob::EndPoint2D> & corner_vector){


  const int   x_bins = h_wire_data->GetNbinsX();
  const float x_min  = h_wire_data->GetXaxis()->GetBinLowEdge(1);
  const float x_max  = h_wire_data->GetXaxis()->GetBinUpEdge(x_bins);

  const int   y_bins = h_wire_data->GetNbinsY();
  const float y_min  = h_wire_data->GetYaxis()->GetBinLowEdge(1);
  const float y_max  = h_wire_data->GetYaxis()->GetBinUpEdge(y_bins);

  const int converted_y_bins = y_bins/fConversion_bins_per_input_y;
  const int converted_x_bins = x_bins/fConversion_bins_per_input_x;

  std::stringstream conversion_name;  conversion_name  << "h_conversion_"   << view << "_" << run_number << "_" << event_number;
  std::stringstream dx_name;          dx_name          << "h_derivative_x_" << view << "_" << run_number << "_" << event_number;
  std::stringstream dy_name;          dy_name          << "h_derivative_y_" << view << "_" << run_number << "_" << event_number;
  std::stringstream cornerScore_name; cornerScore_name << "h_cornerScore_"  << view << "_" << run_number << "_" << event_number;
  std::stringstream maxSuppress_name; maxSuppress_name << "h_maxSuppress_"  << view << "_" << run_number << "_" << event_number;

  TH2F h_conversion  ((conversion_name.str()).c_str(),
		      "Image Conversion Histogram",
		      converted_x_bins,x_min,x_max,
		      converted_y_bins,y_min,y_max);
  TH2F h_derivative_x((dx_name.str()).c_str(),
		      "Partial Derivatives (x)",
		      converted_x_bins,x_min,x_max,
		      converted_y_bins,y_min,y_max);
  TH2F h_derivative_y((dy_name.str()).c_str(),
		      "Partial Derivatives (y)",
		      converted_x_bins,x_min,x_max,
		      converted_y_bins,y_min,y_max);
  TH2D h_cornerScore ((cornerScore_name.str()).c_str(),
		      "Feature Point Corner Score",
		      converted_x_bins,x_min,x_max,
		      converted_y_bins,y_min,y_max);
  TH2D h_maxSuppress ((maxSuppress_name.str()).c_str(),
		      "Corner Points (Maximum Suppressed)",
		      converted_x_bins,x_min,x_max,
		      converted_y_bins,y_min,y_max);
  
  create_image_histo(h_wire_data,&h_conversion);  
  create_derivative_histograms(&h_conversion,&h_derivative_x,&h_derivative_y);
  create_cornerScore_histogram(&h_derivative_x,&h_derivative_y,&h_cornerScore);

  std::vector<recob::EndPoint2D> corner_vector_tmp;
  perform_maximum_suppression(&h_cornerScore,corner_vector_tmp,wireIDs,view,&h_maxSuppress);

  std::stringstream LI_name; LI_name << "h_lineIntegralScore_" << view << "_" << run_number << "_" << event_number;
  TH2F h_lineIntegralScore((LI_name.str()).c_str(),
			   "Line Integral Score",
			   x_bins,x_min,x_max,
			   y_bins,y_min,y_max);
  calculate_line_integral_score(h_wire_data,corner_vector_tmp,corner_vector,&h_lineIntegralScore);
    
}


//-----------------------------------------------------------------------------
// Convert to pixel
void cluster::CornerFinderAlg::create_image_histo(TH2F *h_wire_data, TH2F *h_conversion) {
  
  double temp_integral=0;

  const TF2 fConversion_TF2("fConversion_func",fConversion_func.c_str(),-20,20,-20,20);

  for(int ix=1; ix<=h_conversion->GetNbinsX(); ix++){
    for(int iy=1; iy<=h_conversion->GetNbinsY(); iy++){
      
      temp_integral = h_wire_data->Integral(ix,ix,iy,iy);

      if( temp_integral > fConversion_threshold){

	if(fConversion_algorithm.compare("binary")==0)
	  h_conversion->SetBinContent(ix,iy,10*fConversion_threshold);
	else if(fConversion_algorithm.compare("standard")==0)
	  h_conversion->SetBinContent(ix,iy,temp_integral);

	else if(fConversion_algorithm.compare("function")==0){

	  temp_integral = 0;
	  for(int jx=ix-fConversion_func_neighborhood; jx<=ix+fConversion_func_neighborhood; jx++){
	    for(int jy=iy-fConversion_func_neighborhood; jy<=iy+fConversion_func_neighborhood; jy++){
	      temp_integral += h_wire_data->GetBinContent(jx,jy)*fConversion_TF2.Eval(ix-jx,iy-jy);
	    }
	  }
	  h_conversion->SetBinContent(ix,iy,temp_integral);
	}

	else if(fConversion_algorithm.compare("skeleton")==0){
	  
	  if( (temp_integral > h_wire_data->GetBinContent(ix-1,iy) && temp_integral > h_wire_data->GetBinContent(ix+1,iy))
	      || (temp_integral > h_wire_data->GetBinContent(ix,iy-1) && temp_integral > h_wire_data->GetBinContent(ix,iy+1)))
	    h_conversion->SetBinContent(ix,iy,temp_integral);
	  else 
	    h_conversion->SetBinContent(ix,iy,fConversion_threshold);
	}
	else if(fConversion_algorithm.compare("sk_bin")==0){
	  
	  if( (temp_integral > h_wire_data->GetBinContent(ix-1,iy) && temp_integral > h_wire_data->GetBinContent(ix+1,iy))
	      || (temp_integral > h_wire_data->GetBinContent(ix,iy-1) && temp_integral > h_wire_data->GetBinContent(ix,iy+1)))
	    h_conversion->SetBinContent(ix,iy,10*fConversion_threshold);
	  else 
	    h_conversion->SetBinContent(ix,iy,fConversion_threshold);
	}
	else
	  h_conversion->SetBinContent(ix,iy,temp_integral);
      }

      else
	h_conversion->SetBinContent(ix,iy,fConversion_threshold);
      
    }
  }

}

//-----------------------------------------------------------------------------
// Derivative

void cluster::CornerFinderAlg::create_derivative_histograms(TH2F *h_conversion, TH2F *h_derivative_x, TH2F *h_derivative_y){

  const int x_bins = h_conversion->GetNbinsX();
  const int y_bins = h_conversion->GetNbinsY();

  for(int iy=1+fDerivative_neighborhood; iy<=(y_bins-fDerivative_neighborhood); iy++){
    for(int ix=1+fDerivative_neighborhood; ix<=(x_bins-fDerivative_neighborhood); ix++){
      
      if(fDerivative_method.compare("Sobel")==0){
	
	if(fDerivative_neighborhood==1){
	  h_derivative_x->SetBinContent(ix,iy,
					0.5*(h_conversion->GetBinContent(ix+1,iy)-h_conversion->GetBinContent(ix-1,iy))
					+ 0.25*(h_conversion->GetBinContent(ix+1,iy+1)-h_conversion->GetBinContent(ix-1,iy+1))
					+ 0.25*(h_conversion->GetBinContent(ix+1,iy-1)-h_conversion->GetBinContent(ix-1,iy-1)));
	  h_derivative_y->SetBinContent(ix,iy,
					0.5*(h_conversion->GetBinContent(ix,iy+1)-h_conversion->GetBinContent(ix,iy-1))
					+ 0.25*(h_conversion->GetBinContent(ix-1,iy+1)-h_conversion->GetBinContent(ix-1,iy-1))
					+ 0.25*(h_conversion->GetBinContent(ix+1,iy+1)-h_conversion->GetBinContent(ix+1,iy-1)));
	}
	else if(fDerivative_neighborhood==2){
	  h_derivative_x->SetBinContent(ix,iy,
					12*(h_conversion->GetBinContent(ix+1,iy)-h_conversion->GetBinContent(ix-1,iy))
					+ 8*(h_conversion->GetBinContent(ix+1,iy+1)-h_conversion->GetBinContent(ix-1,iy+1))
					+ 8*(h_conversion->GetBinContent(ix+1,iy-1)-h_conversion->GetBinContent(ix-1,iy-1))
					+ 2*(h_conversion->GetBinContent(ix+1,iy+2)-h_conversion->GetBinContent(ix-1,iy+2))
					+ 2*(h_conversion->GetBinContent(ix+1,iy-2)-h_conversion->GetBinContent(ix-1,iy-2))
					  + 6*(h_conversion->GetBinContent(ix+2,iy)-h_conversion->GetBinContent(ix-2,iy))
					+ 4*(h_conversion->GetBinContent(ix+2,iy+1)-h_conversion->GetBinContent(ix-2,iy+1))
					+ 4*(h_conversion->GetBinContent(ix+2,iy-1)-h_conversion->GetBinContent(ix-2,iy-1))
					+ 1*(h_conversion->GetBinContent(ix+2,iy+2)-h_conversion->GetBinContent(ix-2,iy+2))
					  + 1*(h_conversion->GetBinContent(ix+2,iy-2)-h_conversion->GetBinContent(ix-2,iy-2)));
	  h_derivative_y->SetBinContent(ix,iy,
					12*(h_conversion->GetBinContent(ix,iy+1)-h_conversion->GetBinContent(ix,iy-1))
					+ 8*(h_conversion->GetBinContent(ix-1,iy+1)-h_conversion->GetBinContent(ix-1,iy-1))
					+ 8*(h_conversion->GetBinContent(ix+1,iy+1)-h_conversion->GetBinContent(ix+1,iy-1))
					+ 2*(h_conversion->GetBinContent(ix-2,iy+1)-h_conversion->GetBinContent(ix-2,iy-1))
					+ 2*(h_conversion->GetBinContent(ix+2,iy+1)-h_conversion->GetBinContent(ix+2,iy-1))
					+ 6*(h_conversion->GetBinContent(ix,iy+2)-h_conversion->GetBinContent(ix,iy-2))
					+ 4*(h_conversion->GetBinContent(ix-1,iy+2)-h_conversion->GetBinContent(ix-1,iy-2))
					+ 4*(h_conversion->GetBinContent(ix+1,iy+2)-h_conversion->GetBinContent(ix+1,iy-2))
					+ 1*(h_conversion->GetBinContent(ix-2,iy+2)-h_conversion->GetBinContent(ix-2,iy-2))
					+ 1*(h_conversion->GetBinContent(ix+2,iy+2)-h_conversion->GetBinContent(ix+2,iy-2)));
	}
	else{
	  std::cout << "Sobel derivative not supported for neighborhoods > 2." << std::endl;
	  return;
	}
	
	} //end if Sobel
      
      else if(fDerivative_method.compare("local")==0){
	
	if(fDerivative_neighborhood==1){
	  h_derivative_x->SetBinContent(ix,iy,
					(h_conversion->GetBinContent(ix+1,iy)-h_conversion->GetBinContent(ix-1,iy)));
	  h_derivative_y->SetBinContent(ix,iy,
					(h_conversion->GetBinContent(ix,iy+1)-h_conversion->GetBinContent(ix,iy-1)));
	}
	else{
	  std::cout << "Local derivative not yet supported for neighborhoods > 1." << std::endl;
	  return;
	}
	
	} //end if local

      else{
	std::cout << "Bad derivative algorithm! " << fDerivative_method << std::endl;
	return;
      }

    }
  }


}


//-----------------------------------------------------------------------------
// Corner Score

void cluster::CornerFinderAlg::create_cornerScore_histogram(TH2F *h_derivative_x, TH2F *h_derivative_y, TH2D *h_cornerScore){

  const int x_bins = h_derivative_x->GetNbinsX();
  const int y_bins = h_derivative_y->GetNbinsY();
  
  //the structure tensor elements
  double st_xx, st_xy, st_yy;

  for(int iy=1+fCornerScore_neighborhood; iy<=(y_bins-fCornerScore_neighborhood); iy++){
    for(int ix=1+fCornerScore_neighborhood; ix<=(x_bins-fCornerScore_neighborhood); ix++){

      if(ix==1+fCornerScore_neighborhood){
	st_xx=0.; st_xy=0.; st_yy=0.;

	for(int jx=ix-fCornerScore_neighborhood; jx<=ix+fCornerScore_neighborhood; jx++){
	  for(int jy=iy-fCornerScore_neighborhood; jy<=iy+fCornerScore_neighborhood; jy++){
	    
	    st_xx += h_derivative_x->GetBinContent(jx,jy)*h_derivative_x->GetBinContent(jx,jy);
	    st_yy += h_derivative_y->GetBinContent(jx,jy)*h_derivative_y->GetBinContent(jx,jy);
	    st_xy += h_derivative_x->GetBinContent(jx,jy)*h_derivative_y->GetBinContent(jx,jy);
	    
	  }
	}
      }
      
      // we do it this way to reduce computation time
      else{
	for(int jy=iy-fCornerScore_neighborhood; jy<=iy+fCornerScore_neighborhood; jy++){
	  
	  st_xx -= h_derivative_x->GetBinContent(ix-fCornerScore_neighborhood-1,jy)*h_derivative_x->GetBinContent(ix-fCornerScore_neighborhood-1,jy);
	  st_xx += h_derivative_x->GetBinContent(ix+fCornerScore_neighborhood,jy)*h_derivative_x->GetBinContent(ix+fCornerScore_neighborhood,jy);
	  
	  st_yy -= h_derivative_y->GetBinContent(ix-fCornerScore_neighborhood-1,jy)*h_derivative_y->GetBinContent(ix-fCornerScore_neighborhood-1,jy);
	  st_yy += h_derivative_y->GetBinContent(ix+fCornerScore_neighborhood,jy)*h_derivative_y->GetBinContent(ix+fCornerScore_neighborhood,jy);
	  
	  st_xy -= h_derivative_x->GetBinContent(ix-fCornerScore_neighborhood-1,jy)*h_derivative_y->GetBinContent(ix-fCornerScore_neighborhood-1,jy);
	  st_xy += h_derivative_x->GetBinContent(ix+fCornerScore_neighborhood,jy)*h_derivative_y->GetBinContent(ix+fCornerScore_neighborhood,jy);
	}
      }
      
      if( fCornerScore_algorithm.compare("Noble")==0 ) {
	h_cornerScore->SetBinContent(ix,iy,
				     (st_xx*st_yy-st_xy*st_xy) / (st_xx+st_yy + fCornerScore_Noble_epsilon));
      }
	else if( fCornerScore_algorithm.compare("Harris")==0 ) {
	  h_cornerScore->SetBinContent(ix,iy,
				       (st_xx*st_yy-st_xy*st_xy) - ((st_xx+st_yy)*(st_xx+st_yy)*fCornerScore_Harris_kappa));
	}
	else{
	  std::cout << "BAD CORNER ALGORITHM: " << fCornerScore_algorithm << std::endl;
	  return;
	}
      
    } // end for loop over x bins
  } // end for loop over y bins

}


//-----------------------------------------------------------------------------
// Max Supress
size_t cluster::CornerFinderAlg::perform_maximum_suppression(TH2D *h_cornerScore, 
							     std::vector<recob::EndPoint2D> & corner_vector,
							     std::vector<geo::WireID> wireIDs, 
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
	
	float time_tick = 0.5 * (float)((2*iy-1) * fConversion_bins_per_input_y);
	int wire_number = ( (2*ix-1)*fConversion_bins_per_input_x ) / 2;
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


std::vector<float> cluster::CornerFinderAlg::line_integrals(trkf::BezierTrack& TheTrack, size_t Steps, float threshold){

  std::vector<float> fractions(3);

  for(size_t i=0; i!=Steps; ++i)
    {
      double s = float(i)/float(Steps);
      double uvw[3], ticks[3];
      int c=0, t=0;
      TheTrack.GetProjectedPointUVWT(s, uvw, ticks, c, t);

      for(size_t j=0; j!=WireData_histos.size(); ++j)
        {
          int x = WireData_histos.at(j)->GetXaxis()->FindBin(uvw[j]);
          int y = WireData_histos.at(j)->GetYaxis()->FindBin(ticks[j]);

          if( WireData_histos.at(j)->GetBinContent(x,y) > threshold )
            fractions.at(j) += 1.;
        }
    }
  for(size_t j=0; j!=fractions.size(); ++j)
    fractions.at(j) /= float(Steps);

  return fractions;
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


TH2F* cluster::CornerFinderAlg::GetWireDataHist(unsigned int i_plane){

  if(i_plane >= WireData_histos.size()){
    mf::LogWarning("CornerFinderAlg") << "WARNING:  Requested plane does not exist.";
    return NULL;
  }

  return WireData_histos.at(i_plane);
}

TH2F* cluster::CornerFinderAlg::GetConversionHist(unsigned int i_plane){

  if(i_plane >= fConversion_histos.size()){
    mf::LogWarning("CornerFinderAlg") << "WARNING:  Requested plane does not exist.";
    return NULL;
  }

  return fConversion_histos.at(i_plane);
}

TH2F* cluster::CornerFinderAlg::GetDerivativeXHist(unsigned int i_plane){

  if(i_plane >= fDerivativeX_histos.size()){
    mf::LogWarning("CornerFinderAlg") << "WARNING:  Requested plane does not exist.";
    return NULL;
  }

  return fDerivativeX_histos.at(i_plane);
}

TH2F* cluster::CornerFinderAlg::GetDerivativeYHist(unsigned int i_plane){

  if(i_plane >= fDerivativeY_histos.size()){
    mf::LogWarning("CornerFinderAlg") << "WARNING:  Requested plane does not exist.";
    return NULL;
  }

  return fDerivativeY_histos.at(i_plane);
}

TH2D* cluster::CornerFinderAlg::GetCornerScoreHist(unsigned int i_plane){

  if(i_plane >= fCornerScore_histos.size()){
    mf::LogWarning("CornerFinderAlg") << "WARNING:  Requested plane does not exist.";
    return NULL;
  }

  return fCornerScore_histos.at(i_plane);
}

TH2D* cluster::CornerFinderAlg::GetMaxSuppressHist(unsigned int i_plane){

  if(i_plane >= fMaxSuppress_histos.size()){
    mf::LogWarning("CornerFinderAlg") << "WARNING:  Requested plane does not exist.";
    return NULL;
  }

  return fMaxSuppress_histos.at(i_plane);
}
