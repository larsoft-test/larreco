////////////////////////////////////////////////////////////////////////
/// \file  CornerFinderAlgAlg.h
/// \brief algorithm to find feature 2D points
///
/// \author  
////////////////////////////////////////////////////////////////////////

#ifndef CORNERFINDERALG_H
#define CORNERFINDERALG_H

#include "fhiclcpp/ParameterSet.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "TMath.h"
#include "TH2.h"
#include <vector>
#include <string>
#include "RecoBase/EndPoint2D.h"



namespace cluster { //<---Not sure if this is the right namespace

   class CornerFinderAlg {
   
    public:
    
    explicit CornerFinderAlg(fhicl::ParameterSet const& pset); 
    virtual ~CornerFinderAlg();        

    void   reconfigure(fhicl::ParameterSet const& pset);
    
    /// \ put the output of the algorithm here
    /// \ for now I am putting the example of EndPoint
    /// \ but this is wrong and should be changed
    
//    size_t EndPoint(art::PtrVector<recob::Cluster>                 & clusIn, 
//		    std::vector<recob::EndPoint2D>                 & vtxcol,
//		    std::vector< art::PtrVector<recob::Hit> >      & vtxHitsOut,
//		    art::Event                                const& evt,
//		    std::string                               const& label);



    void TakeInRaw(//art::PtrVector<raw::RawDigit>	        & rawhits,
                   //art::PtrVector<recob::Wire>    & wires,
		   art::Event				const&evt);    
    std::vector<recob::EndPoint2D> get_feature_points();
    std::vector<recob::EndPoint2D> get_feature_points_LineIntegralScore();

    

    private:
    
    
    // Need to list the things we will take in from the .fcl file
    
    // Taking in RawData from the event
    
    std::string  fRawDataModuleLabel;
    std::string  fCalDataModuleLabel;
    
    
    float          fConversion_threshold;
    int            fConversion_bins_per_input_x;
    int            fConversion_bins_per_input_y;
    std::string    fDerivative_method;
    int            fDerivative_neighborhood;
    int            fCornerScore_neighborhood;
    std::string    fCornerScore_algorithm;
    int            fMaxSuppress_neighborhood;
    int            fMaxSuppress_threshold;
    float          fIntegral_bin_threshold;
    float          fIntegral_fraction_threshold;
    
    // Making a vector of histograms
    TH2F* RawData_histos[3];
    
     
    void create_image_histo(TH2F *h_wire_data, TH2F *h_conversion);
    void create_derivative_histograms(TH2F *h_conversion, TH2F *h_derivative_x, TH2F *h_derivative_y);
    void create_cornerScore_histogram(TH2F *h_derivative_x, TH2F *h_derivative_y, TH2D *h_cornerScore);
    size_t perform_maximum_suppression(TH2D *h_cornerScore, 
				   std::vector<recob::EndPoint2D> & corner_vector,
				   geo::View_t view, 
				   TH2D *h_maxSuppress);
				   
    float line_integral(TH2F *hist, int x1, float y1, int x2, float y2, float threshold);				   
    
    size_t calculate_line_integral_score( TH2F* h_wire_data, 
					  std::vector<recob::EndPoint2D> const & corner_vector, 
					  std::vector<recob::EndPoint2D> & corner_lineIntegralScore_vector,
					  TH2F* h_lineIntegralScore);

    void attach_feature_points(TH2F *h_wire_data, 
			       geo::View_t view, 
			       std::vector<recob::EndPoint2D>&);
    void attach_feature_points_LineIntegralScore(TH2F *h_wire_data, 
						 geo::View_t view, 
						 std::vector<recob::EndPoint2D>&);
    

     };//<---End of class CornerFinderAlg


}

#endif //CORNERFINDERALG_H
