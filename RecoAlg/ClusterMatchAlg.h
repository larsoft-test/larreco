////////////////////////////////////////////////////////////////////////
// \file ClusterMatchAlg.h
//
// \brief ClusterMatchAlg header file
//
// \author kazuhiro@nevis.columbia.edu
//
////////////////////////////////////////////////////////////////////////

#ifndef CLUSTERMATCHALG_H
#define CLUSTERMATCHALG_H

// ART includes
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/FindManyP.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"

// LArSoft
#include "SimulationBase/MCTruth.h"
#include "SimulationBase/MCParticle.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/SpacePoint.h"
#include "RecoAlg/SpacePointAlg.h"
#include "Utilities/DetectorProperties.h"
#include "Geometry/Geometry.h"
// STL
#include <set>
#include <vector>

// ROOT
#include <TString.h>
#include <TTree.h>

namespace cluster
{
  class ClusterMatchAlg {

  public:

    /// Enum switch for various matching methods
    enum MatchMethod_t { kRoughZ = 0,      ///< Rough-Z comparison method ... see Match_RoughZ() description
			 kRoughT,          ///< Rough-Time comparison method ... see Match_RoughTime() description
			 kSpacePoint,      ///< Use SpacePoint finder algorithm ... see Match_SpacePoint() description
			 kSumCharge,       ///< Use summed charge comparison ... see Match_SumCharge() description
			 kMATCH_METHOD_MAX ///<
    };

    /**
       Local struct data container to store cluster's basic information
       based on hits that consist the cluster. Looping over hit pointer
       occurs when we create art::PtrVector<recob::Hit> from input file.
       All information that is based on hits and is used for cluster-matching
       should be extracted from there to maximize I/O efficiency as looping
       over hits takes time.  
       In other words... all hits related variables should be stored here!
    */
    struct cluster_info {
      
      unsigned short cluster_index; ///< Cluster's index position in the input cluster vector array
      geo::View_t view;             ///< Wire plane ID
      unsigned short wire_max;      ///< Maximum wire number in this cluster
      unsigned short wire_min;      ///< Minimum wire number in this cluster
      double start_time_max;        ///< Maximum "start time" among all hits in this cluster
      double peak_time_max;         ///< Maximum "peak time"  among all hits in this cluster
      double end_time_max;          ///< Maximum "end time"   among all hits in this cluster
      double start_time_min;        ///< Minimum "start time" among all hits in this cluster
      double peak_time_min;         ///< Minimum "peak time"  among all hits in this cluster
      double end_time_min;          ///< Minimum "end time"   among all hits in this cluster
      double sum_charge;            ///< Summed charge among all hits in this cluster
      
      /// Constructor with cluster's index ID
      cluster_info(unsigned short index) {
	cluster_index = index;
	view = geo::kUnknown;
	wire_max = 0;
	wire_min = 0xffff;
	start_time_max = peak_time_max = end_time_max = -1.;
	start_time_min = peak_time_min = end_time_min = 1.e9;
	sum_charge = -1.;
      };
      
      /// Default constructor
      cluster_info(){
	cluster_info(0xffff);
      };
    };
    
  public:

    /// Default constructor with fhicl parameters
    ClusterMatchAlg(fhicl::ParameterSet const& pset);
    //ClusterMatchAlg(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);

    /// Default destructor
    virtual ~ClusterMatchAlg(){};

    /// Method to report the current configuration
    void ReportConfig() const;

    /// Method to specify input cluster's module name (necessary)
    void SetClusterModName(std::string name);

    /// Method to specify input MCTruth's module name (optional)
    void SetMCTruthModName(std::string name)
    {_ModName_MCTruth = name;}

    /** 
	Method to set boolean vector which tells which cluster to ignroe.
	The input vector must have the same length as the cluster pointer vector
	that can be retrieved with the set cluster module name. Index with false
	entry is ignored from cluster vector.
    */
    void SetInputBoolArray(std::vector<bool> input_cluster_index)
    {_input_cluster_index = input_cluster_index;}

    /**
       Method to run matching algorithms for three planes. 
       Event info must be provided prior to this function call through FillEventInfo() 
       function call. If the function is called more than once w/o supplying the new 
       art::Event object, it does not perform any new matching unless a user explicitly 
       calls ClearEventInfo() and then fill event info again though FillEventInfo().
    */
    void MatchThreePlanes();
    
    /// Two plane version of cluster matching method.
    void MatchTwoPlanes();

    /// Method to retrieve matched cluster combinations. The format is [wire_plane][cluster_index]
    std::vector<std::vector<unsigned int> > GetMatchedClusters() const;

    /// Method to fill input data ... to be called before Match function call.
    void FillEventInfo(const art::Event &evt);

    /// Method to clear event-wise information
    void ClearEventInfo();

  protected:

    /// Internal method to fill cluster-wise information
    bool FillClusterInfo(const art::Event &evt);
    /// Internal method to fill MCTruth information when available
    void FillMCInfo(const art::Event &evt);
    /// Internal method, called only once, to fill detector-wise information
    void PrepareDetParams();
    /// Internal method to create output TTree for quality checking of the algorithm
    void PrepareTree();
    /// Internal method to check if the stored event is the same as what was provided previously
    inline bool SameEvent(const art::Event &evt) const
    { return (_run == evt.run() && _subrun == evt.subRun() && _event_id == evt.id().event() ); }
    //
    // Filtering methods
    //

    /**
       Match clusters based on min/max Z boundary information. It checks clusters' overlap 
       along Z spatial coordinate based on 2 input cluster information.
    */
    bool Match_RoughZ(const cluster_info &ci1,  const cluster_info &ci2,
		      const geo::View_t v1,     const geo::View_t v2 ) const;

    /// Checks min/max hit timing among two clusters and make sure there is an overlap
    bool Match_RoughTime(const cluster_info &ci1, const cluster_info &ci2);

    /**
      Checks min/max hit timing among three clusters and make sure there is an overlap.
      If _overlay_tratio_cut is set, then overlapped-time / cluster-timespan fraction
      is compared to the set cut value to claim a match.
    */
    //bool Match_RoughTime(const cluster_info &ci1, const cluster_info &ci2, const cluster_info &ci3);
    
    /**
       Checks the ratio of two clusters' summed charge. If the ratio is within 1 +/- set cut value,
       two clusters are considered to match.
    */
    bool Match_SumCharge(const cluster_info &uc, const cluster_info &vc);

    /**
       Cluster matching using space points. This can be slow as we run SpacePointFinder
       algorithm per cluster pair.  Three clusters (U, V, W) are considerd to "match" if 
       there found N spacepoints using hits in them and N > min_nsps where "min_nsps" is
       the cut value you can set in SetNSpacePointCut() method.
    */
    bool Match_SpacePoint(const size_t uindex, const size_t vindex, const size_t windex=0);

    //
    // Cut parameter values 
    //
    size_t _num_sps_cut;        ///< Number of SpacePoint used to cut in Match_SpacePoint method
    double _overlay_tratio_cut; ///< Minimum overlayed time fraction among two clusters used in Match_RoughTime method
    double _qratio_cut;         ///< Maximum difference among clusters' charge sum used in Match_SumCharge method

    //
    // Resulting data holder for matched cluster indexes
    //
    std::vector<unsigned int> _matched_uclusters_v; ///< U plane matched clusters' index
    std::vector<unsigned int> _matched_vclusters_v; ///< V plane matched clusters' index
    std::vector<unsigned int> _matched_wclusters_v; ///< W plane matched clusters' index

    //
    // Run control variables
    //
    std::vector<bool> _input_cluster_index; ///< Boolean vector to ignore some clusters
    bool _match_methods[kMATCH_METHOD_MAX]; ///< Boolean list for enabled algorithms
    bool _event_var_filled;                 ///< Boolean to keep track of whether the even data is received or not
    bool _debug_mode;                       ///< Boolean to enable debug mode (call all enabled matching methods)
    unsigned int _event_id;  ///< Processed event's counter
    unsigned int _run;       ///< Processed event's run number
    unsigned int _subrun;    ///< Processed event's subrun number

    std::string _ModName_Cluster; ///< Cluster producer's module name
    std::string _ModName_MCTruth; ///< MCTruth producer's module name

    std::vector<art::PtrVector<recob::Hit> > _uhits_v; ///< Local Hit pointer vector container ... U-plane
    std::vector<art::PtrVector<recob::Hit> > _vhits_v; ///< Local Hit pointer vector container ... V-plane
    std::vector<art::PtrVector<recob::Hit> > _whits_v; ///< Local Hit pointer vector container ... W-plane

    std::vector<cluster_info> _ucluster_v; ///< Local cluster data container... U-plane
    std::vector<cluster_info> _vcluster_v; ///< Local cluster data container... V-plane
    std::vector<cluster_info> _wcluster_v; ///< Local cluster data container... W-plane

    trkf::SpacePointAlg* _sps_algo; ///< SpacePointFinder algorithm pointer

    //
    // Quality control parameters to be saved in the TTree
    //
    TTree* _tree;
    double _mc_E;
    double _mc_Px;
    double _mc_Py;
    double _mc_Pz;
    double _mc_Vx;
    double _mc_Vy;
    double _mc_Vz;
    int _pdgid;
    unsigned int _tot_planes;
    unsigned short _tot_u;
    unsigned short _tot_v;
    unsigned short _tot_w;
    unsigned short _tot_pass_qsum;
    unsigned short _tot_pass_t;
    unsigned short _tot_pass_z;
    unsigned short _tot_pass_sps;
    std::vector<uint16_t> _u_nhits_v;
    std::vector<uint16_t> _v_nhits_v;
    std::vector<uint16_t> _w_nhits_v;
    std::vector<uint16_t> _nsps;
    std::vector<double>   _qratio_v;
    //std::vector<double>   _qratio_v2;
    std::vector<double>   _uv_tratio_v;
    std::vector<double>   _vw_tratio_v;
    std::vector<double>   _wu_tratio_v;
  }; // class ClusterMatchAlg
  
} //namespace cluster
#endif
