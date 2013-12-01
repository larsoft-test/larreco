////////////////////////////////////////////////////////////////////////
//
//  ClusterMatchAlg source
//
////////////////////////////////////////////////////////////////////////

#ifndef CLUSTERMATCHALG_CC
#define CLUSTERMATCHALG_CC

#include "ClusterMatchAlg.h"

namespace cluster{

  //##################################################################
  ClusterMatchAlg::ClusterMatchAlg(fhicl::ParameterSet const& pset):
    _ModName_Cluster(""),
    _ModName_MCTruth("")
  //##################################################################
  {    

    _debug_mode         = pset.get<bool>   ("DebugMode");
    _num_sps_cut        = pset.get<size_t> ("CutParam_NumSpacePoint");
    _overlay_tratio_cut = pset.get<double> ("CutParam_OverlayTimeFraction");
    _qratio_cut         = pset.get<double> ("CutParam_SumChargeRatio");
    std::vector<size_t> algo_list = pset.get<std::vector<size_t> > ("MatchAlgoList");

     _sps_algo = new trkf::SpacePointAlg(pset.get<fhicl::ParameterSet>("SpacePointAlg"));
    _tree = 0;
    _event_var_filled = false;

    for(size_t i=0; i<(size_t)(kMATCH_METHOD_MAX); ++i)

      _match_methods[i]=false;

    for(auto const v : algo_list) {

      if(v >= (size_t)(kMATCH_METHOD_MAX))

	mf::LogError("ClusterMatchAlg")<<Form("Invalid algorithm enum: %zu",v);
      
      else _match_methods[v]=true;
    }

    ReportConfig();

    ClearEventInfo();
  }

  //########################################
  void ClusterMatchAlg::ReportConfig() const
  //########################################
  {
    std::ostringstream msg;
    msg
      << std::endl
      << " ClusterMatchAlg Configuration:              " << std::endl
      << "---------------------------------------------" << std::endl;
    msg << " Debug Mode ... " << (_debug_mode ? "enabled!" : "disabled!") << std::endl;
    msg << " RoughZ ....... " << (_match_methods[kRoughZ] ? "enabled!" : "disabled!") << std::endl;
    msg << " RoughT ....... " << (_match_methods[kRoughT] ? "enabled!" : "disabled!") << std::endl;
    msg << " SpacePoint ... " << (_match_methods[kSpacePoint] ? "enabled!" : "disabled!") << std::endl;
    msg << " SumCharge .... " << (_match_methods[kSumCharge] ? "enabled!" : "disabled!") << std::endl;
    msg << std::endl;
    msg
      << " Overlay-Time Fraction Cut : " << _overlay_tratio_cut << std::endl
      << " Charge-Ratio Diff. Cut    : " << _qratio_cut << std::endl
      << " Minimum # of SpacePoint   : " << _num_sps_cut << std::endl
      << std::endl;
    msg
      << "---------------------------------------------" << std::endl;

    mf::LogWarning("ClusterMatchAlg")<<msg.str();

  }

  //##################################################################
  void ClusterMatchAlg::ClearEventInfo()
  //##################################################################
  {
    // Clear input event data holders
    _ucluster_v.clear();
    _vcluster_v.clear();
    _wcluster_v.clear();
    
    _uhits_v.clear();
    _vhits_v.clear();
    _whits_v.clear();

    _event_var_filled = false;

    // Clear result data holders
    _matched_uclusters_v.clear();
    _matched_vclusters_v.clear();
    _matched_wclusters_v.clear();

    /// Run control variables
    _event_id = 0;
    _run = 0;
    _subrun = 0;

    /// QC variables
    _mc_E  = 0;
    _mc_Px = 0;
    _mc_Py = 0;
    _mc_Pz = 0;
    _mc_Vx = 0;
    _mc_Vy = 0;
    _mc_Vz = 0;
    _pdgid = 0;
    _tot_u = 0;
    _tot_v = 0;
    _tot_w = 0;
    _tot_pass_z = 0;
    _tot_pass_t = 0;
    _tot_pass_sps = 0;
    _tot_pass_qsum = 0;
    _qratio_v.clear();
    _uv_tratio_v.clear();
    _vw_tratio_v.clear();
    _wu_tratio_v.clear();
    _u_nhits_v.clear();
    _v_nhits_v.clear();
    _w_nhits_v.clear();
    _nsps.clear();

  }
  
  //##################################################################
  void ClusterMatchAlg::SetClusterModName(std::string name)
  //##################################################################
  {    
    if(name != _ModName_Cluster) {
      _ModName_Cluster=name;
      ClearEventInfo();
    }

  }

  //##################################################################
  void ClusterMatchAlg::PrepareTree()
  //##################################################################
  {
    art::ServiceHandle<art::TFileService> fileService;
    if(!_tree){
      _tree = fileService->make<TTree>("match_tree","");
      _tree->Branch("mc_E",&_mc_E,"mc_E/D");
      _tree->Branch("mc_Px",&_mc_Px,"mc_Px/D");
      _tree->Branch("mc_Py",&_mc_Py,"mc_Py/D");
      _tree->Branch("mc_Pz",&_mc_Pz,"mc_Pz/D");
      _tree->Branch("mc_Vx",&_mc_Vx,"mc_Vx/D");
      _tree->Branch("mc_Vy",&_mc_Vy,"mc_Vy/D");
      _tree->Branch("mc_Vz",&_mc_Vz,"mc_Vz/D");

      _tree->Branch("pdgid",&_pdgid,"pdgid/I");
      _tree->Branch("tot_u",&_tot_u,"tot_u/s");
      _tree->Branch("tot_v",&_tot_v,"tot_v/s");
      _tree->Branch("tot_w",&_tot_w,"tot_w/s");
      _tree->Branch("tot_pass_t",&_tot_pass_t,"tot_pass_t/s");
      _tree->Branch("tot_pass_z",&_tot_pass_z,"tot_pass_z/s");
      _tree->Branch("tot_pass_sps",&_tot_pass_sps,"tot_pass_sps/s");
      _tree->Branch("tot_pass_qsum",&_tot_pass_qsum,"tot_pass_qsum/s");

      _tree->Branch("uv_tratio_v","std::vector<double>",&_uv_tratio_v);
      _tree->Branch("vw_tratio_v","std::vector<double>",&_vw_tratio_v);
      _tree->Branch("wu_tratio_v","std::vector<double>",&_wu_tratio_v);

      _tree->Branch("qratio_v", "std::vector<double>",&_qratio_v);
      _tree->Branch("u_nhits_v", "std::vector<UShort_t>",&_u_nhits_v);
      _tree->Branch("v_nhits_v", "std::vector<UShort_t>",&_v_nhits_v);
      _tree->Branch("w_nhits_v", "std::vector<UShort_t>",&_w_nhits_v);
      _tree->Branch("nsps",  "std::vector<UShort_t>",&_nsps);
    }

  }

  //##########################################################################################
  void ClusterMatchAlg::FillEventInfo(const art::Event &evt)
  //##########################################################################################  
  {
    if(!(_event_id) && !_event_var_filled) {
      // 1st call ever
      PrepareTree();
    }
    else if(_event_var_filled && SameEvent(evt)) return;
    ClearEventInfo();
    FillMCInfo(evt);
    FillClusterInfo(evt);
    _event_id = evt.id().event();
    _run = evt.run();
    _subrun = evt.subRun();
    _event_var_filled = true;
  }

  //##########################################################################################
  void ClusterMatchAlg::FillMCInfo(const art::Event &evt)
  //##########################################################################################  
  {
    if(!_ModName_MCTruth.size()) return;

    std::vector<const simb::MCTruth*> mciArray;

    try {

      evt.getView(_ModName_MCTruth,mciArray);

    }catch (art::Exception const& e) {

      if (e.categoryCode() != art::errors::ProductNotFound ) throw;

    }

    for(size_t i=0; i < mciArray.size(); ++i){
      
      if(i==1) {
	mf::LogWarning("ClusterMatchAlg")<<" Ignoring > 2nd MCTruth in MC generator...";
	break;
      }
      const simb::MCTruth* mci_ptr(mciArray.at(i));

      for(size_t j=0; j < (size_t)(mci_ptr->NParticles()); ++j){

	if(j==1) {
	  mf::LogWarning("ClusterMatchAlg")<<" Ignoring > 2nd MCParticle in MC generator...";
	  break;
	}

	const simb::MCParticle part(mci_ptr->GetParticle(j));
	
	_pdgid = part.PdgCode();
	_mc_E  = part.E();
	_mc_Px = part.Px();
	_mc_Py = part.Py();
	_mc_Pz = part.Pz();
	_mc_Vx = part.Vx();
	_mc_Vy = part.Vy();
	_mc_Vz = part.Vz();
      }
    }
  }


  //##########################################################################################
  void ClusterMatchAlg::FillClusterInfo(const art::Event &evt)
  //##########################################################################################
  {

    //
    // Pull a vector of cluster & hit pointer from input (std::vector<art::Ptr>)
    //
    std::vector<art::Ptr<recob::Cluster> > cluster_ptr_v;
    art::Handle< std::vector<recob::Cluster> > cluster_handle;
    evt.getByLabel(_ModName_Cluster, cluster_handle);
    if(cluster_handle.isValid()) {
      art::fill_ptr_vector(cluster_ptr_v, cluster_handle);
    }

    //
    // Next, we loop over clusters & save art::PtrVector<recob::Hit> per cluster
    // into the local data container. In the loop, we also save some relevant
    // cluster-wise data in a dedicated cluster_info struct defined in the header.
    //

    // Create association 
    art::FindManyP<recob::Hit> hit_m(cluster_handle, evt, _ModName_Cluster);

    // Ask Geo about time-offset among different wire planes ... used to correct timing
    // difference among different wire planes in the following loop.
    art::ServiceHandle<util::DetectorProperties> det_h;
    double time_offset_uplane = det_h->GetXTicksOffset(geo::kU,0,0);
    double time_offset_vplane = det_h->GetXTicksOffset(geo::kV,0,0);
    double time_offset_wplane = det_h->GetXTicksOffset(geo::kW,0,0);

    // Start looping over clusters
    for(size_t i = 0; i < cluster_ptr_v.size(); ++i) {

      const art::Ptr<recob::Cluster> cluster_ptr = cluster_ptr_v.at(i);
      cluster_info ci((unsigned short)(i));
      ci.view = cluster_ptr->View();

      const std::vector<art::Ptr<recob::Hit> > hit_v = hit_m.at(i);
      art::PtrVector<recob::Hit> hit_ptrv;
      hit_ptrv.reserve(hit_v.size());
      
      // Loop over hits in this cluster
      for(auto const hit : hit_v){
	unsigned int wire = hit->WireID().Wire;
	double tstart = hit->StartTime();
	double tpeak = hit->PeakTime();
	double tend = hit->EndTime();
	ci.sum_charge += hit->Charge();

	ci.wire_max = (ci.wire_max < wire) ? wire : ci.wire_max;
	ci.wire_min = (ci.wire_min > wire) ? wire : ci.wire_min;

	ci.start_time_max = ( ci.start_time_max < tstart ) ? tstart : ci.start_time_max;
	ci.peak_time_max  = ( ci.peak_time_max  < tpeak  ) ? tpeak  : ci.peak_time_max;
	ci.end_time_max   = ( ci.end_time_max   < tend   ) ? tend   : ci.end_time_max;

	ci.start_time_min = ( ci.start_time_min > tstart) ? tstart : ci.start_time_min;
	ci.peak_time_min  = ( ci.peak_time_min  > tpeak)  ? tpeak  : ci.peak_time_min;
	ci.end_time_min   = ( ci.end_time_min   > tend)   ? tend   : ci.end_time_min;

	hit_ptrv.push_back(hit);
      }

      // Save created art::PtrVector & cluster_info struct object
      switch(ci.view){
      case geo::kU: 
	_uhits_v.push_back(hit_ptrv);
	ci.start_time_max -= time_offset_uplane;
	ci.peak_time_max  -= time_offset_uplane;
	ci.end_time_max   -= time_offset_uplane;
	ci.start_time_min -= time_offset_uplane;
	ci.peak_time_min  -= time_offset_uplane;
	ci.end_time_min   -= time_offset_uplane;
	_ucluster_v.push_back(ci);
	break;
      case geo::kV: 
	_vhits_v.push_back(hit_ptrv);
	ci.start_time_max -= time_offset_vplane;
	ci.peak_time_max  -= time_offset_vplane;
	ci.end_time_max   -= time_offset_vplane;
	ci.start_time_min -= time_offset_vplane;
	ci.peak_time_min  -= time_offset_vplane;
	ci.end_time_min   -= time_offset_vplane;
	_vcluster_v.push_back(ci);
	break;
      case geo::kW:
	_whits_v.push_back(hit_ptrv);
	ci.start_time_max -= time_offset_wplane;
	ci.peak_time_max  -= time_offset_wplane;
	ci.end_time_max   -= time_offset_wplane;
	ci.start_time_min -= time_offset_wplane;
	ci.peak_time_min  -= time_offset_wplane;
	ci.end_time_min   -= time_offset_wplane;
	_wcluster_v.push_back(ci);
	break;
      default:
	mf::LogError("ClusterMatchAlg")<<Form("Found an invalid plane ID: %d",cluster_ptr->View());
	continue;
      }

      
    }

    mf::LogWarning("ClusterMatchAlg")
      << Form("Found (U,V,W) = (%zu,%zu,%zu) clusters...",
	      _uhits_v.size(),
	      _vhits_v.size(),
	      _whits_v.size())
      << std::endl;
    
    _tot_u = _ucluster_v.size();
    _tot_v = _vcluster_v.size();
    _tot_w = _wcluster_v.size();
    
  }

  //########################################################################################
  bool ClusterMatchAlg::Match_RoughZ(const cluster_info &ci1,  const cluster_info &ci2,
				     const geo::View_t v1,     const geo::View_t v2 ) const
  //########################################################################################
  {
    art::ServiceHandle<geo::Geometry> geo_h;
    double y, z_min, z_max;
    y = z_min = z_max = -1;
    geo_h->IntersectionPoint(ci1.wire_min, ci2.wire_min, v1, v2, 0, 0, y, z_min);
    geo_h->IntersectionPoint(ci1.wire_max, ci2.wire_max, v1, v2, 0, 0, y, z_max);
    return (z_max > z_min);
  }

  //###########################################################################################
  bool ClusterMatchAlg::Match_RoughTime(const cluster_info &ci1, const cluster_info &ci2)
  //###########################################################################################
  {
    //return (!(ci1.end_time_max < ci2.start_time_min || ci2.end_time_max < ci1.start_time_min));
    double time_overlay = std::min(ci1.end_time_max,ci2.end_time_max) - std::max(ci1.start_time_min,ci2.start_time_min);

    //if(time_overlay <= 0 && !_debug_mode) return false;

    double overlay_tratio = time_overlay /  (ci1.end_time_max - ci1.start_time_min + ci2.end_time_max - ci2.start_time_min) * 2.;
    
    if( (ci1.view==geo::kU && ci2.view==geo::kV) || (ci1.view==geo::kV && ci2.view==geo::kU) )
      _uv_tratio_v.push_back(overlay_tratio);
    else if( (ci1.view==geo::kV && ci2.view==geo::kW) || (ci1.view==geo::kW && ci2.view==geo::kV) )
      _vw_tratio_v.push_back(overlay_tratio);
    else if( (ci1.view==geo::kW && ci2.view==geo::kU) || (ci1.view==geo::kU && ci2.view==geo::kW) )
      _wu_tratio_v.push_back(overlay_tratio);

    return (overlay_tratio > _overlay_tratio_cut);
  }

  //###################################################################################
  bool ClusterMatchAlg::Match_SumCharge(const cluster_info &uc, const cluster_info &vc)
  //###################################################################################
  {
    double qratio = (uc.sum_charge)/(vc.sum_charge);
    
    // For quality check log
    _qratio_v.push_back(qratio);

    return ( (1 - _qratio_cut) < qratio && (qratio) < (1 + _qratio_cut) );
  }

  //#####################################################################################################
  bool ClusterMatchAlg::Match_SpacePoint(const size_t uindex, const size_t vindex, const size_t windex)
  //#####################################################################################################
  {
    bool use_wplane(_wcluster_v.size());

    if( uindex >= _ucluster_v.size() ||
	vindex >= _vcluster_v.size() ||
	(use_wplane && (windex >= _wcluster_v.size())) ) {
      
      mf::LogError("ClusterMatchAlg")
	<< std::endl
	<< Form("Requested to cluster-index (U,V,W) = (%zu,%zu,%zu) where max-length is (%zu,%zu,%zu)",
		uindex, vindex, windex, _ucluster_v.size(), _vcluster_v.size(), _wcluster_v.size())
	<< std::endl;
      return false;
    }

    // Define a time range in which hits are used for spacepoint finding ... here "peak time" is the relevant one
    double trange_min = std::min(_ucluster_v.at(uindex).peak_time_min,_vcluster_v.at(vindex).peak_time_min);
    if(use_wplane) trange_min = std::min(trange_min, _wcluster_v.at(windex).peak_time_min);

    double trange_max = std::max(_ucluster_v.at(uindex).peak_time_max,_vcluster_v.at(vindex).peak_time_max);
    if(use_wplane) trange_max = std::max(trange_max,_wcluster_v.at(windex).peak_time_max);
    
    // Space-point algorithm applies additional dT
    trange_min -= _sps_algo->maxDT();
    trange_max += _sps_algo->maxDT();
    
    // Make PtrVector<recob::Hit> for relevant Hits
    art::PtrVector<recob::Hit> hit_group;
    size_t max_size = _uhits_v.at(uindex).size() + _vhits_v.at(vindex).size();
    if(use_wplane) max_size += _whits_v.at(windex).size();
    hit_group.reserve(max_size);

    // Loop over hits in U-plane
    for(auto const hit : _uhits_v.at(uindex)) {
      if(hit->PeakTime() < trange_min) continue;
      if(hit->PeakTime() > trange_max) continue;
      hit_group.push_back(hit);
    }
    // Check if any hit found in this plane
    size_t u_nhits = hit_group.size();
    if( !u_nhits && !_debug_mode) return false;
    
    // Loop over hits in V-plane
    for(auto const hit: _vhits_v.at(vindex)) {
      if(hit->PeakTime() < trange_min) continue;
      if(hit->PeakTime() > trange_max) continue;
      hit_group.push_back(hit);
    }
    // Check if any hit found in this plane
    size_t v_nhits = hit_group.size() - u_nhits;
    if( !(v_nhits) && !_debug_mode) return false;
    
    // Loop over hits in W-plane
    for(auto const hit: _whits_v.at(windex)) {
      if(hit->PeakTime() < trange_min) continue;
      if(hit->PeakTime() > trange_max) continue;
      hit_group.push_back(hit);
    }
    // Check if any hit found in this plane
    size_t w_nhits = hit_group.size() - u_nhits - v_nhits;
    if( !(w_nhits) && use_wplane && !_debug_mode) return false;

    // Run SpacePoint finder algo
    std::vector<recob::SpacePoint> sps;
    if(u_nhits && v_nhits && (w_nhits && use_wplane)) {
      _sps_algo->clearHitMap();
      _sps_algo->makeSpacePoints(hit_group,sps);
    }

    size_t nsps = sps.size();
    _u_nhits_v.push_back(u_nhits);
    _v_nhits_v.push_back(v_nhits);
    _w_nhits_v.push_back(w_nhits);
    _nsps.push_back(nsps);

    if( nsps < _num_sps_cut ) return false;
    return true;
  }

  //#################################################################################
  std::vector<std::vector<unsigned int> > ClusterMatchAlg::GetMatchedClusters() const
  //#################################################################################
  {
    std::vector<std::vector<unsigned int> > result;
    result.push_back(_matched_uclusters_v);
    result.push_back(_matched_vclusters_v);
    result.push_back(_matched_wclusters_v);
    return result;
  }

  //#######################################################################
  void ClusterMatchAlg::MatchTwoPlanes()
  //#######################################################################
  {
    // Return immediately if this method was already called
    if(_matched_uclusters_v.size()) return;

    bool overlay_2d = false;
    for(size_t uci_index=0; uci_index<_ucluster_v.size(); ++uci_index) {

      for(size_t vci_index=0; vci_index<_vcluster_v.size(); ++vci_index) {

	overlay_2d = true;

	// Apply cuts

	// Rough z-position overlay cut
	if(_match_methods[kRoughZ]) {

	  if(Match_RoughZ(_ucluster_v.at(uci_index), _vcluster_v.at(vci_index), geo::kU, geo::kV))
	    _tot_pass_z++;
	  else if(!_debug_mode) continue;
	  else overlay_2d = false;
	}

	// Sum charge cut
	if(_match_methods[kSumCharge]) {
	  
	  if(Match_SumCharge(_ucluster_v.at(uci_index),_vcluster_v.at(vci_index)))
	    _tot_pass_qsum++;
	  else if(!_debug_mode) continue;
	  else overlay_2d = false;
	}

	// Rough time overlap cut
	if(_match_methods[kRoughT]) {
	  
	  if(Match_RoughTime(_ucluster_v.at(uci_index), _vcluster_v.at(vci_index)))
	    _tot_pass_t++;
	  else if(!_debug_mode) continue;
	  else overlay_2d = false;
	}

	// SpacePoint cut
	if(_match_methods[kSpacePoint]) {
	  
	  if(Match_SpacePoint(uci_index, vci_index)) 
	    _tot_pass_sps++;
	  else if(!_debug_mode) continue;
	  else overlay_2d = false;
	}	      
	
	if(overlay_2d) {
	  _matched_uclusters_v.push_back((unsigned int)(_ucluster_v[uci_index].cluster_index));
	  _matched_vclusters_v.push_back((unsigned int)(_vcluster_v[vci_index].cluster_index));
	}
      } // end of ... _vcluster_v loop
    } // end of ... _ucluster_v loop
    mf::LogWarning("ClusterMatchAlg")<<Form("Found %zu matched cluster pairs...",_matched_uclusters_v.size());
    _tree->Fill();
  }

  //#######################################################################
  void ClusterMatchAlg::MatchThreePlanes()
  //#######################################################################
  {
    // Return immediately if this method was already called
    if(_matched_wclusters_v.size()) return;
    // Handle the case MatchTwoPlanes() was already called (tolerate)
    _matched_wclusters_v.clear();
    _matched_uclusters_v.clear();
    _matched_vclusters_v.clear();

    bool overlay_2d=true;
    bool overlay_3d=true;
    // Loop over all possible u-v-w cluster combination
    for(size_t uci_index=0; uci_index<_ucluster_v.size(); ++uci_index) {

      for(size_t vci_index=0; vci_index<_vcluster_v.size(); ++vci_index) {

	// Apply cuts that can be done with U&V planes here
	overlay_2d = true;

	// Rough z-position overlay cut
	if(_match_methods[kRoughZ]) {

	  if(Match_RoughZ(_ucluster_v.at(uci_index), _vcluster_v.at(vci_index), geo::kU, geo::kV))
	    _tot_pass_z++;
	  else if(!_debug_mode) continue;
	  else overlay_2d = false;
	}

	// Sum charge cut
	if(_match_methods[kSumCharge]) {
	  
	  if(Match_SumCharge(_ucluster_v.at(uci_index),_vcluster_v.at(vci_index)))
	    _tot_pass_qsum++;
	  else if(!_debug_mode) continue;
	  else overlay_2d = false;
	}

	for(size_t wci_index=0; wci_index<_wcluster_v.size(); ++wci_index) {

	  overlay_3d = overlay_2d;
	  // Apply cuts that requires 3 planes here

	  // Rough time overlap cut
	  if(_match_methods[kRoughT]) {
	    
	    bool rough_time_match = Match_RoughTime(_ucluster_v.at(uci_index), _vcluster_v.at(vci_index));
	    if(!_debug_mode && !rough_time_match) continue;
	    
	    rough_time_match = (Match_RoughTime(_vcluster_v.at(vci_index), _wcluster_v.at(wci_index)) && rough_time_match);
	    if(!_debug_mode && !rough_time_match) continue;
	    
	    rough_time_match = (Match_RoughTime(_wcluster_v.at(wci_index), _ucluster_v.at(uci_index)) && rough_time_match);

	    overlay_3d = overlay_3d && rough_time_match;
	    if(rough_time_match) _tot_pass_t++;
	    else if(!_debug_mode) continue;
	  }

	  // SpacePoint cut
	  if(_match_methods[kSpacePoint]) {

	    if(Match_SpacePoint(uci_index, vci_index, wci_index))
	      _tot_pass_sps++;
	    else if(!_debug_mode) continue; 
	    else overlay_3d = false;
	  }	      

	  if(overlay_3d) {
	    _matched_uclusters_v.push_back((unsigned int)(_ucluster_v[uci_index].cluster_index));
	    _matched_vclusters_v.push_back((unsigned int)(_vcluster_v[vci_index].cluster_index));
	    _matched_wclusters_v.push_back((unsigned int)(_wcluster_v[wci_index].cluster_index));
	  }
	} // end of ... _wcluster_v loop
      } // end of ... _vcluster_v loop
    } // end of ... _ucluster_v loop
    mf::LogWarning("ClusterMatchAlg")<<Form("Found %zu matched cluster pairs...",_matched_uclusters_v.size());
    _tree->Fill();
  }

} // namespace match

#endif 
