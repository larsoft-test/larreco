#include "RecoAlg/ClusterParamsAlg.h"
#include "TF1.h"
#include "TH2F.h"
#include "TH1F.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "RecoBase/Hit.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Services/Optional/TFileService.h"

cluster::ClusterParamsAlg::ClusterParamsAlg(fhicl::ParameterSet const& pset)
 : fHBAlg(pset.get< fhicl::ParameterSet >("HoughBaseAlg"))
{

/**Get TFileService and define output Histograms*/
art::ServiceHandle<art::TFileService> tfs;
tgxtest=tfs->make<TH2F>(Form("charge hi distrib"),"charge hi distribution per wires",
 			  10,0,10,10,0,10);
linefittest_cm=tfs->make<TF1>(Form("linefittest_cm"),"pol1", 0,10); 


//fChargeCutoffThreshold	=pset.get<   std::vector<double> > ("ChargeCutoffThreshold");



}




void cluster::ClusterParamsAlg::reconfigure(fhicl::ParameterSet const& pset) 
{
  fChargeCutoffThreshold	=pset.get<   std::vector<double> > ("ChargeCutoffThreshold");    
  fSelectBoxSizePar    		=pset.get<   double > ("SelectBoxSizePar");
  fSelectBoxSizePerp		=pset.get<   double > ("SelectBoxSizePerp");
  //fSelectBoxDiff		=pset.get<   double > ("SelectBoxDiff");
  fMinHitListSize		=pset.get<   double > ("MinHitListSize"); //40;
  fOutlierRadius		=pset.get<   double > ("fOutlierRadius"); //5;
  fForceRightGoing 		=pset.get<   bool > ("ForceRightGoing");
  fHBAlg.reconfigure(pset.get< fhicl::ParameterSet >("HoughBaseAlg"));
  
  fWirePitch = geo->WirePitch(0,1,0);    //wire pitch in cm
  fTimeTick=detp->SamplingRate()/1000.; 

  //define conversion constants
  double fDriftVelocity=larp->DriftVelocity(larp->Efield(),larp->Temperature());
  fWiretoCm=fWirePitch;
  fTimetoCm=fTimeTick*fDriftVelocity;
  fWireTimetoCmCm=(fTimeTick*fDriftVelocity)/fWirePitch;
  
 }




//---------------------------------------------------------------------------------------------
//Calculate Minimum and Maximum wires and times of the cluster. And Mean charge
//--------------------------------------------------------------------------------------------

void cluster::ClusterParamsAlg::FindMinMaxWiresTimes(double &MinWire, double &MaxWire,double &MinTime,double &MaxTime,double &MeanCharge,std::vector< art::Ptr < recob::Hit> > hitlist){

  unsigned int wire,tpc, cstat;
  unsigned int plane;

  double mcharge=0;
  
  if(hitlist.size()==0)
    return;
  
   for(std::vector < art::Ptr < recob::Hit > >::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
    double time = (*hitIter)->PeakTime();  
    //time_C -= (presamplings+10.1);
    GetPlaneAndTPC((*hitIter),plane,cstat,tpc,wire);
    
    if(time>MaxTime)
	MaxTime=time;

    if(time<MinTime)
	MinTime=time;
    
    if(wire>MaxWire)
	MaxWire=wire;

    if(wire<MinWire)
	MinWire=wire;

    mcharge+=(*hitIter)->Charge();
    
  }
  
  MeanCharge=mcharge/hitlist.size();
  
  //fTotalCharge=mcharge;
  
  //fClusterPlane=plane;  // time coordinate of last point for each plane
      
 // std::cout << " +++ maximums " << MinWire << " " << MaxWire << " " <<MinTime << " " << MaxTime << std::endl;
//  std::cout << " +++ meancharge " << MeanCharge << std::endl;
}



//---------------------------------------------------------------------------------------------
//Calculate rough axis of the cluster, using hits above the average charge. Calculate measure of verticalness of the cluster and return in goodness.
//--------------------------------------------------------------------------------------------

void cluster::ClusterParamsAlg::Find2DAxisRoughHighCharge(double &lineslope,double &lineintercept,double &goodness,std::vector < art::Ptr < recob::Hit> > hitlist){
 
  double kMinWire=99999,kMinTime=999999,kMaxWire=0,kMaxTime=0,kMeanCharge=0;
  FindMinMaxWiresTimes(kMinWire,kMaxWire,kMinTime,kMaxTime,kMeanCharge,hitlist);
  std::vector< art::Ptr<recob::Hit> > highhitlist=CreateHighHitlist(kMeanCharge,hitlist);
  std::cout << "$$$$$$$ highhitlist size: " << highhitlist.size() << " " << kMeanCharge <<  std::endl;
  Find2DAxisRough(lineslope,lineintercept,goodness,highhitlist);
  
}

//---------------------------------------------------------------------------------------------
//Calculate rough axis of the cluster. Calculate measure of verticalness of the cluster and return in goodness.
//--------------------------------------------------------------------------------------------

void cluster::ClusterParamsAlg::Find2DAxisRough(double &lineslope,double &lineintercept,double &goodness,std::vector < art::Ptr < recob::Hit> > hitlist){
 
  double time;
  unsigned int wire,tpc, cstat;
  unsigned int plane;

 
  if(hitlist.size()<5)
  {
   std::cout << " hitlist too small," <<hitlist.size() << " bailing " << std::endl;
    return;
  }
 // padding of the selected TGraph in wires and time. 
  int wirepad=20;
  int timepad=wirepad*fWiretoCm/fTimetoCm+0.5;
 
  double kMinWire=99999,kMinTime=999999,kMaxWire=0,kMaxTime=0,kMeanCharge=0;
  
  FindMinMaxWiresTimes(kMinWire,kMaxWire,kMinTime,kMaxTime,kMeanCharge,hitlist);
  
  
  //>\todo check that min/max wires are ok
  int nbinsx= (kMaxWire-kMinWire+2*wirepad)*fWiretoCm;  // nbins to have 
  int nbinsy= (kMaxTime-kMinTime+2.*(double)timepad)*fTimetoCm;  // nbins to have 
 
 
  tgxtest->Reset();
  tgxtest->SetBins(nbinsx,((double)kMinWire-(double)wirepad)*fWiretoCm,			((double)kMaxWire+(double)wirepad)*fWiretoCm,nbinsy,
			  (kMinTime-timepad)*fTimetoCm,(kMaxTime+timepad)*fTimetoCm);
    
  linefittest_cm->SetRange(((double)kMinWire-(double)wirepad)*fWiretoCm,										((double)kMaxWire+(double)wirepad)*fWiretoCm); 
//     linefit2_cm->SetRange(((double)kMinWire-(double)wirepad)*fWiretoCm,										((double)kMaxWire+(double)wirepad)*fWiretoCm);
//   }
  
  
  for(std::vector < art::Ptr < recob::Hit > >::const_iterator hitIter = hitlist.begin(); hitIter != 				hitlist.end();  hitIter++){
    
    time =  (*hitIter)->PeakTime();  
    GetPlaneAndTPC((*hitIter),plane,cstat,tpc,wire);
  
    tgxtest->Fill((double)wire*fWiretoCm,
			      time*fTimetoCm,(*hitIter)->Charge());
    
  }
  
 
  tgxtest->Fit(Form("linefittest_cm"),"+QMRNCFrob=0.95");


  
 double kRMS_wire=tgxtest->GetRMS(1);
 double kRMS_time=tgxtest->GetRMS(2);
 double kChisq;
 
  if(linefittest_cm->GetNDF())
    kChisq=linefittest_cm->GetChisquare()/linefittest_cm->GetNDF();
  else
    kChisq=1;
    
  double kCorrelation=tgxtest->GetCorrelationFactor();
  double kCovariance=tgxtest->GetCovariance();
  
  
  lineslope=linefittest_cm->GetParameter(1)/fWireTimetoCmCm;
  lineintercept=linefittest_cm->GetParameter(0)/fTimetoCm;

    if(kRMS_time && (kMaxTime-kMinTime) ){
  
      if(linefittest_cm->GetNDF() && tgxtest->GetCorrelationFactor() && tgxtest->GetCovariance())
	goodness=kRMS_time/kRMS_wire*((kMaxTime-kMinTime)*fTimetoCm)/((kMaxWire-kMinWire)*fWiretoCm)*kChisq/100000*1/kCorrelation*1/kCovariance*10;
      else
	goodness=-1;
    }

 
 //   std::cout << " old fit: " <<  linefittest_cm->GetParameter(1)/fWireTimetoCmCm << " " <<  linefittest_cm->GetParameter(0)/fTimetoCm << std::endl;
 

    
    
}


//---------------------------------------------------------------------------------------------
//Calculate rough axis of the cluster. Interface for version without the slope yet
//--------------------------------------------------------------------------------------------

void   cluster::ClusterParamsAlg::Find2DStartPointsBasic(std::vector< art::Ptr < recob::Hit> > hitlist,double &wire_start,double &time_start,double &wire_end,double &time_end)
{
  
 double lineslope, lineintercept,goodness;

 Find2DAxisRough(lineslope,lineintercept,goodness,hitlist);
 Find2DStartPointsBasic(lineslope,lineintercept,hitlist,wire_start,time_start,wire_end,time_end);
    
}



//---------------------------------------------------------------------------------------------
//Calculate rough axis of the cluster. Interface for version without the slope yet
//--------------------------------------------------------------------------------------------

void   cluster::ClusterParamsAlg::Find2DStartPointsHighCharge(std::vector< art::Ptr < recob::Hit> > hitlist,double &wire_start,double &time_start,double &wire_end,double &time_end)
{
  
 double lineslope, lineintercept,goodness;

 Find2DAxisRoughHighCharge(lineslope,lineintercept,goodness,hitlist);
 Find2DStartPointsBasic(lineslope,lineintercept,hitlist,wire_start,time_start,wire_end,time_end);
    
}



//---------------------------------------------------------------------------------------------
//Calculate rough axis of the cluster. Calculate measure of verticalness of the cluster and return in goodness.
//--------------------------------------------------------------------------------------------


void   cluster::ClusterParamsAlg::Find2DStartPointsBasic(double lineslope,double lineintercept,std::vector< art::Ptr < recob::Hit> > hitlist,double &wire_start,double &time_start,double &wire_end,double &time_end)
{
  
  //  double time;
  unsigned int wire,plane,tpc,cstat;
  double a,c;
 
  GetPlaneAndTPC((*hitlist.begin()),plane,cstat,tpc,wire);
   
  //get paramters of the straight line fit. (and rescale them to cm/cm)
 
   a=lineslope*fWireTimetoCmCm;
   c=lineintercept*fTimetoCm;
   
  //get slope of lines orthogonal to those found crossing the shower.
  double aprim=0;
  if(a){
    aprim=-1./a*fWireTimetoCmCm*fWireTimetoCmCm;
  }
 else
   aprim = -999999.;
  
 
 //std::cout << "in starting points basic, line slope, intercept, a,c: " << lineslope << " " << lineintercept << " " << a << " " << c << std::endl;
 // find extreme intercepts. For the time being we don't care which one is the start point and which one is the endpoint.
  
  double extreme_intercept_high,extreme_intercept_low;
  
  FindExtremeIntercepts(hitlist,aprim,extreme_intercept_high,extreme_intercept_low);
  double extreme_intercept_start=-999999;
  double extreme_intercept_end=999999;
  
 // int multiplier=1;   // +1 for positive angles, -1 for negative angles. to compensate that we are looking for either the highest (omega >0 ) or lowest (omega<0) intercept.
  
  if(a>=0) {  // for the time being assuming forward going showers
   //  multiplier=1;
     extreme_intercept_end=extreme_intercept_high;
     extreme_intercept_start=extreme_intercept_low;
  }
  else if(a<0){
   // multiplier=-1;
    extreme_intercept_start=extreme_intercept_high;
    extreme_intercept_end=extreme_intercept_low;
  }
  
  std::cout << " extreme intercepts: " << extreme_intercept_start << " " << extreme_intercept_end << std::endl;

  
  double wire_online_end,wire_online_begin,time_online_end,time_online_begin;
  
  gser.GetPointOnLineWSlopes(a,c,extreme_intercept_end,wire_online_end,time_online_end);
  gser.GetPointOnLineWSlopes(a,c,extreme_intercept_start,wire_online_begin,time_online_begin);
  
// std::cout << " ---------:::::::: wire_online begin point " << wire_online_begin << " " << time_online_begin << std::endl;
// std::cout << " ----------:::::::: wire_online end point " << wire_online_end << " " << time_online_end << std::endl;
  //calculate the first and last cluster points on the through line:
  
  art::Ptr<recob::Hit> startHit=FindClosestHit(hitlist, wire_online_begin,time_online_begin);
  art::Ptr<recob::Hit> endHit=FindClosestHit(hitlist, wire_online_end,time_online_end);
  
  GetPlaneAndTPC(startHit,plane,cstat,tpc,wire);
  wire_start=wire;
  time_start=startHit->PeakTime();
  
  GetPlaneAndTPC(endHit,plane,cstat,tpc,wire);
  wire_end=wire;
  time_end=endHit->PeakTime();
  
 
//  CalculateAxisParameters(nClust,hitlist,wire_start,time_start,wire_end,time_end);
}



//---------------------------------------------------------------------------------------------
// Find the trunk of the cluster and use it to determine new rought start/end points.
// uses, rough strarting points found earlier and the slope of the line found earlier
// returns preliminary shower direction //interface that executes all the stuff before it
//--------------------------------------------------------------------------------------------
int cluster::ClusterParamsAlg::FindTrunk(std::vector < art::Ptr < recob::Hit> > hitlist,double &wstn,double &tstn,double &wendn,double &tendn)
{
 double lineslope, lineintercept,goodness;
  
 Find2DAxisRough(lineslope,lineintercept,goodness,hitlist); 
 Find2DStartPointsBasic(lineslope,lineintercept, hitlist,wstn,tstn,wendn,tendn); 
 return FindTrunk(hitlist,wstn,tstn,wendn,tendn,lineslope,lineintercept); 
  
}



//---------------------------------------------------------------------------------------------
// Find the trunk of the cluster and use it to determine new rought start/end points.
// uses, rough strarting points found earlier and the slope of the line found earlier
// returns preliminary shower direction //interface that needs earlier input
//--------------------------------------------------------------------------------------------

int cluster::ClusterParamsAlg::FindTrunk(std::vector < art::Ptr < recob::Hit> > hitlist,
					 double &wstn,
					 double &tstn,
					 double &wendn,
					 double &tendn,
					 double lineslope,
					 double lineintercept)
  {
  int fDirection=0;
  unsigned int currplane=999;
  double fTotalCharge;
  
  
  if(hitlist.size()==0)
    return 0;

  double at=lineslope*fWireTimetoCmCm;
  double ct=lineintercept*fTimetoCm;

  //get slope of lines orthogonal to those found crossing the shower.
  double aprimt=0;
  if(at){
    aprimt=-1./at*fWireTimetoCmCm*fWireTimetoCmCm;
  }
  else
    aprimt=-999999.;
  
 
  // find extreme intercepts. For the time being we don't care which one is the start point and which one is the endpoint.
  
  double extreme_intercept_high,extreme_intercept_low;
    
  FindExtremeIntercepts(hitlist,aprimt,extreme_intercept_high,extreme_intercept_low);
    
  double extreme_intercept_start=-999999;
  double extreme_intercept_end=999999;
  
 
  if(at>=0) {  // for the time being assuming forward going showers
     extreme_intercept_end=extreme_intercept_high;
     extreme_intercept_start=extreme_intercept_low;
  }
  else if(at<0){
    extreme_intercept_start=extreme_intercept_high;
    extreme_intercept_end=extreme_intercept_low;
  }
  
//projections of the start and end points onto the rough cluster axis
  double wire_online_end,wire_online_begin,time_online_end,time_online_begin;

  gser.GetPointOnLineWSlopes(at,ct,extreme_intercept_end,wire_online_end,time_online_end);
  gser.GetPointOnLineWSlopes(at,ct,extreme_intercept_start,wire_online_begin,time_online_begin);
  

  double projectedlength=TMath::Sqrt( ((double)wire_online_end-(double)wire_online_begin)*((double)wire_online_end-(double)wire_online_begin)*fWiretoCm*fWiretoCm +  ((double)wire_online_end-(double)wire_online_begin)*((double)wire_online_end-(double)wire_online_begin)*at*at/fWireTimetoCmCm/fWireTimetoCmCm*fTimetoCm*fTimetoCm  );;

 //number of bins for a test histogram
  //int nbins = (hitlist.size() > 1000) ? 100 : hitlist.size()/10;

  hithistinv=new TH1F(Form("hithistinv_ev"),Form("hithistinv_pl"),(int)projectedlength,0.,projectedlength);  
  hithistinv->Clear();
  hithistinv->SetName(Form("hithistinv_pl"));
  
  for(std::vector < art::Ptr < recob::Hit > >::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
    art::Ptr<recob::Hit> theHit = (*hitIter);
    double time = theHit->PeakTime() ;  
    unsigned int wire,cstat, tpc,plane;
    GetPlaneAndTPC(theHit,plane,cstat,tpc, wire);
    
    if(currplane==999) currplane=plane;
    
    double wire_on_line,time_on_line;
    
    gser.GetPointOnLine(at/fWireTimetoCmCm,wire_online_begin,time_online_begin,wire,time,wire_on_line,time_on_line);
    
  //get linear distance from beginning and orthogonal distance from cluster axis
    double linedist=gser.Get2DDistance(wire_on_line,time_on_line,wire_online_begin,time_online_begin);
    double ortdist=gser.Get2DDistance(wire_on_line,time_on_line,wire,time);
    
   ////////////////////////////////////////////////////////////////////// 
   //calculate the weight along the axis, to guess the Shower direction  
   /////////////////////////////////////////////////////////////////////// 
    double weight= (ortdist<1.) ? 10*theHit->Charge() : 5*theHit->Charge()/ortdist;
    hithistinv->Fill(linedist,weight);
    
    fTotalCharge+=(*hitIter)->Charge();  //calculate total charge while we're at it.
        
  }
 
  //Now we have a histogram with the projection of the charge along the axis.
  
    double start_off=0., end_off=((double)wire_online_end-(double)wire_online_begin); 
    
    double fMaxweight=hithistinv->GetMaximum();  
  
   // calculate integrals along the axis: total, from maximum to start and from maximum to end 
   // double fullinteg=hithistinv->Integral();
    double startinteg=hithistinv->Integral(0,hithistinv->GetMaximumBin());
    double endinteg=hithistinv->Integral(hithistinv->GetMaximumBin(),hithistinv->GetNbinsX());
  
    
    if(fMaxweight>fChargeCutoffThreshold[currplane]) //check whether this operation makes sense
    {
      for(int x= hithistinv->GetMaximumBin();x>0;x--)
	{
	  if( hithistinv->GetBinContent(x)<fChargeCutoffThreshold[currplane] || x==1)
	  {
	    start_off=x;
	    if(  x>1 && hithistinv->Integral(0,x)/startinteg <0.01  ){
	      break;
	  }
	}
      }
 
      for(int x=hithistinv->GetMaximumBin();x<hithistinv->GetNbinsX();x++)
	{
	if( hithistinv->GetBinContent(x)<fChargeCutoffThreshold[currplane] || x==hithistinv->GetNbinsX()-1)
	  {
	  end_off=x;
	  if(  x<hithistinv->GetNbinsX()-1  && hithistinv->Integral(x,hithistinv->GetNbinsX())/endinteg <0.01  ){
	  break;
	  }
	}
      }
    }	
    
    
    
 //now create new rough points on axis based on the findings of the histogram above:
  double w_on_end,w_on_begin,t_on_end,t_on_begin;

  //calculate intervals in wire and time
  double dws= start_off*fabs(((double)wire_online_begin-(double)wire_online_end))/projectedlength;
  double dwe= end_off*fabs(((double)wire_online_begin-(double)wire_online_end))/projectedlength;

  double dts=at*dws/fWireTimetoCmCm;
  double dte=at*dwe/fWireTimetoCmCm;

  w_on_begin=wire_online_begin+dws;
  w_on_end=wire_online_begin+dwe;
  t_on_begin=time_online_begin+dts;
  t_on_end=time_online_begin+dte;

  //std::cout << "inside findTrunk, start points: " << w_on_begin << " "<<t_on_begin << " | "<< w_on_end << "  " << t_on_end << std::endl;
  /////// create new, two bin histogram, to test where the front and back are:
   double lowbin=0,highbin=0;
   double relowbin=0,rehighbin=0;
   double side_weight_start_charge=0,side_weight_end_charge=0;
   
  FindDirectionWeights(lineslope,w_on_begin, t_on_begin,w_on_end,t_on_end, hitlist,highbin,lowbin,rehighbin,relowbin); 
  
  FindSideWeights(lineslope,lineintercept,w_on_begin,t_on_begin,1,hitlist,side_weight_start_charge);
  FindSideWeights(lineslope,lineintercept,w_on_end,t_on_end,-1,hitlist,side_weight_end_charge);

  

  
  ///////////////////////////////////////////////////////
  // At this point the start still means most leftwise, end rightwise.
  //////////////////////////////////////////////////////
  
 
//   std::cout << "==== low and highbin " << lowbin << " " << highbin << std::endl;
//   std::cout << "==== reverse low and highbin " << relowbin << " " << rehighbin << std::endl;
//   std::cout << "==== reverse low and highbin hist" << hitreinv2->GetBinContent(1)<< " " << hitreinv2->GetBinContent(2) << std::endl;
//   
//   std::cout << " histo coords: LowE " <<  hitreinv2->GetBinLowEdge(1) << " " <<  hitreinv2->GetBinLowEdge(2) << " " << hitreinv2->GetBinWidth(1) << " " <<  hitreinv2->GetBinWidth(2) << std::endl;
  
  /////////Refine start points to real hits:
  double wstartnorm=0,tstartnorm=0,wendnorm=0,tendnorm=0;
  double wstartwght=0,tstartwght=0,wendwght=0,tendwght=0;
  double mindiststart=999999,mindistend=999999,minwghtdstart=999999,minwghtdend=999999;
  double inv_intercept_start=t_on_begin*fTimetoCm-aprimt*(double)w_on_begin*fWiretoCm;
  double inv_intercept_end=t_on_end*fTimetoCm-aprimt*(double)w_on_end*fWiretoCm; 
  
  for(std::vector < art::Ptr < recob::Hit > >::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
    art::Ptr<recob::Hit> theHit = (*hitIter);
    double time = theHit->PeakTime() ;  
    unsigned int wire,cstat, tpc,plane;
    GetPlaneAndTPC(theHit,plane,cstat,tpc, wire);
    
    double inv_intercept=time*fTimetoCm-aprimt*(double)wire*fWiretoCm;
    
    // we're basically only interested in the sign of these.
    double inter_distance_start = (inv_intercept - inv_intercept_start)*at;   // the a makes the sign work for negative.
    double inter_distance_end = (inv_intercept_end - inv_intercept)*at;   // the a makes the sign work for negative.
  
    if(inter_distance_start<0 || inter_distance_end<0)
      continue;
    
    
    double linedistst=gser.Get2DDistance(wire,time,w_on_begin,t_on_begin);
    double linedistend=gser.Get2DDistance(wire,time,w_on_end,t_on_end);
    double intercept=time*fTimetoCm-at*(double)wire*fWiretoCm;
    double locweight=(intercept-ct);
   
    //start points
    
    if(linedistst<mindiststart)
      {
      wstartnorm=wire; tstartnorm=time;
      mindiststart=linedistst;
      }
      // ony use the weight if it is big enough - otherwise you risk going the wrong way.
    if(linedistst<minwghtdstart && fabs(side_weight_start_charge/fTotalCharge) > 0.003 && side_weight_start_charge*locweight >=0 )
      {
      wstartwght=wire; tstartwght=time;
      minwghtdstart=linedistst;
      }
    
     //end points 
    
    if(linedistend<mindistend)
      {
      wendnorm=wire; tendnorm=time;
      mindistend=linedistend;
      }
     // ony use the weight if it is big enough - otherwise you risk going the wrong way.
   if(linedistend<minwghtdend &&  fabs(side_weight_start_charge/fTotalCharge) > 0.003 && side_weight_end_charge*locweight >=0 )
      {
      wendwght=wire; tendwght=time;
      minwghtdend=linedistend;
      }
 } // end refining of starting points
  
  if(wstartwght==0)
    {
    wstartwght=wstartnorm;
    tstartwght=tstartnorm;
    }
  if(wendwght==0)
    {
    wendwght=wendnorm;
    tendwght=tendnorm;
    }
  
  std::cout << "**** starting points. standard: " << wstartnorm << "," << tstartnorm << " weighted: " << wstartwght <<","<<tstartwght << std::endl;
  std::cout << "**** ending points. standard: " << wendnorm << "," << tendnorm << " weighted: " << wendwght <<","<<tendwght << std::endl;
  
  if(highbin>lowbin)  // cluster is leftward going
  {
    wstn=wendwght; 
    tstn=tendwght;
    wendn=wstartwght;
    tendn=tstartwght;
  
    fDirection=-1;
  }
  else  // cluster is rightward going 
  {
   wstn=wstartwght; 
   tstn=tstartwght;
   wendn=wendwght;
   tendn=tendwght; 
   fDirection=1;
  }
    
  std::cout << "**** Final starting points. standard: " << wstn << "," << tstn << std::endl;
    

return fDirection;

}





void cluster::ClusterParamsAlg::FindDirectionWeights(double lineslope,double w_on_begin,double t_on_begin,double w_on_end,double t_on_end, std::vector < art::Ptr < recob::Hit> > hitlist,double &HiBin,double &LowBin,double &invHiBin,double &invLowBin)
{
  
  double at=lineslope*fWireTimetoCmCm;
  //double ct=lineintercept*fTimetoCm;
  double smallprojlength= TMath::Sqrt( ((double)w_on_end-(double)w_on_begin)*((double)w_on_end-(double)w_on_begin)*fWiretoCm*fWiretoCm +  ((double)w_on_end-(double)w_on_begin)*((double)w_on_end-(double)w_on_begin)*at*at/fWireTimetoCmCm/fWireTimetoCmCm*fTimetoCm*fTimetoCm  );
  
  //std::cout << " smallprojlength: "<< smallprojlength << std::endl;
  
  //get slope of lines orthogonal to those found crossing the shower.
  double aprimt=0;
  if(at){
    aprimt=-1./at*fWireTimetoCmCm*fWireTimetoCmCm;
  }
  else
    aprimt=-999999.;
  
   //histogram for weight, i.e. hits along the axis:
  hitinv2=new TH1F(Form("hitinv2_ev"),Form("hitinv2_ev"),2,0.,smallprojlength);  
  double lowbin=0,highbin=0;

  //histogram for inverse weight:
  hitreinv2=new TH1F(Form("rehitinv2_ev"),Form("rehitinv2_ev"),2,0.,smallprojlength);  
  double relowbin=0,rehighbin=0;

  hitinv2->Clear();
  hitinv2->SetName(Form("hitinv2_ev"));
  hitreinv2->Clear();
  hitreinv2->SetName(Form("rehitinv2_"));
  
  ///////////////////////////////////////////////////////
  // At this point the start still means most leftwise, end rightwise.
  //////////////////////////////////////////////////////
  
  
   //// these are the intercepts of perpendicular lines going through the start and end point.
   //// I will use the perpendicular intercept as a measurement of distance for simplicity.
    
   double inv_intercept_start=t_on_begin*fTimetoCm-aprimt*(double)w_on_begin*fWiretoCm;
   double inv_intercept_end=t_on_end*fTimetoCm-aprimt*(double)w_on_end*fWiretoCm;
    
  // std::cout << "intercepts in weights: " << inv_intercept_start << " " << inv_intercept_end << " " << hitlist.size()  <<std::endl;
 
   ///////// Calculate weights
    for(std::vector < art::Ptr < recob::Hit > >::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
      art::Ptr<recob::Hit> theHit = (*hitIter);
      double time = theHit->PeakTime() ;  
      unsigned int wire,cstat, tpc,plane;
      GetPlaneAndTPC(theHit,plane,cstat,tpc, wire);
    
      double wire_on_line,time_on_line;
    
      gser.GetPointOnLine(at/fWireTimetoCmCm,w_on_begin,t_on_begin,wire,time,wire_on_line,time_on_line);
      double inv_intercept=time*fTimetoCm-aprimt*(double)wire*fWiretoCm;
    
      // we're basically only interested in the sign of these.
      double inter_distance_start = (inv_intercept - inv_intercept_start)*at;   // the a makes the sign work for negative.
      double inter_distance_end = (inv_intercept_end - inv_intercept)*at;   // the a makes the sign work for negative.
  
  
      // startEnd checking
  
      double linedist=gser.Get2DDistance(wire_on_line,time_on_line,w_on_begin,t_on_begin);
      double ortdist=gser.Get2DDistance(wire_on_line,time_on_line,wire,time);
      double weight= (ortdist<1.) ? 10*theHit->Charge() : 5*theHit->Charge()/ortdist;
      double revweight= (ortdist<1.) ? 0.1 : ortdist;

      hitreinv2->Fill(linedist,revweight);
      hitinv2->Fill(linedist,weight);
 
      if(inter_distance_start>0. && linedist<= smallprojlength*0.5 )
      {
      lowbin+=weight;
      relowbin+=revweight;
      }
      else if (linedist>smallprojlength*0.5 && inter_distance_end > 0. ) 
      { highbin+=weight;
	rehighbin+=revweight;
      }
     
 
       
  }  // end weight calculating loop
  
  
  HiBin=highbin;
  LowBin=lowbin;
  invHiBin=rehighbin;
  invLowBin=relowbin;
  
  
}





void cluster::ClusterParamsAlg::FindSideWeights(double lineslope,double lineintercept,double w_on_begin,double t_on_begin, int direction,std::vector < art::Ptr < recob::Hit> > hitlist,double &side_weight_start_charge)
{
 
 
  double at=lineslope*fWireTimetoCmCm;
  double ct=lineintercept*fTimetoCm;
   
  //get slope of lines orthogonal to those found crossing the shower.
  double aprimt=0;
  if(at){
    aprimt=-1./at*fWireTimetoCmCm*fWireTimetoCmCm;
  }
  else
    aprimt=-999999.;
   

   //// these are the intercepts of perpendicular lines going through the start and end point.
   //// I will use the perpendicular intercept as a measurement of distance for simplicity.
    
   double inv_intercept_start=t_on_begin*fTimetoCm-aprimt*(double)w_on_begin*fWiretoCm;
    
   int startctr=0;

   ///////// Calculate weights
    for(std::vector < art::Ptr < recob::Hit > >::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
      art::Ptr<recob::Hit> theHit = (*hitIter);
      double time = theHit->PeakTime() ;  
      unsigned int wire,cstat, tpc,plane;
      GetPlaneAndTPC(theHit,plane,cstat,tpc, wire);
    
      double wire_on_line,time_on_line;
    
      gser.GetPointOnLine(at/fWireTimetoCmCm,w_on_begin,t_on_begin,wire,time,wire_on_line,time_on_line);
      double inv_intercept=time*fTimetoCm-aprimt*(double)wire*fWiretoCm;
    
      // we're basically only interested in the sign of these.
      double inter_distance_start = (inv_intercept - inv_intercept_start)*at*direction;   // the a makes the sign work for negative., the direction, differentiates between points at start and end.
      
      double linedist=gser.Get2DDistance(wire_on_line,time_on_line,w_on_begin,t_on_begin);
    
 
           
     //now calculate the side weight, i.e. on which side of the axis more charge is found (on the start side)
      if( linedist < 10. && inter_distance_start > 0.  )
	{
	// this is the intercept the line would have if it went through the point
	double intercept=time*fTimetoCm-at*(double)wire*fWiretoCm;
    
	side_weight_start_charge+=(intercept-ct)*theHit->Charge();
	startctr++;
     	}
    
     
     
       
  }  // end weight calculating loop
 //////////////////////////////////////////////////
   if(startctr>0)
    {
      side_weight_start_charge/=startctr;
    }
  
        
}





//---------------------------------------------------------------------------------------------
//Create new Hitlist consisting of hits higher than a Threshold and not further away than 
// radius: 
// and the list not smaller then 40: hits:
//--------------------------------------------------------------------------------------------


std::vector< art::Ptr<recob::Hit> > cluster::ClusterParamsAlg::CreateHighHitlist(double threshold,std::vector< art::Ptr<recob::Hit> > hitlist)
{
  
  std::vector< art::Ptr<recob::Hit> > hitlist_high;
    
  //first select a subset of points that is close to the selected start point:
  std::vector < art::Ptr<recob::Hit> > hitlistlocal;
  
  for(unsigned int ix = 0; ix<  hitlist.size();  ix++){
    	art::Ptr<recob::Hit> theHit = hitlist[ix];
	if(theHit->Charge()>threshold)
	  hitlist_high.push_back(theHit);
    }
  
  //subtracting outliers
    
  // this is the erasing loop
 for(unsigned int ix = 0; ix<  hitlist_high.size();  ix++){
   if(hitlist_high.size()<(unsigned int)fMinHitListSize)
	  break;
   art::Ptr<recob::Hit> theHit = hitlist_high[ix];
   double time = theHit->PeakTime() ;  
   unsigned int plane,cstat,tpc,wire;
   GetPlaneAndTPC(theHit,plane,cstat,tpc,wire);
	
   //ok. now I have a sub selection of hits that are close to the perceived start of the shower.
   SelectLocalHitlist(hitlist_high, hitlistlocal, (double)wire,time,fOutlierRadius);
   
  
   if(hitlistlocal.size()<10 && hitlist_high.size() > ix && hitlist_high.size()>0 ){
      hitlist_high.erase(hitlist_high.begin()+ix);
      ix--;
    //  std::cout << " erasing hit @ w,t" << wire << " "<< time << " (start) "  << std::endl;
    }
   hitlistlocal.clear();
	

 } 
  
return hitlist_high;  
}
/////////////////////////////////////////////////////////////////////////
// refine the provided start points by looking for Hough Lines at the start and end of 
// cluster.
//////////////////////////////////////////////////////////////////////////

void cluster::ClusterParamsAlg::RefineStartPointsHough(std::vector< art::Ptr < recob::Hit> > hitlist, double & wire_start,double & time_start, double & wire_end,double & time_end, int &direction)
{
  //parameters to be set later?
  double linearlimit=fSelectBoxSizePar; 
  double ortlimit=fSelectBoxSizePerp;
  int direction_local=1;
  
  std::vector < art::Ptr<recob::Hit> > hitlistlocal_start;
  std::vector < art::Ptr<recob::Hit> > hitlistlocal_end;
  std::vector < art::Ptr<recob::Hit> > hitlistlocal_refined;  //these in principle should only be hits on the line.
  
  ///////////////////////////// find trunk of shower:
  //double wstn=0.,tstn=0.,wendn=0.,tendn=0.;
  
  //FindTrunk(nClust,hitlist,wstn,tstn,wendn,tendn);
  
  /////
  double lineslope=0; //either add as parameter or pick up from Find Rough
  
  //first select a subset of points that is close to the selected start point:
  
 SelectLocalHitlist(hitlist, hitlistlocal_start, wire_start,time_start,linearlimit,ortlimit,lineslope);
 SelectLocalHitlist(hitlist, hitlistlocal_end, wire_end,time_end,linearlimit,ortlimit,lineslope);
 std::cout << "Hough hitlist_local start size " <<  hitlistlocal_start.size() << std::endl;
 std::cout << "Hough hitlist_local end size " <<  hitlistlocal_end.size() << std::endl;
 std::vector < art::PtrVector<recob::Hit> > houghlines; 
 std::vector < art::PtrVector<recob::Hit> > houghlinesend; 
 
 //save sizes of hitlists found
 unsigned int lochitlistsize=hitlistlocal_start.size();
 unsigned int lochitlistendsize=hitlistlocal_end.size();
   
   
//  for(unsigned int ix=0;ix<hitlistlocal_start.size();ix++) {
//    unsigned int w;
//    //unsigned int c = (*hitlistlocal_start[ix]).Wire()->RawDigit()->Channel(); 
//    //geo->ChannelToWire(c,cs,t,p,w);
//    w  = (*hitlistlocal_start[ix]).WireID().Wire;
//   //  std::cout << " local hits for: " << " plane: " << p << " w,t: "  << w<< " "<< (*hitlistlocal_start[ix]).PeakTime() << std::endl;
//    }
  // std::cout<<std::endl<<std::endl;    
    
  unsigned int iplane,cstat,tpc,wire;
  GetPlaneAndTPC((*hitlist.begin()),iplane,cstat,tpc,wire);
// 
    
  size_t numclus = fHBAlg.FastTransform(hitlistlocal_start, houghlines);
  size_t numclusend = fHBAlg.FastTransform(hitlistlocal_end, houghlinesend);
 
 
  std::cout << "found " << numclus << "lines with HoughBaseAlg ";
  if(numclus>0)
   std::cout << houghlines[0].size();
  std::cout << std::endl;
  std::cout << "found " << numclusend << "end lines with HoughBaseAlg ";
  if(numclusend>0)
   std::cout << houghlinesend[0].size();
  std::cout << std::endl; 

  //int nhoughlines=numclus;
  //int nhoughlinesend=numclusend;
 
  int starthoughclus=-1,endhoughclus=-1;  //index of largest houghline
 // int numclusfunc;
 // std::vector < art::PtrVector<recob::Hit> > *houghutil; 
 // numclusfunc=numclus;
 // houghutil=&houghlines;
 
 // assume that direction inserted is correct (wire_start is start), unless found otherwise
 //if(numclusend>=1 && ( (fReLowBin > fReHighBin && fDirection==1 ) || (fReLowBin < fReHighBin && fDirection==-1 ) ) )
  if(numclusend>=1 )
  {
   unsigned int stsize=0; 
   unsigned int endsize=0;
   for(unsigned int i=0;i<numclus;i++)
   {
     if(houghlines[i].size()>stsize)
     {
//       TGraph tt=TGraph(houghlines[i].size());
//       for(int xx=0;xx<houghlines[i].size();xx++) 
// 	{
// 	unsigned int c = (*houghlines[i][xx]).Wire()->RawDigit()->Channel(); 
// 	unsigned int cs,t,p,w;
// 	geom->ChannelToWire(c,cs,t,p,w);
// 	tt.SetPoint(xx,w*fWiretoCm,(*houghlines[i][xx]).PeakTime()*fTimetoCm);
// 	}
// 	
//       tt.Fit("pol1");	
//       double na=tt.GetFunction("pol1")->GetParameter(1)/fWireTimetoCmCm;	 
//       double nc=tt.GetFunction("pol1")->GetParameter(0)/fTimetoCm;	 
//       
//      // take only lines that are at angle roughly similar to the original line 
//       if(fabs(TMath::ATan(na) - TMath::ATan(lineslopetest) ) < 30*TMath::Pi()/180. ){
	stsize= houghlines[i].size();
	starthoughclus=i; 
// 	std::cout << " saving stsize: " << stsize << std::endl;
//       }
     }
   }
   /// end sizes:
    for(unsigned int i=0;i<numclusend;i++)
   {
     if(houghlinesend[i].size()>endsize)
     {
//       TGraph tt=TGraph(houghlinesend[i].size());
//       for(int xx=0;xx<houghlinesend[i].size();xx++) 
// 	{
// 	unsigned int c = (*houghlinesend[i][xx]).Wire()->RawDigit()->Channel(); 
// 	unsigned int cs,t,p,w;
// 	geom->ChannelToWire(c,cs,t,p,w);
//      // tt.SetPoint(xx,w*fWiretoCm,(*houghlinesend[i][xx]).PeakTime()*fTimetoCm);
//       // std::cout << " hough hits for set: "<< i << " plane: " << p << " w,t: "  << w<< " "<< (*houghlinesend[i][xx]).PeakTime() << std::endl;
//       
// 	}
	
//       tt.Fit("pol1");	
//       double na=tt.GetFunction("pol1")->GetParameter(1)/fWireTimetoCmCm;	 
//       double nc=tt.GetFunction("pol1")->GetParameter(0)/fTimetoCm;	 
//       
//       std::cout << "local line fit, end " << na << ","<<nc << " global " << lineslopetest << " " << lineinterctest << std::endl;
//       std::cout << " Arctans: " << TMath::ATan(na) << " " << TMath::ATan(lineslopetest) << std::endl;
//       
//       if(fabs(TMath::ATan(na) - TMath::ATan(lineslopetest) ) < 30*TMath::Pi()/180. ){
	endsize= houghlinesend[i].size();
	endhoughclus=i;
// 	std::cout << " saving endsize: " << endsize << std::endl;
//       }
     }
   }
   
//    for(int i=0;i<numclusend;i++)
//    endsize= houghlinesend[i].size();
  
   
   
   //if(endsize>stsize )
   if(endsize/lochitlistendsize  > stsize/lochitlistsize  )
   {
     direction_local=-1; 
     std::cout << "!!! reversing direction!!! " << direction_local << std::endl;
   //  numclusfunc=numclusend;
   //  houghutil=&houghlinesend; 
   }
   
  }
   
 
 
 
 //std::cout << "new parameters: " << numclusfunc << " " << (*houghutil)[0].size() << std::endl;
 
  double startwire=10000*(1+direction),starttime=-1,endwire=10000*(1-direction),endtime=-1;
 
  //use longest houghline to find new start
  if(endhoughclus!=-1)
  {
    for(unsigned int ix=0;ix<houghlinesend[endhoughclus].size();ix++) {
	unsigned int w;
	//unsigned int c = houghlinesend[endhoughclus][ix]->Wire()->RawDigit()->Channel(); 
	// geo->ChannelToWire(c,cs,t,p,w);
	//GetPlaneAndTPC((*hitlist.begin()),iplane,cstat,tpc,wire);
	w=houghlinesend[endhoughclus][ix]->WireID().Wire; 
	//p=houghlinesend[endhoughclus][ix]->WireID().Plane; 
   
      
      
	/*    std::cout << " hough hits for set: "<< i << " plane: " << p << " w,t: "  << w<< " "<< (*houghlines[i][ix]).PeakTime() << std::endl;
	*/
	//std::cout << " hough hits for set: "<< endhoughclus << " plane: " << p << " w,t: "  << w<< " "<< houghlinesend[endhoughclus][ix]->PeakTime() << std::endl;
   
	if(w*direction>endwire*direction) //multiply by original direction to find most extreme point 
	  {endwire=w;
	  endtime=houghlinesend[endhoughclus][ix]->PeakTime();
	  }
      }
  }
  
  if(starthoughclus!=-1)
  {
    for(unsigned int ix=0;ix<houghlines[starthoughclus].size();ix++) {
	unsigned int w;
	//  unsigned int c = houghlines[starthoughclus][ix]->Wire()->RawDigit()->Channel(); 
	//  geo->ChannelToWire(c,cs,t,p,w);
    
	w=houghlines[starthoughclus][ix]->WireID().Wire; 
	
	/*    std::cout << " hough hits for set: "<< i << " plane: " << p << " w,t: "  << w<< " "<< (*houghlines[i][ix]).PeakTime() << std::endl;  */    
	//std::cout << " hough hits for set: "<<starthoughclus << " plane: " << p << " w,t: "  << w<< " "<< houghlines[starthoughclus][ix]->PeakTime() << std::endl;
   
	if(w*direction<startwire*direction) //multiply by original direction to find most extreme point 
	  {startwire=w;
	  starttime=houghlines[starthoughclus][ix]->PeakTime();
	  }
    }
  }
  
  
 std::cout << " low, high wire,time: " << startwire <<","<<starttime << " | " << endwire<<","<<endtime<<std::endl;
 direction=direction_local;  // returning the newfound direction for a decision outside.
 
}




////////////////////// add override if forcing rightwise:

int cluster::ClusterParamsAlg::DecideClusterDirection(std::vector < art::Ptr<recob::Hit> > hitlist,
				double lineslope,double &wstn,double &tstn,double &wendn,double &tendn)  
{
  int fDirection=1;
  double HiBin,LowBin,invHiBin,invLowBin;
  double wire_start,time_start,wire_end,time_end;
  
  wire_start=wstn;
  wire_end=wendn;
  time_start=tstn;
  time_end=tendn;
  
  //doesn't change wire,time values - assumes it's already the right ones.
  FindDirectionWeights(lineslope,wire_start,time_start,wire_end,time_end,hitlist,HiBin,LowBin,invHiBin,invLowBin); 
  
  
 if(HiBin>LowBin)  // shower is leftward going (contrary to our previous assumptions)
  {
  wire_start=wendn; 
  time_start=tendn;
  wire_end=wstn;
  time_end=tstn;
  fDirection=-1;
  }
  
  //this one does refine the start and end point, although it leaves them in the same order.
  // it will now be up to the algorithm here to decide whether we should override the decision  
  // based on FindDirectionWeights
  int houghdirection=fDirection;
  RefineStartPointsHough(hitlist,wire_start,time_start,wire_end,time_end,houghdirection);
  
  
  //check whether it makes sense to reverse direction based on the findings of RefineStartPointsHough
   if(fDirection!=houghdirection 
     && ( (invLowBin > invHiBin && fDirection==1 ) || (invLowBin < invHiBin && fDirection==-1 ) ))
   {
    //overriding basic with HoughLine Direction, flipping start and end 
        wstn=wire_end;
	tstn=time_end;
	wendn=wire_start;
	tendn=time_start;
     
     fDirection = houghdirection; 
     
   }
  else   // not flipping
   {  // assigning to pick up the hough refined values
      wstn=wire_start;
      tstn=time_start;
      wendn=wire_end;
      tendn=time_end;
   }
    
  
  if(fForceRightGoing)
  {
   if(wstn > wendn  || fDirection==-1 )  //change only if decided it's backwards up to now. otherwise it's fine as is.
   {
     std::cout << "forcing right going shower" << std::endl;
     fDirection=1;
   
     wstn=wire_end;
     tstn=time_end;
      
     wendn=wire_start;
     tendn=time_start;
   }
   
  }
 

 
 return fDirection;
}




void cluster::ClusterParamsAlg::FindExtremeIntercepts(std::vector < art::Ptr<recob::Hit> > hitlist,
			     double perpslope,
			     double &inter_high,
			     double &inter_low)  
{
  
  inter_high=-999999;
  inter_low=999999;

  unsigned int plane,tpc,wire,cstat;
  
  for(std::vector < art::Ptr < recob::Hit > >::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
    art::Ptr<recob::Hit> theHit = (*hitIter);
    double time = theHit->PeakTime() ;  
    GetPlaneAndTPC(theHit,plane,cstat,tpc,wire);
    
    //wire_bar+=wire;
    //time_bar+=time;	
    //nhits++;
    
    double intercept=time*fTimetoCm-perpslope*(double)wire*fWiretoCm;
    
    if(intercept > inter_high ){
      inter_high=intercept;
    
    }
    if(intercept < inter_low ){
      inter_low=intercept;
    
    }  

    

  }   // end of first HitIter loop, at this point we should have the extreme intercepts 

}




void cluster::ClusterParamsAlg::SelectLocalHitlist(std::vector< art::Ptr < recob::Hit> > hitlist, std::vector < art::Ptr<recob::Hit> > &hitlistlocal, double  wire_start,double time_start, double linearlimit,   double ortlimit, double lineslopetest)
{
  
  double locintercept=time_start-wire_start*lineslopetest;
  
  //std::cout << "$$$$$$$$$$$$$ locintercept: " << locintercept << " " << lineinterctest << std::endl;
  
  for(std::vector < art::Ptr < recob::Hit > >::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
    	art::Ptr<recob::Hit> theHit = (*hitIter);
    	double time = theHit->PeakTime() ;  
    	unsigned int plane,cstat,tpc,wire;
	GetPlaneAndTPC(theHit,plane,cstat,tpc,wire);
	
	double wonline=wire,tonline=time;
	//gser.GetPointOnLine(lineslopetest,lineinterctest,wire,time,wonline,tonline);
	gser.GetPointOnLine(lineslopetest,locintercept,wire,time,wonline,tonline);
	
	//calculate linear distance from start point and orthogonal distance from axis
	double lindist=gser.Get2DDistance(wonline,tonline,wire_start,time_start);
	double ortdist=gser.Get2DDistance(wire,time,wonline,tonline);
	
	//std::cout << " w,t: " << wire << " " << time << " ws,ts " << wonline << " "<< tonline <<" "<< lindist << " " << ortdist << std::endl;
	
	if(lindist<linearlimit && ortdist<ortlimit)
	  hitlistlocal.push_back(theHit);
    
    
    }
    
}
///////////////////

void cluster::ClusterParamsAlg::SelectLocalHitlist(std::vector< art::Ptr < recob::Hit> > hitlist, std::vector < art::Ptr<recob::Hit> > &hitlistlocal, double  wire_start,double time_start, double radlimit)
{
  for(std::vector < art::Ptr < recob::Hit > >::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
    	art::Ptr<recob::Hit> theHit = (*hitIter);
    	double time = theHit->PeakTime() ;  
    	unsigned int plane,cstat,tpc,wire;
	GetPlaneAndTPC(theHit,plane,cstat,tpc,wire);
	
	//calculate linear distance from start point and orthogonal distance from axis
	double lindist=gser.Get2DDistance(wire,time,wire_start,time_start);
	
	if(lindist<radlimit )
	  hitlistlocal.push_back(theHit);
    
    
    }
    
}

/*
void cluster::ClusterParamsAlg::FindExtremeIntercepts(std::vector < art::Ptr<recob::Hit> > hitlist,
			     double perpslope,
			     double &inter_high,
			     double &inter_low)  
{
  
  inter_high=-999999;
  inter_low=999999;

  unsigned int plane,tpc,wire,cstat;
  
  for(std::vector < art::Ptr < recob::Hit > >::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
    art::Ptr<recob::Hit> theHit = (*hitIter);
    double time = theHit->PeakTime() ;  
    GetPlaneAndTPC(theHit,plane,cstat,tpc,wire);
    
    //wire_bar+=wire;
    //time_bar+=time;	
    //nhits++;
    
    double intercept=time*fTimetoCm-perpslope*(double)wire*fWiretoCm;
    
    if(intercept > inter_high ){
      inter_high=intercept;
    
    }
    if(intercept < inter_low ){
      inter_low=intercept;
    
    }  

    

  }   // end of first HitIter loop, at this point we should have the extreme intercepts 

}*/




art::Ptr<recob::Hit> cluster::ClusterParamsAlg::FindClosestHit(std::vector < art::Ptr < recob::Hit > > hitlist,
			     double wire_online,
			     double time_online)
{
  
  double min_length_from_start=99999;
  art::Ptr<recob::Hit> nearHit;
   
  unsigned int plane,tpc,wire,cstat;
   
   
  for(std::vector < art::Ptr < recob::Hit > >::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
    art::Ptr<recob::Hit> theHit = (*hitIter);
    double time = theHit->PeakTime() ;  
    GetPlaneAndTPC(theHit,plane,cstat,tpc,wire);
  //  std::cout << " in find closest hit " << std::endl;
    double dist_mod=gser.Get2DDistance(wire_online,time_online,wire,time);
    //TMath::Sqrt( pow(((double)wire_online-(double)wire*fWirePitch),2)+pow((time_online-time*fDriftVelocity*fTimeTick),2) );	

    if(dist_mod<min_length_from_start){
	//wire_start[plane]=wire;
	//time_start[plane]=time;
	nearHit=(*hitIter);
	min_length_from_start=dist_mod;
	}	

  } 
  
return nearHit;    
}









// ******************************* //
int cluster::ClusterParamsAlg::GetPlaneAndTPC(art::Ptr<recob::Hit> a,
						unsigned int &p,
						unsigned int &cs,
						unsigned int &t,
						unsigned int &w)
{
   
    
     w=a->WireID().Wire; 
     p=a->WireID().Plane; 
     t=a->WireID().TPC; 
     cs=a->WireID().Cryostat; 
     
  return 0;
}




