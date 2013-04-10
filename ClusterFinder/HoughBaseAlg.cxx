/////////////////////////////////////////////////////////////////////
///
/// HoughBaseAlg class
///
/// Ben Carls, bcarls@fnal.gov
///
/// The Hough transform employed by fuzzy clustering is a heavily modified variant of the original 
/// Hough line code. It identifies lines in fuzzy clusters (e.g. muon tracks) and splits them off
/// into new clusters. 
///
/// The algorithm is based on the Progressive Probabilistic Hough Transform (PPHT).
/// See the following paper for details:
///
/// J. Matas et al., Robust Detection of Lines Using the Progressive Probabilistic Hough Transform,
/// Computer Vision and Image Understanding, Volume 78, Issue 1, April 2000, Pages 119–137
///
////////////////////////////////////////////////////////////////////////

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

#include <sstream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <vector>
#include <stdint.h>

#include <TF1.h>
#include <TH1D.h>

#include "CLHEP/Random/RandFlat.h"

#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
 
#include "Filters/ChannelFilter.h"
#include "ClusterFinder/HoughBaseAlg.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/AssociationUtil.h"

// Define parameters that will tell us if we are doing a normal Hough line merge
// or a shower Hough line merge
static const int iMergeShowerPart = 0;
static const int iMergeShower     = 0;
static const int iMergeNormal     = 1;

//------------------------------------------------------------------------------
cluster::HoughBaseAlg::HoughBaseAlg(fhicl::ParameterSet const& pset)
{
  this->reconfigure(pset);
}

//------------------------------------------------------------------------------
cluster::HoughBaseAlg::~HoughBaseAlg()
{
}

//------------------------------------------------------------------------------
void cluster::HoughBaseAlg::reconfigure(fhicl::ParameterSet const& pset)
{ 
  fMaxLines                       = pset.get< int    >("MaxLines"                      );
  fMinHits                        = pset.get< int    >("MinHits"                       );
  fSaveAccumulator                = pset.get< int    >("SaveAccumulator"               );
  fNumAngleCells                  = pset.get< int    >("NumAngleCells"                 );
  fMaxDistance                    = pset.get< double >("MaxDistance"                   );
  fMaxSlope                       = pset.get< double >("MaxSlope"                      );
  fRhoZeroOutRange                = pset.get< int    >("RhoZeroOutRange"               );
  fThetaZeroOutRange              = pset.get< int    >("ThetaZeroOutRange"             );
  fRhoResolutionFactor            = pset.get< int    >("RhoResolutionFactor"           );
  fPerCluster                     = pset.get< int    >("HitsPerCluster"                );
  fMissedHits                     = pset.get< int    >("MissedHits"                    );
  fDoHoughLineMerge               = pset.get< double >("DoHoughLineMerge"              );
  fHoughLineMergeAngle            = pset.get< double >("HoughLineMergeAngle"           );
  fHoughLineMergeCutoff           = pset.get< double >("HoughLineMergeCutoff"          );
  fDoShowerPartHoughLineMerge     = pset.get< double >("DoShowerPartHoughLineMerge"    );
  fShowerPartHoughLineMergeAngle  = pset.get< double >("ShowerPartHoughLineMergeAngle" );
  fShowerPartHoughLineMergeCutoff = pset.get< double >("ShowerPartHoughLineMergeCutoff");
  fDoShowerHoughLineMerge         = pset.get< double >("DoShowerHoughLineMerge"        );
  fShowerHoughLineMergeAngle      = pset.get< double >("ShowerHoughLineMergeAngle"     );
  fShowerHoughLineMergeCutoff     = pset.get< double >("ShowerHoughLineMergeCutoff"    );
  fLineIsolationCut               = pset.get< double >("LineIsolationCut"              );
  fChargeAsymmetryCut             = pset.get< double >("ChargeAsymmetryCut"            );
  fSigmaChargeAsymmetryCut        = pset.get< double >("SigmaChargeAsymmetryCut"       );
  fMergeShowerLikenessCut         = pset.get< double >("MergeShowerLikenessCut"        );
  return;
}

//------------------------------------------------------------------------------
cluster::HoughTransform::HoughTransform()
{  
}


//------------------------------------------------------------------------------
size_t cluster::HoughBaseAlg::Transform(std::vector<art::Ptr<recob::Hit> > const& hits,
					std::vector<unsigned int>                 *fpointId_to_clusterId,
					unsigned int                               clusterId, // The id of the cluster we are examining
					int                                       *nClusters,
					std::vector<unsigned int>                  corners
					)
{

  int nClustersTemp = *nClusters;
  
  HoughTransform c;
  std::vector< art::Ptr<recob::Hit> > hit;

  art::ServiceHandle<geo::Geometry> geom;
  art::ServiceHandle<util::LArProperties> larprop;
  art::ServiceHandle<util::DetectorProperties> detprop;
  filter::ChannelFilter chanFilt;

  uint32_t     channel = hits[0]->Wire()->RawDigit()->Channel();
  unsigned int wire    = 0;
  unsigned int wireMax = 0;
  std::vector<lineSlope> linesFound;
  linesFound.clear();


  //mf::LogInfo("HoughBaseAlg") << "nClusters is: " << *nClusters;


  //size_t cinctr = 0;
  //geo::View_t    view = geom->Cryostat(cs).TPC(t).Plane(p).View();
  geo::SigType_t sigt = geom->SignalType(channel);
  std::vector<int> skip;  
  hit.clear();
  
  //factor to make x and y scale the same units
  double xyScale  = .001*larprop->DriftVelocity(larprop->Efield(),larprop->Temperature());
  xyScale        *= detprop->SamplingRate()/geom->WirePitch(0, 1, 
							    hits[0]->WireID().Plane,
							    hits[0]->WireID().TPC,
							    hits[0]->WireID().Cryostat);

  const double pos[3] = {0., 0.0, 0.};
  double posWorld0[3] = {0.};
  double posWorld1[3] = {0.};
  const geo::WireGeo& wireGeom = geom->Plane(0).Wire(0);
  const geo::WireGeo& wire1Geom = geom->Plane(0).Wire(1);
  wireGeom.LocalToWorld(pos, posWorld0);
  wire1Geom.LocalToWorld(pos, posWorld1);
  double wire_dist = posWorld0[1]- posWorld1[1];
  double tickToDist = larprop->DriftVelocity(larprop->Efield(),larprop->Temperature());
  tickToDist *= 1.e-3 * detprop->SamplingRate(); // 1e-3 is conversion of 1/us to 1/ns


  mf::LogInfo("HoughBaseAlg") << "xyScale: " << xyScale;
  
  int x, y;
  //unsigned int channel, plane, wire, tpc, cstat;
  //there must be a better way to find which plane a cluster comes from
  int dx = geom->Cryostat(hits[0]->WireID().Cryostat).TPC(hits[0]->WireID().TPC).Plane(hits[0]->WireID().Plane).Nwires();//number of wires 
  int dy = hits[0]->Wire()->NSignal();//number of time samples. 
  skip.clear();
  skip.resize(hits.size());
  std::vector<int> listofxmax;
  std::vector<int> listofymax;  
  std::vector<int> hitsTemp;        //indecies ofcandidate hits
  std::vector<int> sequenceHolder; //channels of hits in list
  std::vector<int> currentHits;    //working vector of hits 
  std::vector<int> lastHits;       //best list of hits
  std::vector<art::Ptr<recob::Hit> > clusterHits;
  double indcolscaling = 0.;       //a parameter to account for the different 
        			   ////characteristic hit width of induction and collection plane
  double centerofmassx = 0;
  double centerofmassy = 0;
  double denom = 0; 
  double intercept=0.;
  double slope = 0.;
  //this array keeps track of the hits that have already been associated with a line. 
  int xMax = 0;
  int yMax = 0;
  double rho;
  double theta; 
  int accDx(0), accDy(0);





  /// Outline of PPHT, J. Matas et. al. 
  /// ---------------------------------------
  /// 
  ///LOOP over hits, picking a random one
  ///  Enter the point into the accumulator
  ///  IF it is already in the accumulator or part of a line, skip it
  ///  Store it in a vector of points that have been chosen
  ///
  ///  Find max value in accumulator; IF above threshold, create a line
  ///    Subtract points in line from accumulator
  ///    
  ///
  ///END LOOP over hits, picking a random one
  ///

  ///Init specifies the size of the two-dimensional accumulator 
  ///(based on the arguments, number of wires and number of time samples). 
  c.Init(dx,dy,fRhoResolutionFactor,fNumAngleCells);
  /// Adds all of the hits to the accumulator
  mf::LogInfo("HoughBaseAlg") << "Beginning PPHT";

  c.GetAccumSize(accDy, accDx);

  /// count is how many points are left to randomly insert
  unsigned int count = hits.size();
  std::vector<unsigned int> accumPoints;
  accumPoints.resize(hits.size());
  int nAccum = 0;

  /// Get the random number generator
  art::ServiceHandle<art::RandomNumberGenerator> rng;
  CLHEP::HepRandomEngine & engine = rng -> getEngine();
  CLHEP::RandFlat flat(engine);

  for( ; count > 0; count--){
  

    /// The random hit we are examining
    ///unsigned int randInd = rand() % hits.size();
    unsigned int randInd = (unsigned int)(flat.fire()*hits.size());

    /// If the point isn't in the current fuzzy cluster, skip it
    if(fpointId_to_clusterId->at(randInd) != clusterId)
      continue;

    /// Skip if it's already in a line
    if(skip[randInd]==1)
      continue;

    /// If we have already accumulated the point, skip it
    if(accumPoints[randInd])
      continue;
    accumPoints[randInd]=1;
   
    /// zeroes out the neighborhood of all previous lines  
    for(auto listofxmaxItr = listofxmax.begin(); listofxmaxItr != listofxmax.end(); ++listofxmaxItr) {
      int yClearStart = listofymax[listofxmaxItr-listofxmax.begin()] - fRhoZeroOutRange;
      if (yClearStart < 0) yClearStart = 0;
      
      int yClearEnd = listofymax[listofxmaxItr-listofxmax.begin()] + fRhoZeroOutRange;
      if (yClearEnd >= accDy) yClearEnd = accDy - 1;
      
      int xClearStart = *listofxmaxItr - fThetaZeroOutRange;
      if (xClearStart < 0) xClearStart = 0;
      
      int xClearEnd = *listofxmaxItr + fThetaZeroOutRange;
      if (xClearEnd >= accDx) xClearEnd = accDx - 1;
      
      for (y = yClearStart; y <= yClearEnd; ++y){
        for (x = xClearStart; x <= xClearEnd; ++x){
          c.SetCell(y,x,0);
        }
      }
    }/// end loop over size of listxmax
  
    /// Find the weightiest cell in the accumulator.
    int maxCell = 0;
    xMax = 0;
    yMax = 0;
    uint32_t channel = hits[randInd]->Wire()->RawDigit()->Channel();
    wireMax = hits[randInd]->WireID().Wire;

    /// Add the randomly selected point to the accumulator
    maxCell = c.AddPointReturnMax(wireMax, (int)(hits[randInd]->PeakTime()), &yMax, &xMax, fMinHits);
    nAccum++; 

    //mf::LogVerbatim("HoughBaseAlg") << "cout: " << count << " maxCell: " << maxCell << std::endl;
    //mf::LogVerbatim("HoughBaseAlg") << "xMax: " << xMax << " yMax: " << yMax << std::endl;

    /// The threshold calculation, see http://www.via.cornell.edu/ece547/projects/Hough/webpage/statistics.html
    /// accDx is the number of rho bins,m_rowLength
    //TF1 *threshGaus = new TF1("threshGaus","(1/([0]*std::sqrt(2*TMath::Pi())))*exp(-0.5*pow(((x-[1])/[0]),2))");
    //double sigma = std::sqrt(((double)nAccum/(double)accDx)*(1-1/(double)accDx));
    //double mean = (double)nAccum/(double)accDx;
    //threshGaus->SetParameter(0,sigma);
    //threshGaus->SetParameter(1,mean);
    //mf::LogVerbatim("HoughBaseAlg") << "threshGaus mean: " << mean << " sigma: " << sigma << " accDx: " << accDx << std::endl;
    //mf::LogVerbatim("HoughBaseAlg") << "nAccum: " << nAccum << std::endl;
    //mf::LogVerbatim("HoughBaseAlg") << "threshGaus integral range: " << mean-2*sigma << " to " << maxCell << std::endl;
    //mf::LogVerbatim("HoughBaseAlg") << "threshGaus integral: " << threshGaus->Integral(mean-2*sigma,maxCell) << std::endl;
    //mf::LogVerbatim("HoughBaseAlg") << "threshGaus integral: " << threshGaus->Integral(0,maxCell) << std::endl;


    /// The threshold calculation using a Poisson distribution instead
    //double poisProbSum = 0;
    //for(int j = 0; j <= maxCell; j++){
    //double poisProb = TMath::Poisson(j,mean);
    //poisProbSum+=poisProb;
    //mf::LogVerbatim("HoughBaseAlg") << "Poisson: " << poisProb << std::endl;
    //}
    //mf::LogVerbatim("HoughBaseAlg") << "Poisson prob sum: " << poisProbSum << std::endl;
    //mf::LogVerbatim("HoughBaseAlg") << "Probability it is higher: " << 1-poisProbSum << std::endl;

    // Continue if the probability of finding a point, (1-poisProbSum) is the probability of finding a 
    // value of maxCell higher than what it currently is
    //if( (1-poisProbSum) > 1e-13)
    //continue;


    // The threshold calculation using a Poisson distribution instead
    //double binomProbSum = 0;
    //for(int j = 0; j <= maxCell; j++){
    //double binomProb = TMath::BinomialI(1/(double)accDx,nAccum,j);
    //binomProbSum+=binomProb;
    //mf::LogVerbatim("HoughBaseAlg") << "BinomialI: " << binomProb << std::endl;
    //}
    //mf::LogVerbatim("HoughBaseAlg") << "BinomialI prob sum: " << binomProbSum << std::endl;
    //mf::LogVerbatim("HoughBaseAlg") << "Probability it is higher: " << 1-binomProbSum << std::endl;

    // Continue if the probability of finding a point, (1-poisProbSum) is the probability of finding a 
    // value of maxCell higher than what it currently is
    //if( (1-binomProbSum) > 1e-13)
    //continue;





    /// Continue if the biggest maximum for the randomly selected point is smaller than fMinHits
    if (maxCell < fMinHits) 
      continue;


    /// Find the center of mass of the 3x3 cell system (with maxCell at the center). 
    denom = centerofmassx = centerofmassy = 0;
  
    if(xMax > 0 && yMax > 0 && xMax + 1 < accDx && yMax + 1 < accDy){  
      for(int i = -1; i < 2; ++i){
        for(int j = -1; j < 2; ++j){
          denom += c.GetCell(yMax+i,xMax+j);
          centerofmassx += j*c.GetCell(yMax+i,xMax+j);
          centerofmassy += i*c.GetCell(yMax+i,xMax+j);
        }
      }
      centerofmassx /= denom;
      centerofmassy /= denom;      
    }
    else  centerofmassx = centerofmassy = 0;
  
    ///fill the list of cells that have already been found
    listofxmax.push_back(xMax);
    listofymax.push_back(yMax);

    // Find the lines equation
    c.GetEquation(yMax+centerofmassy, xMax+centerofmassx, rho, theta);
    slope = -1./tan(theta);    
    intercept = (rho/sin(theta));
    double distance;
    /// \todo: the collection plane's characteristic hit width's are, 
    /// \todo: on average, about 5 time samples wider than the induction plane's. 
    /// \todo: this is hard-coded for now.
    if(sigt == geo::kInduction)
      indcolscaling = 5.;
    else
      indcolscaling = 0.;
    indcolscaling = 0;

    if(!isinf(slope) && !isnan(slope)){
      sequenceHolder.clear();
      unsigned int fMaxWire = 0;
      int iMaxWire = 0;
      unsigned int fMinWire = 99999999;
      int iMinWire = -1;
      float background = 0;
      hitsTemp.clear();
      for(auto hitsItr = hits.cbegin(); hitsItr != hits.cend(); ++hitsItr){
        wire = (*hitsItr)->WireID().Wire;
        if(fpointId_to_clusterId->at(hitsItr - hits.begin()) != clusterId)
          continue;
        channel = (*hitsItr)->Wire()->RawDigit()->Channel();
        distance = (std::abs((*hitsItr)->PeakTime()-slope*(double)((*hitsItr)->WireID().Wire)-intercept)/(std::sqrt(pow(xyScale*slope,2)+1)));
        if(distance < fMaxDistance+(((*hitsItr)->EndTime()-(*hitsItr)->StartTime())/2.)+indcolscaling && skip[hitsItr-hits.begin()]!=1){
          hitsTemp.push_back(hitsItr-hits.begin());
          sequenceHolder.push_back(channel);
        }
      }/// end iterator over hits
   
      if(hitsTemp.size() < 5) continue;
      currentHits.clear();  
      lastHits.clear();
      int j; 
      currentHits.push_back(0);
      for(auto sequenceHolderItr = sequenceHolder.begin(); sequenceHolderItr+1 != sequenceHolder.end(); ++sequenceHolderItr) {
        j = 1;
        while((chanFilt.BadChannel(*sequenceHolderItr+j)) == true) j++;
        if(sequenceHolder[sequenceHolderItr-sequenceHolder.begin()+1]-sequenceHolder[sequenceHolderItr-sequenceHolder.begin()] <= j + fMissedHits) currentHits.push_back(sequenceHolderItr-sequenceHolder.begin()+1);
        else if(currentHits.size() > lastHits.size()) {
          lastHits = currentHits;
          currentHits.clear();
        }
        else currentHits.clear();
      } 
      if(currentHits.size() > lastHits.size()) lastHits = currentHits;
      clusterHits.clear();    

      if(lastHits.size() < 5) continue;

      // Find minimum and maximum wires in the Hough line
      fMaxWire = 0;
      iMaxWire = 0;
      fMinWire = 99999999;
      iMinWire = -1;
      for(std::vector<int>::iterator lastHitsItr = lastHits.begin(); lastHitsItr != lastHits.end(); ++lastHitsItr) {
        if(fpointId_to_clusterId->at(hitsTemp[(*lastHitsItr)]) != clusterId)
          continue;
        wire = hits[hitsTemp[(*lastHitsItr)]]->WireID().Wire;
        if(wire < fMinWire){
          fMinWire = wire;
          iMinWire = hitsTemp[(*lastHitsItr)];
        }
        if(wire > fMaxWire){
          fMaxWire = wire;
          iMaxWire = hitsTemp[(*lastHitsItr)];
        }
      }

      // Loop over hits to sum up background around prospective Hough line
      for(auto hitsItr = hits.cbegin(); hitsItr != hits.cend(); ++hitsItr){
        distance = (std::abs((*hitsItr)->PeakTime()-slope*(double)((*hitsItr)->WireID().Wire)-intercept)/(std::sqrt(pow(xyScale*slope,2)+1)));

        // Sum up background hits, use smart distance
        if(fMaxDistance+(((*hitsItr)->EndTime()-(*hitsItr)->StartTime())/2.)+indcolscaling < distance 
           && distance < 100*(fMaxDistance+(((*hitsItr)->EndTime()-(*hitsItr)->StartTime())/2.)+indcolscaling) 
	   ){
          double peakTimePerpMin=-(1/slope)*(double)(wire)+hits[iMinWire]->PeakTime()+(1/slope)*(fMinWire);
          double peakTimePerpMax=-(1/slope)*(double)(wire)+hits[iMaxWire]->PeakTime()+(1/slope)*(fMaxWire);
          if((-1/slope) > 0 && hits[iMinWire]->PeakTime() < peakTimePerpMin && hits[iMaxWire]->PeakTime() > peakTimePerpMax)
            background++;
          if((-1/slope) < 0 && hits[iMinWire]->PeakTime() > peakTimePerpMin && hits[iMaxWire]->PeakTime() < peakTimePerpMax)
            background++;
        }
      }// end loop over hits

      const double pos[3] = {0., 0.0, 0.};
      double posWorld0[3] = {0.};
      double posWorld1[3] = {0.};
      const geo::WireGeo& wireGeom = geom->Plane(0).Wire(0);
      const geo::WireGeo& wire1Geom = geom->Plane(0).Wire(1);
      wireGeom.LocalToWorld(pos, posWorld0);
      wire1Geom.LocalToWorld(pos, posWorld1);
      double wire_dist = posWorld0[1]- posWorld1[1];
      double tickToDist = larprop->DriftVelocity(larprop->Efield(),larprop->Temperature());
      tickToDist *= 1.e-3 * detprop->SamplingRate(); // 1e-3 is conversion of 1/us to 1/ns
      if(std::abs(slope)>fMaxSlope ){
        continue;
      }

      //std::cout << std::endl;
      //std::cout  << "Found line!" << std::endl
                                       //<< "Slope: " << slope << std::endl
                                       //<< "Number of hits: " << lastHits.size() << std::endl
                                       //<< "Wire: " << fMinWire << " Peak time: " 
                                       //<< hits[iMinWire]->PeakTime() << std::endl
                                       //<< "Wire: " << fMaxWire << " Peak time: " 
                                       //<< hits[iMaxWire]->PeakTime() << std::endl;


      // Add new line to list of lines
      double totalQ = 0;
      std::vector<std::pair<double,double> > pHit;
      std::vector<std::pair<double,double> > pHitChargeSigma;
      nClustersTemp++;
      ///std::cout << "nClusters: " << *nClusters << std::endl;
      for(auto lastHitsItr = lastHits.begin(); lastHitsItr != lastHits.end(); ++lastHitsItr) {
        fpointId_to_clusterId->at(hitsTemp[(*lastHitsItr)]) = nClustersTemp-1;
        clusterHits.push_back(hits[hitsTemp[(*lastHitsItr)]]);
        totalQ += clusterHits.back()->Charge();
        wire = hits[hitsTemp[(*lastHitsItr)]]->WireID().Wire;

        pHit.push_back(std::make_pair((hits[hitsTemp[(*lastHitsItr)]]->Wire()->RawDigit()->Channel())*wire_dist,
              ((hits[hitsTemp[(*lastHitsItr)]]->StartTime()+hits[hitsTemp[(*lastHitsItr)]]->EndTime())/2.)*tickToDist));
        pHitChargeSigma.push_back(std::make_pair(clusterHits.back()->Charge(),
              clusterHits.back()->SigmaCharge()));
        skip[hitsTemp[(*lastHitsItr)]]=1;

        /// Subtract points from the accumulator that have already been used
        if(accumPoints[hitsTemp[(*lastHitsItr)]]) 
          c.SubtractPoint(wire, (int)(hits[hitsTemp[(*lastHitsItr)]]->PeakTime()));
        
        if(wire < fMinWire){
          fMinWire = wire;
          iMinWire = hitsTemp[(*lastHitsItr)];
        }
        if(wire > fMaxWire){
          fMaxWire = wire;
          iMaxWire = hitsTemp[(*lastHitsItr)];
        }
      }
      double pCornerMin[2];
      pCornerMin[0] = (hits[iMinWire]->Wire()->RawDigit()->Channel())*wire_dist;
      pCornerMin[1] = ((hits[iMinWire]->StartTime()+hits[iMinWire]->EndTime())/2.)*tickToDist;
      double pCornerMax[2];
      pCornerMax[0] = (hits[iMaxWire]->Wire()->RawDigit()->Channel())*wire_dist;
      pCornerMax[1] = ((hits[iMaxWire]->StartTime()+hits[iMaxWire]->EndTime())/2.)*tickToDist;

      ///std::cout << std::endl;
      ///std::cout << "pCornerMin[0]: " << pCornerMin[0] << " pCornerMin[1]: " << pCornerMin[1] << std::endl;
      ///std::cout << "pCornerMax[0]: " << pCornerMax[0] << " pCornerMax[1]: " << pCornerMax[1] << std::endl;
      linesFound.push_back(lineSlope(nClustersTemp-1,
            slope,
            intercept,
            pCornerMin[0],
            pCornerMin[1],
            pCornerMax[0],
            pCornerMax[1],
            iMinWire,
            iMaxWire,
            pHit,
            pHitChargeSigma));
      ///std::cout << "isolation: " << background/lastHits.size() << std::endl;
       
    }/// end if !isnan

    if(linesFound.size()>(unsigned int)fMaxLines)
      break;

  }/// end loop over hits



  /// Loop over hits to sum up background around the Hough lines
  for(auto linesFoundItr = linesFound.begin(); linesFoundItr < linesFound.end(); linesFoundItr++){
    float totalBkgDist = 0;
    float totalBkgDistCharge = 0;
    int   totalBkg = 0;
    float background = 0;
    for(auto hitsItr = hits.cbegin(); hitsItr != hits.cend(); ++hitsItr){
      /// Veto the hit if it already belongs to a line
      if(fpointId_to_clusterId->at(hitsItr-hits.cbegin()) != clusterId)
        continue;
      float distance = (TMath::Abs((*hitsItr)->PeakTime()-linesFoundItr->clusterSlope*(double)((*hitsItr)->WireID().Wire)-linesFoundItr->clusterIntercept)/(std::sqrt(pow(xyScale*linesFoundItr->clusterSlope,2)+1)));
      /// Sum up background hits, use smart distance
      double peakTimePerpMin=-(1/linesFoundItr->clusterSlope)*(double)((*hitsItr)->WireID().Wire)+hits[linesFoundItr->iMinWire]->PeakTime()+(1/linesFoundItr->clusterSlope)*(hits[linesFoundItr->iMinWire]->WireID().Wire);
      double peakTimePerpMax=-(1/linesFoundItr->clusterSlope)*(double)((*hitsItr)->WireID().Wire)+hits[linesFoundItr->iMaxWire]->PeakTime()+(1/linesFoundItr->clusterSlope)*(hits[linesFoundItr->iMaxWire]->WireID().Wire);
      ///if(fMaxDistance+(((*hitsItr)->EndTime()-(*hitsItr)->StartTime())/2.)+indcolscaling < distance
         ///&& distance < 10*(fMaxDistance+(((*hitsItr)->EndTime()-(*hitsItr)->StartTime())/2.)+indcolscaling)){
      if(distance < 10*(fMaxDistance+(((*hitsItr)->EndTime()-(*hitsItr)->StartTime())/2.)+indcolscaling)){
          totalBkgDist+=distance;
          totalBkgDistCharge+=distance/(*hitsItr)->Charge();
          totalBkg++;
      }
      if(fMaxDistance+(((*hitsItr)->EndTime()-(*hitsItr)->StartTime())/2.)+indcolscaling < distance 
         && distance < 50*(fMaxDistance+(((*hitsItr)->EndTime()-(*hitsItr)->StartTime())/2.)+indcolscaling) 
         ){
        if((-1/slope) > 0 && (*hitsItr)->PeakTime() < peakTimePerpMin && (*hitsItr)->PeakTime() > peakTimePerpMax){
          background++;
        }
        if((-1/slope) < 0 && (*hitsItr)->PeakTime() > peakTimePerpMin && (*hitsItr)->PeakTime() < peakTimePerpMax){
          background++;
        }
      }
    }/// end loop over hits
    //std::cout << std::endl << "total BkgDist: " << totalBkgDist << std::endl;
    //std::cout << std::endl << "total BkgDistCharge: " << totalBkgDistCharge << std::endl;
    linesFoundItr->showerLikeness = totalBkgDistCharge;
    linesFoundItr->isolation = background/linesFoundItr->pHit.size();
  }/// end loop over lines found

  /// Do a merge based on distances between line segments instead of endpoints
  std::map<int,int> mNLineMerges;
  if(fDoShowerHoughLineMerge)     mergeHoughLinesBySegment(0,&linesFound,xyScale,iMergeShower,&mNLineMerges);
  if(fDoShowerPartHoughLineMerge) mergeHoughLinesBySegment(0,&linesFound,xyScale,iMergeShowerPart,&mNLineMerges);
  if(fDoHoughLineMerge)           mergeHoughLinesBySegment(0,&linesFound,xyScale,iMergeNormal,&mNLineMerges);

  /// Average the slopes and y-intercepts of Hough lines in shower like regions
  /// map (cluster number ( slope, y-int))
  std::map<int,std::pair<double,double> > showerLines; /// map (cluster number ( slope, y-int))
  std::map<int,std::pair<int,int> > showerMinMaxWires; /// map (cluster number ( iMinWire, iMaxWire))
  std::map<int,int> linesCount; /// map (cluster number, number of lines)
  for(auto linesFoundItr = linesFound.begin(); linesFoundItr != linesFound.end(); linesFoundItr++){
    if(linesFoundItr->showerLikeness<10)
      continue;
    int lineSize = linesFoundItr->pHit.size();
    double lineSlope = linesFoundItr->clusterSlope;
    double lineIntercept = linesFoundItr->clusterIntercept;
    if(!linesCount.count(linesFoundItr->clusterNumber)){
      showerLines.insert(std::make_pair(linesFoundItr->clusterNumber,std::make_pair(lineSlope*lineSize,lineIntercept*lineSize)));
      showerMinMaxWires.insert(std::make_pair(linesFoundItr->clusterNumber,std::make_pair(linesFoundItr->iMinWire,linesFoundItr->iMaxWire)));
      linesCount.insert(std::make_pair(linesFoundItr->clusterNumber,lineSize));
      continue;
    }
    if(linesFoundItr->iMinWire < showerMinMaxWires.at(linesFoundItr->clusterNumber).first)
      showerMinMaxWires.at(linesFoundItr->clusterNumber).first=linesFoundItr->iMinWire;
    if(linesFoundItr->iMaxWire > showerMinMaxWires.at(linesFoundItr->clusterNumber).second)
      showerMinMaxWires.at(linesFoundItr->clusterNumber).second=linesFoundItr->iMaxWire;
    showerLines.at(linesFoundItr->clusterNumber).first+=lineSlope*lineSize;
    showerLines.at(linesFoundItr->clusterNumber).second+=lineIntercept*lineSize;
    linesCount.at(linesFoundItr->clusterNumber)+=lineSize;
  }
  /// Now actually do the averaging and find line direction
  std::map<int,std::pair<double,double> > showerLinesAverage; // map (cluster number ( slope, y-int))
  for(auto linesCountItr = linesCount.cbegin(); linesCountItr != linesCount.cend(); linesCountItr++){
    showerLines.at(linesCountItr->first).first/=(linesCountItr->second);
    showerLines.at(linesCountItr->first).second/=(linesCountItr->second);
    double averageSlope = showerLines.at(linesCountItr->first).first;
    double averageInt = showerLines.at(linesCountItr->first).second;
    int midWire = hits[(showerMinMaxWires.at(linesCountItr->first).first)]->WireID().Wire/2 + 
      hits[(showerMinMaxWires.at(linesCountItr->first).second)]->WireID().Wire/2 + 
      (hits[(showerMinMaxWires.at(linesCountItr->first).first)]->WireID().Wire & hits[showerMinMaxWires.at(linesCountItr->first).second]->WireID().Wire & 1);
    double midPeakTime = averageSlope*midWire+averageInt;
    showerLinesAverage.insert(std::make_pair(linesCountItr->first,std::make_pair(showerLines.at(linesCountItr->first).first,showerLines.at(linesCountItr->first).second)));
    //std::cout << "slope: " << showerLines.at(linesCountItr->first).first <<
      //" intercept: " << showerLines.at(linesCountItr->first).second 
      //<< " min wire: " << hits[showerMinMaxWires.at(linesCountItr->first).first]->WireID().Wire
      //<< " max wire: " << hits[showerMinMaxWires.at(linesCountItr->first).second]->WireID().Wire
      //<< " mid wire: " << midWire
      //<< " midPeakTime: " << midPeakTime
      //<< std::endl;
    double directionBin1=0;
    double directionBin2=0;
    for(auto hitsItr = hits.cbegin(); hitsItr != hits.cend(); ++hitsItr){
      // Veto the hit if it already belongs to a line
      //if(fpointId_to_clusterId->at(hitsItr-hits.cbegin()) != clusterId)
        //continue;
      float distance = (TMath::Abs((*hitsItr)->PeakTime()-averageSlope*(double)((*hitsItr)->WireID().Wire)-averageInt)/(std::sqrt(pow(xyScale*averageSlope,2)+1)));
      //if(fMaxDistance+(((*hitsItr)->EndTime()-(*hitsItr)->StartTime())/2.)+indcolscaling < distance 
         //&& distance < 10*(fMaxDistance+(((*hitsItr)->EndTime()-(*hitsItr)->StartTime())/2.)+indcolscaling) 
         //){
      if(distance < 10000*(fMaxDistance+(((*hitsItr)->EndTime()-(*hitsItr)->StartTime())/2.)+indcolscaling)){
        double peakTimePerpMid=-(1/averageSlope)*(double)((*hitsItr)->WireID().Wire)+midPeakTime+(1/averageSlope)*(midWire);
        if((-1/slope) < 0 && (*hitsItr)->PeakTime() < peakTimePerpMid)
          directionBin1+=distance/(*hitsItr)->Charge();
        if((-1/slope) > 0 && (*hitsItr)->PeakTime() > peakTimePerpMid)
          directionBin1+=distance/(*hitsItr)->Charge();
        if((-1/slope) > 0 && (*hitsItr)->PeakTime() < peakTimePerpMid)
          directionBin2+=distance/(*hitsItr)->Charge();
        if((-1/slope) < 0 && (*hitsItr)->PeakTime() > peakTimePerpMid)
          directionBin2+=distance/(*hitsItr)->Charge();
      }
    }
    //std::cout << "directionBin1: " << directionBin1 << " directionBin2: " << directionBin2 << std::endl;
    //int showerDirection = 0;
    //if(directionBin1<directionBin2){
      //showerDirection = 1;
      //std::cout << "I'm probably moving right" << std::endl;
    //}
    //else{
      //showerDirection = -1;
      //std::cout << "I'm probably moving left" << std::endl;
    //}
    int smallUnmergedLines = 0;
    int toTheLeft = 0;
    int toTheRight = 0;
    //int showerDirection = 0; // +1 is to the right, -1 is to the left
    for(auto linesFoundItr = linesFound.begin(); linesFoundItr < linesFound.end(); linesFoundItr++){
      if(linesFoundItr->merged)
        continue;
      double segmentDistance = std::min(std::abs(hits[linesFoundItr->iMaxWire]->PeakTime()-(showerLines.at(linesCountItr->first).first*hits[linesFoundItr->iMaxWire]->WireID().Wire+showerLines.at(linesCountItr->first).second)),std::abs(hits[linesFoundItr->iMaxWire]->PeakTime()-(showerLines.at(linesCountItr->first).first*hits[linesFoundItr->iMaxWire]->WireID().Wire+showerLines.at(linesCountItr->first).second)));
      //std::cout << "segmentDistance: " << segmentDistance << std::endl;
      //if(segmentDistance < 100)
        //linesFoundItr->clusterNumber = linesCountItr->first;
      if(segmentDistance < 100)
        smallUnmergedLines++;
      if(hits[linesFoundItr->iMinWire]->WireID().Wire > 
          (hits[showerMinMaxWires.at(linesCountItr->first).first]->WireID().Wire + 
           hits[showerMinMaxWires.at(linesCountItr->first).second]->WireID().Wire))
        toTheRight++;
      if(hits[linesFoundItr->iMaxWire]->WireID().Wire <= 
          (hits[showerMinMaxWires.at(linesCountItr->first).first]->WireID().Wire + 
           hits[showerMinMaxWires.at(linesCountItr->first).second]->WireID().Wire))
        toTheLeft++;
      //std::cout << "hits[linesFoundItr->iMinWire]->WireID().Wire " << hits[linesFoundItr->iMinWire]->WireID().Wire << std::endl;
      //std::cout << "hits[linesFoundItr->iMaxWire]->WireID().Wire " << hits[linesFoundItr->iMaxWire]->WireID().Wire << std::endl;
      //std::cout << "hits[showerMinMaxWires.at(linesCountItr->first).first]->WireID().Wire " << hits[showerMinMaxWires.at(linesCountItr->first).first]->WireID().Wire << std::endl;
    }
    //std::cout << "smallUnmergedLines: " << smallUnmergedLines << std::endl;
    //std::cout << "toTheRight: " << toTheRight << std::endl;
    //std::cout << "toTheLeft: " << toTheLeft << std::endl;
    //if(toTheRight > toTheLeft){
      ////showerDirection = 1;
      //std::cout << "I'm probably moving right" << std::endl;
    //}
    //else if(toTheRight < toTheLeft){
      ////showerDirection = -1;
      //std::cout << "I'm probably moving left" << std::endl;
    //}
    //for(auto linesFoundItr = linesFound.begin(); linesFoundItr < linesFound.end(); linesFoundItr++){
      //if(linesFoundItr->merged)
        //continue;
      //double segmentDistance = std::min(std::abs(hits[linesFoundItr->iMaxWire]->PeakTime()-(showerLines.at(linesCountItr->first).first*hits[linesFoundItr->iMaxWire]->WireID().Wire+showerLines.at(linesCountItr->first).second)),std::abs(hits[linesFoundItr->iMaxWire]->PeakTime()-(showerLines.at(linesCountItr->first).first*hits[linesFoundItr->iMaxWire]->WireID().Wire+showerLines.at(linesCountItr->first).second)));
      //if(segmentDistance >= 100)
        //continue;
      //if(showerDirection > 0)
        //if(hits[linesFoundItr->iMinWire]->WireID().Wire > 
            //(hits[showerMinMaxWires.at(linesCountItr->first).first]->WireID().Wire + 
             //hits[showerMinMaxWires.at(linesCountItr->first).second]->WireID().Wire)/2)
          //linesFoundItr->clusterNumber = linesCountItr->first;
      //if(showerDirection < 0)
        //if(hits[linesFoundItr->iMaxWire]->WireID().Wire <=
            //(hits[showerMinMaxWires.at(linesCountItr->first).first]->WireID().Wire + 
             //hits[showerMinMaxWires.at(linesCountItr->first).second]->WireID().Wire)/2)
          //linesFoundItr->clusterNumber = linesCountItr->first;
    //}
  }







  for( auto itr=mNLineMerges.begin(); itr!=mNLineMerges.end(); itr++){
    mf::LogVerbatim("HoughBaseAlg") << "cluster num: " 
                                    << (*itr).first 
                                    << " number of Hough line merges: " 
                                    << (*itr).second << std::endl;
  }

  // Reassign the merged lines
  for(auto fpointId_to_clusterIdItr = fpointId_to_clusterId->begin(); fpointId_to_clusterIdItr != fpointId_to_clusterId->end(); ++fpointId_to_clusterIdItr){
    if(*fpointId_to_clusterIdItr == clusterId)
      continue;
    for(auto linesFoundItr = linesFound.begin(); linesFoundItr < linesFound.end(); linesFoundItr++){
      if(*fpointId_to_clusterIdItr == linesFoundItr->oldClusterNumber)
        *fpointId_to_clusterIdItr = linesFoundItr->clusterNumber;
    }
  }

  // Find sizes of all merged lines combined
  // For linesFoundSizes, key is cluster number and size is the mapped value
  //std::map<int,float> linesFoundSizes;
  //for(unsigned int i = 0; i < linesFound.size(); i++){
    ////if(linesFoundSizes.find(linesFound[i].clusterNumber) == linesFoundSizes.end())
    //if(!linesFoundSizes.count(linesFound[i].clusterNumber))
      //linesFoundSizes[linesFound[i].clusterNumber] = std::sqrt( pow(linesFound[i].pMin0-linesFound[i].pMax0,2)+pow(linesFound[i].pMin1-linesFound[i].pMax1,2));
    //else 
      //linesFoundSizes[linesFound[i].clusterNumber]+= std::sqrt( pow(linesFound[i].pMin0-linesFound[i].pMax0,2)+pow(linesFound[i].pMin1-linesFound[i].pMax1,2));
    //mf::LogInfo("HoughBaseAlg") << "i: " << i 
				//<< " linesFound[i].clusterNumber: " 
				//<<  linesFound[i].clusterNumber 
				//<< " linesFound[i] size: " 
				//<< std::sqrt( pow(linesFound[i].pMin0-linesFound[i].pMax0,2)
					 //+pow(linesFound[i].pMin1-linesFound[i].pMax1,2)) 
				//<< " size: " 
				//<< linesFoundSizes[linesFound[i].clusterNumber];
  //}

  
  // Merge remains of original fuzzy cluster into nearest Hough line, 
  // assumes the Hough line is a segment
  for(auto hitsItr = hits.cbegin(); hitsItr != hits.cend(); ++hitsItr){
    if(fpointId_to_clusterId->at(hitsItr-hits.begin()) != clusterId)
      continue;
      double p0 = ((*hitsItr)->Wire()->RawDigit()->Channel())*wire_dist;
      double p1 = (((*hitsItr)->StartTime()+(*hitsItr)->EndTime())/2.)*tickToDist;
      double minDistance = 10000;
    for(std::vector<lineSlope>::iterator linesFoundItr = linesFound.begin(); linesFoundItr < linesFound.end(); linesFoundItr++){
      double distance = PointSegmentDistance( p0, p1, linesFoundItr->pMin0, linesFoundItr->pMin1, linesFoundItr->pMax0, linesFoundItr->pMax1);
      //distance/=std::sqrt( pow(linesFound[k].pMin0-linesFound[k].pMax0,2)+pow(linesFound[k].pMin1-linesFound[k].pMax1,2));
      //distance/=linesFoundSizes[linesFound[k].clusterNumber];
      if(distance < minDistance){
        fpointId_to_clusterId->at(hitsItr-hits.begin()) = linesFoundItr->clusterNumber;
        minDistance = distance;
      }
    }
  } 

  hit.clear();
  lastHits.clear();

  listofxmax.clear();
  listofymax.clear();

  // saves a bitmap image of the accumulator (useful for debugging), 
  // with scaling based on the maximum cell value
  if(fSaveAccumulator){   
    unsigned char *outPix = new unsigned char [accDx*accDy];
    //finds the maximum cell in the accumulator for image scaling
    int cell, pix = 0, maxCell = 0;
    for (int y = 0; y < accDy; ++y){ 
      for (int x = 0; x < accDx; ++x){
        cell = c.GetCell(y,x);
        if (cell > maxCell) maxCell = cell;
      }
    }
    for (int y = 0; y < accDy; ++y){
      for (int x = 0; x < accDx; ++x){ 
        //scales the pixel weights based on the maximum cell value     
        if(maxCell > 0)
          pix = (int)((1500*c.GetCell(y,x))/maxCell);
        outPix[y*accDx + x] = pix;
      }
    }
                
    HLSSaveBMPFile("houghaccum.bmp", outPix, accDx, accDy);
    delete [] outPix;
  }// end if saving accumulator

  *nClusters=nClustersTemp;

  return 1; 
}


// Merges based on the distance between line segments
void cluster::HoughBaseAlg::mergeHoughLinesBySegment(unsigned int clusIndexStart,
						     std::vector<lineSlope> *linesFound,
						     double xyScale,
                                                     int mergeStyle,
                                                     std::map<int,int> *mNLineMerges)
{

  // If we have zero or one Hough lines, move on 
  if(linesFound->size() == 0 || linesFound->size() == 1)
    return;

  // If a merge happened, trigger this function again to look at this cluster again!
  bool lineMerged = false;

  // If we reach the last Hough line, move on 
  if(linesFound->size() == clusIndexStart+1)
    return;

  // Min to merge
  std::vector<unsigned int> toMerge; 
  std::vector<double> mergeSlope;
  std::vector<double> mergeTheta;

  // Check if segments are close enough
  for(auto linesFoundClusIndStItr = linesFound->begin(); linesFoundClusIndStItr != linesFound->end(); linesFoundClusIndStItr++){
    if(linesFound->at(clusIndexStart).clusterNumber != linesFoundClusIndStItr->clusterNumber)
      continue;
    for(auto linesFoundToMergeItr = linesFound->begin(); linesFoundToMergeItr != linesFound->end(); linesFoundToMergeItr++){
      if(linesFound->at(clusIndexStart).clusterNumber == linesFoundToMergeItr->clusterNumber)
        continue;
      double segmentDistance = HoughLineDistance(linesFoundClusIndStItr->pMin0,linesFoundClusIndStItr->pMin1,
                                                 linesFoundClusIndStItr->pMax0,linesFoundClusIndStItr->pMax1, 
						 linesFoundToMergeItr->pMin0,linesFoundToMergeItr->pMin1,
                                                 linesFoundToMergeItr->pMax0,linesFoundToMergeItr->pMax1);
      if( (segmentDistance<fHoughLineMergeCutoff && mergeStyle == iMergeNormal) 
          || (segmentDistance<fShowerHoughLineMergeCutoff && mergeStyle == iMergeShower)
          || (segmentDistance<fShowerPartHoughLineMergeCutoff && mergeStyle == iMergeShowerPart))
	{
	  toMerge.push_back(linesFoundToMergeItr-linesFound->begin());
	  mergeSlope.push_back(linesFoundClusIndStItr->clusterSlope*xyScale);
	}
    }
  }

  mergeTheta.resize(toMerge.size());

  // Find the angle between the slopes
  for(auto mergeThetaItr = mergeTheta.begin(); mergeThetaItr != mergeTheta.end(); mergeThetaItr++){
    double toMergeSlope =  linesFound->at(toMerge[mergeThetaItr-mergeTheta.begin()]).clusterSlope*xyScale;
    mergeTheta[mergeThetaItr-mergeTheta.begin()] = atan(std::abs(( toMergeSlope - mergeSlope[mergeThetaItr-mergeTheta.begin()])/(1 + toMergeSlope*mergeSlope[mergeThetaItr-mergeTheta.begin()] )))*(180/TMath::Pi());
  }

  // Perform the merge
  for(auto toMergeItr = toMerge.begin(); toMergeItr != toMerge.end(); toMergeItr++){
    if( (mergeTheta[toMergeItr-toMerge.begin()] < fHoughLineMergeAngle && mergeStyle == iMergeNormal) 
        || (mergeTheta[toMergeItr-toMerge.begin()] < fShowerHoughLineMergeAngle && mergeStyle == iMergeShower)
        || (mergeTheta[toMergeItr-toMerge.begin()] < fShowerPartHoughLineMergeAngle && mergeStyle == iMergeShowerPart)){
      
      // First check averages of charge and sigma charge for hits in lines closest to each other
      int closestToMerge=-1;
      int closestClusIndexStart=-1;
      double closestDistance=999999;
      for (auto toMergePHitItr = linesFound->at(*toMergeItr).pHit.begin(); toMergePHitItr != linesFound->at(*toMergeItr).pHit.end(); toMergePHitItr++) {
        for (auto clusIndStPHitItr = linesFound->at(clusIndexStart).pHit.begin(); clusIndStPHitItr != linesFound->at(clusIndexStart).pHit.end(); clusIndStPHitItr++) {
          double distance = std::sqrt(pow(clusIndStPHitItr->first-(*toMergePHitItr).first,2)+
                    pow(clusIndStPHitItr->second-toMergePHitItr->second,2));
          if(distance < closestDistance){
            closestDistance = distance;
            closestToMerge=toMergePHitItr-linesFound->at(*toMergeItr).pHit.begin();
            closestClusIndexStart=clusIndStPHitItr-linesFound->at(clusIndexStart).pHit.begin();
          }
        }
      }

      // Find four points closest to closestToMerge on the toMerge[i] line
      int closestToMerge1=-1;
      int closestToMerge2=-1;
      int closestToMerge3=-1;
      int closestToMerge4=-1;
      double closestToMergeDist1=999999;
      double closestToMergeDist2=999999;
      double closestToMergeDist3=999999;
      double closestToMergeDist4=999999;
      for (auto toMergePHitItr = linesFound->at(*toMergeItr).pHit.begin(); toMergePHitItr != linesFound->at(*toMergeItr).pHit.end(); toMergePHitItr++) {
        if(closestToMerge==toMergePHitItr-linesFound->at(*toMergeItr).pHit.begin())
          continue;
        double distance = std::sqrt(pow(toMergePHitItr->first-linesFound->at(*toMergeItr).pHit[closestToMerge].first,2)+
                  pow(toMergePHitItr->second-linesFound->at(*toMergeItr).pHit[closestToMerge].second,2));
        if(closestToMergeDist1 > distance){
          closestToMergeDist4 = closestToMergeDist3;
          closestToMerge4 = closestToMerge3;
          closestToMergeDist3 = closestToMergeDist2;
          closestToMerge3 = closestToMerge2;
          closestToMergeDist2 = closestToMergeDist1;
          closestToMerge2 = closestToMerge1;
          closestToMergeDist1 = distance;
          closestToMerge1 = toMergePHitItr-linesFound->at(*toMergeItr).pHit.begin();
        }
        else if(closestToMergeDist2 > distance){
          closestToMergeDist4 = closestToMergeDist3;
          closestToMerge4 = closestToMerge3;
          closestToMergeDist3 = closestToMergeDist2;
          closestToMerge3 = closestToMerge2;
          closestToMergeDist2 = distance;
          closestToMerge2 = toMergePHitItr-linesFound->at(*toMergeItr).pHit.begin();
        }
        else if(closestToMergeDist3 > distance){
          closestToMergeDist4 = closestToMergeDist3;
          closestToMerge4 = closestToMerge3;
          closestToMergeDist3 = distance;
          closestToMerge3 = toMergePHitItr-linesFound->at(*toMergeItr).pHit.begin();
        }
        else if(closestToMergeDist4 > distance){
          closestToMergeDist4 = distance;
          closestToMerge4 = toMergePHitItr-linesFound->at(*toMergeItr).pHit.begin();
        }
      }
      //std::cout << std::endl;
      //std::cout << "first toMerge: " << linesFound->at(*toMergeItr).pHit[closestToMerge].first  << " " << linesFound->at(*toMergeItr).pHit[closestToMerge].second << std::endl;
      //std::cout << "second toMerge: "<< linesFound->at(*toMergeItr).pHit[closestToMerge1].first << " " << linesFound->at(*toMergeItr).pHit[closestToMerge1].second << std::endl;
      //std::cout << "third toMerge: " << linesFound->at(*toMergeItr).pHit[closestToMerge2].first << " " << linesFound->at(*toMergeItr).pHit[closestToMerge2].second << std::endl;
      //std::cout << "fourth toMerge: "<< linesFound->at(*toMergeItr).pHit[closestToMerge3].first << " " << linesFound->at(*toMergeItr).pHit[closestToMerge3].second << std::endl;

      //std::cout << "first toMerge charge: "  << linesFound->at(*toMergeItr).pHitChargeSigma[closestToMerge].first << std::endl;
      //std::cout << "second toMerge charge: " << linesFound->at(*toMergeItr).pHitChargeSigma[closestToMerge1].first << std::endl;
      //std::cout << "third toMerge charge: "  << linesFound->at(*toMergeItr).pHitChargeSigma[closestToMerge2].first << std::endl;
      //std::cout << "fourth toMerge charge: " << linesFound->at(*toMergeItr).pHitChargeSigma[closestToMerge3].first << std::endl;
      //std::cout << "fifth toMerge charge: "  << linesFound->at(*toMergeItr).pHitChargeSigma[closestToMerge4].first << std::endl;


      // Find three points closest to closestToMerge on the clusIndexStart line
      int closestClusIndexStart1=-1;
      int closestClusIndexStart2=-1;
      int closestClusIndexStart3=-1;
      int closestClusIndexStart4=-1;
      double closestClusIndexStartDist1=999999;
      double closestClusIndexStartDist2=999999;
      double closestClusIndexStartDist3=999999;
      double closestClusIndexStartDist4=999999;
      for (auto clusIndStPHitItr = linesFound->at(clusIndexStart).pHit.begin(); clusIndStPHitItr != linesFound->at(clusIndexStart).pHit.end(); clusIndStPHitItr++) {
        if(closestClusIndexStart==clusIndStPHitItr-linesFound->at(clusIndexStart).pHit.begin())
          continue;
        double distance = std::sqrt(pow(clusIndStPHitItr->first-linesFound->at(clusIndexStart).pHit[closestClusIndexStart].first,2)+
                  pow(clusIndStPHitItr->second-linesFound->at(clusIndexStart).pHit[closestClusIndexStart].second,2));
        if(closestClusIndexStartDist1 > distance){
          closestClusIndexStartDist4 = closestClusIndexStartDist3;
          closestClusIndexStart4 = closestClusIndexStart3;
          closestClusIndexStartDist3 = closestClusIndexStartDist2;
          closestClusIndexStart3 = closestClusIndexStart2;
          closestClusIndexStartDist2 = closestClusIndexStartDist1;
          closestClusIndexStart2 = closestClusIndexStart1;
          closestClusIndexStartDist1 = distance;
          closestClusIndexStart1 = clusIndStPHitItr-linesFound->at(clusIndexStart).pHit.begin();
        }
        else if(closestClusIndexStartDist2 > distance){
          closestClusIndexStartDist4 = closestClusIndexStartDist3;
          closestClusIndexStart4 = closestClusIndexStart3;
          closestClusIndexStartDist3 = closestClusIndexStartDist2;
          closestClusIndexStart3 = closestClusIndexStart2;
          closestClusIndexStartDist2 = distance;
          closestClusIndexStart2 = clusIndStPHitItr-linesFound->at(clusIndexStart).pHit.begin();
        }
        else if(closestClusIndexStartDist3 > distance){
          closestClusIndexStartDist4 = closestClusIndexStartDist3;
          closestClusIndexStart4 = closestClusIndexStart3;
          closestClusIndexStartDist3 = distance;
          closestClusIndexStart3 = clusIndStPHitItr-linesFound->at(clusIndexStart).pHit.begin();
        }
        else if(closestClusIndexStartDist4 > distance){
          closestClusIndexStartDist4 = distance;
          closestClusIndexStart4 = clusIndStPHitItr-linesFound->at(clusIndexStart).pHit.begin();
        }
      }
      //std::cout << "first clusIndex: " << linesFound->at(clusIndexStart).pHit[closestClusIndexStart].first << " " << linesFound->at(clusIndexStart).pHit[closestClusIndexStart].second << std::endl;
      //std::cout << "second clusIndex: " << linesFound->at(clusIndexStart).pHit[closestClusIndexStart1].first << " " << linesFound->at(clusIndexStart).pHit[closestClusIndexStart1].second << std::endl;
      //std::cout << "third clusIndex: " << linesFound->at(clusIndexStart).pHit[closestClusIndexStart2].first << " " << linesFound->at(clusIndexStart).pHit[closestClusIndexStart2].second << std::endl;
      //std::cout << "fourth clusIndex: " << linesFound->at(clusIndexStart).pHit[closestClusIndexStart3].first << " " << linesFound->at(clusIndexStart).pHit[closestClusIndexStart3].second << std::endl;

      //std::cout << "first clusIndex charge: "  << linesFound->at(clusIndexStart).pHitChargeSigma[closestClusIndexStart].first << std::endl;
      //std::cout << "second clusIndex charge: " << linesFound->at(clusIndexStart).pHitChargeSigma[closestClusIndexStart1].first << std::endl;
      //std::cout << "third clusIndex charge: "  << linesFound->at(clusIndexStart).pHitChargeSigma[closestClusIndexStart2].first << std::endl;
      //std::cout << "fourth clusIndex charge: " << linesFound->at(clusIndexStart).pHitChargeSigma[closestClusIndexStart3].first << std::endl;
      //std::cout << "fifth clusIndex charge: "  << linesFound->at(clusIndexStart).pHitChargeSigma[closestClusIndexStart4].first << std::endl;

      //Determine average charge and average sigma charge for each line around closest hits
      double toMergeAveCharge = linesFound->at(*toMergeItr).pHitChargeSigma[closestToMerge].first+
        linesFound->at(*toMergeItr).pHitChargeSigma[closestToMerge1].first+
        linesFound->at(*toMergeItr).pHitChargeSigma[closestToMerge2].first+
        linesFound->at(*toMergeItr).pHitChargeSigma[closestToMerge3].first+
        linesFound->at(*toMergeItr).pHitChargeSigma[closestToMerge4].first;
      toMergeAveCharge/=5;
      double clusIndexStartAveCharge = linesFound->at(clusIndexStart).pHitChargeSigma[closestClusIndexStart].first+
        linesFound->at(clusIndexStart).pHitChargeSigma[closestClusIndexStart1].first+
        linesFound->at(clusIndexStart).pHitChargeSigma[closestClusIndexStart2].first+
        linesFound->at(clusIndexStart).pHitChargeSigma[closestClusIndexStart3].first+
        linesFound->at(clusIndexStart).pHitChargeSigma[closestClusIndexStart4].first;
      clusIndexStartAveCharge/=5;
      double toMergeAveSigmaCharge = linesFound->at(*toMergeItr).pHitChargeSigma[closestToMerge].second+
        linesFound->at(*toMergeItr).pHitChargeSigma[closestToMerge1].second+
        linesFound->at(*toMergeItr).pHitChargeSigma[closestToMerge2].second+
        linesFound->at(*toMergeItr).pHitChargeSigma[closestToMerge3].second+
        linesFound->at(*toMergeItr).pHitChargeSigma[closestToMerge4].second;
      toMergeAveSigmaCharge/=5;
      double clusIndexStartAveSigmaCharge = linesFound->at(clusIndexStart).pHitChargeSigma[closestClusIndexStart].second+
        linesFound->at(clusIndexStart).pHitChargeSigma[closestClusIndexStart1].second+
        linesFound->at(clusIndexStart).pHitChargeSigma[closestClusIndexStart2].second+
        linesFound->at(clusIndexStart).pHitChargeSigma[closestClusIndexStart3].second+
        linesFound->at(clusIndexStart).pHitChargeSigma[closestClusIndexStart4].second;
      clusIndexStartAveSigmaCharge/=5;

      double chargeAsymmetry = std::abs(toMergeAveCharge-clusIndexStartAveCharge)/(toMergeAveCharge+clusIndexStartAveCharge);
      double sigmaChargeAsymmetry = std::abs(toMergeAveSigmaCharge-clusIndexStartAveSigmaCharge)/(toMergeAveSigmaCharge+clusIndexStartAveSigmaCharge);

      //std::cout << "toMergeAveCharge: " << toMergeAveCharge << std::endl;
      //std::cout << "clusIndexStartAveCharge: " << clusIndexStartAveCharge << std::endl;
      //std::cout << "toMergeAveSigmaCharge: " << toMergeAveSigmaCharge << std::endl;
      //std::cout << "clusIndexStartAveSigmaCharge: " << clusIndexStartAveSigmaCharge << std::endl;
      //std::cout << "charge asymmetry: " << chargeAsymmetry << std::endl;
      //std::cout << "sigma charge asymmetry: " << sigmaChargeAsymmetry << std::endl;



      if(chargeAsymmetry > fChargeAsymmetryCut &&
          mergeStyle == iMergeNormal)
        continue;

      if(sigmaChargeAsymmetry > fSigmaChargeAsymmetryCut &&
          mergeStyle == iMergeNormal)
        continue;
     
      // Check if both lines is in region that looks showerlike
      if(mergeStyle == iMergeShower){
        if(!(linesFound->at(*toMergeItr).showerLikeness>fMergeShowerLikenessCut) || !(linesFound->at(clusIndexStart).showerLikeness>fMergeShowerLikenessCut))
          continue;
      }

      // Check if one (and exactly one) of the lines is in region that looks showerlike
      //if(mergeStyle == iMergeShowerPart){
        //if((linesFound->at(*toMergeItr).showerLike && !linesFound->at(clusIndexStart).showerLike)
            //|| (!linesFound->at(*toMergeItr).showerLike && linesFound->at(clusIndexStart).showerLike))
          //continue;
      //}

      lineMerged = true;
      linesFound->at(clusIndexStart).merged = true;
      linesFound->at(*toMergeItr).merged = true;

      //for(unsigned int j = 0; j < linesFound->size(); j++){
      for(auto linesFoundItr = linesFound->begin(); linesFoundItr != linesFound->end(); linesFoundItr++){
        if((unsigned int)(*toMergeItr) == linesFoundItr-linesFound->begin())
          continue;

        if(linesFoundItr->clusterNumber == linesFound->at(*toMergeItr).clusterNumber)
          linesFoundItr->clusterNumber = linesFound->at(clusIndexStart).clusterNumber;
      }
      linesFound->at(*toMergeItr).clusterNumber = linesFound->at(clusIndexStart).clusterNumber;
     
      // Count up how many merges we've had
      if(mNLineMerges->count(linesFound->at(clusIndexStart).clusterNumber)==0){
        mNLineMerges->insert(std::pair<int,int>(linesFound->at(clusIndexStart).clusterNumber,1));
      } else
      mNLineMerges->at(linesFound->at(clusIndexStart).clusterNumber)++;

    }  
  }

  if(lineMerged)
    mergeHoughLinesBySegment(clusIndexStart,linesFound,xyScale,mergeStyle,mNLineMerges);
  else
    mergeHoughLinesBySegment(clusIndexStart+1,linesFound,xyScale,mergeStyle,mNLineMerges);
  
  return;

}




//------------------------------------------------------------------------------
double cluster::HoughBaseAlg::HoughLineDistance(double p0MinLine1, 
						double p1MinLine1, 
						double p0MaxLine1, 
						double p1MaxLine1, 
						double p0MinLine2, 
						double p1MinLine2, 
						double p0MaxLine2, 
						double p1MaxLine2)
{
  //distance between two segments in the plane:
  //  one segment is (x11, y11) to (x12, y12) or (p0MinLine1, p1MinLine1) to (p0MaxLine1, p1MaxLine1)
  //  the other is   (x21, y21) to (x22, y22) or (p0MinLine2, p1MinLine2) to (p0MaxLine2, p1MaxLine2)
  double x11 = p0MinLine1; 
  double y11 = p1MinLine1; 
  double x12 = p0MaxLine1; 
  double y12 = p1MaxLine1; 
  double x21 = p0MinLine2; 
  double y21 = p1MinLine2; 
  double x22 = p0MaxLine2; 
  double y22 = p1MaxLine2; 

  if(HoughLineIntersect(x11, y11, x12, y12, x21, y21, x22, y22)) return 0;
  // try each of the 4 vertices w/the other segment
  std::vector<double> distances;
  distances.push_back(PointSegmentDistance(x11, y11, x21, y21, x22, y22));
  distances.push_back(PointSegmentDistance(x12, y12, x21, y21, x22, y22));
  distances.push_back(PointSegmentDistance(x21, y21, x11, y11, x12, y12));
  distances.push_back(PointSegmentDistance(x22, y22, x11, y11, x12, y12));

  double minDistance = 999999; 
  for(unsigned int j = 0; j < distances.size(); j++){
    if (distances[j] < minDistance)
      minDistance = distances[j];
  }
  
  return minDistance;

}




//------------------------------------------------------------------------------
bool cluster::HoughBaseAlg::HoughLineIntersect(double x11,
					       double  y11,
					       double  x12,
					       double  y12,
					       double  x21,
					       double  y21,
					       double  x22,
					       double  y22)
{
  //whether two segments in the plane intersect:
  //one segment is (x11, y11) to (x12, y12)
  //the other is   (x21, y21) to (x22, y22)
  
  double dx1 = x12 - x11; // x2-x1
  double dy1 = y12 - y11; // y2-y1
  double dx2 = x22 - x21; // x4-x3
  double dy2 = y22 - y21; // y4-y3
  //double delta = dx2*dy1 - dy2*dx1; // (x4-x3)(y2-y1) - (y4-y3)(x2-x1)
  double delta = dy2*dx1 - dx2*dy1; // (y4-y3)(x2-x1) - (x4-x3)(y2-y1) 
  if (delta == 0) return false;  // parallel segments

  double t = (dx2*(y11 - y21) + dy2*(x21 - x11)) / delta; // ua
  double s = (dx1*(y11 - y21) + dy1*(x21 - x11)) / delta; // ub
  
  return (0 <= s && s <= 1 && 0 <= t && t <= 1);

}



//------------------------------------------------------------------------------
double cluster::HoughBaseAlg::PointSegmentDistance(double px,
						   double  py,
						   double  x1,
						   double  y1,
						   double  x2,
						   double  y2)
{
  double dx = x2 - x1;
  double dy = y2 - y1;
  if ( dx == 0 && dy == 0 )  // the segment's just a point
    return std::sqrt( pow(px - x1,2) + pow(py - y1,2));

  // Calculate the t that minimizes the distance.
  double t = ((px - x1)*dx + (py - y1)*dy) / (dx*dx + dy*dy);

  // See if this represents one of the segment's
  // end points or a point in the middle.
  if(t < 0){
    dx = px - x1;
    dy = py - y1;
  }
  else if(t > 1) {
    dx = px - x2;
    dy = py - y2;
  }
  else if(0 <= t && t <= 1) {
    double near_x = x1 + t * dx;
    double near_y = y1 + t * dy;
    dx = px - near_x;
    dy = py - near_y;
  }

  return std::sqrt(dx*dx + dy*dy);

}

//------------------------------------------------------------------------------
cluster::HoughTransform::~HoughTransform()
{  
}

//------------------------------------------------------------------------------
bool cluster::HoughTransform::AddPoint(int x, int y)
{
  if (x>m_dx || y>m_dy || x<0.0 || y<0.0)
    return false;
  return DoAddPoint(x, y);
}


//------------------------------------------------------------------------------
int cluster::HoughTransform::AddPointReturnMax(int x, 
					       int y, 
					       int *yMax, 
					       int *xMax, 
					       int minHits)
{
  if (x>m_dx || y>m_dy || x<0.0 || y<0.0)
    return 0;
  return DoAddPointReturnMax(x, y, yMax, xMax, minHits);
}



//------------------------------------------------------------------------------
bool cluster::HoughTransform::SubtractPoint(int x, int y)
{
  if (x>m_dx || y>m_dy || x<0.0 || y<0.0)
    return false;
  return DoSubtractPoint(x, y);
}


//------------------------------------------------------------------------------
void cluster::HoughTransform::Init(int dx, 
				   int dy, 
				   int rhores,
				   int numACells)
{
  m_numAngleCells=numACells;
  m_rhoResolutionFactor = rhores;
  m_accum.clear();
  m_accum.resize(m_numAngleCells);
  m_numAccumulated = 0;   
  //   m_cosTable.clear();
  //   m_sinTable.clear();
  m_cosTable.resize(m_numAngleCells);
  m_sinTable.resize(m_numAngleCells);
  //if (dx == m_dx && dy == m_dy)
  //return;
  m_dx = dx;
  m_dy = dy;
  m_rowLength = (int)(m_rhoResolutionFactor*2. * std::sqrt(dx*dx + dy*dy));
  
  int angleIndex;
  double a, angleStep = TMath::Pi()/m_numAngleCells;
  for (a=0.0, angleIndex = 0; angleIndex < m_numAngleCells; ++angleIndex){
    m_cosTable[angleIndex] = cos(a);
    m_sinTable[angleIndex] = sin(a);
    a += angleStep;
  }
}


//------------------------------------------------------------------------------
int cluster::HoughTransform::GetMax(int &xmax, 
				    int &ymax)
{
  int maxVal = -1;
  for(unsigned int i = 0; i < m_accum.size(); i++){
    for(auto rhoIter=m_accum[i].begin(); rhoIter!=m_accum[i].end(); ++rhoIter){
      if((*rhoIter).second > maxVal) {
	maxVal = (*rhoIter).second;
	xmax = i;
	ymax = (*rhoIter).first;
      }
    }
  }

  return maxVal;
}


//------------------------------------------------------------------------------
int cluster::HoughTransform::DoAddPointReturnMax(int x, 
						 int y, 
						 int *ymax, 
						 int *xmax, 
						 int minHits)
{
  int distCenter = (int)(m_rowLength/2.);
 
  // prime the lastDist variable so our linear fill works below
  int lastDist = (int)(distCenter+(m_rhoResolutionFactor*(m_cosTable[0]*x + m_sinTable[0]*y)));

  //int max_val = minHits-1;
  int max_val = 0;

  // loop through all angles a from 0 to 180 degrees
  for (unsigned int a=1; a<m_cosTable.size(); ++a){
    // Calculate the basic line equation dist = cos(a)*x + sin(a)*y.
    // Shift to center of row to cover negative values
    int dist = (int)(distCenter+(m_rhoResolutionFactor*(m_cosTable[a]*x + m_sinTable[a]*y)));
    // sanity check to make sure we stay within our row
    if (dist >= 0 && dist<m_rowLength){
      if(lastDist==dist){
	m_accum[a][lastDist]++;
        if( max_val < m_accum[a][lastDist]){
          max_val = m_accum[a][lastDist];
          *xmax = lastDist;
          *ymax = a;
        }
        //mf::LogVerbatim("HoughBaseAlg") << "First, a: " << a << " lastDist: " << lastDist << std::endl;
      }
      else{
	// fill in all values in row a, not just a single cell
	int stepDir = dist>lastDist ? 1 : -1;
	for (int cell=lastDist; cell!=dist; cell+=stepDir){   
	  m_accum[a][cell]++;//maybe add weight of hit here?
          // Note, m_accum is a vector of associative containers, "a" calls the vector element, "cell" is the container key, and the ++ iterates the value correspoding to the key
          if(max_val < m_accum[a][cell]){
            max_val = m_accum[a][cell];
            *xmax = cell;
            *ymax = a;
          }
	}      
      }
    }      
    lastDist = dist;
  }
  m_numAccumulated++;
  
  //mf::LogVerbatim("HoughBaseAlg") << "Add point says xmax: " << *xmax << " ymax: " << *ymax << std::endl;

  return max_val;
}


//------------------------------------------------------------------------------
bool cluster::HoughTransform::DoAddPoint(int x, int y)
{
  int distCenter = (int)(m_rowLength/2.);
 
  // prime the lastDist variable so our linear fill works below
  int lastDist = (int)(distCenter+(m_rhoResolutionFactor*(m_cosTable[0]*x + m_sinTable[0]*y)));

  // loop through all angles a from 0 to 180 degrees
  for (unsigned int a=1; a<m_cosTable.size(); ++a){
    // Calculate the basic line equation dist = cos(a)*x + sin(a)*y.
    // Shift to center of row to cover negative values
    int dist = (int)(distCenter+(m_rhoResolutionFactor*(m_cosTable[a]*x + m_sinTable[a]*y)));
    // sanity check to make sure we stay within our row
    if (dist >= 0 && dist<m_rowLength){
      if(lastDist==dist){
	m_accum[a][lastDist]++;
        //mf::LogVerbatim("HoughBaseAlg") << "First, a: " << a << " lastDist: " << lastDist << std::endl;
      }
      else{
	// fill in all values in row a, not just a single cell
	int stepDir = dist>lastDist ? 1 : -1;
	for (int cell=lastDist; cell!=dist; cell+=stepDir){   
	  m_accum[a][cell]++;//maybe add weight of hit here?
          // Note, m_accum is a vector of associative containers, "a" calls the vector element, "cell" is the container key, and the ++ iterates the value correspoding to the key
	}      
      }
    }      
    lastDist = dist;
  }
  m_numAccumulated++;

  return true;
}

//------------------------------------------------------------------------------
bool cluster::HoughTransform::DoSubtractPoint(int x, int y)
{
  int distCenter = (int)(m_rowLength/2.);
 
  // prime the lastDist variable so our linear fill works below
  int lastDist = (int)(distCenter+(m_rhoResolutionFactor*(m_cosTable[0]*x + m_sinTable[0]*y)));
  // loop through all angles a from 0 to 180 degrees
  for (unsigned int a=1; a<m_cosTable.size(); ++a){
    // Calculate the basic line equation dist = cos(a)*x + sin(a)*y.
    // Shift to center of row to cover negative values
    int dist = (int)(distCenter+(m_rhoResolutionFactor*(m_cosTable[a]*x + m_sinTable[a]*y)));
    // sanity check to make sure we stay within our row
    if (dist >= 0 && dist<m_rowLength){
      if(lastDist==dist)
	m_accum[a][lastDist]--;
      else{
	// fill in all values in row a, not just a single cell
	int stepDir = dist>lastDist ? 1 : -1;
	for (int cell=lastDist; cell!=dist; cell+=stepDir){   
	  m_accum[a][cell]--;//maybe add weight of hit here?
          // Note, m_accum is a vector of associative containers, "a" calls the vector element, "cell" is the container key, and the -- iterates the value correspoding to the key
	}      
      }
    }      
    lastDist = dist;
  }
  m_numAccumulated--;

  return true;
}

//------------------------------------------------------------------------------
//this method saves a BMP image of the Hough Accumulator, which can be viewed with gimp
void cluster::HoughBaseAlg::HLSSaveBMPFile(const char *fileName, unsigned char *pix, int dx, int dy)
{
  ofstream bmpFile(fileName, std::ios::binary);
  bmpFile.write("B", 1);
  bmpFile.write("M", 1);
  int bitsOffset = 54 +256*4; 
  int size = bitsOffset + dx*dy; //header plus 256 entry LUT plus pixels
  bmpFile.write((const char *)&size, 4);
  int reserved = 0;
  bmpFile.write((const char *)&reserved, 4);
  bmpFile.write((const char *)&bitsOffset, 4);
  int bmiSize = 40;
  bmpFile.write((const char *)&bmiSize, 4);
  bmpFile.write((const char *)&dx, 4);
  bmpFile.write((const char *)&dy, 4);
  short planes = 1;
  bmpFile.write((const char *)&planes, 2);
  short bitCount = 8;
  bmpFile.write((const char *)&bitCount, 2);
  int i, temp = 0;
  for (i=0; i<6; i++)
    bmpFile.write((const char *)&temp, 4);  // zero out optional color info
  // write a linear LUT
  char lutEntry[4]; // blue,green,red
  lutEntry[3] = 0;  // reserved part
  for (i=0; i<256; i++)
    {
      lutEntry[0] = lutEntry[1] = lutEntry[2] = i;
      bmpFile.write(lutEntry, sizeof lutEntry);
    }
  // write the actual pixels
  bmpFile.write((const char *)pix, dx*dy);
}
 

//------------------------------------------------------------------------------
size_t cluster::HoughBaseAlg::FastTransform(std::vector<art::Ptr<recob::Cluster> >         & clusIn,
					    std::vector<recob::Cluster>                    & ccol,  
					    std::vector< art::PtrVector<recob::Hit> >      & clusHitsOut,
					    art::Event                                const& evt,
					    std::string                               const& label)
{
  std::vector<int> skip;  

  art::FindManyP<recob::Hit> fmh(clusIn, evt, label);

  art::ServiceHandle<geo::Geometry> geom;
  art::ServiceHandle<util::LArProperties> larprop;
  art::ServiceHandle<util::DetectorProperties> detprop;
  filter::ChannelFilter chanFilt;
  HoughTransform c;

  // Get the random number generator
  art::ServiceHandle<art::RandomNumberGenerator> rng;
  CLHEP::HepRandomEngine & engine = rng -> getEngine();
  CLHEP::RandFlat flat(engine);

  extern void SaveBMPFile(const char *f, unsigned char *pix, int dxx, int dyy);
  std::vector< art::Ptr<recob::Hit> > hit;

  for(size_t cs = 0; cs < geom->Ncryostats(); ++cs){
    for(size_t t = 0; t < geom->Cryostat(cs).NTPC(); ++t){
      for(unsigned int p = 0; p < geom->Cryostat(cs).TPC(t).Nplanes(); ++p) {
	art::PtrVector<recob::Cluster>::const_iterator clusterIter = clusIn.begin();
	int clusterID = 0;//the unique ID of the cluster


	geo::View_t    view = geom->Cryostat(cs).TPC(t).Plane(p).View();
	geo::SigType_t sigt = geom->Cryostat(cs).TPC(t).Plane(p).SignalType();

	size_t cinctr = 0;
	while(clusterIter != clusIn.end()) {
	  hit.clear();
	  if(fPerCluster){
	    if((*clusterIter)->View() == view) hit = fmh.at(cinctr);
	  }
	  else{   
	    while(clusterIter != clusIn.end()){
	      if( (*clusterIter)->View() == view ){

		hit = fmh.at(cinctr);
	      }// end if cluster is in correct view
	      clusterIter++;
	      ++cinctr;
	    }//end loop over clusters
	  }//end if not fPerCluster
	  if(hit.size() == 0){ 
	    if(fPerCluster){
	      clusterIter++;
	      ++cinctr;
	    }
	    continue;
	  }
	  //factor to make x and y scale the same units
	  double xyScale  = .001*larprop->DriftVelocity(larprop->Efield(),larprop->Temperature());
	  xyScale        *= detprop->SamplingRate()/geom->WirePitch(0,1,p,t,cs);
	  
	  int x, y;
	  int dx = geom->Cryostat(cs).TPC(t).Plane(p).Nwires();//number of wires 
	  int dy = hit[0]->Wire()->NSignal();//number of time samples. 
	  skip.clear();
	  skip.resize(hit.size());
	  std::vector<int> listofxmax;
	  std::vector<int> listofymax;  
	  std::vector<int> hitTemp;        //indecies ofcandidate hits
	  std::vector<int> sequenceHolder; //channels of hits in list
	  std::vector<int> currentHits;    //working vector of hits 
	  std::vector<int> lastHits;       //best list of hits
	  art::PtrVector<recob::Hit> clusterHits;
	  double indcolscaling = 0.;       //a parameter to account for the different 
	                                   //characteristic hit width of induction and collection plane
	  double centerofmassx = 0;
	  double centerofmassy = 0;
	  double denom = 0; 
	  double intercept=0.;
	  double slope = 0.;
	  //this array keeps track of the hits that have already been associated with a line. 
	  int xMax = 0;
	  int yMax = 0;
	  double rho;
	  double theta; 
	  int accDx(0), accDy(0);


	  //Init specifies the size of the two-dimensional accumulator 
	  //(based on the arguments, number of wires and number of time samples). 
	  //adds all of the hits (that have not yet been associated with a line) to the accumulator
	  c.Init(dx,dy,fRhoResolutionFactor,fNumAngleCells);


          // count is how many points are left to randomly insert
          unsigned int count = hit.size();
          std::vector<unsigned int> accumPoints;
          accumPoints.resize(hit.size());
          int nAccum = 0;
          unsigned int nLinesFound = 0;

          for( ; count > 0; count--){
	   

            // The random hit we are examining
            unsigned int randInd = (unsigned int)(flat.fire()*hit.size());
        
            // Skip if it's already in a line
            if(skip[randInd]==1)
              continue;
        
            // If we have already accumulated the point, skip it
            if(accumPoints[randInd])
              continue;
            accumPoints[randInd]=1;
    
            // zeroes out the neighborhood of all previous lines  
            for(unsigned int i = 0; i < listofxmax.size(); ++i){
              int yClearStart = listofymax[i] - fRhoZeroOutRange;
              if (yClearStart < 0) yClearStart = 0;
              
              int yClearEnd = listofymax[i] + fRhoZeroOutRange;
              if (yClearEnd >= accDy) yClearEnd = accDy - 1;
              
              int xClearStart = listofxmax[i] - fThetaZeroOutRange;
              if (xClearStart < 0) xClearStart = 0;
              
              int xClearEnd = listofxmax[i] + fThetaZeroOutRange;
              if (xClearEnd >= accDx) xClearEnd = accDx - 1;
              
              for (y = yClearStart; y <= yClearEnd; ++y){
                for (x = xClearStart; x <= xClearEnd; ++x){
                  c.SetCell(y,x,0);
                }
              }
            }// end loop over size of listxmax
  
            //find the weightiest cell in the accumulator.
            int maxCell = 0;
            unsigned int wireMax = hit[randInd]->WireID().Wire;
            xMax = 0;
            yMax = 0;

            // Add the randomly selected point to the accumulator
            maxCell = c.AddPointReturnMax(wireMax, (int)(hit[randInd]->PeakTime()), &yMax, &xMax, fMinHits);
            nAccum++; 

            // Continue if the biggest maximum for the randomly selected point is smaller than fMinHits
            if (maxCell < fMinHits) 
              continue;

            // Find the center of mass of the 3x3 cell system (with maxCell at the center). 
            denom = centerofmassx = centerofmassy = 0;
          
            if(xMax > 0 && yMax > 0 && xMax + 1 < accDx && yMax + 1 < accDy){  
              for(int i = -1; i < 2; ++i){
                for(int j = -1; j < 2; ++j){
                  denom += c.GetCell(yMax+i,xMax+j);
                  centerofmassx += j*c.GetCell(yMax+i,xMax+j);
                  centerofmassy += i*c.GetCell(yMax+i,xMax+j);
                }
              }
              centerofmassx /= denom;
              centerofmassy /= denom;      
            }
            else  centerofmassx = centerofmassy = 0;
          
            //fill the list of cells that have already been found
            listofxmax.push_back(xMax);
            listofymax.push_back(yMax);
        
            // Find the lines equation
            c.GetEquation(yMax+centerofmassy, xMax+centerofmassx, rho, theta);
            //c.GetEquation(yMax, xMax, rho, theta);
            slope = -1./tan(theta);    
            intercept = (rho/sin(theta));
            //mf::LogVerbatim("HoughBaseAlg") << std::endl;
            //mf::LogVerbatim("HoughBaseAlg") << "slope: " << slope << " intercept: " << intercept << std::endl; 
            //mf::LogInfo("HoughBaseAlg") << "slope: " << slope << " intercept: " << intercept;
            double distance;
            /// \todo: the collection plane's characteristic hit width's are, 
            /// \todo: on average, about 5 time samples wider than the induction plane's. 
            /// \todo: this is hard-coded for now.
            if(sigt == geo::kInduction)
              indcolscaling = 5.;
            else
              indcolscaling = 0.;
            // What is this?   
            indcolscaling = 0;






	    if(!isinf(slope) && !isnan(slope)){
	      sequenceHolder.clear();
	      hitTemp.clear();
	      for(size_t i = 0; i < hit.size(); ++i){
		distance = (TMath::Abs(hit[i]->PeakTime()-slope*(double)(hit[i]->WireID().Wire)-intercept)/(std::sqrt(pow(xyScale*slope,2)+1)));
		
		if(distance < fMaxDistance+((hit[i]->EndTime()-hit[i]->StartTime())/2.)+indcolscaling  && skip[i]!=1){
		  hitTemp.push_back(i);
		  sequenceHolder.push_back(hit[i]->Channel());
		}
		
	      }// end loop over hits
	      
	      if(hitTemp.size() < 2) continue;
	      currentHits.clear();  
	      lastHits.clear();
	      int j; 
	      currentHits.push_back(0);
	      for(size_t i = 0; i + 1 < sequenceHolder.size(); ++i){  
		j = 1;
		while((chanFilt.BadChannel(sequenceHolder[i]+j)) == true) j++;
		if(sequenceHolder[i+1]-sequenceHolder[i] <= j + fMissedHits) currentHits.push_back(i+1);
		else if(currentHits.size() > lastHits.size()) {
		  lastHits = currentHits;
		  currentHits.clear();
		}
		else currentHits.clear();
	      } 
	      
	      if(currentHits.size() > lastHits.size()) lastHits = currentHits;
	      clusterHits.clear();    
	      double totalQ = 0.;
              const double pos[3] = {0., 0.0, 0.};
              double posWorld0[3] = {0.};
              double posWorld1[3] = {0.};
              const geo::WireGeo& wireGeom = geom->Plane(0).Wire(0);
              const geo::WireGeo& wire1Geom = geom->Plane(0).Wire(1);
              wireGeom.LocalToWorld(pos, posWorld0);
              wire1Geom.LocalToWorld(pos, posWorld1);
              double tickToDist = larprop->DriftVelocity(larprop->Efield(),larprop->Temperature());
              tickToDist *= 1.e-3 * detprop->SamplingRate(); // 1e-3 is conversion of 1/us to 1/ns

              if(lastHits.size() < 5) continue;


	      for(size_t i = 0; i < lastHits.size(); ++i) {
		clusterHits.push_back(hit[hitTemp[lastHits[i]]]);
		totalQ += clusterHits.back()->Charge();
		skip[hitTemp[lastHits[i]]]=1;
	      } 
	      //protection against very steep uncorrelated hits
	      if(std::abs(slope)>fMaxSlope 
		 && std::abs((*clusterHits.begin())->Wire()->RawDigit()->Channel()-
			       clusterHits[clusterHits.size()-1]->Wire()->RawDigit()->Channel())>=0
		 )
		continue;
	     


	      unsigned int sw = (*clusterHits.begin())->WireID().Wire;
	      unsigned int ew = (*(clusterHits.end()-1))->WireID().Wire;
	      
	      recob::Cluster cluster(sw, 0.,
				     (*clusterHits.begin())->PeakTime(), 0.,
				     ew, 0., 
				     (clusterHits[clusterHits.size()-1])->PeakTime(), 0.,
				     slope, 0., 
				     -999., 0., 
				     totalQ,
				     geom->View((*clusterHits.begin())->Channel()),
				     clusterID);	      
	      
	      ++clusterID;
	      ccol.push_back(cluster);
	      clusHitsOut.push_back(clusterHits);

	      //Turn off hit sharing. T. Yang 9/14/12
	      //	      //allow double assignment of first and last hits
	      //	      for(size_t i = 0; i < lastHits.size(); ++i){ 
	      //		if(skip[hitTemp[lastHits[i]]] ==1){
	      //		  channel = hit[hitTemp[lastHits[i]]]->Wire()->RawDigit()->Channel();  
	      //		  if( channel == sc || channel == ec) skip[i] = 0;
	      //		}
	      //	      }
              
	    }// end if !isnan
   
            nLinesFound++;

            if(nLinesFound>(unsigned int)fMaxLines)
              break;

	    
	  }// end loop over hits
	  
	  // saves a bitmap image of the accumulator (useful for debugging), 
	  // with scaling based on the maximum cell value
	  if(fSaveAccumulator){   
	    unsigned char *outPix = new unsigned char [accDx*accDy];
	    //finds the maximum cell in the accumulator for image scaling
	    int cell, pix = 0, maxCell = 0;
	    for (y = 0; y < accDy; ++y){ 
	      for (x = 0; x < accDx; ++x){
		cell = c.GetCell(y,x);
		if (cell > maxCell) maxCell = cell;
	      }
	    }
	    for (y = 0; y < accDy; ++y){
	      for (x = 0; x < accDx; ++x){ 
		//scales the pixel weights based on the maximum cell value     
		if(maxCell > 0)
		  pix = (int)((1500*c.GetCell(y,x))/maxCell);
		outPix[y*accDx + x] = pix;
	      }
	    }
	    	    
	    SaveBMPFile("houghaccum.bmp", outPix, accDx, accDy);
	    delete [] outPix;
	  }// end if saving accumulator
	  
	  hit.clear();
	  lastHits.clear();
	  if(clusterIter != clusIn.end()){
	    clusterIter++;
	    ++cinctr;
	  }
	  listofxmax.clear();
	  listofymax.clear();
	}//end loop over clusters

      }//end loop over planes
    }// end loop over tpcs
  }// end loop over cryostats

  return ccol.size(); 






}



//------------------------------------------------------------------------------
size_t cluster::HoughBaseAlg::Transform(std::vector< art::Ptr<recob::Hit> > const& hits,
					double                                   & slope,
					double                                   & intercept)
{
  HoughTransform c;

  art::ServiceHandle<geo::Geometry> geom;
  int dx = geom->Nwires(0);               //number of wires 
  int dy = hits[0]->Wire()->NSignal();//number of time samples. 

  c.Init(dx,dy,fRhoResolutionFactor,fNumAngleCells);

  for(unsigned int i=0;i < hits.size(); ++i){
    c.AddPoint(hits[i]->WireID().Wire, (int)(hits[i]->PeakTime()));
  }// end loop over hits

  //gets the actual two-dimensional size of the accumulator
  int accDx = 0;
  int accDy = 0;
  c.GetAccumSize(accDy, accDx);

  //find the weightiest cell in the accumulator.
  int xMax = 0;
  int yMax = 0;
  c.GetMax(yMax,xMax);

  //find the center of mass of the 3x3 cell system (with maxCell at the center). 
  double centerofmassx = 0.;
  double centerofmassy = 0.;
  double denom         = 0.;
    
  if(xMax > 0 && yMax > 0 && xMax+1 < accDx && yMax+1 < accDy){  
    for(int i = -1; i < 2; ++i){
      for(int j = -1; j < 2; ++j){
	denom         += c.GetCell(yMax+i,xMax+j);
	centerofmassx += j*c.GetCell(yMax+i,xMax+j);
	centerofmassy += i*c.GetCell(yMax+i,xMax+j);
      }
    }
    centerofmassx /= denom;
    centerofmassy /= denom;      
  }
  else  centerofmassx = centerofmassy = 0;

  double rho   = 0.;
  double theta = 0.;
  c.GetEquation(yMax+centerofmassy, xMax+centerofmassx, rho, theta);
  slope     = -1./tan(theta);    
  intercept = rho/sin(theta);
  
  ///\todo could eventually refine this method to throw out hits that are 
  ///\todo far from the hough line and refine the fit

  return hits.size();
}







