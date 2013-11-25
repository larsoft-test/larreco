////////////////////////////////////////////////////////////////////////
//
// fuzzyClusterAlg.cxx
//
// Ben Carls, bcarls@fnal.gov
//
// This code looks for clusters using a Gustafson-Kessel variant on fuzzy c-means algorithm. The
// clusters are then examined by the HoughBaseAlg to identify Hough lines
// which can then be split off into their own clusters. See the webpage below
// for more information on the fuzzy clustering algorithm.
//
//http://homes.di.unimi.it/~valenti/SlideCorsi/Bioinformatica05/Fuzzy-Clustering-lecture-Babuska.pdf
//http://biosoft.kaist.ac.kr/BISL_homepage/publication/p20050006.pdf
////////////////////////////////////////////////////////////////////////


#include <boost/bind.hpp>

//Framework includes:
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


#include "fhiclcpp/ParameterSet.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "CLHEP/Random/RandFlat.h"
#include "Filters/ChannelFilter.h"

#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "RecoAlg/fuzzyClusterAlg.h"
#include "RecoBase/Hit.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"
#include "Utilities/AssociationUtil.h"


#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>

// Define parameters that will tell us if we are doing a normal Hough line merge
// or a shower Hough line merge
static const int iMergeShower          = 0;
static const int iMergeNormal          = 1;
static const int iMergeShowerIntercept = 2;
static const int iMergeChargeAsymAngle = 3;


namespace cluster{
  const unsigned int kNO_CLUSTER    = UINT_MAX;
  const unsigned int kNOISE_CLUSTER = UINT_MAX-1;
}

//----------------------------------------------------------
// fuzzyClusterAlg stuff
//----------------------------------------------------------
cluster::fuzzyClusterAlg::fuzzyClusterAlg(fhicl::ParameterSet const& pset) 
   : fHBAlg(pset.get< fhicl::ParameterSet >("HoughBaseAlg"))
   , fDBScan(pset.get< fhicl::ParameterSet >("DBScanAlg"))
{
 this->reconfigure(pset); 
}

//----------------------------------------------------------
cluster::fuzzyClusterAlg::~fuzzyClusterAlg()
{
}

//----------------------------------------------------------
void cluster::fuzzyClusterAlg::reconfigure(fhicl::ParameterSet const& p)
{
  fFuzzifier                      = p.get< float  >("Fuzzifier");
  fMaxNumClusters                 = p.get< int    >("MaxNumClusters");
  nIterations                     = p.get< int    >("Iterations");
  fDoFuzzyRemnantMerge            = p.get< int    >("DoFuzzyRemnantMerge"            );
  fFuzzyRemnantMergeCutoff        = p.get< double >("FuzzyRemnantMergeCutoff"        );
  fDoHoughLineMerge               = p.get< int    >("DoHoughLineMerge"               );
  fHoughLineMergeAngle            = p.get< double >("HoughLineMergeAngle"            );
  fHoughLineMergeCutoff           = p.get< double >("HoughLineMergeCutoff"           );
  fDoChargeAsymAngleMerge         = p.get< int    >("DoChargeAsymAngleMerge"         );
  fChargeAsymAngleCut             = p.get< double >("ChargeAsymAngleCut"             );
  fSigmaChargeAsymAngleCut        = p.get< double >("SigmaChargeAsymAngleCut"        );
  fChargeAsymAngleCutoff          = p.get< double >("ChargeAsymAngleCutoff"          );
  fDoShowerHoughLineMerge         = p.get< int    >("DoShowerHoughLineMerge"         );
  fShowerHoughLineMergeAngle      = p.get< double >("ShowerHoughLineMergeAngle"      );
  fShowerHoughLineMergeCutoff     = p.get< double >("ShowerHoughLineMergeCutoff"     );
  fDoShowerHoughLineInterceptMerge= p.get< int    >("DoShowerHoughLineInterceptMerge");
  fShowerLikenessCut              = p.get< double >("ShowerLikenessCut"              );
  fShowerWidthAngle               = p.get< double >("ShowerWidthAngle"               );
  fHBAlg.reconfigure(p.get< fhicl::ParameterSet >("HoughBaseAlg"));
  fDBScan.reconfigure(p.get< fhicl::ParameterSet >("DBScanAlg"));
}

//----------------------------------------------------------
void cluster::fuzzyClusterAlg::InitFuzzy(std::vector<art::Ptr<recob::Hit> >& allhits, 
					 std::set<uint32_t>                  badChannels)
{
  // clear all the data member vectors for the new set of hits
  fps.clear();
  fpointId_to_clusterId.clear();
  fnoise.clear();
  fvisited.clear();
  fsim.clear();
  fsim2.clear();
  fsim3.clear();
  fclusters.clear();
  fWirePitch.clear();

  fBadChannels = badChannels;
  fBadWireSum.clear();

  // Clear the matrix that stores the points needed for fuzzy clustering
  fpsMat.Clear();
  fpsMembership.Clear();
  fpsNewMembership.Clear();
  fpsCentroids.Clear(); 

  //------------------------------------------------------------------
  // Determine spacing between wires (different for each detector)
  ///get 2 first wires and find their spacing (wire_dist)

  art::ServiceHandle<util::LArProperties> larp;
  art::ServiceHandle<util::DetectorProperties> detp;
  
  // Collect the hits in a useful form,
  // and take note of the maximum time width
  fMaxWidth=0.0;
  //fpsMat = TMatrixT<double>(allhits.size(),2);
  //fpsMembership = TMatrixT<double>(iNumClusters, allhits.size());
  fpsMat.ResizeTo(allhits.size(),2);
  for (unsigned int j = 0; j < allhits.size(); ++j){
    int dims = 3;//our point is defined by 3 elements:wire#,center of the hit, and the hit width
    std::vector<float> p(dims);
        
    float tickToDist = larp->DriftVelocity(larp->Efield(),larp->Temperature());
    tickToDist *= 1.e-3 * detp->SamplingRate(); // 1e-3 is conversion of 1/us to 1/ns
    p[0] = (allhits[j]->Channel())*fGeom->WirePitch(fGeom->View(allhits[j]->Channel()));
    p[1] = ((allhits[j]->StartTime()+allhits[j]->EndTime()  )/2.)*tickToDist;
    p[2] = (allhits[j]->EndTime()  -allhits[j]->StartTime())*tickToDist;   //width of a hit in cm

    // check on the maximum width condition
    if ( p[2] > fMaxWidth ) fMaxWidth = p[2];
    
    fps.push_back(p);

    // Store hits in the matrix needed for fuzzy clustering
    fpsMat(j,0) = p[0];
    fpsMat(j,1) = p[1];

  }

  mf::LogInfo("fuzzyCluster") << "InitFuzzy: hits vector size is " << fps.size();

  return;
}


//----------------------------------------------------------
void cluster::fuzzyClusterAlg::computeCentroids(int k)
{
  // fpsCentroids are the weighted centers of the clusters.
  // We multiply fpsMembership by fpsMat to find fpsCentroids, then
  // divide by the normalization (sum of the weights).

  int iNumClusters = k;

  //fpsMembership.Print();
  //fpsMat.Print();
  float normalizationFactor;


  // Determine the elements of u^m_ij
  TMatrixT<float> Uji_m(iNumClusters, fpsMat.GetNrows());
  // For each cluster
  for ( int j = 0; j < iNumClusters; j++)
    // For each hit
    for ( int i = 0; i < fpsMat.GetNrows(); i++)
      // Determine Uji_m
      Uji_m(j,i) = pow(fpsMembership(j,i), fFuzzifier); 


  // Now find sum^N_i=1 u^m_ij*x_i
  // For each cluster
  for ( int j = 0; j < iNumClusters; j++)
    // For each hit
    for ( int i = 0; i < fpsMat.GetNrows(); i++)
      // For each dimension
      for ( int f = 0; f < 2; f++)
        fpsCentroids(j,f) += Uji_m(j,i)*fpsMat(i,f);


  // Divide centroids by the normalization (sum of weights)
  // For each cluster
  for ( int j = 0; j < iNumClusters; j++){
    normalizationFactor = 0;
    // For each hit
    for ( int i = 0; i < fpsMat.GetNrows(); i++)
      normalizationFactor += Uji_m(j,i);
    // For each dimension
    for ( int f = 0; f < 2; f++)
      fpsCentroids(j,f) /= normalizationFactor;
  }


}
//----------------------------------------------------------
void cluster::fuzzyClusterAlg::computeCentroids2(int k)
{
  // Centroids are defined by c_j = (sum^N_i=1 u^m_ij*x_i)/(sum^N_i=1 u^m_ij)

  int iNumClusters = k;
  TMatrixT<float> Uji_m(iNumClusters, fpsMat.GetNrows());
  float normalizationFactor;

  // Zero the centroid matrix
  for ( int j = 0; j < fpsCentroids.GetNrows(); j++)
    for ( int i = 0; i < fpsCentroids.GetNcols(); i++)
      fpsCentroids(j,i) = 0;

  // Determine the elements of u^m_ij
  // For each cluster
  for ( int j = 0; j < iNumClusters; j++)
    // For each hit
    for ( int i = 0; i < fpsMat.GetNrows(); i++)
      // Determine Uji_m
      Uji_m(j,i) = pow ( fpsMembership(j,i), fFuzzifier); 

  // Now find sum^N_i=1 u^m_ij*x_i
  // For each cluster
  for ( int j = 0; j < iNumClusters; j++)
    // For each hit
    for ( int i = 0; i < fpsMat.GetNrows(); i++)
      // For each dimension
      for ( int f = 0; f < 2; f++)
        fpsCentroids(j,f) += Uji_m(j,i)*fpsMat(i,f);

  // Divide centroids by the normalization (sum of weights)
  // For each cluster
  for ( int j = 0; j < iNumClusters; j++){
    normalizationFactor = 0;
    // For each hit
    for ( int i = 0; i < fpsMat.GetNrows(); i++)
      normalizationFactor += Uji_m(j,i);
    // For each dimension
    for ( int f = 0; f < 2; f++)
      fpsCentroids(j,f) /= normalizationFactor;

  }
}
//----------------------------------------------------------
bool cluster::fuzzyClusterAlg::updateMembership(int k)
{
  // We are updating the membership of the data points based on the centroids
  // determined. This is determined using
  // u_ij = 1/( sum^C_k=1 ( ||x_i-c_j||/||x_i-c_k||)^(2/(m-1)))
  
  //fpsMat.Print();
  //fpsCentroids.Print();
  //fpsMembership.Print();
  
  int iNumClusters = k;
  TMatrixT<float> mNormOneXiMinusCj(fpsMat.GetNrows(),iNumClusters);
  std::vector< TMatrixT<float> > clusterCovarianceMats(iNumClusters);
  std::vector<double> clusterRadii(iNumClusters);

  // Determine the elements of u^m_ij
  TMatrixT<float> Uji_m(iNumClusters, fpsMat.GetNrows());
  // For each cluster
  for ( int j = 0; j < iNumClusters; j++)
    // For each hit
    for ( int i = 0; i < fpsMat.GetNrows(); i++)
      // Determine Uji_m
      Uji_m(j,i) = pow(fpsMembership(j,i), fFuzzifier); 



  TVectorT<float> fpsMat_row(2);
  TVectorT<float> fpsCentroids_row(2);
  TVectorT<float> fpsMat_col(2);
  TVectorT<float> fpsCentroids_col(2);
  TVectorT<float> fpsMat_row_t(2);
  TVectorT<float> fpsCentroids_row_t(2);
  TMatrixT<float> fpsDistances(iNumClusters,fpsMat.GetNrows());
  // Calculate the covariance matrix
  // For each clusters
  for ( int j = 0; j < iNumClusters; ++j){
    fpsCentroids_row = TMatrixFRow(fpsCentroids,j);
    TMatrixT<float> clusCovarianceMat(2,2);
    float Uji_m_sum = 0;
    //For each hit
    for ( int i = 0; i < fpsMat.GetNrows(); ++i){
      fpsMat_row = TMatrixFRow(fpsMat,i);
      TMatrixT<float> fpsMatMinusCent_col(1,2);
      TMatrixT<float> fpsMatMinusCent_row(2,1);
      fpsMatMinusCent_col(0,0)=fpsMat_row(0)-fpsCentroids_row(0);
      fpsMatMinusCent_col(0,1)=fpsMat_row(1)-fpsCentroids_row(1);
      fpsMatMinusCent_row(0,0)=fpsMat_row(0)-fpsCentroids_row(0);
      fpsMatMinusCent_row(1,0)=fpsMat_row(1)-fpsCentroids_row(1);
      //std::cout << "fpsMat" << std::endl;
      //fpsMat_row.Print();
      //std::cout << "fpsCentroids" << std::endl;
      //fpsCentroids_row.Print();
      clusCovarianceMat += Uji_m(j,i)*fpsMatMinusCent_row*fpsMatMinusCent_col;
      Uji_m_sum+=Uji_m(j,i);
    }
    //std::cout << "Uji_m_sum: " << Uji_m_sum << std::endl;
    //clusCovarianceMat.Print();
    fBeta = std::min((float)iNumClusters,fBeta+1);
    clusCovarianceMat=clusCovarianceMat*(1/Uji_m_sum);
    //clusCovarianceMat.Print();
    clusterCovarianceMats[j].ResizeTo(clusCovarianceMat);
    clusterCovarianceMats[j] = clusCovarianceMat;
    try
    {
      clusterRadii[j] = fBeta*pow(clusCovarianceMat.Determinant(),0.25)/(double)iNumClusters;
    }
    catch(...){
      mf::LogVerbatim("fuzzyCluster") << "updateMembership: Covariance matrix is singular 1";
      //clusterRadii[j] = 0;
      continue;
    }
    //std::cout << clusterRadii[j] << std::endl;
  }

  //std::cout << "New set" << std::endl;
  //for(auto clusterCovarianceMatsInt = clusterCovarianceMats.begin(); clusterCovarianceMatsInt != clusterCovarianceMats.end();  clusterCovarianceMatsInt++){
    //clusterCovarianceMatsInt->Print();
    //TMatrixDEigen clusCovarianceMatsEigen(*clusterCovarianceMatsInt);
    //clusCovarianceMatsEigen.GetEigenValues().Print();
  //}

  for ( int j = 0; j < iNumClusters; ++j){
    TMatrixT<float> clusCovarianceMatInv = clusterCovarianceMats[j];
    try
    {
      clusCovarianceMatInv.Invert();
    }
    catch(...){
      mf::LogVerbatim("fuzzyCluster") << "updateMembership: Covariance matrix is singular 2";
    }
    fpsCentroids_row = TMatrixFRow(fpsCentroids,j);
    for ( int i = 0; i < fpsMat.GetNrows(); ++i){
      fpsMat_row = TMatrixFRow(fpsMat,i);
      TMatrixT<float> fpsMatMinusCent_row(1,2);
      TMatrixT<float> fpsMatMinusCent_col(2,1);
      fpsMatMinusCent_row(0,0)=fpsMat_row(0)-fpsCentroids_row(0);
      fpsMatMinusCent_row(0,1)=fpsMat_row(1)-fpsCentroids_row(1);
      fpsMatMinusCent_col(0,0)=fpsMat_row(0)-fpsCentroids_row(0);
      fpsMatMinusCent_col(1,0)=fpsMat_row(1)-fpsCentroids_row(1);
      //fpsMatMinusCent_row.Print();
      //fpsMatMinusCent_col.Print();
      //clusCovarianceMatInv.Print();
      TMatrixT<float> tempDistanceSquared = (fpsMatMinusCent_row*(clusCovarianceMatInv*fpsMatMinusCent_col));
      //std::cout << clusCovarianceMatInv.Determinant() << " " << std::endl;
      //std::cout << tempDistanceSquared(0,0) << " " << sqrt(clusCovarianceMatInv.Determinant()) << std::endl;
      //std::cout << tempDistanceSquared(0,0)/sqrt(clusCovarianceMatInv.Determinant()) << " " << pow(clusterRadii[j],2) << std::endl;
      try
      {
        //fpsDistances(j,i) = std::sqrt(std::max((double)0,tempDistanceSquared(0,0)/sqrt(clusterCovarianceMats[j].Determinant()) - pow(clusterRadii[j],2)));
        fpsDistances(j,i) = std::sqrt(std::max((double)0,tempDistanceSquared(0,0)/sqrt(clusCovarianceMatInv.Determinant()) - pow(clusterRadii[j],2)));
      }
      catch(...){
        mf::LogVerbatim("fuzzyCluster") << "updateMembership: Covariance matrix is singular 2";
        fpsDistances(j,i) = 999999;
      }
      //std::cout << "fpsDistances(j,i): " << fpsDistances(j,i) << std::endl;
    }
  }

  ////Look at the eigenvalues
  //std::cout << "New set" << std::endl;
  //for ( int j = 0; j < iNumClusters; j++){
    //TVector eigenvalues(2);
    //TMatrix eigenvectors(2,2);
    //try
    //{
      //eigenvectors = clusterCovarianceMats[j].EigenVectors(eigenvalues);
    //}
    //catch(...){
      //continue;
    //}
    //eigenvalues.Print();
    //std::cout << std::sqrt(eigenvalues[0]) << " " << std::sqrt(eigenvalues[1]) << " " <<  std::sqrt(eigenvalues[0])/std::sqrt(eigenvalues[1])  << std::endl;
  //}

  //Determine the new elements of u_ij
  //std::cout << "fpsDistances: " << std::endl;
  //fpsDistances.Print();
  fpsNewMembership.ResizeTo(fpsMembership);
  float fCoeff;
  // For each hit
  for ( int i = 0; i < fpsMat.GetNrows(); i++){
    // Does the hit have nonzero distance to all clusters?
    int clusWithZeroDistance = 0;
    for ( int j = 0; j < iNumClusters; j++){
      if(fpsDistances(j,i) == 0)
        clusWithZeroDistance++;
    }
    // For each cluster
    for ( int j = 0; j < iNumClusters; j++){
      //if(clusWithZeroDistance==0){
        //std::cout << "no clusters with zero distance" << std::endl;
      //}
      //if(clusWithZeroDistance>0){
        //std::cout << "clusters with zero distance" << std::endl;
      //}
      if(clusWithZeroDistance == 0){ 
        fCoeff = 0;
        // For each cluster
        for ( int k = 0; k < iNumClusters; k++){
          //std::cout << fpsDistances(j,i) << " " << fpsDistances(k,i) << std::endl;
          fCoeff += pow (( fpsDistances(j,i)/fpsDistances(k,i)),2/(fFuzzifier - 1));
        }
        //std::cout << "fCoeff: " << fCoeff << std::endl;
        fpsNewMembership(j,i) = 1/fCoeff;
      }
      if(clusWithZeroDistance>0){
        if(fpsDistances(j,i) > 0){
          fpsNewMembership(j,i) = 0;
        }
        if(fpsDistances(j,i) == 0){
          fpsNewMembership(j,i) = 1/(double)clusWithZeroDistance;
          //fpsNewMembership(j,i) = 1/(double)iNumClusters;
        }
      }
    }
  }
  //fpsNewMembership.Print();

  if(!canStop()){
    fpsMembership = fpsNewMembership;
    return false;
  }
  return true;

}  




//----------------------------------------------------------------
/////////////////////////////////////////////////////////////////
// This is the algorithm that finds clusters:
//
//  Ben Carls' implementation of fuzzyClusterAlg as much like examples as possible
void cluster::fuzzyClusterAlg::run_fuzzy_cluster(std::vector<art::Ptr<recob::Hit> >& allhits) {

  int nMaxClusters = fMaxNumClusters;
  if(fpsMat.GetNrows() < nMaxClusters)
    nMaxClusters = fpsMat.GetNrows();
  
  std::vector<TMatrixF> fpsMembershipStore;
  fpsMembershipStore.clear();
  fpsMembershipStore.resize(nMaxClusters);

  fpointId_to_clusterId.resize(fps.size(), kNO_CLUSTER); // Not zero as before!
  fnoise.resize(fps.size(), false);
  fvisited.resize(fps.size(), false);

  art::ServiceHandle<geo::Geometry> geom;
  art::ServiceHandle<util::LArProperties> larprop;
  art::ServiceHandle<util::DetectorProperties> detprop;

  if(allhits.size()==0)
    return;
  
  //factor to make x and y scale the same units
  uint32_t channel    = allhits[0]->Wire()->RawDigit()->Channel();
  double   wirePitch  = geom->WirePitch(geom->View(channel));
  double   xyScale    = .001*larprop->DriftVelocity(larprop->Efield(),larprop->Temperature());
  xyScale            *= detprop->SamplingRate()/wirePitch;
  double   wire_dist  = wirePitch;
  double   tickToDist = larprop->DriftVelocity(larprop->Efield(),larprop->Temperature());
  tickToDist         *= 1.e-3 * detprop->SamplingRate(); // 1e-3 is conversion of 1/us to 1/ns
  
  float indcolscaling = 0.; // a parameter to account for the different 
        		    // characteristic hit width of induction and collection plane
  /// \todo: the collection plane's characteristic hit width's are, 
  /// \todo: on average, about 5 time samples wider than the induction plane's. 
  /// \todo: this is hard-coded for now.
  geo::SigType_t sigt = geom->SignalType(channel);
  if(sigt == geo::kInduction)
    indcolscaling = 0.;
  else
    indcolscaling = 1.;
  
  //fpsMat.Print();


  int k = nMaxClusters;
  if (k > fpsMat.GetNrows() || k <= 0)
    return;

  int i = 0;
  fpsMembership.ResizeTo(k, fps.size());
  fpsNewMembership.ResizeTo(k, fps.size());
  fpsMembershipStore[k-1].ResizeTo(k, fps.size());
  fpsCentroids.ResizeTo(k,2);


  /// Get the random number generator
  art::ServiceHandle<art::RandomNumberGenerator> rng;
  CLHEP::HepRandomEngine & engine = rng -> getEngine();
  CLHEP::RandFlat flat1(engine);

  //Randomize membership for each hit for fuzzy
  float normalizationFactor;
  for( int i = 0; i < fpsMat.GetNrows(); i++){
    normalizationFactor = 0;
    for( int j = 0; j < k; j++)
      normalizationFactor += fpsMembership(j,i) = flat1.fire(); //(rand() / (RAND_MAX + 0.0));
    for( int j = 0; j < k; j++)
      fpsMembership(j,i) /= normalizationFactor;
  }

  // Compute initial centroids
  computeCentroids(k);


  //for (size_t pid = 0; pid < fps.size(); pid++){
    //mf::LogInfo("fuzzyCluster") << pid ;
    //for (int l = 0; l < k; l++)
      //mf::LogInfo("fuzzyCluster") << l  << fpsMembership(l,pid) ;
  //} 

  // Run iterations of the clustering algorithm
  fBeta = 1;
  while(!updateMembership(k)){
    //std::cout << "k: " << k << std::endl;
    if(k == 1)
      break;
    if(mergeClusters()){
      k--;
      //fpsNewMembership.ResizeTo(fpsMembership);
      //fpsNewMembership = fpsMembership;
    }
    computeCentroids(k);
    if(i + 1 == nIterations)
      break;
    i++;
  }
  
  //for (size_t pid = 0; pid < fps.size(); pid++){
    //mf::LogInfo("fuzzyCluster") << pid ;
    //for (int l = 0; l < k; l++)
      //mf::LogInfo("fuzzyCluster") << l  << fpsMembership(l,pid) ;
  //} 


  fpsMembershipFinal.ResizeTo(k, fps.size());
  fpsMembershipFinal = fpsMembership; 

  int nClusters = 0;
  if(k > 0) nClusters = fpsMembershipFinal.GetNrows();
  unsigned int cid = fpsMembershipFinal.GetNrows();
  //int nClusters = iMinXBClusterNum;
  //unsigned int cid = iMinXBClusterNum;
  //mf::LogInfo("fuzzyCluster") << "Number of clusters found after merging: " << nClusters   ;
  //int nClusters = 2;
  //unsigned int cid = 2;

  //mf::LogInfo("fuzzyCluster") << iMinXBClusterNum  << nClusters ;
  //std::cout << "nClusters: " << nClusters  << std::endl;
  for (size_t pid = 0; pid < fps.size(); pid++){
    //mf::LogInfo("fuzzyCluster") << pid ;
    int iCluster = kNO_CLUSTER;
    float maxClusMembership = -1;
    // not already visited
    if (fpointId_to_clusterId[pid] == kNO_CLUSTER) {
      for (int i = 0; i <= nClusters-1; i++){
        //mf::LogInfo("fuzzyCluster") << i  << fpsMembershipStore[nClusters-1](i,pid) ;
        if ( fpsMembershipFinal(i,pid) > maxClusMembership ) {
          maxClusMembership = fpsMembershipFinal(i,pid); 
          iCluster = i;
        }
      }
    } // if (!visited
    fpointId_to_clusterId[pid] = iCluster;
  } // for

  // Run EndPointAlg over hits to see if 
  //mf::LogVerbatim("fuzzyCluster") << "New plane: " ;
  std::vector<unsigned int> corners;
  corners.clear();
  // nClustersTEmp is how many fuzzy clusters we originally found, it is not how many hough lines we found
  int nClustersTemp = nClusters;

 
  
  // Loop over clusters with the Hough line finder to break the clusters up further
  // list of lines
  std::vector<lineSlope> linesFound;
  if(nClustersTemp > 0)
    for (unsigned int i = 0; i <= (unsigned int)nClustersTemp-1; i++){
      fHBAlg.Transform(allhits, &fpointId_to_clusterId, i, &nClusters, corners, &linesFound);
    }

  // Determine the shower likeness of lines
  for(auto linesFoundItr = linesFound.begin(); linesFoundItr < linesFound.end(); linesFoundItr++){
    float totalBkgDist = 0;
    float totalBkgDistCharge = 0;
    int   totalBkg = 0;
    float fMaxDistance = 0.1;
    for(auto hitsItr = allhits.cbegin(); hitsItr != allhits.cend(); ++hitsItr){
      /// Veto the hit if it already belongs to a line
      //if(fpointId_to_clusterId->at(hitsItr-hits.cbegin()) != clusterId)
        //continue;
      float distance = (TMath::Abs((*hitsItr)->PeakTime()-linesFoundItr->clusterSlope*(double)((*hitsItr)->WireID().Wire)-linesFoundItr->clusterIntercept)/(std::sqrt(pow(xyScale*linesFoundItr->clusterSlope,2)+1)));
      /// Sum up background hits, use smart distance
      double peakTimePerpMin=-(1/linesFoundItr->clusterSlope)*(double)((*hitsItr)->WireID().Wire)+allhits[linesFoundItr->iMinWire]->PeakTime()+(1/linesFoundItr->clusterSlope)*(allhits[linesFoundItr->iMinWire]->WireID().Wire);
      double peakTimePerpMax=-(1/linesFoundItr->clusterSlope)*(double)((*hitsItr)->WireID().Wire)+allhits[linesFoundItr->iMaxWire]->PeakTime()+(1/linesFoundItr->clusterSlope)*(allhits[linesFoundItr->iMaxWire]->WireID().Wire);
      if(distance > 1*(fMaxDistance+(((*hitsItr)->EndTime()-(*hitsItr)->StartTime())/2.)+indcolscaling)
         && distance < 100*(fMaxDistance+(((*hitsItr)->EndTime()-(*hitsItr)->StartTime())/2.)+indcolscaling)){
        if((linesFoundItr->clusterSlope < 0 && (*hitsItr)->PeakTime() < peakTimePerpMin && (*hitsItr)->PeakTime() > peakTimePerpMax)
            || (linesFoundItr->clusterSlope > 0 && (*hitsItr)->PeakTime() > peakTimePerpMin && (*hitsItr)->PeakTime() < peakTimePerpMax)){
          totalBkgDist+=distance;
          totalBkgDistCharge+=distance/(*hitsItr)->Charge();
          totalBkg++;
        }
      }
    }/// end loop over hits
    linesFoundItr->showerLikeness = totalBkgDistCharge/(double)linesFoundItr->hits.size();
    std::cout << "showerLikeness: " << totalBkgDistCharge/(double)linesFoundItr->hits.size() << std::endl;
  }/// end loop over lines found

  if(fDoShowerHoughLineMerge)  mergeHoughLinesBySegment(0,&linesFound,xyScale,iMergeShower,wire_dist,tickToDist);
  if(fDoShowerHoughLineInterceptMerge)     mergeHoughLinesBySegment(0,&linesFound,xyScale,iMergeShowerIntercept,wire_dist,tickToDist);
  if(fDoChargeAsymAngleMerge)  mergeHoughLinesBySegment(0,&linesFound,xyScale,iMergeChargeAsymAngle,wire_dist,tickToDist);

  // Accumulate the merged lines
  std::map<int,mergedLines> mergedLinesMap; 
  for(auto linesFoundItr = linesFound.begin(); linesFoundItr != linesFound.end(); linesFoundItr++){
    if(!linesFoundItr->merged)
      continue;
    if(!mergedLinesMap.count(linesFoundItr->clusterNumber))
      mergedLinesMap.insert(std::make_pair(linesFoundItr->clusterNumber,mergedLines(linesFoundItr->totalQ,linesFoundItr->pMin0,linesFoundItr->pMin1,linesFoundItr->pMax0,linesFoundItr->pMax1,linesFoundItr->clusterNumber,linesFoundItr->showerLikeness)));
    else {
      // Check if new line has minimum that is lower
      if(linesFoundItr->pMin0 < mergedLinesMap.at(linesFoundItr->clusterNumber).pMin0){
        mergedLinesMap.at(linesFoundItr->clusterNumber).pMin0=linesFoundItr->pMin0;
        mergedLinesMap.at(linesFoundItr->clusterNumber).pMin1=linesFoundItr->pMin1;
      }
      if(linesFoundItr->pMax0 > mergedLinesMap.at(linesFoundItr->clusterNumber).pMax0){
        mergedLinesMap.at(linesFoundItr->clusterNumber).pMax0=linesFoundItr->pMax0;
        mergedLinesMap.at(linesFoundItr->clusterNumber).pMax1=linesFoundItr->pMax1;
      }
      mergedLinesMap.at(linesFoundItr->clusterNumber).totalQ+=linesFoundItr->totalQ;
      //mergedLinesMap.at(linesFoundItr->clusterNumber).showerLikeness+=linesFoundItr->showerLikeness;
      if(linesFoundItr->showerLikeness > mergedLinesMap.at(linesFoundItr->clusterNumber).showerLikeness){
        mergedLinesMap.at(linesFoundItr->clusterNumber).showerLikeness=linesFoundItr->showerLikeness;
      }
    }
  }

  /// Average the slopes and y-intercepts of Hough lines in shower like regions
  std::vector<showerLine> showerLines; 
  for(auto linesFoundItr = linesFound.begin(); linesFoundItr != linesFound.end(); linesFoundItr++){
    if(linesFoundItr->showerLikeness<fShowerLikenessCut && !linesFoundItr->merged)
      continue;
    if(linesFoundItr->merged)
      if(mergedLinesMap.at(linesFoundItr->clusterNumber).showerLikeness < fShowerLikenessCut)
        continue;
    int lineSize = linesFoundItr->hits.size();
    double lineSlope = linesFoundItr->clusterSlope;
    double lineIntercept = linesFoundItr->clusterIntercept;
    //std::cout << "line slopey: " << lineSlope << " intercepty: " << linesFoundItr->clusterIntercept << " line size: " << lineSize  << std::endl;
    // Check if showerLines contains a line with this cluster number already
    bool alreadyInserted=false;
    int insertedIndex=-999999;
    for(auto showerLinesItr = showerLines.begin(); showerLinesItr!=showerLines.end();showerLinesItr++){
      if(showerLinesItr->clusterNumber == linesFoundItr->clusterNumber){
        alreadyInserted=true;
        insertedIndex=showerLinesItr-showerLines.begin();
        break;
      }
    }
    if(!alreadyInserted){
      showerLines.push_back(showerLine(linesFoundItr->clusterSlope*lineSize,
              linesFoundItr->clusterIntercept*lineSize,
              linesFoundItr->clusterNumber,
              linesFoundItr->clusterNumber,
              linesFoundItr->hits.size(),
              linesFoundItr->iMinWire,
              linesFoundItr->iMaxWire));
      //std::cout << "min wire: " << hits[linesFoundItr->iMinWire]->WireID().Wire << std::endl;
      //std::cout << "max wire: " << hits[linesFoundItr->iMaxWire]->WireID().Wire << std::endl;
      continue;
    }
    if(allhits[linesFoundItr->iMinWire]->WireID().Wire < allhits[showerLines.at(insertedIndex).iMinWire]->WireID().Wire)
      showerLines.at(insertedIndex).iMinWire=linesFoundItr->iMinWire;
    if(allhits[linesFoundItr->iMaxWire]->WireID().Wire > allhits[showerLines.at(insertedIndex).iMaxWire]->WireID().Wire)
      showerLines.at(insertedIndex).iMaxWire=linesFoundItr->iMaxWire;
    showerLines.at(insertedIndex).slope+=lineSlope*lineSize;
    showerLines.at(insertedIndex).intercept+=lineIntercept*lineSize;
    showerLines.at(insertedIndex).lineSize+=lineSize;
    //std::cout << "min wire: " << hits[showerLines[insertedIndex].iMinWire]->WireID().Wire << std::endl;
    //std::cout << "max wire: " << hits[showerLines[insertedIndex].iMaxWire]->WireID().Wire << std::endl;
  }
  // Sort the averaged lines by size
  std::sort(showerLines.begin(),showerLines.end());
  std::reverse(showerLines.begin(),showerLines.end());

  /// Now actually do the averaging and find line direction
  std::vector<showerLine> showerLinesAverage; // 
  for(auto showerLinesItr = showerLines.begin(); showerLinesItr!= showerLines.end(); showerLinesItr++){

    if(showerLinesItr->showerMerged)
      continue;

    // If the line was already merged, skip it to prevent showers from splintering off
    if(showerLinesItr->clusterNumber != showerLinesItr->oldClusterNumber)
      continue;

    //showerLines.at(showerLinesItr->first).slope/=(showerLinesItr->second.lineSize);
    //showerLines.at(showerLinesItr->first).intercept/=(showerLinesItr->second.lineSize);
    //std::cout << showerLines.at(showerLinesItr-showerLines.begin()).slope << std::endl;
    //std::cout << showerLines.at(showerLinesItr-showerLines.begin()).intercept << std::endl;
    double averageSlope = showerLines.at(showerLinesItr-showerLines.begin()).slope/showerLinesItr->lineSize;
    double averageInt = showerLines.at(showerLinesItr-showerLines.begin()).intercept/showerLinesItr->lineSize;
    int midWire = allhits[showerLinesItr->iMinWire]->WireID().Wire/2 + 
      allhits[showerLinesItr->iMaxWire]->WireID().Wire/2 + 
      (allhits[showerLinesItr->iMinWire]->WireID().Wire & allhits[showerLinesItr->iMaxWire]->WireID().Wire & 1);
    double midPeakTime = averageSlope*midWire+averageInt;
    showerLinesAverage.push_back(showerLine(averageSlope,
              averageInt,
              showerLinesItr->clusterNumber,
              showerLinesItr->clusterNumber,
              showerLinesItr->lineSize,
              showerLinesItr->iMinWire,
              showerLinesItr->iMaxWire));
    mf::LogVerbatim("fuzzyClusterAlg") << "slope: " << averageSlope 
				       << " intercept: " << averageInt
				       << " min wire: " << allhits[showerLinesItr->iMinWire]->WireID().Wire
				       << " max wire: " << allhits[showerLinesItr->iMaxWire]->WireID().Wire
				       << " mid wire: " << midWire
				       << " midPeakTime: " << midPeakTime
				       << " line size: " << showerLinesItr->lineSize;

  }


  // Reassign the merged lines
  for(auto fpointId_to_clusterIdItr = fpointId_to_clusterId.begin(); fpointId_to_clusterIdItr != fpointId_to_clusterId.end(); ++fpointId_to_clusterIdItr){
    //if(*fpointId_to_clusterIdItr == clusterId)
      //continue;
    for(auto linesFoundItr = linesFound.begin(); linesFoundItr < linesFound.end(); linesFoundItr++){
      if(*fpointId_to_clusterIdItr == (unsigned int)linesFoundItr->oldClusterNumber)
        *fpointId_to_clusterIdItr = linesFoundItr->clusterNumber;
    }
  }


  // Find sizes of all merged lines combined
  // For linesFoundSizes, key is cluster number and size is the mapped value
  std::map<int,float> linesFoundSizes;
  for(auto linesFoundItr = linesFound.begin(); linesFoundItr < linesFound.end(); linesFoundItr++){
    if(!linesFoundSizes.count(linesFoundItr->clusterNumber))
      linesFoundSizes[linesFoundItr->clusterNumber] = std::sqrt( pow(linesFoundItr->pMin0-linesFoundItr->pMax0,2)+pow(linesFoundItr->pMin1-linesFoundItr->pMax1,2));
    else 
      linesFoundSizes[linesFoundItr->clusterNumber]+= std::sqrt( pow(linesFoundItr->pMin0-linesFoundItr->pMax0,2)+pow(linesFoundItr->pMin1-linesFoundItr->pMax1,2));
  }
  
  art::PtrVector<recob::Hit> unclusteredhits;
  std::vector<unsigned int> unclusteredhitsToallhits;
  int nDBClusters = 0;
  if(fDoFuzzyRemnantMerge){
    for(auto allhitsItr = allhits.cbegin(); allhitsItr != allhits.cend(); ++allhitsItr){
      bool unclustered = true;
      // nClusters is the number of fuzzy clusters we found, we only assign hits to lines here
      // if they are not already part of hough lines
      if(fpointId_to_clusterId.at(allhitsItr-allhits.begin()) >= (unsigned int) nClustersTemp){
        unclustered = false;
        continue;
      }
      double p0 = ((*allhitsItr)->Wire()->RawDigit()->Channel())*wire_dist;
      double p1 = (((*allhitsItr)->StartTime()+(*allhitsItr)->EndTime())/2.)*tickToDist;
      double minDistance = 999999;
      for(auto linesFoundItr = linesFound.begin(); linesFoundItr < linesFound.end(); linesFoundItr++){
        
        double distance = PointSegmentDistance( p0, p1, linesFoundItr->pMin0, linesFoundItr->pMin1, linesFoundItr->pMax0, linesFoundItr->pMax1);
        // Is the point behind or ahead of the line?
        //if(linesFoundItr->pMin0 > p0){
           //float linesFoundSlope = (linesFoundItr->pMax1 - linesFoundItr->pMin1)/(linesFoundItr->pMax0 - linesFoundItr->pMin0);
           //float pMinHitSlope = (p1 - linesFoundItr->pMin1)/(p0 - linesFoundItr->pMin0);
           //float slopeAngle = atan(std::abs((linesFoundSlope - pMinHitSlope)/(1 + linesFoundSlope*pMinHitSlope)))*(180/TMath::Pi());
           //if(distance < 10 && slopeAngle < 10){
             //fpointId_to_clusterId.at(allhitsItr-allhits.begin()) = linesFoundItr->clusterNumber;
             //unclustered = false;
             //break;
           //}
        //}
        //if (linesFoundItr->pMax0 < p0){
           //float linesFoundSlope = (linesFoundItr->pMax1 - linesFoundItr->pMin1)/(linesFoundItr->pMax0 - linesFoundItr->pMin0);
           //float pMaxHitSlope = (linesFoundItr->pMax1-p1)/(linesFoundItr->pMax0-p0);
           //float slopeAngle = atan(std::abs((linesFoundSlope - pMaxHitSlope)/(1 + linesFoundSlope*pMaxHitSlope)))*(180/TMath::Pi());
           //if(distance < 10 && slopeAngle < 10){
             //fpointId_to_clusterId.at(allhitsItr-allhits.begin()) = linesFoundItr->clusterNumber;
             //unclustered = false;
             //break;
           //}
        //}
        
        // If the line does not look showerlike, skip it
        if(linesFoundItr->showerLikeness<fShowerLikenessCut)
          continue;

        if(distance > fFuzzyRemnantMergeCutoff)
          continue;

        distance/=pow(linesFoundSizes[linesFoundItr->clusterNumber],1/4);
        if(distance < minDistance){
          fpointId_to_clusterId.at(allhitsItr-allhits.begin()) = linesFoundItr->clusterNumber;
          minDistance = distance;
          unclustered = false;
        }
      }
      if(unclustered){
        unclusteredhitsToallhits.push_back(allhitsItr-allhits.begin());
        unclusteredhits.push_back(*allhitsItr);
      }
      
    }

    // Setup DBSCAN for noise and extra hits
    // Start by getting the ChannelFilter
    filter::ChannelFilter chanFilt;
    fDBScan.InitScan(unclusteredhits, chanFilt.SetOfBadChannels());
    fDBScan.run_cluster();

    nDBClusters = fDBScan.fclusters.size();
    for(size_t j = 0; j < fDBScan.fpointId_to_clusterId.size(); ++j){	  
      if (fDBScan.fpointId_to_clusterId[j]== kNO_CLUSTER || fDBScan.fpointId_to_clusterId[j]==kNOISE_CLUSTER) {
      // This shouldn't happen...all points should be clasified by now!
        fpointId_to_clusterId.at(unclusteredhitsToallhits[j]) = kNOISE_CLUSTER;
      } 
      else {
        fpointId_to_clusterId.at(unclusteredhitsToallhits[j]) = fDBScan.fpointId_to_clusterId[j] + nClusters;
      }
    }
  }

  cid = nClusters + nDBClusters;
  
  //mf::LogInfo("fuzzyCluster") << "cid: " << cid ;

  //for(size_t j = 0; j < fpointId_to_clusterId.size(); ++j)
    //mf::LogInfo("fuzzyCluster") << "fpointId_to_clusterId[j]: " << fpointId_to_clusterId[j] << " j: " << j ;


  // Construct clusters, count noise, etc..
  int noise = 0;
  fclusters.resize(cid);
  for(size_t y = 0; y < fpointId_to_clusterId.size(); ++y){
    if (fpointId_to_clusterId[y] == kNO_CLUSTER) {
      // This shouldn't happen...all points should be clasified by now!
      mf::LogWarning("fuzzyCluster") << "Unclassified point!";
    } 
    else if (fpointId_to_clusterId[y]==kNOISE_CLUSTER) {
      ++noise;
    } 
    else {
      unsigned int c = fpointId_to_clusterId[y];
      if (c >= cid) {
	mf::LogWarning("fuzzyCluster") << "Point in cluster " << c 
			      << " when only " << cid 
			      << " clusters were found [0-" << cid-1
				 << "]";
      }
      fclusters[c].push_back(y);
    }
  }  
  mf::LogInfo("fuzzyCluster") << "DWM (R*-tree): Found " 
			   << cid << " clusters...";
  for (unsigned int c = 0; c < cid; ++c){
    mf::LogVerbatim("fuzzyCluster") << "\t" << "Cluster " << c << ":\t" 
			     << fclusters[c].size();
  }
  mf::LogVerbatim("fuzzyCluster") << "\t" << "...and " << noise << " noise points.";
}

//----------------------------------------------------------
bool cluster::fuzzyClusterAlg::mergeClusters()
{
  // If we only have one cluster, move on 
  if(fpsMembership.GetNrows() == 1)
    return false;

  //for (size_t pid = 0; pid < fps.size(); pid++){
    //mf::LogInfo("fuzzyCluster") << pid  << fpsMembershipCrisp[pid]  << fpsMat(pid,0)  << fpsMat(pid,1) ;
  //}
  
  float largestSij=0;
  float largesti=-1;
  float largestj=-1;
  for( int i = 0; i < fpsMembership.GetNrows()-1; i++){
    for( int j = i+1; j < fpsMembership.GetNrows(); j++){
      //std::cout << i << " " << j << std::endl;
      float sumMinUikUjk = 0;
      float sumUik = 0;
      float sumUjk = 0;
      for( int k = 0; k < fpsMembership.GetNcols(); k++){
        sumMinUikUjk+=std::min(fpsMembership(i,k),fpsMembership(j,k));
        sumUik += fpsMembership(i,k);
        sumUjk += fpsMembership(j,k);
      }
      float Sij = sumMinUikUjk/std::min(sumUik,sumUjk);
      if(Sij > largestSij){
        largestSij = Sij;
        largesti=i;
        largestj=j;
      }
    }
  }
  
  //std::cout << "largest Sij: " << largestSij << std::endl;
  if( largestSij > 1/(fpsMembership.GetNrows()-1)){
    mf::LogVerbatim("fuzzyCluster") << "You're going to Merge!";
  }
  else
    return false;

  if(largesti > -1 && largestj > -1){
    TMatrixF fpsMembershipTemp(fpsMembership.GetNrows()-1, fpsMat.GetNrows());
    TMatrixFRow(fpsMembership,largesti) += TMatrixFRow(fpsMembership,largestj);
    for( int j = 0; j < fpsMembership.GetNrows()-1; j++){
      //std::cout << j << " " << largestj << " " << largesti << std::endl;
      if(j < largestj)  TMatrixFRow(fpsMembershipTemp,j) = TMatrixFRow(fpsMembership,j);
      if(j >= largestj) TMatrixFRow(fpsMembershipTemp,j) = TMatrixFRow(fpsMembership,j+1);
    }
    fpsMembership.ResizeTo(fpsMembershipTemp);
    fpsMembership=fpsMembershipTemp;
  }
  
  return true;

}

//--------------------------------------------------------------------
// Merges based on the distance between line segments
void cluster::fuzzyClusterAlg::mergeHoughLinesBySegment(unsigned int clusIndexStart,
							std::vector<lineSlope> *linesFound,
							double xyScale,
							int mergeStyle,
							double wire_dist,
							double tickToDist)
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
  for(auto linesFoundToMergeItr = linesFound->begin(); linesFoundToMergeItr != linesFound->end(); linesFoundToMergeItr++){
    if(linesFound->at(clusIndexStart).clusterNumber == linesFoundToMergeItr->clusterNumber)
      continue;
    double segmentDistance = HoughLineDistance(linesFound->at(clusIndexStart).pMin0,linesFound->at(clusIndexStart).pMin1,
                                               linesFound->at(clusIndexStart).pMax0,linesFound->at(clusIndexStart).pMax1, 
      					       linesFoundToMergeItr->pMin0,linesFoundToMergeItr->pMin1,
                                               linesFoundToMergeItr->pMax0,linesFoundToMergeItr->pMax1);
    if( (segmentDistance<fHoughLineMergeCutoff && mergeStyle == iMergeNormal) 
        || (segmentDistance<fShowerHoughLineMergeCutoff && (mergeStyle == iMergeShower || mergeStyle == iMergeShowerIntercept))
        || (segmentDistance<fChargeAsymAngleCutoff && mergeStyle == iMergeChargeAsymAngle))
    {
      //std::cout << std::endl;
      //std::cout << linesFoundClusIndStItr->minWire << " " << linesFoundClusIndStItr->maxWire << std::endl;
      //std::cout << linesFoundToMergeItr->minWire << " " << linesFoundToMergeItr->maxWire << std::endl;
      //std::cout << segmentDistance << std::endl;
      toMerge.push_back(linesFoundToMergeItr-linesFound->begin());
      mergeSlope.push_back(linesFound->at(clusIndexStart).clusterSlope*xyScale);
    }
  }

  mergeTheta.resize(toMerge.size());

  // Find the angle between the slopes
  for(auto mergeThetaItr = mergeTheta.begin(); mergeThetaItr != mergeTheta.end(); mergeThetaItr++){
    double toMergeSlope = linesFound->at(toMerge[mergeThetaItr-mergeTheta.begin()]).clusterSlope*xyScale;
    mergeTheta[mergeThetaItr-mergeTheta.begin()] = atan(std::abs(( toMergeSlope - mergeSlope[mergeThetaItr-mergeTheta.begin()])/(1 + toMergeSlope*mergeSlope[mergeThetaItr-mergeTheta.begin()] )))*(180/TMath::Pi());
    //std::cout << std::endl;
    //std::cout << "toMergeSlope: " << toMergeSlope/xyScale<< " mergeSlope[clusIndexStart]: " << mergeSlope[mergeThetaItr-mergeTheta.begin()]/tickToDist << std::endl;
    //std::cout << "mergeTheta: " << mergeTheta[mergeThetaItr-mergeTheta.begin()] << std::endl;
  }

  // Perform the merge
  for(auto toMergeItr = toMerge.begin(); toMergeItr != toMerge.end(); toMergeItr++){
    if( (mergeTheta[toMergeItr-toMerge.begin()] < fHoughLineMergeAngle && mergeStyle == iMergeNormal) 
        || (mergeTheta[toMergeItr-toMerge.begin()] < fShowerHoughLineMergeAngle && mergeStyle == iMergeShower)
        || mergeStyle == iMergeShowerIntercept
        || mergeStyle == iMergeChargeAsymAngle){
     

      // First check averages of charge and sigma charge for hits in lines closest to each other
      int closestToMerge=-1;
      int closestClusIndexStart=-1;
      double closestDistance=999999;
      for (auto toMergeHitItr = linesFound->at(*toMergeItr).hits.begin(); toMergeHitItr != linesFound->at(*toMergeItr).hits.end(); toMergeHitItr++) {
        for (auto clusIndStHitItr = linesFound->at(clusIndexStart).hits.begin(); clusIndStHitItr != linesFound->at(clusIndexStart).hits.end(); clusIndStHitItr++) {
          //double distance = std::sqrt(pow(clusIndStHitItr->first-(*toMergeHitItr).first,2)+
                    //pow(clusIndStHitItr->second-toMergeHitItr->second,2));
          double distance = DistanceBetweenHits(*clusIndStHitItr,
                                                *toMergeHitItr,
                                                wire_dist,
                                                tickToDist);
          if(distance < closestDistance){
            closestDistance = distance;
            closestToMerge=toMergeHitItr-linesFound->at(*toMergeItr).hits.begin();
            closestClusIndexStart=clusIndStHitItr-linesFound->at(clusIndexStart).hits.begin();
          }
        }
      }

      // Find up to 9 more points closest to closestToMerge on the toMerge[i] line
      // check if it's closer, insert, delete
      std::vector<std::pair<int,double> > closestToMergeDist;
      for (auto toMergeHitItr = linesFound->at(*toMergeItr).hits.begin(); toMergeHitItr != linesFound->at(*toMergeItr).hits.end(); toMergeHitItr++) {
        if(closestToMerge==toMergeHitItr-linesFound->at(*toMergeItr).hits.begin())
          continue;
          double distance = DistanceBetweenHits(linesFound->at(clusIndexStart).hits[closestClusIndexStart],
                                                *toMergeHitItr,
                                                wire_dist,
                                                tickToDist);

        bool foundCloser = false;
        for(auto closestToMergeDistItr = closestToMergeDist.begin(); closestToMergeDistItr != closestToMergeDist.end(); closestToMergeDistItr++) {
          if(closestToMergeDistItr->second > distance){
            foundCloser = true;
            break;
          }
        }
        if(foundCloser 
            || closestToMergeDist.size() < linesFound->at(*toMergeItr).hits.size()-1
            || closestToMergeDist.size() < 9){
          closestToMergeDist.push_back(std::make_pair(toMergeHitItr-linesFound->at(*toMergeItr).hits.begin(),distance));
          std::sort(closestToMergeDist.begin(), closestToMergeDist.end(), boost::bind(&std::pair<int,double>::second,_1) < boost::bind(&std::pair<int,double>::second,_2));
        }
        if(closestToMergeDist.size() > linesFound->at(*toMergeItr).hits.size()-1 ||
          closestToMergeDist.size() > 9)
          closestToMergeDist.erase(closestToMergeDist.end());
      }
      //for(auto closestToMergeDistItr = closestToMergeDist.begin(); closestToMergeDistItr != closestToMergeDist.end();
        //closestToMergeDistItr++) 
        //std::cout << closestToMergeDistItr->first << " " << closestToMergeDistItr->second << std::endl;



      // Find up to 9 more points closest to closestToMerge on the clusIndexStart line
      std::vector<std::pair<int,double> > closestClusIndexStartDist;
      for (auto clusIndexStartHitItr = linesFound->at(clusIndexStart).hits.begin(); clusIndexStartHitItr != linesFound->at(clusIndexStart).hits.end(); clusIndexStartHitItr++) {
        if(closestClusIndexStart==clusIndexStartHitItr-linesFound->at(clusIndexStart).hits.begin())
          continue;

        double distance = DistanceBetweenHits(*clusIndexStartHitItr,
                                              linesFound->at(*toMergeItr).hits[closestToMerge],
                                              wire_dist,
                                              tickToDist);

        bool foundCloser = false;
        for(auto closestClusIndexStartDistItr = closestClusIndexStartDist.begin(); closestClusIndexStartDistItr != closestClusIndexStartDist.end(); closestClusIndexStartDistItr++) {
          if(closestClusIndexStartDistItr->second > distance){
            foundCloser = true;
            break;
          }
        }
        if(foundCloser 
            || closestClusIndexStartDist.size() < linesFound->at(clusIndexStart).hits.size()-1
            || closestClusIndexStartDist.size() < 9){
          closestClusIndexStartDist.push_back(std::make_pair(clusIndexStartHitItr-linesFound->at(clusIndexStart).hits.begin(),distance));
          std::sort(closestClusIndexStartDist.begin(), closestClusIndexStartDist.end(), boost::bind(&std::pair<int,double>::second,_1) < boost::bind(&std::pair<int,double>::second,_2));
        }
        if(closestClusIndexStartDist.size() > linesFound->at(clusIndexStart).hits.size()-1 ||
          closestClusIndexStartDist.size() > 9)
          closestClusIndexStartDist.erase(closestClusIndexStartDist.end());
     }



      double toMergeAveCharge = linesFound->at(*toMergeItr).hits[closestToMerge]->Charge();
      double toMergeAveSigmaCharge = linesFound->at(*toMergeItr).hits[closestToMerge]->SigmaCharge();
      for(auto closestToMergeDistItr = closestToMergeDist.begin(); closestToMergeDistItr != closestToMergeDist.end(); closestToMergeDistItr++) {
        toMergeAveCharge+=linesFound->at(*toMergeItr).hits[closestToMergeDistItr->first]->Charge();
        toMergeAveSigmaCharge+=linesFound->at(*toMergeItr).hits[closestToMergeDistItr->first]->SigmaCharge();
      }
      double clusIndexStartAveCharge = linesFound->at(clusIndexStart).hits[closestClusIndexStart]->Charge();
      double clusIndexStartAveSigmaCharge = linesFound->at(clusIndexStart).hits[closestClusIndexStart]->SigmaCharge();
      for(auto closestClusIndexStartDistItr = closestClusIndexStartDist.begin(); closestClusIndexStartDistItr != closestClusIndexStartDist.end(); closestClusIndexStartDistItr++) {
        clusIndexStartAveCharge+=linesFound->at(clusIndexStart).hits[closestClusIndexStartDistItr->first]->Charge();
        clusIndexStartAveSigmaCharge+=linesFound->at(clusIndexStart).hits[closestClusIndexStartDistItr->first]->SigmaCharge();
      }



      double chargeAsymmetry = std::abs(toMergeAveCharge-clusIndexStartAveCharge)/(toMergeAveCharge+clusIndexStartAveCharge);
      double sigmaChargeAsymmetry = std::abs(toMergeAveSigmaCharge-clusIndexStartAveSigmaCharge)/(toMergeAveSigmaCharge+clusIndexStartAveSigmaCharge);
      double chargeAsymmetrySinAngle = chargeAsymmetry*pow(std::fabs(sin(mergeTheta[toMergeItr-toMerge.begin()]*TMath::Pi()/180)),1);
      double sigmaChargeAsymmetrySinAngle = sigmaChargeAsymmetry*pow(std::fabs(sin(mergeTheta[toMergeItr-toMerge.begin()]*TMath::Pi()/180)),1);

      //std::cout << std::endl;
      //std::cout << chargeAsymmetry*pow(std::fabs(sin(mergeTheta[toMergeItr-toMerge.begin()]*TMath::Pi()/180)),1) << std::endl;
      //std::cout << sigmaChargeAsymmetry*pow(std::fabs(sin(mergeTheta[toMergeItr-toMerge.begin()]*TMath::Pi()/180)),1) << std::endl;

      if(chargeAsymmetrySinAngle > fChargeAsymAngleCut &&
          mergeStyle == iMergeChargeAsymAngle)
        continue;

      if(sigmaChargeAsymmetrySinAngle > fSigmaChargeAsymAngleCut  &&
          mergeStyle == iMergeChargeAsymAngle)
        continue;

      //double lineLengthAsymm = std::abs((double)linesFound->at(clusIndexStart).hits.size() - (double)linesFound->at(*toMergeItr).hits.size())/((double)linesFound->at(clusIndexStart).hits.size() + (double)linesFound->at(*toMergeItr).hits.size());

      //std::cout << "Length Asymm: " << lineLengthAsymm << std::endl;
      //std::cout << linesFound->at(clusIndexStart).hits.size() << std::endl;

      //if(lineLengthAsymm > 0.75) 
        //continue; 

      // Veto the merge if the lines are not colinear 
      if(mergeStyle == iMergeNormal || mergeStyle == iMergeChargeAsymAngle) {
        // Find where the lines are closest
  
        //distance between two segments in the plane:
        //  one segment is (x11, y11) to (x12, y12) or (p0MinLine1, p1MinLine1) to (p0MaxLine1, p1MaxLine1)
        //  the other is   (x21, y21) to (x22, y22) or (p0MinLine2, p1MinLine2) to (p0MaxLine2, p1MaxLine2)
        double x11 = linesFound->at(*toMergeItr).pMin0; 
        double y11 = linesFound->at(*toMergeItr).pMin1; 
        double x12 = linesFound->at(*toMergeItr).pMax0; 
        double y12 = linesFound->at(*toMergeItr).pMax1; 
        double x21 = linesFound->at(clusIndexStart).pMin0; 
        double y21 = linesFound->at(clusIndexStart).pMin1; 
        double x22 = linesFound->at(clusIndexStart).pMax0; 
        double y22 = linesFound->at(clusIndexStart).pMax1; 
        std::vector<double> distances;

        // Compare toMergerItr min with clusIndexStart min, if this is the min distance, lines are not colinear, merge is vetoed
        distances.push_back(std::sqrt(pow(x11-x21,2) + pow(y11-y21,2)));
        // Compare toMergerItr min with clusIndexStart max
        distances.push_back(std::sqrt(pow(x11-x22,2) + pow(y11-y22,2)));
        // Compare toMergerItr max with clusIndexStart min
        distances.push_back(std::sqrt(pow(x12-x21,2) + pow(y12-y21,2)));
        // Compare toMergerItr max with clusIndexStart max, if this is the min distance, lines are not colinear, merge is vetoed
        distances.push_back(std::sqrt(pow(x12-x22,2) + pow(y12-y22,2)));

        double minDistance = 999999; 
        int minDistanceIndex = -1;
        for(unsigned int j = 0; j < distances.size(); j++){
          if (distances[j] < minDistance){
            minDistance = distances[j];
            minDistanceIndex = j;
          }
        }

        if(minDistanceIndex  == 0 || minDistanceIndex  == 3)
          continue;

      }
      //std::cout << linesFound->at(*toMergeItr).clusterSlope << " " << linesFound->at(clusIndexStart).clusterSlope << std::endl;



      // Check if both lines is in region that looks showerlike
      // Or merge if the distance between the lines is zero and one looks showerlike
//       double segmentDistance = HoughLineDistance(linesFound->at(clusIndexStart).pMin0,
// 						 linesFound->at(clusIndexStart).pMin1,
//                                                  linesFound->at(clusIndexStart).pMax0,
// 						 linesFound->at(clusIndexStart).pMax1, 
//                                                  linesFound->at(*toMergeItr).pMin0,
// 						 linesFound->at(*toMergeItr).pMin1,
//                                                  linesFound->at(*toMergeItr).pMax0,
// 						 linesFound->at(*toMergeItr).pMax1);
      //std::cout << "segmentDistance: " << segmentDistance << std::endl;
      //std::cout << linesFound->at(*toMergeItr).showerLikeness << " " << linesFound->at(clusIndexStart).showerLikeness << std::endl;

      
     
      // If doing a shower merge, only merge if the Hough lines look showerlike
      if(mergeStyle == iMergeShower){
        if(!(linesFound->at(*toMergeItr).showerLikeness>fShowerLikenessCut) || !(linesFound->at(clusIndexStart).showerLikeness>fShowerLikenessCut))
          continue;
      }
      
      
      //// If not doing a shower merge, only merge if one of the Hough lines doesn't look showerlike
      //if(mergeStyle == iMergeChargeAsymAngle){
        //if(linesFound->at(*toMergeItr).showerLikeness>fShowerLikenessCut && linesFound->at(clusIndexStart).showerLikeness>fShowerLikenessCut)
          //continue;
      //}

      
      if(mergeStyle == iMergeShowerIntercept){
        if((linesFound->at(*toMergeItr).showerLikeness<fShowerLikenessCut) && (linesFound->at(clusIndexStart).showerLikeness<fShowerLikenessCut))
          continue;
        double x11 = linesFound->at(*toMergeItr).pMin0; 
        double y11 = linesFound->at(*toMergeItr).pMin1; 
        double x12 = linesFound->at(*toMergeItr).pMax0; 
        double y12 = linesFound->at(*toMergeItr).pMax1; 
        double x21 = linesFound->at(clusIndexStart).pMin0; 
        double y21 = linesFound->at(clusIndexStart).pMin1; 
        double x22 = linesFound->at(clusIndexStart).pMax0; 
        double y22 = linesFound->at(clusIndexStart).pMax1; 
        if(HoughLineIntersect(x11, y11, x12, y12, x21, y21, x22, y22) == 0)
          continue;
      }

      //std::cout << "Merging" << std::endl;
      lineMerged = true;
      linesFound->at(clusIndexStart).merged = true;
      linesFound->at(*toMergeItr).merged = true;

      // For loop over all lines found to reassign lines to clusIndexStart that already belonged to toMerge 
      for(auto linesFoundItr = linesFound->begin(); linesFoundItr != linesFound->end(); linesFoundItr++){
        if((unsigned int)(*toMergeItr) == linesFoundItr-linesFound->begin())
          continue;

        if(linesFoundItr->clusterNumber == linesFound->at(*toMergeItr).clusterNumber){
          linesFoundItr->clusterNumber = linesFound->at(clusIndexStart).clusterNumber;
          //std::cout << "linesFoundItr slope: " << linesFoundItr->clusterSlope << " clusIndexStart slope: " << linesFound->at(clusIndexStart).clusterSlope << std::endl;

        }
      }
      linesFound->at(*toMergeItr).clusterNumber = linesFound->at(clusIndexStart).clusterNumber;
      //std::cout << "main linesFoundItr slope: " << linesFound->at(*toMergeItr).clusterSlope << " clusIndexStart slope: " << linesFound->at(clusIndexStart).clusterSlope << std::endl;
     
    }  
  }

                                                
  if(lineMerged)
    mergeHoughLinesBySegment(clusIndexStart,linesFound,xyScale,mergeStyle,wire_dist,tickToDist);
  else
    mergeHoughLinesBySegment(clusIndexStart+1,linesFound,xyScale,mergeStyle,wire_dist,tickToDist);
  
  return;

}


//------------------------------------------------------------------------------
double cluster::fuzzyClusterAlg::DistanceBetweenHits(art::Ptr<recob::Hit> hit0,
						     art::Ptr<recob::Hit> hit1,
						     double wire_dist,
						     double tickToDist)
{
  double pHit0[2];
  pHit0[0] = (hit0->Wire()->RawDigit()->Channel())*wire_dist;
  pHit0[1] = ((hit0->StartTime()+hit0->EndTime())/2.)*tickToDist;
  double pHit1[2];
  pHit1[0] = (hit1->Wire()->RawDigit()->Channel())*wire_dist;
  pHit1[1] = ((hit1->StartTime()+hit1->EndTime())/2.)*tickToDist;

  return std::sqrt( pow(pHit0[0] - pHit1[0],2) + pow(pHit0[1] - pHit1[1],2));

}

//------------------------------------------------------------------------------
double cluster::fuzzyClusterAlg::HoughLineDistance(double p0MinLine1, 
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
bool cluster::fuzzyClusterAlg::HoughLineIntersect(double x11,
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
double cluster::fuzzyClusterAlg::PointSegmentDistance(double px,
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
