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


//Framework includes:
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "ClusterFinder/fuzzyClusterAlg.h"
#include "RecoBase/Hit.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>

namespace cluster{
  const unsigned int kNO_CLUSTER    = UINT_MAX;
  const unsigned int kNOISE_CLUSTER = UINT_MAX-1;
}

//----------------------------------------------------------
// fuzzyClusterAlg stuff
//----------------------------------------------------------
cluster::fuzzyClusterAlg::fuzzyClusterAlg(fhicl::ParameterSet const& pset) 
   : fHBAlg(pset.get< fhicl::ParameterSet >("HoughBaseAlg"))
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
  fEps            = p.get< double >("eps"   );
  fEps2           = p.get< double >("epstwo");
  fMinPts         = p.get< int    >("minPts");
  fClusterMethod  = p.get< int    >("Method");
  fDistanceMetric = p.get< int    >("Metric");
  fFuzzifier      = p.get< double >("Fuzzifier");
  nMaxClusters    = p.get< int    >("MaxClusters");
  nIterations     = p.get< int    >("Iterations");
  fMergeCutoff    = p.get< double >("MergeCutoff");
  fHBAlg.reconfigure(p.get< fhicl::ParameterSet >("HoughBaseAlg"));
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
  
  // Collect the bad wire list into a useful form
  if (fClusterMethod) { // Using the R*-tree
    fBadWireSum.resize(fGeom->Nchannels());
    unsigned int count=0;
    for (unsigned int i=0; i<fBadWireSum.size(); ++i) {
      count += fBadChannels.count(i);
      fBadWireSum[i] = count;
    }
  }

  // Collect the hits in a useful form,
  // and take note of the maximum time width
  fMaxWidth=0.0;
  //fpsMat = TMatrixT<double>(allhits.size(),2);
  //fpsMembership = TMatrixT<double>(iNumClusters, allhits.size());
  fpsMat.ResizeTo(allhits.size(),2);
  for (unsigned int j = 0; j < allhits.size(); ++j){
    int dims = 3;//our point is defined by 3 elements:wire#,center of the hit, and the hit width
    std::vector<double> p(dims);
        
    double tickToDist = larp->DriftVelocity(larp->Efield(),larp->Temperature());
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

  fpsCentroids = fpsMembership*fpsMat;

  //float sumUk[iNumClusters] = {0};
  std::vector<float> sumUk;
  sumUk.resize(iNumClusters);

  // Zero sumUk
  //for ( int i = 0; i < iNumClusters; i++){
    //mf::LogInfo("fuzzyCluster") << sumUk[i] ;
    //sumUk[i]=0;
  //} 

  for ( int i = 0; i < iNumClusters; i++)
    for ( int j = 0; j < fpsMat.GetNrows(); j++)
      sumUk[i] += fpsMembership(i,j);

  for ( int i = 0; i < iNumClusters; i++)
    for ( int f = 0; f < 2; f++)
      fpsCentroids(i,f) /= sumUk[i];

}
//----------------------------------------------------------
void cluster::fuzzyClusterAlg::computeCentroids2(int k)
{
  // Centroids are defined by c_j = (sum^N_i=1 u^m_ij*x_i)/(sum^N_i=1 u^m_ij)

  int iNumClusters = k;
  TMatrixT<double> Uji_m(iNumClusters, fpsMat.GetNrows());
  double normalizationFactor;

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

  int iNumClusters = k;
  TMatrixT<double> mNormOneXiMinusCj(fpsMat.GetNrows(),iNumClusters);
  std::vector<TMatrixT<double> > clusterCovarianceMats;


  // Determine the elements of u^m_ij
  TMatrixT<double> Uji_m(iNumClusters, fpsMat.GetNrows());
  // For each cluster
  for ( int j = 0; j < iNumClusters; j++)
    // For each hit
    for ( int i = 0; i < fpsMat.GetNrows(); i++)
      // Determine Uji_m
      Uji_m(j,i) = pow ( fpsMembership(j,i), fFuzzifier); 


  TVectorD fpsMat_row(2);
  TVectorD fpsCentroids_row(2);
  TVectorD fpsMat_col(2);
  TVectorD fpsCentroids_col(2);
  TVectorD fpsMat_row_t(2);
  TVectorD fpsCentroids_row_t(2);
  TMatrixT<double> fpsDistances(iNumClusters,fpsMat.GetNrows());
  // Calculate the covariance matrix
  // For each clusters
  for ( int j = 0; j < iNumClusters; j++){
    fpsCentroids_row = TMatrixDRow(fpsCentroids,j);
    TMatrixT<double> clusCovarianceMat(2,2);
    double Uji_m_sum = 0;
    //For each hit
    for ( int i = 0; i < fpsMat.GetNrows(); i++){
      fpsMat_row = TMatrixDRow(fpsMat,i);
      TMatrixT<double> fpsMatMinusCent_col(1,2);
      TMatrixT<double> fpsMatMinusCent_row(2,1);
      fpsMatMinusCent_col(0,0)=fpsMat_row(0)-fpsCentroids_row(0);
      fpsMatMinusCent_col(0,1)=fpsMat_row(1)-fpsCentroids_row(1);
      fpsMatMinusCent_row(0,0)=fpsMat_row(0)-fpsCentroids_row(0);
      fpsMatMinusCent_row(1,0)=fpsMat_row(1)-fpsCentroids_row(1);
      clusCovarianceMat += Uji_m(j,i)*fpsMatMinusCent_row*fpsMatMinusCent_col;
      Uji_m_sum+=Uji_m(j,i);
    }
    clusCovarianceMat=clusCovarianceMat*(1/Uji_m_sum);
    //clusCovarianceMat.Print();
    clusterCovarianceMats.push_back(clusCovarianceMat);
  }

  double rho = 1;
  for ( int j = 0; j < iNumClusters; j++){
    fpsCentroids_row = TMatrixDRow(fpsCentroids,j);
    for ( int i = 0; i < fpsMat.GetNrows(); i++){
      fpsMat_row = TMatrixDRow(fpsMat,i);
      TMatrixT<double> fpsMatMinusCent_col(1,2);
      TMatrixT<double> fpsMatMinusCent_row(2,1);
      fpsMatMinusCent_col(0,0)=fpsMat_row(0)-fpsCentroids_row(0);
      fpsMatMinusCent_col(0,1)=fpsMat_row(1)-fpsCentroids_row(1);
      fpsMatMinusCent_row(0,0)=fpsMat_row(0)-fpsCentroids_row(0);
      fpsMatMinusCent_row(1,0)=fpsMat_row(1)-fpsCentroids_row(1);
      //TMatrixT<double> clusCovarianceMatInv = clusterCovarianceMats.at(j).Invert();
      TMatrixT<double> clusCovarianceMatInv = clusterCovarianceMats.at(j);
      //clusCovarianceMatInv.Print();
      clusCovarianceMatInv.Invert();
      //clusCovarianceMatInv.Print();
      //TMatrixT<double> tempDistanceSquared = rho*std::sqrt(clusterCovarianceMats.at(j).Determinant())*(fpsMatMinusCent_col*(clusCovarianceMatInv*fpsMatMinusCent_row));
      TMatrixT<double> tempDistanceSquared = rho*pow(clusterCovarianceMats.at(j).Determinant(),0.75)*(fpsMatMinusCent_col*(clusCovarianceMatInv*fpsMatMinusCent_row));
      fpsDistances(j,i) = std::sqrt(tempDistanceSquared(0,0));
      //tempDistanceSquared.Print();
    }
  }
  //fpsDistances.Print();

  //Determine the new elements of u_ij
  double fCoeff;
  // For each hit
  for ( int i = 0; i < fpsMat.GetNrows(); i++)
    // For each cluster
    for ( int j = 0; j < iNumClusters; j++){
      fCoeff = 0;
      // For each cluster
      for ( int k = 0; k < iNumClusters; k++)
        fCoeff += pow (( fpsDistances(j,i)/fpsDistances(k,i)),2/(fFuzzifier - 1));
      fpsNewMembership(j,i) = 1/fCoeff;
    }



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

  int iMinXBClusterNum = 0;
  std::vector<float> fXieBeniIndices;
  fXieBeniIndices.clear();
  fXieBeniIndices.resize(nMaxClusters);
  std::vector<TMatrixT<double> > fpsMembershipStore;
  fpsMembershipStore.clear();
  fpsMembershipStore.resize(nMaxClusters);

  fpointId_to_clusterId.resize(fps.size(), kNO_CLUSTER); // Not zero as before!
  fnoise.resize(fps.size(), false);
  fvisited.resize(fps.size(), false);

  //fpsMat.Print();

  for( int k = 2; k <= nMaxClusters; k++){ 

    if (k > fpsMat.GetNrows())
      break;

    int i = 0;
    fpsMembership.ResizeTo(k, fps.size());
    fpsNewMembership.ResizeTo(k, fps.size());
    fpsMembershipStore[k-1].ResizeTo(k, fps.size());
    fpsCentroids.ResizeTo(k,2);

    //Randomize membership for each hit for fuzzy
    double normalizationFactor;
    for( int i = 0; i < fpsMat.GetNrows(); i++){
      normalizationFactor = 0;
      for( int j = 0; j < k; j++)
        normalizationFactor += fpsMembership(j,i) = (rand() / (RAND_MAX + 0.0));
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
    while(!updateMembership(k) && i++ < nIterations)
      computeCentroids2(k);
  
    //for (size_t pid = 0; pid < fps.size(); pid++){
      //mf::LogInfo("fuzzyCluster") << pid ;
      //for (int l = 0; l < k; l++)
        //mf::LogInfo("fuzzyCluster") << l  << fpsMembership(l,pid) ;
    //} 


    //// Determine Xie-Beni index to quantify how well clustering worked
    float fXieBeniNumer = 0; 
    float fXieBeniDenom; 

    // Begin with numerator
    TMatrixT<double> Uji_m(k, fpsMat.GetNrows());
    // Determine the elements of u^m_ij
    // For each cluster
    for ( int j = 0; j < k; j++)
      // For each hit
      for ( int i = 0; i < fpsMat.GetNrows(); i++)
        // Determine Uji_m
        Uji_m(j,i) = pow ( fpsMembership(j,i), fFuzzifier); 

    TMatrixT<double> mNormOneXiMinusCj(fpsMat.GetNrows(),k);

    // Zero the matrix that will store all the norm values
    for ( int i = 0; i < fpsMat.GetNrows(); i++)
      for ( int j = 0; j < k; j++)
        mNormOneXiMinusCj(i,j) = 0;

    // Calculate the norm matrix x_i - c_j
    for ( int j = 0; j < k; j++)
      for ( int i = 0; i < fpsMat.GetNrows(); i++)
        for ( int f = 0; f < 2; f++)
          // x_i - c_j                      x_i               c_j
          mNormOneXiMinusCj(i,j) += pow( fpsMat(i,f) - fpsCentroids(j,f),2);
          
    // Find the square roots of the elements of the norm matrix
    for ( int i = 0; i < fpsMat.GetNrows(); i++)
      for ( int j = 0; j < k; j++)
        mNormOneXiMinusCj(i,j) = std::sqrt(mNormOneXiMinusCj(i,j));


    for ( int l = 0; l < fpsMat.GetNrows(); l++) 
      for ( int j = 0; j < k; j++)
        fXieBeniNumer += Uji_m(j,l)*mNormOneXiMinusCj(l,j);

    //std::cout << "centroids!" << std::endl;
    //fpsCentroids.Print();
    //fpsMembership.Print();

    // Continue with denominator
    float fMinCentDist = 999999999;
    //float fMinCentDist = 100;
    for ( int j = 0; j < k - 1; j++)
      for ( int l = j + 1; l < k; l++){
        //mf::LogInfo("fuzzyCluster") << fpsCentroids(j,0) << fpsCentroids(l,0) <<fpsCentroids(j,1)<<" "<<fpsCentroids(l,1);
        float fCentDist = std::sqrt( pow(fpsCentroids(j,0)-fpsCentroids(l,0),2)+pow(fpsCentroids(j,1)-fpsCentroids(l,1),2));
        if(fCentDist < fMinCentDist)
          fMinCentDist = fCentDist;
      }
    fXieBeniDenom = fpsMat.GetNrows()*fMinCentDist;



    float fXieBeniIndex = fXieBeniNumer/fXieBeniDenom; 
    fXieBeniIndices[k-1]=fXieBeniIndex;
    fpsMembershipStore[k-1]=fpsMembership;
    mf::LogVerbatim("fuzzyCluster") << "Number of clusters: " << k ;
    mf::LogVerbatim("fuzzyCluster") << "Xie-Beni numerator: " << fXieBeniNumer ;
    mf::LogVerbatim("fuzzyCluster") << "Xie-Beni denominator: " << fXieBeniDenom ;
    mf::LogVerbatim("fuzzyCluster") << "Xie-Beni index: " << fXieBeniIndex ;
  }


  // Find the smallest Xie-Beni index
  double fXieBeniIndexMin=999999999;
  for( unsigned int j = 1; j < fXieBeniIndices.size(); j++){ 
    //mf::LogInfo("fuzzyCluster") << fXieBeniIndices[j]  << fXieBeniIndexMin ;
    if (j+1 > (unsigned int)fpsMat.GetNrows())
      break;
    if(fXieBeniIndices[j] < fXieBeniIndexMin){
      fXieBeniIndexMin = fXieBeniIndices[j];
      iMinXBClusterNum = j+1;    
    }
  }


  //mf::LogInfo("fuzzyCluster") << "Number of clusters initially found: " << iMinXBClusterNum   ;

  //iMinXBClusterNum = nMaxClusters;
  if(iMinXBClusterNum != 0){
    // Check if any clusters can be merged, most likely yes
    fpsMembershipFinal.ResizeTo(iMinXBClusterNum, fps.size());
    fpsMembershipFinal = fpsMembershipStore[iMinXBClusterNum-1]; 
    mergeClusters(0); 
  }










  // End clustering algorithm



 

  int nClusters = 0;
  if(iMinXBClusterNum > 0) nClusters = fpsMembershipFinal.GetNrows();
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
  mf::LogVerbatim("fuzzyCluster") << "New plane: " ;
  std::vector<unsigned int> corners;
  corners.clear();
  int nClustersTemp = nClusters;

  
  // Loop over clusters with the Hough line finder to break the clusters up further
  if(nClustersTemp > 0)
    for (unsigned int i = 0; i <= (unsigned int)nClustersTemp-1; i++){
      fHBAlg.Transform(allhits, &fpointId_to_clusterId, i, &nClusters, corners);
    }
  cid = nClusters;
  
  mf::LogInfo("fuzzyCluster") << "cid: " << cid ;

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















// Experimental version which actually does the merging, it is recursive!
//----------------------------------------------------------
void cluster::fuzzyClusterAlg::mergeClusters(int clusIndexStart)
{
  // If we only have one cluster, move on 
  if(fpsMembershipFinal.GetNrows() == 1)
    return;

  // If we reach the last cluster, move on 
  if(fpsMembershipFinal.GetNrows() == clusIndexStart+1)
    return;

  // Determine crip memberships
  std::vector<int>        fpsMembershipCrisp; 
  fpsMembershipCrisp.resize(fps.size());
  for (size_t pid = 0; pid < fps.size(); pid++){
    float maxClusMembership = -1;
    int iCluster = kNO_CLUSTER;
    for (int i = 0; i < fpsMembershipFinal.GetNrows(); i++){
      //mf::LogInfo("fuzzyCluster") << i  << fpsMembershipStore[nClusters-1](i,pid) ;
      if ( fpsMembershipFinal(i,pid) > maxClusMembership ) {
        maxClusMembership = fpsMembershipFinal(i,pid); 
        iCluster = i;
      }
    }
    fpsMembershipCrisp[pid] = iCluster;
  } 

  //for (size_t pid = 0; pid < fps.size(); pid++){
    //mf::LogInfo("fuzzyCluster") << pid  << fpsMembershipCrisp[pid]  << fpsMat(pid,0)  << fpsMat(pid,1) ;
  //}
  
  // Check if clusters can be merged
  int i = clusIndexStart;
  int clusToMerge = -1;
  for( int j = i + 1; j < fpsMembershipFinal.GetNrows(); j++){
    // Compute centroids for clusters i and j
    float centroidi0 = 0; 
    float centroidi1 = 0; 
    int nCentroidi = 0; 
    float centroidj0 = 0; 
    float centroidj1 = 0; 
    int nCentroidj = 0; 
    for (size_t pid = 0; pid < fps.size(); pid++){
      if ( i == fpsMembershipCrisp[pid]){
        centroidi0 += fpsMat(pid,0); 
        centroidi1 += fpsMat(pid,1); 
        nCentroidi++;
      }
      if ( j == fpsMembershipCrisp[pid]){
        centroidj0 += fpsMat(pid,0); 
        centroidj1 += fpsMat(pid,1); 
        nCentroidj++;
      }
    }
    centroidi0/=nCentroidi;
    centroidi1/=nCentroidi;
    centroidj0/=nCentroidj;
    centroidj1/=nCentroidj;
   
    double sumMinUikUjk = 0;
    double sumUik = 0;
    double sumUjk = 0;
    for( int k = 0; k < fpsMembershipFinal.GetNcols(); k++){
      sumMinUikUjk+=std::min(fpsMembershipFinal(i,k),fpsMembershipFinal(i,k));
      sumUik += fpsMembershipFinal(i,k);
      sumUjk += fpsMembershipFinal(j,k);
    }
    double Sij = sumMinUikUjk/std::min(sumUik,sumUjk);
    //std::cout << "Sij: " << Sij << std::endl;
    //if( sqrt( pow(centroidi0-centroidj0,2)+pow(centroidi1-centroidj1,2)) < fMergeCutoff){
    if( Sij < fMergeCutoff){
      mf::LogVerbatim("fuzzyCluster") << "You're going to Merge!";
      clusToMerge = j;
      break;
    }




    // Find point in cluster i closest to cluster j
    //int xiClosest  = -1;
    //float xiDistanceToCj = 9999999;
    //for(int k = 0; k < fpsMat.GetNrows(); k++){
      //if ( i != fpsMembershipCrisp[k])
        //continue;
      //if ( std::sqrt( pow(fpsMat(k,0)-centroidj0,2) + pow(fpsMat(k,1)-centroidj1 ,2)) < xiDistanceToCj){
        //xiClosest = k;
        //xiDistanceToCj = std::sqrt( pow(fpsMat(k,0)-centroidj0 ,2) + pow(fpsMat(k,1)-centroidj1 ,2));
      //}
    //}
    //// Find point in cluster j closest to cluster i
    //int xjClosest  = -1;
    //float xjDistanceToCi = 9999999;
    //for( int k = 0; k < fpsMat.GetNrows(); k++){
      //if ( j != fpsMembershipCrisp[k])
        //continue;
      //if ( std::sqrt( pow(fpsMat(k,0)-centroidi0 ,2) + pow(fpsMat(k,1)-centroidi1 ,2)) < xjDistanceToCi){
        //xjClosest = k;
        //xjDistanceToCi = std::sqrt( pow(fpsMat(k,0)-centroidi0 ,2) + pow(fpsMat(k,1)-centroidi1 ,2));
      //}
    //}
  


    
    // Find point in cluster i closest to cluster j and vice versa, we are comparing points and not centroids in 
    // this case now
    //int xiClosest  = -1;
    //int xjClosest  = -1;
    //float xiDistanceToXj = 9999999;
    //for(int k = 0; k < fpsMat.GetNrows(); k++){
      //if ( i != fpsMembershipCrisp[k])
        //continue;
      //for(int l = 0; l < fpsMat.GetNrows(); l++){
        //if ( j != fpsMembershipCrisp[l])
          //continue;
        //if ( std::sqrt( pow(fpsMat(k,0)-fpsMat(l,0),2) + pow(fpsMat(k,1)-fpsMat(l,1),2)) < xiDistanceToXj){
          //xiClosest = k;
          //xjClosest = l;
          //xiDistanceToXj = std::sqrt(pow(fpsMat(k,0)-fpsMat(l,0),2) + pow(fpsMat(k,1)-fpsMat(l,1),2));
        //}
      //} 
    //}


    // Compare distance between the two points and merge if they are close enough
    //mf::LogInfo("fuzzyCluster") << i  << j ;
    //mf::LogInfo("fuzzyCluster") << i  << centroidi0  << centroidi1 ;
    //mf::LogInfo("fuzzyCluster") << j  << centroidj0  << centroidj1 ;
    //mf::LogInfo("fuzzyCluster") << "xiClosest " << fpsMat(xiClosest,0)  << fpsMat(xiClosest,1) ;
    //mf::LogInfo("fuzzyCluster") << "xjClosest " << fpsMat(xjClosest,0)  << fpsMat(xjClosest,1) ;
    //mf::LogInfo("fuzzyCluster") << std::sqrt( pow(fpsMat(xiClosest,0)-fpsMat(xjClosest,0) ,2) + pow(fpsMat(xiClosest,1)-fpsMat(xjClosest,1) ,2)) ;
    //if (xiClosest > -1 
        //&& xjClosest > -1 
        //&& std::sqrt( pow(fpsMat(xiClosest,0)-fpsMat(xjClosest,0) ,2) + pow(fpsMat(xiClosest,1)-fpsMat(xjClosest,1) ,2)) < fMergeCutoff){
      //// Perform the merging and deletion of the closest cluster
      //mf::LogVerbatim("fuzzyCluster") << "You're going to Merge!"  ;
      //clusToMerge = j;
      //break;
    //}
   


  } 

  if(clusToMerge > 0){
    TMatrixT<double> fpsMembershipTemp(fpsMembershipFinal.GetNrows()-1, fpsMat.GetNrows());
    TMatrixDRow(fpsMembershipFinal,clusIndexStart) += TMatrixDRow(fpsMembershipFinal,clusToMerge);
    for( int j = 0; j < fpsMembershipFinal.GetNrows()-1; j++){
      if(j < clusToMerge)  TMatrixDRow(fpsMembershipTemp,j) = TMatrixDRow(fpsMembershipFinal,j);
      if(j >= clusToMerge) TMatrixDRow(fpsMembershipTemp,j) = TMatrixDRow(fpsMembershipFinal,j+1);
    }
    fpsMembershipFinal.ResizeTo(fpsMembershipTemp);
    fpsMembershipFinal=fpsMembershipTemp;
    
    mergeClusters(clusIndexStart);
  }
  else
    mergeClusters(clusIndexStart+1);

  return;

}







////----------------------------------------------------------
//void cluster::fuzzyClusterAlg::mergeClusters(int clusIndexStart)
//{
  //// Determine crips memberships
  //std::vector<int>        fpsMembershipCrisp; 
  //fpsMembershipCrisp.resize(fps.size());
  //for (size_t pid = 0; pid < fps.size(); pid++){
    //float maxClusMembership = -1;
    //int iCluster = kNO_CLUSTER;
    //for (int i = 0; i < fpsMembershipFinal.GetNrows(); i++){
      ////mf::LogInfo("fuzzyCluster") << i  << fpsMembershipStore[nClusters-1](i,pid) ;
      //if ( fpsMembershipFinal(i,pid) > maxClusMembership ) {
        //maxClusMembership = fpsMembershipFinal(i,pid); 
        //iCluster = i;
      //}
    //}
    //fpsMembershipCrisp[pid] = iCluster;
  //} 

  //for (size_t pid = 0; pid < fps.size(); pid++){
    //mf::LogInfo("fuzzyCluster") << pid  << fpsMembershipCrisp[pid]  << fpsMat(pid,0)  << fpsMat(pid,1) ;
  //}
  
  //// Check if clusters can be merged
  //for( int i = clusIndexStart; i < fpsMembershipFinal.GetNrows(); i++){
    //for( int j = i + 1; j < fpsMembershipFinal.GetNrows(); j++){
      //// Compute centroids for clusters i and j
      //float centroidi0 = 0; 
      //float centroidi1 = 0; 
      //int nCentroidi = 0; 
      //float centroidj0 = 0; 
      //float centroidj1 = 0; 
      //int nCentroidj = 0; 
      //for (size_t pid = 0; pid < fps.size(); pid++){
        //if ( i == fpsMembershipCrisp[pid]){
          //centroidi0 += fpsMat(pid,0); 
          //centroidi1 += fpsMat(pid,1); 
          //nCentroidi++;
        //}
        //if ( j == fpsMembershipCrisp[pid]){
          //centroidj0 += fpsMat(pid,0); 
          //centroidj1 += fpsMat(pid,1); 
          //nCentroidj++;
        //}
      //}
      //centroidi0/=nCentroidi;
      //centroidi1/=nCentroidi;
      //centroidj0/=nCentroidj;
      //centroidj1/=nCentroidj;
      //// Find point in cluster i closest to cluster j
      //int xiClosest  = -1;
      //float xiDistanceToCj = 9999999;
      //for(int k = 0; k < fpsMat.GetNrows(); k++){
        //if ( i != fpsMembershipCrisp[k])
          //continue;
        //if ( std::sqrt( pow(fpsMat(k,0)-centroidj0,2) + pow(fpsMat(k,1)-centroidj1 ,2)) < xiDistanceToCj){
          //xiClosest = k;
          //xiDistanceToCj = std::sqrt( pow(fpsMat(k,0)-centroidj0 ,2) + pow(fpsMat(k,1)-centroidj1 ,2));
        //}
      //}
      //// Find point in cluster j closest to cluster i
      //int xjClosest  = -1;
      //float xjDistanceToCi = 9999999;
      //for( int k = 0; k < fpsMat.GetNrows(); k++){
        //if ( j != fpsMembershipCrisp[k])
          //continue;
        //if ( std::sqrt( pow(fpsMat(k,0)-centroidi0 ,2) + pow(fpsMat(k,1)-centroidi1 ,2)) < xjDistanceToCi){
          //xjClosest = k;
          //xjDistanceToCi = std::sqrt( pow(fpsMat(k,0)-centroidi0 ,2) + pow(fpsMat(k,1)-centroidi1 ,2));
        //}
      //}
      
      //// Compare distance between the two points and merge if they are close enough
      //float fMergeCutoff = 1;
      //mf::LogInfo("fuzzyCluster") << i  << j ;
      //mf::LogInfo("fuzzyCluster") << i  << centroidi0  << centroidi1 ;
      //mf::LogInfo("fuzzyCluster") << j  << centroidj0  << centroidj1 ;
      //mf::LogInfo("fuzzyCluster") << "xiClosest " << fpsMat(xiClosest,0)  << fpsMat(xiClosest,1) ;
      //mf::LogInfo("fuzzyCluster") << "xjClosest " << fpsMat(xjClosest,0)  << fpsMat(xjClosest,1) ;
      //mf::LogInfo("fuzzyCluster") << std::sqrt( pow(fpsMat(xiClosest,0)-fpsMat(xjClosest,0) ,2) + pow(fpsMat(xiClosest,1)-fpsMat(xjClosest,1) ,2)) ;
      //if (std::sqrt( pow(fpsMat(xiClosest,0)-fpsMat(xjClosest,0) ,2) + pow(fpsMat(xiClosest,1)-fpsMat(xjClosest,1) ,2)) < fMergeCutoff){

        //// Perform the merging and deletion of the closest cluster
        //mf::LogInfo("fuzzyCluster") << "Yo Dawg! You Can Merge!"  ;
        
      
      //}


    //} 
  //}

  //return;

//}




