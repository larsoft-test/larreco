////////////////////////////////////////////////////////////////////////
//
// ClusterFinder base class:
// Supply the basic methods for cluster finding algorithms.
// Users must define cluster finding specifics in inheriting classes.
//
// kinga.partyka@yale.edu
//
// joshua.spitz@yale.edu
////////////////////////////////////////////////////////////////////////
#ifndef CLUSTERFINDER_H
#define CLUSTERFINDER_H

#include "FWCore/Framework/interface/Event.h" 
#include "FWCore/Framework/interface/EDProducer.h"
#include "RecoBase/Wire.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RawData/DAQHeader.h"

#include <vector>
#include <string>


class TH1F;
class TH2F;
///Cluster finding and building 
namespace cluster {

   
  class ClusterFinder : public edm::EDProducer {

  public:
          
    explicit ClusterFinder(edm::ParameterSet const& pset); 
    virtual ~ClusterFinder();
 
    /// write access to event
    void produce(edm::Event& evt, edm::EventSetup const&);
    // analyze() Must go into new module. EC, 5-Oct-2010.
    //    void analyze(const edm::Event& evt);
    void beginJob(edm::EventSetup const&);

  private:
    std::string fHitModuleLabel;
    std::string fClusterFinderModuleLabel;
    
    void Make(edm::Event&);  ///< create new Hits
    void Read(edm::Event&);  ///< read in existing Hits
    void GetData(edm::Event&);   ///< get input data as needed

    std::vector<const recob::Cluster*>  fCluster; ///< const buffer of input clusters

      /* This seems dicey to make fClusterVecCol a member of the class. 
	 But it's got to be for ClusterFinder to be useful as a base class.
	 EC, 6-Oct-2010.
       */


    ///<parameters to set
    int          fMakeClusters;         ///< flag to make hits or read them from a file
    virtual void FindCluster(edm::PtrVector<const recob::Hit> &hit)=0;
    
      	 
  }; // class ClusterFinder

}

#endif // CLUSTERFINDER_H
