#ifndef GAUSHITFINDER_H
#define GAUSHITFINDER_H

#include <string>
#include "art/Framework/Core/EDProducer.h" 
#include "TTree.h"
#include "TH1D.h"
///localizations of energy depositions
namespace hit {
   
  class GausHitFinder : public art::EDProducer {
    
  public:
    
    explicit GausHitFinder(fhicl::ParameterSet const& pset); 
    virtual ~GausHitFinder();
         
    void produce(art::Event& evt); 
    void beginJob(); 
    void endJob(); 
    void reconfigure(fhicl::ParameterSet const& p);                

  private:
    std::string     fCalDataModuleLabel;
    double          fMinSigInd;     ///<Induction signal height threshold 
    double          fMinSigCol;     ///<Collection signal height threshold 
    double          fIndWidth;      ///<Initial width for induction fit
    double          fColWidth;      ///<Initial width for collection fit
    double          fIndMinWidth;   ///<Minimum induction hit width
    double          fColMinWidth;   ///<Minimum collection hit width
    int             fMaxMultiHit;   ///<maximum hits for multi fit 
    int             fAreaMethod;    ///<Type of area calculation  
    std::vector<double> fAreaNorms; ///<factors for converting area to same units as peak height 
    
    	double	WireNumber[100000];
	double 	TotalSignal[100000];
	double 	StartIime;
	double 	StartTimeError;
	double 	EndTime;
	double 	EndTimeError;
	int	NumOfHits;
	double	MeanPosition;
	double 	MeanPosError;
	double	Amp;
	double	AmpError;
	double	Charge;
	double	ChargeError;
	double	FitGoodnes;
		
  protected: 
    
  
  }; // class GausHitFinder


}

#endif // GausHITFINDER_H
