#ifndef HOUGHLINEFINDER_H
#define HOUGHLINEFINDER_H


#include "TMath.h"
#include "TObject.h"
#include <vector>
#include <string>

#include "FWCore/Framework/interface/EDProducer.h"

class TH1F;
class TTree;
namespace cluster {
   
  class HoughLineFinder : public edm::EDProducer {
    
  public:
    
    explicit HoughLineFinder(edm::ParameterSet const& pset); 
    ~HoughLineFinder();
         
    void produce(edm::Event& evt, edm::EventSetup const&);
    void beginJob(edm::EventSetup const&);

     
    class HoughTransform {
    public:

      HoughTransform();
      ~HoughTransform();

  
      void Init(int dx, int dy, int rhoresfact, int numACells);
      bool AddPoint(int x, int y);
      int GetCell(int row, int col) { return m_accum[row][col]; }
      void SetCell(int row, int col, int value) { m_accum[row][col] = value; }
      void IncrementCell(int row, int col)      { m_accum[row][col]++;}
      void GetAccumSize(int &numRows, int &numCols) 
      { 
	numRows = m_accum.size();
	numCols  = m_rowLength;
      }
      int NumAccumulated() { return m_numAccumulated; }
      void GetEquation(double row, double col, double &rho, double &theta)
      {
	theta = (TMath::Pi()*row)/m_numAngleCells;
	rho   = (col - (m_rowLength/2.))/m_rhoResolutionFactor;
      }
      int GetMax(int & xmax, int & ymax);

      
    private:
   
      int m_dx, m_dy;
      std::vector<std::map<int,int> > m_accum;  // column=rho, row=theta
      int m_rowLength;
      int m_numAccumulated;
      int m_rhoResolutionFactor;
      int m_numAngleCells;
      std::vector<double> m_cosTable;
      std::vector<double> m_sinTable;
      bool DoAddPoint(int x, int y);
    }; // class HoughTransform  
     
  private:
    std::string fDBScanModuleLabel;    
    int    fMaxLines;      //Max number of lines that can be found 
    int    fMinHits;       //Min number of hits in the accumulator to consider (number of hits required to be considered a line).
    int    fSaveAccumulator;  //Save bitmap image of accumulator for debugging?
    int    fNumAngleCells;    //Number of angle cells in the accumulator (a measure of the angular resolution of the line finder). If this number is too large than the number of votes that fall into the "correct" bin will be small and consistent with noise.
    double fSmootherSigma;
    double fMaxDistance;
    int fRhoZeroOutRange;
    int fThetaZeroOutRange;
    int fRhoResolutionFactor;
    int fPerCluster;
      
  protected:

    friend class HoughTransform;
  };
  
  
}



#endif // HOUGHLineFINDER_H
