////////////////////////////////////////////////////////////////////////
//
// HoughLineFinder class
//
// joshua.spitz@yale.edu
//
//  This algorithm is designed to find lines (Houghclusters) from clusters found by DBSCAN after deconvolution and hit finding.
//  The algorithm is based on: 
//  Queisser, A. "Computing the Hough Transform", C/C++ Users Journal 21, 12 (Dec. 2003).
//  Niblack, W. and Petkovic, D. On Improving the Accuracy of the Hough Transform", Machine Vision and Applications 3, 87 (1990)  
////////////////////////////////////////////////////////////////////////

#include "HoughLineFinder.h"
extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}
// ROOT includes
#include <TCanvas.h>
#include "TDatabasePDG.h"
#include "TSystem.h"

#include <sstream>
#include <fstream>
#include <math.h>
#include <algorithm>


#include "FWCore/Framework/interface/Event.h" 
#include "FWCore/ParameterSet/interface/ParameterSet.h" 
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h" 
#include "DataFormats/Common/interface/Handle.h" 
#include "DataFormats/Common/interface/Ptr.h" 
#include "DataFormats/Common/interface/PtrVector.h" 
#include "FWCore/Framework/interface/MakerMacros.h" 
#include "FWCore/ServiceRegistry/interface/Service.h" 
#include "FWCore/Services/interface/TFileService.h" 
#include "FWCore/Framework/interface/TFileDirectory.h" 
#include "FWCore/MessageLogger/interface/MessageLogger.h" 
 
#include "RawData/RawDigit.h"
#include "Filters/ChannelFilter.h"
#include "SimulationBase/simbase.h"
#include "RecoBase/recobase.h"
#include "Geometry/geo.h"




cluster::HoughLineFinder::HoughLineFinder(edm::ParameterSet const& pset) : 
  fDBScanModuleLabel       (pset.getParameter< std::string >("DBScanModuleLabel")),
  fMaxLines                (pset.getParameter< int >("MaxLines")),
  fMinHits                 (pset.getParameter< int >("MinHits")),
  fSaveAccumulator         (pset.getParameter< int >("SaveAccumulator")),
  fNumAngleCells           (pset.getParameter< int >("NumAngleCells")),
  fRhoResolutionFactor     (pset.getParameter< int >("RhoResolutionFactor")),
  fSmootherSigma           (pset.getParameter< double >("SmootherSigma")),
  fMaxDistance             (pset.getParameter< double >("MaxDistance")),
  fRhoZeroOutRange         (pset.getParameter< int >("RhoZeroOutRange")),
  fThetaZeroOutRange       (pset.getParameter< int >("ThetaZeroOutRange")),
  fPerCluster              (pset.getParameter< int >("HitsPerCluster"))
{
  produces< std::vector<recob::Cluster> >();
}

cluster::HoughLineFinder::~HoughLineFinder()
{
}

cluster::HoughLineFinder::HoughTransform::HoughTransform()
{  
}

cluster::HoughLineFinder::HoughTransform::~HoughTransform()
{  
}

void cluster::HoughLineFinder::produce(edm::Event& evt, edm::EventSetup const&)
{

  //////////////////////////////////////////////////////
  // here is how to get a collection of objects out of the file
  // and connect it to a edm::Handle
  //////////////////////////////////////////////////////
  // Read in the clusterList object(s).
  edm::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fDBScanModuleLabel,clusterListHandle);
  //Point to a collection of clusters to output.
  std::auto_ptr<std::vector<recob::Cluster> > ccol(new std::vector<recob::Cluster>);

  std::vector<int> skip;  

  edm::Service<geo::Geometry> geom;
  filter::ChannelFilter chanFilt;
  HoughTransform c;
  HoughTransform cc;
  HoughTransform ccc;
  extern void SaveBMPFile(const char *f, unsigned char *pix, int dxx, int dyy);
  edm::PtrVector<recob::Hit> cHits;
  edm::PtrVector<recob::Hit> hit;

  edm::PtrVector<recob::Cluster> clusIn;
  std::cout<<" clusterListHandle->size() "<<clusterListHandle->size()<<std::endl;
  for(unsigned int ii = 0; ii < clusterListHandle->size(); ++ii)
    {
      edm::Ptr<recob::Cluster> cluster(clusterListHandle, ii);
      clusIn.push_back(cluster);
    }

  for(int p = 0; p < geom->Nplanes(); p++) {
     std::cout<<"clusIn.size() "<<clusIn.size()<<std::endl;
    edm::PtrVectorItr<recob::Cluster> clusterIter = clusIn.begin();

    // This is the loop over clusters. The algorithm searches for lines on a (DBSCAN) cluster-by-cluster basis. 
    while(clusterIter!=clusIn.end()) 
      {
 	hit.clear();
 	cHits.clear();
 	if(fPerCluster)
 	{
 	  hit = (*clusterIter)->Hits(p,-1);
 	  std::cout<<"hit size: "<<hit.size()<<std::endl;
 	  }
 	else 
 	  {   
 	    while(clusterIter!=clusIn.end()) 
 	      {
 		cHits = (*clusterIter)->Hits(p,-1);
 		if(cHits.size() > 0)
		  {
		    // hit.insert(hit.end(),cHits.begin(),cHits.end());
		    edm::PtrVectorItr<recob::Hit> hitIter = cHits.begin();
		    while (hitIter!=cHits.end())
		      {
			hit.push_back((*hitIter));
			hitIter++;
		      }
		    clusterIter++;
		  }
 	      } 
 	  }
 	std::cout<<"hit.size() "<<hit.size()<<std::endl;  
 	if(hit.size() == 0) 
 	  { 
 	    if(fPerCluster) clusterIter++;
 	    continue;
 	  }

 	int x, y;
	unsigned int channel, plane, wire;
 	//there must be a better way to find which plane a cluster comes from
 	int dx=geom->Nwires(0);//number of wires 
 	int dy=hit[0]->Wire()->fSignal.size();//number of time samples. 
 	skip.clear();
 	skip.resize(hit.size());
 	int clusterID=0;//the unique ID of the cluster
 	std::vector<int> listofxmax;
 	std::vector<int> listofymax;  
 	std::vector<int> hitTemp;  //indecies ofcandidate hits
 	std::vector<int> sequenceHolder;  //channels of hits in list
 	std::vector<int> currentHits; //working vector of hits 
 	std::vector<int> lastHits;  //best list of hits
 	edm::PtrVector<recob::Hit> clusterHits;
 	double indcolscaling=0.;//a parameter to account for the different characteristic hit width of induction and collection plane
 	double centerofmassx=0;
 	double centerofmassy=0;
 	double denom = 0; 
 	double intercept=0.;
 	double slope=0.;
 	//this array keeps track of the hits that have already been associated with a line. 
 	int xMax=0;
 	int yMax=0;
 	double rho;
 	double theta;
 	double norm=100.;//normalize. This is important since newcellvalue will end up as an int.  
 	int accDx(0), accDy(0);
	
 	for (int linenum = 0; linenum < fMaxLines; linenum++)
 	  { 
 	  
 	  std::cout<<"here"<<std::endl;
 	    //Init specifies the size of the two-dimensional accumulator (based on the arguments, number of wires and number of time samples). 
 	    c.Init(dx,dy,fRhoResolutionFactor,fNumAngleCells);
 	    //initialize the smoothing accumulators as well, one each for the two dimensions of the accumulator
 	    cc.Init(dx,dy,fRhoResolutionFactor,fNumAngleCells);
 	    ccc.Init(dx,dy,fRhoResolutionFactor,fNumAngleCells); 
 	    //adds all of the hits (that have not yet been associated with a line) to the accumulator
 	    for(unsigned int i=0;i < hit.size(); i++)
 	      {
 		channel=hit[i]->Wire()->RawDigit()->Channel();
 		geom->ChannelToWire(channel,plane,wire);
 		if (skip[i] != 1)
 		  c.AddPoint(wire,(int)(hit[i]->CrossingTime()));
 	      }
 	    //gets the actual two-dimensional size of the accumulator
 	    c.GetAccumSize(accDy, accDx);
	    
 	    if(fSmootherSigma>0)  
 	      {  
		
		double Weight[accDx];
		double newcellvalue[accDx];
		for(int i=0;i<=(3*fSmootherSigma);i++)
		  {
		    //find the 2D Gaussian's weights
		    Weight[i]=norm*exp((-(pow(i,2)))/(2.*pow(fSmootherSigma,2))); 
		  }
		//two separate 1D Gaussian blurs are used to form a 2D Gaussian blur.
		//smoothing in x (rho) direction 
		for (y=0; y<accDy; y++)
		  {
		    for (x=0; x<accDx; x++)
		      {
			for(int i=0;i<=(3*fSmootherSigma);i++)
			  {
			    //if the Gaussian falls off the accumulator edge, the accumulator is extended using mirrored cells 
			    if(i==0)
			      newcellvalue[x] = Weight[i]*c.GetCell(y,x);
			    else if(x+i>=accDx && x-i > 0)
			      newcellvalue[x] += 2.*Weight[i]*c.GetCell(y,x-i);
			    else if(x-i<0 && x+i < accDx)
			      newcellvalue[x] += 2.*Weight[i]*c.GetCell(y,x+i);
			    else
			      newcellvalue[x] += Weight[i]*(c.GetCell(y,x+i)+c.GetCell(y,x-i));                             
			  }   
			cc.SetCell(y,x,(int)newcellvalue[x]);
		      }
		  }
		//smoothing in y (theta) direction
		for (y=0; y<accDy; y++)
		  {
		    for (x=0; x<accDx; x++)
		      {
			for(int i=0;i<=(3*fSmootherSigma);i++)
			  {
			    //if the Gaussian falls off the accumulator edge, the accumulator is extended using mirrored cells 
			    if(i==0)
			      newcellvalue[x] = Weight[i]*cc.GetCell(y,x);
			    else if(y+i>=accDy && y-i > 0)
			      newcellvalue[x] += 2.*Weight[i]*cc.GetCell(y-i,x);
			    else if(y-i<0 && y+i < accDy)
			      newcellvalue[x] += 2.*Weight[i]*cc.GetCell(y+i,x);
			    else
			      newcellvalue[x] += Weight[i]*(cc.GetCell(y+i,x)+cc.GetCell(y-i,x));                    
			  } 
			ccc.SetCell(y,x,(int)newcellvalue[x]);
		      }
		  }

 	      } 
 	    else
	      {
		ccc = c;
	      }
 	    // zeroes out the neighborhood of all previous lines  
 	    for(unsigned int i=0;i<listofxmax.size();i++)
 	      {
 		int yClearStart = listofymax[i]-fRhoZeroOutRange;
 		if (yClearStart < 0)
 		  yClearStart = 0;
 		int yClearEnd = listofymax[i]+fRhoZeroOutRange;
 		if (yClearEnd >= accDy)
 		  yClearEnd = accDy-1;
 		int xClearStart = listofxmax[i]-fThetaZeroOutRange;
 		if (xClearStart < 0)
 		  xClearStart = 0;
 		int xClearEnd = listofxmax[i]+fThetaZeroOutRange;
 		if (xClearEnd >= accDx)
 		  xClearEnd = accDx-1;

 		for (y=yClearStart; y<=yClearEnd; y++)
 		  {
 		    for (x=xClearStart; x<=xClearEnd; x++)
 		      {
 			ccc.SetCell(y,x,0);
 		      }
 		  }
 	      }
 	    //find the weightiest cell in the smoothed accumulator.
 	    std::cout<<"here2"<<std::endl;
 	    int maxCell = 0;
 	    xMax = 0;
 	    yMax = 0;
 	    maxCell = ccc.GetMax(yMax,xMax);
 	    // break when biggest maximum is smaller than fMinHits
 	    if ((maxCell < fMinHits*pow(norm,2) && fSmootherSigma>0) || (maxCell < fMinHits && fSmootherSigma==0.)) 
 	      break;
	  
 	    //find the center of mass of the 3x3 cell system (with maxCell at the center). 
 	    denom=centerofmassx=centerofmassy=0;

 	    if(xMax>0 && yMax>0 && xMax+1<accDx && yMax+1<accDy)
 	      {  
 		for(int i = -1; i < 2; i++) 
 		  {
 		    for(int j = -1; j < 2; j++) 
 		      {
 			denom+=ccc.GetCell(yMax+i,xMax+j);
 			centerofmassx+=j*ccc.GetCell(yMax+i,xMax+j);
 			centerofmassy+=i*ccc.GetCell(yMax+i,xMax+j);
 		      }
 		  }
 		centerofmassx/=denom;
 		centerofmassy/=denom;      
 	      }
 	    else
 	      {
 		centerofmassx=centerofmassy=0;
 	      }
 	    //fill the list of cells that have already been found
 	    listofxmax.push_back(xMax);
 	    listofymax.push_back(yMax);
	    
 	    ccc.GetEquation(yMax+centerofmassy, xMax+centerofmassx, rho, theta);
 	    slope=-1./tan(theta);    
 	    intercept=(rho/sin(theta));
	    std::cout<<"here3"<<std::endl;
 	    double distance;
 	    if(p==0)
 	      indcolscaling=5.;
 	    else
 	      indcolscaling=0.;
 	    //the collection plane's characteristic hit width's are, on average, about 5 time samples wider than the induction plane's. this is hard-coded for now.
 	    if(!isinf(slope) && !isnan(slope)){
 	      sequenceHolder.clear();
 	      hitTemp.clear();
 	      for(unsigned int i=0;i<hit.size();i++)
 		{
 		  channel=hit[i]->Wire()->RawDigit()->Channel();
 		  geom->ChannelToWire(channel,plane,wire);
 		  //Note that there are 1/0.0743=13.46 time samples per 4.0 mm (wire pitch in ArgoNeuT), assuming a 1.5 mm/us drift velocity for a 500 V/cm E-field 
 		  distance=(TMath::Abs(hit[i]->CrossingTime()-slope*(double)(wire)-intercept)/(sqrt(pow(.0743*slope,2)+1)));
	
 		  if(distance < fMaxDistance+((hit[i]->EndTime()-hit[i]->StartTime())/2.)+indcolscaling && fPerCluster==1){
 		    hitTemp.push_back(i);
 		    sequenceHolder.push_back(channel);
 		  }
 		  else if(distance < fMaxDistance+((hit[i]->EndTime()-hit[i]->StartTime())/2.)+indcolscaling && fPerCluster==0 && skip[i]!=1){
 		    hitTemp.push_back(i);
 		    sequenceHolder.push_back(channel);
 		  }
 		}
 		std::cout<<"here4"<<std::endl;
 	      if(hitTemp.size() < 2) continue;
 	      currentHits.clear();  
 	      lastHits.clear();
 	      int j; 
 	      currentHits.push_back(0);
 	      for(unsigned int i=0;i+1<sequenceHolder.size();i++) {  
 		j = 1;
 		while((chanFilt.BadChannel(sequenceHolder[i]+j))==true) j++;
 		if(sequenceHolder[i+1]-sequenceHolder[i]<=j+1) currentHits.push_back(i+1);
 		else if(currentHits.size() > lastHits.size()) {
 		  lastHits = currentHits;
 		  currentHits.clear();
 		}
 		else currentHits.clear();
 	      } 
 	      if(currentHits.size() > lastHits.size()) lastHits = currentHits;
 	      clusterHits.clear();    
 	      for(unsigned int i = 0; i < lastHits.size();i++) {
 		clusterHits.push_back(hit[hitTemp[lastHits[i]]]);
 		skip[hitTemp[lastHits[i]]]=1;
 	      } 
 	      //protection against very steep uncorrelated hits
 	      if(TMath::Abs(slope)>75. 
		 && TMath::Abs((*clusterHits.begin())->Wire()->RawDigit()->Channel()-
			       (*clusterHits.end())->Wire()->RawDigit()->Channel())>0
		 )
 		continue;
	      
        std::cout<<"here5"<<std::endl;
 	      recob::Cluster cluster(clusterHits);	      

 	      cluster.SetSlope(slope);
 	      cluster.SetIntercept(intercept);
	      channel = (*clusterHits.begin())->Wire()->RawDigit()->Channel(); 
 	      geom->ChannelToWire(channel,plane,wire);
 	      cluster.SetStartWire(wire);
 	      cluster.SetStartTime((*clusterHits.begin())->CrossingTime());
 	      //cluster->SetStartTime(slope*(double)(wire)+intercept);
 	      channel = (*clusterHits.end())->Wire()->RawDigit()->Channel(); 
 	      geom->ChannelToWire(channel,plane,wire);
 	      cluster.SetEndWire(wire);        
 	      cluster.SetEndTime((*clusterHits.end())->CrossingTime());
 	      //cluster->SetEndTime(slope*(double)(wire)+intercept);
 	      cluster.SetID(clusterID);

 	      clusterID++;
	     
 	      ccol->push_back(cluster);
	     
	    }

	  } 
  
	//saves a bitmap image of the accumulator (useful for debugging), with scaling based on the maximum cell value
	if(fSaveAccumulator)
	  {   
	    unsigned char *outPix = new unsigned char [accDx*accDy];
	    //finds the maximum cell in the accumulator for image scaling
	    int cell, pix=0, maxCell=0;
	    int xmaxx, ymaxx;
	    for (y=0; y<accDy; y++)
	      for (x=0; x<accDx; x++)
		{
		  cell = ccc.GetCell(y,x);
		  if (cell > maxCell){
		    maxCell = cell;
		    xmaxx=x;
		    ymaxx=y;}
		}
	    for (y=0; y<accDy; y++)
	      for (x=0; x<accDx; x++)
		{ 
		  //scales the pixel weights based on the maximum cell value     
		  if(maxCell>0)
		    pix = (int)((1500*ccc.GetCell(y,x))/maxCell);
		  outPix[y*accDx + x] = pix;
		}
	    SaveBMPFile("houghaccum.bmp", outPix, accDx, accDy);
	    delete [] outPix;
	  }
	
	hit.clear();
	if(clusterIter!=clusIn.end()) clusterIter++;

      }
  }

  evt.put(ccol);
 
    

 
}

bool cluster::HoughLineFinder::HoughTransform::AddPoint(int x, int y)
{
  if (x>m_dx || y>m_dy || x<0.0 || y<0.0)
    return false;
  return DoAddPoint(x, y);
}
 
void cluster::HoughLineFinder::HoughTransform::Init(int dx, int dy, int rhores,
						    int numACells)
{
  m_numAngleCells=numACells;
  m_rhoResolutionFactor = rhores;
  m_accum.clear();
  m_accum.resize(m_numAngleCells);
  m_numAccumulated = 0;   
  m_cosTable.resize(m_numAngleCells);
  m_sinTable.resize(m_numAngleCells);
  if (dx == m_dx && dy == m_dy)
    return;
  m_dx = dx;
  m_dy = dy;
  m_rowLength = (int)(m_rhoResolutionFactor*2. * sqrt(dx*dx + dy*dy));
  
  unsigned int angleIndex;
  double a, angleStep = TMath::Pi()/m_numAngleCells;
  for (a=0.0, angleIndex=0; angleIndex<m_cosTable.size(); angleIndex++)
    {
      m_cosTable[angleIndex] = cos(a);
      m_sinTable[angleIndex] = sin(a);
      a += angleStep;
    }
}

int cluster::HoughLineFinder::HoughTransform::GetMax(int &xmax, int &ymax)
{
  std::map<int,int>::iterator rhoIter;
  int maxVal = -1;
  for(unsigned int i = 0; i < m_accum.size(); i++)
    {
      for(rhoIter=m_accum[i].begin(); rhoIter!=m_accum[i].end(); rhoIter++)
	{
          if((*rhoIter).second > maxVal) {
            maxVal = (*rhoIter).second;
            xmax = i;
            ymax = (*rhoIter).first;
          }
	}
    }
  return maxVal;
}

bool cluster::HoughLineFinder::HoughTransform::DoAddPoint(int x, int y)
{
  int distCenter = (int)(m_rowLength/2.);
 
  // prime the lastDist variable so our linear fill works below
  int lastDist = (int)(distCenter+(m_rhoResolutionFactor*(m_cosTable[0]*x + m_sinTable[0]*y))) ;
  // loop through all angles a from 0 to 180 degrees
  for (unsigned int a=1; a<m_cosTable.size(); a++)
    {
      // Calculate the basic line equation dist = cos(a)*x + sin(a)*y.
      // Shift to center of row to cover negative values
      int dist = (int)(distCenter+(m_rhoResolutionFactor*(m_cosTable[a]*x + m_sinTable[a]*y))) ;
      // sanity check to make sure we stay within our row
      if (dist >= 0 && dist<m_rowLength)
	{
	  if(lastDist==dist)
	    m_accum[a][lastDist]++;
	  else
	    {
	      // fill in all values in row a, not just a single cell
	      int stepDir = dist>lastDist ? 1 : -1;
	      int cell;
	      for (cell=lastDist; cell!=dist; cell+=stepDir)
		{   
		  m_accum[a][cell]++;//maybe add weight of hit here?
		}      
	    }
	}      
      lastDist = dist;
    }
  m_numAccumulated++;
  return true;
}

//this method saves a BMP image of the Hough Accumulator, which can be viewed with gimp
void SaveBMPFile(const char *fileName, unsigned char *pix, int dx, int dy)
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



