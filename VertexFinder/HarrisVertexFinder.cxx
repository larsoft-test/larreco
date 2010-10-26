////////////////////////////////////////////////////////////////////////
//
// HarrisVertexFinder class
//
// joshua.spitz@yale.edu
//
//  This algorithm is designed to find (weak) vertices from hits after deconvolution and hit finding. 
//  A weak vertex is a vertex that has been found using a dedicated vertex finding algorithm only. A 
//  strong vertex is a vertex that has been found using a dedicated vertex finding algorithm and matched 
//  to a crossing of two or more HoughLineFinder lines. The VertexMatch module finds strong vertices.
////////////////////////////////////////////////////////////////////////
/// The algorithm is based on:
///C. Harris and M. Stephens (1988). "A combined corner and edge detector". Proceedings of the 4th Alvey 
///Vision Conference. pp. 147-151.
///B. Morgan (2010). "Interest Point Detection for Reconstruction in High Granularity Tracking Detectors". 
///arXiv:1006.3012v1 [physics.ins-det]
//Thanks to B. Morgan of U. of Warwick for comments and suggestions

#include <iostream>

// Framework includes
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

#include "HarrisVertexFinder.h"
extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

#include <sstream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include "TMath.h"


#include "RawData/RawDigit.h"
#include "Filters/ChannelFilter.h"
#include "SimulationBase/simbase.h"
#include "RecoBase/recobase.h"
#include "Geometry/geo.h"

//-----------------------------------------------------------------------------
vertex::HarrisVertexFinder::HarrisVertexFinder(edm::ParameterSet const& pset) :
  fHitsModuleLabel  (pset.getParameter< std::string >("HitsModuleLabel")),
  fTimeBins        (pset.getParameter< int         >("TimeBins")      ),
  fMaxCorners      (pset.getParameter< int         >("MaxCorners")    ),
  fGsigma          (pset.getParameter< double      >("Gsigma")        ),
  fWindow          (pset.getParameter< int         >("Window")        ),
  fThreshold       (pset.getParameter< double      >("Threshold")     ),
  fSaveVertexMap   (pset.getParameter< int         >("SaveVertexMap") )
{
  produces< std::vector<recob::Vertex> >();
}

//-----------------------------------------------------------------------------
vertex::HarrisVertexFinder::~HarrisVertexFinder()
{
}

//-----------------------------------------------------------------------------
double vertex::HarrisVertexFinder::Gaussian(int x, int y, double sigma)
{
  double Norm=1./sqrt(2*TMath::Pi()*pow(sigma,2));
  double value=Norm*exp(-(pow(x,2)+pow(y,2))/(2*pow(sigma,2)));
  return value;
}

//-----------------------------------------------------------------------------
double vertex::HarrisVertexFinder::GaussianDerivativeX(int x,int y)
{
  double Norm=1./(sqrt(2*TMath::Pi())*pow(fGsigma,3));
  double value=Norm*(-x)*exp(-(pow(x,2)+pow(y,2))/(2*pow(fGsigma,2)));
  return value;
}

//-----------------------------------------------------------------------------
double vertex::HarrisVertexFinder::GaussianDerivativeY(int x,int y)
{
  double Norm=1./(sqrt(2*TMath::Pi())*pow(fGsigma,3));
  double value=Norm*(-y)*exp(-(pow(x,2)+pow(y,2))/(2*pow(fGsigma,2)));
  return value;
}

//-----------------------------------------------------------------------------
void vertex::HarrisVertexFinder::produce(edm::Event& evt, edm::EventSetup const&)
{

  edm::Service<geo::Geometry> geom;
  
  extern void SaveBMPFile(const char *f, unsigned char *pix, int dxx, int dyy);
  
  edm::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fHitsModuleLabel,clusterListHandle);
  //Point to a collection of vertices to output.
  std::auto_ptr<std::vector<recob::Vertex> > vtxcol(new std::vector<recob::Vertex>);

  filter::ChannelFilter chanFilt;  
  edm::PtrVector<recob::Hit> cHits;
  edm::PtrVector<recob::Hit> hit;
   
  edm::PtrVector<recob::Cluster> clusIn;
  for(unsigned int ii = 0; ii < clusterListHandle->size(); ++ii)
    {
      edm::Ptr<recob::Cluster> cluster(clusterListHandle, ii);
      clusIn.push_back(cluster);
    }

  int flag=0;
  int windex=0;//the wire index to make sure the vertex finder does not fall off the edge of the hit map
  int tindex=0;//the time index to make sure the vertex finder does not fall off the edge of the hit map
  int n=0; //index of window cell. There are 49 cells in the 7X7 Gaussian and Gaussian derivative windows
  int numberwires;
  double numbertimesamples;
  double MatrixAAsum,MatrixBBsum,MatrixCCsum;
  std::vector<double> Cornerness2;
  //gaussian window definitions. The cell weights are calculated here to help the algorithm's speed
  double  w[49]={0.};
  double wx[49]={0.};
  double wy[49]={0.};
  int ctr = 0;
  for(int i = -3; i < 4; ++i){
    for(int j = 3; j > -4; --j){
      w[ctr] = Gaussian(i, j, fGsigma);
      wx[ctr] = GaussianDerivativeX(i,j);
      wy[ctr] = GaussianDerivativeY(i,j);
      ++ctr;
    }
  }
  
  unsigned int channel,plane,wire,wire2;
  for(int p = 0; p < geom->Nplanes(); p++) 
    {
      edm::PtrVector<recob::Hit> vHits;
      edm::PtrVectorItr<recob::Cluster> clusterIter = clusIn.begin();
      hit.clear();
      cHits.clear();      
      while(clusterIter!= clusIn.end() ) {
	cHits = (*clusterIter)->Hits(p,-1);
	if(cHits.size() > 0)
	  //hit.insert(hit.end(),cHits.begin(),cHits.end());
	  edm::PtrVectorItr<recob::Hit> hitIter = cHits.begin();
	clusterIter++;  
      } 
      
      if(hit.size() == 0) 
        continue;

      numberwires=geom->Nwires(0);
      numbertimesamples=hit[0]->Wire()->fSignal.size();
      double MatrixAsum[numberwires][fTimeBins];
      double MatrixBsum[numberwires][fTimeBins];
      double hit_map[numberwires][fTimeBins];//the map of hits 
      int hit_loc[numberwires][fTimeBins];//the index of the hit that corresponds to the potential corner
      double Cornerness[numberwires][fTimeBins];//the "weight" of a corner
  
      for(int wi=0;wi < numberwires; wi++)
	for(int timebin=0;timebin < fTimeBins; timebin++)
	  {
	    hit_map[wi][timebin]=0.;
	    hit_loc[wi][timebin]=-1;
	    Cornerness[wi][timebin]=0.;
	    MatrixAsum[wi][timebin]=0.;
	    MatrixBsum[wi][timebin]=0.;
	  }
             
      for(unsigned int i=0;i < hit.size(); i++)
	{
	  channel=hit[i]->Wire()->RawDigit()->Channel();
	  geom->ChannelToWire(channel,plane,wire);
	  //pixelization using a Gaussian
	  for(int j=0;j <= (int)(hit[i]->EndTime()-hit[i]->StartTime()+.5); j++)    
	    hit_map[wire][(int)((hit[i]->StartTime()+j)*(fTimeBins/numbertimesamples)+.5)]+=Gaussian((int)(j-((hit[i]->EndTime()-hit[i]->StartTime())/2.)+.5),0,hit[i]->EndTime()-hit[i]->StartTime());      
	}
	
      ////Gaussian derivative convolution  
      for(int wire=1;wire < numberwires-1; wire++)
	for(int timebin=1;timebin < fTimeBins-1; timebin++)
	  {
	    MatrixAsum[wire][timebin]=0.;
	    MatrixBsum[wire][timebin]=0.;
	    n=0;
	    for(int i = -3; i <= 3; i++) 
              {
		windex=0;
                for(int j = -3; j <= 3; j++) 
		  {
		    tindex=0;
		    while(wire+i+windex<0)
		      windex++;
		    while(wire+i+windex>numberwires)
		      windex--;
		    while(timebin+j+tindex<0)
		      tindex++;
		    while(timebin+j+tindex>fTimeBins)
		      tindex--;               
		    MatrixAsum[wire][timebin]+=wx[n]*hit_map[wire+i+windex][timebin+j+tindex];  
		    MatrixBsum[wire][timebin]+=wy[n]*hit_map[wire+i+windex][timebin+j+tindex]; 
		    n++;
		  }
              }
	  }
     
      //calculate the cornerness of each pixel while making sure not to fall off the hit map.
      for(int wire=1;wire < numberwires-1; wire++)
	for(int timebin=1;timebin < fTimeBins-1; timebin++)
	  {     
	    MatrixAAsum=0;
	    MatrixBBsum=0;
	    MatrixCCsum=0;
	    //Gaussian smoothing convolution
	    n=0;

	    for(int i = -3; i <= 3; i++) 
              {
		windex=0;
                for(int j = -3; j <= 3; j++) 
		  {
		    tindex=0;
		    while(wire+i+windex<0)
		      windex++;
		    while(wire+i+windex>numberwires)
		      windex--;
		    while(timebin+j+tindex<0)
		      tindex++;
		    while(timebin+j+tindex>fTimeBins)
		      tindex--;               
		    MatrixAAsum+=w[n]*pow(MatrixAsum[wire+i][timebin+j],2);  
		    MatrixBBsum+=w[n]*pow(MatrixBsum[wire+i][timebin+j],2);                   
		    MatrixCCsum+=w[n]*MatrixAsum[wire+i][timebin+j]*MatrixBsum[wire+i][timebin+j]; 
		    n++;
		  }
              }
       
	    if((MatrixAAsum+MatrixBBsum)>0)		
	      Cornerness[wire][timebin]=(MatrixAAsum*MatrixBBsum-pow(MatrixCCsum,2))/(MatrixAAsum+MatrixBBsum);
	    else
	      Cornerness[wire][timebin]=0;
        
	    if(Cornerness[wire][timebin]>0)
	      {	  
		for(unsigned int i=0;i < hit.size(); i++)
		  {
		    channel=hit[i]->Wire()->RawDigit()->Channel();
		    geom->ChannelToWire(channel,plane,wire2);	 
		    //make sure the vertex candidate coincides with an actual hit.
		    if(wire==wire2 && hit[i]->StartTime()<timebin*(numbertimesamples/fTimeBins) 
		       && hit[i]->EndTime()>timebin*(numbertimesamples/fTimeBins))
		      { 	       
			//this index keeps track of the hit number
			hit_loc[wire][timebin]=i;
			Cornerness2.push_back(Cornerness[wire][timebin]);
			break;
		      } 	        
		  }	     
	      }	    
	  }

      std::sort(Cornerness2.rbegin(),Cornerness2.rend());

      for(int vertexnum=0;vertexnum<fMaxCorners;vertexnum++)
	{
	  flag=0;
	  for(int wire=0;wire < numberwires && flag==0; wire++)
	    for(int timebin=0;timebin < fTimeBins && flag==0; timebin++)
	      {    
		if(Cornerness2.size()>(unsigned int)vertexnum)
		  if(Cornerness[wire][timebin]==Cornerness2[vertexnum] 
		     && Cornerness[wire][timebin]>0. 
		     && hit_loc[wire][timebin]>-1)
		    {
		      flag++;      
		      //thresholding
		      if(Cornerness2.size())
			if(Cornerness[wire][timebin]<(fThreshold*Cornerness2[0]))
			  vertexnum=fMaxCorners;
      
		      vHits.push_back(hit[hit_loc[wire][timebin]]);
		      recob::Vertex vertex(vHits);
		      vertex.SetWire(wire);
		      vertex.SetDriftTime(hit[hit_loc[wire][timebin]]->CrossingTime());
		      //weak vertices are given vertex id=1
		      vertex.SetID(1);
		      vertex.SetStrength(Cornerness[wire][timebin]);	  
		      vtxcol->push_back(vertex);
		      vHits.clear();
		      // non-maximal suppression on a square window. The wire coordinate units are 
		      // converted to time ticks so that the window is truly square. 
		      // Note that there are 1/0.0743=13.46 time samples per 4.0 mm (wire pitch in ArgoNeuT), 
		      // assuming a 1.5 mm/us drift velocity for a 500 V/cm E-field 
		      for(int wireout=wire-(int)((fWindow*(numbertimesamples/fTimeBins)*.0743)+.5);
			  wireout <= wire+(int)((fWindow*(numbertimesamples/fTimeBins)*.0743)+.5) ; wireout++)
			for(int timebinout=timebin-fWindow;timebinout <= timebin+fWindow; timebinout++)
			  if(sqrt(pow(wire-wireout,2)+pow(timebin-timebinout,2))<fWindow)//circular window 
			    Cornerness[wireout][timebinout]=0;	  
		    }     
	      }
	}
    
      Cornerness2.clear();
      hit.clear();
      if(clusterIter!=clusIn.end()) clusterIter++;

      if(p==fSaveVertexMap)
	{ 
	  unsigned char *outPix = new unsigned char [fTimeBins*numberwires];
	  //finds the maximum cell in the map for image scaling
	  int cell, pix=0, maxCell=0;
	  int xmaxx, ymaxx;
	  for (int y=0; y<fTimeBins; y++)
	    for (int x=0; x<numberwires; x++)
	      {
		cell = (int)(hit_map[x][y]*1000);
		if (cell > maxCell){
		  maxCell = cell;
		  xmaxx=x;
		  ymaxx=y;
		}
	      }
       
	  for (int y=0; y<fTimeBins; y++)
	    for (int x=0; x<numberwires; x++)
	      { 
		//scales the pixel weights based on the maximum cell value     
		if(maxCell>0)
		  pix = (int)((1500000*hit_map[x][y])/maxCell);
		outPix[y*numberwires + x] = pix;
	      }
	  SaveBMPFile("harrisvertexmap.bmp", outPix, numberwires, fTimeBins);
	  delete [] outPix;
	}  
    
      evt.put(vtxcol); 
    
    
    }
    

}

//-----------------------------------------------------------------------------
//this method saves a BMP image of the vertex map space, which can be viewed with gimp
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



