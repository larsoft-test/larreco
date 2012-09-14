/////////////////////////////////////////////////////////////////////
//
// HoughBaseAlg class
//
// Ben Carls, bcarls@fnal.gov
//
//
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

#include <TF1.h>

#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
 
#include "ClusterFinder/HoughBaseAlg.h"
#include "Filters/ChannelFilter.h"
#include "RecoBase/recobase.h"
#include "Geometry/geo.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/AssociationUtil.h"

//------------------------------------------------------------------------------
cluster::HoughBaseAlg::HoughBaseAlg(fhicl::ParameterSet const& pset)
{
}

//------------------------------------------------------------------------------
cluster::HoughBaseAlg::~HoughBaseAlg()
{
}

//------------------------------------------------------------------------------
void cluster::HoughBaseAlg::reconfigure(fhicl::ParameterSet const& pset)
{ 
  return;
}

//------------------------------------------------------------------------------
cluster::HoughTransform::HoughTransform()
{  
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
int cluster::HoughTransform::AddPointReturnMax(int x, int y, int *yMax, int *xMax, int minHits)
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
void cluster::HoughTransform::Init(int dx, int dy, int rhores,
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
  m_rowLength = (int)(m_rhoResolutionFactor*2. * sqrt(dx*dx + dy*dy));
  
  int angleIndex;
  double a, angleStep = TMath::Pi()/m_numAngleCells;
  for (a=0.0, angleIndex = 0; angleIndex < m_numAngleCells; ++angleIndex){
    m_cosTable[angleIndex] = cos(a);
    m_sinTable[angleIndex] = sin(a);
    a += angleStep;
  }
}


//------------------------------------------------------------------------------
int cluster::HoughTransform::GetMax(int &xmax, int &ymax)
{
  std::map<int,int>::iterator rhoIter;
  int maxVal = -1;
  for(unsigned int i = 0; i < m_accum.size(); i++){
    for(rhoIter=m_accum[i].begin(); rhoIter!=m_accum[i].end(); ++rhoIter){
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
int cluster::HoughTransform::DoAddPointReturnMax(int x, int y, int *ymax, int *xmax, int minHits)
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
        //std::cout << "First, a: " << a << " lastDist: " << lastDist << std::endl;
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
  
  //std::cout << "Add point says xmax: " << *xmax << " ymax: " << *ymax << std::endl;

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
        //std::cout << "First, a: " << a << " lastDist: " << lastDist << std::endl;
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
 


