////////////////////////////////////////////////////////////////////////
//
// VertexMatch class
//
// joshua.spitz@yale.edu
//
//  This algorithm is designed to match vertices found with a dedicated vertex finder 
//  (HarrisVertexFinder) and those found with the HoughLineFinder. A weak vertex is a 
//  vertex that has been found using a dedicated vertex finding algorithm only. A strong 
//  vertex is a vertex that has been found using a dedicated vertex finding algorithm and 
//  matched to a crossing of two or more HoughLineFinder lines.
////////////////////////////////////////////////////////////////////////

#include "VertexMatch.h"
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


vertex::VertexMatch::VertexMatch(edm::ParameterSet const& pset) : 
  fVertexModuleLabel(pset.getParameter< std::string >("VertexModuleLabel")),
  fHoughModuleLabel (pset.getParameter< std::string >("HoughModuleLabel")),
  fMaxDistance      (pset.getParameter< double      >("MaxDistance"))
{
  produces< std::vector<recob::Vertex> >();
}

vertex::VertexMatch::~VertexMatch()
{
}

bool sort_pred(const std::pair<edm::Ptr<recob::Hit>,double>& left, const std::pair<edm::Ptr<recob::Hit>,double>& right)
{
  return left.first < right.first;
}

bool sort_pred2(const std::pair<edm::Ptr<recob::Hit>,double>& left, const std::pair<edm::Ptr<recob::Hit>,double>& right)
{
  return left.second < right.second;
}

void vertex::VertexMatch::produce(edm::Event& evt, edm::EventSetup const&)
{

  edm::Handle< std::vector<recob::Vertex> > vertexListHandle;
  evt.getByLabel(fVertexModuleLabel,vertexListHandle);
  
  edm::Handle< std::vector<recob::Cluster> > houghListHandle;
  evt.getByLabel(fHoughModuleLabel,houghListHandle);
  
  std::auto_ptr<std::vector<recob::Vertex> > mvertexcol(new std::vector<recob::Vertex>);
  
  edm::Service<geo::Geometry> geom;
  //hits associated with a vertex
  edm::PtrVector<recob::Hit> vHits;
  edm::PtrVector<recob::Hit> vertexhit;
  
  std::vector<double> weakvertexstrength; //strength of weak vertices
  std::vector<double> strongvertexstrength; //strength of strong vertices
  //hits associated with a hough line
  edm::PtrVector<recob::Hit> hHits;
  edm::PtrVector<recob::Hit> houghhit;

  //edm::PtrVector< std::pair<recob::Hit,double> > matchedvertex;//vertices associated with a hough line
  
  std::vector< std::pair<edm::Ptr<recob::Hit>, double> > matchedvertex;

  edm::PtrVector<recob::Hit> strongvertex;
  std::vector< std::pair<edm::Ptr<recob::Hit>, double> > strongestvertex; //the strongest strong vertex
  
  edm::PtrVector<recob::Vertex> vertIn;
  for(unsigned int ii = 0; ii < vertexListHandle->size(); ++ii)
    {
      edm::Ptr<recob::Vertex> vertex(vertexListHandle, ii);
      vertIn.push_back(vertex);
    }
  edm::PtrVector<recob::Cluster> houghIn;
  for(unsigned int ii = 0; ii < houghListHandle->size(); ++ii)
    {
      edm::Ptr<recob::Cluster> cluster(houghListHandle, ii);
      houghIn.push_back(cluster);
    }
  unsigned int channel,plane,wire;
  double slope,intercept,distance;
  double starttime, endtime;
  int startwire, endwire;
  double strength; //the strength of a strong vertex
  for(int p = 0; p < geom->Nplanes(); p++) {
    //create the vector of vertex hits 
    edm::PtrVectorItr<recob::Vertex> vertexIter = vertIn.begin();
    edm::PtrVectorItr<recob::Cluster> houghIter = houghIn.begin();
    while(vertexIter!= vertIn.end() ) 
      {
	// vHits = (*vertexIter)->Hits(p,-1);
	//      if(vHits.size() > 0){
	//  	 vertexhit.insert(vertexhit.end(),vHits.begin(),vHits.end());
	//  	 weakvertexstrength.push_back((*vertexIter)->Strength());
	//  	 }
	vHits = (*vertexIter)->Hits(p);
	if(vHits.size() > 0)
	  {
	    edm::PtrVectorItr<recob::Hit> vertexhitIter = vHits.begin();
	    while (vertexhitIter!=vHits.end())
	      {
		vertexhit.push_back((*vertexhitIter));
		vertexhitIter++;
	      }       
	  }
	vertexIter++;

      } 
    vHits.clear();
      
//     if(vHits.size() == 0)
//       continue;
      
    //loop over vector of hough lines and find the vertex hits that are associated with the hough line(s)
    houghIter = houghIn.begin();  
    while(houghIter!= houghIn.end()) 
      { 
	houghhit.clear();
	hHits.clear();
	plane=-1; 
	distance=-1.;
	//create vector of hits associated with hough line
	hHits = (*houghIter)->Hits(p);
        if(hHits.size() > 0)
	  {

 	    edm::PtrVectorItr<recob::Hit> hitIter = hHits.begin();
	    while (hitIter!=hHits.end())
	      {
		houghhit.push_back((*hitIter));
		hitIter++;
	      }  
	  }
        
        if(houghhit.size()){
	  channel=houghhit[0]->Wire()->RawDigit()->Channel();
	  geom->ChannelToWire(channel,plane,wire);
	}
	if(p==plane)
	  {
	    slope=(*houghIter)->Slope();
	    intercept=(*houghIter)->Intercept();    

	    for(unsigned int i=0;i < vertexhit.size(); i++)
	      {   
		distance=-1;
		channel=vertexhit[i]->Wire()->RawDigit()->Channel();
		geom->ChannelToWire(channel,plane,wire);

		starttime=(*houghIter)->StartTime();
		endtime=(*houghIter)->EndTime();
		startwire=(*houghIter)->StartWire();
		endwire=(*houghIter)->EndWire();
		//require the vertices found with HarrisVertexFinder to match up with the endpoints 
		//(within a window) of a Hough line. A strong vertex matches up with at least two Hough lines. 
		if(((TMath::Abs((int)(wire-startwire))<fMaxDistance*.0743)
		    ||(TMath::Abs((int)(wire-endwire))<fMaxDistance*.0743)
		   )
		   &&((TMath::Abs(vertexhit[i]->CrossingTime()-starttime)<fMaxDistance)
		      ||(TMath::Abs(vertexhit[i]->CrossingTime()-endtime)<fMaxDistance)
		     )
		   )          
		  distance=(TMath::Abs(vertexhit[i]->CrossingTime()-slope*(double)wire-intercept)/(sqrt(pow(.0743*slope,2)+1))); 
   
		if(distance<(fMaxDistance+((vertexhit[i]->EndTime()-vertexhit[i]->StartTime())/2.))&&distance>-1)
		  matchedvertex.push_back(std::pair<edm::Ptr<recob::Hit>,double>(vertexhit[i], weakvertexstrength[i]*sqrt(pow(TMath::Abs(endwire-startwire)*.0743,2)+pow(TMath::Abs(endtime-starttime),2))));
		//ala strongestvertex.push_back(std::pair<edm::PtrVector<recob::Hit>,double>(matchedvertex[i].first,strength));

	      }
	  }
	   
	if(vertexhit.size() == 0 || houghhit.size() == 0) 
	  {
	    houghIter++;
	    continue;
	  }

	if(vertexIter!=vertIn.end()) vertexIter++;
	if(houghIter!=houghIn.end()) houghIter++;    
      }
   
    //sort matchedvertex vector to make it easy to find duplicate entries (strong vertices)
    std::sort(matchedvertex.rbegin(), matchedvertex.rend(),sort_pred);

    // the "strength" of a strong vertex is defined as 
    // (HarrisVertexStrength*LengthofHoughLine)_1+(HarrisVertexStrength*LengthofHoughLine)_2+...
    // ...+(HarrisVertexStrength*LengthofHoughLine)_n, where n is the number of vertices 
    // associated with a Hough Line
  
    for(unsigned int i=0;i < matchedvertex.size(); i++) 
      {
	strength=matchedvertex[i].second;
  
	for(unsigned int n=1;n < matchedvertex.size() && i>=n; n++)
	  if(matchedvertex[i].first==matchedvertex[i-n].first)	
	    strength+=matchedvertex[i-n].second;
  
	strongvertexstrength.push_back(strength);
	//make sure there is more than one Hough Line associated with the vertex 
          
	if(strength>matchedvertex[i].second)
	  strongestvertex.push_back(std::pair<edm::Ptr<recob::Hit>,double>(matchedvertex[i].first,strength));
      }
  
  
    //sort the strength of the strong vertices to find the strongest vertex
    std::sort(strongestvertex.rbegin(), strongestvertex.rend(),sort_pred2);

    for(unsigned int i=0;i < matchedvertex.size(); i++)
      {
	// I think this is grabbing first item in pair, itself a pointer then grabbing first 
	// (.begin()) one of those. EC, 18-Oct-2010.
	channel=(matchedvertex[i].first)->Wire()->RawDigit()->Channel();
	geom->ChannelToWire(channel,plane,wire);

	// strongvertex, despite name, is a hit vector.
	strongvertex.push_back(matchedvertex[i].first);
	recob::Vertex vertex(strongvertex);      
	vertex.SetWireNum(wire);
	vertex.SetStrength(strongvertexstrength[i]);
	//vertex->SetDriftTime(matchedvertex[i].first->CrossingTime());
	vertex.SetDriftTime((matchedvertex[i].first)->CrossingTime());

	//find the strong vertices, those vertices that have been associated with more than one hough line  
	if(matchedvertex[i].first==matchedvertex[i-1].first)	
	  {
	    if(strongvertex[0]==(strongestvertex[0].first)&&strongestvertex.size()>0)
	      vertex.SetID(4);//the strongest strong vertex is given a vertex id=4
	    else
	      vertex.SetID(3);//strong vertices are given vertex id=3
	  } 	  
	else	
	  {
	    vertex.SetID(2);//weak vertices that have been associated with an endpoint of a single Hough line are given vertex id=2
	  } 
	  
	mvertexcol->push_back(vertex);
	strongvertex.clear();
	  	 	  
      }

    strongestvertex.clear();
    matchedvertex.clear();
    vertexhit.clear();  
  }
    

  evt.put(mvertexcol);

}


