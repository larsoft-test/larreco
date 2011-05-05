////////////////////////////////////////////////////////////////////////
//
// MuonFilter class
//
// This event filter can act to identify events with only through-going tracks
//
// pagebri3@msu.edu
//
////////////////////////////////////////////////////////////////////////
extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}
#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include "TMath.h"

//Framework Includes
#include "art/Framework/Core/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Persistency/Common/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Core/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 


//Larsoft Includes
#include "MuonFilter.h"
#include "RecoBase/recobase.h"
#include "Geometry/geo.h"
#include "Geometry/WireGeo.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"

namespace filt{



//-------------------------------------------------
MuonFilter::MuonFilter(fhicl::ParameterSet const & pset) 
{
    this->reconfigure(pset);   
}

//-------------------------------------------------
MuonFilter::~MuonFilter()
{
}

//------------------------------------------------
void MuonFilter::reconfigure(fhicl::ParameterSet p)
{
 fClusterModuleLabel = p.get< std::string  >("ClusterModuleLabel");
 fLineModuleLabel    = p.get< std::string  >("LineModuleLabel");
 fTolerance          = p.get< double       >("Tolerance");
 fDelay              = p.get< double       >("Delay");
 fDCenter            = p.get< double       >("DCenter");
 fMaxIon             = p.get< double       >("MaxIon");
 fIonFactor          = p.get< double       >("IonFactor");
 fCuts               = p.get< std::vector<double> >("Cuts");
 
}

void MuonFilter::beginJob()
{
}

void MuonFilter::endJob()
{
}

//-------------------------------------------------
bool MuonFilter::filter(art::Event &evt)
{ 
  //numEventHist->Fill(0);
  int event = evt.id().event(); 
  art::ServiceHandle<geo::Geometry> geom;
  art::ServiceHandle<util::LArProperties> larprop;
  art::ServiceHandle<util::DetectorProperties> detprop;
  double drift = larprop->DriftVelocity(larprop->Efield(), larprop->Temperature())*detprop->SamplingRate()/1000.0;
  //This code only works comparing 2 planes so for now these are the 
  // last induction plane and collection plane
  int vPlane = geom->Nplanes() -1;
  geo::View_t vView = geom->Plane(vPlane).View();
  int uPlane = vPlane-1;
  geo::View_t uView = geom->Plane(uPlane).View();
  double uNumWires = ((double)(geom->Plane(uPlane).Nwires()-1))/2.0;
  double vNumWires = ((double)(geom->Plane(vPlane).Nwires()-1))/2.0;
  double uTheta = geom->Plane(uPlane).Wire(0).ThetaZ();
  double vTheta = geom->Plane(vPlane).Wire(0).ThetaZ();
  //This assumes that the 2 planes are reflections of one another across z
  double theta = 0.5*(TMath::Pi()-TMath::Abs(uTheta-vTheta));
  double temp[3];
  geom->Plane(uPlane).Wire(20).GetCenter(temp,0.0);
  TVector3 p1(temp);
  geom->Plane(uPlane).Wire(21).GetCenter(temp,0.0);
  TVector3 p2(temp);
  geom->Plane(uPlane).Wire(20).GetCenter(temp,1.0);
  TVector3 p3(temp);
  double spacing = TMath::Abs((((p1-p3).Cross((p1-p3).Cross(p1-p2))).Unit()).Dot(p1-p2));
  //std::cout << "Spacing " << spacing << std::endl;
  art::Handle< std::vector< recob::Cluster > > clustHandle;
  evt.getByLabel(fClusterModuleLabel,clustHandle);
  art::PtrVector<recob::Cluster> clusters;
  for(unsigned int i = 0; i < clustHandle->size(); ++i) {
    art::Ptr<recob::Cluster> prod(clustHandle,i);
    clusters.push_back(prod);
  }
  double indIon(0),colIon(0);
  double indDEDX(0), colDEDX(0);
  unsigned int channel;
  int plane, wire;
  std::map<int,int> indMap;
  std::map<int,int> colMap;
  std::vector<std::pair<int,int> > rLook;
  int matchNum=1;
  std::vector<std::vector<double> > tGoing;
  std::vector<std::vector<double> > matched;
  std::vector<double> pointTemp(6);
  std::pair<int,int> pairTemp;
  double ionSum(0.0);
  for(unsigned int cluster = 0; cluster < clusters.size(); cluster++) {
    ionSum=0.0;
    for(unsigned int hit = 0; hit < clusters[cluster]->Hits().size(); hit++) {
      ionSum+=clusters[cluster]->Hits()[hit]->Charge(true);
    }
    if(clusters[cluster]->View() == uView) indIon+=ionSum;
    else if(clusters[cluster]->View() == vView) colIon+=ionSum;
  }
  std::cout << "Ionizations: " << indIon << " " << colIon << std::endl;
  art::Handle<std::vector<recob::Cluster > > lines; 
  art::PtrVector<recob::Cluster> inductionSegments, collectionSegments;
  evt.getByLabel(fLineModuleLabel,lines);
  art::PtrVector<recob::Cluster> lineVec;
  for(unsigned int i = 0; i < lines->size(); ++i) {
    art::Ptr<recob::Cluster> prod(lines,i);
    lineVec.push_back(prod);
  }
  for(unsigned int cl = 0;cl < clusters.size(); cl++) {
    if (clusters[cl]->Hits().size()>0 && clusters[cl]->View()==uView)
      inductionSegments.push_back(clusters[cl]);
    else if(clusters[cl]->Hits().size()>0 && clusters[cl]->View() == vView) collectionSegments.push_back(clusters[cl]);
  } 
  if(inductionSegments.size() ==0 || collectionSegments.size() == 0) { 
    //resultTable->Fill(event, 1);
    std::cout << "At least one plane with no track"<<std::endl;
  }
  else {  
    double x1,x2,y1,y2,z1,z2,uPos1,vPos1,uPos2,vPos2;
    for(int i = 0; i < inductionSegments.size(); i++) { 
      if(indMap[i]) continue;
      for(int j = 0; j < collectionSegments.size(); j++) {
          if(colMap[j]) continue;	

	  art::Ptr<recob::Cluster>  indSeg = inductionSegments[i];
	  art::Ptr<recob::Cluster>  colSeg = collectionSegments[j];
          
	  double trk1Start = indSeg->StartPos()[1]+fDelay;
	  double trk1End = indSeg->EndPos()[1]+fDelay;
	  double trk2Start =colSeg->StartPos()[1];
	  double trk2End =colSeg->EndPos()[1];
	  
	  uPos1 = ((double)indSeg->StartPos()[0])-uNumWires;
	  uPos2 = ((double)indSeg->EndPos()[0])-uNumWires; 
	  vPos1 = ((double) colSeg->StartPos()[0])-vNumWires;
	  vPos2 = ((double)colSeg->EndPos()[0])-vNumWires;
	  std::cout << "I J " << i <<" " << j <<std::endl;
	  std::cout << "Start/end " << indSeg->StartPos()[0] <<" "<< colSeg->StartPos()[0] <<" "<< indSeg->EndPos()[0] <<" "<< colSeg->EndPos()[0] <<std::endl;
          std::cout<<"U's "<< uPos1 <<" " << uPos2 <<"V's "<< vPos1 <<" " << vPos2 << " times " << trk1End <<" "<< trk2End <<" "<< trk1Start <<" "<< trk2Start << std::endl;
	  if((TMath::Abs(uPos1-vPos1)>51||TMath::Abs(uPos2-vPos2)>51) &&
	     (TMath::Abs(uPos1-vPos2)<=51&&TMath::Abs(uPos2-vPos1)<=51)) {
	    std::cout << "Swapped1" << std::endl;
	    Swap(uPos1,uPos2);
	  }
	  
          if((TMath::Abs(trk1Start-trk2Start) > fTolerance && TMath::Abs(trk1End-trk2End) > fTolerance) &&  (TMath::Abs(trk1Start-trk2End) < fTolerance && TMath::Abs(trk1End-trk2Start) < fTolerance)) {
            Swap(trk1Start,trk1End);
            Swap(uPos1,uPos2);
            std::cout << "Swapped2" << std::endl;
          }
	  std::cout << "Times: " << trk1Start <<" "<< trk2Start <<" "<<trk1End <<" "<<trk2End<<std::endl;
          if((TMath::Abs(trk1Start-trk2Start) < fTolerance && TMath::Abs(trk1End-trk2End) < fTolerance) && (TMath::Abs(uPos1-vPos1) <=53 && TMath::Abs(uPos2-vPos2) <= 53))  {
	    z1 = spacing/(2*TMath::Sin(theta))*(uPos1+ vPos1);
	    z2 = spacing/(2*TMath::Sin(theta))*(uPos2+ vPos2);
	    y1 = spacing/(2*TMath::Cos(theta))*(uPos1-vPos1)+1;
	    y2 = spacing/(2*TMath::Cos(theta))*(uPos2-vPos2)+1;
	    x1 = (trk1Start+trk2Start)/2.0*drift-fDCenter;
	    x2 = (trk1End+trk2End)/2.0*drift-fDCenter;
	    std::cout <<"Match " << matchNum <<" "<< x1 << " "<< y1 << " "<< z1 << " "<< x2 << " "<< y2 << " "<< z2 << std::endl;
	    bool x1edge,x2edge,y1edge, y2edge,z1edge,z2edge;
	    indMap[i]=matchNum;
	    colMap[j]=matchNum;
	    matchNum++;
	    pointTemp[0]=x1;
	    pointTemp[1]=y1;
	    pointTemp[2]=z1;
	    pointTemp[3]=x2;
	    pointTemp[4]=y2;
	    pointTemp[5]=z2;	
	    x1edge =(TMath::Abs(x1) -fCuts[0] > 0);
	    x2edge =(TMath::Abs(x2) -fCuts[0] > 0);
	    y1edge =(TMath::Abs(y1) -fCuts[1] > 0);
	    y2edge =(TMath::Abs(y2) -fCuts[1] > 0);
	    z1edge =(TMath::Abs(z1) -fCuts[2] > 0);
	    z2edge =(TMath::Abs(z2) -fCuts[2] > 0);
	    if((x1edge||y1edge||z1edge) && (x2edge||y2edge||z2edge)) 
	      {
                tGoing.push_back(pointTemp);
		std::cout << "outside   Removed induction ion: ";          
		for(unsigned int h = 0; h < indSeg->Hits().size(); h++)                {
		  std::cout<<indSeg->Hits()[h]->Charge(true)<< " ";
                  indIon -= indSeg->Hits()[h]->Charge(true);
		}
		std::cout << std::endl <<"Removed collection ion: ";
                for(unsigned int h = 0; h < colSeg->Hits().size(); h++)
		  { 
		    std::cout <<colSeg->Hits()[h]->Charge(true)<< " ";
		    colIon -= colSeg->Hits()[h]->Charge(true);
		  }
		std::cout <<std::endl<<"Ionization outside track I/C: " << indIon << " "<<colIon<< std::endl;	 
	      }
            else if((x1edge || y1edge || z1edge) && !(x2edge || y2edge || z2edge)   && (z2-z1)>1.2) 
	      {
                tGoing.push_back(pointTemp);
		std::cout << "stopping   Removed induction ion: ";          
		for(unsigned int h = 0; h < indSeg->Hits().size(); h++)                {
		  std::cout<<indSeg->Hits()[h]->Charge(true)<< " ";
                  indIon -= indSeg->Hits()[h]->Charge(true);
		}
		std::cout << std::endl <<"Removed collection ion: ";
                for(unsigned int h = 0; h < colSeg->Hits().size(); h++)
		  { 
		    std::cout <<colSeg->Hits()[h]->Charge(true)<< " ";
		    colIon -= colSeg->Hits()[h]->Charge(true);
		  }
		std::cout <<std::endl<<"Ionization outside track I/C: " << indIon << " "<<colIon<< std::endl;
	      }
	    else {
	      pairTemp = std::make_pair(i,j);
	      std::cout << "rLook matchnum " <<matchNum << " "<<i << " "<<j <<std::endl; 
	      rLook.push_back(pairTemp);
	      matched.push_back(pointTemp);  
            }
            break;  //advances i, makes j=0;           	 
	  }
      }
    }
  }
  double distance=0;
  for(int i = 0; i < tGoing.size(); i++)for(int j = 0; j < matched.size();j++){
      std::cout << tGoing.size() <<" " << matched.size() <<" "<<i <<" "<<j<< std::endl;
      //test if one is contained within the other in the z-direction 
      if((tGoing[i][2] <= matched[j][2]) && 
	 (tGoing[i][5] >= matched[j][5])) { 
	TVector3 a1(&tGoing[i][0]);
	TVector3 a2(&tGoing[i][3]);
	TVector3 b1(&matched[j][0]);
	distance = TMath::Abs((((a1-a2).Cross((a1-a2).Cross(a1-b1))).Unit()).Dot(a1-b1));
        std::cout <<"distance "<< distance<< std::endl;
        if(distance < 6){ 
	  std::cout <<"Removing delta ion "<<rLook.size()<<" "<<rLook[j].first<<" "<<matchNum<< std::endl;
	  art::PtrVector<recob::Hit > temp = inductionSegments[rLook[j].first]->Hits();
          for(int h = 0; h < temp.size();h++)	 
	    indIon -= temp[h]->Charge(true);
          temp = collectionSegments[rLook[j].second]->Hits();
	  for(int h = 0; h < temp.size();h++)	   
	    colIon -= temp[h]->Charge(true);   
	}
      }
    } 
  std::cout <<"indIon "<<indIon*fIonFactor <<" colIon " << colIon << std::endl;
  if((indIon*fIonFactor > fMaxIon) && (colIon > fMaxIon)) {
    //resultTable->Fill(0);
    //numEventHist->Fill(1);
    //totIonSelHist->Fill(indIon, colIon);
    //seldEdXHist->Fill(indDEDX, colDEDX);
    return true;
  }
  //totIonRejHist->Fill(indIon, colIon);
  //rejdEdXHist->Fill(indDEDX, colDEDX);
  return false;
}
 
void MuonFilter::Swap(double & x, double &y)  {

  double temp;
  temp = x;
  x=y;
  y=temp;

  return;
}

}
