////////////////////////////////////////////////////////////////////////
//
// AxisRecoAlg class
//
// Andrzej.Szelc@yale.edu
//
//  This algorithm is designed to reconstruct the 3D angle of a track or shower based on the 2D angles of planes. 
//  An expanded version of the formulas from Ornella and Maddalena.
////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <vector>
extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}
#include <sstream>
#include <fstream>
#include <math.h>
#include <algorithm>

#include "TMath.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 
// Framework includes
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "ShowerFinder/AxisRecoAlg.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/AssociationUtil.h"
#include "RecoBase/recobase.h"
#include "Geometry/geo.h"

//-----------------------------------------------------------------------------
shwf::AxisRecoAlg::AxisRecoAlg() 
{
  //this->reconfigure(pset);
}

//-----------------------------------------------------------------------------
shwf::AxisRecoAlg::~AxisRecoAlg()
{
}

//-----------------------------------------------------------------------------
void shwf::AxisRecoAlg::reconfigure(fhicl::ParameterSet const& p)
{

}

//-----------------------------------------------------------------------------
// omega0 and omega1 are calculated as:
//  angle based on distances in wires and time - rescaled to cm.
// tan(angle)*fMean_wire_pitch/(ftimetick*fdriftvelocity);
// as those calculated with Get2Dangle
int shwf::AxisRecoAlg::Get3DaxisN(int iplane0,int iplane1,double omega0, double omega1,double &phi,double &theta){


 double l(0),m(0),n(0);
 std::vector< double > angle;  // this should be geometry nplanes. 
 angle.resize(3);
 
 // pretend collection and induction planes. 
 // "Collection" is the plane with the vertical angle equal to zero. 
 // If both are non zero collection is the one with the negative angle. 
 unsigned int Cplane=0,Iplane=1;   
 //then angleC and angleI are the respective angles to vertical in these planes and slopeC, slopeI are the tangents of the track.
 double angleC,angleI,slopeC,slopeI;
 // don't know how to reconstruct these yet, so exit with error.
 if(omega0==0 || omega1==0){
   phi=-900;
   theta=-900;
   mf::LogWarning("AxisRecoAlg") << " event parallell in one of the planes, exiting ";
   return 0;
 }
 
 //////////insert check for existence of planes.
 
 angle[iplane0]=geom->Plane(iplane0).Wire(0).ThetaZ(false)-TMath::Pi()/2.; // wire angle with respect to 
 angle[iplane1]=geom->Plane(iplane1).Wire(0).ThetaZ(false)-TMath::Pi()/2.; // wire angle with respect to
 
 
 if(angle[iplane0] == 0){   
   // first plane is at 0 degrees
   Cplane=iplane0;
   Iplane=iplane1;
   slopeC=omega0;
   slopeI=omega1;
   std::cout << "+++++ new calc first test case 1 angle[0]==0 "<< std::endl;
 }
 else if(angle[iplane1] == 0){  
   // second plane is at 0 degrees
   Cplane = iplane1;
   Iplane = iplane0;
   slopeC=omega1;
   slopeI=omega0;
   std::cout << "+++++ new calc first test case 2 angle[1]==0 "<< std::endl;
 }
 else if(angle[iplane0] != 0 && angle[iplane1] != 0){
   //both planes are at non zero degree - find the one with deg<0
   if(angle[iplane1] < angle[iplane0]){
     Cplane = iplane1;
     Iplane = iplane0;
     slopeC=omega1;
     slopeI=omega0;
     std::cout << "+++++ new calc first test case 3 angle[1]< angle[0] "<< std::endl;
   }
   else if(angle[iplane1] > angle[iplane0]){
     Cplane = iplane0;
     Iplane = iplane1;
     slopeC=omega0;
     slopeI=omega1;
     std::cout << "+++++ new calc first test case 4 angle[1]> angle[0] "<< std::endl;
   }
   else{
     //throw error - same plane.
     return -1;
   }	

 }
 angleC=angle[Cplane];
 angleI=angle[Iplane];
 
 //0 -1 factor depending on if one of the planes is vertical.
 bool nfact = !(angle[Cplane]);

 l = 1;
 
 m = (1/(2*sin(angleI)))*((cos(angleI)/(slopeC*cos(angleC)))-(1/slopeI) 
				 +nfact*(  cos(angleI)/slopeC-1/slopeI  )     );
 
 n = (1/(2*cos(angleC)))*((1/slopeC)+(1/slopeI) +nfact*((1/slopeC)-(1/slopeI)));
 
 // Direction angles
 phi =  atan(n/l);
 theta = acos( m / (sqrt(pow(l,2)+pow(m,2)+pow(n,2)) ) ) ;
 
 
 // solve the ambiguities due to tangent periodicity
 float Phi = phi > 0. ? (TMath::Pi()/2)-phi : fabs(phi)-(TMath::Pi()/2) ; 
 float Theta = 0;
 if(Phi < 0)Theta = (TMath::Pi()/2)-theta;
 if(Phi > 0)Theta = theta-(TMath::Pi()/2);
 
 phi=Phi*180/TMath::Pi();
 theta=Theta*180/TMath::Pi();
 return 0;
 
}



////////////////////////////
double shwf::AxisRecoAlg::Get2Dangle(double wireend,double wirestart,double timeend,double timestart){

return Get2Dangle(wireend-wirestart,timeend-timestart);
  
}

////////////////////////////
double shwf::AxisRecoAlg::Get2Dangle(double dwire,double dtime){

  double fwirepitch = geom->WirePitch(0,1,0);
  double ftimetick=detp->SamplingRate()/1000.; 
  double fdriftvelocity=larp->DriftVelocity(larp->Efield(),larp->Temperature());
  double BC,AC;
  double omega;
 
  BC = ((double)dwire)*fwirepitch; // in cm
  AC = ((double)dtime)*ftimetick*fdriftvelocity; //in cm 
  omega = asin(  AC/sqrt(pow(AC,2)+pow(BC,2)) );
 
  if(BC<0)  // for the time being. Will check if it works for AC<0
    { 
    if(AC!=0)
      omega= AC/fabs(AC)*TMath::Pi()-omega;  //subtract when negative, add when positive
    else    
      omega=TMath::Pi();
    } 
  
  //return omega;
  return tan(omega)*fwirepitch/(ftimetick*fdriftvelocity);

}

