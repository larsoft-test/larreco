/* Copyright 2008-2009, Technische Universitaet Muenchen,
   Authors: Christian Hoeppner & Sebastian Neubert

   This file is part of GENFIT.

   GENFIT is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   GENFIT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with GENFIT.  If not, see <http://www.gnu.org/licenses/>.
*/

#include"GFMaterialEffects.h"
#include<iostream>
#include "stdlib.h"

#include"math.h"

#include"GFAbsEnergyLoss.h"
#include"GFEnergyLossBetheBloch.h"
#include"GFEnergyLossBrems.h"
#include"GFEnergyLossCoulomb.h"
#include "GFException.h"

#include"GFGeoMatManager.h"


genf::GFMaterialEffects::~GFMaterialEffects(){
  for(unsigned int i=0;i<fEnergyLoss.size();++i) delete fEnergyLoss.at(i);
  delete geoMatManager;
}

genf::GFMaterialEffects::GFMaterialEffects(){
  // append all EnergyLoss classes
  fEnergyLoss.push_back(new GFEnergyLossBetheBloch());
  fEnergyLoss.push_back(new GFEnergyLossBrems());
  fEnergyLoss.push_back(new GFEnergyLossCoulomb());
  
  geoMatManager = new GFGeoMatManager(); // handles material parameters and geometry - GFGeoMatManager uses TGeoMaterial and TGeoManager
}

double genf::GFMaterialEffects::effects(const std::vector<TVector3>& points, 
                                  const std::vector<double>& pointPaths, 
                                  const double& mom,
                                  const int& pdg,
                                  const bool& doNoise,
                                        TMatrixT<Double_t>* noise,
                                  const TMatrixT<Double_t>* jacobian,
                                  const TVector3* directionBefore, 
                                  const TVector3* directionAfter){

  assert(points.size()==pointPaths.size());

  double momLoss=0.;
                     
  for(unsigned int i=1;i<points.size();++i){
    TVector3 p1=points.at(i-1);
    TVector3 p2=points.at(i);
    TVector3 dir=p2-p1;
    double dist=dir.Mag();
    double realPath = pointPaths.at(i);
    
    if (dist > 1.E-8) { // do material effects only if distance is not too small
      dir*=1./dist; //normalize dir

      double X(0.);
      double matDensity, matZ, matA, radiationLength, mEE;
      double step;
      
      geoMatManager->initTrack(p1.X(),p1.Y(),p1.Z(),dir.X(),dir.Y(),dir.Z());

      while(X<dist){
        
        geoMatManager->getMaterialParameters(matDensity, matZ, matA, radiationLength, mEE); 

        step = geoMatManager->stepOrNextBoundary(dist-X);
        
        // Loop over EnergyLoss classes
        if(matZ>1.E-3){
          for(unsigned int j=0;j<fEnergyLoss.size();++j){
            momLoss += realPath/dist*fEnergyLoss.at(j)->energyLoss(step,
                                                                   mom,
                                                                   pdg,              
                                                                   matDensity,
                                                                   matZ,
                                                                   matA,
                                                                   radiationLength,
                                                                   mEE,
                                                                   doNoise,
                                                                   noise,
                                                                   jacobian,
                                                                   directionBefore,
                                                                   directionAfter);
          }
        }
        X += step;
      }
    }
  }
  //std::cout << "GFMaterialEffects::effects(): momLoss " << momLoss << std::endl;
  return momLoss;
}





double genf::GFMaterialEffects::stepper(const double& maxDist,
                                  const double& posx,
                                  const double& posy,
                                  const double& posz,
                                  const double& dirx,
                                  const double& diry,
                                  const double& dirz,
                                  const double& mom,
                                  const int& pdg){

  static const double maxPloss = .005; // maximum relative momentum loss allowed

  geoMatManager->initTrack(posx,posy,posz,dirx,diry,dirz);

  double X(0.);
  double dP = 0.;
  double momLoss = 0.;
  double matDensity, matZ, matA, radiationLength, mEE;
  double step;

  while(X<maxDist){
    
    geoMatManager->getMaterialParameters(matDensity, matZ, matA, radiationLength, mEE); 
    // I need to revisit this. EC, 7-Jan-2011.
    // You know what? F*ck it. Just force this to be LAr.... is what I *could/will* say here ....

    matDensity = 1.40; matZ = 18.0; matA = 39.95; radiationLength=13.947; mEE=188.0;
    if (floor(matZ)!=18.0 and floor(matA)!=39.0) // Non-LAr!
      {
	std::cout << "GFMaterialEffects: matDensity, matZ, matA, radiationLength, mEE" << matDensity<<", "<< matZ<<", "<< matA<<", "<< radiationLength<<", "<< mEE<< std::endl;
	GFException exc("GFMaterialEffects: You have somehow stepped outside the TPC!!!! This is probably very bad.",__LINE__,__FILE__);
	throw exc;
      }
    else
      {
	//	std::cout << "GFMaterialEffects: matDensity, matZ, matA, radiationLength, mEE" << matDensity<<", "<< matZ<<", "<< matA<<", "<< radiationLength<<", "<< mEE<< std::endl;
      }
    //
    step = geoMatManager->stepOrNextBoundary(maxDist-X);
    
    // Loop over EnergyLoss classes
    if(matZ>1.E-3){
      for(unsigned int j=0;j<fEnergyLoss.size();++j){
        momLoss += fEnergyLoss.at(j)->energyLoss(step,
                                                 mom,
                                                 pdg,              
                                                 matDensity,
                                                 matZ,
                                                 matA,
                                                 radiationLength,
                                                 mEE);  
      }
    }
    
    if(dP + momLoss > mom*maxPloss){
      double fraction = (mom*maxPloss-dP)/momLoss;
      assert(fraction>0.&&fraction<1.);
      dP+=fraction*momLoss;
      X+=fraction*step;
      break;
    }
    
    dP += momLoss;
    X += step;
  }

  return X;                 
}

//ClassImp(GFMaterialEffects)

