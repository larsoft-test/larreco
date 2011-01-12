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

/** @addtogroup RKTrackRep
 * @{
 */
 
#ifndef GFMATERIALEFFECTS_H
#define GFMATERIALEFFECTS_H

#include<iostream>
#include"TObject.h"
#include<vector>
#include"TVector3.h"
#include"TMatrixT.h"

#include"GFAbsEnergyLoss.h"
#include"GFGeoMatManager.h"

  
/** @brief  Handles energy loss classes. Contains stepper and energy loss/noise matrix calculation
 * 
 *
 *  @author Christian H&ouml;ppner (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Sebastian Neubert  (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Johannes Rauch  (Technische Universit&auml;t M&uuml;nchen, author)
 * 
 *  This class handles the different energy loss classes that inherit from GFAbsEnergyLoss. 
 *  It provides functionality to limit the stepsize of an extrapolation in order not to
 *  exceed a specified maximum momentum loss. After propagation, the energy loss 
 *  for the given length and (optionally) the noise matrix can be calculated.
 *  
 *  
 */

namespace genf {

class GFMaterialEffects : public TObject{ 
 public:
    
  GFMaterialEffects();
  virtual ~GFMaterialEffects();

  //! Calculates energy loss in the travelled path, optional calculation of noise matrix
  double effects(const std::vector<TVector3>& points, 
                 const std::vector<double>& pointPaths, 
                 const double& mom,
                 const int& pdg,
                 const bool& doNoise = false,
                       TMatrixT<Double_t>* noise = NULL,
                 const TMatrixT<Double_t>* jacobian = NULL,
                 const TVector3* directionBefore = NULL, 
                 const TVector3* directionAfter = NULL);

  //! Returns maximum length so that a specified momentum loss will not be exceeded
  /**  The stepper returns the maximum length that the particle may travel, so that a specified relative momentum loss will not be exceeded.
  */
  double stepper(const double& maxDist,
                 const double& posx,
                 const double& posy,
                 const double& posz,
                 const double& dirx,
                 const double& diry,
                 const double& dirz,
                 const double& mom,
                 const int& pdg);
  double stepper(const double& maxDist,
                 const TVector3& pos, 
                 const TVector3& dir,
                 const double& mom,
                 const int& pdg){
    return stepper(maxDist, pos.X(),pos.Y(),pos.Z(),dir.X(),dir.Y(),dir.Z(),mom,pdg);};

 private:
  //! contains energy loss classes
  std::vector<GFAbsEnergyLoss*> fEnergyLoss;
  //! interface to material and geometry
  GFGeoMatManager *geoMatManager; 
 
  // public:
  //lassDef(GFMaterialEffects,1)

};

} // end namespace

#endif

/* @} **/

