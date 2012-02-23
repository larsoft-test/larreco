/* Copyright 2008-2010, Technische Universitaet Muenchen,
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
#include "GFSpacepointHitPolicy.h"

#include "assert.h"

#include "TMath.h"

#include "GFAbsRecoHit.h"

const std::string genf::GFSpacepointHitPolicy::fPolicyName = "GFSpacepointHitPolicy";

TMatrixT<Double_t> 
genf::GFSpacepointHitPolicy::hitCoord(GFAbsRecoHit* hit,const GFDetPlane& plane)
{
  TMatrixT<Double_t> returnMat(2,1);


  TMatrixT<Double_t> _D(3,1);
  TVector3 _U;
  TVector3 _V;

  _D[0][0] = (plane.getO())[0];
  _D[1][0] = (plane.getO())[1];
  _D[2][0] = (plane.getO())[2];

  _D *= -1.; 
  _D += hit->getRawHitCoord();
  //now the vector _D points from the origin of the plane to the hit point
  // This is TVector3(zero) for us. EC.

  _U = plane.getU();
  _V = plane.getV();


  returnMat[0][0] = _D[0][0] * _U[0] + _D[1][0] * _U[1] + _D[2][0] * _U[2];
  returnMat[1][0] = _D[0][0] * _V[0] + _D[1][0] * _V[1] + _D[2][0] * _V[2];
  //std::cout << "hitCoord="<<std::endl;
  //returnMat.Print();
  return returnMat;
}


TMatrixT<Double_t> 
genf::GFSpacepointHitPolicy::hitCoord(GFAbsRecoHit* hit,const GFDetPlane& plane,const GFDetPlane& planePrev)
{
  //  TMatrixT<Double_t> returnMat(2,1);
  TMatrixT<Double_t> returnMat(5,1);

  TMatrixT<Double_t> _D(3,1);
  TVector3 _U;
  TVector3 _V;
  TVector3 _Norm;
  TVector3 _NormPrev;

  _D[0][0] = (plane.getO())[0];
  _D[1][0] = (plane.getO())[1];
  _D[2][0] = (plane.getO())[2];

  _D *= -1.; 
  _D += hit->getRawHitCoord();
  //now the vector _D points from the origin of the plane to the hit point
  // This is TVector3(zero) for us. EC.

  _U = plane.getU();
  _V = plane.getV();
  _Norm = plane.getNormal();
  _NormPrev = planePrev.getNormal();

  returnMat[0][0] = _NormPrev.Angle(_Norm); //  I will overwrite in GFKalman.
  returnMat[1][0] = 0.; // certainly wrong, but I will overwrite in GFKalman.
  returnMat[2][0] = 0.; // ibid
  returnMat[3][0] = _D[0][0] * _U[0] + _D[1][0] * _U[1] + _D[2][0] * _U[2];
  returnMat[4][0] = _D[0][0] * _V[0] + _D[1][0] * _V[1] + _D[2][0] * _V[2];
  //std::cout << "hitCoord="<<std::endl;
  //returnMat.Print();
  return returnMat;
}


TMatrixT<Double_t> 
genf::GFSpacepointHitPolicy::hitCov(GFAbsRecoHit* hit,const GFDetPlane& plane)
{
  TVector3 _U;
  TVector3 _V;

  _U = plane.getU();
  _V = plane.getV();

  TMatrixT<Double_t> rawCov = hit->getRawHitCov();

  // This Jacobian concerns the transfrom from planar to detector coords.
  TMatrixT<Double_t> jac(3,2);
  
  // jac = dF_i/dx_j = s_unitvec * t_untivec, with s=u,v and t=x,y,z
  // jac = dF_i/dx_j = s_unitvec * t_untivec, with s=|q|/p,du/dw, dv/dw, u,v and t=th, x,y,z
  
  jac[0][0] = _U[0];
  jac[1][0] = _U[1];
  jac[2][0] = _U[2];
  jac[0][1] = _V[0];
  jac[1][1] = _V[1];
  jac[2][1] = _V[2];
  
  TMatrixT<Double_t> jac_orig = jac;
  TMatrixT<Double_t> jac_t = jac.T();

  TMatrixT<Double_t> result=jac_t * (rawCov * jac_orig);
  //std::cout << "hitCov="<<std::endl;
  //result.Print();
  return  result;
}

TMatrixT<Double_t> 
genf::GFSpacepointHitPolicy::hitCov(GFAbsRecoHit* hit,const GFDetPlane& plane, const GFDetPlane& planePrev, const TMatrixT<Double_t>& state)
{
  TVector3 _U;
  TVector3 _V;

  _U = plane.getU();
  _V = plane.getV();

  // Eventually getRawHitCov() will pull in Herbs 
  TMatrixT<Double_t> rawCov = hit->getRawHitCov();
  TVector3 tmpRawCov (rawCov[0][0],rawCov[1][1],rawCov[2][2]);
  rawCov.ResizeTo(7,7);
  rawCov[0][0] = tmpRawCov[0]; // x
  rawCov[1][1] = tmpRawCov[1]; // y
  rawCov[2][2] = tmpRawCov[2]; // z
  rawCov[3][3] = pow(0.2,2.0); // Unit Px
  rawCov[4][4] = pow(0.2,2.0); // Unit Py
  rawCov[5][5] = pow(0.2,2.0); // Unit Pz
  rawCov[6][6] = pow(0.1,2.0); // theta. 0.3/3mm, say.

  // This Jacobian concerns the transfrom from planar to detector coords.
  //   TMatrixT<Double_t> jac(3,2);
  TMatrixT<Double_t> jac(7,5); // X,Y,Z,UX,UY,UZ,Theta in detector coords
  
  // jac = dF_i/dx_j = s_unitvec * t_unitvec, with s=|q|/p,du/dw, dv/dw, u,v and t=th, x,y,z

  //Double_t C = 0.0136*sqrt(0.3/14); // for now. EC, 8-Jan-2012.
  // More correctly :
  Double_t dist(0.3); // place holder!!
  Double_t C;
  Double_t p (2.5); // place holder!!!
  p = abs(1/state[0][0]);
  Double_t mass = 0.104; // CHANGE ME!
  Double_t mom = fabs(1.0/state[0][0]);
  Double_t beta = mom/sqrt(mass*mass+mom*mom);
  dist = (plane.getO()-planePrev.getO()).Mag();
  C = 0.0136/beta*sqrt(dist/14.0)*(1+0.038*log(dist/14.0)); 
  TVector3 _D (plane.getNormal());

  TVector3 u=plane.getU();
  TVector3 v=plane.getV();
  TVector3 w=u.Cross(v);

  // Below line has been done, with code change in GFKalman that has updated
  // the plane orientation by now.
  //  TVector3 pTilde = 1.0 * (w + state[1][0] * u + state[2][0] * v);
  TVector3 pTilde = w;
  double pTildeMag = pTilde.Mag();

  jac.Zero();
  jac[6][0] = 1.0; // C; // 1/C;

  jac[0][3] = _U[0];
  jac[1][3] = _U[1];
  jac[2][3] = _U[2];

  jac[0][4] = _V[0];
  jac[1][4] = _V[1];
  jac[2][4] = _V[2];

  // cnp'd from RKTrackRep.cxx, line ~496
  // da{x,y,z}/du'
  jac[3][1] = 1.0/pTildeMag*(u.X()-pTilde.X()/(pTildeMag*pTildeMag)*u*pTilde);
  jac[4][1] = 1.0/pTildeMag*(u.Y()-pTilde.Y()/(pTildeMag*pTildeMag)*u*pTilde);
  jac[5][1] = 1.0/pTildeMag*(u.Z()-pTilde.Z()/(pTildeMag*pTildeMag)*u*pTilde);
  // da{x,y,z}/dv'
  jac[3][2] = 1.0/pTildeMag*(v.X()-pTilde.X()/(pTildeMag*pTildeMag)*v*pTilde);
  jac[4][2] = 1.0/pTildeMag*(v.Y()-pTilde.Y()/(pTildeMag*pTildeMag)*v*pTilde);
  jac[5][2] = 1.0/pTildeMag*(v.Z()-pTilde.Z()/(pTildeMag*pTildeMag)*v*pTilde);

  TMatrixT<Double_t> jac_orig = jac;
  TMatrixT<Double_t> jac_t = jac.T();

  TMatrixT<Double_t> result=jac_t * (rawCov * jac_orig);
  //std::cout << "hitCov="<<std::endl;
  //result.Print();
  return  result;
}

const genf::GFDetPlane&
genf::GFSpacepointHitPolicy::detPlane(GFAbsRecoHit* hit, GFAbsTrackRep* rep)
{
  TMatrixT<Double_t> rawcoord = hit->getRawHitCoord();
  TVector3 point(rawcoord[0][0],rawcoord[1][0],rawcoord[2][0]);

  TVector3 poca,dirInPoca;
  rep->extrapolateToPoint(point,poca,dirInPoca);

  fPlane.setO(point);
  fPlane.setNormal(dirInPoca);

  return fPlane;
}

//ClassImp(GFSpacepointHitPolicy)
