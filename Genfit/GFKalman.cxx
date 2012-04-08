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
#include "GFKalman.h"

#include "assert.h"
#include <iostream>
#include <sstream>
#include <iomanip>

#include "TMath.h"
#include "TRandom.h"
#include "TDatabasePDG.h"

#include "GFTrack.h"
#include "GFAbsRecoHit.h"
#include "GFAbsTrackRep.h"
#include "GFTrackCand.h"
#include "GFException.h"

#include "cetlib/exception.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 
  
#define COVEXC "cov_is_zero"


genf::GFKalman::GFKalman():fInitialDirection(1),fNumIt(3),fBlowUpFactor(50.),fMomLow(-100.0),fMomHigh(100.0)
{
  art::ServiceHandle<art::TFileService> tfs;
  //fIhitvUpdate = tfs->make<TH1D>("GenfitUpdatevIhit", ";Updatesv Ihit;", 3301, -0.5, 3300.5);
  //fUpdate = tfs->make<TH1D>("GenfitUpdate", ";Updates;",100, -0.1, +0.1);
}

genf::GFKalman::~GFKalman(){;}

void genf::GFKalman::processTrack(GFTrack* trk){
  int direction=fInitialDirection;
  assert(direction==1 || direction==-1);
  //  trk->clearGFBookkeeping();
  trk->clearRepAtHit();
  /*
  int nreps=trk->getNumReps();
  for(int i=0; i<nreps; ++i){ 
    trk->getBK(i)->setNhits(trk->getNumHits());
    trk->getBK(i)->bookMatrices("fPredSt");
    trk->getBK(i)->bookMatrices("fPredCov");
    trk->getBK(i)->bookMatrices("bPredSt");
    trk->getBK(i)->bookMatrices("bPredCov");
    trk->getBK(i)->bookGFDetPlanes("f");
    trk->getBK(i)->bookGFDetPlanes("b");
    trk->getBK(i)->bookNumbers("fChi2");
    trk->getBK(i)->bookNumbers("bChi2");
  }
  */

  /*why is there a factor of two here (in the for statement)?
    Because we consider one full iteration to be one back and
    one forth fitting pass */
  for(int ipass=0; ipass<2*fNumIt; ipass++){
    //    std::cout << "GFKalman:: numreps, numhits, ipass, direction"<< trk->getNumReps()<<", "<< trk->getNumHits()<<", "<< ipass<<", " <<direction << std::endl;
    //    if(floor(ipass/2)==fNumIt && direction==1) blowUpCovsDiag(trk);
    if(ipass>0) blowUpCovs(trk);
    if(direction==1){
      trk->setNextHitToFit(0);
    }
    else {
      trk->setNextHitToFit(trk->getNumHits()-1);
    }


    try
      {
	fittingPass(trk,direction);
      }
    catch(GFException &)
      {
	throw cet::exception("GFKalman.cxx: ") << " Line " << __LINE__ << ", " << __FILE__ << " ...Rethrow. \n";
      }
	
    //save first and last plane,state&cov after the fitting pass
    if(direction==1){//forward at last hit
      int nreps=trk->getNumReps();
      for(int i=0; i<nreps; ++i){
	trk->getTrackRep(i)->
	  setLastPlane( trk->getTrackRep(i)->getReferencePlane() );
	trk->getTrackRep(i)->
	  setLastState( trk->getTrackRep(i)->getState() );
	trk->getTrackRep(i)->
	  setLastCov( trk->getTrackRep(i)->getCov() );
      }
    }
    else{//backward at first hit
      int nreps=trk->getNumReps();
      for(int i=0; i<nreps; ++i){
	trk->getTrackRep(i)->
	  setFirstPlane( trk->getTrackRep(i)->getReferencePlane() );
	trk->getTrackRep(i)->
	  setFirstState( trk->getTrackRep(i)->getState() );
	trk->getTrackRep(i)->
	  setFirstCov( trk->getTrackRep(i)->getCov() );
      }
    }

    //switch direction of fitting and also inside all the reps
    if(direction==1) direction=-1;
    else direction=1;
    switchDirection(trk);

  }
  return;
}

void
genf::GFKalman::switchDirection(GFTrack* trk){
  int nreps=trk->getNumReps();
  for(int i=0; i<nreps; ++i){
    trk->getTrackRep(i)->switchDirection();
  }
}

void genf::GFKalman::blowUpCovsDiag(GFTrack* trk){
  int nreps=trk->getNumReps();
  for(int irep=0; irep<nreps; ++irep){
    GFAbsTrackRep* arep=trk->getTrackRep(irep);
    //dont do it for already compromsied reps, since they wont be fitted anyway
    if(arep->getStatusFlag()==0) { 
      TMatrixT<Double_t> cov = arep->getCov();
      for(int i=0;i<cov.GetNrows();++i){
	for(int j=0;j<cov.GetNcols();++j){
	  if(i!=j){//off diagonal
	    cov[i][j]=0.;
	  }
	  else{//diagonal
	    cov[i][j] = cov[i][j] * fBlowUpFactor;
	    //	    cov[0][0] = 0.1;
	  }
	}
      }
      arep->setCov(cov);
    }
  }  
}
void genf::GFKalman::blowUpCovs(GFTrack* trk){
  int nreps=trk->getNumReps();
  for(int irep=0; irep<nreps; ++irep){
    GFAbsTrackRep* arep=trk->getTrackRep(irep);
    //dont do it for already compromsied reps, since they wont be fitted anyway
    if(arep->getStatusFlag()==0) { 
      TMatrixT<Double_t> cov = arep->getCov();
      for(int i=0;i<cov.GetNrows();++i){
	for(int j=0;j<cov.GetNcols();++j){
	  //if(i!=j){//off diagonal
	  //  cov[i][j]=0.;
	  //}
	  //else{//diagonal
	  cov[i][j] = cov[i][j] * fBlowUpFactor;
	  //}
	}
      }
      arep->setCov(cov);
    }
  }  
}

void
genf::GFKalman::fittingPass(GFTrack* trk, int direction){
  //loop over hits

  unsigned int nhits=trk->getNumHits();
  unsigned int starthit=trk->getNextHitToFit();

  int nreps=trk->getNumReps();
  int ihit=(int)starthit;

  //clear chi2 sum and ndf sum in track reps
  for(int i=0;i<nreps;++i){
    GFAbsTrackRep* arep=trk->getTrackRep(i);
    arep->setChiSqu(0.);
    arep->setNDF(0);
  }

  //clear failedHits and outliers
  trk->getBK()->clearFailedHits();

  while((ihit<(int)nhits && direction==1) || (ihit>-1 && direction==-1)){
    //    GFAbsRecoHit* ahit=trk->getHit(ihit);
    // loop over reps
    for(int irep=0; irep<nreps; ++irep){
	  GFAbsTrackRep* arep=trk->getTrackRep(irep);
	  if(arep->getStatusFlag()==0) { 
	    try {
	      processHit(trk,ihit,irep,direction);
	    }
	    catch(GFException& e) {
	      trk->addFailedHit(irep,ihit);
	      std::cerr << e.what();
	      e.info();
	      if(e.isFatal()) {
		arep->setStatusFlag(1);
		continue; // go to next rep immediately
	      }
	    }	
	  }
    }// end loop over reps
    ihit+=direction;
  }// end loop over hits
  trk->setNextHitToFit(ihit-2*direction);
  //trk->printGFBookkeeping();
}

double genf::GFKalman::chi2Increment(const TMatrixT<Double_t>& r,const TMatrixT<Double_t>& H,
			     const TMatrixT<Double_t>& cov,const TMatrixT<Double_t>& V){

  // residuals covariances:R=(V - HCH^T)
  TMatrixT<Double_t> R(V);
  TMatrixT<Double_t> covsum1(cov,TMatrixT<Double_t>::kMultTranspose,H);
  TMatrixT<Double_t> covsum(H,TMatrixT<Double_t>::kMult,covsum1);

  R-=covsum;

  // chisq= r^TR^(-1)r
  double det=0.;
  TMatrixT<Double_t> Rsave(R);
  R.SetTol(1.0e-30); // to avoid inversion problem, EC, 8-Aug-2011. Was23, 9-Jan-2012.

  try
    {
      R.Invert(&det);
    }
  catch (cet::exception &)
    {
      GFException e("Kalman Chi2Increment: R is not even invertable. But keep plowing on ... ",__LINE__,__FILE__);
      //e.setFatal();
      //throw e;
    }
  if(TMath::IsNaN(det)) {
    GFException e("Kalman Chi2Increment: det of covsum is nan",__LINE__,__FILE__);
    e.setFatal();
    throw e;
  }
  TMatrixT<Double_t> residTranspose(r);
  residTranspose.T();
  TMatrixT<Double_t> chisq=residTranspose*(R*r);
  assert(chisq.GetNoElements()==1);

  if(TMath::IsNaN(chisq[0][0])){
	GFException exc("chi2 is nan",__LINE__,__FILE__);
	exc.setFatal();
	std::vector<double> numbers;
	numbers.push_back(det);
	exc.setNumbers("det",numbers);
	std::vector< TMatrixT<Double_t> > matrices;
	matrices.push_back(r);
	matrices.push_back(V);
	matrices.push_back(Rsave);
	matrices.push_back(R);
	matrices.push_back(cov);
	exc.setMatrices("r, V, Rsave, R, cov",matrices);
    throw exc;
  }

  return chisq[0][0];
}


double
genf::GFKalman::getChi2Hit(GFAbsRecoHit* hit, GFAbsTrackRep* rep)
{
  // get prototypes for matrices
  int repDim=rep->getDim();
  TMatrixT<Double_t> state(repDim,1);
  TMatrixT<Double_t> cov(repDim,repDim);;
  GFDetPlane pl=hit->getDetPlane(rep);
  rep->extrapolate(pl,state,cov);


  TMatrixT<Double_t> H = hit->getHMatrix(rep);
  // get hit covariances  
  TMatrixT<Double_t> V=hit->getHitCov(pl);
  
  TMatrixT<Double_t> r=hit->residualVector(rep,state,pl);
  assert(r.GetNrows()>0);

  r[0][0] = fabs(r[0][0]);
  //this is where chi2 is calculated
  double chi2 = chi2Increment(r,H,cov,V);

  return chi2/r.GetNrows();
}

  
void
genf::GFKalman::processHit(GFTrack* tr, int ihit, int irep,int direction){
  GFAbsRecoHit* hit = tr->getHit(ihit);
  GFAbsTrackRep* rep = tr->getTrackRep(irep);
  const GFTrackCand& cand = tr->getCand();

  // get prototypes for matrices
  int repDim = rep->getDim();
  TMatrixT<Double_t> state(repDim,1);
  TMatrixT<Double_t> cov(repDim,repDim);
  GFDetPlane pl, plPrev;
  unsigned int nhits=tr->getNumHits();
  int phit=ihit;
  static TMatrixT<Double_t> oldState(5,1);
  static std::vector<TVector3> pointsPrev;
  unsigned int indLkBk(3);

  if (direction>0 && ihit>0)
    {
      phit = ihit - 1;
    }
  if (direction<0 && ihit<((int)nhits-1))
    {
      phit = ihit + 1;
    }
  GFAbsRecoHit* hitPrev = tr->getHit(phit);

  /* do an extrapolation, if the trackrep irep is not given
   * at this ihit position. This will usually be the case, but
   * not if the fit turns around
   */
  //  std::cout << "GFKalman::ProcessHit(): ihit is " << ihit << std::endl;
  //  std::cout << "GFKalman::ProcessHit(): direction is "  << direction << std::endl;
  //rep->Print();


  Double_t thetaPlanes(0.0);
  if(ihit!=tr->getRepAtHit(irep)){
    //std::cout << "not same" << std::endl;
    // get the (virtual) detector plane. This call itself calls 
    // extrapolateToPoint(), which calls Extrap() with just 2 args and 
    // default 3rd arg, which propagates track through
    // material to find the next plane. But it in fact seems 
    // to just return this pl and plPrev, as desired. Something in Extrap()
    // kicks it out... even though I can see it walking through material ...
    // A mystery.

    pl=hit->getDetPlane(rep);
    plPrev=hitPrev->getDetPlane(rep);

      // std::cout << "GFKalman::ProcessHit(): hit is ... " << ihit << std::endl;
      //hit->Print();
      //std::cout << "GFKalman::ProcessHit(): plane is ... " <<  std::endl;
      //pl.Print();

    //do the extrapolation. This calls Extrap(), which wraps up RKutta and 
    // GFMaterials calls. state is intialized, pl is unchanged. Latter behavior
    // I will alter below.
    try{
      rep->extrapolate(pl,state,cov);
    }
    catch (cet::exception &)
      {
	throw cet::exception("GFKalman.cxx: ") << " Line " << __LINE__ << ", " << __FILE__ << " ...Rethrow. \n";
      }

  }
  else{
    //std::cout << "same" << std::endl;
    //    return;

    pl = rep->getReferencePlane();
    plPrev = hitPrev->getDetPlane(rep);
    state = rep->getState();
    cov = rep->getCov();
  }
  // Update plane. This is a code change from Genfit out of box.
  // state has accumulated the material/magnetic field changes since last
  // step. Here, we make the *plane* at this step know about that bend.
  TVector3 u(pl.getU());
  TVector3 v(pl.getV());
  TVector3 wold(u.Cross(v));
  Double_t sign(1.0);
  if ((direction==-1) && (ihit==nhits-1)) sign = -1.0;

  TVector3 pTilde = direction * (wold + state[1][0] * u + state[2][0] * v);
  TVector3 w(pTilde.Unit());
  // Find angle/vector through which we rotate. Use it subsequently.
  TVector3 rot(wold.Cross(w));
  Double_t ang(TMath::ACos(w*wold));
  ang = TMath::Min(ang,0.95*TMath::Pi()/2.0);
  /*
  if ((ihit==0 && direction==1) ||
      (ihit==nhits-1 && direction==-1) ||
      ihit!=phit
      ) // Don't do it for non-turnaround endpoints.
    {
      u.Rotate(ang,rot);
      v.Rotate(ang,rot);
      
      pl.setNormal(w);
      pl.setU(u.Unit());
      pl.setV(v.Unit());
    }
  */
  if(cov[0][0]<1.E-50){
        std::cout<<"GFKalman::processHit() 0. Calling Exception."<<std::endl;
	GFException exc(COVEXC,__LINE__,__FILE__);
	cov.Print();
        std::cout<<"GFKalman::processHit() 1. About to throw GFException."<<std::endl;
	throw exc;
        std::cout<<"GFKalman::processHit() 2. Back from GFException."<<std::endl;
  }
  
  /*
  if(direction==1){
    tr->getBK(irep)->setMatrix("fPredSt",ihit,state);
    tr->getBK(irep)->setMatrix("fPredCov",ihit,cov);
    tr->getBK(irep)->setGFDetPlane("f",ihit,pl);
  }
  else{
    tr->getBK(irep)->setMatrix("bPredSt",ihit,state);
    tr->getBK(irep)->setMatrix("bPredCov",ihit,cov);
    tr->getBK(irep)->setGFDetPlane("b",ihit,pl);
  }
  */

  //  TMatrixT<Double_t> origcov=rep->getCov();
  int pdg = -13; // cand above is not useful. Grr. How should I do this?
  // pdg = rep->getPDG();
  TParticlePDG * part = TDatabasePDG::Instance()->GetParticle(pdg);
  Double_t mass = part->Mass();

  Double_t dist = (pl.getO()-plPrev.getO()).Mag();
  //  Double_t mass = 0.104; // CHANGE ME!
  Double_t mom = fabs(1.0/state[0][0]);
  Double_t beta = mom/sqrt(mass*mass+mom*mom);
  const Double_t lowerLim(0.01);
  if (isnan(dist) || dist<=0.0) dist=lowerLim; // don't allow 0s here.
  if (isnan(beta) || beta<0.04) beta=0.04;
  TMatrixT<Double_t> H=hit->getHMatrix(rep,beta,dist);

  // get hit covariances  
  //  TMatrixT<Double_t> V=hit->getHitCov(pl);

// will force huge error on p, and thus insensitivity to p.
  TMatrixT<Double_t> V=hit->getHitCov(pl,plPrev,state,mass);
  tr->setMeasuredCov(ihit,V);
  // calculate kalman gain ------------------------------
  TMatrixT<Double_t> Gain(calcGain(cov,V,H));

  //std::cout << "GFKalman:: processHits(), state is  " << std::endl;
  //rep->getState().Print();

  TMatrixT<Double_t> res=hit->residualVector(rep,state,pl,plPrev);
  // calculate update -----------------------------------


  //  TVector3 pointer((pl.getO()-plPrev.getO()).Unit());

  TMatrixT<Double_t> rawcoord = hit->getRawHitCoord();
  TVector3 point(rawcoord[0][0],rawcoord[1][0],rawcoord[2][0]);
  TMatrixT<Double_t> prevrawcoord = hitPrev->getRawHitCoord();
  TVector3 pointPrev(prevrawcoord[0][0],prevrawcoord[1][0],prevrawcoord[2][0]);
  pointsPrev.push_back(pointPrev);
  TMatrixT<Double_t> Hnew(H);
  if (ihit==phit || ihit==0) pointsPrev.clear();
  /*
  if ((unsigned int)pointsPrev.size() >= indLkBk) 
    {
      pointPrev = pointsPrev[pointsPrev.size()-indLkBk];
      // recalculate the expected angular deflection
      Hnew = hit->getHMatrix(rep,beta,fabs((pointPrev-point).Mag()));
      Gain = calcGain(cov,V,Hnew); // need to do this with new Hnew
    }
  */
  TVector3 pointer((point-pointPrev).Unit());
  static TVector3 pointerPrev(pointer);
  if (ihit==0   ) {pointer[0] = 0.0;pointer[1] = 0.0;pointer[2] = 1.0;}
  if (ihit==phit) {pointer[0] = 0.0;pointer[1] = 0.0;pointer[2] = -1.0;}
  double thetaMeas = TMath::Min(fabs(pointer.Angle(pointerPrev)),0.95*TMath::Pi()/2.0);
  // Below line introduced because it's not true we predict the angle to be
  // precisely ang. If we'd taken this quantity from our transport result
  // (with measurement errors in it) we'd expect a smear like this on ang.
  // EC, 22-Feb-2012. 
  // Error is taken on th^2
  ang = (H*state)[0][0];
  // This 5 below makes some the difference. EC, 24-Mar-2012.
  thetaPlanes = fabs(ang*ang);// + 5.*fabs(gRandom->Gaus(0.0,V[0][0]));
  thetaPlanes = TMath::Min(sqrt(thetaPlanes),0.95*TMath::Pi()/2.0);

  Double_t dtheta =  thetaMeas - thetaPlanes; // was fabs(res[0][0]). EC, 26-Jan-2012
  if (ihit==phit || ihit==0) dtheta = 0.0;

  oldState = state;

  res[0][0] = dtheta;
  // The muon was propagated through material until its poca
  // to this current hit was found. At that point its plane is defined, 
  // in which du/dw=dv/dw=0 by construction. But, there's a deviation
  // in those two quantities measured in the *previous* plane which we want.
  TVector3 uPrev(plPrev.getU());
  TVector3 vPrev(plPrev.getV());
  TVector3 wPrev(u.Cross(v));
  // We want the newest momentum vector's du,v/dw defined in the plPrev plane. 
  // That is the predicted du,v/dw. We will subtract from the actual.
  res[1][0] = (pointer*uPrev)/(pointer*wPrev) - (w*uPrev)/(w*wPrev);
  res[2][0] = (pointer*vPrev)/(pointer*wPrev) - (w*vPrev)/(w*wPrev);

  TMatrixT<Double_t> update=Gain*res;
  //fUpdate->Fill((Double_t)(update[0][0]));
  //fIhitvUpdate->Fill((Double_t)ihit,(Double_t)(update[0][0]));

  // Don't let crazy ass spacepoints which pull the fit all over the place
  // be allowed to update the state. Uxe 10xRMS from fUpdate.
  // Don't allow fits above 15 GeV, or so. Else, 1/p will likely 
  // migrate across zero. Similarly, don't allow tiny momentum.

  if (fabs(update[0][0])>1.0e-1 ) update[0][0] = 0.0;
  state+=update; // prediction overwritten!      
  cov-=Gain*(Hnew*cov);

  if (fabs(1.0/state[0][0])<fMomLow) state[0][0] = 1.0/fMomLow*fabs(state[0][0])/state[0][0];
  if (fabs(1.0/state[0][0])>fMomHigh) state[0][0] = 1.0/fMomHigh*fabs(state[0][0])/state[0][0];


  // Let's also calculate the "filtered" plane, and pointer here.
  // Use those to calculate filtered chisq.
  GFDetPlane plFilt(pl);
  TVector3 uf(pl.getU());
  TVector3 vf(pl.getV());
  TVector3 wf(uf.Cross(vf));
  TVector3 Of(pl.getO());
  // direction *
  Double_t sign2(1.0);
  if (direction==-1 && ihit==nhits-1) sign2 = -1.0;

  TVector3 pf = direction*(wf + state[1][0] * uf + state[2][0] * vf);
  TVector3 pposf = Of + state[3][0] * uf + state[4][0] * vf;
  Double_t angf = TMath::Min(fabs(pf.Angle(wf)),0.95*TMath::Pi()/2.0);
  TVector3 rotf(wf.Cross(pf.Unit()));
  
  if ((ihit==0 && direction==1) ||
      (ihit==nhits-1 && direction==-1) ||
      ihit!=phit
      ) // Don't do it for non-turnaround endpoints.
    {
      uf.Rotate(angf,rotf);
      vf.Rotate(angf,rotf);
      wf = uf.Cross(vf);
      plFilt.setU(uf.Unit());
      plFilt.setV(vf.Unit());
      //plFilt.setO(pposf);
      plFilt.setNormal(pf.Unit());
    }
  
  // calculate filtered chisq from filtered residuals
  TMatrixT<Double_t> r=hit->residualVector(rep,state,plFilt,plPrev);
  dtheta = thetaMeas - TMath::Min(fabs(wold.Angle(plFilt.getNormal())),0.95*TMath::Pi()/2.);
  r[0][0] = dtheta;
  r[1][0] = (pointer*uPrev)/(pointer*wPrev) - (wf*uPrev)/(wf*wPrev);
  r[2][0] = (pointer*vPrev)/(pointer*wPrev) - (wf*vPrev)/(wf*wPrev);

  double chi2 = chi2Increment(r,Hnew,cov,V);
  int ndf = r.GetNrows();
  rep->addChiSqu( chi2 );
  rep->addNDF( ndf );
  pointerPrev = pointer;

  /*
  if(direction==1){
    tr->getBK(irep)->setNumber("fChi2",ihit,chi2/ndf);
  }
  else{
    tr->getBK(irep)->setNumber("bChi2",ihit,chi2/ndf);
  }
  */

  // if we survive until here: update TrackRep
  //rep->setState(state);
  //rep->setCov(cov);
  //rep->setReferencePlane(pl);

  // Since I've tilted the planes, d(u,v)/dw=0 by construction.
  state[1][0] = 0.0;
  state[2][0] = 0.0;
  rep->setData(state,pl,&cov); // This is where fState,
                               // fRefPlane and fCov are updated.
  tr->setRepAtHit(irep,ihit);

}


TMatrixT<Double_t>
genf::GFKalman::calcGain(const TMatrixT<Double_t>& cov, 
			 const TMatrixT<Double_t>& HitCov,
			 const TMatrixT<Double_t>& H){

  // calculate covsum (V + HCH^T)

  // Comment next 3 out for normal running.
  //cov.Print();
  //H.Print();
  //HitCov.Print();
  TMatrixT<Double_t> covsum1(cov,TMatrixT<Double_t>::kMultTranspose,H);
  //  covsum1.Print();
  TMatrixT<Double_t> covsum(H,TMatrixT<Double_t>::kMult,covsum1);

  //covsum.Print();

  covsum+=HitCov;
  //covsum = (covsum,TMatrixT<Double_t>::kPlus,HitCov);
  // covsum.Print();

  // invert
  double det=0;
  covsum.SetTol(1.0e-23); // to avoid inversion problem, EC, 8-Aug-2011.
  covsum.Invert(&det);
  //    std::cout << "GFKalman:: calGain(), det is  "<< det << std::endl;
  if(TMath::IsNaN(det)) {
    GFException e("Kalman Gain: det of covsum is nan",__LINE__,__FILE__);
    e.setFatal();
    throw e;
  }

  if(det==0){
	GFException exc("cannot invert covsum in Kalman Gain - det=0",
						__LINE__,__FILE__);
	exc.setFatal();
	std::vector< TMatrixT<Double_t> > matrices;
	matrices.push_back(cov);
	matrices.push_back(HitCov);
	matrices.push_back(covsum1);
	matrices.push_back(covsum);
	exc.setMatrices("cov, HitCov, covsum1 and covsum",matrices);
    throw exc;

  }

  // calculate gain
  TMatrixT<Double_t> gain1(H,TMatrixT<Double_t>::kTransposeMult,covsum);
  TMatrixT<Double_t> gain(cov,TMatrixT<Double_t>::kMult,gain1);

  return gain;
}

