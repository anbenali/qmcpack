//////////////////////////////////////////////////////////////////
// (c) Copyright 2006- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCWaveFunctions/AGPDeterminant.h"
#include "Numerics/DeterminantOperators.h"
#include "Numerics/OhmmsBlas.h"
#include "Numerics/MatrixOperators.h"

namespace qmcplusplus {

  AGPDeterminant::AGPDeterminant(BasisSetType* bs): GeminalBasis(bs), NumPtcls(0){
  }
  AGPDeterminant::~AGPDeterminant() {}

  void AGPDeterminant::resize(int nup, int ndown) {
    BasisSize=GeminalBasis->size();

    if(NumPtcls == 0) { //use Numptcls to ensure resizing only once
      Lambda.resize(BasisSize,BasisSize);
      if(nup>ndown) {
        LambdaUP.resize(nup-ndown,BasisSize);
      }

      Nup=nup;
      Ndown=ndown;
      NumPtcls=nup+ndown;

      psiM.resize(nup,nup);
      psiM_temp.resize(nup,nup);
      psiU.resize(nup);
      psiD.resize(ndown);
      phiT.resize(NumPtcls,BasisSize);
      phiTv.resize(BasisSize);

      WorkSpace.resize(nup);
      Pivot.resize(nup);
    }

    app_log() << "  AGPDetermiant::resize checking size nup, ndown, basis " << nup << " " << ndown 
      << " " << BasisSize << endl;
  }

  void AGPDeterminant::resetTargetParticleSet(ParticleSet& P) { 
    GeminalBasis->resetTargetParticleSet(P);
  }

  /** Calculate the value of the Dirac determinant for particles
   *@param P input configuration containing N particles
   *@param G a vector containing N gradients
   *@param L a vector containing N laplacians
   *@return the value of the determinant
   *
   *\f$ (first,first+nel). \f$  Add the gradient and laplacian 
   *contribution of the determinant to G(radient) and L(aplacian)
   *for local energy calculations.
   */ 
  AGPDeterminant::ValueType
  AGPDeterminant::evaluate(ParticleSet& P, 
      ParticleSet::ParticleGradient_t& G, 
      ParticleSet::ParticleLaplacian_t& L){

    GeminalBasis->evaluate(P);

    /* evaluate psi_up(iat)= \sum_{j} C_{ij} \phi_j^{u}(r_{iat}) 
     * psi_down(iat-Nup) =  \sum_{j} C_{ij} \phi_j^{d}(r_{iat})
     */
    MatrixOperators::product(GeminalBasis->Y, Lambda, phiT);

    //psiM=0.0;
    for(int u=0; u<Nup; u++) {
      //paired block
      for(int d=0, jat=Nup; d<Ndown; d++,jat++) {
        psiM(d,u) = BLAS::dot(BasisSize,phiT[u],GeminalBasis->y(jat));
      }
      //unpaired block Ndown x unpaired
      for(int d=Ndown,unpaired=0; d<Nup; d++,unpaired++) {
        psiM(d,u) = BLAS::dot(BasisSize,LambdaUP[unpaired],GeminalBasis->y(u));
      }
    }

    CurrentDet = Invert(psiM.data(),Nup,Nup,WorkSpace.data(),Pivot.data());

    for(int iat=0; iat<Nup; iat++) {
      GradType rv;
      ValueType lap=0;
      int jat=Nup;
      for(int d=0; d<Ndown; d++,jat++) {
        ValueType dfac=psiM(iat,d);
        rv += dfac*dot(phiT[jat],GeminalBasis->dy(iat),BasisSize);
        lap += dfac*dot(phiT[jat],GeminalBasis->d2y(iat),BasisSize);
      }
      for(int d=Ndown,unpaired=0; d<Nup; d++,unpaired++) {
        ValueType dfac=psiM(iat,d);
        rv += dfac*dot(LambdaUP[unpaired],GeminalBasis->dy(iat),BasisSize);
        lap += dfac*dot(LambdaUP[unpaired],GeminalBasis->d2y(iat),BasisSize);
      }
      G(iat) += rv;
      L(iat) += lap-dot(rv,rv);
    }

    for(int jat=Nup,d=0; jat<NumPtcls; jat++,d++) {
      GradType rv;
      ValueType lap=0;
      for(int u=0; u<Nup; u++) {
        ValueType dfac=psiM(u,d);
        rv += dfac*dot(phiT[u],GeminalBasis->dy(jat),BasisSize);
        lap += dfac*dot(phiT[u],GeminalBasis->d2y(jat),BasisSize);
      }
      G(jat) += rv;
      L(jat) += lap-dot(rv,rv);
    }

    return CurrentDet;
  }

  void
  AGPDeterminant::evaluateLogAndStore(ParticleSet& P) {

    GeminalBasis->evaluate(P);

    /* evaluate psi_up(iat)= \sum_{j} C_{ij} \phi_j^{u}(r_{iat}) 
     * psi_down(iat-Nup) =  \sum_{j} C_{ij} \phi_j^{d}(r_{iat})
     */
    MatrixOperators::product(GeminalBasis->Y, Lambda, phiT);

    for(int u=0; u<Nup; u++) {
      //paired block
      for(int d=0, jat=Nup; d<Ndown; d++,jat++) {
        psiM(d,u) = BLAS::dot(BasisSize,phiT[u],GeminalBasis->y(jat));
      }
      //unpaired block Ndown x unpaired
      for(int d=Ndown,unpaired=0; d<Nup; d++,unpaired++) {
        psiM(d,u) = BLAS::dot(BasisSize,LambdaUP[unpaired],GeminalBasis->y(u));
      }
    }

    CurrentDet = Invert(psiM.data(),Nup,Nup,WorkSpace.data(),Pivot.data());

    for(int iat=0; iat<Nup; iat++) {
      for(int d=0,jat=Nup; d<Ndown; d++,jat++) {
        dpsiU(iat,d)=dot(phiT[jat],GeminalBasis->dy(iat),BasisSize);
        d2psiU(iat,d)=dot(phiT[jat],GeminalBasis->d2y(iat),BasisSize);
      }
      for(int d=Ndown,unpaired=0; d<Nup; d++,unpaired++) {
        dpsiU(iat,d)=dot(LambdaUP[unpaired],GeminalBasis->dy(iat),BasisSize);
        d2psiU(iat,d)=dot(LambdaUP[unpaired],GeminalBasis->d2y(iat),BasisSize);
      }
      GradType rv=dot(psiM[iat],dpsiU[iat],Nup);
      ValueType lap=dot(psiM[iat],d2psiU[iat],Nup);
      myG[iat]=rv;
      myL[iat]=lap-dot(rv,rv);
    }

    for(int jat=Nup,d=0; jat<NumPtcls; jat++,d++) {
      GradType rv;
      ValueType lap=0;
      for(int u=0; u<Nup; u++) {
        ValueType dfac=psiM(u,d);
        rv += dfac*(dpsiD(d,u)=dot(phiT[u],GeminalBasis->dy(jat),BasisSize));
        lap += dfac*(d2psiD(d,u)=dot(phiT[u],GeminalBasis->d2y(jat),BasisSize));
      }
      myG[jat]=rv;
      myL[jat]=lap-dot(rv,rv);
    }

    dY = GeminalBasis->dY;
    d2Y = GeminalBasis->d2Y;
  }

  AGPDeterminant::ValueType 
  AGPDeterminant::registerData(ParticleSet& P, PooledData<RealType>& buf) {

    if(myG.size() == 0) {
      myG.resize(NumPtcls);
      myL.resize(NumPtcls);
      myG_temp.resize(NumPtcls);
      myL_temp.resize(NumPtcls);
      dY.resize(NumPtcls,BasisSize);
      d2Y.resize(NumPtcls,BasisSize);
      dpsiU.resize(Nup,Nup);
      dpsiD.resize(Ndown,Nup);
      d2psiU.resize(Nup,Nup);
      d2psiD.resize(Ndown,Nup);
      dpsiUv.resize(Nup);
      d2psiUv.resize(Nup);
      dpsiDv.resize(Nup);
      d2psiDv.resize(Nup);

      FirstAddressOfdVU=&(dpsiU(0,0)[0]);
      LastAddressOfdVU=FirstAddressOfdVU+DIM*Nup*Nup;
      FirstAddressOfdVD=&(dpsiD(0,0)[0]);
      LastAddressOfdVD=FirstAddressOfdVD+DIM*Ndown*Nup;
      FirstAddressOfdY=&(dY(0,0)[0]);
      LastAddressOfdY=FirstAddressOfdY+DIM*NumPtcls*BasisSize;
      FirstAddressOfG = &myG[0][0];
      LastAddressOfG = FirstAddressOfG + DIM*NumPtcls;
    }

    evaluateLogAndStore(P);

    P.G += myG;
    P.L += myL;

    //copy psiM to temporary
    psiM_temp = psiM;

    if(UseBuffer) {  //add the data: determinant, inverse, gradient and laplacians
      buf.add(CurrentDet);
      buf.add(psiM.begin(),psiM.end());
      buf.add(phiT.begin(),phiT.end());
      buf.add(d2psiU.begin(),d2psiU.end());
      buf.add(d2psiD.begin(),d2psiD.end());
      buf.add(FirstAddressOfdVU,LastAddressOfdVU);
      buf.add(FirstAddressOfdVD,LastAddressOfdVD);
      buf.add(d2Y.begin(),d2Y.end());
      buf.add(FirstAddressOfdY,LastAddressOfdY);
      buf.add(FirstAddressOfG,LastAddressOfG);
      buf.add(myL.begin(), myL.end());
    }

    SignValue = (CurrentDet<0.0)?-1.0:1.0;
    return LogValue = std::log(abs(CurrentDet));
  }

  AGPDeterminant::ValueType 
  AGPDeterminant::updateBuffer(ParticleSet& P, PooledData<RealType>& buf) {

    evaluateLogAndStore(P);

    P.G += myG;
    P.L += myL;

    if(UseBuffer) {
      buf.put(CurrentDet);
      buf.put(psiM.begin(),psiM.end());
      buf.put(phiT.begin(),phiT.end());
      buf.put(d2psiU.begin(),d2psiU.end());
      buf.put(d2psiD.begin(),d2psiD.end());
      buf.put(FirstAddressOfdVU,LastAddressOfdVU);
      buf.put(FirstAddressOfdVD,LastAddressOfdVD);
      buf.put(d2Y.begin(),d2Y.end());
      buf.put(FirstAddressOfdY,LastAddressOfdY);
      buf.put(FirstAddressOfG,LastAddressOfG);
      buf.put(myL.begin(), myL.end());
    }

    return CurrentDet;
  }

  void AGPDeterminant::copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf) {
    if(UseBuffer) {
      buf.get(CurrentDet);
      buf.get(psiM.begin(),psiM.end());
      buf.get(phiT.begin(),phiT.end());
      buf.get(d2psiU.begin(),d2psiU.end());
      buf.get(d2psiD.begin(),d2psiD.end());
      buf.get(FirstAddressOfdVU,LastAddressOfdVU);
      buf.get(FirstAddressOfdVD,LastAddressOfdVD);
      buf.get(d2Y.begin(),d2Y.end());
      buf.get(FirstAddressOfdY,LastAddressOfdY);
      buf.get(FirstAddressOfG,LastAddressOfG);
      buf.get(myL.begin(), myL.end());
      //copy current inverse of the determinant
      psiM_temp = psiM;
    }
  }


    /** return the ratio only for the  iat-th partcle move
     * @param P current configuration
     * @param iat the particle thas is being moved
     */
    AGPDeterminant::ValueType 
    AGPDeterminant::ratio(ParticleSet& P, int iat) {

      std::copy(phiT[iat],phiT[iat]+BasisSize,phiTv.begin());
      GeminalBasis->evaluate(P,iat);
      //BLAS::gemv(Lambda.rows(),Lambda.cols(), Lambda.data(), GeminalBasis->y(0), phiT[iat]);

      const ValueType* restrict y_ptr=GeminalBasis->y(0);
      if(iat<Nup) {
        for(int d=0,jat=Nup; d<Ndown; d++,jat++) {
          psiU[d]=BLAS::dot(BasisSize,y_ptr,phiT[jat]);
        }
        //unpaired block Ndown x unpaired
        for(int d=Ndown,unpaired=0; d<Nup; d++,unpaired++) {
          psiU[d] = BLAS::dot(BasisSize,LambdaUP[unpaired],y_ptr);
        }
        return DetRatio(psiM, psiU.data(),iat);
      } else {
        for(int u=0; u<Nup; u++) {
          psiD[u]=BLAS::dot(BasisSize,y_ptr,phiT[u]);
        }
        return DetRatioTranspose(psiM, psiD.data(),iat-Nup);
      }
    }

    /** return the ratio
     * @param P current configuration
     * @param iat particle whose position is moved
     * @param dG differential Gradients
     * @param dL differential Laplacians
     *
     * Data member *_temp contain the data assuming that the move is accepted
     * and are used to evaluate differential Gradients and Laplacians.
     */
    AGPDeterminant::ValueType 
    AGPDeterminant::ratio(ParticleSet& P, int iat,
		    ParticleSet::ParticleGradient_t& dG, 
		    ParticleSet::ParticleLaplacian_t& dL) {

      //copy the iat-row to temporary vectors, restore when rejected
      std::copy(phiT[iat],phiT[iat]+BasisSize,phiTv.begin());

      GeminalBasis->evaluateAll(P,iat);

      BLAS::gemv(Lambda.rows(),Lambda.cols(), Lambda.data(), GeminalBasis->y(0), phiT[iat]);

      if(iat<Nup)  
        ratioUp(P,iat);
      else
        ratioDown(P,iat);

      for(int kat=0; kat<Nup; kat++) {
        GradType rv=dot(psiM_temp[kat],dpsiU[kat],Nup);
        ValueType lap=dot(psiM_temp[kat],d2psiU[kat],Nup);
        lap -= dot(rv,rv);
        dG[kat] += (rv-myG[kat]); myG_temp[kat]=rv;
        dL[kat] += (lap-myL[kat]); myL_temp[kat]=lap;
      }

      for(int jat=Nup,d=0; jat<NumPtcls; jat++,d++) {
        GradType rv;
        ValueType lap=0;
        for(int u=0; u<Nup; u++) {
          ValueType dfac=psiM_temp(u,d);
          rv += dfac*dpsiD(d,u);
          lap += dfac*d2psiD(d,u);
        }
        lap -= dot(rv,rv);
        dG[jat] +=  (rv-myG[jat]); myG_temp[jat]=rv;
        dL[jat] +=  (lap-myL[jat]); myL_temp[jat]=lap;
      }

      return curRatio;
    }

    void
    AGPDeterminant::ratioUp(ParticleSet& P, int iat) {
      const ValueType* restrict y_ptr=GeminalBasis->y(0);
      for(int d=0,jat=Nup; d<Ndown; d++,jat++) {
        psiU[d]=BLAS::dot(BasisSize,y_ptr,phiT[jat]);
      }
      //unpaired block Ndown x unpaired
      for(int d=Ndown,unpaired=0; d<Nup; d++,unpaired++) {
        psiU[d] = BLAS::dot(BasisSize,LambdaUP[unpaired],y_ptr);
      }

      curRatio = DetRatio(psiM_temp, psiU.data(),iat);
      DetUpdate(psiM_temp,psiU,workV1,workV2,iat,curRatio);

      std::copy(dpsiU[iat],dpsiU[iat]+Nup,dpsiUv.begin());
      std::copy(d2psiU[iat],d2psiU[iat]+Nup,d2psiUv.begin());

      const GradType* restrict  dy_ptr = GeminalBasis->dy(0);
      const ValueType* restrict d2y_ptr = GeminalBasis->d2y(0);
      for(int d=0, jat=Nup; d<Ndown; d++,jat++) {
        dpsiU(iat,d)=dot(phiT[jat],dy_ptr,BasisSize);
        d2psiU(iat,d)=dot(phiT[jat],d2y_ptr,BasisSize);
      }
      for(int d=Ndown,unpaired=0; d<Nup; d++,unpaired++) {
        dpsiU(iat,d)=dot(LambdaUP[unpaired],dy_ptr,BasisSize);
        d2psiU(iat,d)=dot(LambdaUP[unpaired],d2y_ptr,BasisSize);
      }

      for(int jat=Nup,d=0; jat<NumPtcls; jat++,d++) {
        dpsiDv[d]=dpsiD(d,iat);
        dpsiD(d,iat)=dot(phiT[iat],dY[jat],BasisSize);

        d2psiDv[d]=d2psiD(d,iat);
        d2psiD(d,iat)=dot(phiT[iat],d2Y[jat],BasisSize);
      }
    }

    void
    AGPDeterminant::ratioDown(ParticleSet& P, int iat) {
      const ValueType* restrict y_ptr=GeminalBasis->y(0);
      int d=iat-Nup;
      for(int u=0; u<Nup; u++) {
        psiD[u]=BLAS::dot(BasisSize,y_ptr,phiT[u]);
      }

      curRatio = DetRatioTranspose(psiM_temp, psiD.data(),d);
      DetUpdateTranspose(psiM_temp,psiD,workV1,workV2,d,curRatio);

      std::copy(dpsiD[d],dpsiD[d]+Nup,dpsiDv.begin());
      std::copy(d2psiD[d],d2psiD[d]+Nup,d2psiDv.begin());

      const GradType* restrict dy_ptr = GeminalBasis->dy(0);
      const ValueType* restrict d2y_ptr = GeminalBasis->d2y(0);
      for(int u=0; u<Nup; u++) {
        dpsiD(d,u)=dot(phiT[u],dy_ptr,BasisSize);
        d2psiD(d,u)=dot(phiT[u],d2y_ptr,BasisSize);
      }

      for(int kat=0; kat<Nup; kat++) {
        dpsiUv[kat]=dpsiU(kat,d);
        dpsiU(kat,d)=dot(phiT[iat],dY[kat],BasisSize);

        d2psiUv[kat]=d2psiU(kat,d);
        d2psiU(kat,d)=dot(phiT[iat],d2Y[kat],BasisSize);
      }
    }


    /** move was accepted, update the real container
     */
    void AGPDeterminant::acceptMove(ParticleSet& P, int iat) {
      CurrentDet *= curRatio;
      myG = myG_temp;
      myL = myL_temp;
      psiM = psiM_temp;
      std::copy(GeminalBasis->dy(0),GeminalBasis->dy(0)+BasisSize,dY[iat]);
      std::copy(GeminalBasis->d2y(0),GeminalBasis->d2y(0)+BasisSize,d2Y[iat]);
      curRatio=1.0;
    }

    /** move was rejected. copy the real container to the temporary to move on
     */
    void AGPDeterminant::restore(int iat) {
      psiM_temp = psiM;
      std::copy(phiTv.begin(), phiTv.end(),phiT[iat]);
      if(iat<Nup) {
        std::copy(dpsiUv.begin(), dpsiUv.end(),dpsiU[iat]);
        std::copy(d2psiUv.begin(), d2psiUv.end(),d2psiU[iat]);
        for(int d=0; d<Ndown; d++) {
          dpsiD(d,iat)=dpsiDv[d];
          d2psiD(d,iat)=d2psiDv[d];
        }
      } else {
        int d=iat-Nup;
        std::copy(dpsiDv.begin(), dpsiDv.end(),dpsiD[d]);
        std::copy(d2psiDv.begin(), d2psiDv.end(),d2psiD[d]);
        for(int kat=0; kat<Nup; kat++) {
          dpsiU(kat,d)=dpsiUv[kat];
          d2psiU(kat,d) = d2psiUv[kat];
        }
      }
      curRatio=1.0;
    }
    
    void AGPDeterminant::update(ParticleSet& P, 
		ParticleSet::ParticleGradient_t& dG, 
		ParticleSet::ParticleLaplacian_t& dL,
		int iat) {
      cout << "What is going on here?" << endl;
    }

  AGPDeterminant::ValueType 
  AGPDeterminant::evaluate(ParticleSet& P, PooledData<RealType>& buf) {
    if(UseBuffer) {
      buf.put(CurrentDet);
      buf.put(psiM.begin(),psiM.end());
      buf.put(phiT.begin(),phiT.end());
      buf.put(d2psiU.begin(),d2psiU.end());
      buf.put(d2psiD.begin(),d2psiD.end());
      buf.put(FirstAddressOfdVU,LastAddressOfdVU);
      buf.put(FirstAddressOfdVD,LastAddressOfdVD);
      buf.put(d2Y.begin(),d2Y.end());
      buf.put(FirstAddressOfdY,LastAddressOfdY);
      buf.put(FirstAddressOfG,LastAddressOfG);
      buf.put(myL.begin(), myL.end());
    }
    return CurrentDet;
  }

  void AGPDeterminant::resizeByWalkers(int nwalkers) {
    //don't know what to do here
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
