//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "SlaterDet.h"
#include "Message/Communicate.h"

namespace qmcplusplus
{
// for return types
using PsiValueType = WaveFunctionComponent::PsiValueType;

SlaterDet::SlaterDet(ParticleSet& targetPtcl, const std::string& class_name) : WaveFunctionComponent(class_name)
{
  Optimizable  = false;
  is_fermionic = true;

  Last.resize(targetPtcl.groups());
  for (int i = 0; i < Last.size(); ++i)
    Last[i] = targetPtcl.last(i) - 1;

  Dets.resize(targetPtcl.groups());
}

///destructor
SlaterDet::~SlaterDet()
{
  ///clean up SPOSet
}

///add a new DiracDeterminant to the list of determinants
void SlaterDet::add(Determinant_t* det, int ispin)
{
  if (Dets[ispin] != nullptr)
  {
    APP_ABORT("SlaterDet::add(Determinant_t* det, int ispin) is alreaded instantiated.");
  }
  else
    Dets[ispin].reset(det);
  Optimizable = Optimizable || det->Optimizable;
}

void SlaterDet::checkInVariables(opt_variables_type& active)
{
  myVars.clear();
  if (Optimizable)
    for (int i = 0; i < Dets.size(); i++)
    {
      Dets[i]->checkInVariables(active);
      Dets[i]->checkInVariables(myVars);
    }
}

void SlaterDet::checkOutVariables(const opt_variables_type& active)
{
  myVars.clear();
  if (Optimizable)
    for (int i = 0; i < Dets.size(); i++)
    {
      Dets[i]->checkOutVariables(active);
      myVars.insertFrom(Dets[i]->myVars);
    }
  myVars.getIndex(active);
}

///reset all the Dirac determinants, Optimizable is true
void SlaterDet::resetParameters(const opt_variables_type& active)
{
  if (Optimizable)
    for (int i = 0; i < Dets.size(); i++)
      Dets[i]->resetParameters(active);
}

void SlaterDet::reportStatus(std::ostream& os) {}

PsiValueType SlaterDet::ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
{
  return Dets[getDetID(iat)]->ratioGrad(P, iat, grad_iat);
}
void SlaterDet::ratioGradAsync(ParticleSet& P, int iat, PsiValueType& ratio, GradType& grad_iat)
{
  Dets[getDetID(iat)]->ratioGradAsync(P, iat, ratio, grad_iat);
}

PsiValueType SlaterDet::ratioGradWithSpin(ParticleSet& P, int iat, GradType& grad_iat, ComplexType& spingrad_iat)
{
  return Dets[getDetID(iat)]->ratioGradWithSpin(P, iat, grad_iat, spingrad_iat);
}

void SlaterDet::mw_ratioGrad(const RefVector<WaveFunctionComponent>& wfc_list,
                             const RefVector<ParticleSet>& P_list,
                             int iat,
                             std::vector<PsiValueType>& ratios,
                             std::vector<GradType>& grad_now)
{
  const int det_id = getDetID(iat);
  Dets[det_id]->mw_ratioGrad(extract_DetRef_list(wfc_list, det_id), P_list, iat, ratios, grad_now);
}

void SlaterDet::mw_ratioGradAsync(const RefVector<WaveFunctionComponent>& wfc_list,
                                  const RefVector<ParticleSet>& P_list,
                                  int iat,
                                  std::vector<PsiValueType>& ratios,
                                  std::vector<GradType>& grad_now)
{
  const int det_id = getDetID(iat);
  Dets[det_id]->mw_ratioGradAsync(extract_DetRef_list(wfc_list, det_id), P_list, iat, ratios, grad_now);
}

void SlaterDet::evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios)
{
  for (int i = 0; i < Dets.size(); ++i)
    Dets[i]->evaluateRatiosAlltoOne(P, ratios);
}

SlaterDet::LogValueType SlaterDet::evaluateLog(ParticleSet& P,
                                               ParticleSet::ParticleGradient_t& G,
                                               ParticleSet::ParticleLaplacian_t& L)
{
  LogValue = 0.0;
  for (int i = 0; i < Dets.size(); ++i)
    LogValue += Dets[i]->evaluateLog(P, G, L);
  return LogValue;
}

void SlaterDet::mw_evaluateLog(const RefVector<WaveFunctionComponent>& WFC_list,
                               const RefVector<ParticleSet>& P_list,
                               const RefVector<ParticleSet::ParticleGradient_t>& G_list,
                               const RefVector<ParticleSet::ParticleLaplacian_t>& L_list)
{
  constexpr LogValueType czero(0);

  for (WaveFunctionComponent& wfc : WFC_list)
    wfc.LogValue = czero;

  for (int i = 0; i < Dets.size(); ++i)
  {
    const auto Det_list(extract_DetRef_list(WFC_list, i));
    Dets[i]->mw_evaluateLog(Det_list, P_list, G_list, L_list);
    for (int iw = 0; iw < WFC_list.size(); iw++)
      WFC_list[iw].get().LogValue += Det_list[iw].get().LogValue;
  }
}

SlaterDet::LogValueType SlaterDet::evaluateGL(ParticleSet& P,
                                              ParticleSet::ParticleGradient_t& G,
                                              ParticleSet::ParticleLaplacian_t& L,
                                              bool from_scratch)
{
  LogValue = 0.0;
  for (int i = 0; i < Dets.size(); ++i)
    LogValue += Dets[i]->evaluateGL(P, G, L, from_scratch);
  return LogValue;
}

void SlaterDet::mw_evaluateGL(const RefVector<WaveFunctionComponent>& WFC_list,
                              const RefVector<ParticleSet>& P_list,
                              const RefVector<ParticleSet::ParticleGradient_t>& G_list,
                              const RefVector<ParticleSet::ParticleLaplacian_t>& L_list,
                              bool fromscratch)
{
  constexpr LogValueType czero(0);

  for (WaveFunctionComponent& wfc : WFC_list)
    wfc.LogValue = czero;

  for (int i = 0; i < Dets.size(); ++i)
  {
    const auto Det_list(extract_DetRef_list(WFC_list, i));
    Dets[i]->mw_evaluateGL(Det_list, P_list, G_list, L_list, fromscratch);
    for (int iw = 0; iw < WFC_list.size(); iw++)
      WFC_list[iw].get().LogValue += Det_list[iw].get().LogValue;
  }
}

void SlaterDet::recompute(ParticleSet& P)
{
  for (int i = 0; i < Dets.size(); ++i)
    Dets[i]->recompute(P);
}

void SlaterDet::evaluateHessian(ParticleSet& P, HessVector_t& grad_grad_psi)
{
  grad_grad_psi.resize(P.getTotalNum());
  HessVector_t tmp;
  tmp.resize(P.getTotalNum());
  for (int i = 0; i < Dets.size(); ++i)
  {
    tmp = 0;
    Dets[i]->evaluateHessian(P, tmp);
    //  app_log()<<"squee ----- "<<i<< std::endl;
    //  app_log()<<"grad_grad_psi = "<<grad_grad_psi<< std::endl;
    //  app_log()<<"tmp = "<<tmp<< std::endl<< std::endl;
    grad_grad_psi += tmp;
  }
}

void SlaterDet::createResource(ResourceCollection& collection)
{
  for (int i = 0; i < Dets.size(); ++i)
    Dets[i]->createResource(collection);
}

void SlaterDet::acquireResource(ResourceCollection& collection)
{
  for (int i = 0; i < Dets.size(); ++i)
    Dets[i]->acquireResource(collection);
}

void SlaterDet::releaseResource(ResourceCollection& collection)
{
  for (int i = 0; i < Dets.size(); ++i)
    Dets[i]->releaseResource(collection);
}

void SlaterDet::registerData(ParticleSet& P, WFBufferType& buf)
{
  DEBUG_PSIBUFFER(" SlaterDet::registerData ", buf.current());
  for (int i = 0; i < Dets.size(); ++i)
    Dets[i]->registerData(P, buf);
  DEBUG_PSIBUFFER(" SlaterDet::registerData ", buf.current());
}

SlaterDet::LogValueType SlaterDet::updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch)
{
  DEBUG_PSIBUFFER(" SlaterDet::updateBuffer ", buf.current());
  LogValue = 0.0;
  for (int i = 0; i < Dets.size(); ++i)
    LogValue += Dets[i]->updateBuffer(P, buf, fromscratch);
  DEBUG_PSIBUFFER(" SlaterDet::updateBuffer ", buf.current());
  return LogValue;
}

void SlaterDet::copyFromBuffer(ParticleSet& P, WFBufferType& buf)
{
  DEBUG_PSIBUFFER(" SlaterDet::copyFromBuffer ", buf.current());
  for (int i = 0; i < Dets.size(); i++)
    Dets[i]->copyFromBuffer(P, buf);
  DEBUG_PSIBUFFER(" SlaterDet::copyFromBuffer ", buf.current());
}

WaveFunctionComponentPtr SlaterDet::makeClone(ParticleSet& tqp) const
{
  SlaterDet* myclone   = new SlaterDet(tqp);
  myclone->Optimizable = Optimizable;
  for (int i = 0; i < Dets.size(); ++i)
  {
    Determinant_t* newD = Dets[i]->makeCopy(std::unique_ptr<SPOSet>(Dets[i]->getPhi()->makeClone()));
    myclone->add(newD, i);
  }
  return myclone;
}

} // namespace qmcplusplus
