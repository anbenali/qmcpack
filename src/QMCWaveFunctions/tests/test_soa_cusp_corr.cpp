
//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "Configuration.h"

#include "Message/Communicate.h"
#include "Numerics/OneDimGridBase.h"
#include "Particle/DistanceTableData.h"
#include "ParticleIO/XMLParticleIO.h"
#include "Numerics/GaussianBasisSet.h"

#include "QMCWaveFunctions/lcao/LCAOrbitalSet.h"
#include "QMCWaveFunctions/lcao/CuspCorrection.h"

#include "QMCWaveFunctions/SPOSetBuilderFactory.h"

namespace qmcplusplus
{
TEST_CASE("readCuspInfo", "[wavefunction]")
{
  OHMMS::Controller->initialize(0, NULL);
  Communicate* c = OHMMS::Controller;

  typedef OneDimGridBase<double> GridType;

  Matrix<CuspCorrectionParameters> info;
  int num_center       = 3;
  int orbital_set_size = 7;
  info.resize(num_center, orbital_set_size);

  bool okay = readCuspInfo("hcn_downdet.cuspInfo.xml", "downdet", orbital_set_size, info);
  REQUIRE(okay);

  // N
  REQUIRE(info(0, 0).redo == Approx(0.0));                   // redo
  REQUIRE(info(0, 0).C == Approx(0.0));                      // C
  REQUIRE(info(0, 0).sg == Approx(1.0));                     // sg
  REQUIRE(info(0, 0).Rc == Approx(0.0769130700800000));      // rc
  REQUIRE(info(0, 0).alpha[0] == Approx(2.29508580995773));  // a1
  REQUIRE(info(0, 0).alpha[1] == Approx(-7.00028778782666)); // a2
  REQUIRE(info(0, 0).alpha[2] == Approx(0.834942828252775)); // a3
  REQUIRE(info(0, 0).alpha[3] == Approx(-4.61597420905980)); // a4
  REQUIRE(info(0, 0).alpha[4] == Approx(31.6558091872316));  // a5

  // Spot check a few values from these centers
  // C
  REQUIRE(info(0, 6).C == Approx(0.0));        // C
  REQUIRE(info(0, 6).alpha[4] == Approx(0.0)); // a5

  // H
  REQUIRE(info(2, 4).alpha[4] == Approx(-404.733151049101)); // a5
}


TEST_CASE("applyCuspInfo", "[wavefunction]")
{
  OHMMS::Controller->initialize(0, NULL);
  Communicate* c = OHMMS::Controller;

  Libxml2Document doc;
  bool okay = doc.parse("hcn.structure.xml");
  REQUIRE(okay);
  xmlNodePtr root = doc.getRoot();
  Tensor<int, 3> tmat;
  tmat(0, 0) = 1;
  tmat(1, 1) = 1;
  tmat(2, 2) = 1;

  ParticleSet ions;
  XMLParticleParser parse_ions(ions, tmat);
  OhmmsXPathObject particleset_ion("//particleset[@name='ion0']", doc.getXPathContext());
  REQUIRE(particleset_ion.size() == 1);
  parse_ions.put(particleset_ion[0]);

  REQUIRE(ions.groups() == 3);
  REQUIRE(ions.R.size() == 3);
  ions.update();

  ParticleSet elec;
  XMLParticleParser parse_elec(elec, tmat);
  OhmmsXPathObject particleset_elec("//particleset[@name='e']", doc.getXPathContext());
  REQUIRE(particleset_elec.size() == 1);
  parse_elec.put(particleset_elec[0]);

  REQUIRE(elec.groups() == 2);
  REQUIRE(elec.R.size() == 14);

  elec.R = 0.0;

  elec.addTable(ions, DT_SOA);
  elec.update();

  Libxml2Document doc2;
  okay = doc2.parse("hcn.wfnoj.xml");
  REQUIRE(okay);
  xmlNodePtr root2 = doc2.getRoot();

  TrialWaveFunction psi(c);

  WaveFunctionComponentBuilder::PtclPoolType particle_set_map;
  particle_set_map["e"]    = &elec;
  particle_set_map["ion0"] = &ions;

  SPOSetBuilderFactory bf(elec, psi, particle_set_map);

  OhmmsXPathObject MO_base("//determinantset", doc2.getXPathContext());
  REQUIRE(MO_base.size() == 1);

  SPOSetBuilder* bb = bf.createSPOSetBuilder(MO_base[0]);
  REQUIRE(bb != NULL);

  OhmmsXPathObject slater_base("//determinant", doc2.getXPathContext());
  bb->loadBasisSetFromXML(MO_base[0]);
  SPOSet* sposet = bb->createSPOSet(slater_base[0]);

  LCAOrbitalSet* lcob = dynamic_cast<LCAOrbitalSet*>(sposet);
  REQUIRE(lcob != NULL);


  LCAOrbitalSet phi = LCAOrbitalSet(lcob->myBasisSet);
  phi.setOrbitalSetSize(lcob->OrbitalSetSize);
  phi.BasisSetSize = lcob->BasisSetSize;
  phi.setIdentity(false);

  LCAOrbitalSet eta = LCAOrbitalSet(lcob->myBasisSet);
  eta.setOrbitalSetSize(lcob->OrbitalSetSize);
  eta.BasisSetSize = lcob->BasisSetSize;
  eta.setIdentity(false);

  *(eta.C) = *(lcob->C);
  *(phi.C) = *(lcob->C);


  int num_center = 3;
  std::vector<bool> corrCenter(num_center, "true");

  // N is first atom
  int center_idx = 0;

  typedef QMCTraits::RealType RealType;

  splitPhiEta(center_idx, corrCenter, phi, eta);

  // 1S orbital on N
  CHECK((*phi.C)(0, 0) == Approx(1.00180500));
  CHECK((*eta.C)(0, 0) == Approx(0.0));

  int orbital_set_size = 7;
  Matrix<CuspCorrectionParameters> info;
  info.resize(num_center, orbital_set_size);
  okay = readCuspInfo("hcn_downdet.cuspInfo.xml", "downdet", orbital_set_size, info);

  REQUIRE(okay);
  Vector<double> xgrid;
  Vector<double> rad_orb;
  int ngrid = 10;
  xgrid.resize(ngrid);
  double dx = 0.12/ngrid;
  for (int i = 0; i < ngrid; i++)
  {
    xgrid[i] = dx * (i + 1);
  }

  rad_orb.resize(ngrid);

  int mo_idx = 0;
  computeRadialPhiBar(&elec, &ions, mo_idx, center_idx, &phi, xgrid, rad_orb, info(center_idx, mo_idx));

  // Comparisons generated from gen_cusp_corr.py
  //  Center  0  MO 0 rc =  0.07691307008
  REQUIRE(rad_orb[0] == Approx(9.1266186340)); // x = 0.012
  REQUIRE(rad_orb[1] == Approx(8.3939106599)); // x = 0.024
  REQUIRE(rad_orb[2] == Approx(7.7213972780)); // x = 0.036
  REQUIRE(rad_orb[3] == Approx(7.1039662640)); // x = 0.048
  REQUIRE(rad_orb[4] == Approx(6.5370601478)); // x = 0.06
  REQUIRE(rad_orb[5] == Approx(6.0165935481)); // x = 0.072
  REQUIRE(rad_orb[6] == Approx(5.5390213984)); // x = 0.084
  REQUIRE(rad_orb[7] == Approx(5.1023814795)); // x = 0.096
  REQUIRE(rad_orb[8] == Approx(4.7033287383)); // x = 0.108
  REQUIRE(rad_orb[9] == Approx(4.3370522377)); // x = 0.12


  mo_idx = 1;
  computeRadialPhiBar(&elec, &ions, mo_idx, center_idx, &phi, xgrid, rad_orb, info(center_idx, mo_idx));

  //  Center  0  MO 1 rc =  0.060909477888
  REQUIRE(rad_orb[0] == Approx(-0.0099816961)); // x = 0.012
  REQUIRE(rad_orb[1] == Approx(-0.0092950723)); // x = 0.024
  REQUIRE(rad_orb[2] == Approx(-0.0086498844)); // x = 0.036
  REQUIRE(rad_orb[3] == Approx(-0.0080440071)); // x = 0.048
  REQUIRE(rad_orb[4] == Approx(-0.0074778482)); // x = 0.06
  REQUIRE(rad_orb[5] == Approx(-0.0069529708)); // x = 0.072
  REQUIRE(rad_orb[6] == Approx(-0.0064707256)); // x = 0.084
  REQUIRE(rad_orb[7] == Approx(-0.0060313791)); // x = 0.096
  REQUIRE(rad_orb[8] == Approx(-0.0056312867)); // x = 0.108
  REQUIRE(rad_orb[9] == Approx(-0.0052652668)); // x = 0.12


  // Reset the MO matrices for another center

  *(eta.C) = *(lcob->C);
  *(phi.C) = *(lcob->C);


  // C is second atom
  center_idx = 1;
  splitPhiEta(center_idx, corrCenter, phi, eta);

  // 1S orbital on N
  CHECK((*phi.C)(0, 0) == Approx(0.0));
  CHECK((*eta.C)(0, 0) == Approx(1.00180500));

  mo_idx = 0;
  computeRadialPhiBar(&elec, &ions, mo_idx, center_idx, &phi, xgrid, rad_orb, info(center_idx, mo_idx));

  //  Center  1  MO 0 rc =  0.105
  REQUIRE(rad_orb[0] == Approx(0.0017535517)); // x = 0.012
  REQUIRE(rad_orb[1] == Approx(0.0016496533)); // x = 0.024
  REQUIRE(rad_orb[2] == Approx(0.0015544835)); // x = 0.036
  REQUIRE(rad_orb[3] == Approx(0.0014678130)); // x = 0.048
  REQUIRE(rad_orb[4] == Approx(0.0013891000)); // x = 0.06
  REQUIRE(rad_orb[5] == Approx(0.0013175785)); // x = 0.072
  REQUIRE(rad_orb[6] == Approx(0.0012523246)); // x = 0.084
  REQUIRE(rad_orb[7] == Approx(0.0011923038)); // x = 0.096
  REQUIRE(rad_orb[8] == Approx(0.0011364095)); // x = 0.108
  REQUIRE(rad_orb[9] == Approx(0.0010837868)); // x = 0.12


  removeSTypeOrbitals(corrCenter, *lcob);

  CHECK((*lcob->C)(0, 0) == Approx(0.0));
  CHECK((*lcob->C)(0, 1) == Approx(0.0));
  CHECK((*lcob->C)(0, 2) == Approx(0.0));
  CHECK((*lcob->C)(0, 3) != 0.0);

  SPOSetBuilderFactory::clear();
}

TEST_CASE("HCN MO with cusp", "[wavefunction]")
{
  OHMMS::Controller->initialize(0, NULL);
  Communicate* c = OHMMS::Controller;

  Libxml2Document doc;
  bool okay = doc.parse("hcn.structure.xml");
  REQUIRE(okay);
  xmlNodePtr root = doc.getRoot();
  Tensor<int, 3> tmat;
  tmat(0, 0) = 1;
  tmat(1, 1) = 1;
  tmat(2, 2) = 1;

  ParticleSet ions;
  XMLParticleParser parse_ions(ions, tmat);
  OhmmsXPathObject particleset_ion("//particleset[@name='ion0']", doc.getXPathContext());
  REQUIRE(particleset_ion.size() == 1);
  parse_ions.put(particleset_ion[0]);

  REQUIRE(ions.groups() == 3);
  REQUIRE(ions.R.size() == 3);
  ions.update();

  ParticleSet elec;
  XMLParticleParser parse_elec(elec, tmat);
  OhmmsXPathObject particleset_elec("//particleset[@name='e']", doc.getXPathContext());
  REQUIRE(particleset_elec.size() == 1);
  parse_elec.put(particleset_elec[0]);

  REQUIRE(elec.groups() == 2);
  REQUIRE(elec.R.size() == 14);

  elec.R = 0.0;

  elec.addTable(ions, DT_SOA);
  elec.update();

  Libxml2Document doc2;
  okay = doc2.parse("hcn.wfnoj.xml");
  REQUIRE(okay);
  xmlNodePtr root2 = doc2.getRoot();

  TrialWaveFunction psi(c);

  WaveFunctionComponentBuilder::PtclPoolType particle_set_map;
  particle_set_map["e"]    = &elec;
  particle_set_map["ion0"] = &ions;

  SPOSetBuilderFactory bf(elec, psi, particle_set_map);

  OhmmsXPathObject MO_base("//determinantset", doc2.getXPathContext());
  REQUIRE(MO_base.size() == 1);

  xmlSetProp(MO_base[0], (const xmlChar*)"cuspCorrection", (const xmlChar*)"yes");

  SPOSetBuilder* bb = bf.createSPOSetBuilder(MO_base[0]);
  REQUIRE(bb != NULL);

  OhmmsXPathObject slater_base("//determinant", doc2.getXPathContext());
  bb->loadBasisSetFromXML(MO_base[0]);
  SPOSet* sposet = bb->createSPOSet(slater_base[0]);


  SPOSet::ValueVector_t values;
  SPOSet::GradVector_t dpsi;
  SPOSet::ValueVector_t d2psi;
  values.resize(7);
  dpsi.resize(7);
  d2psi.resize(7);

  elec.R = 0.0;
  elec.update();
  ParticleSet::SingleParticlePos_t newpos;
  elec.makeMove(0, newpos);

  sposet->evaluate(elec, 0, values);

  // Values from gen_cusp_corr.py
  REQUIRE(values[0] == Approx(0.00945227));
  REQUIRE(values[1] == Approx(0.0200836));
  REQUIRE(values[2] == Approx(0.416375));
  REQUIRE(values[3] == Approx(-0.0885443));
  REQUIRE(values[4] == Approx(0.273159));
  REQUIRE(values[5] == Approx(0));
  REQUIRE(values[6] == Approx(0));

  // Put electron near N atom
  elec.R[0][0] = -1.09;
  elec.update();
  elec.makeMove(0, newpos);

  values = 0.0;
  sposet->evaluate(elec, 0, values);
  //std::cout << "values = " << values << std::endl;
  // Values from gen_cusp_corr.py
  REQUIRE(values[0] == Approx(9.5150713253));
  REQUIRE(values[1] == Approx(-0.0086731542));
  REQUIRE(values[2] == Approx(-1.6426151116));
  REQUIRE(values[3] == Approx(0.6569242017));
  REQUIRE(values[4] == Approx(0.9775522176));
  REQUIRE(values[5] == Approx(0.0000000000));
  REQUIRE(values[6] == Approx(0.0000000000));


  values = 0.0;
  sposet->evaluate(elec, 0, values, dpsi, d2psi);

  //std::cout << "values = " << values << std::endl;
  //std::cout << "dpsi = " << dpsi << std::endl;
  //std::cout << "d2psi = " << d2psi << std::endl;

  // Values from gen_cusp_corr.py
  REQUIRE(values[0] == Approx(9.5150713253));
  REQUIRE(values[1] == Approx(-0.0086731542));
  REQUIRE(values[2] == Approx(-1.6426151116));
  REQUIRE(values[3] == Approx(0.6569242017));
  REQUIRE(values[4] == Approx(0.9775522176));
  REQUIRE(values[5] == Approx(0.0000000000));
  REQUIRE(values[6] == Approx(0.0000000000));

  REQUIRE(dpsi[0][0] == Approx(-66.5007223213));
  REQUIRE(dpsi[0][1] == Approx(0.0000000000));
  REQUIRE(dpsi[0][2] == Approx(0.0000000000));
  REQUIRE(d2psi[0] == Approx(-21540.9990552510));

  REQUIRE(values[1] == Approx(-0.0086731542));
  REQUIRE(dpsi[1][0] == Approx(0.0616909346));
  REQUIRE(dpsi[1][1] == Approx(0.0000000000));
  REQUIRE(dpsi[1][2] == Approx(0.0000000000));
  REQUIRE(d2psi[1] == Approx(19.8720529007));


  SPOSet::ValueMatrix_t all_values;
  SPOSet::GradMatrix_t all_grad;
  SPOSet::ValueMatrix_t all_lap;
  all_values.resize(7, 7);
  all_grad.resize(7, 7);
  all_lap.resize(7, 7);


  sposet->evaluate_notranspose(elec, 0, 7, all_values, all_grad, all_lap);

  // Values from gen_cusp_corr.py
  REQUIRE(values[0] == Approx(9.5150713253));

  REQUIRE(all_values[0][0] == Approx(9.5150713253));
  REQUIRE(all_grad[0][0][0] == Approx(-66.5007223213));
  REQUIRE(all_grad[0][0][1] == Approx(0.0000000000));
  REQUIRE(all_grad[0][0][2] == Approx(0.0000000000));
  REQUIRE(all_lap[0][0] == Approx(-21540.9990552510));

  REQUIRE(all_values[0][1] == Approx(-0.0086731542));
  REQUIRE(all_grad[0][1][0] == Approx(0.0616909346));
  REQUIRE(all_grad[0][1][1] == Approx(0.0000000000));
  REQUIRE(all_grad[0][1][2] == Approx(0.0000000000));
  REQUIRE(all_lap[0][1] == Approx(19.8720529007));


  // Test the makeClone method
  SPOSet* sposet_clone = sposet->makeClone();

  sposet_clone->evaluate_notranspose(elec, 0, 7, all_values, all_grad, all_lap);

  // Values from gen_cusp_corr.py
  REQUIRE(values[0] == Approx(9.5150713253));

  REQUIRE(all_values[0][0] == Approx(9.5150713253));
  REQUIRE(all_grad[0][0][0] == Approx(-66.5007223213));
  REQUIRE(all_grad[0][0][1] == Approx(0.0000000000));
  REQUIRE(all_grad[0][0][2] == Approx(0.0000000000));
  REQUIRE(all_lap[0][0] == Approx(-21540.9990552510));

  REQUIRE(all_values[0][1] == Approx(-0.0086731542));
  REQUIRE(all_grad[0][1][0] == Approx(0.0616909346));
  REQUIRE(all_grad[0][1][1] == Approx(0.0000000000));
  REQUIRE(all_grad[0][1][2] == Approx(0.0000000000));
  REQUIRE(all_lap[0][1] == Approx(19.8720529007));

  SPOSetBuilderFactory::clear();
}

// Test case with multiple atoms of the same type
TEST_CASE("Ethanol MO with cusp", "[wavefunction]")
{
  OHMMS::Controller->initialize(0, NULL);
  Communicate* c = OHMMS::Controller;

  Libxml2Document doc;
  bool okay = doc.parse("ethanol.structure.xml");
  REQUIRE(okay);
  xmlNodePtr root = doc.getRoot();
  Tensor<int, 3> tmat;
  tmat(0, 0) = 1;
  tmat(1, 1) = 1;
  tmat(2, 2) = 1;

  ParticleSet ions;
  XMLParticleParser parse_ions(ions, tmat);
  OhmmsXPathObject particleset_ion("//particleset[@name='ion0']", doc.getXPathContext());
  REQUIRE(particleset_ion.size() == 1);
  parse_ions.put(particleset_ion[0]);

  REQUIRE(ions.groups() == 3);
  REQUIRE(ions.R.size() == 9);
  ions.update();

  ParticleSet elec;
  XMLParticleParser parse_elec(elec, tmat);
  OhmmsXPathObject particleset_elec("//particleset[@name='e']", doc.getXPathContext());
  REQUIRE(particleset_elec.size() == 1);
  parse_elec.put(particleset_elec[0]);

  REQUIRE(elec.groups() == 2);
  REQUIRE(elec.R.size() == 26);

  elec.R = 0.0;

  elec.addTable(ions, DT_SOA);
  elec.update();

  Libxml2Document doc2;
  okay = doc2.parse("ethanol.wfnoj.xml");
  REQUIRE(okay);
  xmlNodePtr root2 = doc2.getRoot();

  TrialWaveFunction psi(c);

  WaveFunctionComponentBuilder::PtclPoolType particle_set_map;
  particle_set_map["e"]    = &elec;
  particle_set_map["ion0"] = &ions;

  SPOSetBuilderFactory bf(elec, psi, particle_set_map);

  OhmmsXPathObject MO_base("//determinantset", doc2.getXPathContext());
  REQUIRE(MO_base.size() == 1);

  xmlSetProp(MO_base[0], (const xmlChar*)"cuspCorrection", (const xmlChar*)"yes");

  SPOSetBuilder* bb = bf.createSPOSetBuilder(MO_base[0]);
  REQUIRE(bb != NULL);

  OhmmsXPathObject slater_base("//determinant", doc2.getXPathContext());
  bb->loadBasisSetFromXML(MO_base[0]);
  SPOSet* sposet = bb->createSPOSet(slater_base[0]);


  SPOSet::ValueVector_t values;
  SPOSet::GradVector_t dpsi;
  SPOSet::ValueVector_t d2psi;
  values.resize(13);
  dpsi.resize(13);
  d2psi.resize(13);

  elec.R = 0.0;
  // Put electron near O atom
  elec.R[0][0] = -2.10;
  elec.R[0][1] = 0.50;

  elec.update();
  ParticleSet::SingleParticlePos_t newpos;
  elec.makeMove(0, newpos);

  sposet->evaluate(elec, 0, values);

  // Values from gen_cusp_corr.py
  REQUIRE(values[0] == Approx(4.3617329704));
  REQUIRE(values[1] == Approx(0.0014119853));
  REQUIRE(values[2] == Approx(0.0001156461));
  REQUIRE(values[3] == Approx(-0.6722670611));
  REQUIRE(values[4] == Approx(0.2762949842));
  REQUIRE(values[5] == Approx(0.2198735778));
  REQUIRE(values[6] == Approx(0.0659454461));
  REQUIRE(values[7] == Approx(0.2952071056));
  REQUIRE(values[8] == Approx(0.0322071389));
  REQUIRE(values[9] == Approx(0.0877981239));
  REQUIRE(values[10] == Approx(-0.2151873873));
  REQUIRE(values[11] == Approx(0.4250074750));
  REQUIRE(values[12] == Approx(0.0767950823));

  sposet->evaluate(elec, 0, values, dpsi, d2psi);

  REQUIRE(values[0] == Approx(4.3617329704));
  REQUIRE(values[1] == Approx(0.0014119853));
  REQUIRE(values[2] == Approx(0.0001156461));
  REQUIRE(values[3] == Approx(-0.6722670611));
  REQUIRE(values[4] == Approx(0.2762949842));
  REQUIRE(values[5] == Approx(0.2198735778));
  REQUIRE(values[6] == Approx(0.0659454461));
  REQUIRE(values[7] == Approx(0.2952071056));
  REQUIRE(values[8] == Approx(0.0322071389));
  REQUIRE(values[9] == Approx(0.0877981239));
  REQUIRE(values[10] == Approx(-0.2151873873));
  REQUIRE(values[11] == Approx(0.4250074750));
  REQUIRE(values[12] == Approx(0.0767950823));

  REQUIRE(dpsi[0][0] == Approx(-27.2844138432));
  REQUIRE(dpsi[0][1] == Approx(15.9958208598));
  REQUIRE(dpsi[0][2] == Approx(0.0195317131));
  REQUIRE(d2psi[0] == Approx(-293.2869628790));

  REQUIRE(dpsi[12][0] == Approx(1.7548511775));
  REQUIRE(dpsi[12][1] == Approx(2.2759333828));
  REQUIRE(dpsi[12][2] == Approx(-1.4878277937));
  REQUIRE(d2psi[12] == Approx(-4.3399821309));


  SPOSet::ValueMatrix_t all_values;
  SPOSet::GradMatrix_t all_grad;
  SPOSet::ValueMatrix_t all_lap;
  all_values.resize(13, 13);
  all_grad.resize(13, 13);
  all_lap.resize(13, 13);

  sposet->evaluate_notranspose(elec, 0, 7, all_values, all_grad, all_lap);

  REQUIRE(all_values[0][0] == Approx(4.3617329704));
  REQUIRE(all_grad[0][0][0] == Approx(-27.2844138432));
  REQUIRE(all_grad[0][0][1] == Approx(15.9958208598));
  REQUIRE(all_grad[0][0][2] == Approx(0.0195317131));
  REQUIRE(all_lap[0][0] == Approx(-293.2869628790));

  REQUIRE(all_values[0][11] == Approx(0.4250074750));
  REQUIRE(all_grad[0][11][0] == Approx(-0.3947036210));
  REQUIRE(all_grad[0][11][1] == Approx(0.9883840215));
  REQUIRE(all_grad[0][11][2] == Approx(1.7863218842));
  REQUIRE(all_lap[0][11] == Approx(-33.5202249813));

  SPOSetBuilderFactory::clear();
}

TEST_CASE("CuspCorrection He", "[wavefunction]")
{
  OHMMS::Controller->initialize(0, NULL);
  Communicate *c = OHMMS::Controller;

  ParticleSet elec;
  std::vector<int> agroup(2);
  agroup[0] = 1;
  agroup[1] = 1;
  elec.setName("e");
  elec.create(agroup);
  elec.R[0] = 0.0;

  SpeciesSet &tspecies = elec.getSpeciesSet();
  int upIdx = tspecies.addSpecies("u");
  int downIdx = tspecies.addSpecies("d");
  int massIdx = tspecies.addAttribute("mass");
  tspecies(massIdx, upIdx) = 1.0;
  tspecies(massIdx, downIdx) = 1.0;

  ParticleSet ions;
  ions.setName("ion0");
  ions.create(1);
  ions.R[0] = 0.0;
  SpeciesSet &ispecies = ions.getSpeciesSet();
  int heIdx = ispecies.addSpecies("He");
  ions.update();

  elec.addTable(ions,DT_SOA);
  elec.update();
  Libxml2Document doc;

  bool okay = doc.parse("he_sto3g.wfj.xml");
  REQUIRE(okay);
  xmlNodePtr root = doc.getRoot();

  TrialWaveFunction psi(c);

  WaveFunctionComponentBuilder::PtclPoolType particle_set_map;
  particle_set_map["e"]    = &elec;
  particle_set_map["ion0"] = &ions;

  SPOSetBuilderFactory bf(elec, psi, particle_set_map);

  OhmmsXPathObject MO_base("//determinantset", doc.getXPathContext());
  REQUIRE(MO_base.size() == 1);

  SPOSetBuilder* bb = bf.createSPOSetBuilder(MO_base[0]);
  REQUIRE(bb != NULL);

  OhmmsXPathObject slater_base("//determinant", doc.getXPathContext());
  bb->loadBasisSetFromXML(MO_base[0]);
  SPOSet* sposet = bb->createSPOSet(slater_base[0]);

  LCAOrbitalSet* lcob = dynamic_cast<LCAOrbitalSet*>(sposet);
  REQUIRE(lcob != NULL);

  typedef QMCTraits::RealType RealType;
  RealType rc = 0.1;
  int npts = 10;


  LCAOrbitalSet phi = LCAOrbitalSet(lcob->myBasisSet);
  phi.setOrbitalSetSize(lcob->OrbitalSetSize);
  phi.BasisSetSize = lcob->BasisSetSize;
  phi.setIdentity(false);

  LCAOrbitalSet eta = LCAOrbitalSet(lcob->myBasisSet);
  eta.setOrbitalSetSize(lcob->OrbitalSetSize);
  eta.BasisSetSize = lcob->BasisSetSize;
  eta.setIdentity(false);

  *(eta.C) = *(lcob->C);
  *(phi.C) = *(lcob->C);


  int num_center = 1;
  int center_idx = 0;
  std::vector<bool> corrCenter(num_center, "true");
  splitPhiEta(center_idx, corrCenter, phi, eta);

  OneMolecularOrbital phiMO(&elec, &ions, &phi);
  int orb_idx = 0;
  phiMO.changeOrbital(center_idx, orb_idx);

  OneMolecularOrbital etaMO(&elec, &ions, &eta);
  etaMO.changeOrbital(center_idx, orb_idx);

  //RealType phiAtRc = phiMO.phi(rc);
  //std::cout << "Phi at Rc = " << phiAtRc << std::endl;

  RealType dx = rc*1.2/npts;
  ValueVector_t pos(npts);
  ValueVector_t ELideal(npts);
  for (int i = 0; i < npts; i++) {
    pos[i] = (i+1.0)*dx;
  }

  RealType Z = 2.0;


  ValueType valRc;
  GradType gradRc;
  ValueType lapRc;
  phiMO.phi_vgl(rc, valRc, gradRc, lapRc);
   
  RealType C = 0.0;
  RealType valAtZero = phiMO.phi(0.0);
  RealType eta0 = etaMO.phi(0.0); 
  REQUIRE(eta0 == Approx(0.0)); // For He

  TinyVector<ValueType, 5> X;
  evalX(valRc, gradRc, lapRc, rc, Z, C, valAtZero, eta0, X);

  // From gen_cusp_corr.py
  REQUIRE(X[0] == Approx(-0.033436891110336));
  REQUIRE(X[1] == Approx(-0.653568722769692));
  REQUIRE(X[2] == Approx(-5.819488164002633));
  REQUIRE(X[3] == Approx(-2.000000000000000));
  REQUIRE(X[4] == Approx(-0.000396345019839));

  //TinyVector<ValueType, 5> alpha;
  CuspCorrectionParameters data;
  data.Rc = rc;
  data.C = C;
  X2alpha(X, rc, data.alpha);

  // From gen_cusp_corr.py
  REQUIRE(data.alpha[0] == Approx(-0.000396345019839));
  REQUIRE(data.alpha[1] == Approx(-2.000000000000000));
  REQUIRE(data.alpha[2] == Approx(56.659413909100188));
  REQUIRE(data.alpha[3] == Approx(-599.993590267020409));
  REQUIRE(data.alpha[4] == Approx(2003.589050855219512));

  CuspCorrection cusp(data);
  RealType Zeff = getZeff(Z, eta0, phiMO.phi(0.0));
  ValueVector_t ELorig(npts);
  RealType ELorigAtRc = getOriginalLocalEnergy(pos, Zeff, rc, phiMO, ELorig);
  
  // Original local energy
  // From gen_cusp_corr.py
  REQUIRE(ELorig[0] == Approx(-156.654088753559449));
  REQUIRE(ELorig[1] == Approx(-73.346068180623860));
  REQUIRE(ELorig[2] == Approx(-45.610385939854496));
  REQUIRE(ELorig[3] == Approx(-31.780236703094037));
  REQUIRE(ELorig[4] == Approx(-23.522092887496903));
  REQUIRE(ELorig[5] == Approx(-18.057926774366479));
  REQUIRE(ELorig[6] == Approx(-14.196956436578184));
  REQUIRE(ELorig[7] == Approx(-11.343582162638119));
  REQUIRE(ELorig[8] == Approx(-9.166698588588746));
  REQUIRE(ELorig[9] == Approx(-7.467419605354912));

  getIdealLocalEnergy(pos, Z, rc, ELorigAtRc, ELideal);

  //std::cout << "EL ideal = " << ELideal << std::endl;

  REQUIRE(ELorigAtRc== Approx(-10.5545686903018));
  // Ideal local energy
  // From gen_cusp_corr.py
  REQUIRE(ELideal[0] == Approx(-10.634967820121256));
  REQUIRE(ELideal[1] == Approx(-10.630023370426210));
  REQUIRE(ELideal[2] == Approx(-10.622438262014274));
  REQUIRE(ELideal[3] == Approx(-10.612683107096672));
  REQUIRE(ELideal[4] == Approx(-10.601175522784974));
  REQUIRE(ELideal[5] == Approx(-10.588284366044373));
  REQUIRE(ELideal[6] == Approx(-10.574333734041067));
  REQUIRE(ELideal[7] == Approx(-10.559606738039809));
  REQUIRE(ELideal[8] == Approx(-10.544349058872488));
  REQUIRE(ELideal[9] == Approx(-10.528772291863625));


  ValueVector_t ELcurr(npts);
  Zeff = getZeff(Z, eta0, cusp.phiBar(0.0, phiMO));
  getCurrentLocalEnergy(pos, Zeff, rc, ELorigAtRc, cusp, phiMO, ELcurr);
 // Current local energy
  // From gen_cusp_corr.py
  REQUIRE(ELcurr[0] == Approx(-130.055946501151169));
  REQUIRE(ELcurr[1] == Approx(-95.141127120655000));
  REQUIRE(ELcurr[2] == Approx(-66.353414984296620));
  REQUIRE(ELcurr[3] == Approx(-43.358705539629348));
  REQUIRE(ELcurr[4] == Approx(-26.111020098618731));
  REQUIRE(ELcurr[5] == Approx(-14.663411862083322));
  REQUIRE(ELcurr[6] == Approx(-9.047916153231112));
  REQUIRE(ELcurr[7] == Approx(-9.224544860293122));
  REQUIRE(ELcurr[8] == Approx(-9.166698588588746));
  REQUIRE(ELcurr[9] == Approx(-7.467419605354912));



  RealType chi2 = getELchi2(ELcurr, ELideal);
  //REQUIRE(chi2 == Approx(26452.9857955927));
  REQUIRE(chi2 == Approx(25854.2846426018)); 

#if 0
  MinimizePhiAtZero minPhi0(ELcurr, ELideal, cusp);

  RealType start_phi0 = phiMO.phi(0.0);
  Bracket_min_t<RealType> bracket = bracket_minimum(minPhi0, start_phi0);

  auto min_res = find_mininum(minPhi0, bracket);
#endif

  cusp.cparam.Rc = rc;
  chi2 = minimizeForPhiAtZero(cusp, phiMO, Z, eta0, pos, ELcurr, ELideal);
  RealType phi0 = cusp.phiBar(0.0, phiMO);
  std::cout << "phi0 = " << phi0 << std::endl;
  REQUIRE(phi0 == Approx(1.10489791512241));


  std::cout << "alpha = " << cusp.cparam.alpha << std::endl;
  // From gen_cusp_corr.py
  // brents method
  //REQUIRE(cusp.cparam.alpha[0] == Approx(0.099752931388349));
  //REQUIRE(cusp.cparam.alpha[1] == Approx(-2.000000000000000));
  //REQUIRE(cusp.cparam.alpha[2] == Approx(-3.430151935813066));
  //REQUIRE(cusp.cparam.alpha[3] == Approx(201.200620998489711));
  //REQUIRE(cusp.cparam.alpha[4] == Approx(-1000.889241390443090));
  double eps = 0.0003;

  // golden section search
  CHECK(cusp.cparam.alpha[0] == Approx(0.099752947906715));
  CHECK(cusp.cparam.alpha[1] == Approx(-2.000000000000000));
  CHECK(cusp.cparam.alpha[2] == Approx(-3.430161846832382).epsilon(eps));
  CHECK(cusp.cparam.alpha[3] == Approx(201.200753145414012).epsilon(eps));
  CHECK(cusp.cparam.alpha[4] == Approx(-1000.889736941408955).epsilon(eps));


#if 0

  std::cout << "Rc = 0.2" << std::endl;
  cusp.cparam.Rc = 0.2;

  dx = cusp.cparam.Rc*1.2/npts;
  for (int i = 0; i < npts; i++) {
    pos[i] = (i+1.0)*dx;
  }

  phi0 = minimizeForPhiAtZero(cusp, Z, eta0, pos, ELcurr, ELideal);
  std::cout << "phi0 = " << phi0 << std::endl;
  std::cout << "alpha = " << cusp.cparam.alpha << std::endl;
  CHECK(phi0 == Approx(1.21370813860125));
  CHECK(cusp.cparam.alpha[0] == Approx(0.193680245429731));
  CHECK(cusp.cparam.alpha[1] == Approx(-2.000000000000000));
  CHECK(cusp.cparam.alpha[2] == Approx(-2.457636715462080).epsilon(eps));
  CHECK(cusp.cparam.alpha[3] == Approx(44.150681525995104).epsilon(eps));
  CHECK(cusp.cparam.alpha[4] == Approx(-110.516657997161062).epsilon(eps));

  std::cout << "Rc = 0.3" << std::endl;
  cusp.cparam.Rc = 0.3;
  dx = cusp.cparam.Rc*1.2/npts;
  for (int i = 0; i < npts; i++) {
    pos[i] = (i+1.0)*dx;
  }
  phi0 = minimizeForPhiAtZero(cusp, Z, eta0, pos, ELcurr, ELideal);
  std::cout << "phi0 = " << phi0 << std::endl;
  std::cout << "alpha = " << cusp.cparam.alpha << std::endl;

  CHECK(phi0 == Approx(1.30756352742608));
  CHECK(cusp.cparam.alpha[0] == Approx(0.268165507600785));
  CHECK(cusp.cparam.alpha[1] == Approx(-2.000000000000000));
  CHECK(cusp.cparam.alpha[2] == Approx(-1.287424611600981).epsilon(eps));
  CHECK(cusp.cparam.alpha[3] == Approx(13.267353922388702).epsilon(eps));
  CHECK(cusp.cparam.alpha[4] == Approx(-22.573580351003290).epsilon(eps));


  std::cout << "Rc = 0.4" << std::endl;
  cusp.cparam.Rc = 0.4;
  dx = cusp.cparam.Rc*1.2/npts;
  for (int i = 0; i < npts; i++) {
    pos[i] = (i+1.0)*dx;
  }
  phi0 = minimizeForPhiAtZero(cusp, Z, eta0, pos, ELcurr, ELideal);
  std::cout << "phi0 = " << phi0 << std::endl;
  std::cout << "alpha = " << cusp.cparam.alpha << std::endl;

  CHECK(phi0 == Approx(1.36744420208896));


  CHECK(cusp.cparam.alpha[0] == Approx(0.312943445940581));
  CHECK(cusp.cparam.alpha[1] == Approx(-2.000000000000000));
  CHECK(cusp.cparam.alpha[2] == Approx(-0.266042849512209).epsilon(eps));
  CHECK(cusp.cparam.alpha[3] == Approx(2.771373038096821).epsilon(eps));
  CHECK(cusp.cparam.alpha[4] == Approx(-3.700333116423872).epsilon(eps));
#endif

  std::cout << "Rc = 0.5" << std::endl;
  cusp.cparam.Rc = 0.5;
  dx = cusp.cparam.Rc*1.2/npts;
  for (int i = 0; i < npts; i++) {
    pos[i] = (i+1.0)*dx;
  }

  ELorigAtRc = getOriginalLocalEnergy(pos, Zeff, cusp.cparam.Rc, phiMO, ELorig);
  getIdealLocalEnergy(pos, Z, cusp.cparam.Rc, ELorigAtRc, ELideal);
#if 0
  std::cout << "EL idea = " << std::endl;
  for (int i = 0; i < npts; i++) {
    std::cout << i << " " << ELideal[i] << std::endl;
  }
  std::cout << std::endl;
#endif
  

  phi0 = minimizeForPhiAtZero(cusp, phiMO, Z, eta0, pos, ELcurr, ELideal);
  std::cout << "phi0 = " << phi0 << std::endl;
  std::cout << "alpha = " << cusp.cparam.alpha << std::endl;

  CHECK(phi0 == Approx(1.38615192165457));
  CHECK(cusp.cparam.alpha[0] == Approx(0.326531501603760));
  CHECK(cusp.cparam.alpha[1] == Approx(-2.000000000000000));
  CHECK(cusp.cparam.alpha[2] == Approx(0.346373691083105));
  CHECK(cusp.cparam.alpha[3] == Approx(-0.685391898685964).epsilon(eps));
  CHECK(cusp.cparam.alpha[4] == Approx(0.649896548160093).epsilon(eps));


  minimizeForRc(cusp, phiMO, Z, rc, rc, eta0, pos, ELcurr, ELideal);
  std::cout << "alpha = " << cusp.cparam.alpha << std::endl;
  std::cout << "rc = " << cusp.cparam.Rc << std::endl;

  //CHECK(cusp.cparam.Rc == Approx());

  SPOSetBuilderFactory::clear();

}




} // namespace qmcplusplus
