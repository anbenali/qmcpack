
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
  for (int i = 0; i < ngrid; i++)
  {
    xgrid[i] = 0.012 * (i + 1);
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

} // namespace qmcplusplus