
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

#include "QMCWaveFunctions/MolecularOrbitals/LocalizedBasisSet.h"
#include "QMCWaveFunctions/MolecularOrbitals/SphericalBasisSet.h"
#include "QMCWaveFunctions/MolecularOrbitals/NGOBuilder.h"
#include "QMCWaveFunctions/MolecularOrbitals/CuspCorr.h"
#include "QMCWaveFunctions/MolecularOrbitals/LCOrbitalSetWithCorrection.h"

#include "QMCWaveFunctions/SPOSetBuilderFactory.h"

namespace qmcplusplus
{

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


  #ifdef ENABLE_SOA
    elec.addTable(ions,DT_SOA);
  #else
    elec.addTable(ions,DT_AOS);
  #endif
    elec.update();


  Libxml2Document doc;
  bool okay = doc.parse("he_sto3g.wfj.xml");
  REQUIRE(okay);
  xmlNodePtr root = doc.getRoot();

  TrialWaveFunction psi(c);

  WaveFunctionComponentBuilder::PtclPoolType particle_set_map;
  particle_set_map["e"] = &elec;
  particle_set_map["ion0"] = &ions;


  SPOSetBuilderFactory bf(elec, psi, particle_set_map);

  OhmmsXPathObject MO_base("//determinantset", doc.getXPathContext());
  REQUIRE(MO_base.size() == 1);

  SPOSetBuilder *bb = bf.createSPOSetBuilder(MO_base[0]);
  REQUIRE(bb != NULL);

  OhmmsXPathObject slater_base("//determinant", doc.getXPathContext());
  bb->loadBasisSetFromXML(MO_base[0]);
  SPOSet *sposet = bb->createSPOSet(slater_base[0]);


  typedef OneDimGridBase<double> GridType;
  typedef LCOrbitalSet<LocalizedBasisSet<SphericalBasisSet<NGOrbital, GridType>>, false> OrbType;
  OrbType *lcob = dynamic_cast<OrbType *>(sposet);
  REQUIRE(lcob != NULL);


  typedef CuspCorr<LocalizedBasisSet<SphericalBasisSet<NGOrbital, GridType>>>  CuspCorrType;
  typedef OrbitalSetTraits<double>::ValueVector_t ValueVector_t;

  double rc = 0.1;
  int npts = 10;
  CuspCorrType cusp(rc, npts, &elec, &ions);

  typedef NGOBuilder::CenteredOrbitalType COT;
  OrbType bs_phi(lcob->myBasisSet);
  bs_phi.setOrbitalSetSize(lcob->OrbitalSetSize);
  bs_phi.BasisSetSize = lcob->BasisSetSize;
  bs_phi.setIdentity(false);

  *(bs_phi.C) = *(lcob->C);

  OrbType bs_eta(lcob->myBasisSet);
  bs_eta.setOrbitalSetSize(lcob->OrbitalSetSize);
  bs_eta.BasisSetSize = lcob->BasisSetSize;
  bs_eta.setIdentity(false);
  *(bs_eta.C) = *(lcob->C);
  // For He in a minimal basis, there are only s-type orbitals
  (*bs_eta.C)(0,0) = 0.0;

  cusp.setPhiAndEta(&bs_phi, &bs_eta);

  cusp.curCenter = 0;
  cusp.curOrb = 0;
  cusp.Z = 2.0;
  cusp.computeValAtZero();

  ValueVector_t X(5);
  cusp.evalX(X);

  //std::cout << "X = " << X << std::endl;

  // From gen_cusp_corr.py
  REQUIRE(X[0] == Approx(-0.033436891110336));
  REQUIRE(X[1] == Approx(-0.653568722769692));
  REQUIRE(X[2] == Approx(-5.819488164002633));
  REQUIRE(X[3] == Approx(-2.000000000000000));
  REQUIRE(X[4] == Approx(-0.000396345019839));

  cusp.X2alpha(X);

  // From gen_cusp_corr.py
  REQUIRE(cusp.alpha[0] == Approx(-0.000396345019839));
  REQUIRE(cusp.alpha[1] == Approx(-2.000000000000000));
  REQUIRE(cusp.alpha[2] == Approx(56.659413909100188));
  REQUIRE(cusp.alpha[3] == Approx(-599.993590267020409));
  REQUIRE(cusp.alpha[4] == Approx(2003.589050855219512));


  cusp.allocateELspace();

  // compute original EL
  cusp.ELorigAtRc = cusp.getELorig();

  // compute current EL
  cusp.getELcurr();

  std::cout << "ELorigAtRc = " << cusp.ELorigAtRc << std::endl;
  // compute ideal EL
  cusp.getELideal(cusp.ELorigAtRc);

  // Grid for local energy evaluations
  // From gen_cusp_corr.py
  REQUIRE(cusp.pos[0] == Approx(0.012000000000000));
  REQUIRE(cusp.pos[1] == Approx(0.024000000000000));
  REQUIRE(cusp.pos[2] == Approx(0.036000000000000));
  REQUIRE(cusp.pos[3] == Approx(0.048000000000000));
  REQUIRE(cusp.pos[4] == Approx(0.060000000000000));
  REQUIRE(cusp.pos[5] == Approx(0.072000000000000));
  REQUIRE(cusp.pos[6] == Approx(0.084000000000000));
  REQUIRE(cusp.pos[7] == Approx(0.096000000000000));
  REQUIRE(cusp.pos[8] == Approx(0.108000000000000));
  REQUIRE(cusp.pos[9] == Approx(0.120000000000000));

  // Original local energy
  // From gen_cusp_corr.py
  REQUIRE(cusp.ELorig[0] == Approx(-156.654088753559449));
  REQUIRE(cusp.ELorig[1] == Approx(-73.346068180623860));
  REQUIRE(cusp.ELorig[2] == Approx(-45.610385939854496));
  REQUIRE(cusp.ELorig[3] == Approx(-31.780236703094037));
  REQUIRE(cusp.ELorig[4] == Approx(-23.522092887496903));
  REQUIRE(cusp.ELorig[5] == Approx(-18.057926774366479));
  REQUIRE(cusp.ELorig[6] == Approx(-14.196956436578184));
  REQUIRE(cusp.ELorig[7] == Approx(-11.343582162638119));
  REQUIRE(cusp.ELorig[8] == Approx(-9.166698588588746));
  REQUIRE(cusp.ELorig[9] == Approx(-7.467419605354912));


 // Current local energy
  // From gen_cusp_corr.py
  REQUIRE(cusp.ELcurr[0] == Approx(-130.055946501151169));
  REQUIRE(cusp.ELcurr[1] == Approx(-95.141127120655000));
  REQUIRE(cusp.ELcurr[2] == Approx(-66.353414984296620));
  REQUIRE(cusp.ELcurr[3] == Approx(-43.358705539629348));
  REQUIRE(cusp.ELcurr[4] == Approx(-26.111020098618731));
  REQUIRE(cusp.ELcurr[5] == Approx(-14.663411862083322));
  REQUIRE(cusp.ELcurr[6] == Approx(-9.047916153231112));
  REQUIRE(cusp.ELcurr[7] == Approx(-9.224544860293122));
  REQUIRE(cusp.ELcurr[8] == Approx(-9.166698588588746));
  REQUIRE(cusp.ELcurr[9] == Approx(-7.467419605354912));


  REQUIRE(cusp.ELorigAtRc== Approx(-10.5545686903018));
  // Ideal local energy
  // From gen_cusp_corr.py
  REQUIRE(cusp.ELideal[0] == Approx(-10.634967820121256));
  REQUIRE(cusp.ELideal[1] == Approx(-10.630023370426210));
  REQUIRE(cusp.ELideal[2] == Approx(-10.622438262014274));
  REQUIRE(cusp.ELideal[3] == Approx(-10.612683107096672));
  REQUIRE(cusp.ELideal[4] == Approx(-10.601175522784974));
  REQUIRE(cusp.ELideal[5] == Approx(-10.588284366044373));
  REQUIRE(cusp.ELideal[6] == Approx(-10.574333734041067));
  REQUIRE(cusp.ELideal[7] == Approx(-10.559606738039809));
  REQUIRE(cusp.ELideal[8] == Approx(-10.544349058872488));
  REQUIRE(cusp.ELideal[9] == Approx(-10.528772291863625));


  double chi2 = cusp.getchi2();
  REQUIRE(chi2 == Approx(25854.2846426019));


  CuspCorrectionParameters data;
  Vector<double> xgrid;
  Vector<double> rad_orb;
  xgrid.resize(1);
  xgrid[0] = 0.012;
  rad_orb.resize(1);

  rc = 0.5;
  cusp.Rc = rc;
  cusp.ELorigAtRc = cusp.getELorig();
  cusp.getELcurr();
  cusp.getELideal(cusp.ELorigAtRc);
  std::cout << "ELorigAtRc = " << cusp.ELorigAtRc << std::endl;
  std::cout << "Beta0 = " << cusp.beta0 << std::endl;
  std::cout << "pos = " << cusp.pos << std::endl;
  std::cout << "ELideal at rc=0.5" << std::endl << cusp.ELideal << std::endl;

  chi2 = cusp.execute(0, 0, 2.0, &bs_phi, &bs_eta, xgrid, rad_orb, "none", rc, data);
  std::cout << "Final chi2 = " << chi2 << std::endl;
  std::cout << "data alpha = " << data.alpha << std::endl;
  std::cout << "data sg = " << data.sg << std::endl;
  std::cout << "data C = " << data.C << std::endl;
  std::cout << "data Rc = " << data.Rc << std::endl;
  std::cout <<  std::endl;

#if 1
  cusp.Rc_max = rc;
  cusp.executeWithRCLoop(0, 0, 2.0, &bs_phi, &bs_eta, xgrid, rad_orb, "none", rc, data);
  std::cout << "After rc loop" << std::endl;
  std::cout << "Final chi2 = " << chi2 << std::endl;
  std::cout << "data alpha = " << data.alpha << std::endl;
  std::cout << "data sg = " << data.sg << std::endl;
  std::cout << "data Rc = " << data.Rc << std::endl;
#endif

  SPOSetBuilderFactory::clear();

}

TEST_CASE("readCuspInfo", "[wavefunction]")
{
  OHMMS::Controller->initialize(0, NULL);
  Communicate *c = OHMMS::Controller;

  typedef OneDimGridBase<double> GridType;
  typedef LCOrbitalSetWithCorrection<LocalizedBasisSet<SphericalBasisSet<NGOrbital, GridType>>, false> OrbType;

  OrbType orb;

  Matrix<CuspCorrectionParameters> info;
  int num_center = 3;
  int orbital_set_size = 7;
  info.resize(num_center, orbital_set_size);
  orb.cuspInfoFile = "hcn_downdet.cuspInfo.xml";
  orb.objectName = "downdet";
  orb.setOrbitalSetSize(orbital_set_size);
  bool okay = orb.readCuspInfo(info);
  REQUIRE(okay);

  // N
  REQUIRE(info(0,0).redo == Approx(0.0));
  REQUIRE(info(0,0).C == Approx(0.0));
  REQUIRE(info(0,0).sg == Approx(1.0));
  REQUIRE(info(0,0).Rc == Approx(0.0769130700800000));
  REQUIRE(info(0,0).alpha[0] == Approx(2.29508580995773));
  REQUIRE(info(0,0).alpha[1] == Approx(-7.00028778782666));
  REQUIRE(info(0,0).alpha[2] == Approx(0.834942828252775));
  REQUIRE(info(0,0).alpha[3] == Approx(-4.61597420905980));
  REQUIRE(info(0,0).alpha[4] == Approx(31.6558091872316));


  // Spot check a few values from these centers
  // C
  REQUIRE(info(0,6).C == Approx(0.0));
  REQUIRE(info(0,6).alpha[4] == Approx(0.0));

  // H
  REQUIRE(info(2,4).alpha[4] == Approx(-404.733151049101));

}

}
