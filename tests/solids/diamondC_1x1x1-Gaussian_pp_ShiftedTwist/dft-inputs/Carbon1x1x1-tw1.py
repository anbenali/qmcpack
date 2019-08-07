#!/usr/bin/env python
import numpy
import h5py
from pyscf.pbc import gto, scf, dft
from mpi4pyscf.pbc import df as mpidf

cell = gto.Cell()
cell.a             = '''
                    2.0 0 0  
                    0 2.0 0  
                    0 0 2.0 
                     '''
cell.dimension     = 3
cell.basis         = 'bfd-vdz'
cell.ecp           = 'bfd'
cell.unit          = 'B'
cell.a = numpy.eye(3)*3.5668
cell.atom = '''C     0.      0.      0.    
              C     0.8917  0.8917  0.8917
              C     1.7834  1.7834  0.    
              C     2.6751  2.6751  0.8917
              C     1.7834  0.      1.7834
              C     2.6751  0.8917  2.6751
              C     0.      1.7834  1.7834
              C     0.8917  2.6751  2.6751'''
cell.drop_exponent = 0.15
cell.verbose       = 5
cell.charge        = 0
cell.spin          = 0
cell.build()

sp_twist=[0.333,0.333,0.333]
twist = numpy.asarray([0.333,0.333,0.333]) / 1.0
kmesh=[1,1,1]
kpts = cell.make_kpts((1,1,1), with_gamma_point=False,  wrap_around=True, scaled_center=twist)

mf = scf.KRHF(cell,kpts)
mf.with_df = mpidf.FFTDF(cell,kpts)
#mf.with_df_cderi_to_save = 'df_ints.h5'   # new
#mf.with_df.build()
#mf.xc = 'pbe,pbe'

mf.with_df._cderi = 'df_ints.h5'
mf.exxdiv = 'ewald'


mf.chkfile = "tw1_ds30.dump"                         # store checkpoint file in scf.dump
#dm = mf.from_chk('tw1_ds30.dump')
#mf = scf.addons.smearing2_(mf, sigma=.0001, method='fermi', dS=30)
mf.max_cycle = 200




e_scf=mf.kernel()

ener = open('e_scf','w')
ener.write('%s\n' % (e_scf))
print('e_scf',e_scf)
ener.close()

title="S1-twist1"
from PyscfToQmcpack import savetoqmcpack
savetoqmcpack(cell,mf,title=title,kmesh=kmesh,kpts=kpts,sp_twist=sp_twist)

