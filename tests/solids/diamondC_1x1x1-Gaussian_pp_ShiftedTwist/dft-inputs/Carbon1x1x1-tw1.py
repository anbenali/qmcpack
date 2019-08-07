#!/usr/bin/env python
import numpy
import h5py
from pyscf.pbc import gto, scf, dft
from pyscf.pbc import df 

cell = gto.Cell()
cell.a             = '''
         3.37316115       3.37316115       0.00000000
         0.00000000       3.37316115       3.37316115
         3.37316115       0.00000000       3.37316115'''
cell.atom = '''  
   C        0.00000000       0.00000000       0.00000000
   C        1.686580575      1.686580575      1.686580575 
            ''' 
cell.basis         = 'bfd-vdz'
cell.ecp           = 'bfd'
cell.unit          = 'B'
cell.drop_exponent = 0.15
cell.verbose       = 5
cell.charge        = 0
cell.spin          = 0
cell.build()

sp_twist=[0.333,0.333,0.333]
twist = numpy.asarray([0.333,0.333,0.333]) / 1.0
kmesh=[1,1,1]
kpts = cell.make_kpts((1,1,1), with_gamma_point=False,  wrap_around=True, scaled_center=twist)

print kpts
#axes = numpy.eye(3)*3.5668
#kaxes = 2*numpy.pi*numpy.linalg.inv(axes).T
#print (numpy.dot(kpts,axes)/numpy.pi/2)[13]
#exit()

mf = scf.KRHF(cell,kpts)
#mf.xc = 'pbe,pbe'

mf.exxdiv = 'ewald'

mf.max_cycle = 200




e_scf=mf.kernel()

ener = open('e_scf','w')
ener.write('%s\n' % (e_scf))
print('e_scf',e_scf)
ener.close()

title="S1-twist1"
from PyscfToQmcpack import savetoqmcpack
savetoqmcpack(cell,mf,title=title,kmesh=kmesh,kpts=kpts,sp_twist=kpts)

