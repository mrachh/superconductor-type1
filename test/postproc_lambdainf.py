from numpy import *
from numpy.linalg import *
from pylab import * 
from scipy.io import FortranFile



def get_bbpcomp(idzk):

    fname = '/mnt/home/mrachh/ceph/superconductor-type1-data/stell-data/statj_soln_20_60_norder8_ibg0_idzk'+str(idzk)+'.dat'
    
    print(fname)
    f = FortranFile(fname,'r')
    niter = int(f.read_ints(int32))
    rres = f.read_reals(float)
    soln = f.read_reals(float)
    npts = int(2400*45)
    bjmcomp = f.read_reals(float).reshape((3,npts),order="F")
    bbmcomp = f.read_reals(float).reshape((3,npts),order="F")
    bbpcomp = f.read_reals(float).reshape((3,npts),order="F")
    f.close()
    return bbpcomp

def get_bbpcomp_lambdainf():

    fname = '/mnt/home/mrachh/ceph/superconductor-type1-data/stell-data/statj_soln_lambdainf_20_60_norder8.dat'
    print(fname)
    f = FortranFile(fname,'r')
    niter = int(f.read_ints(int32))
    rres = f.read_reals(float)
    soln = f.read_reals(float)
    npts = int(2400*45)
    bbpcomp = f.read_reals(float).reshape((3,npts),order="F")
    f.close()
    return bbpcomp

npts = int(2400*45)
fname2 = 'triaskel_stellinfo.dat' 
f2 = FortranFile(fname2,'r')
uvpts = f2.read_reals(float).reshape((2,npts),order="F")
f2.close()


bbpcomp4 = get_bbpcomp(4)
bbpcomp = bbpcomp4
figure(1)
clf()
tricontourf(uvpts[0,:],uvpts[1,:],bbpcomp[0,:])
savefig('lambda0.pdf',bbox_inches='tight')

bbpcomp5 = get_bbpcomp(5)
bbpcomp = bbpcomp5

bbpcomp6 = get_bbpcomp(6)
bbpcomp = bbpcomp6
figure(2)
clf()
tricontourf(uvpts[0,:],uvpts[1,:],bbpcomp[0,:])
savefig('lambda2.pdf',bbox_inches='tight')

bbpcomp7 = get_bbpcomp(7)
bbpcomp = bbpcomp7


bbpcomp8 = get_bbpcomp(8)
bbpcomp = bbpcomp8
figure(3)
clf()
tricontourf(uvpts[0,:],uvpts[1,:],bbpcomp[0,:])
savefig('lambda4.pdf',bbox_inches='tight')

bbpcomp9 = get_bbpcomp(9)
bbpcomp = bbpcomp4
tricontourf(uvpts[0,:],uvpts[1,:],bbpcomp[0,:])


bbpcomp10 = get_bbpcomp(10)
bbpcomp = bbpcomp10
figure(4)
tricontourf(uvpts[0,:],uvpts[1,:],bbpcomp[0,:])
savefig('lambda6.pdf',bbox_inches='tight')

bbpcomp_ex = get_bbpcomp_lambdainf()
figure(5)
tricontourf(uvpts[0,:],uvpts[1,:],bbpcomp[0,:])
savefig('lambdainf.pdf',bbox_inches='tight')


errs = zeros(7)
errs[0] = norm(bbpcomp4-bbpcomp_ex)/norm(bbpcomp_ex) 
errs[1] = norm(bbpcomp5-bbpcomp_ex)/norm(bbpcomp_ex) 
errs[2] = norm(bbpcomp6-bbpcomp_ex)/norm(bbpcomp_ex) 
errs[3] = norm(bbpcomp7-bbpcomp_ex)/norm(bbpcomp_ex) 
errs[4] = norm(bbpcomp8-bbpcomp_ex)/norm(bbpcomp_ex) 
errs[5] = norm(bbpcomp9-bbpcomp_ex)/norm(bbpcomp_ex) 
errs[6] = norm(bbpcomp10-bbpcomp_ex)/norm(bbpcomp_ex) 

ii = 2.0**(-arange(7))

figure(6)
clf()
loglog(ii,errs,'k.',markersize=10)
xlim((0.01,1.2))
ylim((0.01,1.2))
savefig('lambdaconv.pdf',bbox_inches='tight')
show()
