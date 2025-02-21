from numpy import *
from numpy.linalg import *
from pylab import * 
from scipy.io import FortranFile



def get_bbpcomp(idzk, ifb):

    if ifb:
        fname = './static_curr_results/ap_0_bm_1_iref_4_idzk_'+str(idzk)+'/bbpcomp.dat'
    else: 
        fname = './static_curr_results/ap_1_bm_0_iref_4_idzk_'+str(idzk)+'/bbpcomp.dat'
    npts = int(1024*2*45)
        
    bbpcomp = loadtxt(fname) 
    return bbpcomp.transpose()

def get_bbpcomp_lambdainf(idzk, ifb):

    if ifb:
        fname = '/mnt/home/mrachh/ceph/superconductor-type1-data/thinshell-data/statj_soln_lambdainf_32_16_norder8_ifb1.dat'
    else:
        fname = '/mnt/home/mrachh/ceph/superconductor-type1-data/thinshell-data/statj_soln_lambdainf_32_16_norder8_ifb0.dat'
    npts = int(1024*45)
    print(fname)
    f = FortranFile(fname,'r')
    ipars = f.read_ints(int32).reshape((2,1),order="F")
    iref = int(f.read_ints(int32))
    ifb = int(f.read_ints(int32))
    tsolve = f.read_reals(float)
    niter = int(f.read_ints(int32))
    rres = f.read_reals(float)
    soln = f.read_reals(float)
    bbpcomp = f.read_reals(float).reshape((3,npts),order="F")
    
    f.close()
    return bbpcomp


ifb = 0

errs_b = zeros(7)
errs_a = zeros(7)

fname = './outer-torus-wts_iref4.dat'
f = FortranFile(fname, 'r')
wts_a = f.read_reals(float)
f.close()

fname = './inner-torus-wts_iref4.dat'
f = FortranFile(fname, 'r')
wts_b = f.read_reals(float)
f.close()

for i in arange(7):
    bbpcomp_inf = get_bbpcomp_lambdainf(i+4, ifb)
    npts = shape(bbpcomp_inf)[1]
    bbpcomp = get_bbpcomp(i+4,ifb)
    bbpcomp = bbpcomp[:,0:npts]

    bbpcomp = multiply(bbpcomp, sqrt(wts_a))
    bbpcomp_inf = multiply(bbpcomp_inf, sqrt(wts_a))
    errs_a[i] = norm(bbpcomp_inf - bbpcomp)/norm(bbpcomp)
    

ifb = 1

for i in arange(7):
    bbpcomp_inf = get_bbpcomp_lambdainf(i+4, ifb)
    npts = shape(bbpcomp_inf)[1]
    bbpcomp = get_bbpcomp(i+4,ifb)
    bbpcomp = bbpcomp[:,npts::]

    bbpcomp = multiply(bbpcomp, sqrt(wts_b))
    bbpcomp_inf = multiply(bbpcomp_inf, sqrt(wts_b))
    errs_b[i] = norm(bbpcomp_inf - bbpcomp)/norm(bbpcomp)
    
    

ii = 2.0**(-arange(7))

figure(6)
clf()
loglog(ii,errs_a,'k.',markersize=10)
xlim((0.01,1.2))
ylim((0.01,1.2))
savefig('lambdaconv_thinshell_a.pdf',bbox_inches='tight')
show()

figure(7)
clf()
loglog(ii,errs_b,'k.',markersize=10)
xlim((0.01,1.2))
#ylim((0.01,1.2))
savefig('lambdaconv_thinshell_b.pdf',bbox_inches='tight')
show()
