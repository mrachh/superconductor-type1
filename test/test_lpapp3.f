      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: wts(:),rsigma(:)
      integer ipars(2)
      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      real *8 dpars(2)

      integer, allocatable :: norders(:),ixyzs(:),iptype(:)

      real *8 xyz_out(3),xyz_in(3)
      real *8, allocatable :: ffform(:,:,:),ffformex(:,:,:)
      real *8, allocatable :: ffforminv(:,:,:),ffformexinv(:,:,:)
      real *8, allocatable :: curv(:)
      real *8, allocatable :: spqp(:),grads0qm(:,:)
      real *8, allocatable :: s0laps0qp(:),s0laps0qm(:)

      complex *16, allocatable :: rhs(:)
      real *8, allocatable :: sigma(:), pot(:),pot1(:),rrhs(:)
      real *8, allocatable :: unm(:,:),xnm(:,:)
      real *8, allocatable :: drhs(:,:)
      real *8, allocatable :: sigma2(:)
      real *8, allocatable :: sigma0(:),pot0(:)
      real *8, allocatable :: errs(:)

      real *8, allocatable :: wnear(:)
      real *8, allocatable :: targs(:,:)
      real *8, allocatable :: cms(:,:),rads(:),rad_near(:)
      real *8, allocatable :: srcover(:,:),wover(:)
      real *8, allocatable :: f(:,:),df(:,:,:)
      real *8, allocatable :: w1(:,:),w2(:,:)
      real *8, allocatable :: blm(:,:),bmm(:,:)

      integer, allocatable :: col_ptr(:),row_ind(:)
      integer, allocatable :: ixyzso(:),novers(:)
      integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)

      complex *16 zalpha,zbeta,zgamma,zdelta,zeta,zteta,zk,ztetap
      complex *16 ztetam
      real *8 dalpha,dbeta,dgamma,ddelta,deta,dteta
      real *8 wtmp1(3)
      complex *16 fjvals(0:100),fhvals(0:100),fjder(0:100),fhder(0:100) 



      real *8 thet,phi,eps_gmres
      complex * 16 zpars(3)
      integer numit,niter
      character *100 title,dirname
      character *300 fname

      real *8, allocatable :: w(:,:)

      logical isout0,isout1

      complex *16 ztmp,ima

      data ima/(0.0d0,1.0d0)/


      call prini(6,13)

      done = 1
      pi = atan(done)*4


      zk = 4.4d0+ima*0.0d0
      zpars(1) = zk 
      zpars(2) = -ima*zk
      zpars(3) = 2.0d0

      
      xyz_in(1) = 0.11d0
      xyz_in(2) = 0.0d-5
      xyz_in(3) = 0.37d0

      xyz_out(1) = -3.5d0
      xyz_out(2) = 3.1d0
      xyz_out(3) = 20.1d0

      igeomtype = 1
      ipars(1) = 1 
      npatches12*(4**ipars(1))

      norder = 5 
      npols = (norder+1)*(norder+2)/2

      npts = npatches*npols
      allocate(srcvals(12,npts),srccoefs(9,npts))
      ifplot = 0

      call setup_geom(igeomtype,norder,npatches,ipars, 
     1       srcvals,srccoefs,ifplot,fname)

      allocate(norders(npatches),ixyzs(npatches+1),iptype(npatches))

      do i=1,npatches
        norders(i) = norder
        ixyzs(i) = 1 +(i-1)*npols
        iptype(i) = 1
      enddo

      print *, 'npts=',npts

      ixyzs(npatches+1) = 1+npols*npatches
      allocate(wts(npts))
      call get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,wts)


      allocate(sigma(6*npts),rhs(npts),pot(6*npts),rrhs(npts))
      allocate(ffform(2,2,npts))
      allocate(curv(npts))
      call get_mean_curvature(npatches,norders,ixyzs,iptype,npts, 
     1   srccoefs,srcvals,curv)

c
c       define rhs to be one of the ynm's
c
      nn = 3
      mm = 1
      nmax = nn
      allocate(w(0:nmax,0:nmax))
      call l3getsph(nmax,mm,nn,12,srcvals,rhs,npts,w)
c
c  set the value of dzk
c
      dzk = 0.19d0
cc      dzk = 1.0d0


      njh = 5
      ifder = 1
      rscale = 1.0d0
      zk = ima*dzk
      call prin2('zk=*',zk,2)
      call besseljs3d(njh,zk,rscale,fjvals,ifder,fjder)
      call h3dall(njh,zk,rscale,fhvals,ifder,fhder)

      zalpha = ima*((nn+1)*fjvals(nn)*fhvals(nn-1)+
     1   nn*fjvals(nn+1)*fhvals(nn)-zk*fjvals(nn+1)*fhvals(nn-1))
      dalpha = real(zalpha)
      call prin2('zalpha=*',zalpha,2)
      
      zbeta = fjvals(nn)*fhvals(nn)*zk*zk/ima/zk
      dbeta = real(zbeta)
      call prin2('zbeta=*',zbeta,2)
      
      zgamma = ima*zk*fjvals(nn)*fhvals(nn)*sqrt(nn*(nn+1.0d0))
      dgamma = real(zgamma)
      call prin2('zgamma=*',zgamma,2)
      
      zdelta = ima*((fjvals(nn)+zk*fjder(nn))*fhvals(nn)*zk + 
     1  (fhvals(nn)+zk*fhder(nn))*fjvals(nn)*zk)/2.0d0
      ddelta = real(zdelta)
      call prin2('zdelta=*',zdelta,2)
      
      zeta = ima*sqrt(nn*(nn+1.0d0))*(fjvals(nn)*fhvals(nn-1)-
     1   fjvals(nn+1)*fhvals(nn))
      deta = real(zeta)
      
      call prin2('zeta=*',zeta,2)
      
      zteta = ima*zk**2/2*(fjvals(nn)*fhder(nn)+fjder(nn)*fhvals(nn))
      call prin2('zteta=*',zteta,2)
      dteta = real(zteta)


      ztetap = ima*zk**2*fjvals(nn)*fhder(nn)
      dtetap = real(ztetap)
      call prin2('ztetap=*',ztetap,2)

      ztetam = ima*zk**2*fjder(nn)*fhvals(nn)
      dtetam = real(ztetam)
      call prin2('ztetam=*',ztetam,2)


      
      

c
c  set rcrhom to 1 for now for testing purposes
c

      rcrhom = 0
      rcrhop = 0
      rcrmum = 0

      rcqm = 0
      rcqp = 0
      rcrm = 0

      rcrhom = 2.9d0
      rcrhop = -0.17d0
      rcrmum = 1.3d0

      rcqm = 0.7d0
      rcqp = 1.0d0
      rcrm = -1.1d0

      do i=1,npts
        rrhs(i) = real(rhs(i))

        sigma(i) = rcrhom*rrhs(i)
        sigma(npts+i) = rcrhop*rrhs(i)
        sigma(2*npts+i) = rcrmum*rrhs(i)
        sigma(3*npts+i) = rcqm*rrhs(i)
        sigma(4*npts+i) = rcqp*rrhs(i)
        sigma(5*npts+i) = rcrm*rrhs(i)
      enddo

      allocate(unm(3,npts),xnm(3,npts))
      call surf_grad(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,rrhs,unm)
      do i=1,npts
        unm(1:3,i) = unm(1:3,i)/sqrt(nn*(nn+1.0d0))
        call cross_prod3d(srcvals(10,i),unm(1,i),xnm(1,i))
      enddo


      eps = 0.51d-5

      allocate(ipatch_id(npts),uvs_targ(2,npts))
      call get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts,
     1  ipatch_id,uvs_targ)
      
c
c       precompute near quadrature correction
c
c
      iptype_avg = floor(sum(iptype)/(npatches+0.0d0))
      norder_avg = floor(sum(norders)/(npatches+0.0d0))

      call get_rfacs(norder_avg,iptype_avg,rfac,rfac0)


      allocate(cms(3,npatches),rads(npatches),rad_near(npatches))

      call get_centroid_rads(npatches,norders,ixyzs,iptype,npts, 
     1     srccoefs,cms,rads)

C$OMP PARALLEL DO DEFAULT(SHARED) 
      do i=1,npatches
        rad_near(i) = rads(i)*rfac
      enddo
C$OMP END PARALLEL DO     

      ntarg = npts
      allocate(targs(3,npts))
C$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,npts 
        targs(1,i) = srcvals(1,i)
        targs(2,i) = srcvals(2,i)
        targs(3,i) = srcvals(3,i)
      enddo
C$OMP END PARALLEL DO      

c
c    find near quadrature correction interactions
c
      call findnearmem(cms,npatches,rad_near,3,targs,npts,nnz)

      allocate(row_ptr(npts+1),col_ind(nnz))
      
      call findnear(cms,npatches,rad_near,3,targs,npts,row_ptr, 
     1        col_ind)

      allocate(iquad(nnz+1)) 
      call get_iquad_rsc(npatches,ixyzs,ntarg,nnz,row_ptr,col_ind,
     1         iquad)



c
c    estimate oversampling for far-field, and oversample geometry
c

      ikerorder = -1
      allocate(novers(npatches),ixyzso(npatches+1))

      zpars = 0
      ndtarg = 3
      call get_far_order(eps,npatches,norders,ixyzs,iptype,cms,
     1    rads,npts,srccoefs,ndtarg,npts,targs,ikerorder,zpars,
     2    nnz,row_ptr,col_ind,rfac,novers,ixyzso)

      npts_over = ixyzso(npatches+1)-1

      allocate(srcover(12,npts_over),wover(npts_over))

      call oversample_geom(npatches,norders,ixyzs,iptype,npts, 
     1   srccoefs,srcvals,novers,ixyzso,npts_over,srcover)

      call get_qwts(npatches,novers,ixyzso,iptype,npts_over,
     1        srcover,wover)


c
c   compute near quadrature correction
c
      nquad = iquad(nnz+1)-1
      allocate(wnear(10*nquad))
      
C$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,10*nquad
        wnear(i) = 0
      enddo
C$OMP END PARALLEL DO    

      call prinf('finished generating near field info*',i,0)
      call prinf('finished generating far field orders*',i,0)
      call prinf('npts_over=*',npts_over,1)
      call prin2('eps=*',eps,1)

      iquadtype = 1

c      goto 1111
      call getnearquad_statj_gendeb(npatches,norders,
     1      ixyzs,iptype,npts,srccoefs,srcvals,
     1      eps,dzk,iquadtype,nnz,row_ptr,col_ind,
     1      iquad,rfac0,nquad,wnear)
      call prinf('finished generating near quadrature correction*',i,0)
 1111 continue      

      call prinf('entering layer potential eval*',i,0)
      call prinf('npts=*',npts,1)

      allocate(spqp(npts),grads0qm(3,npts),s0laps0qp(npts))
      allocate(s0laps0qm(npts))


      call statj_gendebproc_qpqm(npatches,norders,ixyzs,iptype,
     1  npts,srccoefs,srcvals,eps,nnz,row_ptr,col_ind,iquad,nquad,
     2  wnear,sigma,novers,npts_over,ixyzso,srcover,wover,curv,
     3  spqp,grads0qm,s0laps0qp,s0laps0qm)
      

      print *, "done with computation"

      do i=1,npts
        s0laps0qp(i) = s0laps0qp(i) - sigma(4*npts+i)/4
        s0laps0qm(i) = s0laps0qm(i) - sigma(3*npts+i)/4
      enddo

      erraqm = 0
      erraqp = 0
      ra = 0
      rqm  = -rcqm*(nn+0.0d0)*(nn+1.0d0)/(2*nn+1.0d0)**2
      rqp  = -rcqp*(nn+0.0d0)*(nn+1.0d0)/(2*nn+1.0d0)**2
      do i=1,npts
        erraqp = erraqp + (s0laps0qp(i) - rqp*rrhs(i))**2*wts(i)
        erraqm = erraqm + (s0laps0qm(i) - rqm*rrhs(i))**2*wts(i)
        rfac = s0laps0qp(i)/rqp/rrhs(i)
        if(i.lt.5) print *, s0laps0qp(i),rqp*rrhs(i),rfac
        ra = ra + rrhs(i)**2*wts(i)
      enddo
      erraqp = sqrt(erraqp/ra)
      erraqm = sqrt(erraqm/ra)
      call prin2('error in s0laps0qp=*',erraqp,1)
      call prin2('error in s0laps0qm=*',erraqm,1)
c
c  test error in spqp
c
      erra = 0
      ra = 0
      rspqp = -1.0d0/2/(2*nn+1.0d0)*rcqp
      do i=1,npts
        erra = erra + (rspqp*rrhs(i)-spqp(i))**2*wts(i)
        ra = ra + rrhs(i)**2*wts(i)
      enddo
      erra = sqrt(erra/ra)
      call prin2('error in spqp=*',erra,1)
c
c  test error in grads0qm
c

      rspqm = -1.0d0/2/(2*nn+1.0d0)*rcqm
      rugrads0 = sqrt((nn+0.0d0)*(nn+1.0d0))/(2*nn+1.0d0)*rcqm
      print *, 'rugrads0=',rugrads0
      erra = 0
      ra = 0
      do i=1,npts
        wtmp1(1:3) = rrhs(i)*rspqm*srcvals(10:12,i) + 
     1     unm(1:3,i)*rugrads0
        erra = erra + (wtmp1(1) - grads0qm(1,i))**2*wts(i)
        erra = erra + (wtmp1(2) - grads0qm(2,i))**2*wts(i)
        erra = erra + (wtmp1(3) - grads0qm(3,i))**2*wts(i)
        ra = ra + rrhs(i)**2*wts(i)
        if(i.lt.5) print *, wtmp1(1),grads0qm(1,i),
     1     wtmp1(1)/grads0qm(1,i)
        
      enddo

      erra = sqrt(erra/ra)
      call prin2('error in grads0qm=*',erra,1)




      stop
      end




      subroutine setup_geom(igeomtype,norder,npatches,ipars, 
     1    srcvals,srccoefs,ifplot,fname)
      implicit real *8 (a-h,o-z)
      integer igeomtype,norder,npatches,ipars(*),ifplot
      character (len=*) fname
      real *8 srcvals(12,*), srccoefs(9,*)
      real *8, allocatable :: uvs(:,:),umatr(:,:),vmatr(:,:),wts(:)

      real *8, pointer :: ptr1,ptr2,ptr3,ptr4
      integer, pointer :: iptr1,iptr2,iptr3,iptr4
      real *8, target :: p1(10),p2(10),p3(10),p4(10)
      real *8, allocatable, target :: triaskel(:,:,:)
      real *8, allocatable, target :: deltas(:,:)
      integer, allocatable :: isides(:)
      integer, target :: nmax,mmax

      procedure (), pointer :: xtri_geometry


      external xtri_stell_eval,xtri_sphere_eval
      
      npols = (norder+1)*(norder+2)/2
      allocate(uvs(2,npols),umatr(npols,npols),vmatr(npols,npols))
      allocate(wts(npols))

      call vioreanu_simplex_quad(norder,npols,uvs,umatr,vmatr,wts)

      if(igeomtype.eq.1) then
        itype = 2
        allocate(triaskel(3,3,npatches))
        allocate(isides(npatches))
        npmax = npatches
        ntri = 0
        call xtri_platonic(itype, ipars(1), npmax, ntri, 
     1      triaskel, isides)

        xtri_geometry => xtri_sphere_eval
        ptr1 => triaskel(1,1,1)
        ptr2 => p2(1)
        ptr3 => p3(1)
        ptr4 => p4(1)


        if(ifplot.eq.1) then
           call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2, 
     1         ptr3,ptr4, norder,'Triangulated surface of the sphere')
        endif


        call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,ptr3,ptr4,
     1     npols,uvs,umatr,srcvals,srccoefs)
      endif

      if(igeomtype.eq.2) then
        done = 1
        pi = atan(done)*4
        umin = 0
        umax = 2*pi
        vmin = 0
        vmax = 2*pi
        allocate(triaskel(3,3,npatches))
        nover = 0
        call xtri_rectmesh_ani(umin,umax,vmin,vmax,ipars(1),ipars(2),
     1     nover,npatches,npatches,triaskel)

        mmax = 2
        nmax = 1
        xtri_geometry => xtri_stell_eval

        allocate(deltas(-1:mmax,-1:nmax))
        deltas(-1,-1) = 0.17d0
        deltas(0,-1) = 0
        deltas(1,-1) = 0
        deltas(2,-1) = 0

        deltas(-1,0) = 0.11d0
        deltas(0,0) = 1
        deltas(1,0) = 4.5d0
        deltas(2,0) = -0.25d0

        deltas(-1,1) = 0
        deltas(0,1) = 0.07d0
        deltas(1,1) = 0
        deltas(2,1) = -0.45d0

        ptr1 => triaskel(1,1,1)
        ptr2 => deltas(-1,-1)
        iptr3 => mmax
        iptr4 => nmax

        if(ifplot.eq.1) then
           call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2, 
     1         iptr3,iptr4, norder,
     2         'Triangulated surface of the stellarator')
        endif

        call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,iptr3,iptr4,
     1     npols,uvs,umatr,srcvals,srccoefs)
      endif
      
      return  
      end


      subroutine test_exterior_pt(npatches,norder,npts,srcvals,
     1   srccoefs,wts,xyzout,isout)
c
c
c  this subroutine tests whether the pt xyzin, is
c  in the exterior of a surface, and also estimates the error
c  in representing e^{ir/2}/r and \grad e^{ir/2}/r \cdot n
c  centered at the interior point. Whether a point 
c  is in the interior or not is tested using Gauss' 
c  identity for the flux due to a point charge
c
c
c  input:
c    npatches - integer
c       number of patches
c    norder - integer
c       order of discretization
c    npts - integer
c       total number of discretization points on the surface
c    srccoefs - real *8 (9,npts)
c       koornwinder expansion coefficients of geometry info
c    xyzout -  real *8 (3)
c       point to be tested
c
c  output: 
c    isout - boolean
c      whether the target is in the interior or not
c

      implicit none
      integer npatches,norder,npts,npols
      real *8 srccoefs(9,npts),srcvals(12,npts),xyzout(3),wts(npts)
      real *8 tmp(3)
      real *8 dpars,done,pi
      real *8, allocatable :: rsurf(:),err_p(:,:) 
      integer ipars,norderhead,nd
      complex *16, allocatable :: sigma_coefs(:,:), sigma_vals(:,:)
      complex *16 zk,val

      integer ipatch,j,i
      real *8 ra,ds
      logical isout

      done = 1
      pi = atan(done)*4

      npols = (norder+1)*(norder+2)/2


      zk = 0

      ra = 0



      do ipatch=1,npatches
        do j=1,npols
          i = (ipatch-1)*npols + j
          call h3d_sprime(xyzout,srcvals(1,i),dpars,zk,ipars,val)
          call cross_prod3d(srcvals(4,i),srcvals(7,i),tmp)
          ds = sqrt(tmp(1)**2 + tmp(2)**2 + tmp(3)**2)
          ra = ra + real(val)*wts(i)
        enddo
      enddo

      if(abs(ra+4*pi).le.1.0d-3) isout = .false.
      if(abs(ra).le.1.0d-3) isout = .true.

      return
      end

   



      subroutine l3getsph(nmax,mm,nn,ndx,xyzs,ynms,npts,ynm)
      implicit real *8 (a-h,o-z)
      real *8 :: xyzs(ndx,npts)
      complex *16 ynms(npts),ima
      real *8 rat1(10000),rat2(10000)
      real *8 ynm(0:nmax,0:nmax)
      data ima/(0.0d0,1.0d0)/
  
      call ylgndrini(nmax, rat1, rat2)
  
      do i=1,npts
        x=xyzs(1,i)
        y=xyzs(2,i)
        z=xyzs(3,i)
        r=sqrt(x**2+y**2+z**2)
        call cart2polar(xyzs(1,i),r,theta,phi)
        ctheta = cos(theta)
        call ylgndrf(nmax, ctheta, ynm, rat1, rat2)
        ynms(i) = ynm(nn,abs(mm))*exp(ima*mm*phi)        
      enddo
       
      return
      end



