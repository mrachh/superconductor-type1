      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: wts(:),rsigma(:)
      integer ipars(2)

      integer, allocatable :: norders(:),ixyzs(:),iptype(:)

      real *8 xyz_out(3),xyz_in(3)
      real *8, allocatable :: ffform(:,:,:),ffformex(:,:,:)
      real *8, allocatable :: ffforminv(:,:,:),ffformexinv(:,:,:)

      complex *16, allocatable :: sigma(:),rhs(:)
      real *8, allocatable :: drhs(:),drhsuv(:,:),rhsdiv(:),rhsdiv2(:)
      real *8, allocatable :: drhs_cart(:,:)
      real *8, allocatable :: drhs_cart2(:,:)
      real *8, allocatable :: errs(:)
      real *8 thet,phi,eps_gmres
      real *8, allocatable :: avals(:,:),bvals(:,:)
      real *8, allocatable :: awts(:),bwts(:)
      real *8, allocatable :: auv(:,:),buv(:,:)
      real *8, allocatable :: hvec(:,:,:)
      real *8, allocatable :: hinterpa(:,:,:),hinterpb(:,:,:)
      real *8, allocatable :: rrhs1(:),rrhs2(:)
      real *8 raint(2), rbint(2),rhvec(2)
      integer, allocatable :: apatches(:),bpatches(:)
      complex * 16 zpars(3)
      integer numit,niter
      character *100 title,dirname
      character *300 fname

      integer ipatch_id
      real *8 wtmp1(3),wtmp2(3)
      real *8 uvs_targ(2)
      real *8, allocatable :: w(:,:)

      logical isout0,isout1

      complex *16 pot,potex,ztmp,ima

      data ima/(0.0d0,1.0d0)/


      call prini(6,13)

      done = 1
      pi = atan(done)*4



      igeomtype = 4
      ipars(1) = 3*2
      ipars(2) = 9*2
      npatches = 2*ipars(1)*ipars(2)

      if(igeomtype.eq.2) fname = 'stell.vtk'
      if(igeomtype.eq.4) fname = 'torus.vtk'
      if(igeomtype.eq.5) fname = 'torus2.vtk'

      norder = 8 
      npols = (norder+1)*(norder+2)/2

      npts = npatches*npols
      allocate(srcvals(12,npts),srccoefs(9,npts))
      ifplot = 1

      call setup_geom(igeomtype,norder,npatches,ipars, 
     1       srcvals,srccoefs,ifplot,fname)

      allocate(norders(npatches),ixyzs(npatches+1),iptype(npatches))

      do i=1,npatches
        norders(i) = norder
        ixyzs(i) = 1 +(i-1)*npols
        iptype(i) = 1
      enddo

      ixyzs(npatches+1) = 1+npols*npatches
      allocate(wts(npts))
      call get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,wts)


      m = 16
      na = ipars(2)*m
      nb = ipars(1)*m
      allocate(avals(9,na),awts(na),auv(2,na),apatches(na))
      allocate(bvals(9,nb),bwts(nb),buv(2,nb),bpatches(nb))
      call get_ab_cycles_torusparam(npatches,norders,ixyzs,iptype,
     1   npts,srccoefs,srcvals,ipars,m,na,avals,awts,apatches,auv,
     2   nb,bvals,bwts,bpatches,buv)


      rvecb = 0
      rscb = 0

      rveca = 0
      rsca = 0
      
      do i=1,nb
        rvecb = rvecb + (bvals(3,i)*bvals(4,i) + 
     1        bvals(1,i)*bvals(5,i) + bvals(2,i)*bvals(6,i))*bwts(i)
        ds = sqrt(bvals(4,i)**2 + bvals(5,i)**2 + bvals(6,i)**2)
        rscb = rscb + bwts(i)*ds
      enddo

      do i=1,na
        rveca = rveca + (avals(3,i)*avals(4,i) + 
     1        avals(1,i)*avals(5,i) + avals(2,i)*avals(6,i))*awts(i)
        ds = sqrt(avals(4,i)**2 + avals(5,i)**2 + avals(6,i)**2)
        rsca = rsca + awts(i)*ds
      enddo

      allocate(hvec(3,npts,2))

      ra1 = 0
      rb1 = 0
      ra2 = 0
      rb2 = 0
      do i=1,npts
        rr1 = srcvals(1,i)**2 + srcvals(2,i)**2
        hvec(1,i,1) = -srcvals(2,i)/rr1
        hvec(2,i,1) = srcvals(1,i)/rr1
        hvec(3,i,1) = 0
        call cross_prod3d(srcvals(10,i),hvec(1,i,1),hvec(1,i,2))
      enddo

      allocate(hinterpa(3,na,2),hinterpb(3,nb,2))


      call fun_surf_interp(3,npatches,norders,ixyzs,iptype,
     1   npts,hvec(1,1,1),na,apatches,auv,hinterpa(1,1,1))
      call fun_surf_interp(3,npatches,norders,ixyzs,iptype,npts,
     1   hvec(1,1,2),na,apatches,auv,hinterpa(1,1,2))

      call fun_surf_interp(3,npatches,norders,ixyzs,iptype,npts,
     1   hvec(1,1,1),nb,bpatches,buv,hinterpb(1,1,1)) 
      call fun_surf_interp(3,npatches,norders,ixyzs,iptype,npts,
     1   hvec(1,1,2),nb,bpatches,buv,hinterpb(1,1,2))
      
      do j=1,2
        raint(j) = 0
        do i=1,na
          raint(j) = raint(j) + (hinterpa(1,i,j)*avals(4,i) + 
     1      hinterpa(2,i,j)*avals(5,i) + hinterpa(3,i,j)*avals(6,i))*
     2      awts(i)
        enddo
      enddo

      call prin2('raint=*',raint,2)
      do j=1,2
        rbint(j) = 0
        do i=1,nb
          rbint(j) = rbint(j) + (hinterpb(1,i,j)*bvals(4,i) + 
     1      hinterpb(2,i,j)*bvals(5,i) + hinterpb(3,i,j)*bvals(6,i))*
     2      bwts(i)
        enddo
      enddo

      call prin2('rbint=*',rbint,2)
      do i=1,npts
        hvec(1:3,i,1) = hvec(1:3,i,1)/rbint(1)
        hvec(1:3,i,2) = hvec(1:3,i,2)/raint(2)
      enddo


      call fun_surf_interp(3,npatches,norders,ixyzs,iptype,
     1   npts,hvec(1,1,1),na,apatches,auv,hinterpa(1,1,1))
      call fun_surf_interp(3,npatches,norders,ixyzs,iptype,npts,
     1   hvec(1,1,2),na,apatches,auv,hinterpa(1,1,2))

      call fun_surf_interp(3,npatches,norders,ixyzs,iptype,npts,
     1   hvec(1,1,1),nb,bpatches,buv,hinterpb(1,1,1)) 
      call fun_surf_interp(3,npatches,norders,ixyzs,iptype,npts,
     1   hvec(1,1,2),nb,bpatches,buv,hinterpb(1,1,2))
      

      do j=1,2
        raint(j) = 0
        do i=1,na
          raint(j) = raint(j) + (hinterpa(1,i,j)*avals(4,i) + 
     1      hinterpa(2,i,j)*avals(5,i) + hinterpa(3,i,j)*avals(6,i))*
     2      awts(i)
        enddo
      enddo

      call prin2('raint renorm=*',raint,2)
      do j=1,2
        rbint(j) = 0
        do i=1,nb
          rbint(j) = rbint(j) + (hinterpb(1,i,j)*bvals(4,i) + 
     1      hinterpb(2,i,j)*bvals(5,i) + hinterpb(3,i,j)*bvals(6,i))*
     2      bwts(i)
        enddo
      enddo
      call prin2('rbint renorm=*',rbint,2)

c
c  now compute wedge product of the harmonic fields with themselves
c
c
      do j=1,2
        rhvec(j) = 0
        if(j.eq.1) jj=2
        if(j.eq.2) jj=1
        do i=1,npts
          call cross_prod3d(hvec(1,i,j),hvec(1,i,jj),wtmp1)
          rhvec(j) = rhvec(j) + (srcvals(10,i)*wtmp1(1) + 
     1       srcvals(11,i)*wtmp1(2) + srcvals(12,i)*wtmp1(3))*wts(i)
        enddo
      enddo
      call prin2('rhvec=*',rhvec,2)

      rprod = raint(2)*rbint(1)
      erra = abs(rprod - rhvec(2))

      ra = 0
      do i=1,na
        ra = ra + awts(i)*sqrt(avals(4,i)**2 + avals(5,i)**2 +
     1      avals(6,i)**2)

      enddo
      call prin2('length of a cycle=*',ra,1)
      rb = 0
      do i=1,nb
        rb = rb + bwts(i)*sqrt(bvals(4,i)**2 + bvals(5,i)**2 + 
     1    bvals(6,i)**2)
      enddo

      call prin2('length of b cycle=*',rb,1)



      return
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


      external xtri_stell_eval,xtri_sphere_eval,xtri_wtorus_eval
      
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
        vmin = 2*pi
        vmax = 0
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

      if(igeomtype.eq.3) then
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
        call prinf('npatches=*',npatches,1)
         
        p1(1) = 1
        p1(2) = 2
        p1(3) = 0.25d0

        p2(1) = 1.0d0
        p2(2) = 1.0d0
        p2(3) = 1.0d0

c
c         numberof oscillations
c
        p4(1) = 3.0d0


        ptr1 => triaskel(1,1,1)
        ptr2 => p1(1)
        ptr3 => p2(1)
        ptr4 => p4(1)
        xtri_geometry => xtri_wtorus_eval
        if(ifplot.eq.1) then
           call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2, 
     1         ptr3,ptr4, norder,
     2         'Triangulated surface of the wtorus')
        endif


        call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,ptr3,ptr4,
     1     npols,uvs,umatr,srcvals,srccoefs)
      endif
      
      if(igeomtype.eq.4) then
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
        call prinf('npatches=*',npatches,1)
         
        p1(1) = 1.0d0
        p1(2) = 1.75d0
        p1(3) = 0.25d0

        p2(1) = 1.0d0
        p2(2) = 1.0d0
        p2(3) = 1.0d0

c
c         number of oscillations
c
        p4(1) = 0.0d0


        ptr1 => triaskel(1,1,1)
        ptr2 => p1(1)
        ptr3 => p2(1)
        ptr4 => p4(1)
        xtri_geometry => xtri_wtorus_eval
        if(ifplot.eq.1) then
           call xtri_vtk_surf(fname,npatches,xtri_geometry, ptr1,ptr2, 
     1         ptr3,ptr4, norder,
     2         'Triangulated surface of the torus')
        endif


        call getgeominfo(npatches,xtri_geometry,ptr1,ptr2,ptr3,ptr4,
     1     npols,uvs,umatr,srcvals,srccoefs)
      endif

      if(igeomtype.eq.5) then
        done = 1
        pi = atan(done)*4
        umin = 0
        umax = 2*pi
        vmin = 2*pi
        vmax = 0
        allocate(triaskel(3,3,npatches))
        nover = 0
        call xtri_rectmesh_ani(umin,umax,vmin,vmax,ipars(1),ipars(2),
     1     nover,npatches,npatches,triaskel)

        mmax = 2
        nmax = 1
        xtri_geometry => xtri_stell_eval

        allocate(deltas(-1:mmax,-1:nmax))
        deltas(-1,-1) = 0.17d0*0
        deltas(0,-1) = 0
        deltas(1,-1) = 0
        deltas(2,-1) = 0

        deltas(-1,0) = 0.11d0*0
        deltas(0,0) = 1*0
        deltas(1,0) = 2.0d0
        deltas(2,0) = -0.25d0*0

        deltas(-1,1) = 0
        deltas(0,1) = 1.0d0
        deltas(1,1) = 0
        deltas(2,1) = -0.45d0*0


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
      integer ipars,norderhead,nd

      integer ipatch,j,i
      real *8 ra,ds,r,dx,dy,dz,val
      logical isout

      done = 1
      pi = atan(done)*4

      npols = (norder+1)*(norder+2)/2
      call prin2('xyzout=*',xyzout,3)


      ra = 0
      do ipatch=1,npatches
        do j=1,npols
          i = (ipatch-1)*npols + j
          dx = xyzout(1) - srcvals(1,i)
          dy = xyzout(2) - srcvals(2,i)
          dz = xyzout(3) - srcvals(3,i)

          r = sqrt(dx**2 + dy**2 + dz**2)
          
          val = dx*srcvals(10,i) + dy*srcvals(11,i) + dz*srcvals(12,i)
          val = val/r**3 
          ra = ra + val*wts(i)
        enddo
      enddo
      print *, "ra=",ra
cc      call prin2('ra=*',ra,1)

      if(abs(ra+4*pi).le.2*pi) isout = .false.
      if(abs(ra+4*pi).gt.2*pi) isout = .true.

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






