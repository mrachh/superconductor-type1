      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: wts(:),rsigma(:)
      integer ipars(2)
      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      real *8 dpars(2)

      integer, allocatable :: norders(:),ixyzs(:),iptype(:)

      real *8 xyz_out(3),xyz_in(3),dxyz(3)

      complex *16, allocatable :: zbbm(:,:),zbbp(:,:),zjm(:,:)
      complex *16, allocatable :: zjm0(:,:),zjm1(:)
      real *8 bbpex(3),bbpc(3)
      real *8, allocatable :: bbm(:,:),bbp(:,:),jm(:,:)
      real *8, allocatable :: ptmp(:)
      real *8, allocatable :: rhs(:,:),soln(:,:)
      real *8, allocatable :: errs(:)

      real *8 thet,phi,eps_gmres
      real *8 df2(3),rinttmp(6),cf2(1)
      complex *16 vf1(3),vf2(3),w1(3),w2(3)
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

      
      xyz_in(1) = 0.11d0
      xyz_in(2) = 0.01d-5
      xyz_in(3) = 0.07d0

      xyz_out(1) = -3.5d0
      xyz_out(2) = 3.1d0
      xyz_out(3) = 5.1d0

      igeomtype = 1
      ipars(1) = 2 
      npatches=12*(4**ipars(1))

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

      eps = 0.51d-5

      allocate(rhs(npts,6),soln(npts,6))
      allocate(zbbm(3,npts),zbbp(3,npts),zjm(3,npts))
      allocate(bbm(3,npts),bbp(3,npts),jm(3,npts),ptmp(npts))

      vf1(1) = 1.1d0
      vf1(2) = -0.2d0
      vf1(3) = 0.33d0
      
      vf2(1) = 0.61d0
      vf2(2) = 0.19d0
      vf2(3) = -0.71d0

      df2 = vf2
      cf2 = 1.1d0

      ifinit = 0
      dzk = 1.0d0
      ztmp = dzk*ima
      ifinit = 0
      
      do i=1,npts
        bbm(1:3,i) = 0
        jm(1:3,i) = 0
      enddo


      allocate(zjm0(3,npts),zjm1(npts))
      ztmp = 0.0d0

      thresh = 1.0d-16

      call prin2('df2=*',df2,3)
      bbp = 0
      call prin2('xyz_in=*',xyz_in,3)

cc      call l3ddirectdg(1,xyz_in,df2,1,srcvals(1:3,1:npts),npts,ptmp,
cc     1  bbp,thresh)
      call l3ddirectcg(1,xyz_in,cf2,1,srcvals(1:3,1:npts),npts,ptmp,
     1  bbp,thresh)
      call prin2('ptmp=*',ptmp,24)

cc      do i=1,npts
cc        dx = srcvals(1,i) - xyz_in(1)
cc        dy = srcvals(2,i) - xyz_in(2)
cc        dz = srcvals(3,i) - xyz_in(3)
cc        dxyz(1) = dx
cc        dxyz(2) = dy
cc        dxyz(3) = dz
cc        r = sqrt(dx**2 + dy**2 + dz**2)
cc        bbp(1:3,i) = dxyz(1:3)/r**3
cc      enddo

      
      call prin2('bbp=*',bbp,24)

      do i=1,npts
        rhs(i,1:3) = bbm(1:3,i)/dzk - bbp(1:3,i) 
        rhs(i,4:6) = jm(1:3,i)
      enddo

      numit = 100
      allocate(errs(numit+1))
      eps_gmres = 1.0d-6
      call statj_gendeb_solver(npatches,norders,ixyzs,iptype,
     1  npts,srccoefs,srcvals,eps,dzk,numit,rhs,eps_gmres,
     2  niter,errs,rres,soln)
      call prin2('rres=*',rres,1)

      rinttmp(1:6) = 0
      do i=1,npts
        rinttmp(1:6) = rinttmp(1:6) + soln(i,1:6)*wts(i)
      enddo
      call prin2('integral of densities for final solution=*',rinttmp,6)

c
c   test solution at an exterior point
c

      bbpc(1:3) = 0
      bbpex(1:3) = 0
      ra = 0

      do i=1,npts
        dx = xyz_out(1) - srcvals(1,i)
        dy = xyz_out(2) - srcvals(2,i)
        dz = xyz_out(3) - srcvals(3,i)
        
        r = sqrt(dx**2 + dy**2 + dz**2)

        bbpc(1) = bbpc(1) - dx/r**3*soln(i,5)*wts(i)
        bbpc(2) = bbpc(2) - dy/r**3*soln(i,5)*wts(i)
        bbpc(3) = bbpc(3) - dz/r**3*soln(i,5)*wts(i)
        ra = ra + wts(i)
      enddo
      bbpc(1:3) = bbpc(1:3)/4/pi
      print *, "ra=",ra

cc      call l3ddirectdg(1,xyz_in,df2,1,xyz_out,1,ptmp,
cc     1  bbpex,thresh)
      call l3ddirectcg(1,xyz_in,cf2,1,xyz_out,1,ptmp,
     1  bbpex,thresh)

      erra = abs(bbpc(1)-bbpex(1)) + abs(bbpc(2)-bbpex(2)) + 
     1   (bbpc(3)-bbpex(3))
      ra = abs(bbpex(1)) + abs(bbpex(2)) + abs(bbpex(3))
      call prin2('bbpex=*',bbpex,3)
      call prin2('bbpc=*',bbpc,3)
      print *, bbpex(1)/bbpc(1)
      print *, bbpex(2)/bbpc(2)
      print *, bbpex(3)/bbpc(3)
      call prin2('absolute error in exterior b field=*',erra/ra,1)


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



