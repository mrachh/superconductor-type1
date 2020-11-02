
      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: wts(:)
      character *100 fname
      integer ipars(2)

      integer, allocatable :: norders(:),ixyzs(:),iptype(:)
      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      integer, allocatable :: novers(:),ixyzso(:)
      integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)
      real *8, allocatable :: srcover(:,:),wover(:)

      real *8, allocatable :: hvec(:,:)
      real *8, allocatable :: rhs(:)
      real *8, allocatable :: hvecdiv(:)

      real *8, allocatable :: cms(:,:),rads(:),rad_near(:)

      real *8, allocatable :: wnear(:)
      real *8 xyz_out(3),xyz_in(3)
      complex *16 zpars


      call prini(6,13)

      done = 1
      pi = atan(done)*4

      igeomtype = 4

      ipars(1) = 8*2
      ipars(2) = 2*2

      npatches = 2*ipars(1)*ipars(2)
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

      ixyzs(npatches+1) = 1+npols*npatches
      allocate(wts(npts))
      call get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,wts)

      ra = 0
      do i=1,npts
        ra = ra + wts(i)
      enddo
      call prin2('surface area of torus=*',ra,1)

      allocate(hvec(3,npts),rhs(3*npts),hvecdiv(npts))
      do i=1,npts
        rr = srcvals(1,i)**2 + srcvals(2,i)**2
        hvec(1,i) = -srcvals(2,i)/rr
        hvec(2,i) = srcvals(1,i)/rr
        hvec(3,i) = 0.0d0
      enddo

      call prin2('hvec=*',hvec,2*3*npols)
      call prin2('srcvals=*',srcvals(1:3,1:npts),3*npols)

      call surf_div(npatches,norders,ixyzs,iptype,npts,srccoefs,
     1  srcvals,hvec,hvecdiv)
      
      ra = 0
      do i=1,npts
        ra = ra + hvecdiv(i)**2*wts(i)
      enddo
      ra = sqrt(ra)
      call prin2('error in surface divergence of hvec=*',ra,1)

      allocate(cms(3,npatches),rads(npatches),rad_near(npatches))
      call get_centroid_rads(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,cms,rads)
      
      call get_rfacs(norder,1,rfac,rfac0)
      call prin2('rfac=*',rfac,1)

      do i=1,npatches
        rad_near(i) = rads(i)*rfac
      enddo
      
      call findnearmem(cms,npatches,rad_near,12,srcvals,npts,nnz)
      call prinf('nnz=*',nnz,1)
      allocate(row_ptr(npts+1),col_ind(nnz),iquad(nnz+1))
      call findnear(cms,npatches,rad_near,12,srcvals,npts,row_ptr,
     1  col_ind)
      
      call get_iquad_rsc(npatches,ixyzs,npts,nnz,row_ptr,col_ind,iquad)
      nquad = iquad(nnz+1)-1


c
c  oversample geometry
c
      allocate(ixyzso(npatches+1),novers(npatches))
      eps = 0.51d-5
      ikerorder = 0
      zpars = 0.0d0
      call get_far_order(eps,npatches,norders,ixyzs,iptype,cms,
     1 rads,npts,srccoefs,12,npts,srcvals,ikerorder,zpars,nnz,
     2 row_ptr,col_ind,rfac,novers,ixyzso)
      call prinf('novers=*',novers,20)
      call prinf('ixyzso=*',ixyzso,20)

      nptso = ixyzso(npatches+1) - 1
      allocate(srcover(12,nptso),wover(nptso))
      call oversample_geom(npatches,norders,ixyzs,iptype,npts,
     1 srccoefs,srcvals,novers,ixyzso,nptso,srcover)

      call get_qwts(npatches,novers,ixyzso,iptype,npts,srcover,wover)


      dzk = 1.1d0
      iquadtype = 1
      allocate(wnear(10*nquad))
      wnear = 0
      call getnearquad_statj_gendeb(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,eps,dzk,iquadtype,nnz,row_ptr,col_ind,iquad,
     2  rfac0,nquad,wnear)

      hvecdiv = 0
      do i=1,npts
        rhs(i) = hvec(1,i)
        rhs(i+npts) = hvec(2,i)
        rhs(i+2*npts) = hvec(3,i)
      enddo

      call prinf('nptso=*',nptso,1)
      call lpcomp_divs0tan_addsub(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,eps,nnz,row_ptr,col_ind,iquad,nquad,
     2  wnear(2*nquad+1),rhs,novers,nptso,ixyzso,srcover,wover,
     3  hvecdiv) 
      
      ra = 0
      do i=1,npts
        ra = ra + hvecdiv(i)**2*wts(i)
      enddo
      ra = sqrt(ra)
      call prin2('error in surface divs0tan=*',ra,1)


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

        p2(1) = 1.2d0
        p2(2) = 1.0d0
        p2(3) = 1.7d0

c
c         numberof oscillations
c
        p4(1) = 5.0d0


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
        p1(2) = 3.0d0
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
          call h3d_sprime(xyzout,12,srcvals(1,i),0,dpars,1,zk,0,ipars,
     1       val)

          call cross_prod3d(srcvals(4,i),srcvals(7,i),tmp)
          ds = sqrt(tmp(1)**2 + tmp(2)**2 + tmp(3)**2)
          ra = ra + real(val)*wts(i)
        enddo
      enddo

      if(abs(ra+4*pi).le.1.0d-1) isout = .false.
      if(abs(ra).le.1.0d-1) isout = .true.

      return
      end

   




