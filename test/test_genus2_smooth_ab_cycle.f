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
      real *8, allocatable :: hvec1(:,:),hvec2(:,:)
      real *8, allocatable :: srcinterpa(:,:),srcinterpb(:,:)
      real *8, allocatable :: hinterp1a(:,:),hinterp2a(:,:)
      real *8, allocatable :: hinterp1b(:,:),hinterp2b(:,:)
      real *8, allocatable :: hinterp1aex(:,:),hinterp2aex(:,:)
      real *8, allocatable :: hinterp1bex(:,:),hinterp2bex(:,:)
      real *8, allocatable :: rrhs1(:),rrhs2(:)
      integer, allocatable :: apatches(:),bpatches(:)
      integer, allocatable :: iaxyzs(:),ibxyzs(:)
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


      norder = 2 
      npols = (norder+1)*(norder+2)/2
      fname = '../../fmm3dbie/geometries/genus_2_o08_r03.go3'
      call open_gov3_geometry_mem(fname,npatches,npts)
      allocate(srcvals(12,npts),srccoefs(9,npts),wts(npts))
      allocate(norders(npatches),ixyzs(npatches+1),iptype(npatches))
      call open_gov3_geometry(fname,npatches,norders,ixyzs,iptype,
     1   npts,srcvals,srccoefs,wts)

      call surf_vtk_plot(npatches,norders,ixyzs,iptype,npts,srccoefs,
     1   srcvals,'genus2.vtk','a')

      fname = '2torus_acycle_ref3.dat'
      open(unit=33,file=trim(fname))
      read(33,*) ngenus
      read(33,*) ne1
      read(33,*) ne2

      m = 20
      na = (ne1+ne2)*m
      allocate(avals(9,na),awts(na),apatches(na),auv(2,na))
      allocate(iaxyzs(3))
      close(33)
      call get_cycle_readfile(npatches,norders,ixyzs,iptype,npts,
     1   srccoefs,srcvals,ngenus,m,na,fname,avals,awts,apatches,auv,
     2   iaxyzs)
      
      na1 = iaxyzs(2)-iaxyzs(1)  
      call vtk_curv_plot(na1,9,avals,'acycle1.vtk','a')
      na2 = iaxyzs(3) - iaxyzs(2)
      call vtk_curv_plot(na2,9,avals(1,iaxyzs(2)),'acycle2.vtk','a')


      fname = '2torus_bcycle_ref3.dat'
      open(unit=33,file=trim(fname))
      read(33,*) ngenus
      read(33,*) ne1
      read(33,*) ne2

      m = 20
      nb = (ne1+ne2)*m
      allocate(bvals(9,nb),bwts(nb),bpatches(nb),buv(2,nb))
      allocate(ibxyzs(3))
      close(33)
      call get_cycle_readfile(npatches,norders,ixyzs,iptype,npts,
     1   srccoefs,srcvals,ngenus,m,nb,fname,bvals,bwts,bpatches,buv,
     2   ibxyzs)
      
      nb1 = ibxyzs(2)-ibxyzs(1)  
      call vtk_curv_plot(nb1,9,bvals,'bcycle1.vtk','a')
      nb2 = ibxyzs(3) - ibxyzs(2)
      call vtk_curv_plot(nb2,9,bvals(1,ibxyzs(2)),'bcycle2.vtk','a')
      stop


      rvecb = 0
      rscb = 0

      rveca = 0
      rsca = 0 

      do igen=1,ngenus
        rscb = 0
        rvecb = 0
        do i=ibxyzs(igen),ibxyzs(igen+1)-1
          rvecb = rvecb + (bvals(3,i)*bvals(4,i) + 
     1        bvals(1,i)*bvals(5,i) + bvals(2,i)*bvals(6,i))*bwts(i)
          ds = sqrt(bvals(4,i)**2 + bvals(5,i)**2 + bvals(6,i)**2)
          rscb = rscb + bwts(i)*ds
        enddo
        print *, igen,rscb
        print *, igen,rvecb
      enddo
      
      do igen=1,ngenus
        rsca = 0
        rveca = 0
        do i=iaxyzs(igen),iaxyzs(igen+1)-1
          rveca = rveca + (avals(3,i)*avals(4,i) + 
     1        avals(1,i)*avals(5,i) + avals(2,i)*avals(6,i))*awts(i)
          ds = sqrt(avals(4,i)**2 + avals(5,i)**2 + avals(6,i)**2)
          rsca = rsca + awts(i)*ds
        enddo
        print *, igen,rsca
        print *, igen,rveca
      enddo
      stop


      allocate(hvec1(3,npts),hvec2(3,npts))

      ra1 = 0
      rb1 = 0
      ra2 = 0
      rb2 = 0
      do i=1,npts
        rr1 = srcvals(1,i)**2 + srcvals(2,i)**2
        hvec1(1:3,i) = srcvals(4:6,i)/rr1
        call cross_prod3d(srcvals(10,i),hvec1(1,i),hvec2(1,i))
      enddo

      allocate(hinterp1a(3,na),hinterp2a(3,na))
      allocate(hinterp1aex(3,na),hinterp2aex(3,na))

      allocate(hinterp1b(3,nb),hinterp2b(3,nb))
      allocate(hinterp1bex(3,nb),hinterp2bex(3,nb))


      call fun_surf_interp(3,npatches,norders,ixyzs,iptype,npts,hvec1,
     1   na,apatches,auv,hinterp1a)
      call fun_surf_interp(3,npatches,norders,ixyzs,iptype,npts,hvec2,
     1   na,apatches,auv,hinterp2a)

      call fun_surf_interp(3,npatches,norders,ixyzs,iptype,npts,hvec1,
     1   nb,bpatches,buv,hinterp1b) 
      call fun_surf_interp(3,npatches,norders,ixyzs,iptype,npts,hvec2,
     1   nb,bpatches,buv,hinterp2b)

      allocate(srcinterpa(12,na),srcinterpb(12,nb))
      call geom_coefs_interp(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,na,apatches,auv,srcinterpa)

c
c  test normals
c
      ra = 0
      erra = 0
      do i=1,na
        ra = ra + srcinterpa(10,i)**2*awts(i)
        ra = ra + srcinterpa(11,i)**2*awts(i)
        ra = ra + srcinterpa(12,i)**2*awts(i)
        erra = erra + (srcinterpa(10,i)-avals(7,i))**2*awts(i)
        erra = erra + (srcinterpa(11,i)-avals(8,i))**2*awts(i)
        erra = erra + (srcinterpa(12,i)-avals(9,i))**2*awts(i)
      enddo


      erra = sqrt(erra/ra)
      call prin2('error in normals on acycles=*',erra,1)

cc      call prin2('srcinterpa=*',srcinterpa(1:3,1:12),36)
cc      call prin2('avals=*',avals(1:3,1:12),36)
      
      call geom_coefs_interp(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,nb,bpatches,buv,srcinterpb) 
c
c  test normals
c
      ra = 0
      erra = 0
      do i=1,nb
        ra = ra + srcinterpb(10,i)**2*bwts(i)
        ra = ra + srcinterpb(11,i)**2*bwts(i)
        ra = ra + srcinterpb(12,i)**2*bwts(i)
        erra = erra + (srcinterpb(10,i)-bvals(7,i))**2*bwts(i)
        erra = erra + (srcinterpb(11,i)-bvals(8,i))**2*bwts(i)
        erra = erra + (srcinterpb(12,i)-bvals(9,i))**2*bwts(i)
      enddo

      erra = sqrt(erra/ra)
      call prin2('error in normals on bcycles=*',erra,1)
      
c
c  compute exact values of the interpolant
c

      ra1 = 0
      erra1 = 0
      ra2 = 0
      erra2 = 0
      do i=1,na 
        rr1 = srcinterpa(1,i)**2 + srcinterpa(2,i)**2 
        hinterp1aex(1:3,i) = srcinterpa(4:6,i)/rr1
        call cross_prod3d(srcinterpa(10,i),hinterp1aex(1,i),
     1     hinterp2aex(1,i))
        ds = sqrt(avals(4,i)**2 + avals(5,i)**2 + avals(6,i)**2)
        ra1 = ra1 + (hinterp1aex(1,i)**2 + hinterp1aex(2,i)**2 + 
     1      hinterp1aex(3,i)**2)*awts(i)*ds
        erra1 = erra1 + (hinterp1aex(1,i)-hinterp1a(1,i))**2*awts(i)*ds
        erra1 = erra1 + (hinterp1aex(2,i)-hinterp1a(2,i))**2*awts(i)*ds
        erra1 = erra1 + (hinterp1aex(3,i)-hinterp1a(3,i))**2*awts(i)*ds

        ra2 = ra2 + (hinterp2aex(1,i)**2 + hinterp2aex(2,i)**2 + 
     1      hinterp2aex(3,i)**2)*awts(i)*ds
        erra2 = erra2 + (hinterp2aex(1,i)-hinterp2a(1,i))**2*awts(i)*ds
        erra2 = erra2 + (hinterp2aex(2,i)-hinterp2a(2,i))**2*awts(i)*ds
        erra2 = erra2 + (hinterp2aex(3,i)-hinterp2a(3,i))**2*awts(i)*ds
      enddo

      erra1 = sqrt(erra1/ra1)
      erra2 = sqrt(erra2/ra2)
      call prin2('error in interpolated density on acycle=*',erra1,1)
      call prin2('error in interpolated density2 on acycle=*',erra2,1)

      stop


      ra1 = 0
      erra1 = 0
      ra2 = 0
      erra2 = 0
      do i=1,nb 
        rr1 = srcinterpb(1,i)**2 + srcinterpb(2,i)**2 
        hinterp1bex(1:3,i) = srcinterpb(4:6,i)/rr1
        call cross_prod3d(srcinterpb(10,i),hinterp1bex(1,i),
     1     hinterp2bex(1,i))
        ds = sqrt(bvals(4,i)**2 + bvals(5,i)**2 + bvals(6,i)**2)
        ra1 = ra1 + (hinterp1bex(1,i)**2 + hinterp1bex(2,i)**2 + 
     1      hinterp1bex(3,i)**2)*bwts(i)*ds
        erra1 = erra1 + (hinterp1bex(1,i)-hinterp1b(1,i))**2*bwts(i)*ds
        erra1 = erra1 + (hinterp1bex(2,i)-hinterp1b(2,i))**2*bwts(i)*ds
        erra1 = erra1 + (hinterp1bex(3,i)-hinterp1b(3,i))**2*bwts(i)*ds

        ra2 = ra2 + (hinterp2bex(1,i)**2 + hinterp2bex(2,i)**2 + 
     1      hinterp2bex(3,i)**2)*bwts(i)*ds
        erra2 = erra2 + (hinterp2bex(1,i)-hinterp2b(1,i))**2*bwts(i)*ds
        erra2 = erra2 + (hinterp2bex(2,i)-hinterp2b(2,i))**2*bwts(i)*ds
        erra2 = erra2 + (hinterp2bex(3,i)-hinterp2b(3,i))**2*bwts(i)*ds
      enddo

      erra1 = sqrt(erra1/ra1)
      erra2 = sqrt(erra2/ra2)
      call prin2('error in interpolated density on bcycle=*',erra1,1)
      call prin2('error in interpolated density2 on bcycle=*',erra2,1)

c
c
c   compute all 4 integrals
c

      ra1 = 0
      ra2 = 0
      do i=1,na
        ra1 = ra1 + (hinterp1aex(1,i)*avals(4,i) +
     1     hinterp1aex(2,i)*avals(5,i) + hinterp1aex(3,i)*avals(6,i))*
     2     awts(i)
        ra2 = ra2 + (hinterp2aex(1,i)*avals(4,i) +
     1     hinterp2aex(2,i)*avals(5,i) + hinterp2aex(3,i)*avals(6,i))*
     2     awts(i)
      enddo

      call prin2('ra1=*',ra1,1)
      call prin2('ra2=*',ra2,1)

      rb1 = 0
      rb2 = 0
      do i=1,nb
        rb1 = rb1 + (hinterp1bex(1,i)*bvals(4,i) +
     1     hinterp1bex(2,i)*bvals(5,i) + hinterp1bex(3,i)*bvals(6,i))*
     2     bwts(i)
        rb2 = rb2 + (hinterp2bex(1,i)*bvals(4,i) +
     1     hinterp2bex(2,i)*bvals(5,i) + hinterp2bex(3,i)*bvals(6,i))*
     2     bwts(i)      
      enddo

      call prin2('rb1=*',rb1,1)
      call prin2('rb2=*',rb2,1)

c
c  compute integral of n\times first density on a cycle
c
c  
      ra = 0
      do i=1,na
        call cross_prod3d(srcinterpa(10,i),hinterp1aex(1,i),wtmp1)
        ra = ra + (wtmp1(1)*avals(4,i) +
     1     wtmp1(2)*avals(5,i)+wtmp1(3)*avals(6,i))*awts(i)
      enddo
      call prin2('integral of n cross density in vplus=*',ra,1)
      

      allocate(rrhs1(npts),rrhs2(npts))

      call surf_div(npatches,norders,ixyzs,iptype,npts,srccoefs,
     1  srcvals,hvec1,rrhs1)

      ra = 0
      erra = 0
      do i=1,npts
        ra = ra + abs(hvec1(1,i)**2 + hvec1(2,i)**2 +
     1     hvec1(3,i)**2)*wts(i)
        erra = erra + rrhs1(i)**2*wts(i)
      enddo
      erra = sqrt(erra/ra)
      ra = sqrt(ra)
      call prin2('error in surf div=*',erra,1)
      call prin2('integral of first harmonic field=*',ra,1)

      call surf_div(npatches,norders,ixyzs,iptype,npts,srccoefs,
     1  srcvals,hvec1,rrhs2)
      
      ra = 0
      erra = 0
      do i=1,npts
        ra = ra + abs(hvec2(1,i)**2 + hvec2(2,i)**2 +
     1     hvec2(3,i)**2)*wts(i)
        erra = erra + rrhs2(i)**2*wts(i)
      enddo
      erra = sqrt(erra/ra)
      ra = sqrt(ra)
      call prin2('error in n cross surf div=*',erra,1)
      call prin2('integral of second harmonic field=*',ra,1)


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






