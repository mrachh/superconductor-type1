      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: wts(:),rsigma(:),rpot(:)
      integer ipars(2)

      integer, allocatable :: norders(:),ixyzs(:),iptype(:)
      
      real *8 xyz_in_targ(3,100),xyz_out_targ(3,100)
      real *8 xyz_out_src(3),xyz_in_src(3)
      real *8 vtmp1(3),bbpc(3),bbpex(3)
      real *8 bbmc(3),bbmex(3),bjmc(3),bjmex(3)
      real *8, allocatable :: ffform(:,:,:),ffformex(:,:,:)
      real *8, allocatable :: ffforminv(:,:,:),ffformexinv(:,:,:)
      real *8, allocatable :: curv(:)

      complex *16, allocatable :: zbbm(:,:),zbbp(:,:),zjm(:,:)
      real *8, allocatable :: bbm(:,:),bbp(:,:),bjm(:,:)
      real *8, allocatable :: bbp_a(:,:),bbp_b(:,:)
      real *8, allocatable :: bbm_a(:,:),bbm_b(:,:)
      real *8, allocatable :: rhs(:),soln(:)
      real *8, allocatable :: errs(:)
      real *8, allocatable :: bbpcomp(:,:),bjmcomp(:,:),bbmcomp(:,:)
      real *8, allocatable :: bjmcomp_targ(:,:),bbmcomp_targ(:,:)
      real *8, allocatable :: bbpcomp_targ(:,:)
      real *8, allocatable :: bjm_targ(:,:),bbm_targ(:,:)
      real *8, allocatable :: bbp_targ(:,:)
      complex *16, allocatable :: zjm_targ(:,:),zbbm_targ(:,:)
      real *8 thet,phi,eps_gmres
      real *8, allocatable :: avals(:,:),bvals(:,:)
      real *8, allocatable :: awts(:),bwts(:)
      real *8, allocatable :: auv(:,:),buv(:,:)
      real *8, allocatable :: hvecs(:,:,:),bbphvecs(:,:,:)
      real *8, allocatable :: rhstmp(:),outtmp(:)
      real *8, allocatable :: hvecs_div(:,:),hvecs_div2(:)
      real *8, allocatable :: hvecs_a(:,:,:),bbphvecs_a(:,:,:)
      real *8, allocatable :: hvecs_b(:,:,:),bbphvecs_b(:,:,:)
      integer, allocatable :: apatches(:),bpatches(:)
      integer iaxyzs(2),ibxyzs(2)
      real *8 vf2(3,2),dpars(3),cf2(2),rzf2(3),dpars2(2)
      real *8 xyz_start(3), dxyz(3)
      complex *16 zf2(3),zk
      complex * 16 zpars(3)
      integer numit,niter,ndims(3)
      character *100 title,dirname
      character *1000 fname,fname1

      integer, allocatable :: ipatch_id(:),ipatch_id_targ(:)
      real *8, allocatable :: uvs_targ(:,:),uvs_src(:,:)
      real *8, allocatable :: wnear(:),wnear_targ(:)
      real *8, allocatable :: sources(:,:),targs(:,:)
      real *8, allocatable :: cms(:,:),rads(:),rad_near(:)
      integer, allocatable :: iquad(:),row_ptr(:),col_ind(:)
      integer, allocatable :: iquad_targ(:),row_ptr_targ(:)
      integer, allocatable :: col_ind_targ(:)
      integer, allocatable :: isout(:)


      integer, allocatable :: novers(:),ixyzso(:)
      real *8, allocatable :: srcover(:,:),wover(:)
      real *8, allocatable :: ptmp(:)
      real *8, allocatable :: laps02rhom(:),laps02rhop(:)
      real *8, allocatable :: laps02mum(:),blm0(:,:),bmm0(:,:)
      real *8, allocatable :: blm(:,:),bmm(:,:)
      real *8, allocatable :: wtmp1(:,:),wtmp2(:,:),wtmp3(:,:)
      real *8, allocatable :: wtmp4(:,:)
      real *8 rinttmp(6)

      logical isout0,isout1

      complex *16 pot,potex,ztmp,ima

      data ima/(0.0d0,1.0d0)/


      call prini(6,13)

      done = 1
      pi = atan(done)*4
c
c
c  ntarg test must be less than 100
c

      ntargtest = 10

      ibg = 1
c  
c  The penetration depth is 10^(-idzk) for idzk=0,1,2
c  but for idzk =3, penetration depth is 1/30
c
c  idzk = 4:10, 2^(4-idzk)
c
      idzk = 4

c
c  Note for this code, igeomtype can only be 1,2,3,4,6 
c
c  igeomtype = 1 => sphere
c  igeomtype = 2 => stellartor
c  igeomtype = 3 => wiggly torus
c  igeomtype = 4 => torus
c  igeomtype = 6 => thin shell tori
c
      igeomtype = 4
      iref = 2
      if(igeomtype.eq.1) then
        ipars(1) = iref
        npatches = 12*(4**(ipars(1)))
        fname = 'sphere.vtk'
        dirname = '/Users/mrachh/git/' // 
     1     'superconductor-type1-data/sphere-data/'
        ngenus = 0
        xyz_in_src(1) = 0.201d0
        xyz_in_src(2) = 0.102d0
        xyz_in_src(3) = 0.011d0

        xyz_out_src(1) = -1.2d0
        xyz_out_src(2) = 1.1d0
        xyz_out_src(3) = 1.7d0
      endif

      if(igeomtype.eq.2) then
        ipars(1) = 5*2**(iref)
        ipars(2) = 15*2**(iref)
        dirname = '/mnt/home/mrachh/ceph/' // 
     1     'superconductor-type1-data/stell-data/'
c        ipars(1) = 30
c        ipars(2) = 90
        npatches = 2*ipars(1)*ipars(2)

        xyz_in_src(1) = -4.503d0
        xyz_in_src(2) = -0.01d0
        xyz_in_src(3) = 0.001d0

        xyz_out_src(1) = -3.5d0
        xyz_out_src(2) = 7.7d0
        xyz_out_src(3) = 5.7d0


        fname = 'stell.vtk'
        ngenus = 1
      endif

      if(igeomtype.eq.3) then
        ipars(1) = 8*2**(iref)
        ipars(2) = 4*2**(iref)
        npatches = 2*ipars(1)*ipars(2)
        ngenus = 1
        
        fname = 'wtorus.vtk'
        dirname = '/mnt/home/mrachh/ceph/' // 
     1     'superconductor-type1-data/wtorus-data/'
        xyz_in_src(1) = -2.001d0
        xyz_in_src(2) = 0.002d0
        xyz_in_src(3) = 0.001d0

        uu = 0.91d0*pi
        vv = 0.17d0*2*pi
        rr = 1.33d0 
        xyz_out_src(1) = (rr*cos(uu) + 2+0.25d0*cos(3*vv))*cos(vv)
        xyz_out_src(2) = (rr*cos(uu) + 2+0.25d0*cos(3*vv))*sin(vv)
        xyz_out_src(3) = rr*sin(uu)

      endif

      if(igeomtype.eq.4) then
        ipars(1) = 4*2**(iref)
        ipars(2) = 2*2**(iref)
        npatches = 2*ipars(1)*ipars(2)

        dirname = '/Users/mrachh/git/' // 
     1     'superconductor-type1-data/torus-data/'
        fname = 'torus.vtk'
        ngenus = 1

        xyz_in_src(1) = 2.001d0
        xyz_in_src(2) = 0.002d0
        xyz_in_src(3) = 0.001d0

        uu = 0.75d0*2*pi
        vv = 0.22d0*2*pi
        rr = 1.37d0 
        xyz_out_src(1) = (rr*cos(uu) + 2)*cos(vv)
        xyz_out_src(2) = (rr*cos(uu) + 2)*sin(vv)
        xyz_out_src(3) = rr*sin(uu)
      endif

      if(igeomtype.eq.6) then
        ipars(1) = 4*2**(iref)
        ipars(2) = 2*2**(iref)
        npatches = 4*ipars(1)*ipars(2)

        dirname = '/Users/mrachh/git/' // 
     1     'superconductor-type1-data/thin-shell-torus-data/'
        fname = 'thin-torus.vtk'
        ngenus = 1

        xyz_in_src(1) = 2.001d0
        xyz_in_src(2) = 0.002d0
        xyz_in_src(3) = 0.001d0

        uu = 0.75d0*2*pi
        vv = 0.22d0*2*pi
        rr = 1.37d0 
        xyz_out_src(1) = (rr*cos(uu) + 2)*cos(vv)
        xyz_out_src(2) = (rr*cos(uu) + 2)*sin(vv)
        xyz_out_src(3) = rr*sin(uu)
      endif

      norder = 6
      npols = (norder+1)*(norder+2)/2

      npts = npatches*npols
      call prinf('npatches=*',npatches,1)

      allocate(srcvals(12,npts),srccoefs(9,npts))
      ifplot = 1


      call setup_geom(igeomtype,norder,npatches,ipars, 
     1       srcvals,srccoefs,ifplot,fname)
cc      call prin2('srcvals=*',srcvals,12*npts)
cc      call prin2('srccoefs=*',srccoefs,9*npts)

      allocate(norders(npatches),ixyzs(npatches+1),iptype(npatches))

      do i=1,npatches
        norders(i) = norder
        ixyzs(i) = 1 +(i-1)*npols
        iptype(i) = 1
      enddo

      ixyzs(npatches+1) = 1+npols*npatches

      print *, "npatches=",npatches
      print *, "npts=",npts
      call prin2('srcvals=*',srcvals(1,npts/2+1), 24)
      call prin2('srccoefs=*',srccoefs(1,npts/2+1), 24)

      call surf_vtk_plot_vec(npatches, norders, ixyzs, iptype, 
     1  npts, srccoefs, srcvals, srcvals(10:12,:),
     2   'thin-shell-normals.vtk','a')


      allocate(wts(npts))

      call get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,wts)


      allocate(cms(3,npatches),rads(npatches),rad_near(npatches))

      call get_centroid_rads(npatches,norders,ixyzs,iptype,npts, 
     1     srccoefs,cms,rads)
      call prin2('cms=*',cms,24)
      call prin2('rads=*',rads,12)

      if(igeomtype.eq.1) then
        call get_sphere_testtarg(ntargtest,xyz_in_targ,xyz_out_targ)
      endif

      if(igeomtype.eq.2) then
       htest = 0.5d0*2**(iref) 
       call get_stell_testtarg(ntargtest,xyz_in_targ,xyz_out_targ,
     1    npatches,npts,cms,rads,srcvals,wts,htest)
      endif

      if(igeomtype.eq.3) then
        call get_wtorus_testtarg(ntargtest,xyz_in_targ,xyz_out_targ)
      endif

      if(igeomtype.eq.4) then
        call get_torus_testtarg(ntargtest,xyz_in_targ,xyz_out_targ)
      endif

      call prin2('xyz_in_targ=*',xyz_in_targ,3*ntargtest)
      call prin2('xyz_out_targ=*',xyz_out_targ,3*ntargtest)


      print *, " "
      print *, " "
      print *, "Testing interior targets"
      do i=1,ntargtest
        isout0 = .false.
        call test_exterior_pt(npts,srcvals,wts,xyz_in_targ(1,i),isout0)
          print *, i,isout0
      enddo
      print *, "=============="
      print *, " "
      print *, " "
      print *, "Testing exterior targets"

      do i=1,ntargtest
        isout0 = .false.
        call test_exterior_pt(npts,srcvals,wts,xyz_out_targ(1,i),isout0)
          print *, i,isout0
      enddo
      print *, "=============="
c
c  set file names for reading or writing harmonic vector fields
c
c c
c       write(fname,'(a,a,i3.3,a,i3.3,a,i1,a)') trim(dirname),
c      1     'hvecs_',ipars(1),
c      1    '_',ipars(2),'_',norder,'_1.dat'
c       if(igeomtype.ge.2)
c      1   open(unit=78,file=trim(fname),form='unformatted')
c       write(fname,'(a,a,i3.3,a,i3.3,a,i1,a)') trim(dirname),
c      1     'hvecs_',ipars(1),
c      1     '_',ipars(2),'_',norder,'_2.dat'
c       print *, fname 
c       if(igeomtype.ge.2) 
c      1   open(unit=79,file=trim(fname),form='unformatted')

c
c   set a and b cycle params
c

      if(igeomtype.ge.2) then
        m = 40
        na = ipars(2)*m
        nb = ipars(1)*m
        allocate(avals(9,na),awts(na),auv(2,na),apatches(na))
        allocate(bvals(9,nb),bwts(nb),buv(2,nb),bpatches(nb))
        allocate(bbphvecs_a(3,na,2),bbphvecs_b(3,nb,2))
        allocate(hvecs_a(3,na,2),hvecs_b(3,nb,2))
        call get_ab_cycles_torusparam(npatches,norders,ixyzs,iptype,
     1     npts,srccoefs,srcvals,ipars,m,na,avals,awts,apatches,auv,
     2     nb,bvals,bwts,bpatches,buv)
      endif

      if(igeomtype.eq.1) then
        na = 0
        nb = 0
        ntmp = 10
        allocate(avals(9,ntmp),awts(ntmp),auv(2,ntmp),apatches(ntmp))
        allocate(bvals(9,ntmp),bwts(ntmp),buv(2,ntmp),bpatches(ntmp))
        allocate(bbphvecs_a(3,ntmp,2),bbphvecs_b(3,ntmp,2))
        allocate(hvecs_a(3,ntmp,2),hvecs_b(3,ntmp,2))
        avals = 0
        awts = 0
        auv = 0
        apatches = 0
        bvals = 0
        bwts = 0
        buv = 0
        bpatches = 0
      endif
c
c
c   compute harmonic vector fields
c
     
      print *, "starting harmonic vector field computation"


      allocate(hvecs(3,npts,2),bbphvecs(3,npts,2),hvecs_div(npts,2))
      allocate(hvecs_div2(npts))
      allocate(rhstmp(npts*3),outtmp(npts*3))
      bbphvecs = 0

      if(igeomtype.eq.4.or.igeomtype.eq.5) then
        do i=1,npts
          rr1 = srcvals(1,i)**2 + srcvals(2,i)**2
          hvecs(1,i,1) = -srcvals(2,i)/rr1
          hvecs(2,i,1) = srcvals(1,i)/rr1
          hvecs(3,i,1) = 0 
          call cross_prod3d(srcvals(10,i),hvecs(1:3,i,1),
     1       hvecs(1:3,i,2))
        enddo
        call surf_div(npatches,norders,ixyzs,iptype,npts, 
     1   srccoefs,srcvals,hvecs(1,1,1),hvecs_div2)
        errest = 0
        do i=1,npts
          errest = errest + hvecs_div2(i)**2*wts(i)
        enddo
        errest = sqrt(errest)
        call prin2('errest=*',errest,1)
      endif

      if(igeomtype.eq.2.or.igeomtype.eq.3) then

        ifread = 1
        ifwrite = 0 
        if(ifread.eq.0) then
          eps = 0.51d-7
          call cpu_time(tt1)
C$         tt1 = omp_get_wtime()         
          print *, "Here1"
          call get_harm_vec_field(npatches,norders,ixyzs,iptype, 
     1      npts,srccoefs,srcvals,wts,eps,hvecs(1,1,1),errest)
          call cpu_time(tt2)
C$         tt2 = omp_get_wtime()          
          call prin2('errest=*',errest,1)
          print *, "here2"
          print *, "time taken = ",tt2-tt1
          do i=1,npts
            call cross_prod3d(srcvals(10,i),hvecs(1,i,1),hvecs(1,i,2))
          enddo
        else
          print *, "here1"
          read(78) hvecs(1:3,1:npts,1)
          read(79) hvecs(1:3,1:npts,2)

          ra1 = 0
          ra2 = 0
          do i=1,npts
            ra1 = ra1 + (hvecs(1,i,1)**2 + hvecs(2,i,1)**2 + 
     1         hvecs(3,i,1)**2)*wts(i)
            ra2 = ra2 + (hvecs(1,i,2)**2 + hvecs(2,i,2)**2 + 
     1         hvecs(3,i,2)**2)*wts(i)
          enddo
          ra1 = sqrt(ra1)
          ra2 = sqrt(ra2)
          do i=1,npts
            hvecs(1:3,i,1) = hvecs(1:3,i,1)/ra1
            hvecs(1:3,i,2) = hvecs(1:3,i,2)/ra2
          enddo
          ra1 = 0
          ra2 = 0
          do i=1,npts
            ra1 = ra1 + (hvecs(1,i,1)**2 + hvecs(2,i,1)**2 + 
     1         hvecs(3,i,1)**2)*wts(i)
            ra2 = ra2 + (hvecs(1,i,2)**2 + hvecs(2,i,2)**2 + 
     1         hvecs(3,i,2)**2)*wts(i)
          enddo
          ra1 = sqrt(ra1)
          ra2 = sqrt(ra2)
          call prin2('ra1=*',ra1,1)
          call prin2('ra2=*',ra2,1)
          
          close(78)
          close(79)
          call surf_div(npatches,norders,ixyzs,iptype,npts, 
     1     srccoefs,srcvals,hvecs(1,1,1),hvecs_div2)
          errest = 0
          do i=1,npts
            errest = errest + hvecs_div2(i)**2*wts(i)
          enddo
          errest = sqrt(errest)
          call prin2('errest=*',errest,1)
        endif
        if(ifwrite.eq.1) then
          write(78) hvecs(1:3,1:npts,1)
          write(79) hvecs(1:3,1:npts,2)
          close(78)
          close(79)
        endif
      endif


      iaxyzs(1) = 1
      iaxyzs(2) = na+1

      ibxyzs(1) = 1
      ibxyzs(2) = nb+1

c
c
      vf2(1:3,1) = 1.0d0 
      cf2(1) = 1*1.0d0

      thet1 = 0.7d0*2*pi
      thet2 = 0.33d0*2*pi
      thet3 = 0.6d0*2*pi
      zf2(1) = cos(thet1) 
      zf2(2) = cos(thet2) 
      zf2(3) = cos(thet3)

      ra = sqrt(abs(zf2(1))**2 + abs(zf2(2))**2 + abs(zf2(3))**2)
      zf2(1:3) = zf2(1:3)/ra
      rzf2(1:3) = real(zf2(1:3))


      thresh = 1.0d-16
      allocate(bbp(3,npts),ptmp(npts))
      allocate(bjm(3,npts),bbm(3,npts))
      allocate(zjm(3,npts),zbbm(3,npts))
      ptmp = 0
      bbp = 0
      bjm = 0
      bbm = 0

      ntarg = npts
      allocate(sources(3,npts))
C$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,npts 
        sources(1,i) = srcvals(1,i)
        sources(2,i) = srcvals(2,i)
        sources(3,i) = srcvals(3,i)
      enddo
C$OMP END PARALLEL DO      

      call prin2('thresh=*',thresh,1)
      call prin2('xyz_in_src=*',xyz_in_src,3)
      call prin2('cf2=*',cf2,1)
      call prin2('vf2=*',vf2,3)
      call l3ddirectcdg(1,xyz_in_src,cf2,vf2,1,sources,npts,
     1 ptmp,bbp,thresh)
      call prin2('bbp=*',bbp,24)
      dzk = 10**(idzk)

      if(idzk.eq.3) dzk = 30.0d0

      if(idzk.ge.4) dzk = 2**(idzk-4)

      zk = ima*dzk
      init = 0

      call fieldsED(zk,xyz_out_src,srcvals,npts,zjm,zbbm,zf2,init)
      bjm = real(zjm)
      bbm = real(zbbm)
      call prin2('bjm=*',bjm,24)
      call prin2('bbm=*',bbm,24)
      
      call prin2('ptmp=*',ptmp,24)
      allocate(rhs(6*npts+4*ngenus),soln(6*npts+4*ngenus))
      rhs = 0
      soln = 0
      do i=1,npts
        rhs(i) = bbm(1,i)/dzk-bbp(1,i)
        rhs(i+npts) = bbm(2,i)/dzk-bbp(2,i)
        rhs(i+2*npts) = bbm(3,i)/dzk-bbp(3,i)
        rhs(i+3*npts) = bjm(1,i)
        rhs(i+4*npts) = bjm(2,i)
        rhs(i+5*npts) = bjm(3,i)
      enddo

      if(igeomtype.ge.2) then

        allocate(bbp_a(3,na),bbp_b(3,nb),bbm_a(3,na),bbm_b(3,nb))
        call fun_surf_interp(3,npatches,norders,ixyzs,iptype,npts,
     1     bbp,na,apatches,auv,bbp_a)
        call fun_surf_interp(3,npatches,norders,ixyzs,iptype,npts,
     1     bbp,nb,bpatches,buv,bbp_b)

        call fun_surf_interp(3,npatches,norders,ixyzs,iptype,npts,
     1     bbm,na,apatches,auv,bbm_a)
        call fun_surf_interp(3,npatches,norders,ixyzs,iptype,npts,
     1     bbm,nb,bpatches,buv,bbm_b)

        rra = 0 
        do i=1,na
          rhs(6*npts+1) = rhs(6*npts+1) + (bbm_a(1,i)*avals(4,i) + 
     1      bbm_a(2,i)*avals(5,i) + bbm_a(3,i)*avals(6,i))*awts(i)
          rhs(6*npts+3) = rhs(6*npts+3) + (bbp_a(1,i)*avals(4,i) + 
     1      bbp_a(2,i)*avals(5,i) + bbp_a(3,i)*avals(6,i))*awts(i)
          rra = rra + 
     1       sqrt(avals(4,i)**2 + avals(5,i)**2 + avals(6,i)**2)*awts(i)
        enddo

        rrb = 0
        do i=1,nb
          rhs(6*npts+2) = rhs(6*npts+2) + (bbm_b(1,i)*bvals(4,i) + 
     1      bbm_b(2,i)*bvals(5,i) + bbm_b(3,i)*bvals(6,i))*bwts(i)
          rhs(6*npts+4) = rhs(6*npts+4) + (bbp_b(1,i)*bvals(4,i) + 
     1      bbp_b(2,i)*bvals(5,i) + bbp_b(3,i)*bvals(6,i))*bwts(i)
          rrb = rrb + 
     1      sqrt(bvals(4,i)**2 + bvals(5,i)**2 + bvals(6,i)**2)*bwts(i)
        enddo
      endif




      eps = 0.51d-8

      allocate(ipatch_id(npts),uvs_src(2,npts))
      call get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts,
     1  ipatch_id,uvs_src)
      call prinf('ipatch_id=*',ipatch_id,20)
      call prin2('uvs_src=*',uvs_src,48)
      
c
c       precompute near quadrature correction
c
c
      iptype_avg = floor(sum(iptype)/(npatches+0.0d0))
      norder_avg = floor(sum(norders)/(npatches+0.0d0))

      call get_rfacs(norder_avg,iptype_avg,rfac,rfac0)
      call prin2('rfac=*',rfac,1)
      call prin2('rfac0=*',rfac0,1)

C$OMP PARALLEL DO DEFAULT(SHARED) 
      do i=1,npatches
        rad_near(i) = rads(i)*rfac
      enddo
C$OMP END PARALLEL DO     

      call prin2('rad_near=*',rad_near,12)
c
c    find near quadrature correction interactions
c
      call findnearmem(cms,npatches,rad_near,3,sources,npts,nnz)

      allocate(row_ptr(npts+1),col_ind(nnz))
      
      call findnear(cms,npatches,rad_near,3,sources,npts,row_ptr, 
     1        col_ind)

      allocate(iquad(nnz+1)) 
      call get_iquad_rsc(npatches,ixyzs,ntarg,nnz,row_ptr,col_ind,
     1         iquad)

c
c    estimate oversampling for far-field, and oversample geometry
c

      ikerorder = 0
      allocate(novers(npatches),ixyzso(npatches+1))

      zpars = 0
      ndtarg = 3
      call get_far_order(eps,npatches,norders,ixyzs,iptype,cms,
     1    rads,npts,srccoefs,ndtarg,npts,sources,ikerorder,zpars,
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

      if(ibg.eq.1) then
        rbeta = 0.0d0
        rgamma = 0.0d0
      endif

      if(ibg.eq.2) then
        rbeta = 1.0d0
        rgamma = 0.0d0
      endif

      if(ibg.eq.3) then
        rbeta = 0.0d0
        rgamma = -1.0d0*dzk
      endif

      if(ibg.eq.4) then
        rbeta = 1.0d0
        rgamma = 1.0d0
      endif

      dpars(1) = dzk
      dpars(2) = rbeta
      dpars(3) = rgamma
      wnear = 0

c      goto 1111
      epsquad = 0.51d-7
      call cpu_time(t1)
C$      t1 = omp_get_wtime()      
      call getnearquad_statj_gendeb(npatches,norders,
     1      ixyzs,iptype,npts,srccoefs,srcvals,
     1      epsquad,dzk,iquadtype,nnz,row_ptr,col_ind,
     1      iquad,rfac0,nquad,wnear)
 1111 continue     
      call prinf('finished generating near quadrature correction*',i,0)
      call cpu_time(t2)
C$      t2 = omp_get_wtime()      
      call prin2('time taken to generate near quadrature=*',t2-t1,1)

      call prinf('entering layer potential eval*',i,0)
      call prinf('npts=*',npts,1)

c
c
c  compute bbphvecs
c
      hvecs_div = 0
      do igen=1,2*ngenus
        do i=1,npts
          rhstmp(i) = hvecs(1,i,igen) 
          rhstmp(i+npts) = hvecs(2,i,igen) 
          rhstmp(i+2*npts) = hvecs(3,i,igen) 
        enddo
        call prin2('rhstmp=*',rhstmp,24)
        outtmp = 0

        call lpcomp_s0curl_addsub(npatches,norders,ixyzs,iptype,npts,
     1    srccoefs,srcvals,eps,nnz,row_ptr,col_ind,iquad,nquad,
     2    wnear(2*nquad+1),rhstmp,novers,npts_over,ixyzso,
     2    srcover,wover,outtmp)
        do i=1,npts
          bbphvecs(1,i,igen) = outtmp(i)
          bbphvecs(2,i,igen) = outtmp(i+npts)
          bbphvecs(3,i,igen) = outtmp(i+2*npts)
        enddo
        hvecs_div = 0
        call lpcomp_divs0tan_addsub(npatches,norders,ixyzs,iptype,
     1   npts,srccoefs,srcvals,eps,nnz,row_ptr,col_ind,iquad,nquad,
     2   wnear(2*nquad+1),rhstmp,novers,npts_over,ixyzso,
     3   srcover,wover,hvecs_div(1:npts,igen))
        do j=1,npts
          call cross_prod3d(srcvals(10,j),hvecs(1,j,igen),vtmp1)
          bbphvecs(1:3,j,igen) = bbphvecs(1:3,j,igen) - vtmp1(1:3)/2
        enddo
      enddo
      print *, "here"

      numit = 300
      allocate(errs(numit+1))
      call prinf('ngenus=*',ngenus,1)
      eps_gmres = 1.0d-8
      if(ngenus.ge.1)  call prin2('rhs_projs=*',rhs(6*npts+1),4)
      print *, "here"


      do i=1,npts 
        write (776,'(3(2x,e11.5))') bbphvecs(1,i,1),
     1      bbphvecs(2,i,1), bbphvecs(3,i,1)
      enddo

      do i=1,npts
        write (777,'(3(2x,e11.5))') bbphvecs(1,i,2),
     1      bbphvecs(2,i,2), bbphvecs(3,i,2)
      enddo


      call cpu_time(t1)
C$       t1 = omp_get_wtime()      

      call statj_gendeb_solver(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,eps,dpars,ngenus,hvecs,bbphvecs,na,iaxyzs,
     2  apatches,auv,avals,awts,nb,ibxyzs,bpatches,buv,bvals,bwts,
     3  numit,rhs,eps_gmres,niter,errs,rres,soln)
      call cpu_time(t2)
C$       t2 = omp_get_wtime()      
      
      call prin2('solve time=*',t2-t1,1)
      
      
      call prin2('projs=*',soln(6*npts+1),4)



      allocate(bbpcomp(3,npts),bjmcomp(3,npts),bbmcomp(3,npts))

      print *, "here"
      call prin2('eps=*',eps,1)
      call prin2('soln=*',soln,24)
      call prinf('ngenus=*',ngenus,1)
      call prin2('dpars=*',dpars,3)


      call lpcomp_statj_gendeb_postproc(npatches,norders,ixyzs,
     1  iptype,npts,srccoefs,srcvals,eps,dpars,nnz,row_ptr,col_ind,
     2  iquad,nquad,wnear,ngenus,hvecs,bbphvecs,soln,novers,npts_over,
     3  ixyzso,srcover,wover,bjmcomp,bbmcomp,bbpcomp)
      call prin2('bjmcomp=*',bjmcomp,24)
      call prin2('bbmcomp=*',bbmcomp,24)
      call prin2('bbpcomp=*',bbpcomp,24)

      errj = 0
      errbm = 0
      errbp = 0
      rsurfintl2 = 0
      ra = 0
      ra2 = 0
      do i=1,npts
        do j=1,3
          errj  = errj + (bjmcomp(j,i) - bjm(j,i))**2*wts(i)
          errbm = errbm + (bbmcomp(j,i) - bbm(j,i))**2*wts(i)
          errbp = errbp + (bbpcomp(j,i) - bbp(j,i))**2*wts(i)
          ra2 = ra2 + bbp(j,i)**2*wts(i)
          ra2 = ra2 + bjm(j,i)**2*wts(i)
          ra2 = ra2 + bbm(j,i)**2*wts(i)
        enddo
        do j=1,6
          ra = ra + rhs((j-1)*npts+i)**2*wts(i)
        enddo
        do j=4,6
          rsurfintl2 = rsurfintl2 + soln((j-1)*npts+i)**2*wts(i)
        enddo
      enddo
      rsurfintl2 = sqrt(rsurfintl2)
      errbdry = sqrt((errj + errbm + errbp))/rsurfintl2



      errj = sqrt(errj)/rsurfintl2
      errbm = sqrt(errbm)/rsurfintl2
      errbp = sqrt(errbp)/rsurfintl2
      call prin2('error in current=*',errj,1)
      call prin2('error in interior magnetic field=*',errbm,1)
      call prin2('error in exterior magnetic field=*',errbp,1)
      call prin2('rsurfintl2 =*',rsurfintl2,1)


c
c  test solution at interior and exterior points  
c
      rinttmp = 0
      rsurf = 0
      do i=1,npts
        rsurf = rsurf + wts(i)
        do j=1,6
          rinttmp(j) = rinttmp(j) + soln(i+(j-1)*npts)*wts(i)
        enddo
      enddo
      rinttmp(1:6) = rinttmp(1:6)/rsurf
      call prin2('average integral of densities=*',rinttmp,6)
      
      allocate(curv(npts),ffforminv(2,2,npts))
      call get_mean_curvature(npatches,norders,ixyzs,iptype,npts, 
     1   srccoefs,srcvals,curv)
      call get_inv_first_fundamental_form(npatches,norders,ixyzs, 
     1   iptype,npts,srccoefs,srcvals,ffforminv)
      
      allocate(wtmp1(3,npts),wtmp2(3,npts),wtmp3(3,npts),wtmp4(3,npts))

      do i=1,npts
         wtmp1(1:3,i) = ffforminv(1,1,i)*srcvals(4:6,i) + 
     1      ffforminv(1,2,i)*srcvals(7:9,i)
         wtmp2(1:3,i) = ffforminv(2,1,i)*srcvals(4:6,i) + 
     1      ffforminv(2,2,i)*srcvals(7:9,i)
         call cross_prod3d(srcvals(10,i),wtmp1(1:3,i),wtmp3(1:3,i))
         call cross_prod3d(srcvals(10,i),wtmp2(1:3,i),wtmp4(1:3,i))
      enddo

      allocate(laps02rhom(npts),laps02rhop(npts),laps02mum(npts))
      allocate(blm0(3,npts),bmm0(3,npts),bmm(3,npts),blm(3,npts))

      
      call statj_gendebproc_rhomrhopmum(npatches,norders,ixyzs, 
     1  iptype,npts,srccoefs,srcvals,eps,nnz,row_ptr,col_ind,iquad, 
     2  nquad,wnear,soln,novers,npts_over,ixyzso,srcover,wover,curv, 
     3  wtmp1,wtmp2,wtmp3,wtmp4,dzk,rbeta,rgamma,laps02rhom,laps02rhop, 
     4  laps02mum,blm,bmm)

c
c  add in contribution of harmonic vector fields
c
      do i=1,npts
        do igen=1,2*ngenus
          blm(1:3,i) = blm(1:3,i) + soln(6*npts+igen)*hvecs(1:3,i,igen)
        enddo
      enddo

c
c  error at exterior targets
c
      erra = 0
      ra = 0
      do j=1,ntargtest
        bbpc(1:3) = 0
        bbpex(1:3) = 0
        ptmp = 0
        call l3ddirectcdg(1,xyz_in_src,cf2,vf2,1,xyz_out_targ(1,j),
     1     1,ptmp,bbpex,thresh)

        do i=1,npts
          dx = xyz_out_targ(1,j) - srcvals(1,i)
          dy = xyz_out_targ(2,j) - srcvals(2,i)
          dz = xyz_out_targ(3,j) - srcvals(3,i)
        
          sig = soln(4*npts+i)
        
          r = sqrt(dx**2 + dy**2 + dz**2)

          bbpc(1) = bbpc(1) - dx/r**3*sig*wts(i)
          bbpc(2) = bbpc(2) - dy/r**3*sig*wts(i)
          bbpc(3) = bbpc(3) - dz/r**3*sig*wts(i)

          do igen=1,2*ngenus
        
            bbpc(1) = bbpc(1) +1.0d0/r**3*wts(i)*soln(6*npts+2+igen)*
     1        (dy*hvecs(3,i,igen)-dz*hvecs(2,i,igen))
            bbpc(2) = bbpc(2) +1.0d0/r**3*wts(i)*soln(6*npts+2+igen)*
     1        (dz*hvecs(1,i,igen)-dx*hvecs(3,i,igen))
            bbpc(3) = bbpc(3) +1.0d0/r**3*wts(i)*soln(6*npts+2+igen)*
     1        (dx*hvecs(2,i,igen)-dy*hvecs(1,i,igen))
          enddo
        enddo
        bbpc(1:3) = bbpc(1:3)/4/pi

        erra = erra + (bbpc(1)-bbpex(1))**2 + (bbpc(2)-bbpex(2))**2 + 
     1   (bbpc(3)-bbpex(3))**2
        ra = ra + abs(bbpex(1))**2 + abs(bbpex(2))**2 + abs(bbpex(3))**2
      enddo


      errextsq = erra
      rextsq = ra
      errpt = sqrt(erra)/rsurfintl2
      call prin2('relativd error in exterior b field=*',errpt,1)

c
c  error in interior targets
c

      call prin2('rzf2=*',rzf2,3)
      call prin2('xyz_in_targ=*',xyz_in_targ,3)
      do j=1,ntargtest
        bbmc(1:3) = 0
        bbmex(1:3) = 0
        bjmc(1:3) = 0
        bjmex(1:3) = 0
        ptmp = 0
        dx = xyz_in_targ(1,j)-xyz_out_src(1)
        dy = xyz_in_targ(2,j)-xyz_out_src(2) 
        dz = xyz_in_targ(3,j)-xyz_out_src(3) 
        r = sqrt(dx**2 + dy**2 + dz**2)
        r1 = (-dzk/r**2-1/r**3)
        r2 = (dzk**2/r**3 + 3*dzk/r**4 + 3/r**5)
        au1 = (-dzk**2/r+r1)*exp(-dzk*r)/dzk
        au2 = dx*rzf2(1) + dy*rzf2(2) + dz*rzf2(3)
        au2 = au2*r2*exp(-dzk*r)/dzk
        bjmex(1) = rzf2(1)*au1 + dx*au2
        bjmex(2) = rzf2(2)*au1 + dy*au2
        bjmex(3) = rzf2(3)*au1 + dz*au2

        bbmex(1) = (rzf2(3)*dy - rzf2(2)*dz)*r1*exp(-dzk*r) 
        bbmex(2) = (rzf2(1)*dz - rzf2(3)*dx)*r1*exp(-dzk*r) 
        bbmex(3) = (rzf2(2)*dx - rzf2(1)*dy)*r1*exp(-dzk*r) 
        

        do i=1,npts
          dx = xyz_in_targ(1,j) - srcvals(1,i)
          dy = xyz_in_targ(2,j) - srcvals(2,i)
          dz = xyz_in_targ(3,j) - srcvals(3,i)
        
          r = sqrt(dx**2 + dy**2 + dz**2)
          r1 = (-dzk/r**2-1/r**3)
          gr1 = exp(-dzk*r)/r
          gr2 = r1*exp(-dzk*r)

          bjmc(1) = bjmc(1) + (-dzk*gr1*blm(1,i) -
     1       gr2*dx*soln(5*npts+i) -
     2       (bmm(3,i)*dy-bmm(2,i)*dz)*gr2)*wts(i)
          bjmc(2) = bjmc(2) + (-dzk*gr1*blm(2,i) -
     1       gr2*dy*soln(5*npts+i) -
     2       (bmm(1,i)*dz-bmm(3,i)*dx)*gr2)*wts(i)
          bjmc(3) = bjmc(3) + (-dzk*gr1*blm(3,i) -
     1       gr2*dz*soln(5*npts+i) -
     2       (bmm(2,i)*dx-bmm(1,i)*dy)*gr2)*wts(i)

          bbmc(1) = bbmc(1) + (-dzk*gr1*bmm(1,i) +
     1       gr2*dx*soln(3*npts+i) +
     2       (blm(3,i)*dy-blm(2,i)*dz)*gr2)*wts(i)
          bbmc(2) = bbmc(2) + (-dzk*gr1*bmm(2,i) +
     1       gr2*dy*soln(3*npts+i) +
     2       (blm(1,i)*dz-blm(3,i)*dx)*gr2)*wts(i)
          bbmc(3) = bbmc(3) + (-dzk*gr1*bmm(3,i) +
     1       gr2*dz*soln(3*npts+i) +
     2       (blm(2,i)*dx-blm(1,i)*dy)*gr2)*wts(i)
        enddo
        bbmc(1:3) = bbmc(1:3)/4/pi
        bjmc(1:3) = bjmc(1:3)/4/pi
        erra = erra + abs(bbmc(1)-bbmex(1))**2 + 
     1   abs(bbmc(2)-bbmex(2))**2 + 
     1   (bbmc(3)-bbmex(3))**2
        ra = ra + abs(bbmex(1))**2 + abs(bbmex(2))**2 + abs(bbmex(3))**2

        erra = erra + (bjmc(1)-bjmex(1))**2 + (bjmc(2)-bjmex(2))**2 + 
     1   (bjmc(3)-bjmex(3))**2
        ra = ra + abs(bjmex(1))**2 + abs(bjmex(2))**2 + abs(bjmex(3))**2
      enddo
      call prin2('erra=*',erra,1)
      call prin2('ra=*',ra,1)
      errvol_abs = sqrt(erra)
      errvol = sqrt(erra)/rsurfintl2
      call prin2('l2 rel error at interior + exterior targets=*',
     1  errvol,1)
      call prin2('l2 abs error at interior + exterior targets=*',
     1  errvol_abs,1)
      stop

      open(unit=81,file='res_mar17_2022.txt',access='append')
c
c  igeomtype,norder,npatches,niter,eps,eps_gmres,errbdry,errvol
c  

 1211 format(2x,i1,2x,i2,2x,i6,2x,i3,2(2x,i2),4(2x,e11.5)) 
      write(81,1211) igeomtype,norder,npatches,niter,ibg,idzk,
     1   eps,eps_gmres,errbdry,errvol
      close(81)

       
      if(igeomtype.eq.1) then
        write(fname,'(a,a,i2.2,a,i2.2,a,i1,a,i1,a,i1,a)') 
     1    trim(dirname),'analytic_soln_',ipars(1),'_',ipars(1),
     2    '_norder',norder,'_ibg',ibg,'_idzk',idzk,'.dat'
      else
        write(fname,'(a,a,i2.2,a,i2.2,a,i1,a,i1,a,i1,a)') 
     1    trim(dirname),'analytic_soln_',ipars(1),'_',ipars(2),
     2    '_norder',norder,'_ibg',ibg,'_idzk',idzk,'.dat'
      endif

 1233 format(3(2x,e22.16))
      open(unit=82,file=trim(fname))
      write(82,1233) xyz_in_src(1),xyz_in_src(2),xyz_in_src(3)
      write(82,1233) xyz_out_src(1),xyz_out_src(2),xyz_out_src(3)
      write(82,*) cf2(1)
      write(82,1233) vf2(1,1),vf2(2,1),vf2(3,1)
      write(82,1233) real(zf2(1)),real(zf2(2)),real(zf2(3))
      write(82,1233) imag(zf2(1)),imag(zf2(2)),imag(zf2(3))

      write(82,*) niter
      write(82,*) rres
      do i=1,6*npts+4*ngenus
        write(82,*) soln(i)
      enddo

      do i=1,npts
        write(82,1233) bjmcomp(1,i),bjmcomp(2,i),bjmcomp(3,i)
      enddo

      do i=1,npts
        write(82,1233) bbmcomp(1,i),bbmcomp(2,i),bbmcomp(3,i)
      enddo

      do i=1,npts
        write(82,1233) bbpcomp(1,i),bbpcomp(2,i),bbpcomp(3,i)
      enddo

      do i=1,niter
        write(82,*) errs(i)
      enddo
      close(82)


c
c
c      generate grid of targets
c
      ntarg = 0
      if(igeomtype.ge.2) then
        nlat = 101

        if(igeomtype.ge.3) then
          xyz_start(1) = 2.0d0
          xyz_start(2) = -2.0d0
          xyz_start(3) = -2.0d0

          dxyz(1) = 1.0d0
          dxyz(2) = 4.0d0/(nlat-1.0d0)
          dxyz(3) = 4.0d0/(nlat-1.0d0)

          ntarg = nlat*nlat
          allocate(targs(3,ntarg))

          ndims(1) = 1
          ndims(2) = nlat
          ndims(3) = nlat

          do i=1,nlat
            do j=1,nlat
              ipt = (i-1)*nlat + j
              targs(1,ipt) = 2.0d0
              targs(2,ipt) = xyz_start(2) + dxyz(2)*(j-1)
              targs(3,ipt) = xyz_start(3) + dxyz(3)*(i-1)
            enddo
          enddo

        endif

        if(igeomtype.eq.2) then
          xyz_start(1) = -5.9d0
          xyz_start(2) = 0.0d0
          xyz_start(3) = -1.5d0

          dxyz(1) = 3.0d0/(nlat-1.0d0)
          dxyz(2) = 1.0d0
          dxyz(3) = 3.0d0/(nlat-1.0d0)

          ntarg = nlat*nlat
          allocate(targs(3,ntarg))

          ndims(1) = nlat
          ndims(2) = 1
          ndims(3) = nlat

          do i=1,nlat
            do j=1,nlat
              ipt = (i-1)*nlat + j
              targs(1,ipt) = xyz_start(1) + dxyz(1)*(j-1)
              targs(2,ipt) = xyz_start(2) 
              targs(3,ipt) = xyz_start(3) + dxyz(3)*(i-1)
            enddo
          enddo

        endif

        allocate(uvs_targ(2,ntarg),ipatch_id_targ(ntarg))
        allocate(rpot(ntarg),isout(ntarg))


        do ipt=1,ntarg
          isout(ipt) = 0.0d0
          ipatch_id_targ(ipt) = -1
          uvs_targ(1,ipt) = 0
          uvs_targ(2,ipt) = 0
        enddo


      endif

      allocate(rsigma(npts))

      do i=1,npts
        rsigma(i) = 1.0d0
      enddo

      if(ntarg.eq.0) then
         allocate(targs(3,1), uvs_targ(2,1))
         allocate(ipatch_id_targ(1),rpot(1))
         allocate(isout(1))
      endif

      dpars2(1) = 0.0d0
      dpars2(2) = 1.0d0
      call prinf('ntarg=*',ntarg,1)

      do i=1,ntarg
        rpot(i) = 0
      enddo
      call prin2('eps=*',eps,1)
cc      call prinf('ipatch_id_targ=*',ipatch_id_targ,ntarg)
cc      call prin2('uvs_targ=*',uvs_targ,2*ntarg)
      call prinf('npatches=*',npatches,1)
      call prinf('norders=*',norders,10)
      call prinf('ixyzs=*',ixyzs,10)
      call prin2('dpars2=*',dpars2,2)
      ndtarg = 3
      
      call lap_comb_dir_eval(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,ndtarg,ntarg,targs,ipatch_id_targ,uvs_targ,eps,
     2  dpars2,rsigma,rpot)
      do i=1,ntarg
        if(abs(rpot(i)).le.1.0d-1) isout(i) = 1
      enddo
cc      call prinf('isout=*',isout,ntarg)


c
c    find near quadrature correction interactions
c
      nnz_targ = 0 
      call findnearmem(cms,npatches,rad_near,3,targs,ntarg,nnz_targ)
      call prinf('nnz_targ=*',nnz_targ,1)

      allocate(row_ptr_targ(ntarg+1),col_ind_targ(nnz_targ))
      
      call findnear(cms,npatches,rad_near,3,targs,ntarg,row_ptr_targ, 
     1        col_ind_targ)

      allocate(iquad_targ(nnz_targ+1)) 
      call get_iquad_rsc(npatches,ixyzs,ntarg,nnz_targ,
     1      row_ptr_targ,col_ind_targ,iquad_targ)
      nquad_targ = iquad_targ(nnz_targ+1)-1
      allocate(wnear_targ(7*nquad_targ))
      call getnearquad_statj_gendeb_vol(npatches,norders,ixyzs,iptype,
     1   npts,srccoefs,srcvals,3,ntarg,targs,ipatch_id_targ,uvs_targ,
     2   eps,dpars,iquadtype,nnz_targ,row_ptr_targ,col_ind_targ,
     3   iquad_targ,rfac0,nquad_targ,wnear_targ)

c
c  compute exact field at targets
c
      zk = ima*dzk
      allocate(zjm_targ(3,ntarg),zbbm_targ(3,ntarg))
      allocate(bjm_targ(3,ntarg),bbm_targ(3,ntarg))
      allocate(bbp_targ(3,ntarg))
      bbp_targ = 0
      do i=1,ntarg
       call l3ddirectcdg(1,xyz_in_src,cf2,vf2,1,targs(1,i),
     1     1,ptmp,bbp_targ(1,i),thresh)
      enddo


      
      init = 0
      call fieldsED2(zk,xyz_out_src,3,targs,ntarg,zjm_targ,
     1   zbbm_targ,zf2,init)
      bbm_targ = real(zbbm_targ)
      bjm_targ = real(zjm_targ)

      allocate(bjmcomp_targ(3,ntarg),bbmcomp_targ(3,ntarg))
      allocate(bbpcomp_targ(3,ntarg))
      call prinf('ntarg=*',ntarg,1)
      call lpcomp_statj_gendeb_postproc_vol(npatches,norders,ixyzs,
     1  iptype,npts,srccoefs,srcvals,eps,dpars,nnz,row_ptr,col_ind,
     2  iquad,nquad,wnear,ngenus,hvecs,bbphvecs,soln,novers,npts_over,
     3  ixyzso,srcover,wover,ntarg,targs,nnz_targ,row_ptr_targ,
     4  col_ind_targ,iquad_targ,nquad_targ,wnear_targ,bjmcomp_targ,
     5  bbmcomp_targ,bbpcomp_targ)

      errtargin = 0
      rtargin = 0
      ntargin = 0
      errtargout = 0
      rtargout =0
      do i=1,ntarg
        if(isout(i).eq.0) then
          ntargin = ntargin + 1
          do j=1,3
            rtargin = rtargin + (bjm_targ(j,i))**2
            rtargin = rtargin + (bbm_targ(j,i))**2
            errtargin =errtargin + (bjm_targ(j,i)-bjmcomp_targ(j,i))**2
            errtargin =errtargin + (bbm_targ(j,i)-bbmcomp_targ(j,i))**2
          enddo
        endif
        if(isout(i).eq.1) then
          bjmcomp_targ(1:3,i) = 0
          bbmcomp_targ(1:3,i) = 0
          do j=1,3
            rtargout = rtargout + bbp_targ(j,i)**2
            errtargout=errtargout+(bbp_targ(j,i)-bbpcomp_targ(j,i))**2
          enddo
        endif
      enddo
      errtargin_avg = sqrt(errtargin/(ntargin+1.0d0))
      errtargin_rel = sqrt(errtargin)/rsurfintl2
      call prin2('bjm_targ=*',bjm_targ,24)
      call prin2('bjmcomp_targ=*',bjmcomp_targ,24)
      call prin2('bbm_targ=*',bbm_targ,24)
      call prin2('bbmcomp_targ=*',bbmcomp_targ,24)

      call prin2('errtargin rel=*',errtargin_rel,1)
      call prin2('errtargin avg=*',errtargin_avg,1)
      errtargout_rel = sqrt(errtargout)/rsurfintl2
      call prin2('errtargout rel=*',errtargout_rel,1)


     
      if(igeomtype.ge.2) then
        write(fname,'(a,a,i2.2,a,i2.2,a,i1,a,i1,a,i1,a)') 
     1    trim(dirname),'bjm_',ipars(1),'_',ipars(2),'_norder',norder,
     1    '_ibg',ibg,'_idzk',idzk,'.vtk'
        write(fname1,'(a,a,i2.2,a,i2.2,a,i1,a,i1,a,i1,a)') 
     1    trim(dirname),'bbm_',ipars(1),'_',ipars(2),'_norder',norder,
     1    '_ibg',ibg,'_idzk',idzk,'.vtk'
        print *, "fname=",fname
        print *, "fname1=",fname1
        call surf_vtk_plot_vec(npatches,norders,ixyzs,iptype,npts,
     1    srccoefs,srcvals,bjmcomp,trim(fname),'a')
        call surf_vtk_plot_vec(npatches,norders,ixyzs,iptype,npts,
     1    srccoefs,srcvals,bbmcomp,trim(fname1),'a')

        write(fname,'(a,a,i2.2,a,i2.2,a,i1,a,i1,a,i1,a)') 
     1    trim(dirname),'bjm_',ipars(1),'_',ipars(2),'_norder',norder,
     1    '_ibg',ibg,'_idzk',idzk,'_plane1.vtk'
        write(fname1,'(a,a,i2.2,a,i2.2,a,i1,a,i1,a,i1,a)') 
     1    trim(dirname),'bbm_',ipars(1),'_',ipars(2),'_norder',norder,
     1    '_ibg',ibg,'_idzk',idzk,'_plane1.vtk'
        print *, "fname=",fname
        print *, "fname1=",fname1


       call vtk_write_plane_vec(ndims,ntarg,xyz_start,dxyz,bjmcomp_targ,
     1     'abc',trim(fname))
       call vtk_write_plane_vec(ndims,ntarg,xyz_start,dxyz,bbmcomp_targ,
     1   'abc',trim(fname1))
      endif


      return
      end







      subroutine setup_geom(igeomtype,norder,npatches,ipars, 
     1    srcvals,srccoefs,ifplot,fname)
      implicit real *8 (a-h,o-z)
      integer igeomtype,norder,npatches,ipars(*),ifplot
      character (len=*) fname
      character (len=1000) fname2
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
      
      if(igeomtype.eq.6) then
        npatches0 = npatches/2
        npatchestmp = 0 
        done = 1
        pi = atan(done)*4
        allocate(triaskel(3,3,npatches0))


        umin = 0
        umax = 2*pi
        vmin = 0
        vmax = 2*pi
        nover = 0
        call xtri_rectmesh_ani(umin,umax,vmin,vmax,ipars(1),ipars(2),
     1     nover,npatches0,npatchestmp,triaskel)
        call prinf('npatches=*',npatches,1)

        ll = len(trim(fname))
        fname2 = trim(fname(1:ll-4)//'-1.vtk')        
         
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
           call xtri_vtk_surf(fname2,npatches0,xtri_geometry,
     1         ptr1,ptr2,ptr3,ptr4,norder,
     2         'Triangulated surface of the torus')
        endif

        call getgeominfo(npatches0,xtri_geometry,ptr1,ptr2,ptr3,ptr4,
     1     npols,uvs,umatr,srcvals,srccoefs)


        umin = 0
        umax = 2*pi
        vmin = 2*pi
        vmax = 0
        nover = 0
        call xtri_rectmesh_ani(umin,umax,vmin,vmax,ipars(1),ipars(2),
     1     nover,npatches0,npatchestmp,triaskel)
        call prinf('npatches=*',npatches,1)

        ll = len(trim(fname))
        fname2 = trim(fname(1:ll-4)//'-2.vtk')        
         
        p1(1) = 0.5d0
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
           call xtri_vtk_surf(fname2,npatches0,xtri_geometry,
     1         ptr1,ptr2,ptr3,ptr4,norder,
     2         'Triangulated surface of the torus')
        endif

        npts0 = npatches0*npols
        print *, "npts0=*",npts0
        print *, "npatches0=*",npatches0
        

        call getgeominfo(npatches0,xtri_geometry,ptr1,ptr2,ptr3,ptr4,
     1     npols,uvs,umatr,srcvals(1,npts0+1),srccoefs(1,npts0+1))


      endif

      
      return  
      end
c
c
c
c
c
      subroutine get_stell_testtarg(ntargtest,xyz_in,xyz_out,npatches,
     1    npts,cms,rads,srcvals,wts,htest)
      implicit real *8 (a-h,o-z)
      real *8 triainfo(3,3),cms(3,npatches),rads(npatches),htest
      real *8 xyz(3),dxyzduv(3,2),rnorm(3),xtest(3)
      real *8 xyz_in(3,ntargtest),xyz_out(3,ntargtest)
      real *8 srcvals(12,npts),wts(npts)
      real *8, allocatable :: deltas(:,:)
      logical isout
      integer iseed
   
      done = 1
      pi = atan(done)*4

      triainfo(1:3,1:3) = 0
      triainfo(1,2) = 2*pi
      triainfo(2,3) = 2*pi

      mmax = 2
      nmax = 1

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

      ii = 0
      ictr = 0
      do while (ii<ntargtest)
         ictr = ictr+1
         u = hkrand(0)
         v = hkrand(0)
         r = 0.1d0 + 0.9d0*hkrand(0)
         call xtri_stell_eval(1,u,v,xyz,dxyzduv,triainfo,deltas,
     1      mmax,nmax)
         call cross_prod3d(dxyzduv(1,1),dxyzduv(1,2),rnorm)
         ds = sqrt(rnorm(1)**2 + rnorm(2)**2 + rnorm(3)**2)
         rnorm(1:3) = rnorm(1:3)/ds
         
         xtest(1:3) = xyz(1:3) + rnorm(1:3)*r
         do ipatch=1,npatches
           rr = (xtest(1) - cms(1,ipatch))**2 + 
     1        (xtest(2)-cms(2,ipatch))**2 + 
     1        (xtest(3)-cms(3,ipatch))**2
           if(rr.lt.(htest*rads(ipatch))**2) then
             goto 1111
           endif
         enddo
         isout = .false.
         call test_exterior_pt(npts,srcvals,wts,xtest,isout)
         if(.not.isout) then
           ii = ii + 1
           xyz_in(1:3,ii) = xtest(1:3)
         endif
 1111    continue
      enddo
      call prinf('ictr in=*',ictr,1)


      ii = 0
      ictr = 0
      do while (ii<ntargtest)
         ictr = ictr+1
         u = hkrand(0)
         v = hkrand(0)
         r = 1.5d0 + 10.0d0*hkrand(0)
         call xtri_stell_eval(1,u,v,xyz,dxyzduv,triainfo,deltas,
     1      mmax,nmax)
         call cross_prod3d(dxyzduv(1,1),dxyzduv(1,2),rnorm)
         ds = sqrt(rnorm(1)**2 + rnorm(2)**2 + rnorm(3)**2)
         rnorm(1:3) = rnorm(1:3)/ds
         
         xtest(1:3) = xyz(1:3) - rnorm(1:3)*r
         do ipatch=1,npatches
           rr = (xtest(1) - cms(1,ipatch))**2 + 
     1        (xtest(2)-cms(2,ipatch))**2 + 
     1        (xtest(3)-cms(3,ipatch))**2
           if(rr.lt.(htest*rads(ipatch))**2) then
             goto 1211
           endif
         enddo
         isout = .false.
         call test_exterior_pt(npts,srcvals,wts,xtest,isout)
         if(isout) then
           ii = ii + 1
           xyz_out(1:3,ii) = xtest(1:3)
         endif
 1211    continue
      enddo
      call prinf('ictr out=*',ictr,1)

      return
      end
c
c
c
c
c
      subroutine get_sphere_testtarg(n,xyz_in,xyz_out)
      implicit real *8 (a-h,o-z)
      real *8 xyz_in(3,n),xyz_out(3,n)

      done = 1
      pi = atan(done)*4

      do i=1,n
        rr = hkrand(0)*0.67d0
        thet = hkrand(0)*pi
        phi = hkrand(0)*2*pi
        xyz_in(1,i) = rr*sin(thet)*cos(phi)
        xyz_in(2,i) = rr*sin(thet)*sin(phi)
        xyz_in(3,i) = rr*cos(thet)
      enddo

      do i=1,n
        rr = 1.33d0 + hkrand(0)*0.67d0
        thet = hkrand(0)*pi
        phi = hkrand(0)*2*pi
        xyz_out(1,i) = rr*sin(thet)*cos(phi)
        xyz_out(2,i) = rr*sin(thet)*sin(phi)
        xyz_out(3,i) = rr*cos(thet)
      enddo

      return
      end
c
c
c
c
c

      subroutine get_torus_testtarg(n,xyz_in,xyz_out)
      implicit real *8 (a-h,o-z)
      real *8 xyz_in(3,n),xyz_out(3,n)

      done = 1
      pi = atan(done)*4

      do i=1,n
        uu = hkrand(0)*2*pi
        vv = hkrand(0)*2*pi
        rr = 0.67d0*hkrand(0)
        xyz_in(1,i) = (rr*cos(uu) + 2)*cos(vv)
        xyz_in(2,i) = (rr*cos(uu) + 2)*sin(vv)
        xyz_in(3,i) = rr*sin(uu)
      enddo

      do i=1,n
        uu = hkrand(0)*2*pi
        vv = hkrand(0)*2*pi
        rr = 0.67d0*hkrand(0) + 1.33d0
        xyz_out(1,i) = (rr*cos(uu) + 2)*cos(vv)
        xyz_out(2,i) = (rr*cos(uu) + 2)*sin(vv)
        xyz_out(3,i) = rr*sin(uu)
      enddo

      return
      end
c
c
c
c
c

      subroutine get_wtorus_testtarg(n,xyz_in,xyz_out)
      implicit real *8 (a-h,o-z)
      real *8 xyz_in(3,n),xyz_out(3,n)

      done = 1
      pi = atan(done)*4

      do i=1,n
        uu = hkrand(0)*2*pi
        vv = hkrand(0)*2*pi
        rr = 0.67d0*hkrand(0)
        xyz_in(1,i) = (rr*cos(uu) + 2+0.25d0*cos(3*vv))*cos(vv)
        xyz_in(2,i) = (rr*cos(uu) + 2+0.25d0*cos(3*vv))*sin(vv)
        xyz_in(3,i) = rr*sin(uu)
      enddo

      do i=1,n
        uu = hkrand(0)*2*pi
        vv = hkrand(0)*2*pi
        rr = 0.67d0*hkrand(0) + 1.33d0
        xyz_out(1,i) = (rr*cos(uu) + 2+0.25d0*cos(3*vv))*cos(vv)
        xyz_out(2,i) = (rr*cos(uu) + 2+0.25d0*cos(3*vv))*sin(vv)
        xyz_out(3,i) = rr*sin(uu)
      enddo

      return
      end
c
c
c
c
c



      subroutine test_exterior_pt(npts,srcvals,wts,xyzout,isout)
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
c    npts - integer
c       total number of discretization points on the surface
c    srvals - real *8 (12,npts)
c       source geometry info
c    xyzout -  real *8 (3)
c       point to be tested
c
c  output: 
c    isout - boolean
c      whether the target is in the interior or not
c

      implicit none
      integer npts
      real *8 srcvals(12,npts),xyzout(3),wts(npts)
      real *8 tmp(3)
      real *8 dpars,done,pi
      integer ipars,norderhead,nd

      integer ipatch,j,i
      real *8 ra,ds,r,dx,dy,dz,val
      logical isout

      done = 1
      pi = atan(done)*4


      ra = 0
      do i=1,npts
        dx = xyzout(1) - srcvals(1,i)
        dy = xyzout(2) - srcvals(2,i)
        dz = xyzout(3) - srcvals(3,i)

        r = sqrt(dx**2 + dy**2 + dz**2)
          
        val = dx*srcvals(10,i) + dy*srcvals(11,i) + dz*srcvals(12,i)
        val = val/r**3 
        ra = ra + val*wts(i)
      enddo

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






