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


      real *8, allocatable :: bbp(:,:)
      real *8, allocatable :: bbp_a(:,:),bbp_b(:,:)
      real *8, allocatable :: rhs(:),soln(:)
      real *8, allocatable :: errs(:)
      real *8, allocatable :: bbpcomp(:,:)
      real *8, allocatable :: bbpcomp_targ(:,:)
      real *8, allocatable :: bbp_targ(:,:)
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

c
c  Note for this code, igeomtype can only be 1,2,3,4 
c
c  igeomtype = 1 => sphere
c  igeomtype = 2 => stellartor
c  igeomtype = 3 => wiggly torus
c  igeomtype = 4 => torus
c
      igeomtype = 2
      iref = 2
      if(igeomtype.eq.1) then
        ipars(1) = iref
        npatches = 12*(4**(ipars(1)))
        fname = 'sphere.vtk'
        dirname = '/mnt/home/mrachh/ceph/' // 
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
        dirname = '/mnt/home/mrachh/ceph/' // 
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

      norder = 8
      npols = (norder+1)*(norder+2)/2

      npts = npatches*npols
      call prinf('npatches=*',npatches,1)

      allocate(srcvals(12,npts),srccoefs(9,npts))
      ifplot = 1

      call setup_geom(igeomtype,norder,npatches,ipars, 
     1       srcvals,srccoefs,ifplot,fname)

      allocate(norders(npatches),ixyzs(npatches+1),iptype(npatches))

      if(igeomtype.eq.2) call write_triaskelpts(norder,npatches,ipars)
      stop

      do i=1,npatches
        norders(i) = norder
        ixyzs(i) = 1 +(i-1)*npols
        iptype(i) = 1
      enddo

      ixyzs(npatches+1) = 1+npols*npatches
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
c
      if(igeomtype.ge.2) then
        write(fname,'(a,a,i3.3,a,i3.3,a,i1,a)') trim(dirname),
     1       'hvecs_',ipars(1),
     1      '_',ipars(2),'_',norder,'_1.dat'
          open(unit=78,file=trim(fname),form='unformatted')
        write(fname,'(a,a,i3.3,a,i3.3,a,i1,a)') trim(dirname),
     1       'hvecs_',ipars(1),
     1       '_',ipars(2),'_',norder,'_2.dat'
        print *, fname 
        open(unit=79,file=trim(fname),form='unformatted')
      endif


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


      thresh = 1.0d-16
      allocate(bbp(3,npts),ptmp(npts))
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

      allocate(rhs(npts+2*ngenus),soln(npts+2*ngenus))
      rhs = 0
      soln = 0
      allocate(bbp_a(3,na),bbp_b(3,nb))
      if(igeomtype.eq.2.or.igeomtype.eq.5) then
        rhs(npts+2) = 1
      endif

      if(igeomtype.eq.4) then
        rhs(npts+1) = 1
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
      allocate(wnear(3*nquad))
      wnear = 0
      dpars = 0

      


      t1 = 0
      t2 = 0
cc      goto 1111
      call cpu_time(t1)
C$      t1 = omp_get_wtime()      
      call getnearquad_statj_gendeb_grads0(npatches,norders,
     1      ixyzs,iptype,npts,srccoefs,srcvals,
     1      eps,dzk,iquadtype,nnz,row_ptr,col_ind,
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
     2    wnear,rhstmp,novers,npts_over,ixyzso,
     2    srcover,wover,outtmp)
        do i=1,npts
          bbphvecs(1,i,igen) = outtmp(i)
          bbphvecs(2,i,igen) = outtmp(i+npts)
          bbphvecs(3,i,igen) = outtmp(i+2*npts)
        enddo
        do j=1,npts
          call cross_prod3d(srcvals(10,j),hvecs(1,j,igen),vtmp1)
          bbphvecs(1:3,j,igen) = bbphvecs(1:3,j,igen) - vtmp1(1:3)/2
        enddo
      enddo
      print *, "finished computing bbphvecs"

      numit = 300
      allocate(errs(numit+1))
      call prinf('ngenus=*',ngenus,1)
      eps_gmres = 1.0d-8
      if(ngenus.ge.1)  call prin2('rhs_projs=*',rhs(npts+1),2)
      print *, "here"


      call cpu_time(t1)
C$       t1 = omp_get_wtime()      

      call statj_gendeb_lambdainf_solver(npatches,norders,ixyzs,iptype,
     1  npts,srccoefs,srcvals,eps,dpars,ngenus,hvecs,bbphvecs,na,iaxyzs,
     2  apatches,auv,avals,awts,nb,ibxyzs,bpatches,buv,bvals,bwts,
     3  numit,rhs,eps_gmres,niter,errs,rres,soln)
      call cpu_time(t2)
C$       t2 = omp_get_wtime()      
      
      call prin2('solve time=*',t2-t1,1)
      
      if(ngenus.ge.1) call prin2('projs=*',soln(npts+1),2)


      allocate(bbpcomp(3,npts))

      print *, "here"
      call prin2('eps=*',eps,1)
      call prin2('soln=*',soln,24)
      call prinf('ngenus=*',ngenus,1)
      call prin2('dpars=*',dpars,3)


      call lpcomp_statj_gendeb_lambdainf_postproc(npatches,norders,
     1  ixyzs,iptype,npts,srccoefs,srcvals,eps,dpars,nnz,row_ptr,
     2  col_ind,iquad,nquad,wnear,ngenus,hvecs,bbphvecs,soln,
     3  novers,npts_over,ixyzso,srcover,wover,bbpcomp)
      call prin2('bbpcomp=*',bbpcomp,24)

       
      write(fname,'(a,a,i2.2,a,i2.2,a,i1,a)') 
     1    trim(dirname),'statj_soln_lambdainf_',ipars(1),'_',ipars(2),
     2    '_norder',norder,'.dat'


 1233 format(3(2x,e22.16))
      open(unit=80,file=trim(fname),form='unformatted')

      write(80) niter
      write(80) rres
      write(80) soln
      write(80) bbpcomp
      write(80) errs

      close(80)




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






      subroutine write_triaskelpts(norder,npatches,ipars)
      implicit real *8 (a-h,o-z)
      real *8, allocatable :: triaskel(:,:,:)
      real *8, allocatable :: uvs(:,:)
      real *8, allocatable :: uvpts(:,:)
      integer ipars(2)

      npols = (norder+1)*(norder+2)/2
      npts = npatches*npols

      allocate(uvs(2,npols),uvpts(2,npts))

      call get_vioreanu_nodes(norder,npols,uvs)

      done = 1.0d0
      pi = atan(done)*4



      allocate(triaskel(3,3,npatches))
      umin = 0
      umax = 2*pi
      vmin = 2*pi
      vmax = 0
      nover = 0
      call xtri_rectmesh_ani(umin,umax,vmin,vmax,ipars(1),ipars(2),
     1     nover,npatches,npatches,triaskel)
      
      do itri=1,npatches 
        do j=1,npols
          ipt = (itri-1)*npols + j

          x0=triaskel(1,1,itri)
          y0=triaskel(2,1,itri)
          z0=triaskel(3,1,itri)

          x1=triaskel(1,2,itri)
          y1=triaskel(2,2,itri)
          z1=triaskel(3,2,itri)

          x2=triaskel(1,3,itri)
          y2=triaskel(2,3,itri)
          z2=triaskel(3,3,itri)

          uvpts(1,ipt) = x0 + uvs(1,j)*(x1-x0) + uvs(2,j)*(x2-x0) 
          uvpts(2,ipt) = y0 + uvs(1,j)*(y1-y0) + uvs(2,j)*(y2-y0) 

        enddo
      enddo

      open(unit=334,file='triaskel_stellinfo.dat',form='unformatted')
      write(334) uvpts
      close(334)


      return
      end
