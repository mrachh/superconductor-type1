      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: wts(:),rsigma(:),rpot(:)
      integer ipars(2)

      integer, allocatable :: norders(:),ixyzs(:),iptype(:)

      real *8 xyz_out(3),xyz_in(3,2)
      real *8 vtmp1(3),bbpc(3),bbpex(3)
      real *8, allocatable :: ffform(:,:,:),ffformex(:,:,:)
      real *8, allocatable :: ffforminv(:,:,:),ffformexinv(:,:,:)
      real *8, allocatable :: curv(:)

      complex *16, allocatable :: zbbm(:,:),zbbp(:,:),zjm(:,:)
      real *8, allocatable :: bbm(:,:),bbp(:,:),bjm(:,:)
      real *8, allocatable :: bbp_a(:,:),bbp_b(:,:)
      real *8, allocatable :: rhs(:),soln(:)
      real *8, allocatable :: errs(:)
      real *8, allocatable :: bbpcomp(:,:),bjmcomp(:,:),bbmcomp(:,:)
      real *8, allocatable :: bjmcomp_targ(:,:),bbmcomp_targ(:,:)
      real *8, allocatable :: bbpcomp_targ(:,:)
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
      real *8 vf2(3,2),dpars(3),cf2(2),dpars2(2)
      real *8 xyz_start(3), dxyz(3)
      complex * 16 zpars(3)
      integer numit,niter,ndims(3)
      character *100 title,dirname
      character *1000 fname,fname1,fname2,fname3,fname4,fname5


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



      ibg = 0
      idzk = 10
      igeomtype = 2 
      if(igeomtype.eq.4) then
        ipars(1) = 4*4
        ipars(2) = 2*4
        npatches = 2*ipars(1)*ipars(2)
        dirname = '/mnt/home/mrachh/ceph/' // 
     1     'superconductor-type1-data/torus-data/'
        fname = 'torus.vtk'

        xyz_in(1,1) = 2.001d0
        xyz_in(2,1) = 0.002d0
        xyz_in(3,1) = 0.001d0

        xyz_in(1,2) = 2.003d0
        xyz_in(2,2) = -0.01d0
        xyz_in(3,2) = 0.001d0

        xyz_out(1) = -3.5d0
        xyz_out(2) = 7.1d0
        xyz_out(3) = 5.7d0

        ngenus = 1
      endif


      if(igeomtype.eq.2) then
        ipars(1) = 20
        ipars(2) = 60
        npatches = 2*ipars(1)*ipars(2)
        
        fname = 'stellarator.vtk'
        dirname = '/mnt/home/mrachh/ceph/' // 
     1     'superconductor-type1-data/stell-data/'

        xyz_in(1,1) = -4.501d0
        xyz_in(2,1) = 0.002d0
        xyz_in(3,1) = 0.001d0

        xyz_in(1,2) = -4.503d0
        xyz_in(2,2) = -0.01d0
        xyz_in(3,2) = 0.001d0

        xyz_out(1) = -3.5d0
        xyz_out(2) = 10.1d0
        xyz_out(3) = 5.7d0

        ngenus = 1

      endif

      if(igeomtype.eq.1) then
        ipars(1) = 1
        npatches = 12*(4**(ipars(1)))
        fname = 'sphere.vtk'
        dirname = '/mnt/home/mrachh/ceph/' // 
     1     'superconductor-type1-data/sphere-data/'

        xyz_in(1,1) = 0.201d0
        xyz_in(2,1) = 0.102d0
        xyz_in(3,1) = 0.011d0

        xyz_in(1,2) = -0.3003d0
        xyz_in(2,2) = -0.01d0
        xyz_in(3,2) = 0.001d0

        xyz_out(1) = -3.5d0
        xyz_out(2) = 7.1d0
        xyz_out(3) = 5.7d0

        ngenus = 0
      endif



      norder = 8 
      npols = (norder+1)*(norder+2)/2

      npts = npatches*npols
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
      allocate(wts(npts))

      call get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,wts)
      call test_exterior_pt(npatches,norders,npts,
     1  srcvals,srccoefs,wts,xyz_in,isout0)
      print *, "isout=",isout0

      call test_exterior_pt(npatches,norders,npts,
     1  srcvals,srccoefs,wts,xyz_out,isout1)
      print *, "isout=",isout1

c
c  set file names for reading or writing harmonic vector fields
c
c
      write(fname,'(a,a,i3.3,a,i3.3,a,i1,a)') trim(dirname),
     1     'hvecs_',ipars(1),
     1    '_',ipars(2),'_',norder,'_1.dat'
        open(unit=78,file=trim(fname),form='unformatted')
      write(fname,'(a,a,i3.3,a,i3.3,a,i1,a)') trim(dirname),
     1     'hvecs_',ipars(1),
     1     '_',ipars(2),'_',norder,'_2.dat'
      print *, fname 
      open(unit=79,file=trim(fname),form='unformatted')

c
c  set a and b cycle params
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
          call get_harm_vec_field(npatches,norders,ixyzs,iptype, 
     1      npts,srccoefs,srcvals,wts,eps,hvecs(1,1,1),errest)
          call prin2('errest=*',errest,1)
          do i=1,npts
            call cross_prod3d(srcvals(10,i),hvecs(1,i,1),hvecs(1,i,2))
          enddo
        else
          read(78) hvecs(1:3,1:npts,1)
          read(79) hvecs(1:3,1:npts,2)
          close(78)
          close(79)
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
c   compute the boundary data, for now assume that only external
c   B field is applied and that the interior fields are 0 
c
      vf2(1,1) = hkrand(0)
      vf2(2,1) = hkrand(0)
      vf2(3,1) = hkrand(0)

      vf2(1:3,1) = vf2(1:3,1)*1.0d0*1 
      vf2(1:3,2) = 0
      cf2(1) = 1*0
      cf2(2) = -1*0 

      

      thresh = 1.0d-16
      allocate(bbp(3,npts),ptmp(npts))
      allocate(bjm(3,npts),bbm(3,npts))
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

      call l3ddirectcdg(1,xyz_in,cf2,vf2,2,sources,npts,
     1 ptmp,bbp,thresh)
      call prin2('bbp=*',bbp,24)
      
      call prin2('ptmp=*',ptmp,24)
      allocate(rhs(6*npts+4*ngenus),soln(6*npts+4*ngenus))
      rhs = 0
      soln = 0
      dzk = 10**(idzk)
      if(idzk.eq.3) dzk = 30.0d0
      if(idzk.ge.4) dzk = 2**(idzk-4)
      print *, "dzk=",dzk
cc      call surf_vtk_plot_scalar(npatches,norders,ixyzs,iptype,npts, 
cc     1   srccoefs,srcvals,rhs(2*npts+1),'rhs3-torus.vtk','a')

      if(igeomtype.eq.4) then

        allocate(bbp_a(3,na),bbp_b(3,nb))
        rhs(6*npts+1) = 1 
        rhs(6*npts+3) = 1
      endif
c
c  The a and b cycles for the stellarator discretization are
c  flipped
c
      if(igeomtype.eq.2.or.igeomtype.eq.5) then
        allocate(bbp_a(3,na),bbp_b(3,nb))
        rhs(6*npts+2) = 1
        rhs(6*npts+4) = 1
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


      allocate(cms(3,npatches),rads(npatches),rad_near(npatches))

      call get_centroid_rads(npatches,norders,ixyzs,iptype,npts, 
     1     srccoefs,cms,rads)
      call prin2('cms=*',cms,24)
      call prin2('rads=*',rads,12)

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
        rgamma = 1.0d0
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
      call getnearquad_statj_gendeb(npatches,norders,
     1      ixyzs,iptype,npts,srccoefs,srcvals,
     1      eps,dzk,iquadtype,nnz,row_ptr,col_ind,
     1      iquad,rfac0,nquad,wnear)
 1111 continue     
      call prinf('finished generating near quadrature correction*',i,0)

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

      numit = 330
      allocate(errs(numit+1))
      call prinf('ngenus=*',ngenus,1)
      eps_gmres = 1.0d-8
      ifreadsol = 0
 1233 format(3(2x,e22.16))

      write(fname,'(a,a,i2.2,a,i2.2,a,i1,a,i1,a,i2.2,a)') 
     1    trim(dirname),'statj_soln_',ipars(1),'_',ipars(2),
     2    '_norder',norder,'_ibg',ibg,'_idzk',idzk,'.dat'

      write(fname1,'(a,a,i2.2,a,i2.2,a,i1,a,i1,a,i2.2,a)') 
     1    trim(dirname),'statj_bjm_',ipars(1),'_',ipars(2),
     2    '_norder',norder,'_ibg',ibg,'_idzk',idzk,'.vtk'

      write(fname2,'(a,a,i2.2,a,i2.2,a,i1,a,i1,a,i2.2,a)') 
     1    trim(dirname),'statj_bbm_',ipars(1),'_',ipars(2),
     2    '_norder',norder,'_ibg',ibg,'_idzk',idzk,'.vtk'

      write(fname3,'(a,a,i2.2,a,i2.2,a,i1,a,i1,a,i2.2,a)') 
     1    trim(dirname),'statj_bbp_',ipars(1),'_',ipars(2),
     2    '_norder',norder,'_ibg',ibg,'_idzk',idzk,'.vtk'

      write(fname4,'(a,a,i2.2,a,i2.2,a,i1,a,i1,a,i2.2,a)') 
     1    trim(dirname),'statj_bjm_',ipars(1),'_',ipars(2),
     2    '_norder',norder,'_ibg',ibg,'_idzk',idzk,'_plane2.vtk'

      write(fname5,'(a,a,i2.2,a,i2.2,a,i1,a,i1,a,i2.2,a)') 
     1    trim(dirname),'statj_bbm_',ipars(1),'_',ipars(2),
     2    '_norder',norder,'_ibg',ibg,'_idzk',idzk,'_plane2.vtk'


      allocate(bbpcomp(3,npts),bjmcomp(3,npts),bbmcomp(3,npts))
      if(ifreadsol.eq.0) then
        call statj_gendeb_solver(npatches,norders,ixyzs,iptype,npts,
     1    srccoefs,srcvals,eps,dpars,ngenus,hvecs,bbphvecs,na,iaxyzs,
     2    apatches,auv,avals,awts,nb,ibxyzs,bpatches,buv,bvals,bwts,
     3    numit,rhs,eps_gmres,niter,errs,rres,soln)
      
        print *, "here"
        call prin2('eps=*',eps,1)
        call prin2('soln=*',soln,24)
        call prinf('ngenus=*',ngenus,1)
        call prin2('dpars=*',dpars,3)

        call lpcomp_statj_gendeb_postproc(npatches,norders,ixyzs,
     1    iptype,npts,srccoefs,srcvals,eps,dpars,nnz,row_ptr,col_ind,
     2    iquad,nquad,wnear,ngenus,hvecs,bbphvecs,soln,novers,npts_over,
     3    ixyzso,srcover,wover,bjmcomp,bbmcomp,bbpcomp)

        open(unit=80,file=fname,form='unformatted')
        write(80) niter
        write(80) rres
        write(80) soln
        write(80) bjmcomp
        write(80) bbmcomp
        write(80) bbpcomp
        write(80) errs(1:niter)

        close(80)
      else
        open(unit=80,file=fname,form='unformatted')
        read(80) niter
        read(80) rres
        read(80) soln
        read(80) bjmcomp
        read(80) bbmcomp
        read(80) bbpcomp
        read(80) errs(1:niter)

      endif
      call prin2('errs=*',errs,niter)
      
     
      call surf_vtk_plot_vec(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,bjmcomp,trim(fname1),'a')
      call surf_vtk_plot_vec(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,bbmcomp,trim(fname2),'a')
      call surf_vtk_plot_vec(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,bbpcomp,trim(fname3),'a')

c
c
c      generate grid of targets
c
      ntarg = 0
      if(igeomtype.ge.2) then
        nlat = 1001

        if(igeomtype.eq.4.or.igeomtype.eq.5) then
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

          print *, "Here"
        endif

        allocate(uvs_targ(2,ntarg),ipatch_id_targ(ntarg))
        allocate(rsigma(npts),rpot(ntarg),isout(ntarg))

        do i=1,npts
          rsigma(i) = 1.0d0
        enddo

        do ipt=1,ntarg
          isout(ipt) = 0.0d0
          ipatch_id_targ(ipt) = -1
          uvs_targ(1,ipt) = 0
          uvs_targ(2,ipt) = 0
        enddo


      endif

      dpars2(1) = 0.0d0
      dpars2(2) = 1.0d0
      call prinf('ntarg=*',ntarg,1)

      do i=1,ntarg
        rpot(i) = 0
      enddo
      call prin2('eps=*',eps,1)
      call prinf('npatches=*',npatches,1)
      call prinf('norders=*',norders,10)
      call prinf('ixyzs=*',ixyzs,10)
      call prin2('dpars2=*',dpars2,2)
      ndtarg = 3
      
      call lpcomp_lap_comb_dir(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,ndtarg,ntarg,targs,ipatch_id_targ,uvs_targ,eps,
     2  dpars2,rsigma,rpot)
      do i=1,ntarg
        if(abs(rpot(i)).le.1.0d-1) isout(i) = 1
      enddo


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

      allocate(bjmcomp_targ(3,ntarg),bbmcomp_targ(3,ntarg))
      allocate(bbpcomp_targ(3,ntarg))
      call prinf('ntarg=*',ntarg,1)
      call lpcomp_statj_gendeb_postproc_vol(npatches,norders,ixyzs,
     1  iptype,npts,srccoefs,srcvals,eps,dpars,nnz,row_ptr,col_ind,
     2  iquad,nquad,wnear,ngenus,hvecs,bbphvecs,soln,novers,npts_over,
     3  ixyzso,srcover,wover,ntarg,targs,nnz_targ,row_ptr_targ,
     4  col_ind_targ,iquad_targ,nquad_targ,wnear_targ,bjmcomp_targ,
     5  bbmcomp_targ,bbpcomp_targ)

      do i=1,ntarg
        if(isout(i).eq.1) then
          bjmcomp_targ(1:3,i) = 0
          bbmcomp_targ(1:3,i) = 0
        endif
      enddo


      call vtk_write_plane_vec(ndims,ntarg,xyz_start,dxyz,bjmcomp_targ,
     1   'abc',trim(fname4))
      call vtk_write_plane_vec(ndims,ntarg,xyz_start,dxyz,bbmcomp_targ,
     1   'abc',trim(fname5))


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

      if(abs(ra+4*pi).le.1.0d-1) isout = .false.
      if(abs(ra).le.1.0d-1) isout = .true.

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






