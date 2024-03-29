      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: wts(:),rsigma(:)
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
      real *8 vf2(3,2),dpars(3),cf2(2)
      complex * 16 zpars(3)
      integer numit,niter
      character *100 title,dirname
      character *300 fname

      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      real *8, allocatable :: wnear(:)
      real *8, allocatable :: targs(:,:)
      real *8, allocatable :: cms(:,:),rads(:),rad_near(:)
      integer, allocatable :: iquad(:),row_ptr(:),col_ind(:)
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



      ibg = 3

      igeomtype = 2
      if(igeomtype.eq.5) then
        ipars(1) = 20
        ipars(2) = 60
        npatches = 2*ipars(1)*ipars(2)
        fname = 'torus2.vtk'

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

      if(igeomtype.eq.4) then
        ipars(1) = 4*8
        ipars(2) = 2*8
        npatches = 2*ipars(1)*ipars(2)
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

      if(igeomtype.eq.3) then
        ipars(1) = 8*8
        ipars(2) = 8*4
        npatches = 2*ipars(1)*ipars(2)
        
        fname = 'wtorus.vtk'

        xyz_in(1,1) = -2.001d0
        xyz_in(2,1) = 0.002d0
        xyz_in(3,1) = 0.001d0

        xyz_in(1,2) = -2.003d0
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
      allocate(wts(npts))

      isout = .false.
      call get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,wts)
      call test_exterior_pt(npatches,norders,npts,
     1  srcvals,srccoefs,wts,xyz_in,isout0)
      print *, "isout=",isout0

      call test_exterior_pt(npatches,norders,npts,
     1  srcvals,srccoefs,wts,xyz_out,isout1)
      print *, "isout=",isout1


      if(igeomtype.eq.2) then
        write(fname,'(a,i2.2,a,i2.2,a,i1,a)') 'stell_hvecs_',ipars(1),
     1    '_',ipars(2),'_',norder,'_1.dat'
        open(unit=78,file=trim(fname))
        write(fname,'(a,i2.2,a,i2.2,a,i1,a)') 'stell_hvecs_',ipars(1),
     1     '_',ipars(2),'_',norder,'_2.dat'
        print *, fname 
        open(unit=79,file=trim(fname))
      endif
      if(igeomtype.eq.4) then
        write(fname,'(a,i2.2,a,i2.2,a,i1,a)') 'torus_hvecs_',ipars(1),
     1     '_',ipars(2),'_',norder,'_1.dat'
        open(unit=78,file=trim(fname))
        write(fname,'(a,i2.2,a,i2.2,a,i1,a)') 'torus_hvecs_',ipars(1),
     1     '_',ipars(2),'_',norder,'_2.dat'
        open(unit=79,file=trim(fname))
      endif

      if(igeomtype.eq.5) then
        write(fname,'(a,i2.2,a,i2.2,a,i1,a)') 'torus2_hvecs_',ipars(1),
     1     '_',ipars(2),'_',norder,'_1.dat'
        open(unit=78,file=trim(fname))
        write(fname,'(a,i2.2,a,i2.2,a,i1,a)') 'torus2_hvecs_',ipars(1),
     1     '_',ipars(2),'_',norder,'_2.dat'
        open(unit=79,file=trim(fname))
      endif


      if(igeomtype.eq.3) then
        write(fname,'(a,i2.2,a,i2.2,a,i1,a)') 'wtorus_hvecs_',ipars(1),
     1     '_',ipars(2),'_',norder,'_1.dat'
        open(unit=78,file=trim(fname))
        write(fname,'(a,i2.2,a,i2.2,a,i1,a)') 'wtorus_hvecs_',ipars(1),
     1     '_',ipars(2),'_',norder,'_2.dat'
        open(unit=79,file=trim(fname))
      endif



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
          eps = 1.0d-7
          call get_harm_vec_field(npatches,norders,ixyzs,iptype, 
     1      npts,srccoefs,srcvals,wts,eps,hvecs(1,1,1),errest)
          call prin2('errest=*',errest,1)
          do i=1,npts
            call cross_prod3d(srcvals(10,i),hvecs(1,i,1),hvecs(1,i,2))
          enddo
        else
          do i=1,npts
            read(78,*) hvecs(1,i,1),hvecs(2,i,1),hvecs(3,i,1)
            read(79,*) hvecs(1,i,2),hvecs(2,i,2),hvecs(3,i,2)
          enddo
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
          do i=1,npts
            write(78,*) hvecs(1,i,1),hvecs(2,i,1),hvecs(3,i,1)
            write(79,*) hvecs(1,i,2),hvecs(2,i,2),hvecs(3,i,2)
          enddo
          close(78)
          close(79)
        endif
      endif

c
c
c   compute the boundary data, for now assume that only external
c   B field is applied and that the interior fields are 0 
c
      vf2(1:3,1) = 1 
      vf2(1:3,2) = 0
      cf2(1) = 1*1
      cf2(2) = -1*0 

      

      thresh = 1.0d-16
      allocate(bbp(3,npts),ptmp(npts))
      allocate(bjm(3,npts),bbm(3,npts))
      ptmp = 0
      bbp = 0
      bjm = 0
      bbm = 0

      ntarg = npts
      allocate(targs(3,npts))
C$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,npts 
        targs(1,i) = srcvals(1,i)
        targs(2,i) = srcvals(2,i)
        targs(3,i) = srcvals(3,i)
      enddo
C$OMP END PARALLEL DO      

      call l3ddirectcdg(1,xyz_in,cf2,vf2,2,targs,npts,
     1 ptmp,bbp,thresh)
      call prin2('bbp=*',bbp,24)
      
      call prin2('ptmp=*',ptmp,24)
      allocate(rhs(6*npts+4*ngenus),soln(6*npts+4*ngenus))
      rhs = 0
      soln = 0
      dzk = 1.0d0
      do i=1,npts
        rhs(i) = bbm(1,i)/dzk-bbp(1,i)
        rhs(i+npts) = bbm(2,i)/dzk-bbp(2,i)
        rhs(i+2*npts) = bbm(3,i)/dzk-bbp(3,i)
        rhs(i+3*npts) = bjm(1,i)
        rhs(i+4*npts) = bjm(2,i)
        rhs(i+5*npts) = bjm(3,i)
      enddo

      if(igeomtype.ge.2) then
        print *, "Here"

        allocate(bbp_a(3,na),bbp_b(3,nb))
        call fun_surf_interp(3,npatches,norders,ixyzs,iptype,npts,
     1     bbp,na,apatches,auv,bbp_a)
        call fun_surf_interp(3,npatches,norders,ixyzs,iptype,npts,
     1     bbp,nb,bpatches,buv,bbp_b)

        rra = 0 
        do i=1,na
          rhs(6*npts+3) = rhs(6*npts+3) + (bbp_a(1,i)*avals(4,i) + 
     1      bbp_a(2,i)*avals(5,i) + bbp_a(3,i)*avals(6,i))*awts(i)
          rra = rra + 
     1       sqrt(avals(4,i)**2 + avals(5,i)**2 + avals(6,i)**2)*awts(i)
        enddo

        rrb = 0
        do i=1,nb
          rhs(6*npts+4) = rhs(6*npts+4) + (bbp_b(1,i)*bvals(4,i) + 
     1      bbp_b(2,i)*bvals(5,i) + bbp_b(3,i)*bvals(6,i))*bwts(i)
          rrb = rrb + 
     1      sqrt(bvals(4,i)**2 + bvals(5,i)**2 + bvals(6,i)**2)*bwts(i)
        enddo
        rhs(6*npts+3) = rhs(6*npts+3)*0 
        rhs(6*npts+4) = rhs(6*npts+4)*0
        call prin2('proj3=*',rhs(6*npts+3),1)
        call prin2('proj4=*',rhs(6*npts+4),1)
      endif




      eps = 0.51d-7

      allocate(ipatch_id(npts),uvs_targ(2,npts))
      call get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts,
     1  ipatch_id,uvs_targ)
      call prinf('ipatch_id=*',ipatch_id,20)
      call prin2('uvs_targ=*',uvs_targ,48)
      
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

      ikerorder = 0
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

      numit = 150
      allocate(errs(numit+1))
      call prinf('ngenus=*',ngenus,1)
      eps_gmres = 1.0d-8
      call prin2('rhs_projs=*',rhs(6*npts+1),4)
      print *, "here"
      call statj_gendeb_solver(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,eps,dpars,ngenus,hvecs,bbphvecs,na,apatches,
     2  auv,avals,awts,nb,bpatches,buv,bvals,bwts,numit,rhs,eps_gmres,
     3  niter,errs,rres,soln)
      
      call prin2('projs=*',soln(6*npts+1),4)
      bbpc(1:3) = 0
      bbpex(1:3) = 0
      ptmp = 0
      call l3ddirectcdg(1,xyz_in,cf2,vf2,2,xyz_out,1,ptmp,bbpex,thresh)
      ra = 0

      do i=1,npts
        dx = xyz_out(1) - srcvals(1,i)
        dy = xyz_out(2) - srcvals(2,i)
        dz = xyz_out(3) - srcvals(3,i)
        
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
        ra = ra + wts(i)
      enddo

      bbpc(1:3) = bbpc(1:3)/4/pi
      print *, "ra=",ra

      erra = abs(bbpc(1)-bbpex(1)) + abs(bbpc(2)-bbpex(2)) + 
     1   (bbpc(3)-bbpex(3))
      ra = abs(bbpex(1)) + abs(bbpex(2)) + abs(bbpex(3))
      call prin2('bbpex=*',bbpex,3)
      call prin2('bbpc=*',bbpc,3)
      print *, bbpex(1)/bbpc(1)
      print *, bbpex(2)/bbpc(2)
      print *, bbpex(3)/bbpc(3)
      errpt = erra/ra
      call prin2('relativd error in exterior b field=*',erra/ra,1)


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
      ra = 0
      do i=1,npts
        do j=1,3
          errj  = errj + (bjmcomp(j,i) - bjm(j,i))**2*wts(i)
          errbm = errbm + (bbmcomp(j,i) - bbm(j,i))**2*wts(i)
          errbp = errbp + (bbpcomp(j,i) - bbp(j,i))**2*wts(i)
        enddo
        do j=1,6
          ra = ra + rhs((j-1)*npts+i)**2*wts(i)
        enddo
      enddo
      errj = sqrt(errj/ra)
      errbm = sqrt(errbm/ra)
      errbp = sqrt(errbp/ra)
      call prin2('error in current=*',errj,1)
      call prin2('error in interior magnetic field=*',errbm,1)
      call prin2('error in exterior magnetic field=*',errbp,1)

      if(igeomtype.eq.2) then
      write(fname,'(a,i2.2,a,i2.2,a,i1,a,i1,a)') 'stell8_soln_',
     1    ipars(1),'_',ipars(2),'_',norder,'_',ibg,'.dat'
      endif
      if(igeomtype.eq.4) then
      write(fname,'(a,i2.2,a,i2.2,a,i1,a,i1,a)') 'torus_soln_',ipars(1),
     1    '_',ipars(2),'_',norder,'_',ibg,'.dat'
      endif

      if(igeomtype.eq.3) then
      write(fname,'(a,i2.2,a,i2.2,a,i1,a,i1,a)') 'wtorus_soln_',
     1    ipars(1),'_',ipars(2),'_',norder,'_',ibg,'.dat'
      endif
      if(igeomtype.eq.5) then
      write(fname,'(a,i2.2,a,i2.2,a,i1,a,i1,a)') 'torus2_soln_',
     1    ipars(1),'_',ipars(2),'_',norder,'_',ibg,'.dat'
      endif


      open(unit=80,file=fname)
      write(80,*) niter
      write(80,*) rres
 1230 format(5(2x,e11.5))     
      write(80,1230) errest,errpt,errj,errbm,errbp
      do i=1,6*npts+4*ngenus
        write(80,*) soln(i)
      enddo

 1233 format(6(2x,e22.16))
      do i=1,npts
        write(80,1233) bjmcomp(1,i),bjmcomp(2,i),bjmcomp(3,i),bjm(1,i),
     1    bjm(2,i),bjm(3,i)
      enddo

      do i=1,npts
        write(80,1233) bbmcomp(1,i),bbmcomp(2,i),bbmcomp(3,i),bbm(1,i),
     1    bbm(2,i),bbm(3,i)
      enddo

      do i=1,npts
        write(80,1233) bbpcomp(1,i),bbpcomp(2,i),bbpcomp(3,i),bbp(1,i),
     1    bbp(2,i),bbp(3,i)
      enddo

      do i=1,niter
        write(80,*) errs(i)
      enddo
      close(80)

      stop

      

c
c  compute integral of densities
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

      
      dzk = 0
      call statj_gendebproc_rhomrhopmum(npatches,norders,ixyzs, 
     1  iptype,npts,srccoefs,srcvals,eps,nnz,row_ptr,col_ind,iquad, 
     2  nquad,wnear,soln,novers,npts_over,ixyzso,srcover,wover,curv, 
     3  wtmp1,wtmp2,wtmp3,wtmp4,dzk,rbeta,rgamma,laps02rhom,laps02rhop, 
     4  laps02mum,blm0,bmm0)
      
      call surf_vtk_plot_vec(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,bmm,'bmm-torus-110-reft.vtk','a')
      call surf_vtk_plot_vec(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,blm,'blm-torus-110-reft.vtk','a')
      call surf_vtk_plot_vec(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,bmm0,'bmm0-torus-110-reft.vtk','a')
      call surf_vtk_plot_vec(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,blm0,'blm0-torus-110-reft.vtk','a')
      

      call surf_vtk_plot_scalar(npatches,norders,ixyzs,iptype,npts,
     1   srccoefs,srcvals,soln(3*npts+1),'rqm-torus-110-reft.vtk','a')

      call surf_vtk_plot_scalar(npatches,norders,ixyzs,iptype,npts,
     1   srccoefs,srcvals,soln(4*npts+1),'rqp-torus-110-reft.vtk','a')

      call surf_vtk_plot_scalar(npatches,norders,ixyzs,iptype,npts,
     1   srccoefs,srcvals,soln(5*npts+1),'rrm-torus-110-reft.vtk','a')

      call surf_vtk_plot_scalar(npatches,norders,ixyzs,iptype,npts,
     1   srccoefs,srcvals,soln(1),'rrhom-torus-110-reft.vtk','a')

      call surf_vtk_plot_scalar(npatches,norders,ixyzs,iptype,npts,
     1   srccoefs,srcvals,soln(npts+1),'rrhop-torus-110-reft.vtk','a')

      call surf_vtk_plot_scalar(npatches,norders,ixyzs,iptype,npts,
     1   srccoefs,srcvals,soln(2*npts+1),'rmum-torus-110-reft.vtk','a')

      call surf_quadratic_msh_vtk_plot(npatches,norders,ixyzs,iptype,
     1  npts,srccoefs,srcvals,'torus-msh.vtk','a')


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






