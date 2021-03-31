      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: wts(:),rsigma(:)
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
      real *8 vf2(3,2),dpars(3),cf2(2),rzf2(3)
      complex *16 zf2(3),zk
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
c
c
c  ntarg test must be less than 100
c

      ntargtest = 10

      ibg = 3

c
c  Note for this code, igeomtype can only be 1,3,or, 4 
c
      igeomtype = 3
      if(igeomtype.eq.4) then
        ipars(1) = 4*4
        ipars(2) = 2*4
        npatches = 2*ipars(1)*ipars(2)
        fname = 'torus.vtk'

        xyz_in_src(1) = 2.001d0
        xyz_in_src(2) = 0.002d0
        xyz_in_src(3) = 0.001d0

        uu = 0.75d0*2*pi
        vv = 0.22d0*2*pi
        rr = 1.37d0 
        xyz_out_src(1) = (rr*cos(uu) + 2)*cos(vv)
        xyz_out_src(2) = (rr*cos(uu) + 2)*sin(vv)
        xyz_out_src(3) = rr*sin(uu)

        do i=1,ntargtest
          uu = hkrand(0)*2*pi
          vv = hkrand(0)*2*pi
          rr = 0.67d0*hkrand(0)
          xyz_in_targ(1,i) = (rr*cos(uu) + 2)*cos(vv)
          xyz_in_targ(2,i) = (rr*cos(uu) + 2)*sin(vv)
          xyz_in_targ(3,i) = rr*sin(uu)
        enddo

        do i=1,ntargtest
          uu = hkrand(0)*2*pi
          vv = hkrand(0)*2*pi
          rr = 0.67d0*hkrand(0) + 1.33d0
          xyz_out_targ(1,i) = (rr*cos(uu) + 2)*cos(vv)
          xyz_out_targ(2,i) = (rr*cos(uu) + 2)*sin(vv)
          xyz_out_targ(3,i) = rr*sin(uu)
        enddo


        ngenus = 1

      endif

      if(igeomtype.eq.3) then
        ipars(1) = 8*8
        ipars(2) = 8*4
        npatches = 2*ipars(1)*ipars(2)
        
        fname = 'wtorus.vtk'

        xyz_in_src(1) = -2.001d0
        xyz_in_src(2) = 0.002d0
        xyz_in_src(3) = 0.001d0

        uu = 0.91d0*pi
        vv = 0.17d0*2*pi
        rr = 1.33d0 
        xyz_out_src(1) = (rr*cos(uu) + 2+0.25d0*cos(3*vv))*cos(vv)
        xyz_out_src(2) = (rr*cos(uu) + 2+0.25d0*cos(3*vv))*sin(vv)
        xyz_out_src(3) = rr*sin(uu)


        do i=1,ntargtest
          uu = hkrand(0)*2*pi
          vv = hkrand(0)*2*pi
          rr = 0.67d0*hkrand(0)
          xyz_in_targ(1,i) = (rr*cos(uu) + 2+0.25d0*cos(3*vv))*cos(vv)
          xyz_in_targ(2,i) = (rr*cos(uu) + 2+0.25d0*cos(3*vv))*sin(vv)
          xyz_in_targ(3,i) = rr*sin(uu)
        enddo

        do i=1,ntargtest
          uu = hkrand(0)*2*pi
          vv = hkrand(0)*2*pi
          rr = 0.67d0*hkrand(0) + 1.33d0
          xyz_out_targ(1,i) = (rr*cos(uu) + 2+0.25d0*cos(3*vv))*cos(vv)
          xyz_out_targ(2,i) = (rr*cos(uu) + 2+0.25d0*cos(3*vv))*sin(vv)
          xyz_out_targ(3,i) = rr*sin(uu)
        enddo


        ngenus = 1

      endif

      if(igeomtype.eq.1) then
        ipars(1) = 3
        npatches = 12*(4**(ipars(1)))
        fname = 'sphere.vtk'

        xyz_in_src(1) = 0.201d0
        xyz_in_src(2) = 0.102d0
        xyz_in_src(3) = 0.011d0

        xyz_out_src(1) = -1.2d0
        xyz_out_src(2) = 1.1d0
        xyz_out_src(3) = 1.7d0

        do i=1,ntargtest
          rr = hkrand(0)*0.67d0
          thet = hkrand(0)*pi
          phi = hkrand(0)*2*pi
          xyz_in_targ(1,i) = rr*sin(thet)*cos(phi)
          xyz_in_targ(2,i) = rr*sin(thet)*sin(phi)
          xyz_in_targ(3,i) = rr*cos(thet)
        enddo

        do i=1,ntargtest
          rr = 1.33d0 + hkrand(0)*0.67d0
          thet = hkrand(0)*pi
          phi = hkrand(0)*2*pi
          xyz_out_targ(1,i) = rr*sin(thet)*cos(phi)
          xyz_out_targ(2,i) = rr*sin(thet)*sin(phi)
          xyz_out_targ(3,i) = rr*cos(thet)
        enddo


        ngenus = 0
      endif



      norder = 5 
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

      call get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,wts)
      do i=1,ntargtest
        isout0 = .false.
        call test_exterior_pt(npatches,norders,npts,
     1    srcvals,srccoefs,wts,xyz_in_targ(1,i),isout0)
          print *, i,isout0
      enddo

      do i=1,ntargtest
        isout0 = .false.
        call test_exterior_pt(npatches,norders,npts,
     1    srcvals,srccoefs,wts,xyz_out_targ(1,i),isout0)
          print *, i,isout0
      enddo

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

        ifread = 0
        ifwrite = 1
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
      vf2(1:3,1) = 0 
      vf2(1:3,2) = 0
      cf2(1) = 1*1.0d0
      cf2(2) = -1*0 

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
      allocate(targs(3,npts))
C$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,npts 
        targs(1,i) = srcvals(1,i)
        targs(2,i) = srcvals(2,i)
        targs(3,i) = srcvals(3,i)
      enddo
C$OMP END PARALLEL DO      

      call prin2('thresh=*',thresh,1)
      call prin2('xyz_in_src=*',xyz_in_src,3)
      call prin2('cf2=*',cf2,1)
      call prin2('vf2=*',vf2,3)
      call l3ddirectcdg(1,xyz_in_src,cf2,vf2,1,targs,npts,
     1 ptmp,bbp,thresh)
      call prin2('bbp=*',bbp,24)
      dzk = 1.0d0

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
      if(ngenus.ge.1)  call prin2('rhs_projs=*',rhs(6*npts+1),4)
      print *, "here"
      call statj_gendeb_solver(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,eps,dpars,ngenus,hvecs,bbphvecs,na,apatches,
     2  auv,avals,awts,nb,bpatches,buv,bvals,bwts,numit,rhs,eps_gmres,
     3  niter,errs,rres,soln)
      
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
      enddo
      errbdry = sqrt((errj + errbm + errbp)/ra2)



      errj = sqrt(errj/ra)
      errbm = sqrt(errbm/ra)
      errbp = sqrt(errbp/ra)
      call prin2('error in current=*',errj,1)
      call prin2('error in interior magnetic field=*',errbm,1)
      call prin2('error in exterior magnetic field=*',errbp,1)


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

        erra = abs(bbpc(1)-bbpex(1))**2 + abs(bbpc(2)-bbpex(2))**2 + 
     1   (bbpc(3)-bbpex(3))**2
        ra = abs(bbpex(1))**2 + abs(bbpex(2))**2 + abs(bbpex(3))**2
      enddo


      errextsq = erra
      rextsq = ra
      errpt = sqrt(erra/ra)
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
        erra = abs(bbmc(1)-bbmex(1))**2 + abs(bbmc(2)-bbmex(2))**2 + 
     1   (bbmc(3)-bbmex(3))**2
        ra = abs(bbmex(1))**2 + abs(bbmex(2))**2 + abs(bbmex(3))**2

        erra = abs(bjmc(1)-bjmex(1))**2 + abs(bjmc(2)-bjmex(2))**2 + 
     1   (bjmc(3)-bjmex(3))**2
        ra = abs(bjmex(1))**2 + abs(bjmex(2))**2 + abs(bjmex(3))**2
      enddo
      erra = sqrt(erra/ra)
      errvol = erra
      call prin2('l2 error at interior + exterior targets=*',erra,1)
      open(unit=81,file='paper_res.txt',access='append')
c
c  igeomtype,norder,npatches,niter,eps,eps_gmres,errbdry,errvol
c  

 1211 format(2x,i1,2x,i2,2x,i6,2x,i2,4(2x,e11.5)) 
      write(81,1211) igeomtype,norder,npatches,niter,eps,eps_gmres,
     1   errbdry,errvol


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






