      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: wts(:),rsigma(:)
      integer ipars(2)

      integer, allocatable :: norders(:),ixyzs(:),iptype(:)

      real *8 xyz_out(3),xyz_in(3,2)
      real *8 vtmp1(3),bbpc(3),bbpex(3)
      real *8, allocatable :: ffform(:,:,:),ffformex(:,:,:)
      real *8, allocatable :: ffforminv(:,:,:),ffformexinv(:,:,:)

      complex *16, allocatable :: zbbm(:,:),zbbp(:,:),zjm(:,:)
      real *8, allocatable :: bbm(:,:),bbp(:,:),jm(:,:)
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
      real *8, allocatable :: hvecs_div(:,:)
      real *8, allocatable :: hvecs_a(:,:,:),bbphvecs_a(:,:,:)
      real *8, allocatable :: hvecs_b(:,:,:),bbphvecs_b(:,:,:)
      integer, allocatable :: apatches(:),bpatches(:)
      real *8 vf2(3,2),dpars(2),cf2(2)
      complex * 16 zpars(3)
      integer numit,niter
      character *100 title,dirname
      character *300 fname

      real *8 wtmp1(3),wtmp2(3)
      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:)
      real *8, allocatable :: wnear(:)
      real *8, allocatable :: targs(:,:)
      real *8, allocatable :: cms(:,:),rads(:),rad_near(:)
      integer, allocatable :: iquad(:),row_ptr(:),col_ind(:)
      integer, allocatable :: novers(:),ixyzso(:)
      real *8, allocatable :: srcover(:,:),wover(:)
      real *8, allocatable :: ptmp(:)
      real *8 rinttmp(6)

      logical isout0,isout1

      complex *16 pot,potex,ztmp,ima

      data ima/(0.0d0,1.0d0)/


      call prini(6,13)

      done = 1
      pi = atan(done)*4



      igeomtype = 4
      ipars(1) = 4*4
      ipars(2) = 4*4
      npatches = 2*ipars(1)*ipars(2)
      fname = 'torus.vtk'


      xyz_in(1,1) = 3.24d0
      xyz_in(2,1) = 0.002d0
      xyz_in(3,1) = 0.001d0

      xyz_in(1,2) = 3.26d0
      xyz_in(2,2) = -0.01d0
      xyz_in(3,2) = 0.001d0

      xyz_out(1) = -3.5d0
      xyz_out(2) = 7.1d0
      xyz_out(3) = 5.7d0


      norder = 5 
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
 
      m = 32
      na = ipars(2)*m
      nb = ipars(1)*m
      allocate(avals(9,na),awts(na),auv(2,na),apatches(na))
      allocate(bvals(9,nb),bwts(nb),buv(2,nb),bpatches(nb))
      call get_ab_cycles_torusparam(npatches,norders,ixyzs,iptype,
     1   npts,srccoefs,srcvals,ipars,m,na,avals,awts,apatches,auv,
     2   nb,bvals,bwts,bpatches,buv)
      call prin2('bvals=*',bvals(6,1:nb),nb)

      allocate(hvecs(3,npts,2),bbphvecs(3,npts,2),hvecs_div(npts,2))
      allocate(rhstmp(npts*3),outtmp(npts*3))
      bbphvecs = 0

      do i=1,npts
        rr1 = srcvals(1,i)**2 + srcvals(2,i)**2
        hvecs(1,i,1) = -srcvals(2,i)/rr1
        hvecs(2,i,1) = srcvals(1,i)/rr1
        hvecs(3,i,1) = 0 
        call cross_prod3d(srcvals(10,i),hvecs(1:3,i,1),hvecs(1:3,i,2))
      enddo
      call surf_div(npatches,norders,ixyzs,iptype,npts,srccoefs,
     1  srcvals,hvecs(1,1,1),hvecs_div(1,1))
      ra = 0
      do i=1,npts
        ra = ra + hvecs_div(i,1)**2*wts(i)
      enddo
      ra = sqrt(ra) 
      call prin2('error in surf div hvecs_div=*',ra,1)
c
c
c   compute the boundary data, for now assume that only external
c   B field is applied and that the interior fields are 0 
c
      vf2(1,1) = hkrand(0)
      vf2(2,1) = hkrand(0)
      vf2(3,1) = hkrand(0)

      vf2(1:3,1) = vf2(1:3,1)*1.0d1*0 
      vf2(1:3,2) = 0
      cf2(1) = 1
      cf2(2) = -1*0 

      

      thresh = 1.0d-16
      allocate(bbp(3,npts),ptmp(npts))
      ptmp = 0
      bbp = 0
      call l3ddirectcdg(1,xyz_in,cf2,vf2,2,srcvals(1:3,1:npts),npts,
     1 ptmp,bbp,thresh)
      
      call prin2('ptmp=*',ptmp,24)
      allocate(rhs(6*npts+4),soln(6*npts+4))
      rhs = 0
      soln = 0
      do i=1,npts
        rhs(i) = -bbp(1,i)
        rhs(i+npts) = -bbp(2,i)
        rhs(i+2*npts) = -bbp(3,i)
      enddo

      allocate(bbp_a(3,na),bbp_b(3,nb))
      call fun_surf_interp(3,npatches,norders,ixyzs,iptype,npts,
     1   bbp,na,apatches,auv,bbp_a)
      call fun_surf_interp(3,npatches,norders,ixyzs,iptype,npts,
     1   bbp,nb,bpatches,buv,bbp_b)

      
      do i=1,na
        rhs(6*npts+3) = rhs(6*npts+3) + (bbp_a(1,i)*avals(4,i) + 
     1    bbp_a(2,i)*avals(5,i) + bbp_a(3,i)*avals(6,i))*awts(i)
      enddo

      do i=1,nb
        rhs(6*npts+4) = rhs(6*npts+4) + (bbp_b(1,i)*bvals(4,i) + 
     1    bbp_b(2,i)*bvals(5,i) + bbp_b(3,i)*bvals(6,i))*bwts(i)
      enddo
      call prin2('proj3=*',rhs(6*npts+3),1)
      call prin2('proj4=*',rhs(6*npts+4),1)
      rhs(6*npts+3) = 0 
      rhs(6*npts+4) = 0

c
c   test hvecs by computing the surface divergence at a point
c   in the interior and exterior for the first density
c
      rint = 0
      rext = 0
      do i=1,npts
        dx = xyz_in(1,1) - srcvals(1,i)
        dy = xyz_in(2,1) - srcvals(2,i)
        dz = xyz_in(3,1) - srcvals(3,i)
        
        rr = sqrt(dx**2 + dy**2 + dz**2)
        rint = rint + dx*hvecs(1,i,1)*wts(i)/rr**3
        rint = rint + dy*hvecs(2,i,1)*wts(i)/rr**3
        rint = rint + dz*hvecs(3,i,1)*wts(i)/rr**3

        dx = xyz_out(1) - srcvals(1,i)
        dy = xyz_out(2) - srcvals(2,i)
        dz = xyz_out(3) - srcvals(3,i)
        
        rr = sqrt(dx**2 + dy**2 + dz**2)
        rext = rext + dx*hvecs(1,i,1)*wts(i)/rr**3
        rext = rext + dy*hvecs(2,i,1)*wts(i)/rr**3
        rext = rext + dz*hvecs(3,i,1)*wts(i)/rr**3
      enddo

      call prin2('rint=*',rint,1)
      call prin2('rext=*',rext,1)



c
c
c  compute bbphvecs
c


      eps = 0.51d-5

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
      call prin2('rads=*',rads,24)

C$OMP PARALLEL DO DEFAULT(SHARED) 
      do i=1,npatches
        rad_near(i) = rads(i)*rfac
      enddo
C$OMP END PARALLEL DO     

      call prin2('rad_near=*',rad_near,24)
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
      dzk = 1.0d0
      wnear = 0


      call getnearquad_statj_gendeb(npatches,norders,
     1      ixyzs,iptype,npts,srccoefs,srcvals,
     1      eps,dzk,iquadtype,nnz,row_ptr,col_ind,
     1      iquad,rfac0,nquad,wnear)
      call prinf('finished generating near quadrature correction*',i,0)

      call prinf('entering layer potential eval*',i,0)
      call prinf('npts=*',npts,1)

      ra = 0
      do i=1,npts
        ra = ra + wts(i)
      enddo
      call prin2('surface area of torus=*',ra,1)

      ra2 = 0
      do i=1,npts_over
        ra2 = ra2 + wover(i)
      enddo
      call prin2('surface area of torus with oversampled points=*',
     1   ra2,1)
      call prin2('error=*',abs(ra2-ra),1)

      ngenus = 1
      hvecs_div = 0
      do igen=1,2*ngenus
        do i=1,npts
          rhstmp(i) = hvecs(1,i,igen) 
          rhstmp(i+npts) = hvecs(2,i,igen) 
          rhstmp(i+2*npts) = hvecs(3,i,igen) 
        enddo
        call prin2('rhstmp=*',rhstmp,24)
        outtmp = 0
cc        call prinf('nnz=*',nnz,1)
cc        call prinf('row_ptr=*',row_ptr,20)
cc        call prinf('col_ind=*',col_ind,100)

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
cc        print *, "Here"
        do j=1,npts
          call cross_prod3d(srcvals(10,j),hvecs(1,j,igen),vtmp1)
          bbphvecs(1:3,j,igen) = bbphvecs(1:3,j,igen) - vtmp1(1:3)/2
          write(47+i,'(3(2x,e11.5))') bbphvecs(1,j,igen),
     1       bbphvecs(2,j,igen),bbphvecs(3,j,igen)
          write(49+i,'(3(2x,e11.5))') hvecs(1,j,igen),
     1       hvecs(2,j,igen),hvecs(3,j,igen)
          write(51+i,'(3(2x,e11.5))') vtmp1(1),vtmp1(2),vtmp1(3)
        enddo

        ra = 0
        do j=1,npts
          ra = ra + hvecs_div(j,igen)**2*wts(j)
        enddo
cc        call prin2('hvecs_div=*',hvecs_div(1,igen),24)
        ra = sqrt(ra)
        call prin2('error in div s0(hvec)=*',ra,1)
      enddo

c
c
c  Check that integral of bbphvecs on a and b cycles
c  must agree with -\int_{} n\times hvec
c
c
      allocate(bbphvecs_a(3,na,2),bbphvecs_b(3,nb,2))
      allocate(hvecs_a(3,na,2),hvecs_b(3,nb,2))
      do igen=1,2*ngenus
        call fun_surf_interp(3,npatches,norders,ixyzs,iptype,npts,
     1   hvecs(1,1,igen),na,apatches,auv,hvecs_a(1,1,igen))
        call fun_surf_interp(3,npatches,norders,ixyzs,iptype,npts,
     1   bbphvecs(1,1,igen),na,apatches,auv,bbphvecs_a(1,1,igen))
        call fun_surf_interp(3,npatches,norders,ixyzs,iptype,npts,
     1   hvecs(1,1,igen),nb,bpatches,buv,hvecs_b(1,1,igen))
        call fun_surf_interp(3,npatches,norders,ixyzs,iptype,npts,
     1   bbphvecs(1,1,igen),nb,bpatches,buv,bbphvecs_b(1,1,igen))

        ra1 = 0
        ra2 = 0
        ra3 = 0
        ra4 = 0
        do i=1,na
          call cross_prod3d(avals(7,i),hvecs_a(1,i,igen),vtmp1)
          ra1 = ra1 - (vtmp1(1)*avals(4,i)+vtmp1(2)*avals(5,i)+
     1       vtmp1(3)*avals(6,i))*awts(i)
          ra2 = ra2 + (bbphvecs_a(1,i,igen)*avals(4,i)+
     1      bbphvecs_a(2,i,igen)*avals(5,i) +
     2      bbphvecs_a(3,i,igen)*avals(6,i))*awts(i)
        enddo

        do i=1,nb
          call cross_prod3d(bvals(7,i),hvecs_b(1,i,igen),vtmp1)
          ra3 = ra3 - (vtmp1(1)*bvals(4,i)+vtmp1(2)*bvals(5,i)+
     1       vtmp1(3)*bvals(6,i))*bwts(i)
          ra4 = ra4 + (bbphvecs_b(1,i,igen)*bvals(4,i)+
     1      bbphvecs_b(2,i,igen)*bvals(5,i) +
     2      bbphvecs_b(3,i,igen)*bvals(6,i))*bwts(i)
        enddo
        call prinf('igen=*',igen,1)
        call prin2('ra1=*',ra1,1)
        call prin2('ra2=*',ra2,1)
        call prin2('ra3=*',ra3,1)
        call prin2('ra4=*',ra4,1)
      enddo
      numit = 200
      allocate(errs(numit+1))
      call prinf('ngenus=*',ngenus,1)
      eps_gmres = 1.0d-9
      call statj_gendeb_solver(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,eps,dzk,ngenus,hvecs,bbphvecs,na,apatches,auv,
     2  avals,awts,nb,bpatches,buv,bvals,bwts,numit,rhs,eps_gmres,
     3  niter,errs,rres,soln)
      

      call prin2('projs=*',soln(6*npts+1),4)
      bbpc(1:3) = 0
      bbpex(1:3) = 0
      call l3ddirectcdg(1,xyz_in,cf2,vf2,2,xyz_out,1,ptmp,bbpex,thresh)

      c0 = soln(6*npts+3)
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
        
          bbpc(1) = bbpc(1) -1.0d0/r**3*wts(i)*soln(6*npts+2+igen)*
     1        (dy*hvecs(3,i,igen)-dz*hvecs(2,i,igen))
          bbpc(2) = bbpc(2) -1.0d0/r**3*wts(i)*soln(6*npts+2+igen)*
     1        (dz*hvecs(1,i,igen)-dx*hvecs(3,i,igen))
          bbpc(3) = bbpc(3) -1.0d0/r**3*wts(i)*soln(6*npts+2+igen)*
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
      call prin2('relativd error in exterior b field=*',erra/ra,1)

      dpars(1) = 1.0d0
      dpars(2) = 0.0d0

      rinttmp(1:6) = 0
      do i=1,npts
        do j=1,6
          rinttmp(j) = rinttmp(j) + soln((j-1)*npts+i)*wts(i)
        enddo
      enddo
      call prin2('integral of densities=*',rinttmp,6)

      ptmp = 0
      call prin2('dpars=*',dpars,2)
      call prin2('eps=*',eps,1)

cc      call lpcomp_lap_comb_dir_addsub(npatches,norders,ixyzs,iptype,
cc     1    npts,srccoefs,srcvals,12,npts,srcvals,eps,dpars,nnz,row_ptr,
cc     2    col_ind,iquad,nquad,wnear,soln(2*npts+1),
cc     3    novers,npts_over,ixyzso,srcover,wover,ptmp)

      call prin2('soln=*',soln(2*npts+1),24)
      call lpcomp_lap_comb_dir(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,12,npts,srcvals,ipatch_id,uvs_Targ,eps,dpars,
     2  soln(2*npts+1),ptmp)
      call prin2('ptmp=*',ptmp,24)
      rint = 0
      do i=1,npts
        rint = rint + ptmp(i)*wts(i)
      enddo
      print *, "rint=",rint
      call prin2('integral of S of third density=*',rint,1)
       

      ptmp = 0
      call prin2('soln=*',soln,24)
      call lpcomp_lap_comb_dir(npatches,norders,ixyzs,iptype,npts,
     1  srccoefs,srcvals,12,npts,srcvals,ipatch_id,uvs_Targ,eps,dpars,
     2  soln,ptmp)
      call prin2('ptmp=*',ptmp,24)
      rint = 0
      do i=1,npts
        rint = rint + ptmp(i)*wts(i)
      enddo
      call prin2('integral of S of first density=*',rint,1)

      do i=1,npts
        write(70,'(6(2x,e11.5))') soln(i),soln(npts+i),soln(2*npts+i),
     1    soln(3*npts+i),soln(4*npts+i),soln(5*npts+i)

      enddo
      stop
      
      allocate(bbpcomp(3,npts),bjmcomp(3,npts),bbmcomp(3,npts))
      ngenus = 1
      
cc      call lpcomp_statj_gendeb_postproc_addsub(npatches,norders,
cc     1  ixyzs,iptype,npts,srccoefs,srcvals,eps,dzk,nnz,row_ptr,
cc     2  col_ind,iquad,nquad,wnear,ngenus,hvecs,bbphvecs,soln,
cc     3  novers,npts_over,ixyzso,srcover,wover,bjmcomp,bbmcomp,
cc     4  bbpcomp)

      errbm = 0
      ra = 0.0d0
      do i=1,npts
        do j=1,6
          ra = ra + rhs(i+(j-1)*npts)**2*wts(i)
        enddo
      enddo
      do i =1,npts
        errbm = errbm + bbmcomp(1,i)**2*wts(i) 
        errbm = errbm + bbmcomp(2,i)**2*wts(i) 
        errbm = errbm + bbmcomp(3,i)**2*wts(i) 
      enddo
      errbm = sqrt(errbm/ra)
      call prin2('error in bm=*',errbm,1)

       
      errjm = 0
      do i =1,npts
        errjm = errjm + bjmcomp(1,i)**2*wts(i) 
        errjm = errjm + bjmcomp(2,i)**2*wts(i) 
        errjm = errjm + bjmcomp(3,i)**2*wts(i) 
      enddo
      errjm = sqrt(errjm/ra)
      call prin2('error in jm=*',errjm,1)


     

      

      errbp = 0
      rbp = 0
      do i =1,npts
        rbp = rbp + bbpcomp(1,i)**2*wts(i) 
        rbp = rbp + bbpcomp(2,i)**2*wts(i) 
        rbp = rbp + bbpcomp(3,i)**2*wts(i) 
        errbp = errbp + (bbpcomp(1,i)-bbp(1,i))**2*wts(i) 
        errbp = errbp + (bbpcomp(2,i)-bbp(2,i))**2*wts(i) 
        errbp = errbp + (bbpcomp(3,i)-bbp(3,i))**2*wts(i) 
      enddo
      errbp = sqrt(errbp/ra)
      call prin2('error in bp=*',errbp,1)

      call prin2('bbp=*',bbp,24)
      call prin2('bbpcomp=*',bbpcomp,24)
      call prin2('ra=*',ra,1)





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






