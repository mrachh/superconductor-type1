      implicit real *8 (a-h,o-z)
      
      do idzk = 4, 10
        do idatatype = 1,2
          call solve_thinshell(idzk, idatatype)
        enddo
      enddo

      stop
      end


      subroutine solve_thinshell(idzk, idatatype)
      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: wts(:),rsigma(:),rpot(:)
      integer ipars(2)
      integer, allocatable :: norders(:),ixyzs(:),iptype(:)

CCC   for thin shell only 
      real *8, allocatable :: srcvals1(:,:),srccoefs1(:,:)
      real *8, allocatable :: srcvals2(:,:),srccoefs2(:,:)
      integer, allocatable :: norders1(:),ixyzs1(:),iptype1(:)
      integer, allocatable :: norders2(:),ixyzs2(:),iptype2(:)
      
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
      real *8, allocatable :: rhs(:),soln(:),soln1(:),soln2(:)
      real *8, allocatable :: errs(:)
      real *8, allocatable :: errp(:),errp1(:),errp2(:)
      real *8  errm,errm1,errm2
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
      integer iaxyzs(3),ibxyzs(3)
      real *8 vf2(3,2),dpars(3),cf2(2),rzf2(3),dpars2(2)
      real *8 xyz_start(3), dxyz(3)
      complex *16 zf2(3),zk
      complex * 16 zpars(3)
      integer numit,niter,ndims(3)
      character *100 title,dirname
      character *1000 fname,fname1,fname2,fname3,fname4 

CCC   for thin shell only 
      real *8, allocatable :: avals1(:,:),bvals1(:,:)
      real *8, allocatable :: awts1(:),bwts1(:)
      real *8, allocatable :: auv1(:,:),buv1(:,:)
      real *8, allocatable :: avals2(:,:),bvals2(:,:)
      real *8, allocatable :: awts2(:),bwts2(:)
      real *8, allocatable :: auv2(:,:),buv2(:,:)
      integer, allocatable :: apatches1(:),bpatches1(:)
      integer, allocatable :: apatches2(:),bpatches2(:)
      complex * 16 zpars1(3), zpars2(3)
      real *8, allocatable :: rhstmp1(:),outtmp1(:)
      real *8, allocatable :: rhstmp2(:),outtmp2(:)

      real *8, allocatable :: bbp1_a(:,:),bbp2_a(:,:)
      real *8, allocatable :: bbp1_b(:,:),bbp2_b(:,:)
      real *8, allocatable :: bbm1_a(:,:),bbm2_a(:,:)
      real *8, allocatable :: bbm1_b(:,:),bbm2_b(:,:)


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

CCC     for thin-shell only 
      real *8, allocatable :: cms1(:,:),rads1(:),rad_near1(:)
      real *8, allocatable :: cms2(:,:),rads2(:),rad_near2(:)
      real *8, allocatable :: sources1(:,:),targs1(:,:)
      real *8, allocatable :: sources2(:,:),targs2(:,:)
      real *8, allocatable :: wnear1(:),wnear2(:)

CCC     for thin shell only 
      integer, allocatable :: iquad1(:),row_ptr1(:),col_ind1(:)
      integer, allocatable :: iquad2(:),row_ptr2(:),col_ind2(:)
      integer, allocatable :: novers1(:),ixyzso1(:)
      integer, allocatable :: novers2(:),ixyzso2(:),ixyzso_2(:)
      real *8, allocatable :: srcover1(:,:),wover1(:)
      real *8, allocatable :: srcover2(:,:),wover2(:)
      

CCC     for thin shell only 
      real *8, allocatable :: hvecs_shell(:,:,:,:)
      real *8, allocatable :: bbphvecs_shell(:,:,:,:)
      real *8, allocatable :: hvecs_div_shell(:,:,:)
      real *8, allocatable :: hvecs_div2_shell(:,:)
      real *8, allocatable :: wts1(:),wts2(:)


      real *8, allocatable :: hvecs1(:,:,:),hvecs2(:,:,:)
      real *8, allocatable :: bbphvecs1(:,:,:),bbphvecs2(:,:,:)
      real *8, allocatable :: hvecs1_div(:,:),hvecs2_div(:,:)
      real *8, allocatable :: hvecs1_div2(:),hvecs2_div2(:)
CCC 

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

      ntargtest = 1
      ibg = 1
c  
c  The penetration depth is 10^(-idzk) for idzk=0,1,2
c  but for idzk =3, penetration depth is 1/30
c
c  idzk = 4:10, 2^(4-idzk)
c
      dzk = 10**(idzk)
      if(idzk.eq.3) dzk = 30.0d0
      if(idzk.ge.4) dzk = 2**(idzk-4)
c      dzk = 0.1d0
      zk = ima*dzk

c
c  Note for this code, igeomtype can only be 6 
c
c  igeomtype = 6 => thin shell tori
c
      igeomtype = 6
c      iref = 1
      iref = 4


      ipars(1) = 2*2**(iref)
      ipars(2) = 1*2**(iref)
      npatches = 4*ipars(1)*ipars(2)

      dirname = '/mnt/home/mrachh/ceph/' // 
     1     'superconductor-type1-data/thinshell-data/'
      fname = 'thin-torus.vtk'
      ngenus = 1

      xyz_in_src(1) = 2.751d0
      xyz_in_src(2) = 0.002d0
      xyz_in_src(3) = 0.001d0 

      uu = 0.75d0*2*pi
      vv = 0.22d0*2*pi
      rr = 1.37d0 
      xyz_out_src(1) = (rr*cos(uu) + 2)*cos(vv)
      xyz_out_src(2) = (rr*cos(uu) + 2)*sin(vv)
      xyz_out_src(3) = rr*sin(uu) + 10.0d0

      norder = 8
      npols = (norder+1)*(norder+2)/2

      npts = npatches*npols
      call prinf('npatches=*',npatches,1)

      allocate(srcvals(12,npts),srccoefs(9,npts))
      ifplot = 1


      call setup_geom(igeomtype,norder,npatches,ipars, 
     1       srcvals,srccoefs,ifplot,fname)




c
cc      call prin2('srcvals=*',srcvals,12*npts)
cc      call prin2('srccoefs=*',srccoefs,9*npts)

      allocate(norders(npatches),ixyzs(npatches+1),iptype(npatches))

      do i=1,npatches
        norders(i) = norder
        ixyzs(i) = 1 +(i-1)*npols
        iptype(i) = 1
      enddo

      ixyzs(npatches+1) = 1+npols*npatches

c       
      npatches1 = npatches/2 
      npatches2 = npatches-npatches1
      npts1 =  npatches1*npols
      npts2 = npts-npts1 

      allocate(srcvals1(12,npts1),srccoefs1(9,npts1))
      allocate(srcvals2(12,npts2),srccoefs2(9,npts2))

      srcvals1 = srcvals(1:12,1:npts1)
      srccoefs1 = srccoefs(1:9,1:npts1)
      srcvals2 = srcvals(1:12,(npts1+1):npts)
      srccoefs2 = srccoefs(1:9,(npts1+1):npts)

      allocate(norders1(npatches1),ixyzs1(npatches1+1),
     1          iptype1(npatches1))

      do i=1,npatches1
        norders1(i) = norder
        ixyzs1(i) = 1 +(i-1)*npols
        iptype1(i) = 1
      enddo

      ixyzs1(npatches1+1) = 1+npols*npatches1


      call surf_vtk_plot_vec(npatches1,norders1,ixyzs1,
     1       iptype1,npts1,srccoefs1,srcvals1,srcvals1(10:12,:),
     2      'surface1-normals.vtk','a')


            
      allocate(norders2(npatches2),ixyzs2(npatches2+1),
     1        iptype2(npatches2))
        
      do i=1,npatches2
        norders2(i) = norder
        ixyzs2(i) = 1 +(i-1)*npols
        iptype2(i) = 1
      enddo 

      ixyzs2(npatches2+1) = 1+npols*npatches2

      call surf_vtk_plot_vec(npatches2,norders2,ixyzs2,
     1      iptype2,npts2,srccoefs2,srcvals2,srcvals2(10:12,:),
     2     'surface2-normals.vtk','a')
           
    
      allocate(wts(npts))
      call get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,wts)


c
c   set a and b cycle params
c

      npatches1 = npatches/2
      call prinf('npatches1=*',npatches1,1)

      m = 40
      na = 2*ipars(2)*m
      nb = 2*ipars(1)*m
      npts1 = npatches1*npols
      npts2 = npts-npts1

      na1 = na/2 
      na2 = na-na1 

      nb1 = nb/2 
      nb2 = nb-nb1

      allocate(avals1(9,na1),awts1(na1),auv1(2,na1),
     1            apatches1(na1))
      allocate(bvals1(9,nb1),bwts1(nb1),buv1(2,nb1),
     1             bpatches1(nb1))
      print *, "npatches1=",npatches1
      call prinf('nb1=*',nb1,1)

      call get_ab_cycles_torusparam(npatches1,norders1,
     1         ixyzs1,iptype1,npts1,srccoefs1,srcvals1,
     2         ipars,m,na1,avals1,awts1,apatches1,auv1,
     3         nb1,bvals1,bwts1,bpatches1,buv1) 
      call vtk_curv_plot(na1,9,avals1,'acycle1.vtk','a1')
      call vtk_curv_plot(nb1,9,bvals1,'bcycle1.vtk','b1')



      allocate(avals2(9,na2),awts2(na2),auv2(2,na2),
     1             apatches2(na2))
      allocate(bvals2(9,nb2),bwts2(nb2),buv2(2,nb2),
     1             bpatches2(nb2))

      call get_ab_cycles_torusparam(npatches2,norders2,
     1         ixyzs2,iptype2,npts2,srccoefs2,srcvals2,
     2         ipars,m,na2,avals2,awts2,apatches2,auv2,
     3         nb2,bvals2,bwts2,bpatches2,buv2) 
      call vtk_curv_plot(na2,9,avals2,'acycle2.vtk','a2')
      call vtk_curv_plot(nb2,9,bvals2,'bcycle2.vtk','b2')


      allocate(avals(9,na),awts(na),auv(2,na),
     1           apatches(na))
      allocate(bvals(9,nb),bwts(nb),buv(2,nb),
     1           bpatches(nb))
      avals(1:9,1:na1) = avals1
      awts(1:na1) = awts1
      auv(1:2,1:na1) = auv1 
      apatches(1:na1) = apatches1
      avals(1:9,(na1+1):na) = avals2
      awts((na1+1):na) = awts2
      auv(1:2,(na1+1):na) = auv2 
      do i=na1+1,na
        apatches(i) = apatches2(i-na1) + npatches1
      enddo
      

      bvals(1:9,1:nb1) = bvals1
      bwts(1:nb1) = bwts1
      buv(1:2,1:nb1) = buv1 
      bpatches(1:nb1) = bpatches1
      bvals(1:9,(nb1+1):nb) = bvals2
      bwts((nb1+1):nb) = bwts2
      buv(1:2,(nb1+1):nb) = buv2
      do i=nb1+1,nb
        bpatches(i) = bpatches2(i-nb1) + npatches1
      enddo


c
c
c   compute harmonic vector fields
c

c
c   hvecs, bbphvecs need to be of size (3,npts,2,2)
c   hvecs_div(npts,2,2) (can be ignored)
c     hvecs(:,:,:,1) -> corresponds to outer torus
c     hvecs(:,:,:,2) -> inner torus
c
c     hvecs(1:3,:,1,:) -> first harmonic vector field (-y/r^2, x/r^2, 0)
c     hvecs(1:3,:,2,:) -> second harmonic vector field
c                         n\times hvecs(:,:,1,:)
c     
c     v_{1}^{+} = hvecs(:,:,1,1)
c     v_{2}^{+} = hvecs(:,:,2,1)
c     
c     v_{1}^{-} = hvecs(:,:,1,2)
c     v_{2}^{-} = hvecs(:,:,2,2)
c
c     bbphvecs = \nabla \times S_{0} [hvecs]
c
c     hvecs_div = \nabla \times S_{0}[n \times n\times hvecs]
c     
c

      call prinf('npts1=*',npts1,1)
      call prinf('npts2=*',npts2,1)

      allocate(hvecs1(3,npts1,2))
      allocate(hvecs2(3,npts2,2))
      allocate(bbphvecs1(3,npts1,2))
      allocate(bbphvecs2(3,npts2,2))
      allocate(hvecs1_div(npts1,2))
      allocate(hvecs2_div(npts2,2))
      allocate(hvecs1_div2(npts1))
      allocate(hvecs2_div2(npts2))


      hvecs1 = 0
      hvecs2 = 0 
      bbphvecs1 = 0 
      bbphvecs2 = 0 
      hvecs1_div = 0 
      hvecs2_div = 0
      hvecs1_div2 = 0 
      hvecs2_div2 = 0 


CCC         surface 1     CCCC 
      do i=1,npts1
        rr1 = srcvals1(1,i)**2 + srcvals1(2,i)**2
        hvecs1(1,i,1) = -srcvals1(2,i)/rr1
        hvecs1(2,i,1) = srcvals1(1,i)/rr1
        hvecs1(3,i,1) = 0 
        call cross_prod3d(srcvals1(10,i),hvecs1(1:3,i,1),
     1       hvecs1(1:3,i,2))
      enddo 

      call prin2('srcvals1=*',srcvals1(1,1),24)
      call prin2('hvecs1=*',hvecs1(1,1,1),24)
  

      allocate(wts1(npts1))
      call get_qwts(npatches1,norders1,ixyzs1,
     1        iptype1,npts1,srcvals1,wts1)

      call surf_div(npatches1,norders1,ixyzs1,iptype1,npts1, 
     1  srccoefs1,srcvals1,hvecs1(1,1,1),hvecs1_div2(1))

      errest = 0
      do i=1,npts1
        errest = errest + hvecs1_div2(i)**2*wts1(i)
      enddo
      errest = sqrt(errest)
      call prin2('errest1=*',errest,1)



CCC         surface 2     CCCC


      do i=1,npts2
        rr1 = srcvals2(1,i)**2 + srcvals2(2,i)**2
        hvecs2(1,i,1) = -srcvals2(2,i)/rr1
        hvecs2(2,i,1) = srcvals2(1,i)/rr1
        hvecs2(3,i,1) = 0 
        call cross_prod3d(srcvals2(10,i),hvecs2(1:3,i,1),
     1     hvecs2(1:3,i,2))
      enddo 


      allocate(wts2(npts2))
      call get_qwts(npatches2,norders2,ixyzs2,
     1   iptype2,npts2,srcvals2,wts2)

    
      call surf_div(npatches2,norders2,ixyzs2,iptype2,npts2, 
     1  srccoefs2,srcvals2,hvecs2(1,1,1),hvecs2_div2(1))

      errest = 0
      do i=1,npts2
        errest = errest + hvecs2_div2(i)**2*wts2(i)
      enddo
      errest = sqrt(errest)
      call prin2('errest2=*',errest,1)


        
      allocate(rhstmp(npts*3),outtmp(npts*3))
      
      iaxyzs(1) = 1
      iaxyzs(2) = na1 + 1
      iaxyzs(3) = na + 1

      ibxyzs(1) = 1
      ibxyzs(2) = nb1 + 1
      ibxyzs(3) = nb + 1


      allocate(bbp(3,npts),ptmp(npts))
      allocate(bjm(3,npts),bbm(3,npts))
      allocate(zjm(3,npts),zbbm(3,npts))

      ntarg = npts
      allocate(sources(3,npts))
C$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,npts 
        sources(1,i) = srcvals(1,i)
        sources(2,i) = srcvals(2,i)
        sources(3,i) = srcvals(3,i)
      enddo
C$OMP END PARALLEL DO      


      allocate(rhs(6*npts+8*ngenus),soln(6*npts+8*ngenus))
      rhs = 0
      soln = 0
      
      if (idatatype.eq.1) then 
        rhs(6*npts+1) = dzk 
        rhs(6*npts+5) = 1 
      else 
        rhs(6*npts+4) = dzk 
        rhs(6*npts+8) = 1 
      endif

      if(idatatype.eq.1) ifb = 1
      if(idatatype.eq.2) ifb = 0
      if(idzk.lt.10) then
          write(fname,'(a,a,i1,a,i2.2,a,i2.2,a,i1,a,i1,a)') 
     1    trim(dirname),'statj_soln_idzk_',idzk,'_',ipars(1),
     2    '_',ipars(2),'_norder',norder,'_ifb',ifb,'.dat'
      else
          write(fname,'(a,a,i2,a,i2.2,a,i2.2,a,i1,a,i1,a)') 
     1    trim(dirname),'statj_soln_idzk_',idzk,'_',ipars(1),
     2    '_',ipars(2),'_norder',norder,'_ifb',ifb,'.dat'
      endif

      print *, trim(fname)


c
c
c  Our ordering: interior a_+, interior a_-, interior b_+, interior b_-
c                exterior a_+, exterior a_-, exterior b_+, exterior b_-
c   


      eps = 0.51d-7

      allocate(ipatch_id(npts),uvs_src(2,npts))
      call get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts,
     1  ipatch_id,uvs_src)
c      call prinf('ipatch_id=*',ipatch_id,20)
c      call prin2('uvs_src=*',uvs_src,48)

      allocate(cms(3,npatches),rads(npatches),rad_near(npatches))
      call get_centroid_rads(npatches,norders,ixyzs,iptype,npts,
     1     srccoefs,cms,rads)
      
      
c
c       precompute near quadrature correction
c
c
c
c  findnearmem + findnear determines which pairs of sources
c  (thought of as patches) and targets require special
c  quadrature. It's job is to populate row_ptr, col_ind,
c  and nnz (number of non-zero entries).
c
c  Input: surface info: npatches, npts, ixyzs, iptype, norders,
c          srcvals, srccoefs, wts
c
c  Output: row_ptr, col_ind, iquad, nnz, novers, ixyzso, srcover
c          wover (intermediate quantities: rad_near)
c
c
c  Currently these things are computed for \Omega^{+} \cup \Omega^{-}
c
c  We need versions of these quantities for just i) \Omega^{+}, ii) \Omega^{-}
c  and iii) \Omega^{+} \cup \Omega^{-}
c

      
        iptype_avg = floor(sum(iptype)/(npatches+0.0d0))
        norder_avg = floor(sum(norders)/(npatches+0.0d0))

        call get_rfacs(norder_avg,iptype_avg,rfac,rfac0)
c        call prin2('rfac=*',rfac,1)
c        call prin2('rfac0=*',rfac0,1)

C$OMP PARALLEL DO DEFAULT(SHARED) 
        do i=1,npatches
         rad_near(i) = rads(i)*rfac
        enddo
C$OMP END PARALLEL DO     

        call prin2('rad_near=*',rad_near,4)
c
c    find near quadrature correction interactions
c 
        call findnearmem(cms,npatches,rad_near,3,sources,npts,nnz)

       allocate(row_ptr(npts+1),col_ind(nnz))
      
        call findnear(cms,npatches,rad_near,3,sources,npts,row_ptr, 
     1        col_ind)

c
c  iquad is indexing array for quadrature corrections
c  sources as points and targets as points
c
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

       rbeta = 0.0d0
       rgamma = 0.0d0

        dpars(1) = dzk
        dpars(2) = rbeta
        dpars(3) = rgamma
        wnear = 0


CCC     thin shell 

CCC     for surface 1 
        iptype_avg1 = floor(sum(iptype1)/(npatches1+0.0d0))
        norder_avg1 = floor(sum(norders1)/(npatches1+0.0d0))

        call get_rfacs(norder_avg1,iptype_avg1,rfac1,rfac01)
c        call prin2('rfac1=*',rfac1,1)
c        call prin2('rfac01=*',rfac01,1)


        allocate(cms1(3,npatches1),rads1(npatches1),
     1       rad_near1(npatches1))
        call get_centroid_rads(npatches1,norders1,ixyzs1,iptype1,npts1, 
     1     srccoefs1,cms1,rads1)
c        call prin2('cms1=*',cms1,24)
c        call prin2('rads1=*',rads1,12)


        do i=1,npatches1
          rad_near1(i) = rads1(i)*rfac1
        enddo
c        call prin2('rad_near1=*',rad_near1,4)
c
c    find near quadrature correction interactions
        allocate(sources1(3,npts1))
        do i=1,npts1 
          sources1(1,i) = srcvals1(1,i)
          sources1(2,i) = srcvals1(2,i)
          sources1(3,i) = srcvals1(3,i)
        enddo
c 
        call findnearmem(cms1,npatches1,rad_near1,3,
     1             sources1,npts1,nnz1)
c        call prinf('nnz1=*',nnz1,1)

        allocate(row_ptr1(npts1+1),col_ind1(nnz1))

        call findnear(cms1,npatches1,rad_near1,3,sources1,
     1        npts1,row_ptr1,col_ind1)

c
c  iquad is indexing array for quadrature corrections
c  sources as points and targets as points
c

        allocate(iquad1(nnz1+1)) 
        ntarg1 = npts1
        call get_iquad_rsc(npatches1,ixyzs1,ntarg1,nnz1,
     1         row_ptr1,col_ind1,iquad1)
c        call prinf('iquad1=*',iquad1,24)

c
c    estimate oversampling for far-field, and oversample geometry
c

        ikerorder = 0
        allocate(novers1(npatches1),ixyzso1(npatches1+1))

        zpars1 = 0
        ndtarg = 3
        call get_far_order(eps,npatches1,norders1,ixyzs1,iptype1,cms1,
     1    rads1,npts1,srccoefs1,ndtarg,npts1,sources1,ikerorder,zpars1,
     2    nnz1,row_ptr1,col_ind1,rfac1,novers1,ixyzso1)

        npts_over1 = ixyzso1(npatches1+1)-1
c        call prinf('npts_over1=*',npts_over1,1)

        allocate(srcover1(12,npts_over1),wover1(npts_over1))

        call oversample_geom(npatches1,norders1,ixyzs1,
     1   iptype1,npts1,srccoefs1,srcvals1,novers1,ixyzso1,
     2   npts_over1,srcover1)

c        call prin2('srcover1=*',srcover1,24)

        call get_qwts(npatches1,novers1,ixyzso1,iptype1,
     1        npts_over1,srcover1,wover1)
c        call prin2('wover1=*',wover1,24)

c
c   compute near quadrature correction
c
        nquad1 = iquad1(nnz1+1)-1
        allocate(wnear1(10*nquad1))
        do i=1,10*nquad1
            wnear1(i) = 0
        enddo

       call prinf('finished generating near field info*',i,0)
       call prinf('finished generating far field orders*',i,0)
c       call prinf('npts_over1=*',npts_over1,1)
       call prin2('eps=*',eps,1)

CCC     for surface 2 
        iptype_avg2 = floor(sum(iptype2)/(npatches2+0.0d0))
        norder_avg2 = floor(sum(norders2)/(npatches2+0.0d0))

        call get_rfacs(norder_avg2,iptype_avg2,rfac2,rfac02)
c        call prin2('rfac2=*',rfac2,1)
c        call prin2('rfac02=*',rfac02,1)


        allocate(cms2(3,npatches2),rads2(npatches2),
     1       rad_near2(npatches2))
        call get_centroid_rads(npatches2,norders2,ixyzs2,iptype2,npts2, 
     1     srccoefs2,cms2,rads2)
c        call prin2('cms2=*',cms2,12)
c        call prin2('rads2=*',rads2,4)


        do i=1,npatches2
          rad_near2(i) = rads2(i)*rfac2
        enddo
c        call prin2('rad_near2=*',rad_near2,4)
c
c    find near quadrature correction interactions
        allocate(sources2(3,npts2))
        do i=1,npts2 
          sources2(1,i) = srcvals2(1,i)
          sources2(2,i) = srcvals2(2,i)
          sources2(3,i) = srcvals2(3,i)
        enddo
c 
        call findnearmem(cms2,npatches2,rad_near2,3,
     1             sources2,npts2,nnz2)
c        call prinf('nnz2=*',nnz2,1)

        allocate(row_ptr2(npts2+1),col_ind2(nnz2))

        call findnear(cms2,npatches2,rad_near2,3,sources2,
     1        npts2,row_ptr2,col_ind2)

c
c  iquad is indexing array for quadrature corrections
c  sources as points and targets as points
c

        allocate(iquad2(nnz2+1)) 
        ntarg2 = npts2
        call get_iquad_rsc(npatches2,ixyzs2,ntarg2,nnz2,
     1         row_ptr2,col_ind2,iquad2)
c        call prinf('iquad2=*',iquad2,24)

c
c    estimate oversampling for far-field, and oversample geometry
c

        ikerorder = 0
        allocate(novers2(npatches2),ixyzso2(npatches2+1))

        zpars2 = 0
        ndtarg = 3
        call get_far_order(eps,npatches2,norders2,ixyzs2,iptype2,cms2,
     1    rads2,npts2,srccoefs2,ndtarg,npts2,sources2,ikerorder,zpars2,
     2    nnz2,row_ptr2,col_ind2,rfac2,novers2,ixyzso2)

        npts_over2 = ixyzso2(npatches2+1)-1
c        call prinf('npts_over2=*',npts_over2,1)

        allocate(srcover2(12,npts_over2),wover2(npts_over2))

        call oversample_geom(npatches2,norders2,ixyzs2,
     1   iptype2,npts2,srccoefs2,srcvals2,novers2,ixyzso2,
     2   npts_over2,srcover2)

c        call prin2('srcover2=*',srcover2,24)

        call get_qwts(npatches2,novers2,ixyzso2,iptype2,
     1        npts_over2,srcover2,wover2)
c        call prin2('wover2=*',wover2,24)
c        print *, "npts_over=",npts_over
c        print *, "npts_over1=",npts_over1
c        print *, "npts_over2=",npts_over2

c
c   compute near quadrature correction
c
        nquad2 = iquad2(nnz2+1)-1
        allocate(wnear2(10*nquad2))
        do i=1,10*nquad2
            wnear2(i) = 0
        enddo




c

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

c      call prin2('wnear=*',wnear,24)


C
C     surface 1 
C
      call getnearquad_statj_gendeb(npatches1,norders1,
     1      ixyzs1,iptype1,npts1,srccoefs1,srcvals1,
     2      epsquad,dzk,iquadtype,nnz1,row_ptr1,col_ind1,
     3      iquad1,rfac01,nquad1,wnear1)

c      call prin2('wnear1=*',wnear1,24)

C
C     surface 2 
C
      call getnearquad_statj_gendeb(npatches2,norders2,
     1      ixyzs2,iptype2,npts2,srccoefs2,srcvals2,
     2      epsquad,dzk,iquadtype,nnz2,row_ptr2,col_ind2,
     3      iquad2,rfac02,nquad2,wnear2)
c


c
c
c  compute bbphvecs
c
c  hvecs(:,:,:,1) use info for \Omega^{+} and compute
c  bbphevcs(:,:,:,1)
c
c  hvecs(:,:,:,2) use info for \Omega^{-} and compute
c  bbphvecs(:,:,:,2)  
c   
      

      allocate(rhstmp1(npts1*3),outtmp1(npts1*3))
      allocate(rhstmp2(npts2*3),outtmp2(npts2*3))

C
CCCC       surface 1 
C
      hvecs1_div = 0 
      hvecs2_div = 0 
      do igen=1,2*ngenus
c        call prin2('hvecs1=*',hvecs1(1,1,igen),24)
        do i=1,npts1 
          rhstmp1(i) = hvecs1(1,i,igen) 
          rhstmp1(i+npts1) = hvecs1(2,i,igen) 
          rhstmp1(i+2*npts1) = hvecs1(3,i,igen) 
        enddo 
        
        outtmp1 = 0
        call lpcomp_s0curl_addsub(npatches1,norders1,ixyzs1,
     1    iptype1,npts1,srccoefs1,srcvals1,eps,nnz1,row_ptr1,
     2    col_ind1,iquad1,nquad1,wnear1(2*nquad1+1),rhstmp1,
     3    novers1,npts_over1,ixyzso1,srcover1,wover1,outtmp1)
          
        do i=1,npts1
          if(igen.eq.1) then
            write(75,'(3(2x,e11.5))') outtmp1(i),
     1           outtmp1(i+npts1), outtmp1(i+2*npts1)
          endif
          bbphvecs1(1,i,igen) = outtmp1(i)
          bbphvecs1(2,i,igen) = outtmp1(i+npts1)
          bbphvecs1(3,i,igen) = outtmp1(i+2*npts1)
        enddo

         
        
        vtmp1 = 0 
        do j=1,npts1
          call cross_prod3d(srcvals1(10,j),hvecs1(1,j,igen),
     1         vtmp1)
          bbphvecs1(1:3,j,igen) = bbphvecs1(1:3,j,igen) 
     1                - vtmp1(1:3)/2
        enddo


c        call prin2('outtmp1=*',outtmp1,24)
c        call prin2('rhstmp1=*',rhstmp1,24)

      enddo 


c      call prin2('bbphvecs1=*',bbphvecs1(1,1,1),24)
c      call prin2('bbphvecs1=*',bbphvecs1(1,1,2),24)
c

C
CCCC       surface 2 
C
      do igen=1,2*ngenus
        do i=1,npts2 
          rhstmp2(i) = hvecs2(1,i,igen) 
          rhstmp2(i+npts2) = hvecs2(2,i,igen) 
          rhstmp2(i+2*npts2) = hvecs2(3,i,igen) 
        enddo 
        
        outtmp2 = 0
        call lpcomp_s0curl_addsub(npatches2,norders2,ixyzs2,
     1    iptype2,npts2,srccoefs2,srcvals2,eps,nnz2,row_ptr2,
     2    col_ind2,iquad2,nquad2,wnear2(2*nquad2+1),rhstmp2,
     3    novers2,npts_over2,ixyzso2,srcover2,wover2,outtmp2)
          
        do i=1,npts2
          bbphvecs2(1,i,igen) = outtmp2(i)
          bbphvecs2(2,i,igen) = outtmp2(i+npts2)
          bbphvecs2(3,i,igen) = outtmp2(i+2*npts2)
        enddo 
        
        vtmp1 = 0 
        do j=1,npts2
          call cross_prod3d(srcvals2(10,j),hvecs2(1,j,igen),
     1         vtmp1)
          bbphvecs2(1:3,j,igen) = bbphvecs2(1:3,j,igen) 
     1                - vtmp1(1:3)/2
        enddo

c        call prin2('outtmp2=*',outtmp2,24)
c        call prin2('rhstmp2=*',rhstmp2,24)
      enddo 



      print *, "here"

      numit = 300
      allocate(errs(numit+1))
      call prinf('ngenus=*',ngenus,1)
      eps_gmres = 1.0d-8
      if(ngenus.ge.1)  call prin2('rhs_projs=*',rhs(6*npts+1),8)
      print *, "here"

c
c
c  Yuguan: Make sure calling sequence is consistent
c
c

      call cpu_time(t1)
C$       t1 = omp_get_wtime()      
      call statj_gendeb_solver_thinshell_guru(npatches,npatches1,
     1 norders,ixyzs,iptype,npts,srccoefs,srcvals,eps,dpars,hvecs1,
     2 hvecs2,bbphvecs1,bbphvecs2,na,na1,iaxyzs,apatches,auv,avals,
     3 awts,nb,nb1,ibxyzs,bpatches,buv,bvals,bwts,nnz,row_ptr,col_ind,
     4 iquad,nquad,wnear,nnz1,npts1,row_ptr1,col_ind1,iquad1,nquad1,
     5 wnear1,nnz2,npts2,row_ptr2,col_ind2,iquad2,nquad2,wnear2,rfac,
     6 numit,rhs,eps_gmres,niter,errs,rres,soln,cms,rads)
      call cpu_time(t2)
C$       t2 = omp_get_wtime()      
      tsolve = t2-t1 
      call prin2('solve time=*',t2-t1,1)      
      
      call prin2('projs=*',soln(6*npts+1),8)




      allocate(bbpcomp(3,npts),bjmcomp(3,npts),bbmcomp(3,npts))

      print *, "here"
      call prin2('eps=*',eps,1)
      call prin2('soln=*',soln,24)
      call prinf('ngenus=*',ngenus,1)
      call prin2('dpars=*',dpars,3)


c        call lpcomp_statj_gendeb_postproc(npatches,norders,ixyzs,
c     1    iptype,npts,srccoefs,srcvals,eps,dpars,nnz,row_ptr,col_ind,
c     2    iquad,nquad,wnear,ngenus,hvecs,bbphvecs,soln,novers,npts_over,
c     3    ixyzso,srcover,wover,bjmcomp,bbmcomp,bbpcomp)


      call lpcomp_statj_gendeb_postproc_thinshell(npatches,
     1   npatches1,norders,ixyzs,iptype,npts,srccoefs,
     1   srcvals,eps,dpars,nnz,row_ptr,col_ind,iquad,nquad,wnear,
     1   npts1,nnz1,row_ptr1,col_ind1,iquad1,nquad1,wnear1,npts2,
     1   nnz2,row_ptr2,col_ind2,iquad2,nquad2,wnear2,ngenus,hvecs1,
     1   hvecs2,bbphvecs1,bbphvecs2,soln,novers,npts_over,ixyzso,
     1   srcover,wover,bjmcomp,bbmcomp,bbpcomp)

      do i=1,npts
        bbmcomp(1:3,i) = bbmcomp(1:3,i)/dzk
      enddo

      call prin2('bbmcomp=*',bbmcomp,24)
      call prin2('bbpcomp=*',bbpcomp,24)
      open(unit=80,file=trim(fname),form='unformatted')
      write(80) ipars
      write(80) iref
      write(80) ifb
      write(80) tsolve
      write(80) niter
      write(80) rres
      write(80) rhs
      write(80) soln
      write(80) bjmcomp
      write(80) bbmcomp
      write(80) bbpcomp

      close(80)

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






