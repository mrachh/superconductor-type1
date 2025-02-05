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
      real *8 vtmp1(3),bbpc(3),bbpex(3), vtmp2(3)
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
      character *1000 fname,fname1

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
      idzk = 7
      dzk = 10**(idzk)
      if(idzk.eq.3) dzk = 30.0d0
      if(idzk.ge.4) dzk = 2**(idzk-4)
      dzk = 0.1d0
      zk = ima*dzk

c
c  Note for this code, igeomtype can only be 6 
c
c  igeomtype = 6 => thin shell tori
c
      igeomtype = 6
c      iref = 1
      iref = 2


      ipars(1) = 2*2**(iref)
      ipars(2) = 1*2**(iref)
      npatches = 4*ipars(1)*ipars(2)

      dirname = '/Users/mrachh/git/' // 
     1   'superconductor-type1-data/thin-shell-torus-data/'
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

      thet = 0
CCC         surface 1     CCCC 
      do i=1,npts1
        rr1 = srcvals1(1,i)**2 + srcvals1(2,i)**2
        vtmp1(1) = -srcvals1(2,i)/rr1
        vtmp1(2) = srcvals1(1,i)/rr1
        vtmp1(3) = 0 
        call cross_prod3d(srcvals1(10,i),vtmp1,vtmp2)
        hvecs1(1:3,i,1) = cos(thet)*vtmp1(1:3) + 
     1     sin(thet)*vtmp2(1:3)
        hvecs1(1:3,i,2) = -sin(thet)*vtmp1(1:3) + 
     1     cos(thet)*vtmp2(1:3)
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
        vtmp1(1) = -srcvals2(2,i)/rr1
        vtmp1(2) = srcvals2(1,i)/rr1
        vtmp1(3) = 0 
        call cross_prod3d(srcvals2(10,i),vtmp1,vtmp2)
        hvecs2(1:3,i,1) = cos(thet)*vtmp1(1:3) + 
     1     sin(thet)*vtmp2(1:3)
        hvecs2(1:3,i,2) = -sin(thet)*vtmp1(1:3) + 
     1     cos(thet)*vtmp2(1:3)
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
c
c  Compute boundary data for analytic solution test
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


      thresh = 1.0d-14
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

      print *, "npts=",npts
      
      call l3ddirectcdg(1,xyz_out_src,cf2,vf2,1,sources,npts,
     1     ptmp,bbp,thresh)
      do i=1,npts1
        bbp(1:3,i) = 0
      enddo
c      do i=npts1+1,npts 
c        bbp(1,i) = srcvals(1,i)**2+srcvals(2,i)**2-2*srcvals(3,i)**2
c        bbp(2,i) = srcvals(1,i)
c        bbp(3,i) = srcvals(3,i)**3-3*srcvals(1,i)**2*srcvals(3,i)
c      enddo 
      call prin2('bbp=*',bbp(1,npts1+1),24)
      call surf_vtk_plot_vec(npatches1,norders1,ixyzs1,iptype1,npts1,
     1  srccoefs1,srcvals1,bbp,'bbp1.vtk','a')
      call surf_vtk_plot_vec(npatches2,norders2,ixyzs2,iptype2,npts2,
     1  srccoefs2,srcvals2,bbp(1,npts1+1),'bbp2.vtk','a')
      allocate(errp(npatches))
      
      init = 0

      call fieldsED(zk,xyz_out_src,srcvals,npts,zjm,zbbm,zf2,init)
      bjm = real(zjm)
      bbm = real(zbbm)
      errm = 0
      call surf_fun_error(3,npatches,norders,ixyzs,iptype,npts,bbm,
     1  wts,errp,errm)
      call prin2('bjm=*',bjm,24)
      call prin2('bbm=*',bbm,24)
      call prin2('ptmp=*',ptmp,24)

      allocate(rhs(6*npts+8*ngenus),soln(6*npts+8*ngenus))
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
c      bbp(1,:) = srcvals(2,:)
c      bbp(2,:) = 0
c      bbp(3,:) = srcvals(1,:)

      print *, maxval(bbp(1,:))
      print *, maxval(bbp(2,:))
      print *, maxval(bbp(3,:))

c
c
c  Our ordering: interior a_+, interior a_-, interior b_+, interior b_-
c                exterior a_+, exterior a_-, exterior b_+, exterior b_-
c   

      allocate(bbp1_a(3,na1),bbp2_a(3,na2))
      allocate(bbp1_b(3,nb1),bbp2_b(3,nb2))
      allocate(bbm1_a(3,na1),bbm2_a(3,na2))
      allocate(bbm1_b(3,nb1),bbm2_b(3,nb2))


C
CCC     for surface 1 
C
      print *, "npatches1=",npatches1
      print *, "npatches2=",npatches2
      call fun_surf_interp(3,npatches1,norders1,ixyzs1,iptype1,
     1     npts1,bbp,na1,apatches1,auv1,bbp1_a)
      call fun_surf_interp(3,npatches1,norders1,ixyzs1,iptype1,
     1     npts1,bbp,nb1,bpatches1,buv1,bbp1_b)

      call fun_surf_interp(3,npatches1,norders1,ixyzs1,iptype1,
     1     npts1,bbm,na1,apatches1,auv1,bbm1_a)
      call fun_surf_interp(3,npatches1,norders1,ixyzs1,iptype1,
     1     npts1,bbm,nb1,bpatches1,buv1,bbm1_b)

C
CCC     for surface 2
C

      call fun_surf_interp(3,npatches2,norders2,ixyzs2,iptype2,
     1     npts2,bbp(1,npts1+1),na2,apatches2,auv2,bbp2_a)
      call fun_surf_interp(3,npatches2,norders2,ixyzs2,iptype2,
     1     npts2,bbp(1,npts1+1),nb2,bpatches2,buv2,bbp2_b)

      call fun_surf_interp(3,npatches2,norders2,ixyzs2,iptype2,
     1     npts2,bbm(1,npts1+1),na2,apatches2,auv2,bbm2_a)
      call fun_surf_interp(3,npatches2,norders2,ixyzs2,iptype2,
     1     npts2,bbm(1,npts1+1),nb2,bpatches2,buv2,bbm2_b)

      print *,"na=",na
      print *, "na1=",na1
      print *, "na2=",na2
      print *, "nb=",nb
      print *, "nb1=",nb1
      print *, "nb2=",nb2
C
      rra = 0 
      do i=1,na1
        rhs(6*npts+1) = rhs(6*npts+1) + (bbm1_a(1,i)*avals1(4,i) + 
     1  bbm1_a(2,i)*avals1(5,i) + bbm1_a(3,i)*avals1(6,i))*awts1(i)
        rhs(6*npts+5) = rhs(6*npts+5) + (bbp1_a(1,i)*avals1(4,i) + 
     1  bbp1_a(2,i)*avals1(5,i) + bbp1_a(3,i)*avals1(6,i))*awts1(i)
        rra = rra + 
     1  sqrt(avals1(4,i)**2+avals1(5,i)**2+avals1(6,i)**2)*awts1(i)
      enddo
      print *, "rra1 = ",rra
      

      rra = 0
      do i=1,na2
        rhs(6*npts+2) = rhs(6*npts+2) + (bbm2_a(1,i)*avals2(4,i) + 
     1  bbm2_a(2,i)*avals2(5,i) + bbm2_a(3,i)*avals2(6,i))*awts2(i)
        rhs(6*npts+6) = rhs(6*npts+6) + (bbp2_a(1,i)*avals2(4,i) + 
     1  bbp2_a(2,i)*avals2(5,i) + bbp2_a(3,i)*avals2(6,i))*awts2(i)
        rra = rra + 
     1  sqrt(avals2(4,i)**2 + avals2(5,i)**2 + avals2(6,i)**2)*awts2(i)
      enddo
      print *, "rra2=",rra

      rrb = 0
      do i=1,nb1
        rhs(6*npts+3) = rhs(6*npts+3) + (bbm1_b(1,i)*bvals1(4,i) + 
     1  bbm1_b(2,i)*bvals1(5,i) + bbm1_b(3,i)*bvals1(6,i))*bwts1(i)
        rhs(6*npts+7) = rhs(6*npts+7) + (bbp1_b(1,i)*bvals1(4,i) + 
     1  bbp1_b(2,i)*bvals1(5,i) + bbp1_b(3,i)*bvals1(6,i))*bwts1(i)
        rrb = rrb + 
     1  sqrt(bvals1(4,i)**2 + bvals1(5,i)**2 + bvals1(6,i)**2)*bwts1(i)
      enddo
      print *, "rrb1=",rrb

      rrb = 0 
      do i=1,nb2
        rhs(6*npts+4) = rhs(6*npts+4)+(bbm2_b(1,i)*bvals2(4,i)+ 
     1  bbm2_b(2,i)*bvals2(5,i) + bbm2_b(3,i)*bvals2(6,i))*bwts2(i)
        rhs(6*npts+8) = rhs(6*npts+8) + (bbp2_b(1,i)*bvals2(4,i) + 
     1  bbp2_b(2,i)*bvals2(5,i) + bbp2_b(3,i)*bvals2(6,i))*bwts2(i)
        rrb = rrb + 
     1  sqrt(bvals2(4,i)**2 + bvals2(5,i)**2 + bvals2(6,i)**2)*bwts2(i)
      enddo
      print *, "rrb2=",rrb

      call prin2('rhs_projs=*',rhs(6*npts+1),8)
      print *, "estimated error in resolving bbm=",errm
      print *, "rra = ",rra
      print *, "rrb = ", rrb

      

      eps = 0.51d-6

      allocate(ipatch_id(npts),uvs_src(2,npts))
      call get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts,
     1  ipatch_id,uvs_src)
c      call prinf('ipatch_id=*',ipatch_id,20)
c      call prin2('uvs_src=*',uvs_src,48)

      allocate(cms(3,npatches),rads(npatches),rad_near(npatches))
      call get_centroid_rads(npatches,norders,ixyzs,iptype,npts, 
     1     srccoefs,cms,rads)
c      call prin2('cms=*',cms,24)
c      call prin2('rads=*',rads,6)
      
      
      
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

      ifquadread = 1
      if(ifquadread.eq.1) then
        open(unit=99,file='quad_iref2.bin',form='unformatted')
        read(99) wnear
        read(99) wnear1
        read(99) wnear2
        close(99)
      else


        call cpu_time(t1)
C$        t1 = omp_get_wtime()      
        call getnearquad_statj_gendeb(npatches,norders,
     1      ixyzs,iptype,npts,srccoefs,srcvals,
     1      epsquad,dzk,iquadtype,nnz,row_ptr,col_ind,
     1      iquad,rfac0,nquad,wnear)
 1111 continue     
        call prinf('finished generating near quadrature correction*',
     1     i,0)
        call cpu_time(t2)
C$      t2 = omp_get_wtime()      
        call prin2('time taken to generate near quadrature=*',t2-t1,1)


C
C     surface 1 
C
        call getnearquad_statj_gendeb(npatches1,norders1,
     1      ixyzs1,iptype1,npts1,srccoefs1,srcvals1,
     2      epsquad,dzk,iquadtype,nnz1,row_ptr1,col_ind1,
     3      iquad1,rfac01,nquad1,wnear1)


C
C     surface 2 
C
        call getnearquad_statj_gendeb(npatches2,norders2,
     1      ixyzs2,iptype2,npts2,srccoefs2,srcvals2,
     2      epsquad,dzk,iquadtype,nnz2,row_ptr2,col_ind2,
     3      iquad2,rfac02,nquad2,wnear2)
      endif
      
      ifquadwrite = 0
      if(ifquadwrite.eq.1) then
        open(unit=99,file='quad_iref2.bin',form='unformatted')
        write(99) wnear
        write(99) wnear1
        write(99) wnear2
        close(99)
      endif

      call prin2('wnear=*',wnear,24)
      call prin2('wnear1=*',wnear1,24)
      call prin2('wnear2=*',wnear2,24)

c
c      print *, maxval(wnear)
c      print *, maxval(wnear1)
c      print *, maxval(wnear2)
c      print *, "nquad=", nquad
c      print *, "nquad1=",nquad1
c      print *, "nquad2=",nquad2
c      call prinf('ixyzs2=*',ixyzs2,npatches2+1)


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

        call prinf('igen=*',igen,1)
        call prin2('bbphvecs1=*',bbphvecs1(1,1,igen),24)


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
        call prinf('igen=*',igen,1)
        call prin2('bbphvecs2=*',bbphvecs2(1,1,igen),24)

c        call prin2('outtmp2=*',outtmp2,24)
c        call prin2('rhstmp2=*',rhstmp2,24)
      enddo



      print *, "here"

      numit = 300
      allocate(errs(numit+1))
      call prinf('ngenus=*',ngenus,1)
      eps_gmres = 1.0d-6
      if(ngenus.ge.1)  call prin2('rhs_projs=*',rhs(6*npts+1),8)
      print *, "here"

      rhs = 0
      rhs(6*npts+8) = 0

      rsurint1 = 0
      do i=1,npts1
        rsurfint1 = rsurfint1 + wts1(i)
      enddo

      rsurint2 = 0
      do i=1,npts2
        rsurfint2 = rsurfint2 + wts2(i)
      enddo

      do i=1,npts2
        rhs(5*npts+npts1+i) = 1.0d0/rsurfint2
      enddo
      

      soln = 0
      call lpcomp_statj_gendeb_thinshell_addsub(npatches,npatches1,
     1   norders,ixyzs,iptype,npts,srccoefs,srcvals,eps,dpars,nnz,
     2   row_ptr,col_ind,iquad,nquad,wnear,nnz1,npts1,row_ptr1,
     3   col_ind1,iquad1,nquad1,wnear1,nnz2,npts2,row_ptr2,
     4   col_ind2,iquad2,nquad2,wnear2,hvecs1,bbphvecs1,hvecs2,
     5   bbphvecs2,na,na1,iaxyzs,apatches,auv,avals,awts,nb,
     6   nb1,ibxyzs,bpatches,buv,bvals,bwts,rhs,novers,npts_over,
     7   ixyzso,srcover,wover,soln)
      soln = soln + rhs
      call prin2('soln1=*',soln,24)
      call prin2('soln2=*',soln(npts+1),24)
      call prin2('soln3=*',soln(2*npts+1),24)
      call prin2('soln4=*',soln(3*npts+1),24)
      call prin2('soln5=*',soln(4*npts+1),24)
      call prin2('soln6=*',soln(5*npts+1),24)
      call prin2('soln proj=*',soln(6*npts+1),8)

      ra = 0
      do i=1,npts
        do j=0,5
          ra = ra + soln(j*npts+i)**2*wts(i)
        enddo
      enddo

      do i=1,8
        ra = ra + soln(6*npts+i)**2
      enddo
      ra = sqrt(ra)
      call prin2('norm of vector=*',ra,1)

      stop


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
      
      call prin2('solve time=*',t2-t1,1)      
      
      call prin2('projs=*',soln(6*npts+1),8)


      if(1.eq.1) then
        allocate(bbpcomp(3,npts),bjmcomp(3,npts),bbmcomp(3,npts))

        print *, "here"
        call prin2('eps=*',eps,1)
        call prin2('soln=*',soln,24)
        call prinf('ngenus=*',ngenus,1)
        call prin2('dpars=*',dpars,3)

        call prin2('sol comp1 surface 1=*', soln,24)
        call prin2('sol comp1 surface 2=*', soln(npts1+1),24)
        call prin2('sol comp2 surface 1=*', soln(npts+1), 24)
        call prin2('sol comp2 surface 2=*', soln(npts+npts1+1), 24)
        call prin2('sol comp3 surface 1=*', soln(2*npts+1), 24)
        call prin2('sol comp3 surface 2=*', soln(2*npts+npts1+1), 24)
        call prin2('sol comp4 surface 1=*', soln(3*npts+1), 24)
        call prin2('sol comp4 surface 2=*', soln(3*npts+npts1+1), 24)
        call prin2('sol comp5 surface 1=*', soln(4*npts+1), 24)
        call prin2('sol comp5 surface 2=*', soln(4*npts+npts1+1), 24)
        call prin2('sol comp6 surface 1=*', soln(5*npts+1), 24)
        call prin2('sol comp6 surface 2=*', soln(5*npts+npts1+1), 24)

c
c
c       new homework 
c

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
        call prin2('absolute errj=*', sqrt(errj), 1)
        call prin2('absolute errbm=*', sqrt(errbm), 1)
        call prin2('absolute errbp=*', sqrt(errbp), 1)


        errj = sqrt(errj)/rsurfintl2
        errbm = sqrt(errbm)/rsurfintl2
        errbp = sqrt(errbp)/rsurfintl2
        call prin2('error in current=*',errj,1)
        call prin2('error in interior magnetic field=*',errbm,1)
        call prin2('error in exterior magnetic field=*',errbp,1)
        call prin2('rsurfintl2 =*',rsurfintl2,1)
      endif
      call surf_vtk_plot_scalar(npatches1,norders1,ixyzs1,iptype1,npts1,
     1  srccoefs1,srcvals1,soln(3*npts+1),'qm1_iref3_1s.vtk','a')
      call surf_vtk_plot_scalar(npatches2,norders2,ixyzs2,iptype2,npts2,
     1  srccoefs2,srcvals2,soln(3*npts+npts1+1),'qm2_iref3_1s.vtk','a')

      call surf_vtk_plot_scalar(npatches1,norders1,ixyzs1,iptype1,npts1,
     1  srccoefs1,srcvals1,soln(4*npts+1),'qp1_iref3_1s.vtk','a')
      call surf_vtk_plot_scalar(npatches2,norders2,ixyzs2,iptype2,npts2,
     1  srccoefs2,srcvals2,soln(4*npts+npts1+1),'qp2_iref3_1s.vtk','a')


      call surf_vtk_plot_scalar(npatches1,norders1,ixyzs1,iptype1,npts1,
     1  srccoefs1,srcvals1,soln(5*npts+1),'rm1_iref3_1s.vtk','a')
      call surf_vtk_plot_scalar(npatches2,norders2,ixyzs2,iptype2,npts2,
     1  srccoefs2,srcvals2,soln(5*npts+npts1+1),'rm2_iref3_1s.vtk','a')


      rqm1 = 0
      rqp1 = 0
      rrm1 = 0

      rqm2 = 0
      rqp2 = 0
      rrm2 = 0
      do i=1,npts1
        rqm1 = rqm1 + soln(3*npts+i)**2*wts1(i)
        rqp1 = rqp1 + soln(4*npts+i)**2*wts1(i)
        rrm1 = rrm1 + soln(5*npts+i)**2*wts1(i)
      enddo

      do i=1,npts2
        rqm2 = rqm2 + soln(3*npts+npts1+i)**2*wts2(i)
        rqp2 = rqp2 + soln(4*npts+npts1+i)**2*wts2(i)
        rrm2 = rrm2 + soln(5*npts+npts1+i)**2*wts2(i)
      enddo

      rqm1 = sqrt(rqm1)
      rqp1 = sqrt(rqp1)
      rrm1 = sqrt(rrm1)
      
      rqm2 = sqrt(rqm2)
      rqp2 = sqrt(rqp2)
      rrm2 = sqrt(rrm2)

      call prin2('rqm1 = *', rqm1, 1)
      call prin2('rqp1 = *', rqp1, 1)
      call prin2('rrm1 = *', rrm1, 1)

      call prin2('rqm2 = *', rqm2, 1)
      call prin2('rqp2 = *', rqp2, 1)
      call prin2('rrm2 = *', rrm2, 1)

      

      

c     
c     can stop here, expect to get digits 
c 


      rsurfintl2 = 0 
      do i=1,npts 
        do j=0,5 
          rsurfintl2 = rsurfintl2 + soln(j*npts + i)**2*wts(i)
        enddo 
      enddo 
      rsurfintl2 = sqrt(rsurfintl2)


c
c  test solution at interior and exterior points  
c
      rinttmp = 0
      rsurf = 0
      do i=1,npts1
        rsurf = rsurf + wts(i)
        do j=1,6
          rinttmp(j) = rinttmp(j) + soln(i+(j-1)*npts)*wts(i)
        enddo
      enddo
      rinttmp(1:6) = rinttmp(1:6)/rsurf
      call prin2('average integral of densities on surf1=*',rinttmp,6)
      
      rinttmp = 0
      rsurf = 0
      do i=npts1+1,npts
        rsurf = rsurf + wts(i)
        do j=1,6
          rinttmp(j) = rinttmp(j) + soln(i+(j-1)*npts)*wts(i)
        enddo
      enddo
      rinttmp(1:6) = rinttmp(1:6)/rsurf
      call prin2('average integral of densities on surf2=*',rinttmp,6)
      
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
c
c  Yuguna: Mimic whatever is done in the addsub layer potential evaluator
c  routine to make sure blm = \ell^{-} is correct
c
c




      allocate(soln1(6*npts1))
      do i=1,npts1
        do j=0,5 
          soln1(j*npts1+i) = soln(j*npts+i)
        enddo
      enddo




      
      call statj_gendebproc_rhomrhopmum(npatches1,norders1,ixyzs1,
     1 iptype,npts1,srccoefs1,srcvals1,eps,nnz1,row_ptr1,col_ind1, 
     2 iquad1,nquad1,wnear1,soln1,novers1,npts_over1,ixyzso1,srcover1,
     3 wover1,curv,wtmp1,wtmp2,wtmp3,wtmp4,dzk,rbeta,rgamma,laps02rhom,
     4 laps02rhop,laps02mum,blm,bmm)


      allocate(soln2(6*npts2))
      do i=1,npts2
        do j=0,5 
          soln2(j*npts2+i) = soln(j*npts+npts1+i)
        enddo
      enddo


      call statj_gendebproc_rhomrhopmum(npatches2,norders2,ixyzs2,
     1 iptype,npts2,srccoefs2,srcvals2,eps,nnz2,row_ptr2,col_ind2,
     2 iquad2,nquad2,wnear2,soln2,novers2,npts_over2,ixyzso2,srcover2,
     3 wover2,curv(npts1+1), wtmp1(1,npts1+1),wtmp2(1,npts1+1),
     4 wtmp3(1,npts1+1),wtmp4(1,npts1+1),dzk,rbeta,rgamma,
     5 laps02rhom(npts1+1),laps02rhop(npts1+1),laps02mum(npts1+1),
     6 blm(1,npts1+1),bmm(1,npts1+1))
     

c
c  add in contribution of harmonic vector fields
c
      do i=1,npts
        do igen=1,2*ngenus
          if (i.le.npts1) then 
            blm(1:3,i) = blm(1:3,i) + 
     1       soln(6*npts+igen+4)*hvecs1(1:3,i,igen)
          else 
            blm(1:3,i) = blm(1:3,i) + 
     1       soln(6*npts+igen+2+4)*hvecs2(1:3,i-npts1,igen)
          endif 
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
        call l3ddirectcdg(1,xyz_in_src,cf2,vf2,1,xyz_out_src,
     1     1,ptmp,bbpex,thresh)
        bbpex = 0

        do i=1,npts
          dx = xyz_out_src(1) - srcvals(1,i)
          dy = xyz_out_src(2) - srcvals(2,i)
          dz = xyz_out_src(3) - srcvals(3,i)
        
          sig = soln(4*npts+i)
        
          r = sqrt(dx**2 + dy**2 + dz**2)

          bbpc(1) = bbpc(1) - dx/r**3*sig*wts(i)
          bbpc(2) = bbpc(2) - dy/r**3*sig*wts(i)
          bbpc(3) = bbpc(3) - dz/r**3*sig*wts(i)
c
c  Note: For targets outside the outer torus, we add contributions
c  of c_{+} and for targets inside of inner torus, we need
c  to add contributions only from the c_{-} part
c 
c
          if (i.le.npts1) then 
            do igen=1,2*ngenus
        
              bbpc(1) = bbpc(1)+1.0d0/r**3*wts(i)*soln(6*npts+igen)
     1      *(dy*hvecs1(3,i,igen)-dz*hvecs1(2,i,igen))
              bbpc(2) = bbpc(2)+1.0d0/r**3*wts(i)*soln(6*npts+igen)
     1      *(dz*hvecs1(1,i,igen)-dx*hvecs1(3,i,igen))
              bbpc(3) = bbpc(3)+1.0d0/r**3*wts(i)*soln(6*npts+igen)
     1      *(dx*hvecs1(2,i,igen)-dy*hvecs1(1,i,igen))
            enddo

          else 

c            do igen=1,2*ngenus
c        
c              bbpc(1) = bbpc(1) + 1.0d0/r**3*wts(i)*
c     1                    soln(6*npts+igen+2)*
c     2      (dy*hvecs2(3,i-npts1,igen)-dz*hvecs2(2,i-npts1,igen))
c              bbpc(2) = bbpc(2) + 1.0d0/r**3*wts(i)*
c     1                      soln(6*npts+igen+2)*
c     2      (dz*hvecs2(1,i-npts1,igen)-dx*hvecs2(3,i-npts1,igen))
c              bbpc(3) = bbpc(3) + 1.0d0/r**3*wts(i)*
c     1                     soln(6*npts+igen+2)*
c     2      (dx*hvecs2(2,i-npts1,igen)-dy*hvecs2(1,i-npts1,igen))
c          enddo

          endif 

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
      call prin2('xyz_in_src=*',xyz_in_src,3)
      do j=1,ntargtest
        bbmc(1:3) = 0
        bbmex(1:3) = 0
        bjmc(1:3) = 0
        bjmex(1:3) = 0
        ptmp = 0
        dx = xyz_in_src(1)-xyz_out_src(1)
        dy = xyz_in_src(2)-xyz_out_src(2) 
        dz = xyz_in_src(3)-xyz_out_src(3) 
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
          dx = xyz_in_src(1) - srcvals(1,i)
          dy = xyz_in_src(2) - srcvals(2,i)
          dz = xyz_in_src(3) - srcvals(3,i)
        
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


      call prin2('bjmex=*',bjmex,3)
      call prin2('bjm=*',bjmc,3)
      call prin2('bbmex=*',bbmex,3)
      call prin2('bbm=*',bbmc,3)
      call prin2('bbpex=*',bbpex,3)
      call prin2('bbp=*',bbpc,3)

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


      stop 


     
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






