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
      real *8, allocatable :: rhs_test(:), soln_test(:)
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
      character *100 title,dirname, dirname2
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
      real *8 v1(3),v2(3),p(3),rp(3)
      real *8, allocatable :: xcircle(:,:),dxcircle(:,:)

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
      idzk = 10
      dzk = 10**(idzk)
      if(idzk.eq.3) dzk = 30.0d0
      if(idzk.ge.4) dzk = 2**(idzk-4)
      zk = ima*dzk

      ifb = 0

      if(ifb.eq.0) ifa = 1
      if(ifb.eq.1) ifa = 0

c
c  Note for this code, igeomtype can only be 6 
c
c  igeomtype = 6 => thin shell tori
c
      igeomtype = 6
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


            
      allocate(norders2(npatches2),ixyzs2(npatches2+1),
     1        iptype2(npatches2))
        
      do i=1,npatches2
        norders2(i) = norder
        ixyzs2(i) = 1 +(i-1)*npols
        iptype2(i) = 1
      enddo 

      ixyzs2(npatches2+1) = 1+npols*npatches2

      allocate(wts(npts))
      call get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,wts)


      write(fname,'(a,a,i2.2,a,i2.2,a,i1,a,i1,a)') 
     1    trim(dirname),'statj_soln_lambdainf_',ipars(1),'_',ipars(2),
     2    '_norder',norder,'_ifb',ifb,'.dat'

      allocate(bbpcomp(3,npts1),soln(npts1+2))

      print *, "fname=", trim(fname)

      open(unit=80,file=trim(fname),form='unformatted')
 
      read(80) ipars
      read(80) iref
      read(80) ifb
      read(80) tsolve
      read(80) niter
      read(80) rres
      read(80) soln
      read(80) bbpcomp

      close(80)
      iref = iref + 1

      write(dirname2, '(a,i1,a,i1,a,i1,a)') 
     1   'static_curr_results/ap_',
     1   ifa,'_bm_',ifb,'_iref_',iref,'_idzk_0'

      if(ifb.eq.0) then
      
        fname = trim(dirname2)//'/bbp_outer.vtk'
        call surf_vtk_plot_vec(npatches1, norders1, ixyzs1, iptype1, 
     1    npts1, srccoefs1, srcvals1, bbpcomp, trim(fname), 'a')
      else

        fname = trim(dirname2)//'/bbp_inner.vtk'
        call surf_vtk_plot_vec(npatches2, norders2, ixyzs2, iptype2, 
     1    npts2, srccoefs2, srcvals2, bbpcomp, trim(fname),'a')
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






