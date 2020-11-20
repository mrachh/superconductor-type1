      
      implicit real *8 (a-h,o-z)
      real *8, allocatable :: amat(:,:,:),bmat(:,:,:)
      real *8, allocatable :: amatc(:,:,:),bmatc(:,:,:)
      real *8, allocatable :: aeigr(:,:),beigr(:,:)
      real *8, allocatable :: aeigi(:,:),beigi(:,:)
      real *8, allocatable :: aceigr(:,:),bceigr(:,:)
      real *8, allocatable :: aceigi(:,:),bceigi(:,:)
      real *8, allocatable :: work(:),wtmp(:)

      lw = 100000
      allocate(wtmp(10000),work(lw))

      call prini(6,13)


      nmax = 100
      allocate(amat(6,6,0:nmax),bmat(3,3,0:nmax))
      allocate(amatc(6,6,0:nmax),bmatc(3,3,0:nmax))
      allocate(aeigr(6,0:nmax),beigr(3,0:nmax))
      allocate(aeigi(6,0:nmax),beigi(3,0:nmax))
      allocate(aceigr(6,0:nmax),bceigr(3,0:nmax))
      allocate(aceigi(6,0:nmax),bceigi(3,0:nmax))

      rbeta = 1.0d0
      rgamma = 0.0d0

      do i=0,nmax
        call get_ynm_mat(i,rbeta,rgamma,amat(1,1,i),amatc(1,1,i),
     1     bmat(1,1,i),bmatc(1,1,i))
        print *, bmat(2,2,0)
        call prin2('amat=*',amat(1,1,i),36)
        call prin2('bmat=*',bmat(1,1,i),9)
        write(77,'(6(2x,e11.5))') amatc(4:6,2,i),amatc(4:6,5,i)
        write(78,'(3(2x,e11.5))') bmatc(1:3,2,i)
        call dgeev('N','N',6,amat(1,1,i),6,aeigr(1,i),aeigi(1,i),
     1     wtmp,1,wtmp,1,work,lw,info)
        call dgeev('N','N',3,bmat(1,1,i),3,beigr(1,i),beigi(1,i),
     1     wtmp,1,wtmp,1,work,lw,info)
        
        call dgeev('N','N',6,amatc(1,1,i),6,aceigr(1,i),aceigi(1,i),
     1     wtmp,1,wtmp,1,work,lw,info)
        call dgeev('N','N',3,bmatc(1,1,i),3,bceigr(1,i),bceigi(1,i),
     1     wtmp,1,wtmp,1,work,lw,info)
        

        write(34,'(6(2x,e11.5))') aeigr(1:6,i)
        write(35,'(6(2x,e11.5))') aceigr(1:6,i)
        write(36,'(3(2x,e11.5))') beigr(1:3,i)
        write(37,'(3(2x,e11.5))') bceigr(1:3,i)
      enddo


      return
      end


      subroutine get_ynm_mat(nn,rbeta,rgamma,amat,amatc,bmat,bmatc)

      implicit real *8 (a-h,o-z) 
      real *8 amat(6,6),bmat(3,3),amatc(6,6),bmatc(3,3)

      complex *16 zalpha,zbeta,zgamma,zdelta,zeta,zteta,zk,ztetap
      complex *16 ztetam
      real *8 dalpha,dbeta,dgamma,ddelta,deta,dteta
      complex *16 fjvals(0:1000),fhvals(0:1000),
     1  fjder(0:1000),fhder(0:1000) 

      complex * 16 zpars(3)
      integer numit,niter
      character *100 title,dirname
      character *300 fname

      real *8, allocatable :: w(:,:)

      logical isout0,isout1

      complex *16 ztmp,ima

      data ima/(0.0d0,1.0d0)/


      call prini(6,13)

      done = 1
      pi = atan(done)*4



c
c       define rhs to be one of the ynm's
c
      mm = 1
      nmax = nn
c
c  set the value of dzk
c
      dzk = 1.0d0


      njh = nn+5
      ifder = 1
      rscale = 1.0d0
      zk = ima*dzk
      call prin2('zk=*',zk,2)
      call besseljs3d(njh,zk,rscale,fjvals,ifder,fjder)
      call h3dall(njh,zk,rscale,fhvals,ifder,fhder)


      zalpha = ima*((nn+1)*fjvals(nn)*fhvals(nn-1)+
     1   nn*fjvals(nn+1)*fhvals(nn)-zk*fjvals(nn+1)*fhvals(nn-1))
      dalpha = real(zalpha)
cc      call prin2('zalpha=*',zalpha,2)
      
      zbeta = fjvals(nn)*fhvals(nn)*zk*zk/ima/zk
      dbeta = real(zbeta)
cc      call prin2('zbeta=*',zbeta,2)
      
      zgamma = ima*zk*fjvals(nn)*fhvals(nn)*sqrt(nn*(nn+1.0d0))
      dgamma = real(zgamma)
cc      call prin2('zgamma=*',zgamma,2)
      
      zdelta = ima*((fjvals(nn)+zk*fjder(nn))*fhvals(nn)*zk + 
     1  (fhvals(nn)+zk*fhder(nn))*fjvals(nn)*zk)/2.0d0
      ddelta = real(zdelta)
cc      call prin2('zdelta=*',zdelta,2)
      
      zeta = ima*sqrt(nn*(nn+1.0d0))*(fjvals(nn)*fhvals(nn-1)-
     1   fjvals(nn+1)*fhvals(nn))
      deta = real(zeta)
      
cc      call prin2('zeta=*',zeta,2)
      
      zteta = ima*zk**2/2*(fjvals(nn)*fhder(nn)+fjder(nn)*fhvals(nn))
      dteta = real(zteta)
cc      call prin2('zteta=*',zteta,2)


      ztetap = ima*zk**2*fjvals(nn)*fhder(nn)
      dtetap = real(ztetap)
cc      call prin2('ztetap=*',ztetap,2)

      ztetam = ima*zk**2*fjder(nn)*fhvals(nn)
      dtetam = real(ztetam)
cc      call prin2('ztetam=*',ztetam,2)


      rslaps = -(nn+0.0d0)*(nn+1.0d0)/(2*nn+1.0d0)**2
      if(nn.eq.0) rslaps = rslaps + 1.0d0
      rs = 1
      
      rinv1 = -1.0d0/(nn+0.0d0)/(nn+1.0d0)

      if(nn.eq.0) rinv1 = 0.0d0
      amat =  0
      bmat = 0

      amat(1,1) = rslaps
      amat(1,4) = -rs
      amat(2,2) = rslaps
      amat(2,5) = -rs
      amat(3,3) = rslaps
      amat(3,6) = -rs


      amat(4,1) = -dzk*(nn+0.0d0)*(nn+1.0d0)*dalpha/(2*nn+1.0d0)**3
      amat(4,2) = 
     1   -rbeta*(nn+0.0d0)*(nn+1.0d0)*ddelta/(2*nn+1.0d0)**2/dzk
      amat(4,3) = 0
      amat(4,4) = dgamma*sqrt(nn*(nn+1.0d0))/(2*nn+1.0d0)/dzk
      amat(4,5) = -nn*(nn+1.0d0)/(2*nn+1.0d0)**2
      amat(4,6) = 0

      amat(5,1) = -dzk*sqrt(nn*(nn+1.0d0))*deta/(2*nn+1.0d0)**2
      amat(5,2) = rbeta*dbeta*nn*(nn+1.0d0)/(2*nn+1.0d0)**2/dzk 
      amat(5,3) = 0
      amat(5,4) = dtetam/dzk
      amat(5,5) = (nn+1.0d0)/(2*nn+1.0d0)
      amat(5,6) = 0

      amat(6,1) = 0
      amat(6,2) = -dbeta*rgamma*nn*(nn+1.0d0)/(2*nn+1.0d0)**2
      amat(6,3) = dzk**2*sqrt(nn*(nn+1.0d0))/(2*nn+1.0d0)*deta
      amat(6,4) = 0
      amat(6,5) = 0
      amat(6,6) = -dtetam

      bmat(1:3,1:3) = amat(4:6,4:6) + amat(4:6,1:3)*rinv1
      print *, amat(5,5),bmat(2,2)


      amatc = 0
      bmatc = 0
      amatc(1:3,1:6) = -4*amat(1:3,1:6)
      amatc(4,1:6) = (amat(5,1:6) + 2*amat(4,1:6))*dzk
      amatc(5,1:6) = (amat(5,1:6) - 2*amat(4,1:6))
      amatc(6,1:6) = -2*amat(6,1:6)

      bmatc(1,1:3) = (bmat(2,1:3) + 2*bmat(1,1:3))*dzk
      bmatc(2,1:3) = (bmat(2,1:3) - 2*bmat(1,1:3))
      bmatc(3,1:3) = -2*bmat(3,1:3)

      if(nn.eq.0) then
        amatc = 0
        bmatc = 0
        amatc(1,1) = -4
        amatc(2,2) = -4
        amatc(3,3) = -4
        amatc(2,5) = 4
        amatc(4,4) = 0.3d0
        amatc(4,5) = dzk
        amatc(5,4) = 1.3d0
        amatc(5,5) = 1.0d0
        amatc(6,6) = 0.3d0

        bmatc(1,1) = 0.3d0
        bmatc(1,2) = dzk
        bmatc(2,1) = 1.3d0
        bmatc(2,2) = 1.0d0
        bmatc(3,3) = 0.3d0
      endif

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


      external xtri_stell_eval,xtri_sphere_eval
      
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
        vmin = 0
        vmax = 2*pi
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
      real *8, allocatable :: rsurf(:),err_p(:,:) 
      integer ipars,norderhead,nd
      complex *16, allocatable :: sigma_coefs(:,:), sigma_vals(:,:)
      complex *16 zk,val

      integer ipatch,j,i
      real *8 ra,ds
      logical isout

      done = 1
      pi = atan(done)*4

      npols = (norder+1)*(norder+2)/2


      zk = 0

      ra = 0



      do ipatch=1,npatches
        do j=1,npols
          i = (ipatch-1)*npols + j
          call h3d_sprime(xyzout,srcvals(1,i),dpars,zk,ipars,val)
          call cross_prod3d(srcvals(4,i),srcvals(7,i),tmp)
          ds = sqrt(tmp(1)**2 + tmp(2)**2 + tmp(3)**2)
          ra = ra + real(val)*wts(i)
        enddo
      enddo

      if(abs(ra+4*pi).le.1.0d-3) isout = .false.
      if(abs(ra).le.1.0d-3) isout = .true.

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



