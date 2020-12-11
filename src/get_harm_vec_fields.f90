      subroutine get_harm_vec_field(npatches,norders,ixyzs,iptype, &
        npts,srccoefs,srcvals,wts,eps,H,errest)
!
!  This subroutine computes the harmonic vector fields for a genus
!  ngenus object described in the go3 format.
!
!  The harmonic fields returned are orthogonalized using one pass
!  of gram-schmidt
!
!  An estimated error which computes the max of surface divergence 
!  and the n \times surface divergence. 
!
!
!  Input arguments:
! 
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization on each patch 
!    - ixyzs: integer(npatches+1)
!        ixyzs(i) denotes the starting location in srccoefs,
!        and srcvals array corresponding to patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangular patch discretized using RV nodes
!    - npts: integer
!        total number of discretization points on the boundary
!    - srccoefs: real *8 (9,npts)
!        koornwinder expansion coefficients of xyz, dxyz/du,
!        and dxyz/dv on each patch. 
!        For each point 
!          * srccoefs(1:3,i) is xyz info
!          * srccoefs(4:6,i) is dxyz/du info
!          * srccoefs(7:9,i) is dxyz/dv info
!    - srcvals: real *8 (12,npts)
!        xyz(u,v) and derivative info sampled at the 
!        discretization nodes on the surface
!          * srcvals(1:3,i) - xyz info
!          * srcvals(4:6,i) - dxyz/du info
!          * srcvals(7:9,i) - dxyz/dv info
!          * srcvals(10:12,i) - normals info
!    - eps: real *8
!        precision requested
!    - ngenus: integer
!        genus of the geometry
!
!  Output arguments:
!    hvecs: real *8 (3,npts,2*ngenus)
!      the orthonormalized harmonic vector fields, the first ngenus 
!      components orthonormal to each other and the remaining
!      ngenus components are defined to be the n \times the first
!      ngenus components. 
!    
!
      implicit real *8 (a-h,o-z)
      
      real *8 srcvals(12,npts),srccoefs(9,npts),wts(npts)
      integer norders(npatches),ixyzs(npatches+1),iptype(npatches)
      real *8 H(3,npts),errest
      real *8, allocatable :: rsigma(:),nF(:,:)
      integer ipars(2),d(3)
      integer, allocatable :: ipatch_id(:)
      real *8, allocatable :: uvs_targ(:,:),V(:,:),beta(:),alpha(:)
      real *8 dpars(2)

      real *8, allocatable :: ffform(:,:,:),ffformex(:,:,:)
      real *8, allocatable :: ffforminv(:,:,:),ffformexinv(:,:,:)

      real *8, allocatable :: sigma(:), pot(:),pot1(:),rrhs(:)
      real *8, allocatable :: errs(:),rrhs2(:),u(:),rrhs2t(:)

      real *8, allocatable :: wnear(:),rrhs1(:),nFt(:,:)
      real *8, allocatable :: targs(:,:),F(:,:),Ft(:,:)
      real *8, allocatable :: cms(:,:),rads(:),rad_near(:)
      real *8, allocatable :: srcover(:,:),wover(:)
      real *8, allocatable :: xmat(:,:),xtmp(:,:),slmat(:,:)
      real *8, allocatable :: slapmat(:,:), s0mat(:,:),nsgbeta(:,:)

      real *8, allocatable :: rfds(:),sgalpha(:,:),sgbeta(:,:),gu(:,:)
      integer, allocatable :: ifds(:)
      complex *16, allocatable :: zfds(:)


      integer, allocatable :: col_ptr(:),row_ind(:)
      integer, allocatable :: ixyzso(:),novers(:)
      integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)
      integer, allocatable :: seed(:),seed_old(:)


      real *8 thet,phi,eps_gmres,Wg,Wu,rnd,Hnorm,Fnorm
      real *8 rpars(4,3),xyz0(3),dvec(3),vtmp1(3),vtmp2(3)
      complex * 16 zpars(3)
      integer numit,niter,sizer,n_min,n_max,ran_len,ran_int,clock
      character *100 title,dirname
      character *300 fname,fname1,fname2,fname3,fname4,fname5

      real *8, allocatable :: w(:,:),rsc(:)

      logical isout0,isout1

      complex *16 ztmp,ima

      data ima/(0.0d0,1.0d0)/


      done = 1
      pi = atan(done)*4

      
      
      allocate(rrhs2t(npts))
      allocate(sigma(npts),pot(npts),rrhs(npts))
      allocate(ffform(2,2,npts),rrhs2(npts),u(npts))
      allocate(rrhs1(npts),Ft(3,npts),nFt(3,npts))
      allocate(V(3,npts),F(3,npts),nF(3,npts))
      allocate(alpha(npts),beta(npts),nsgbeta(3,npts))
      allocate(sgalpha(3,npts),sgbeta(3,npts))
      allocate(gu(3,npts))

!
!       precompute near quadrature correction
!
!
      iptype_avg = floor(sum(iptype)/(npatches+0.0d0))
      norder_avg = floor(sum(norders)/(npatches+0.0d0))

      call get_rfacs(norder_avg,iptype_avg,rfac,rfac0)


      allocate(cms(3,npatches),rads(npatches),rad_near(npatches))

      call get_centroid_rads(npatches,norders,ixyzs,iptype,npts, &
          srccoefs,cms,rads)

!$OMP PARALLEL DO DEFAULT(SHARED) 
      do i=1,npatches
        rad_near(i) = rads(i)*rfac
      enddo
!$OMP END PARALLEL DO     

      ntarg = npts
      allocate(targs(3,npts))
!$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,npts 
        targs(1,i) = srcvals(1,i)
        targs(2,i) = srcvals(2,i)
        targs(3,i) = srcvals(3,i)
      enddo
!$OMP END PARALLEL DO      

!
!    find near quadrature correction interactions
!
      call findnearmem(cms,npatches,rad_near,3,targs,npts,nnz)

      allocate(row_ptr(npts+1),col_ind(nnz))
      
      call findnear(cms,npatches,rad_near,3,targs,npts,row_ptr, &
             col_ind)

      allocate(iquad(nnz+1)) 
      call get_iquad_rsc(npatches,ixyzs,ntarg,nnz,row_ptr,col_ind, &
           iquad)

      allocate(ipatch_id(npts),uvs_targ(2,npts))
      call get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts, &
       ipatch_id,uvs_targ)
      
!
!    estimate oversampling for far-field, and oversample geometry
!

      ikerorder = 0
      allocate(novers(npatches),ixyzso(npatches+1))

      zpars = 0
      ndtarg = 3
      call get_far_order(eps,npatches,norders,ixyzs,iptype,cms, &
         rads,npts,srccoefs,ndtarg,npts,targs,ikerorder,zpars, &
         nnz,row_ptr,col_ind,rfac,novers,ixyzso)

      npts_over = ixyzso(npatches+1)-1

      allocate(srcover(12,npts_over),wover(npts_over))

      call oversample_geom(npatches,norders,ixyzs,iptype,npts, &
        srccoefs,srcvals,novers,ixyzso,npts_over,srcover)

      call get_qwts(npatches,novers,ixyzso,iptype,npts_over, &
             srcover,wover)


!
!   compute near quadrature correction
!
      nquad = iquad(nnz+1)-1
      allocate(wnear(4*nquad))
      
!$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,4*nquad
        wnear(i) = 0
      enddo
!$OMP END PARALLEL DO    

      call prinf('finished generating near field info*',i,0)
      call prinf('finished generating far field orders*',i,0)
      call prinf('npts_over=*',npts_over,1)
      call prin2('eps=*',eps,1)

      

      iquadtype = 1
      wnear = 0

      call getnearquad_lap_bel2fast(npatches,norders, &
           ixyzs,iptype,npts,srccoefs,srcvals, &
           ipatch_id,uvs_targ,eps,iquadtype,nnz,row_ptr,col_ind, &
           iquad,rfac0,nquad,wnear)
      call prinf('finished generating near quadrature correction*',i,0)

     
      call prin2('starting FAST iterative solve*',i,0)
      numit = 50
      allocate(errs(numit+1))
     
      Wg = 0 
      do i=1,npts
        sigma(i) = 0 
        alpha(i) = 0
        beta(i) = 0
      enddo


      xmin = srcvals(1,1)
      xmax = srcvals(1,1)

      ymin = srcvals(2,1)
      ymax = srcvals(2,1)

      zmin = srcvals(3,1)
      zmax = srcvals(3,1)

      do i=1,npts
        if(srcvals(1,i).gt.xmax) xmax = srcvals(1,i)
        if(srcvals(1,i).lt.xmin) xmin = srcvals(1,i)
        if(srcvals(2,i).gt.ymax) ymax = srcvals(2,i)
        if(srcvals(2,i).lt.ymin) ymin = srcvals(2,i)
        if(srcvals(3,i).gt.zmax) zmax = srcvals(3,i)
        if(srcvals(3,i).lt.zmin) zmin = srcvals(3,i)
      enddo

      xc = (xmin+xmax)/2
      yc = (ymin+ymax)/2
      zc = (zmin+zmax)/2

      xsize = xmax-xmin
      ysize = ymax-ymin
      zsize = zmax-zmin

      do i=1,4
        do j=1,3
          rpars(i,j) = hkrand(0) - 0.5d0
        enddo
      enddo

      xyz0(1) = 10.1d0
      xyz0(2) = 11.2d0
      xyz0(3) = -9.1d0
      dvec(1) = 1.1d0
      dvec(2) = -0.3d0
      dvec(3) = 0.7d0
      do i=1,npts
!!        xx = (srcvals(1,i) - xc)/xsize
!!        yy = (srcvals(2,i) - yc)/ysize
!!        zz = (srcvals(3,i) - zc)/zsize
!!        do j=1,3
!!          V(j,i) = rpars(1,j) + (rpars(2,j)*xx + rpars(3,j)*yy+ &
!!           rpars(4,j)*zz)**3/2
!!        enddo
        vtmp1(1:3) = srcvals(1:3,i) - xyz0(1:3)  
        call cross_prod3d(dvec,vtmp1,vtmp2)
        rr = sqrt(vtmp1(1)**2 + vtmp1(2)**2 + vtmp1(3)**2)
        V(1:3,i) = vtmp2/rr**3
      enddo

      do i=1,npts
        Wg = Wg + wts(i) 
      enddo
 

      print *, 'calling tangential projection'

      do i=1,npts
        call cross_cross_prod3d(srcvals(10,i),srcvals(10,i),V(1,i), &
          F(1,i))
      enddo
   
      ra = 0
      do i=1,npts
        ra=ra+F(1,i)**2*wts(i)
        ra=ra+F(2,i)**2*wts(i)
        ra=ra+F(3,i)**2*wts(i) 
      enddo
      ra = sqrt(ra)
      Fnorm = ra
      rr = 1.0d0/Fnorm 
!
!     scale F to norm 1
!
      do i=1,npts
        F(1,i)=rr*F(1,i)
        F(2,i)=rr*F(2,i)
        F(3,i)=rr*F(3,i) 
      enddo

      do i=1,npts
        call cross_prod3d(srcvals(10,i),F(1,i),nF(1,i))
      enddo



      call surf_div(npatches,norders,ixyzs,iptype,npts, &
        srccoefs,srcvals,F,rrhs1)

      call surf_div(npatches,norders,ixyzs,iptype,npts, &
        srccoefs,srcvals,nF,rrhs2)

      ra1 = 0
      ra2 = 0
      rsurf = 0
      do i=1,npts
         rrhs2(i) = -1.0d0*rrhs2(i) 
         ra1 = ra1 + rrhs1(i)*wts(i)
         ra2 = ra2 + rrhs2(i)*wts(i)
         rsurf = rsurf+wts(i)
      enddo
      call prin2('integral of boundary data1=*',ra1/rsurf,1)
      call prin2('integral of boundary data2=*',ra2/rsurf,1)
      call prin2('rsurf=*',rsurf,1)

      call prin2('rrhs1=*',rrhs1,24)
      call prin2('rrhs2=*',rrhs2,24)




      eps_gmres = 1.0d-7

      call lap_bel_solver2fast(npatches,norders,ixyzs,iptype, &
       npts,srccoefs,srcvals,eps,numit,rrhs1,eps_gmres,niter, &
       errs,rres,sigma) 


      dpars(1) = 1.0d0
      dpars(2) = 0
 

      call lpcomp_lap_comb_dir_addsub(npatches,norders,ixyzs, &
       iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,eps, &
       dpars,nnz,row_ptr,col_ind,iquad,nquad,wnear(0*nquad+1), &
       sigma,novers,npts_over,ixyzso,srcover,wover,alpha)
 
!    Surface integral should be zero
      Wu = 0 
      do i=1,npts
        Wu = Wu + alpha(i)*wts(i) 
      enddo
      call prin2('surface integral of alpha=*',Wu,1)

      erra = 0
      do i=1,npts
        erra=erra+((alpha(i))**2)*wts(i)
      enddo
      erra = sqrt(erra)
      call prin2('norm of surface alpha=*',erra,1)

      sigma = 0
      niter  = 0
      eps_gmres = 1.0d-7
      call lap_bel_solver2fast(npatches,norders,ixyzs,iptype, &
        npts,srccoefs,srcvals,eps,numit,rrhs2,eps_gmres,niter, &
        errs,rres,sigma) 


      dpars(1) = 1.0d0
      dpars(2) = 0
 

      call lpcomp_lap_comb_dir_addsub(npatches,norders,ixyzs, &
       iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,eps, &
       dpars,nnz,row_ptr,col_ind,iquad,nquad,wnear(0*nquad+1), &
       sigma,novers,npts_over,ixyzso,srcover,wover,beta)
 
!    Surface integral should be zero
      Wu = 0 
      do i=1,npts
        Wu = Wu + beta(i)*wts(i) 
      enddo
      call prin2('surface integral of beta=*',Wu,1)

      erra = 0
      ra = 0
      rr = ((1.0d0)/(4*(2*nn+1.0d0)**2)) 
      do i=1,npts
        erra=erra+((beta(i))**2)*wts(i)
      enddo
      erra = sqrt(erra)
      call prin2('norm of surface beta=*',erra,1)


      call surf_grad(npatches,norders,ixyzs,iptype,npts, &
        srccoefs,srcvals,alpha,sgalpha)

      call surf_grad(npatches,norders,ixyzs,iptype,npts, &
        srccoefs,srcvals,beta,sgbeta)

      do i=1,npts
        call cross_prod3d(srcvals(10,i),sgbeta(1,i),nsgbeta(1,i))
      enddo


      ra = 0
!$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,npts 
        H(1,i) = F(1,i)-sgalpha(1,i)-nsgbeta(1,i)
        H(2,i) = F(2,i)-sgalpha(2,i)-nsgbeta(2,i)
        H(3,i) = F(3,i)-sgalpha(3,i)-nsgbeta(3,i)
        ra = ra + (H(1,i)**2 + H(2,i)**2 + H(3,i)**2)*wts(i)
      enddo
!$OMP END PARALLEL DO     
      ra = sqrt(ra)
      call prin2('ra=*',ra,1)

      do i=1,npts
        H(1,i) = H(1,i)/ra
        H(2,i) = H(2,i)/ra
        H(3,i) = H(3,i)/ra
      enddo

      rrhs1 = 0
      call surf_div(npatches,norders,ixyzs,iptype,npts, &
        srccoefs,srcvals,H,rrhs1)
      call prin2('rrhs1=*',rrhs1,24)

      erra = 0
      do i=1,npts
        erra=  erra + ((rrhs1(i))**2)*wts(i)
      enddo
      erra = sqrt(erra)
      errest = erra
      call prin2('harmonic comp error in div =*',erra,1)

      do i=1,npts
        call cross_prod3d(srcvals(10,i),H(1,i),nsgbeta(1,i))
      enddo

      rrhs2 = 0
      call surf_div(npatches,norders,ixyzs,iptype,npts, &
        srccoefs,srcvals,nsgbeta,rrhs2)
      call prin2('rrhs2=*',rrhs2,24)

      erra = 0
      do i=1,npts
        erra=  erra + ((rrhs2(i))**2)*wts(i)
      enddo
      erra = sqrt(erra)
      if(erra.gt.errest) errest = erra
      call prin2('harmonic comp error in curl =*',erra,1)


      end subroutine get_harm_vec_field

