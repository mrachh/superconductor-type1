

      subroutine getnearquad_statj_gendeb(npatches,norders,&
       ixyzs,iptype,npts,srccoefs,srcvals, &
       eps,dpars,iquadtype,nnz,row_ptr,col_ind,&
       iquad,rfac0,nquad,wnear)
!
!  This subroutine generates the near field quadrature
!  for the representations:
!
!  wdnear: 
!    S_{0}, D_{0}, \nabla_{1} S_{0}, \nabla_{2} S_{0}, \nabla_{3} S_{0},
!    S_{0}'' + D_{0}'', S_{ik}, \nabla_{x} S_{ik}, \nabla_{y} S_{ik},
!    \nabla_{z} S_{ik}
!
!
!  wznear:
!    S_{k}, \nabla_{1} S_{k}, \nabla_{2} S_{k}, \nabla_{3} S_{k}
!
!  The quadrature is computed by the following strategy
!  targets within a sphere of radius rfac0*rs
!  of a patch centroid is handled using adaptive integration
!  where rs is the radius of the bounding sphere
!  for the patch
!  
!  All other targets in the near field are handled via
!  oversampled quadrature
!
!  The recommended parameter for rfac0 is 1.25d0
!  
!  NOTES:
!    - wdnear must be of size 10*nquad as 10 different layer
!      potentials are returned
!      * the first kernel is S_{0} 
!      * the second kernel is D_{0} 
!      * the third kernel is \nabla_{x} S_{0}
!      * the fourth kernel is \nabla_{y} S_{0}
!      * the fifth kernel is \nabla_{z} S_{0}
!      * the sixth kernel is S_{0}'' + D_{0}'
!      * the seventh kernel is S_{ik}
!      * the eighth kernel is \nabla_{x} S_{ik}
!      * the ninth kernel is \nabla_{y} S_{ik}
!      * the tenth kernel is \nabla_{z} S_{ik}
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
!    - dpars: real *8 (1)
!        the parameter k for the yukawa kernels
!    - iquadtype: integer
!        quadrature type
!          * iquadtype = 1, use ggq for self + adaptive integration
!            for rest
!    - nnz: integer
!        number of source patch-> target interactions in the near field
!    - row_ptr: integer(npts+1)
!        row_ptr(i) is the pointer
!        to col_ind array where list of relevant source patches
!        for target i start
!    - col_ind: integer (nnz)
!        list of source patches relevant for all targets, sorted
!        by the target number
!    - iquad: integer(nnz+1)
!        location in wnear_ij array where quadrature for col_ind(i)
!        starts for a single kernel. In this case the different kernels
!        are matrix entries are located at (m-1)*nquad+iquad(i), where
!        m is the kernel number
!    - rfac0: integer
!        radius parameter for near field
!    - nquad: integer
!        number of near field entries corresponding to each source target
!        pair. The size of wnear is 4*nquad since there are 4 kernels
!        per source target pair
!
!  Output arguments
!    - wnear: complex *16(4*nquad)
!        The desired near field quadrature
!               
!

      implicit none 
      integer npatches,norders(npatches),npts,nquad
      integer ixyzs(npatches+1),iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps,rfac0
      integer ndtarg,ntarg
      integer iquadtype
      integer nnz,ipars(2)
      real *8 dpars(1)
      integer row_ptr(npts+1),col_ind(nnz),iquad(nnz+1)
      real *8 wnear(10*nquad)
      real *8, allocatable :: ipatch_id(:),uvs_targ(:,:)

      complex *16 alpha,beta,ima,zk,zpars
      integer i,j,ndi,ndd,ndz

      integer ipv

      procedure (), pointer :: fker
      external l3d_slp,l3d_dlp,l3d_sgradx,l3d_sgrady,l3d_sgradz
      external y3d_slp,y3d_sgradx,y3d_sgrady,y3d_sgradz
      external l3d_spp_sum_dp


      ndz=0
      ndd=1
      ndi=0
      ndtarg = 12
      ntarg = npts

      allocate(ipatch_id(npts),uvs_targ(2,npts))
!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        ipatch_id(i) = -1
        uvs_targ(1,i) = 0
        uvs_targ(2,i) = 0
      enddo
!$OMP END PARALLEL DO      

      call get_patch_id_uvs(npatches,norders,ixyzs,iptype,npts,&
        ipatch_id,uvs_targ)

      ipv=0
      fker => l3d_slp
      call dgetnearquad_ggq_guru(npatches,norders,ixyzs,&
        iptype,npts,srccoefs,srcvals,ndtarg,ntarg,srcvals,&
        ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,&
        ndi,ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear)
      print *, "done with kernel 1"

      fker => l3d_dlp
      ipv = 0
      call dgetnearquad_ggq_guru(npatches,norders,ixyzs, &
        iptype,npts,srccoefs,srcvals,ndtarg,ntarg,srcvals, &
        ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars, &
        ndi,ipars,nnz,row_ptr,col_ind,iquad, &
        rfac0,nquad,wnear(nquad+1))
      print *, "done with kernel 2"
 1000 continue      
      fker => l3d_sgradx
      ipv = 1
      call dgetnearquad_ggq_guru(npatches,norders,ixyzs,&
       iptype,npts,srccoefs,srcvals,ndtarg,ntarg,srcvals,&
       ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,&
       ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear(2*nquad+1))

      print *, "done with kernel 3"
      fker => l3d_sgrady
      ipv = 1
      call dgetnearquad_ggq_guru(npatches,norders,ixyzs,&
       iptype,npts,srccoefs,srcvals,ndtarg,ntarg,srcvals,&
       ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,&
       ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear(3*nquad+1))

      print *, "done with kernel 4"
      fker => l3d_sgradz
      ipv = 1
      call dgetnearquad_ggq_guru(npatches,norders,ixyzs,&
       iptype,npts,srccoefs,srcvals,ndtarg,ntarg,srcvals,&
       ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,&
       ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear(4*nquad+1))

      print *, "done with kernel 5"

      fker => l3d_spp_sum_dp
      ipv = 0
      call dgetnearquad_ggq_guru(npatches,norders,ixyzs,&
       iptype,npts,srccoefs,srcvals,ndtarg,ntarg,srcvals,&
       ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,&
       ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear(5*nquad+1))

      print *, "done with kernel 6"
      goto 1111


      fker => y3d_slp
      ipv = 0
      call dgetnearquad_ggq_guru(npatches,norders,ixyzs,&
       iptype,npts,srccoefs,srcvals,ndtarg,ntarg,srcvals,&
       ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,&
       ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear(6*nquad+1))
      print *, "done with kernel 7"
      
      fker => y3d_sgradx
      ipv = 1
      call dgetnearquad_ggq_guru(npatches,norders,ixyzs,&
       iptype,npts,srccoefs,srcvals,ndtarg,ntarg,srcvals,&
       ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,&
       ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear(7*nquad+1))
      print *, "done with kernel 8"

      fker => y3d_sgrady
      ipv = 1
      call dgetnearquad_ggq_guru(npatches,norders,ixyzs,&
       iptype,npts,srccoefs,srcvals,ndtarg,ntarg,srcvals,&
       ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,&
       ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear(8*nquad+1))
      print *, "done with kernel 9"

      fker => y3d_sgradz
      ipv = 1
      call dgetnearquad_ggq_guru(npatches,norders,ixyzs,&
       iptype,npts,srccoefs,srcvals,ndtarg,ntarg,srcvals,&
       ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,&
       ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear(9*nquad+1))
      print *, "done with kernel 10"

 1111 continue       

      return
      end subroutine getnearquad_statj_gendeb
!
!
!
!
!
!


      subroutine statj_gendebproc_qpqm(npatches,norders,ixyzs, &
       iptype,npts,srccoefs,srcvals,eps,nnz,row_ptr,col_ind,iquad, &
       nquad,wnear,sigma,novers,nptso,ixyzso,srcover,whtsover,curv, &
       spqp,grads0qm,s0laps0qp,s0laps0qm)
!
!  This subroutine evaluates the following potentials:
!    spqp = S_{0}' [q^{+}]
!    grads0qm = \nabla S_{0} [qm]
!    s0laps0qp = S_{0} \Delta_{\Gamma} S_{0}[qp]
!    s0laps0qm = S_{0} \Delta_{\Gamma} S_{0}[qm]
!
!  where the near field is precomputed and stored
!  in the row sparse compressed format.
!
!  NOTES: 
!    - on output, the identity terms are not included
!
!  The fmm is used to accelerate the far-field and 
!  near-field interactions are handled via precomputed quadrature
!
!  Using add and subtract - no need to call tree and set fmm parameters
!  can directly call existing fmm library
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
!    - nnz: integer *8
!        number of source patch-> target interactions in the near field
!    - row_ptr: integer(npts+1)
!        row_ptr(i) is the pointer
!        to col_ind array where list of relevant source patches
!        for target i start
!    - col_ind: integer (nnz)
!        list of source patches relevant for all targets, sorted
!        by the target number
!    - iquad: integer(nnz+1)
!        location in wnear_ij array where quadrature for col_ind(i)
!        starts for a single kernel. In this case the different kernels
!        are matrix entries are located at (m-1)*nquad+iquad(i), where
!        m is the kernel number
!    - nquad: integer
!        number of near field entries corresponding to each source target
!        pair. The size of wnear is 4*nquad since there are 4 kernels
!        per source target pair
!    - wnear: real *8(3*nquad)
!        Precomputed near field quadrature
!    - sigma: real *8(6*npts)
!        * sigma(1:npts) is rho^{-} 
!        * sigma(npts+1:2*npts) is rho^{+} 
!        * sigma(2*npts+1:3*npts) is mu^{-} 
!        * sigma(3*npts+1:4*npts) is q^{-} 
!        * sigma(4*npts+1:5*npts) is q^{+}
!        * sigma(5*npts+1:6*npts) is r^{-}
!    - novers: integer(npatches)
!        order of discretization for oversampled sources and
!        density
!    - ixyzso: integer(npatches+1)
!        ixyzso(i) denotes the starting location in srcover,
!        corresponding to patch i
!    - nptso: integer
!        total number of oversampled points
!    - srcover: real *8 (12,nptso)
!        oversampled set of source information
!    - whtsover: real *8 (nptso)
!        smooth quadrature weights at oversampled nodes
!
!  Output arguments:
!    - spqp: real *8 (npts)
!         S_{0}'[q^{+}]
!    - grads0qm: real *8 (3,npts)
!         \nabla S_{0} [q^{-}]
!    - s0laps0qp: real *8 (npts)
!         S_{0} \Delta_{\Gamma} S_{0} [q^{+}]
!    - s0laps0qm: real *8 (npts)
!         S_{0} \Delta_{\Gamma} S_{0} [q^{-}]

!

      implicit none
      integer npatches,norder,npols,npts
      integer ndtarg,ntarg
      integer norders(npatches),ixyzs(npatches+1)
      integer ixyzso(npatches+1),iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps
      integer nnz,row_ptr(npts+1),col_ind(nnz),nquad
      integer iquad(nnz+1)
      real *8 sigma(6*npts)
      real *8 wnear(10*nquad)

      integer novers(npatches)
      integer nover,npolso,nptso
      real *8 srcover(12,nptso),whtsover(nptso)
      real *8 spqp(npts),grads0qm(3,npts),s0laps0qp(npts)
      real *8 s0laps0qm(npts)
      real *8 curv(npts)
      real *8, allocatable :: wts(:)

      real *8 rhom,rhop,rmum,uf,vf,wtmp
      real *8 u1,u2,u3,u4,w1,w2,w3,w4,w5,w0

      real *8, allocatable :: sources(:,:),srctmp(:,:)
      real *8, allocatable :: charges0(:,:),dipvec0(:,:,:)
      real *8, allocatable :: sigmaover(:,:),abc0(:,:),abc1(:,:)
      real *8, allocatable :: abc2(:,:),abc3(:,:),abc4(:,:)
      real *8, allocatable :: pot_aux(:,:),grad_aux(:,:,:)
      real *8, allocatable :: hess_aux(:,:,:)
      real *8, allocatable :: dpottmp(:),dgradtmp(:,:)
      real *8 vtmp1(3),vtmp2(3),vtmp3(3),rncj,errncj

      integer ns,nt
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      complex *16 tmp(10),val,E(4)

      integer i,j,jpatch,jquadstart,jstart,count1,count2
      real *8 radexp,epsfmm

      integer ipars(2)
      real *8 dpars(1),timeinfo(10),t1,t2,omp_get_wtime
      real *8 rqmint,rsurf

      real *8, allocatable :: radsrc(:)
      real *8, allocatable :: srctmp2(:,:)
      real *8, allocatable :: ctmp0(:,:),dtmp0(:,:,:)
      real *8 thresh,ra,erra
      real *8 rr,rmin
      real *8 over4pi
      real *8 rbl(3),rbm(3)
      integer nss,ii,l,npover
      complex *16 ima,ztmp,zpars

      integer nd,ntarg0,nmax
      integer ndd,ndz,ndi

      real *8 ttot,done,pi
      data ima/(0.0d0,1.0d0)/
      data over4pi/0.07957747154594767d0/

      parameter (ntarg0=1)

      ns = nptso
      ntarg = npts
      done = 1
      pi = atan(done)*4

      ifpgh = 0
      ifpghtarg = 2
      spqp = 0
      grads0qm = 0
      s0laps0qp = 0
      s0laps0qm = 0
      allocate(sources(3,ns),srctmp(3,npts))
      allocate(wts(npts))
      call get_qwts(npatches,norders,ixyzs,iptype,npts,&
       srcvals,wts)

      rsurf = 0
      rqmint = 0
      do i=1,npts
        rsurf = rsurf  + wts(i)
        rqmint = rqmint + sigma(3*npts+i)*wts(i)
      enddo
      rqmint = rqmint/rsurf
      call prin2('rsurf=*',rsurf,1)

!
!  estimate max number of sources in the near field of any target
!
!    
      call get_near_corr_max(npts,row_ptr,nnz,col_ind,npatches,ixyzso,&
        nmax)


!
!  extract source and target info

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ns
        sources(1,i) = srcover(1,i)
        sources(2,i) = srcover(2,i)
        sources(3,i) = srcover(3,i)
      enddo
!$OMP END PARALLEL DO      


!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        srctmp(1,i) = srcvals(1,i)
        srctmp(2,i) = srcvals(2,i)
        srctmp(3,i) = srcvals(3,i)
      enddo
!$OMP END PARALLEL DO

      nd = 2
      allocate(sigmaover(nd,ns))
      allocate(abc0(nd,npts),abc1(nd,npts),abc2(nd,npts), &
        abc3(nd,npts),abc4(nd,npts))
! 
!       oversample densities
!
      do i=1,npts
        abc0(1,i) = sigma(3*npts+i) - rqmint
        abc0(2,i) = sigma(4*npts+i)
      enddo

      call oversample_fun_surf(nd,npatches,norders,ixyzs,iptype,& 
         npts,abc0,novers,ixyzso,ns,sigmaover)



!
!  step 1: evaluate D0 (qm,qp), 
!  and D0'(qm,qp)
!  and add the corresponding potentials to their respective blocks
!
!  and S0 (qm,qp), 
!  S0' (qm,qp), and S0''(qm,qp)
!
!  abc1 = D0(qm,qp)
!  abc2 = D0'(qm,qp) + S0''(qm,qp)
!  abc3 = S0'(qm,qp)
!
      allocate(dipvec0(nd,3,ns),charges0(nd,ns))
      allocate(pot_aux(nd,npts),grad_aux(nd,3,npts), &
        hess_aux(nd,6,npts))

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ns
        charges0(1:2,i) = sigmaover(1:2,i)*whtsover(i)*over4pi
        dipvec0(1:2,1,i) = srcover(10,i)*sigmaover(1:2,i)*whtsover(i)* &
          over4pi
        dipvec0(1:2,2,i) = srcover(11,i)*sigmaover(1:2,i)*whtsover(i)* &
          over4pi
        dipvec0(1:2,3,i) = srcover(12,i)*sigmaover(1:2,i)*whtsover(i)* &
          over4pi
      enddo
!$OMP END PARALLEL DO      


      pot_aux = 0
      grad_aux = 0
      hess_aux = 0
      abc1 = 0

      call lfmm3d_t_d_g_vec(nd,eps,ns,sources,dipvec0,npts,srctmp, &
        abc1,grad_aux)

!      print *, "done with first fmm"

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        abc2(1:2,i) = grad_aux(1:2,1,i)*srcvals(10,i) + &
           grad_aux(1:2,2,i)*srcvals(11,i) + &
           grad_aux(1:2,3,i)*srcvals(12,i)
      enddo
!$OMP END PARALLEL DO

!      print *, "finished computing abc2"

      pot_aux = 0
      grad_aux = 0
      hess_aux = 0
      call lfmm3d_t_c_h_vec(nd,eps,ns,sources,charges0,npts,srctmp, &
        pot_aux,grad_aux,hess_aux)

!      print *, "done with second fmm"

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(u1,u2,u3,u4)
      do i=1,npts
        abc3(1:2,i) = grad_aux(1:2,1,i)*srcvals(10,i) + &
           grad_aux(1:2,2,i)*srcvals(11,i) + &
           grad_aux(1:2,3,i)*srcvals(12,i)
        grads0qm(1:3,i) = grad_aux(1,1:3,i)

        abc2(1:2,i) = abc2(1:2,i) + &
                    hess_aux(1:2,1,i)*srcvals(10,i)*srcvals(10,i) + &
                    hess_aux(1:2,2,i)*srcvals(11,i)*srcvals(11,i) + &
                    hess_aux(1:2,3,i)*srcvals(12,i)*srcvals(12,i) + &
                    2*hess_aux(1:2,4,i)*srcvals(11,i)*srcvals(10,i) + &
                    2*hess_aux(1:2,5,i)*srcvals(12,i)*srcvals(10,i) + &
                    2*hess_aux(1:2,6,i)*srcvals(11,i)*srcvals(12,i)
      enddo
!$OMP END PARALLEL DO

!      print *, "finished computing bll,blm,abc2,abc3"
!
!  Add in precorrected quadratures for relevant quantities
!

      call get_fmm_thresh(12,ns,srcover,12,npts,srcvals,thresh)
!
!       Add near field precomputed contribution
!
      call cpu_time(t1)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart) &
!$OMP PRIVATE(jstart,npols,l,wtmp,w1,w2,w3,w4,w5) &
!$OMP PRIVATE(w0)
      do i=1,npts
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch) 
          do l=1,npols

            w0 = wnear(jquadstart+l-1)
            w1 = wnear(nquad+jquadstart+l-1)
            w2 = wnear(2*nquad+jquadstart+l-1)
            w3 = wnear(3*nquad+jquadstart+l-1)
            w4 = wnear(4*nquad+jquadstart+l-1)
            w5 = wnear(5*nquad+jquadstart+l-1)

            grads0qm(1,i) = grads0qm(1,i) + w2*abc0(1,jstart+l-1) 
            grads0qm(2,i) = grads0qm(2,i) + w3*abc0(1,jstart+l-1) 
            grads0qm(3,i) = grads0qm(3,i) + w4*abc0(1,jstart+l-1) 

            abc1(1:2,i) =  abc1(1:2,i) + w1*abc0(1:2,jstart+l-1)
            abc2(1:2,i) =  abc2(1:2,i) + w5*abc0(1:2,jstart+l-1)

            wtmp = w2*srcvals(10,i) + w3*srcvals(11,i) + &
              w4*srcvals(12,i)

            abc3(1:2,i) =  abc3(1:2,i) + wtmp*abc0(1:2,jstart+l-1)
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO


!     Remove near contribution of the FMM
!

      allocate(ctmp0(2,nmax),dtmp0(2,3,nmax))
      allocate(dpottmp(2),dgradtmp(2,3),srctmp2(3,nmax))
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2) &
!$OMP PRIVATE(ctmp0,dtmp0,l,jstart,nss,dpottmp,dgradtmp,rr,u1)
      do i=1,npts
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          do l=ixyzso(jpatch),ixyzso(jpatch+1)-1
            nss = nss+1
            srctmp2(1,nss) = srcover(1,l)
            srctmp2(2,nss) = srcover(2,l)
            srctmp2(3,nss) = srcover(3,l)

            ctmp0(1:2,nss)=charges0(1:2,l)

            dtmp0(1:2,1:3,nss)=dipvec0(1:2,1:3,l)
            
            rr = sqrt((srcover(1,l)-srcvals(1,i))**2 + &
              (srcover(2,l)-srcvals(2,i))**2 + &
              (srcover(3,l)-srcvals(3,i))**2)
            if(rr.lt.thresh) goto 1311
            call l3d_spp_sum_dp(srcover(1,l),12,srcvals(1,i),0,dpars,&
              0,zpars,0,ipars,u1)
            abc2(1:2,i) = abc2(1:2,i) - u1*sigmaover(1:2,l)*whtsover(l)
 1311       continue            
          enddo
        enddo

        dpottmp = 0

        call l3ddirectdp(nd,srctmp2,dtmp0,nss,srctmp(1,i), &
          ntarg0,dpottmp,thresh)

        abc1(1:2,i) = abc1(1:2,i) - dpottmp(1:2)

        dpottmp = 0
        dgradtmp = 0

        call l3ddirectcg(nd,srctmp2,ctmp0,nss,srctmp(1,i), &
          ntarg0,dpottmp,dgradtmp,thresh)
        
        grads0qm(1:3,i) = grads0qm(1:3,i) - dgradtmp(1,1:3)
        
        abc3(1:2,i) = abc3(1:2,i)-(dgradtmp(1:2,1)*srcvals(10,i) + &
           dgradtmp(1:2,2)*srcvals(11,i) + &
           dgradtmp(1:2,3)*srcvals(12,i))
      enddo
!$OMP END PARALLEL DO     
      
      do i=1,npts
        spqp(i) = abc3(2,i)
      enddo

!
! Completed computation of the following quantities
!  abc1 = D0(qm,qp)
!  abc2 = D0'(qm,qp) + S0''(qm,qp)
!  abc3 = S0'(qm,qp)
!  spqp = S0'[qp]

!
!  Now compute D0 of abc1 and store in abc0 since it is no
!  longer used
!
      call oversample_fun_surf(nd,npatches,norders,ixyzs,iptype,& 
          npts,abc1,novers,ixyzso,ns,sigmaover)

!$OMP PARALLEL DO DEFAULT(SHARED) 
      do i=1,ns
        dipvec0(1:2,1,i) = srcover(10,i)*sigmaover(1:2,i)*whtsover(i)* &
          over4pi
        dipvec0(1:2,2,i) = srcover(11,i)*sigmaover(1:2,i)*whtsover(i)* &
          over4pi
        dipvec0(1:2,3,i) = srcover(12,i)*sigmaover(1:2,i)*whtsover(i)* &
          over4pi
      enddo
!$OMP END PARALLEL DO      

      pot_aux = 0
      abc0 = 0

      call lfmm3d_t_d_p_vec(nd,eps,ns,sources,dipvec0,npts,srctmp, &
        abc0)
!
!  Add near quadrature correction
!

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart) &
!$OMP PRIVATE(jstart,npols,l,w1)
      do i=1,npts
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch) 
          do l=1,npols

            w1 = wnear(nquad+jquadstart+l-1)
            abc0(1:2,i) = abc0(1:2,i) + w1*abc1(1:2,jstart+l-1)
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO      

!
! Subtract near contributions computed via fmm
!
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2) &
!$OMP PRIVATE(dtmp0,l,jstart,nss,dpottmp)
      do i=1,npts
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          do l=ixyzso(jpatch),ixyzso(jpatch+1)-1
            nss = nss+1
            srctmp2(1,nss) = srcover(1,l)
            srctmp2(2,nss) = srcover(2,l)
            srctmp2(3,nss) = srcover(3,l)

            dtmp0(1:2,1:3,nss)=dipvec0(1:2,1:3,l)
          enddo
        enddo

        dpottmp = 0

        call l3ddirectdp(nd,srctmp2,dtmp0,nss,srctmp(1,i), &
          ntarg0,dpottmp,thresh)

        abc0(1:2,i) = abc0(1:2,i) - dpottmp(1:2)
      enddo
!$OMP END PARALLEL DO     



!
!  End of computing abc0 = D0 (abc1) = D0^2[qm,qp]
!
!
!
!  The 6 densities now being considered are ordered as follows:
!  
!  abc4(1) = ((D0'+S0'') + 2HS0')[qm]
!     + 1/|\Gamma|\int_{\Gamma} S_{0} [\qm]
!  abc4(2) = ((D0'+S0'') + 2HS0')[qp]
!     + 1/|\Gamma|\int_{\Gamma} S_{0} [\qp]
!
!  We first organize these densities and oversample them
!
!
      do i=1,npts
        abc4(1:2,i) = abc2(1:2,i) + 2*curv(i)*abc3(1:2,i)
      enddo


      sigmaover=  0

      call oversample_fun_surf(nd,npatches,norders,ixyzs,iptype,& 
          npts,abc4,novers,ixyzso,ns,sigmaover)

!
!$OMP PARALLEL DO DEFAULT(SHARED) 
      do i=1,ns
        charges0(1:2,i) = sigmaover(1:2,i)*whtsover(i)*over4pi
      enddo
!$OMP END PARALLEL DO      

      pot_aux = 0

      call lfmm3d_t_c_p_vec(nd,eps,ns,sources,charges0,npts,srctmp, &
        pot_aux)


!
!  Add near quadrature correction
!

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart) &
!$OMP PRIVATE(jstart,npols,l,w1)
      do i=1,npts
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch) 
          do l=1,npols
            w1 = wnear(jquadstart+l-1)
            pot_aux(1:2,i) = pot_aux(1:2,i) + w1*abc4(1:2,jstart+l-1)
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO      

!
! Subtract near contributions computed via fmm
!
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2) &
!$OMP PRIVATE(ctmp0,l,jstart,nss,dpottmp)
      do i=1,npts
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          do l=ixyzso(jpatch),ixyzso(jpatch+1)-1
            nss = nss+1
            srctmp2(1,nss) = srcover(1,l)
            srctmp2(2,nss) = srcover(2,l)
            srctmp2(3,nss) = srcover(3,l)

            ctmp0(1:2,nss)=charges0(1:2,l)
          enddo
        enddo

        dpottmp = 0

        call l3ddirectcp(nd,srctmp2,ctmp0,nss,srctmp(1,i), &
          ntarg0,dpottmp,thresh)

        pot_aux(1:2,i) = pot_aux(1:2,i) - dpottmp(1:2)
      enddo
!$OMP END PARALLEL DO      

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        s0laps0qp(i) = abc0(2,i) - pot_aux(2,i)
        s0laps0qm(i) = abc0(1,i) - pot_aux(1,i)
      enddo
!$OMP END PARALLEL DO


      return
      end subroutine statj_gendebproc_qpqm
      





      subroutine statj_gendebproc_rhomrhopmum(npatches,norders,ixyzs, &
       iptype,npts,srccoefs,srcvals,eps,nnz,row_ptr,col_ind,iquad, &
       nquad,wnear,sigma,novers,nptso,ixyzso,srcover,whtsover,curv, &
       wtmp1,wtmp2,wtmp3,wtmp4,dzk,rbeta,rgamma,laps02rhom,laps02rhop, &
       laps02mum,blm,bmm)
!
!  This subroutine evaluates the following potentials:
!    laps02rhom = \Delta_{\Gamma} S_{0}^2 [\rho^{-}]
!    laps02rhop = \Delta_{\Gamma} S_{0}^2 [\rho^{+}]
!    laps02mum = \Delta_{\Gamma} S_{0}^2 [\mu^{-}]
!    blm = -dzk*(\nabla_{\Gamma} S_{0}^2 [\mu^{-}] + rbeta* n \times
!       \nabla_{\Gamma} S_{0}^2 [\rho^{+}])
!    bmm = dzk*(\nabla_{\Gamma} S_{0}^2 [\rho^{-}] + rgamma* n \times
!       \nabla_{\Gamma} S_{0}^2 [\rho^{+}])
!
!  where the near field is precomputed and stored
!  in the row sparse compressed format.
!
!  NOTES: 
!    - on output, the identity terms are not included
!
!  The fmm is used to accelerate the far-field and 
!  near-field interactions are handled via precomputed quadrature
!
!  Using add and subtract - no need to call tree and set fmm parameters
!  can directly call existing fmm library
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
!    - nnz: integer *8
!        number of source patch-> target interactions in the near field
!    - row_ptr: integer(npts+1)
!        row_ptr(i) is the pointer
!        to col_ind array where list of relevant source patches
!        for target i start
!    - col_ind: integer (nnz)
!        list of source patches relevant for all targets, sorted
!        by the target number
!    - iquad: integer(nnz+1)
!        location in wnear_ij array where quadrature for col_ind(i)
!        starts for a single kernel. In this case the different kernels
!        are matrix entries are located at (m-1)*nquad+iquad(i), where
!        m is the kernel number
!    - nquad: integer
!        number of near field entries corresponding to each source target
!        pair. The size of wnear is 4*nquad since there are 4 kernels
!        per source target pair
!    - wnear: real *8(3*nquad)
!        Precomputed near field quadrature
!    - sigma: real *8(6*npts)
!        * sigma(1:npts) is rho^{-} 
!        * sigma(npts+1:2*npts) is rho^{+} 
!        * sigma(2*npts+1:3*npts) is mu^{-} 
!        * sigma(3*npts+1:4*npts) is q^{-} 
!        * sigma(4*npts+1:5*npts) is q^{+}
!        * sigma(5*npts+1:6*npts) is r^{-}
!    - novers: integer(npatches)
!        order of discretization for oversampled sources and
!        density
!    - ixyzso: integer(npatches+1)
!        ixyzso(i) denotes the starting location in srcover,
!        corresponding to patch i
!    - nptso: integer
!        total number of oversampled points
!    - srcover: real *8 (12,nptso)
!        oversampled set of source information
!    - whtsover: real *8 (nptso)
!        smooth quadrature weights at oversampled nodes
!    - curv: real *8 (npts)
!        mean curvature at the discretization nodes
!    - wtmp1-4: real *8 (3,npts)
!        vectors for computing contribution to 
!        blm and bmm
!    - dzk,rbeta,rgamma: real *8
!       skin depth and parameters for contribution of 
!       q^{+} in the two currents
!
!  Output arguments:
!    - laps02rhom: real *8 (npts)
!         \Delta_{\Gamma} S_{0}^2 [rho^{-}]
!    - laps02rhop: real *8 (npts)
!         \Delta_{\Gamma} S_{0}^2 [rho^{+}]
!    - laps02mum: real *8 (npts)
!         \Delta_{\Gamma} S_{0}^2 [mu^{-}]
!    - blm: real *8 (3,npts)
!         See definition above
!    - bmm: real *8 (3,npts)
!         See definition above

      implicit none
      integer npatches,norder,npols,npts
      integer ndtarg,ntarg
      integer norders(npatches),ixyzs(npatches+1)
      integer ixyzso(npatches+1),iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps
      integer nnz,row_ptr(npts+1),col_ind(nnz),nquad
      integer iquad(nnz+1)
      real *8 sigma(6*npts)
      real *8 wnear(10*nquad)

      integer novers(npatches)
      integer nover,npolso,nptso
      real *8 srcover(12,nptso),whtsover(nptso)
      real *8 curv(npts)
      real *8 wtmp1(3,npts),wtmp2(3,npts),wtmp3(3,npts)
      real *8 wtmp4(3,npts)
      real *8 rbeta,rgamma,dzk
      real *8 laps02rhom(npts),laps02rhop(npts)
      real *8 laps02mum(npts),blm(3,npts),bmm(3,npts)
      real *8, allocatable :: wts(:)

      real *8 rhom,rhop,rmum,uf,vf,wtmp
      real *8 srhop,srhom,srmum
      real *8 u1,u2,u3,u4,u5,u6,w1,w2,w3,w4,w5,w0

      real *8, allocatable :: sources(:,:),srctmp(:,:)
      real *8, allocatable :: charges0(:,:),dipvec0(:,:,:)
      real *8, allocatable :: sigmaover(:,:),abc0(:,:),abc1(:,:)
      real *8, allocatable :: abc2(:,:),abc3(:,:),abc4(:,:)
      real *8, allocatable :: pot_aux(:,:),grad_aux(:,:,:)
      real *8, allocatable :: hess_aux(:,:,:)
      real *8, allocatable :: dpottmp(:),dgradtmp(:,:)
      real *8 vtmp1(3),vtmp2(3),vtmp3(3),rncj,errncj

      integer ns,nt
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      complex *16 tmp(10),val,E(4)

      integer i,j,jpatch,jquadstart,jstart,count1,count2
      real *8 radexp,epsfmm

      integer ipars(2)
      real *8 dpars(1),timeinfo(10),t1,t2,omp_get_wtime
      real *8 rqmint,rsurf

      real *8, allocatable :: radsrc(:)
      real *8, allocatable :: srctmp2(:,:)
      real *8, allocatable :: ctmp0(:,:),dtmp0(:,:,:)
      real *8 thresh,ra,erra
      real *8 rr,rmin
      real *8 over4pi
      real *8 rbl(3),rbm(3)
      integer nss,ii,l,npover
      complex *16 ima,ztmp,zpars

      integer nd,ntarg0,nmax
      integer ndd,ndz,ndi

      real *8 ttot,done,pi
      data ima/(0.0d0,1.0d0)/
      data over4pi/0.07957747154594767d0/

      parameter (ntarg0=1)


      ns = nptso
      ntarg = npts
      done = 1
      pi = atan(done)*4

      ifpgh = 0
      ifpghtarg = 2
      allocate(sources(3,ns),srctmp(3,npts))
      allocate(wts(npts))
      call get_qwts(npatches,norders,ixyzs,iptype,npts,&
       srcvals,wts)

!
!  estimate max number of sources in the near field of any target
!
!    
      call get_near_corr_max(npts,row_ptr,nnz,col_ind,npatches,ixyzso,&
        nmax)


!
!  extract source and target info

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ns
        sources(1,i) = srcover(1,i)
        sources(2,i) = srcover(2,i)
        sources(3,i) = srcover(3,i)
      enddo
!$OMP END PARALLEL DO      


!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        srctmp(1,i) = srcvals(1,i)
        srctmp(2,i) = srcvals(2,i)
        srctmp(3,i) = srcvals(3,i)
      enddo
!$OMP END PARALLEL DO

      nd = 3
      allocate(sigmaover(nd,ns))
      allocate(abc0(nd,npts),abc1(nd,npts),abc2(nd,npts), &
        abc3(nd,npts),abc4(nd,npts))
! 
!       oversample densities
!
      do i=1,npts
        abc0(1,i) = sigma(i) 
        abc0(2,i) = sigma(npts+i)
        abc0(3,i) = sigma(2*npts+i)
      enddo

      call oversample_fun_surf(nd,npatches,norders,ixyzs,iptype,& 
         npts,abc0,novers,ixyzso,ns,sigmaover)



!
!
!  abc1 = S0(rhom,rhop,mum)
!  abc2 = S0'(rhom,,rhop,mum)
!
      allocate(dipvec0(nd,3,ns),charges0(nd,ns))
      allocate(pot_aux(nd,npts),grad_aux(nd,3,npts), &
        hess_aux(nd,6,npts))

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ns
        charges0(1:3,i) = sigmaover(1:3,i)*whtsover(i)*over4pi
      enddo
!$OMP END PARALLEL DO      


      pot_aux = 0
      grad_aux = 0
      abc1 = 0

      call lfmm3d_t_c_g_vec(nd,eps,ns,sources,charges0,npts,srctmp, &
        abc1,grad_aux)

!      print *, "done with first fmm"

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        abc2(1:3,i) = grad_aux(1:3,1,i)*srcvals(10,i) + &
           grad_aux(1:3,2,i)*srcvals(11,i) + &
           grad_aux(1:3,3,i)*srcvals(12,i)
      enddo
!$OMP END PARALLEL DO

!
!  Add in precorrected quadratures for relevant quantities
!

      call get_fmm_thresh(12,ns,srcover,12,npts,srcvals,thresh)
!
!       Add near field precomputed contribution
!
      call cpu_time(t1)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart) &
!$OMP PRIVATE(jstart,npols,l,wtmp,w1,w2,w3,w4,w5) &
!$OMP PRIVATE(w0)
      do i=1,npts
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch) 
          do l=1,npols

            w0 = wnear(jquadstart+l-1)
            w1 = wnear(nquad+jquadstart+l-1)
            w2 = wnear(2*nquad+jquadstart+l-1)
            w3 = wnear(3*nquad+jquadstart+l-1)
            w4 = wnear(4*nquad+jquadstart+l-1)
            w5 = wnear(5*nquad+jquadstart+l-1)


            abc1(1:3,i) =  abc1(1:3,i) + w0*abc0(1:3,jstart+l-1)
            wtmp = w2*srcvals(10,i) + w3*srcvals(11,i) + &
              w4*srcvals(12,i)
            abc2(1:3,i) =  abc2(1:3,i) + wtmp*abc0(1:3,jstart+l-1)
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO


!     Remove near contribution of the FMM
!

      allocate(ctmp0(3,nmax),dtmp0(3,3,nmax))
      allocate(dpottmp(3),dgradtmp(3,3),srctmp2(3,nmax))
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2) &
!$OMP PRIVATE(ctmp0,l,jstart,nss,dpottmp,dgradtmp)
      do i=1,npts
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          do l=ixyzso(jpatch),ixyzso(jpatch+1)-1
            nss = nss+1
            srctmp2(1,nss) = srcover(1,l)
            srctmp2(2,nss) = srcover(2,l)
            srctmp2(3,nss) = srcover(3,l)

            ctmp0(1:3,nss)=charges0(1:3,l)
          enddo
        enddo

        dpottmp = 0
        dgradtmp = 0

        call l3ddirectcg(nd,srctmp2,ctmp0,nss,srctmp(1,i), &
          ntarg0,dpottmp,dgradtmp,thresh)
        
        abc1(1:3,i) = abc1(1:3,i) - dpottmp(1:3)
        abc2(1:3,i) = abc2(1:3,i) - (dgradtmp(1:3,1)*srcvals(10,i) + &
           dgradtmp(1:3,2)*srcvals(11,i) + &
           dgradtmp(1:3,3)*srcvals(12,i))
      enddo
!$OMP END PARALLEL DO    

!
! Completed computation of the following quantities
!
! abc1 = S0[rhom,rhop,mum]
! abc2 = S0'[rhom,rhop,mum]
!
      


!
!  Now compute S0' of abc2 and store in abc0 since it is no
!  longer used
!
      call oversample_fun_surf(nd,npatches,norders,ixyzs,iptype,& 
          npts,abc2,novers,ixyzso,ns,sigmaover)

!$OMP PARALLEL DO DEFAULT(SHARED) 
      do i=1,ns
        charges0(1:3,i) = sigmaover(1:3,i)*whtsover(i)*over4pi
      enddo
!$OMP END PARALLEL DO      

      pot_aux = 0
      abc0 = 0
      grad_aux = 0

      call lfmm3d_t_c_g_vec(nd,eps,ns,sources,charges0,npts,srctmp, &
        pot_aux,grad_aux)

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        abc0(1:3,i) = grad_aux(1:3,1,i)*srcvals(10,i) + &
           grad_aux(1:3,2,i)*srcvals(11,i) + &
           grad_aux(1:3,3,i)*srcvals(12,i)
      enddo
!$OMP END PARALLEL DO

!
!  Add near quadrature correction
!

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart) &
!$OMP PRIVATE(jstart,npols,l,w1,w2,w3,w4,wtmp)
      do i=1,npts
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch) 
          do l=1,npols
            w1 = wnear(nquad+jquadstart+l-1)
            w2 = wnear(2*nquad+jquadstart+l-1)
            w3 = wnear(3*nquad+jquadstart+l-1)
            w4 = wnear(4*nquad+jquadstart+l-1)
            wtmp = w2*srcvals(10,i) + w3*srcvals(11,i) + &
              w4*srcvals(12,i)
            abc0(1:3,i) =  abc0(1:3,i) + wtmp*abc2(1:3,jstart+l-1)
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO      

!
! Subtract near contributions computed via fmm
!
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2) &
!$OMP PRIVATE(ctmp0,l,jstart,nss,dpottmp,dgradtmp)
      do i=1,npts
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          do l=ixyzso(jpatch),ixyzso(jpatch+1)-1
            nss = nss+1
            srctmp2(1,nss) = srcover(1,l)
            srctmp2(2,nss) = srcover(2,l)
            srctmp2(3,nss) = srcover(3,l)

            ctmp0(1:3,nss)=charges0(1:3,l)
          enddo
        enddo

        dpottmp = 0
        dgradtmp = 0

        call l3ddirectcg(nd,srctmp2,ctmp0,nss,srctmp(1,i), &
          ntarg0,dpottmp,dgradtmp,thresh)
        abc0(1:3,i) = abc0(1:3,i)-(dgradtmp(1:3,1)*srcvals(10,i) + &
           dgradtmp(1:3,2)*srcvals(11,i) + &
           dgradtmp(1:3,3)*srcvals(12,i))

      enddo
!$OMP END PARALLEL DO     



!
!  End of computing abc0 = S0' (abc2) = S0'^2[rhom,rhop,mum]
!  abc1 = S0(rhom,rhop,mum)
!
!
!  abc4 = (S'' + D') + 2HS' (abc1)
!
!  
!
!

      sigmaover=  0

      call oversample_fun_surf(nd,npatches,norders,ixyzs,iptype,& 
          npts,abc1,novers,ixyzso,ns,sigmaover)

!
!$OMP PARALLEL DO DEFAULT(SHARED) 
      do i=1,ns
        charges0(1:3,i) = sigmaover(1:3,i)*whtsover(i)*over4pi
        dipvec0(1:3,1,i) = srcover(10,i)*sigmaover(1:3,i)*whtsover(i)* &
          over4pi
        dipvec0(1:3,2,i) = srcover(11,i)*sigmaover(1:3,i)*whtsover(i)* &
          over4pi
        dipvec0(1:3,3,i) = srcover(12,i)*sigmaover(1:3,i)*whtsover(i)* &
          over4pi
      enddo
!$OMP END PARALLEL DO      

      pot_aux = 0
      grad_aux = 0

      call lfmm3d_t_d_g_vec(nd,eps,ns,sources,dipvec0,npts,srctmp, &
        pot_aux, grad_aux)

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        abc4(1:3,i) = grad_aux(1:3,1,i)*srcvals(10,i) + &
           grad_aux(1:3,2,i)*srcvals(11,i) + &
           grad_aux(1:3,3,i)*srcvals(12,i)
      enddo
!$OMP END PARALLEL DO


      pot_aux = 0
      grad_aux = 0
      hess_aux = 0
      call lfmm3d_t_c_h_vec(nd,eps,ns,sources,charges0,npts,srctmp, &
        pot_aux,grad_aux,hess_aux)


!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(u1,u2,u3,u4,u5,u6)
      do i=1,npts
        abc4(1:3,i) = abc4(1:3,i) + (grad_aux(1:3,1,i)*srcvals(10,i) + &
           grad_aux(1:3,2,i)*srcvals(11,i) + &
           grad_aux(1:3,3,i)*srcvals(12,i))*curv(i)*2

        abc4(1:3,i) = abc4(1:3,i) + &
                    hess_aux(1:3,1,i)*srcvals(10,i)*srcvals(10,i) + &
                    hess_aux(1:3,2,i)*srcvals(11,i)*srcvals(11,i) + &
                    hess_aux(1:3,3,i)*srcvals(12,i)*srcvals(12,i) + &
                    2*hess_aux(1:3,4,i)*srcvals(11,i)*srcvals(10,i) + &
                    2*hess_aux(1:3,5,i)*srcvals(12,i)*srcvals(10,i) + &
                    2*hess_aux(1:3,6,i)*srcvals(11,i)*srcvals(12,i)

         call dot_prod3d(grad_aux(1,1:3,i),srcvals(4,i),u1)
         call dot_prod3d(grad_aux(1,1:3,i),srcvals(7,i),u2)

         call dot_prod3d(grad_aux(2,1:3,i),srcvals(4,i),u3)
         call dot_prod3d(grad_aux(2,1:3,i),srcvals(7,i),u4)

         call dot_prod3d(grad_aux(3,1:3,i),srcvals(4,i),u5)
         call dot_prod3d(grad_aux(3,1:3,i),srcvals(7,i),u6)

         blm(1:3,i) = -(u5*wtmp1(1:3,i) + u6*wtmp2(1:3,i))*dzk + &
            rbeta*(u3*wtmp3(1:3,i) + u4*wtmp4(1:3,i))
         bmm(1:3,i) = (u1*wtmp1(1:3,i) + u2*wtmp2(1:3,i))*dzk + &
            rgamma*(u3*wtmp3(1:3,i) + u4*wtmp4(1:3,i))
      enddo
!$OMP END PARALLEL DO

!
!       Add near field precomputed contribution
!
      call cpu_time(t1)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart) &
!$OMP PRIVATE(jstart,npols,l,wtmp,uf,vf,srhop,srhom,srmum,w1,w2,w3,w4,w5) &
!$OMP PRIVATE(w0)
      do i=1,npts
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch) 
          do l=1,npols
            srhom = abc1(1,jstart+l-1)
            srhop = abc1(2,jstart+l-1)
            srmum = abc1(3,jstart+l-1)

            w0 = wnear(jquadstart+l-1)
            w1 = wnear(nquad+jquadstart+l-1)
            w2 = wnear(2*nquad+jquadstart+l-1)
            w3 = wnear(3*nquad+jquadstart+l-1)
            w4 = wnear(4*nquad+jquadstart+l-1)
            w5 = wnear(5*nquad+jquadstart+l-1)

            abc4(1:3,i) =  abc4(1:3,i) + w5*abc1(1:3,jstart+l-1)

            wtmp = w2*srcvals(10,i) + w3*srcvals(11,i) + &
              w4*srcvals(12,i)
            abc4(1:3,i) =  abc4(1:3,i) + &
               wtmp*abc1(1:3,jstart+l-1)*2*curv(i)

            uf = w2*srcvals(4,i) + w3*srcvals(5,i) + &
              w4*srcvals(6,i)
            vf = w2*srcvals(7,i) + w3*srcvals(8,i) + &
              w4*srcvals(9,i)

            blm(1:3,i) = blm(1:3,i) - dzk*(uf*wtmp1(1:3,i) + &
              vf*wtmp2(1:3,i))*srmum + rbeta*(uf*wtmp3(1:3,i) + &
                vf*wtmp4(1:3,i))*srhop

            bmm(1:3,i) = bmm(1:3,i) + dzk*(uf*wtmp1(1:3,i) + &
              vf*wtmp2(1:3,i))*srhom + rgamma*(uf*wtmp3(1:3,i) + &
                vf*wtmp4(1:3,i))*srhop
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2) &
!$OMP PRIVATE(ctmp0,l,jstart,nss,dpottmp,dgradtmp,u1,u2,u3,u4,u5,u6,rr)
      do i=1,npts
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          do l=ixyzso(jpatch),ixyzso(jpatch+1)-1
            nss = nss+1
            srctmp2(1,nss) = srcover(1,l)
            srctmp2(2,nss) = srcover(2,l)
            srctmp2(3,nss) = srcover(3,l)

            ctmp0(1:3,nss)=charges0(1:3,l)

            
            rr = sqrt((srcover(1,l)-srcvals(1,i))**2 + &
              (srcover(2,l)-srcvals(2,i))**2 + &
              (srcover(3,l)-srcvals(3,i))**2)
            if(rr.lt.thresh) goto 1311
            u1 = 0
            call l3d_spp_sum_dp(srcover(1,l),12,srcvals(1,i),0,dpars,&
              0,zpars,0,ipars,u1)
            abc4(1:3,i) = abc4(1:3,i) - u1*sigmaover(1:3,l)*whtsover(l)
 1311       continue            
          enddo
        enddo

        dpottmp = 0
        dgradtmp = 0

        call l3ddirectcg(nd,srctmp2,ctmp0,nss,srctmp(1,i), &
          ntarg0,dpottmp,dgradtmp,thresh)
        
        abc4(1:3,i) = abc4(1:3,i)-(dgradtmp(1:3,1)*srcvals(10,i) + &
           dgradtmp(1:3,2)*srcvals(11,i) + &
           dgradtmp(1:3,3)*srcvals(12,i))*2*curv(i)

         call dot_prod3d(dgradtmp(3,1:3),srcvals(4,i),u5)
         call dot_prod3d(dgradtmp(3,1:3),srcvals(7,i),u6)

         call dot_prod3d(dgradtmp(1,1:3),srcvals(4,i),u1)
         call dot_prod3d(dgradtmp(1,1:3),srcvals(7,i),u2)

         call dot_prod3d(dgradtmp(2,1:3),srcvals(4,i),u3)
         call dot_prod3d(dgradtmp(2,1:3),srcvals(7,i),u4)


         blm(1:3,i) =blm(1:3,i)+(u5*wtmp1(1:3,i) + u6*wtmp2(1:3,i))*dzk 
         blm(1:3,i) =blm(1:3,i)-(u3*wtmp3(1:3,i) + u4*wtmp4(1:3,i))*rbeta

         bmm(1:3,i) =bmm(1:3,i)-(u1*wtmp1(1:3,i) + u2*wtmp2(1:3,i))*dzk 
         bmm(1:3,i) =bmm(1:3,i)-(u3*wtmp3(1:3,i) + u4*wtmp4(1:3,i))*rgamma

      enddo
!$OMP END PARALLEL DO     

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        laps02rhom(i) = -abc4(1,i) + abc0(1,i)
        laps02rhop(i) = -abc4(2,i) + abc0(2,i)
        laps02mum(i) = -abc4(3,i) + abc0(3,i)
      enddo
!$OMP END PARALLEL DO

      call prin2('abc4=*',abc4,24)
      call prin2('abc0=*',abc0,24)


      return
      end subroutine statj_gendebproc_rhomrhopmum
      


