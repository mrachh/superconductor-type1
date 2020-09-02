!  This file contains the following user callable routines
!
!  getnearquad_statj_gendeb - generate all the near quadrature
!  corrections needed for the various layer potentials
!
!

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

      subroutine lpcomp_statj_gendeb_addsub(npatches,norders,ixyzs,&
        iptype,npts,srccoefs,srcvals,eps,dzk,nnz,row_ptr,col_ind, &
        iquad,nquad,wnear,sigma,novers,nptso,ixyzso,srcover,whtsover, &
        pot)

!
!  This subroutine evaluates the sets of potentials to the following
!  representations:
!
!  1. -4 S_{0} \Delta_{\Gamma} S_{0} \rho^- = -4 S_{0} [q^-]
!  2. -4 S_{0} \Delta_{\Gamma} S_{0} \rho^{+} = -4 S_{0} [q^+]
!  3. -4 S_{0} \Delta_{\Gamma} S_{0} \mu^{-} = -4 S_{0} [r^{-}]
!  a. -n \cdot (S_{ik} [m-] + 1/k* \nabla S_{ik}[q^{-}] 
!     +  \nabla \times S_{ik}[\ell^{-}]  = B\cdot n
!  b. \nabla \cdot S_{0} [ \bn \times \bn \times (S_{ik} [m^{-}] + 
!     -1/k*\nabla S_{ik} [q^{-}] -1/k* \nabla \times S_{ik} + 
!      \nabla S_{0}[q^{+}]] = -\nabla \cdot S_{0} [\bn \times \bn \times B]
!  c. m^{-} = k \nabla_{\Gamma} S_{0} [\rho^{-}] 
!     -k n \times \nabla_{\Gamma}S_{0} [\rho^{+}]
!  d. \ell^{-} = -k \nabla_{\Gamma} S_{0}[\mu^{-}]
!  4. (a+b)*2*k
!  5. (a-b)
!  6. 2*(n \cdot (-\nabla S_{ik}[r^{-}] - \nabla \times S_{ik} [m^{-}] -k 
!        S_{ik} [\ell^{-}]) = 0 
!
!
!
!  where the near field is precomputed and stored
!  in the row sparse compressed format.
!
!  NOTES: 
!    - on output, the identity terms are not included
!  Input arguments:
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
!    - dzk: real *8 (1)
!        the parameter k for the yukawa kernels
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
!    - wnear: real *8(10*nquad)
!        Precomputed near field quadrature
!          * the first kernel is S_{0} 
!          * the second kernel is D_{0} 
!          * the third kernel is \nabla_{x} S_{0}
!          * the fourth kernel is \nabla_{y} S_{0}
!          * the fifth kernel is \nabla_{z} S_{0}
!          * the sixth kernel is S_{0}'' + D_{0}'
!          * the seventh kernel is S_{ik}
!          * the eighth kernel is \nabla_{x} S_{ik}
!          * the ninth kernel is \nabla_{y} S_{ik}
!          * the tenth kernel is \nabla_{z} S_{ik}
!    - sigma: real *8(6*npts)
!        * sigma(1:npts) is the density \rho^-
!        * sigma(npts+1:2*npts) is the density \rho^+
!        * sigma(2*npts+1:3*npts) is the density \mu^-
!        * sigma(3*npts+1:4*npts) is the desnity q^-
!        * sigma(4*npts+1:5*npts) is the desnity q^+
!        * sigma(5*npts+1:6*npts) is the density r^-
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
!    - pot: real *8 (6*npts)
!

      implicit none
      integer npatches,norder,npols,npts
      integer ndtarg,ntarg
      integer norders(npatches),ixyzs(npatches+1)
      integer ixyzso(npatches+1),iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps
      real *8 dzk
      complex *16 zpars(5),zpars_tmp(4)
      complex *16 zk0,zk1
      integer nnz,row_ptr(npts+1),col_ind(nnz),nquad
      integer iquad(nnz+1)
      real *8 sigma(6*npts)
      real *8 wnear(10*nquad)

      integer novers(npatches)
      integer nover,npolso,nptso
      real *8 srcover(12,nptso),whtsover(nptso)
      real *8 pot(6*npts)
      real *8, allocatable :: wts(:)

      real *8 rhom,rhop,rmum,uf,vf,wtmp
      real *8 u1,u2,u3,u4,w1,w2,w3,w4,w5

      real *8, allocatable :: sources(:,:),srctmp(:,:)
      real *8, allocatable :: charges0(:,:),dipvec0(:,:,:)
      real *8, allocatable :: sigmaover(:,:)
      real *8, allocatable :: sigmause(:,:)
      real *8, allocatable :: abc1(:,:),abc2(:,:),abc3(:,:)
      real *8, allocatable :: abc4(:,:),abc0(:,:)
      complex *16, allocatable :: zcharges0(:,:)
      real *8, allocatable :: pot_aux(:,:),grad_aux(:,:,:)
      real *8, allocatable :: hess_aux(:,:,:)
      complex *16, allocatable :: zpot_aux(:,:),zgrad_aux(:,:,:)
      real *8, allocatable :: bmm(:,:),blm(:,:)
      real *8, allocatable :: bjm(:,:),bbm(:,:),bbp(:,:)
      real *8, allocatable :: ffforminv(:,:,:),curv(:)
      real *8, allocatable :: wtmp1(:,:),wtmp2(:,:),wtmp3(:,:), &
        wtmp4(:,:)
      real *8, allocatable :: dpottmp(:),dgradtmp(:,:)
      complex *16, allocatable :: zpottmp(:),zgradtmp(:,:)



      
      integer ns,nt
      complex *16 alpha,beta
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      complex *16 tmp(10),val,E(4)

      integer i,j,jpatch,jquadstart,jstart,count1,count2
      complex *16 zdotu,pottmp,gradtmp(3)
      complex *16 ep0,ep1,ep0sq,ep1sq,ep0inv,ep1inv

      real *8 radexp,epsfmm

      integer ipars(2)
      real *8 dpars(1),timeinfo(10),t1,t2,omp_get_wtime

      real *8, allocatable :: radsrc(:)
      real *8, allocatable :: srctmp2(:,:)
      real *8, allocatable :: ctmp0(:,:),dtmp0(:,:,:)
      complex *16, allocatable :: zctmp0(:,:),zdtmp0(:,:,:)
      real *8 thresh,ra,erra
      real *8 rr,rmin
      real *8 over4pi
      real *8 rbl(3),rbm(3)
      integer nss,ii,l,npover
      complex *16 ima,ztmp
      complex *16 jvals(100),hvals(100),fjder(100),fhder(100)
      real *8 errbl,errbm,errjn,rfac,rjn,rscale,rrjn,rtmp
      integer ifder,njh,n

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
      pot = 0
      allocate(sources(3,ns),srctmp(3,npts))
!
!  Compute inverse of first fundamental form and mean curvature 
!
      allocate(curv(npts),ffforminv(2,2,npts))
      call get_mean_curvature(npatches,norders,ixyzs,iptype,npts, &
        srccoefs,srcvals,curv)
      call get_inv_first_fundamental_form(npatches,norders,ixyzs, &
        iptype,npts,srccoefs,srcvals,ffforminv)
      
      allocate(wtmp1(3,npts),wtmp2(3,npts),wtmp3(3,npts),wtmp4(3,npts))

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
         wtmp1(1:3,i) = ffforminv(1,1,i)*srcvals(4:6,i) + &
           ffforminv(1,2,i)*srcvals(7:9,i)
         wtmp2(1:3,i) = ffforminv(2,1,i)*srcvals(4:6,i) + &
           ffforminv(2,2,i)*srcvals(7:9,i)
         call cross_prod3d(srcvals(10,i),wtmp1(1,i),wtmp3(1,i))
         call cross_prod3d(srcvals(10,i),wtmp2(1,i),wtmp4(1,i))
      enddo
!$OMP END PARALLEL DO      



!
!  estimate max number of sources in the near field of any target
!
!    
      call get_near_corr_max(npts,row_ptr,nnz,col_ind,npatches,ixyzso,&
        nmax)

!
!  Allocate various densities
!

      allocate(sigmaover(3,ns))
      allocate(abc0(3,npts),abc1(3,npts),abc2(3,npts), &
        abc3(3,npts),abc4(6,npts))
      allocate(bmm(3,npts),blm(3,npts))
! 
!       oversample densities
!
      do i=1,npts
        abc0(1,i) = sigma(i)
        abc0(2,i) = sigma(npts+i)
        abc0(3,i) = sigma(2*npts+i)
      enddo

      call oversample_fun_surf(3,npatches,norders,ixyzs,iptype,& 
         npts,abc0,novers,ixyzso,ns,sigmaover)


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


!
!  step 1: evaluate D0 (rhom, rhop, mum), and D0'(rhom,rhop,mum)
!  and add the corresponding potentials to their respective blocks
!  also, hold on to the respective blocks in pot_aux_tmp,
!  and grad_aux_tmp
!
!  and S0 (rhom,rhop,mum), S0' (rhom,rhop,mum), and S0'' (rhom, rhop,
!  mum)
!
!  abc1 = D0(rhom,rhop,mum)
!  abc2 = D0'(rhom,rhop,mum) + S0''(rhom,rhop,mum)
!  abc3 = S0'(rhom,rhop,mum)
!  blm = -dzk*\nabla_{\Gamma} S_{0}[\mum] 
!  bmm = dzk*\nabla_{\Gamma}S_{0}[\rhom] &
!     - dzk* (n \times \nabla_{\Gamma} S_{0}[\rhop]
!
      nd = 3
      allocate(dipvec0(nd,3,ns),charges0(nd,ns))
      allocate(pot_aux(nd,npts),grad_aux(nd,3,npts), &
        hess_aux(nd,6,npts))

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
      hess_aux = 0
      abc1 = 0

      call lfmm3d_t_d_g_vec(nd,eps,ns,sources,dipvec0,npts,srctmp, &
        abc1,grad_aux)

      print *, "done with first fmm"

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        abc2(1:3,i) = grad_aux(1:3,1,i)*srcvals(10,i) + &
           grad_aux(1:3,2,i)*srcvals(11,i) + &
           grad_aux(1:3,3,i)*srcvals(12,i)
      enddo
!$OMP END PARALLEL DO

      print *, "finished computing abc2"

      pot_aux = 0
      grad_aux = 0
      hess_aux = 0
      call lfmm3d_t_c_h_vec(nd,eps,ns,sources,charges0,npts,srctmp, &
        pot_aux,grad_aux,hess_aux)

      print *, "done with second fmm"

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(u1,u2,u3,u4)
      do i=1,npts
        abc3(1:3,i) = grad_aux(1:3,1,i)*srcvals(10,i) + &
           grad_aux(1:3,2,i)*srcvals(11,i) + &
           grad_aux(1:3,3,i)*srcvals(12,i)

        abc2(1:3,i) = abc2(1:3,i) + &
                    hess_aux(1:3,1,i)*srcvals(10,i)*srcvals(10,i) + &
                    hess_aux(1:3,2,i)*srcvals(11,i)*srcvals(11,i) + &
                    hess_aux(1:3,3,i)*srcvals(12,i)*srcvals(12,i) + &
                    2*hess_aux(1:3,4,i)*srcvals(11,i)*srcvals(10,i) + &
                    2*hess_aux(1:3,5,i)*srcvals(12,i)*srcvals(10,i) + &
                    2*hess_aux(1:3,6,i)*srcvals(11,i)*srcvals(12,i)

         call dot_prod3d(grad_aux(3,1:3,i),srcvals(4,i),u1)
         call dot_prod3d(grad_aux(3,1:3,i),srcvals(7,i),u2)

         blm(1:3,i) = -(u1*wtmp1(1:3,i) + u2*wtmp2(1:3,i))*dzk
         call dot_prod3d(grad_aux(1,1:3,i),srcvals(4,i),u1)
         call dot_prod3d(grad_aux(1,1:3,i),srcvals(7,i),u2)
         
         call dot_prod3d(grad_aux(2,1:3,i),srcvals(4,i),u3)
         call dot_prod3d(grad_aux(2,1:3,i),srcvals(7,i),u4)

         bmm(1:3,i) = dzk*(u1*wtmp1(1:3,i) + u2*wtmp2(1:3,i) - &
           u3*wtmp3(1:3,i) - u4*wtmp4(1:3,i))
      enddo
!$OMP END PARALLEL DO

      print *, "finished computing bll,blm,abc2,abc3"
!
!  Add in precorrected quadratures for relevant quantities
!

      call get_fmm_thresh(12,ns,srcover,12,npts,srcvals,thresh)
      print *, "thresh=",thresh
!
!       Add near field precomputed contribution
!
      call cpu_time(t1)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart) &
!$OMP PRIVATE(jstart,npols,l,wtmp,uf,vf,rhop,rhom,rmum,w1,w2,w3,w4,w5)
      do i=1,npts
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch) 
          do l=1,npols
            rhom = abc0(1,jstart+l-1)
            rhop = abc0(2,jstart+l-1)
            rmum = abc0(3,jstart+l-1)

            w1 = wnear(nquad+jquadstart+l-1)
            w2 = wnear(2*nquad+jquadstart+l-1)
            w3 = wnear(3*nquad+jquadstart+l-1)
            w4 = wnear(4*nquad+jquadstart+l-1)
            w5 = wnear(5*nquad+jquadstart+l-1)


            abc1(1,i) =  abc1(1,i) + w1*rhom
            abc1(2,i) =  abc1(2,i) + w1*rhop
            abc1(3,i) =  abc1(3,i) + w1*rmum 

            abc2(1,i) =  abc2(1,i) + w5*rhom
            abc2(2,i) =  abc2(2,i) + w5*rhop 
            abc2(3,i) =  abc2(3,i) + w5*rmum 

            wtmp = w2*srcvals(10,i) + w3*srcvals(11,i) + &
              w4*srcvals(12,i)
            uf = w2*srcvals(4,i) + w3*srcvals(5,i) + &
              w4*srcvals(6,i)
            vf = w2*srcvals(7,i) + w3*srcvals(8,i) + &
              w4*srcvals(9,i)

            abc3(1,i) =  abc3(1,i) + wtmp*rhom
            abc3(2,i) =  abc3(2,i) + wtmp*rhop
            abc3(3,i) =  abc3(3,i) + wtmp*rmum

            blm(1:3,i) = blm(1:3,i) - dzk*(uf*wtmp1(1:3,i) + &
              vf*wtmp2(1:3,i))*rmum

            bmm(1:3,i) = bmm(1:3,i) + dzk*(uf*wtmp1(1:3,i) + &
              vf*wtmp2(1:3,i))*rhom
            bmm(1:3,i) = bmm(1:3,i) - dzk*(uf*wtmp3(1:3,i) + &
              vf*wtmp4(1:3,i))*rhop
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

      print *, "finished near computation for abc2,abc3,bmm,blm"

!     Remove near contribution of the FMM
!

      allocate(ctmp0(3,nmax),dtmp0(3,3,nmax))
      allocate(dpottmp(3),dgradtmp(3,3),srctmp2(3,nmax))
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2) &
!$OMP PRIVATE(ctmp0,dtmp0,l,jstart,nss,dpottmp,dgradtmp,u1,u2,u3,u4,rr)
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

            dtmp0(1:3,1:3,nss)=dipvec0(1:3,1:3,l)
            
            rr = sqrt((srcover(1,l)-srcvals(1,i))**2 + &
              (srcover(2,l)-srcvals(2,i))**2 + &
              (srcover(3,l)-srcvals(3,i))**2)
            if(rr.lt.thresh) goto 1311
            call l3d_spp_sum_dp(srcover(1,l),12,srcvals(1,i),0,dpars,&
              0,zpars,0,ipars,u1)
            abc2(1,i) = abc2(1,i) - u1*sigmaover(1,l)*whtsover(l)
            abc2(2,i) = abc2(2,i) - u1*sigmaover(2,l)*whtsover(l)
            abc2(3,i) = abc2(3,i) - u1*sigmaover(3,l)*whtsover(l)
 1311       continue            
          enddo
        enddo

        dpottmp = 0

        call l3ddirectdp(nd,srctmp2,dtmp0,nss,srctmp(1,i), &
          ntarg0,dpottmp,thresh)

        abc1(1:3,i) = abc1(1:3,i) - dpottmp(1:3)

        dpottmp = 0
        dgradtmp = 0

        call l3ddirectcg(nd,srctmp2,ctmp0,nss,srctmp(1,i), &
          ntarg0,dpottmp,dgradtmp,thresh)
        
        abc3(1:3,i) = abc3(1:3,i)-(dgradtmp(1:3,1)*srcvals(10,i) + &
           dgradtmp(1:3,2)*srcvals(11,i) + &
           dgradtmp(1:3,3)*srcvals(12,i))

         call dot_prod3d(dgradtmp(3,1:3),srcvals(4,i),u1)
         call dot_prod3d(dgradtmp(3,1:3),srcvals(7,i),u2)

         blm(1:3,i) =blm(1:3,i)+(u1*wtmp1(1:3,i) + u2*wtmp2(1:3,i))*dzk

         call dot_prod3d(dgradtmp(1,1:3),srcvals(4,i),u1)
         call dot_prod3d(dgradtmp(1,1:3),srcvals(7,i),u2)

         call dot_prod3d(dgradtmp(2,1:3),srcvals(4,i),u3)
         call dot_prod3d(dgradtmp(3,1:3),srcvals(7,i),u4)

         bmm(1:3,i) = bmm(1:3,i) - dzk*(u1*wtmp1(1:3,i) + &
           u2*wtmp2(1:3,i) - u3*wtmp3(1:3,i) - u4*wtmp4(1:3,i))
      enddo
!$OMP END PARALLEL DO     

!
! Completed computation of the following quantities
!  abc1 = D0(rhom,rhop,mum)
!  abc2 = D0'(rhom,rhop,mum) + S0''(rhom,rhop,mum)
!  abc3 = S0'(rhom,rhop,mum)
!  blm = -dzk*\nabla_{\Gamma} S_{0}[\mum] 
!  bmm = dzk*\nabla_{\Gamma}S_{0}[\rhom] &
!     - dzk* (n \times \nabla_{\Gamma} S_{0}[\rhop]
!
!

      print *, "finished all computation for abc1,abc2,abc3,blm,bmm"
!
!  Now compute D0 of abc1 and store in abc0 since it is no
!  longer used
!
      call oversample_fun_surf(nd,npatches,norders,ixyzs,iptype,& 
          npts,abc1,novers,ixyzso,ns,sigmaover)

!$OMP PARALLEL DO DEFAULT(SHARED) 
      do i=1,ns
        dipvec0(1:3,1,i) = srcover(10,i)*sigmaover(1:3,i)*whtsover(i)* &
          over4pi
        dipvec0(1:3,2,i) = srcover(11,i)*sigmaover(1:3,i)*whtsover(i)* &
          over4pi
        dipvec0(1:3,3,i) = srcover(12,i)*sigmaover(1:3,i)*whtsover(i)* &
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
            abc0(1:3,i) = abc0(1:3,i) + w1*abc1(1:3,jstart+l-1)
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

            dtmp0(1:3,1:3,nss)=dipvec0(1:3,1:3,l)
          enddo
        enddo

        dpottmp = 0

        call l3ddirectdp(nd,srctmp2,dtmp0,nss,srctmp(1,i), &
          ntarg0,dpottmp,thresh)

        abc0(1:3,i) = abc0(1:3,i) - dpottmp(1:3)
      enddo
!$OMP END PARALLEL DO      


!
!  End of computing abc0 = D0 (abc1) = D0^2[rhom,rhop,rmum]
!
      deallocate(dipvec0,charges0,pot_aux,grad_aux,hess_aux)
      deallocate(sigmaover,dpottmp,dgradtmp,ctmp0,dtmp0)

      nd = 6
      allocate(charges0(nd,ns),sigmaover(nd,ns))
      allocate(pot_aux(nd,npts),grad_aux(nd,3,npts))
      allocate(dpottmp(nd),dgradtmp(nd,3))
      allocate(ctmp0(nd,nmax))

!
!
!  The 6 densities now being considered are ordered as follows:
!  
!  abc4(1) = ((D0'+S0'') + 2HS0')[\rhom]
!  abc4(2) = ((D0'+S0'') + 2HS0')[\rhop]
!  abc4(3) = ((D0'+S0'') + 2HS0')[\rmum]
!  abc4(4) = q^(-)
!  abc4(5) = q^(+)
!  abc4(6) = r^(-)
!
!  We first organize these densities and oversample them
!
!
      call prin2('abc3=*',abc3,24)
      call prin2('abc2=*',abc2,24)
      call prin2('sigma4=*',sigma(3*npts+1),24)
      call prin2('sigma5=*',sigma(4*npts+1),24)
      call prin2('sigma6=*',sigma(5*npts+1),24)
      do i=1,npts
        abc4(1:3,i) = abc2(1:3,i) + 2*curv(i)*abc3(1:3,i)
        abc4(4,i) = sigma(3*npts+i)
        abc4(5,i) = sigma(4*npts+i)
        abc4(6,i) = sigma(5*npts+i)
      enddo

      call prin2('curv=*',curv,24)
      call prin2('abc4=*',abc4(1:3,1:10),30)
      call prinf('nd=*',nd,1)

      allocate(bbp(3,npts))
      bbp = 0

      sigmaover=  0

      call oversample_fun_surf(nd,npatches,norders,ixyzs,iptype,& 
          npts,abc4,novers,ixyzso,ns,sigmaover)

!
!  Note that we compute the quanitites along with the gradient 
!  since \grad S_{0}[q^{+}] is required in the later representations
!
!$OMP PARALLEL DO DEFAULT(SHARED) 
      do i=1,ns
        charges0(1:6,i) = sigmaover(1:6,i)*whtsover(i)*over4pi
      enddo
!$OMP END PARALLEL DO      

      pot_aux = 0
      grad_aux = 0

      call prinf('nd=*',nd,1)

      call lfmm3d_t_c_g_vec(nd,eps,ns,sources,charges0,npts,srctmp, &
        pot_aux,grad_aux)

!
!
!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        bbp(1:3,i) = grad_aux(5,1:3,i) 
      enddo
!$OMP END PARALLEL DO

!
!  Add near quadrature correction
!

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart) &
!$OMP PRIVATE(jstart,npols,l,w1,w2,w3,w4)
      do i=1,npts
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch) 
          do l=1,npols
            w1 = wnear(jquadstart+l-1)
            w2 = wnear(jquadstart+l-1+2*nquad)
            w3 = wnear(jquadstart+l-1+3*nquad)
            w4 = wnear(jquadstart+l-1+4*nquad)
            pot_aux(1:6,i) = pot_aux(1:6,i) + w1*abc4(1:6,jstart+l-1)
            bbp(1,i) = bbp(1,i) + w2*abc4(5,jstart+l-1)
            bbp(2,i) = bbp(2,i) + w3*abc4(5,jstart+l-1)
            bbp(3,i) = bbp(3,i) + w4*abc4(5,jstart+l-1)
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

            ctmp0(1:6,nss)=charges0(1:6,l)
          enddo
        enddo

        dpottmp = 0
        dgradtmp = 0

        call l3ddirectcg(nd,srctmp2,ctmp0,nss,srctmp(1,i), &
          ntarg0,dpottmp,dgradtmp,thresh)

        pot_aux(1:6,i) = pot_aux(1:6,i) - dpottmp(1:6)
        bbp(1:3,i) = bbp(1:3,i) - dgradtmp(5,1:3)
      enddo
!$OMP END PARALLEL DO      


!
! At this stage we are ready to compile the first three components
! of the output potential
!
      do j=1,3
        do i=1,npts
          pot(i+(j-1)*npts) = -4*(abc0(j,i)-pot_aux(j,i)-pot_aux(j+3,i))
        enddo
      enddo

!
!  At this stage the following quantities have already been
!  computed
!  bpp 
!  blm
!  bmm
!
!  Now we proceed to compute bjm and bbm
!  given by the formulae
!
!  bjm = -dzk*S_{ik} [blm] - \nabla S_{ik} [r^{-}] - \nabla \times S_{ik} [bmm]
!  bbm = -dzk*S_{ik} [bmm] + \nabla S_{ik} [q^{-}] + \nabla \times S_{ik} [bll^{-}]
!
      
      nd = 6
      deallocate(charges0,sigmaover)
      deallocate(pot_aux,grad_aux)
      deallocate(dpottmp,dgradtmp)
      deallocate(ctmp0,abc0)

      nd = 8
      allocate(zcharges0(nd,ns),sigmaover(nd,ns),abc0(nd,npts))

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        abc0(1:3,i) = blm(1:3,i)
        abc0(4:6,i) = bmm(1:3,i)
        abc0(7,i) = sigma(3*npts+i)
        abc0(8,i) = sigma(5*npts+i)
      enddo
!$OMP END PARALLEL DO

      call oversample_fun_surf(nd,npatches,norders,ixyzs,iptype,& 
          npts,abc0,novers,ixyzso,ns,sigmaover)
        
!
!$OMP PARALLEL DO DEFAULT(SHARED) 
      do i=1,ns
        zcharges0(1:8,i) = sigmaover(1:8,i)*whtsover(i)*over4pi
      enddo
!$OMP END PARALLEL DO      

      allocate(zpot_aux(nd,npts),zgrad_aux(nd,3,npts))
      allocate(pot_aux(nd,npts),grad_aux(nd,3,npts))
      zpot_aux = 0
      zgrad_aux = 0

      call prinf('nd=*',nd,1)

      zk0 = dzk*ima
      call prin2('zk0=*',zk0,2)

      print *, "before fmm"

      call hfmm3d_t_c_g_vec(nd,eps,zk0,ns,sources,zcharges0,npts,srctmp, &
        zpot_aux,zgrad_aux)

      print *, "after fmm"


      pot_aux = real(zpot_aux)
      grad_aux = real(zgrad_aux)


!
!  Add near quadrature correction
!

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart) &
!$OMP PRIVATE(jstart,npols,l,w1,w2,w3,w4)
      do i=1,npts
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch) 
          do l=1,npols
            w1 = wnear(jquadstart+l-1+6*nquad)
            w2 = wnear(jquadstart+l-1+7*nquad)
            w3 = wnear(jquadstart+l-1+8*nquad)
            w4 = wnear(jquadstart+l-1+9*nquad)
            pot_aux(1:8,i) = pot_aux(1:8,i) + w1*abc0(1:8,jstart+l-1)
            grad_aux(1:8,1,i) = grad_aux(1:8,1,i) + w2*abc0(1:8,jstart+l-1)
            grad_aux(1:8,2,i) = grad_aux(1:8,2,i) + w3*abc0(1:8,jstart+l-1)
            grad_aux(1:8,3,i) = grad_aux(1:8,3,i) + w4*abc0(1:8,jstart+l-1)
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO     
      
      print *, "after helmholtz near correction"
      print *, "nmax=",nmax
      print *, "nd=",nd


!
! Subtract near contributions computed via fmm
!
      allocate(zpottmp(nd),zgradtmp(nd,3))
      allocate(zctmp0(nd,nmax))
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2) &
!$OMP PRIVATE(zctmp0,l,jstart,nss,zpottmp,zgradtmp)
      do i=1,npts
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          do l=ixyzso(jpatch),ixyzso(jpatch+1)-1
            nss = nss+1
            srctmp2(1,nss) = srcover(1,l)
            srctmp2(2,nss) = srcover(2,l)
            srctmp2(3,nss) = srcover(3,l)

            zctmp0(1:8,nss)=zcharges0(1:8,l)
          enddo
        enddo

        zpottmp = 0
        zgradtmp = 0

        call h3ddirectcg(nd,zk0,srctmp2,zctmp0,nss,srctmp(1,i), &
          ntarg0,zpottmp,zgradtmp,thresh)

        pot_aux(1:8,i) = pot_aux(1:8,i) - real(zpottmp(1:8))
        grad_aux(1:8,1:3,i) = grad_aux(1:8,1:3,i) - real(zgradtmp(1:8,1:3))
      enddo
!$OMP END PARALLEL DO      

      print *, "finished pot eval"
      deallocate(zpottmp,zgradtmp,zpot_aux,zgrad_aux)

      allocate(bjm(3,npts),bbm(3,npts))
!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        bjm(1,i) = -dzk*pot_aux(1,i) & 
            - grad_aux(8,1,i) 
!            - (grad_aux(6,2,i)-grad_aux(5,3,i))
        bjm(2,i) = -dzk*pot_aux(2,i) & 
            - grad_aux(8,2,i) 
!            - (grad_aux(4,3,i)-grad_aux(6,1,i))
        bjm(3,i) = -dzk*pot_aux(3,i) & 
           - grad_aux(8,3,i) 
!           - (grad_aux(5,1,i)-grad_aux(4,2,i))

        bbm(1,i) = -dzk*pot_aux(4,i) + grad_aux(7,1,i) + &
           (grad_aux(3,2,i)-grad_aux(2,3,i))
        bbm(2,i) = -dzk*pot_aux(5,i) + grad_aux(7,2,i) + &
           (grad_aux(1,3,i)-grad_aux(3,1,i))
        bbm(3,i) = -dzk*pot_aux(6,i) + grad_aux(7,3,i) + &
           (grad_aux(2,1,i)-grad_aux(1,2,i))
      enddo
!$OMP END PARALLEL DO

       print *, "done computing bjm,bbm"

!
! End of computing bj-, bb- and bb^+
!
! Only thing left to compute now is \nabla \cdot S_{0} [n \times n \times bb^-]
!

      rr = 0

      allocate(wts(npts))
      call get_qwts(npatches,norders,ixyzs,iptype,npts,srcvals,wts)

      do i=1,npts
        rr = rr + sigma(i)**2*wts(i)
      enddo

!
!  Now test accuracy of bj- \cdot n
!
      njh = 5
      ifder = 1
      rscale = 1.0d0
      call prin2('zk0=*',zk0,2)
      call besseljs3d(njh,zk0,rscale,jvals,ifder,fjder)
      call prin2('jvals=*',jvals,2*njh+2)

      call h3dall(njh,zk0,rscale,hvals,ifder,fhder)
      call prin2('hvals=*',hvals,2*njh+2)

      n = 2

      ztmp = ima*(jvals(3)*hvals(2) - jvals(4)*hvals(3))
      call prin2('ztmp=*',ztmp,2)
      rtmp = real(ztmp)
      rfac = rtmp*dzk**2*n*(n+1.0d0)/(2*n+1.0d0)

      ztmp = -ima*dzk**2*(jvals(3)*fhder(3)+fjder(3)*hvals(3))/2
      call prin2('ztmp=*',ztmp,2)
      rtmp = real(ztmp)
      rfac = rfac +1.76d0*rtmp 
      
      errjn = 0
      rrjn = 0
      do i=1,npts
        rjn = bjm(1,i)*srcvals(10,i) + bjm(2,i)*srcvals(11,i)+&
           bjm(3,i)*srcvals(12,i)
        errjn = errjn + abs(rjn-rfac*sigma(i))**2*wts(i)
        rrjn = rrjn + abs(rjn)**2*wts(i)
      enddo

      errjn = sqrt(errjn/rrjn)
      call prin2('error in normal component of j=*',errjn,1)

      stop

      

      return
      end subroutine lpcomp_statj_gendeb_addsub
!
!
!
!
!
