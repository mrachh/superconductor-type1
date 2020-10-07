!
!  This file contains the following user callable routines:
!    
!
!
!
subroutine get_second_fundamental_form(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,sfform)
!
!------------------------
!  This subroutine computes the second fundamental form at
!  the discretization nodes on the surface.
!
!  The second fundamental form is
!  
!  .. math::
!    
!    \begin{bmatrix} x_{uu} \cdot n & x_{uv} \cdot n \\
!    x_{uv} \cdot n & x_{vv} \cdot n \end{bmatrix}
!
!  Input arguments:
!
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization of each patch
!    - ixyzs: integer(npatches+1)
!        starting location of points on patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangle discretized using RV nodes
!    - npts: integer
!        total number of points on the surface
!    - srccoefs: double precision (9,npts)
!        koornwinder expansion coefs of geometry info
!    - srcvals: double precision (12,npts)
!        xyz, dxyz/du,dxyz/dv, normals at all nodes
!
!  Output arguments:
!
!    - sfform: double precision(2,2,npts)
!        second fundamental form at the discretization nodes
!--------------------------
!
  
  implicit none
  integer, intent(in) :: npatches,norders(npatches)
  integer, intent(in) :: ixyzs(npatches+1),iptype(npatches)
  integer, intent(in) :: npts
  real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts)
  real *8, intent(out) :: sfform(2,2,npts)
  integer i
  real *8 L,M,N
  real *8, allocatable :: dxuv(:,:)
  real *8, allocatable :: dxuv2(:,:,:)


  allocate(dxuv(6,npts))
! Calculating x_{uu}, x_{uv}, x_{uv} stored in dxuv2 
 
  do i=1,npts
    dxuv(1,i) = srcvals(4,i)  
    dxuv(2,i) = srcvals(5,i)  
    dxuv(3,i) = srcvals(6,i)  
    dxuv(4,i) = srcvals(7,i)  
    dxuv(5,i) = srcvals(8,i)  
    dxuv(6,i) = srcvals(9,i)   
  enddo

  allocate(dxuv2(6,2,npts))

  call get_surf_uv_grad(6,npatches,norders,ixyzs,iptype,npts,dxuv,dxuv2)



  do i=1,npts
    L = dxuv2(1,1,i)*srcvals(10,i) + dxuv2(2,1,i)*srcvals(11,i) + dxuv2(3,1,i)*srcvals(12,i) ! Calculation of L, M, N. L = x_uu \cdot n 
    M = dxuv2(1,2,i)*srcvals(10,i) + dxuv2(2,2,i)*srcvals(11,i) + dxuv2(3,2,i)*srcvals(12,i)  ! M = x_uv \cdot n
    N = dxuv2(4,2,i)*srcvals(10,i) + dxuv2(5,2,i)*srcvals(11,i) + dxuv2(6,2,i)*srcvals(12,i)  ! N = x_vv \cdot n
    sfform(1,1,i) = L
    sfform(2,1,i) = M
    sfform(1,2,i) = M
    sfform(2,2,i) = N
  enddo

  return
end subroutine get_second_fundamental_form
!
!
!
subroutine get_mean_curvature(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,mean_curv)
!
!------------------------
!  This subroutine computes the mean curvature at
!  the discretization nodes on the surface.
!
!  The mean curvature is
!  
!  .. math::
!    
!    0.5*Trace(II \cdot I^{-1}) \\
!    
!
!  Input arguments:
!
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization of each patch
!    - ixyzs: integer(npatches+1)
!        starting location of points on patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangle discretized using RV nodes
!    - npts: integer
!        total number of points on the surface
!    - srccoefs: double precision (9,npts)
!        koornwinder expansion coefs of geometry info
!    - srcvals: double precision (12,npts)
!        xyz, dxyz/du,dxyz/dv, normals at all nodes
!
!  Output arguments:
!
!    - mean_curv: double precision(npts)
!        mean curvature at the discretization nodes
!--------------------------
!
  
  implicit none
  integer, intent(in) :: npatches,norders(npatches)
  integer, intent(in) :: ixyzs(npatches+1),iptype(npatches)
  integer, intent(in) :: npts
  real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts)
  real *8, intent(out) :: mean_curv(npts)
  integer i
  real *8, allocatable :: ffform(:,:,:),ffforminv(:,:,:)
  real *8, allocatable :: sfform(:,:,:)

  allocate(ffform(2,2,npts))

  call get_first_fundamental_form(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,ffform)



  allocate(ffforminv(2,2,npts))

  call get_inv_first_fundamental_form(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,ffforminv)

  allocate(sfform(2,2,npts))

  call get_second_fundamental_form(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,sfform)
  
!  print *,"Point on surface:", srcvals(1,3),srcvals(2,3), srcvals(3,3) 
!  print *,"First fundamental form=", ffform(:,:, 3) 
!  print *,"Inverse first fundamental form=", ffforminv(:,:, 3)
!  print *,"Second fundamental form=", sfform(:,:, 3)
 


! Calculating mean curvature 
 
  do i=1,npts
    mean_curv(i) = -0.5*(sfform(1,1,i)*ffforminv(1,1,i) + &
                     sfform(1,2,i)*ffforminv(2,1,i) + &
                     sfform(2,1,i)*ffforminv(1,2,i) + &
                     sfform(2,2,i)*ffforminv(2,2,i))
  enddo
!  print *,"Mean=", mean_curv(3)
 
  return
end subroutine get_mean_curvature

!
!
!
!
!
!
subroutine surf_grad(npatches,norders,ixyzs,iptype,npts, &
  srccoefs,srcvals,fin,gradf)
!
!-----------------------------
!  Compute the surface gradient of scalar function fin  
!
!  Input arguments:
!
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization of each patch
!    - ixyzs: integer(npatches+1)
!        starting location of points on patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangle discretized using RV nodes
!    - npts: integer
!        total number of points on the surface
!    - srccoefs: double precision (9,npts)
!        koornwinder expansion coefs of geometry info
!    - srcvals: double precision (12,npts)
!        xyz, dxyz/du,dxyz/dv, normals at all nodes
!
!    - fin: double precision (npts)
!         vector function on surface
!  Output arguments:
!
!    - gradf: double precision(3,npts)
!        surface gradient 
!        
!-----------------------------
!
!

  implicit none
  integer, intent(in) :: npatches,norders(npatches)
  integer, intent(in) :: ixyzs(npatches+1),iptype(npatches)
  integer, intent(in) :: npts
  real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),fin(npts)
  real *8, intent(out) :: gradf(3,npts)
  real *8, allocatable :: ffforminv(:,:,:)
  real *8, allocatable :: dfuv(:,:,:)
  real *8 E,F1,F2,G,W_sq,a,b,fdu,fdv
  integer i,istart,npols,j,l

  allocate(ffforminv(2,2,npts))



  call get_inv_first_fundamental_form(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,ffforminv)

  allocate(dfuv(1,2,npts))
! Calculating f_{u}, f_{v}} stored in dfuv 
 
  call get_surf_uv_grad(1,npatches,norders,ixyzs,iptype,npts,fin,dfuv)



  do i=1,npts
    E = ffforminv(1,1,i)
    F1 = ffforminv(1,2,i)
    F2 = ffforminv(2,1,i)
    G = ffforminv(2,2,i)
    fdu = dfuv(1,1,i)
    fdv = dfuv(1,2,i)
    a = E*fdu+F1*fdv
    b = F2*fdu+G*fdv

    gradf(1,i) = a*srcvals(4,i) + b*srcvals(7,i)
    gradf(2,i) = a*srcvals(5,i) + b*srcvals(8,i)
    gradf(3,i) = a*srcvals(6,i) + b*srcvals(9,i)


                 
  enddo 

  return
end subroutine surf_grad
!




subroutine surf_grad2(npatches,norders,ixyzs,iptype,npts, &
  srccoefs,srcvals,fin,gradf)
!
!-----------------------------
!  Compute the surface gradient of scalar function fin  
!
!  Input arguments:
!
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization of each patch
!    - ixyzs: integer(npatches+1)
!        starting location of points on patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangle discretized using RV nodes
!    - npts: integer
!        total number of points on the surface
!    - srccoefs: double precision (9,npts)
!        koornwinder expansion coefs of geometry info
!    - srcvals: double precision (12,npts)
!        xyz, dxyz/du,dxyz/dv, normals at all nodes
!
!    - fin: double precision (npts)
!         vector function on surface
!  Output arguments:
!
!    - gradf: double precision(3,npts)
!        surface gradient 
!        
!-----------------------------
!
!

  implicit none
  integer, intent(in) :: npatches,norders(npatches)
  integer, intent(in) :: ixyzs(npatches+1),iptype(npatches)
  integer, intent(in) :: npts
  real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),fin(3,npts)
  real *8, intent(out) :: gradf(3,npts)
  real *8, allocatable :: ffform(:,:,:)
  real *8, allocatable :: dfuv(:,:,:)
  real *8 E,F,G,W_sq,a(3),b(3),fdu,fdv
  integer i,istart,npols,j,l

  allocate(ffform(2,2,npts))



  call get_first_fundamental_form(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,ffform)

  allocate(dfuv(1,2,npts))
! Calculating f_{u}, f_{v}} stored in dfuv 
 
  call get_surf_uv_grad(1,npatches,norders,ixyzs,iptype,npts,fin,dfuv)



  do i=1,npts
    E = ffform(1,1,i)
    F = ffform(1,2,i)
    G = ffform(2,2,i)
    W_sq = E*G - F**2
    fdu = dfuv(1,1,i)
    fdv = dfuv(1,2,i)
    a(1) = (G*srcvals(4,i)-F*srcvals(7,i))/W_sq 
    a(2) = (G*srcvals(5,i)-F*srcvals(8,i))/W_sq 
    a(3) = (G*srcvals(6,i)-F*srcvals(9,i))/W_sq 
    b(1) = (E*srcvals(7,i)-F*srcvals(4,i))/W_sq 
    b(2) = (E*srcvals(8,i)-F*srcvals(5,i))/W_sq 
    b(3) = (E*srcvals(9,i)-F*srcvals(6,i))/W_sq 


    gradf(1,i) = a(1)*fdu+b(1)*fdv
    gradf(2,i) = a(2)*fdu+b(2)*fdv
    gradf(3,i) = a(3)*fdu+b(3)*fdv



                 
  enddo 

  return
end subroutine surf_grad2
!

!
!
!
!
!
subroutine surf_div(npatches,norders,ixyzs,iptype,npts, &
  srccoefs,srcvals,fin,divf)
!
!-----------------------------
!  Compute the surface divergence of vector function fin  
!
!  Input arguments:
!
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization of each patch
!    - ixyzs: integer(npatches+1)
!        starting location of points on patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangle discretized using RV nodes
!    - npts: integer
!        total number of points on the surface
!    - srccoefs: double precision (9,npts)
!        koornwinder expansion coefs of geometry info
!    - srcvals: double precision (12,npts)
!        xyz, dxyz/du,dxyz/dv, normals at all nodes
!
!    - fin: double precision (3,npts)
!         vector function on surface
!  Output arguments:
!
!    - divf: double precision(npts)
!        surface divergence 
!        
!-----------------------------
!
!

  implicit none
  integer, intent(in) :: npatches,norders(npatches)
  integer, intent(in) :: ixyzs(npatches+1),iptype(npatches)
  integer, intent(in) :: npts
  real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts),fin(3,npts)
  real *8, intent(out) :: divf(npts)
  real *8, allocatable :: ffform(:,:,:)
  real *8, allocatable :: dfuv(:,:,:)
  real *8 E,F,G,W_sq,a(3),b(3),fdu(3),fdv(3)
  integer i,istart,npols,j,l

  allocate(ffform(2,2,npts))



  call get_first_fundamental_form(npatches,norders,ixyzs,iptype, &
  npts,srccoefs,srcvals,ffform)

  allocate(dfuv(3,2,npts))
! Calculating f_{u}, f_{v}} stored in dfuv 
 
  call get_surf_uv_grad(3,npatches,norders,ixyzs,iptype,npts,fin,dfuv)



  do i=1,npts
    E = ffform(1,1,i)
    F = ffform(1,2,i)
    G = ffform(2,2,i)
    W_sq = E*G - F**2
    fdu(1) = dfuv(1,1,i)
    fdu(2) = dfuv(2,1,i)
    fdu(3) = dfuv(3,1,i)
    fdv(1) = dfuv(1,2,i)
    fdv(2) = dfuv(2,2,i)
    fdv(3) = dfuv(3,2,i)
    a(1) = (G*fdu(1)-F*fdv(1))/W_sq 
    a(2) = (G*fdu(2)-F*fdv(2))/W_sq 
    a(3) = (G*fdu(3)-F*fdv(3))/W_sq 
    b(1) = (E*fdv(1)-F*fdu(1))/W_sq 
    b(2) = (E*fdv(2)-F*fdu(2))/W_sq 
    b(3) = (E*fdv(3)-F*fdu(3))/W_sq 


    divf(i) = a(1)*srcvals(4,i)+a(2)*srcvals(5,i)+a(3)*srcvals(6,i)+ &
              b(1)*srcvals(7,i)+b(2)*srcvals(8,i)+b(3)*srcvals(9,i)
                 
  enddo 

  return
end subroutine surf_div
!
!
!
!
!
