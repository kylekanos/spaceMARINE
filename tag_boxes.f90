module tag_boxes_module

  use multifab_module
  use bl_error_module

  implicit none 

contains

  subroutine tag_boxes(mf,tagboxes,dx,lev,aux_tag_mf)

    type( multifab)         , intent(in   ) :: mf
    type(lmultifab)         , intent(inout) :: tagboxes
    real(dp_t)              , intent(in   ) :: dx
    integer                 , intent(in   ) :: lev
    type(multifab), optional, intent(in   ) :: aux_tag_mf
    ! aux_tag_mf allows user to pass in additional multifabs for tagging logic
    
    ! local variables
    real(kind = dp_t), pointer :: mfp(:,:,:,:)
    logical          , pointer :: tp(:,:,:,:)
    integer           :: i, lo(get_dim(mf)), hi(get_dim(mf)), ng

    if (present(aux_tag_mf)) then
       call bl_error("tag_boxes.f90: aux_tag_mf passed to tag_boxes without implementation")
    end if

    ng = nghost(mf)

    do i = 1, nfabs(mf)
       mfp => dataptr(mf, i)
       tp  => dataptr(tagboxes, i)
       lo =  lwb(get_box(tagboxes, i))
       hi =  upb(get_box(tagboxes, i))
       select case (get_dim(mf))
         case (2)
            call tag_boxes_2d(tp(:,:,1,1),mfp(:,:,1,:),lo,hi,ng,dx,lev)
         case  (3)
            call tag_boxes_3d(tp(:,:,:,1),mfp(:,:,:,:),lo,hi,ng,dx,lev)
       end select
    end do

  end subroutine tag_boxes

  !> @brief flag the boxes for refinement
  !! uses method of Jiang's MAP code
  subroutine tag_boxes_2d(tagbox,mf,lo,hi,ng,dx,lev)
    integer          , intent(in   ) :: lo(:),hi(:),ng
    logical          , intent(  out) :: tagbox(lo(1)   :,lo(2)   :)
    real(kind = dp_t), intent(in   ) ::     mf(lo(1)-ng:,lo(2)-ng:, 1:)
    real(dp_t)       , intent(in   ) :: dx
    integer          , intent(in   ) :: lev

    ! local variables
    real(dp_t), dimension(ubound(mf,dim=3)) :: dva2, d2va2, dva, thr_va0
    real(dp_t) :: num_tmp, f_weight, s_weight, thr_fl, eps,lv_exp
    integer :: i,j,n,nv

    ! initially say that we do not want to tag any cells for refinement
    tagbox = .false.
    nv = ubound(dva2,dim=1)
    
    ! simple variables...
    num_tmp = sqrt(2d0)*0.5
    f_weight = 0.8d0                   ! how much we depend on 1st derivatives
    s_weight = 1d0 - f_weight          ! how much we depend on 2nd derivatives
    thr_fl = 0.05d0
    eps = 1d-12
    lv_exp  = 2d0**(f_weight*(lev - 1))
    thr_va0 = [(0.1d0, i=1,nv)]
    
    
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
       ! compute first derivative
          do n=1,nv
             dva2(n) = num_tmp*sqrt((mf(i+1,j,n) - mf(i-1,j,n))**2 + (mf(i,j+1,n) - mf(i,j-1,n))**2)
          enddo !- n
       ! compute 2nd derivative
          do n=1,nv
             d2va2(n) = sqrt((mf(i+1,j,n) - 2d0*mf(i,j,n) + mf(i-1,j,n))**2 +&
                             (mf(i,j+1,n) - 2d0*mf(i,j,n) + mf(i,j-1,n))**2)
          enddo !- n
       ! compute ratio of derivatives
          do n=1,nv
             dva(n) = (f_weight * dva2(n) / (abs(mf(i,j,n)) + eps) + &
                       s_weight * d2va2(n) / (dva2(n) + thr_fl * (abs(mf(i,j,n)) + eps))) * lv_exp
          enddo !- n
       ! flag region
          tagbox(i,j) = any(dva > thr_va0)
       enddo !- i
    enddo !- j
    
  end subroutine tag_boxes_2d

  !> @brief flag the boxes for refinement
  !! uses method of Jiang's MAP code
  subroutine tag_boxes_3d(tagbox,mf,lo,hi,ng,dx,lev)

    integer          , intent(in   ) :: lo(:),hi(:),ng
    logical          , intent(  out) :: tagbox(lo(1)   :,lo(2)   :,lo(3)   :)
    real(kind = dp_t), intent(in   ) ::     mf(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:, 1:)
    real(dp_t)       , intent(in   ) :: dx
    integer          , intent(in   ) :: lev

    ! local variables
    real(dp_t), dimension(ubound(mf,dim=4)) :: dva2, d2va2, dva, thr_va0
    real(dp_t) :: num_tmp, f_weight, s_weight, thr_fl, eps, lv_exp
    integer :: i,j,k,n,nv

    ! initially say that we do not want to tag any cells for refinement
    tagbox = .false.
    nv = ubound(dva2,dim=1)
    
    ! simple variables...
    num_tmp = sqrt(2d0)*0.5
    f_weight = 0.8d0
    s_weight = 1d0 - f_weight
    thr_fl = 0.05d0
    eps = 1d-12
    lv_exp  = 2d0**(f_weight*(lev - 1))
    thr_va0 = [(0.1d0, i=1,nv)]
    
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
          ! compute first derivative
             do n=1,nv
                dva2(n) = num_tmp*sqrt((mf(i+1,j,k,n) - mf(i-1,j,k,n))**2 + &
                                       (mf(i,j+1,k,n) - mf(i,j-1,k,n))**2 + &
                                       (mf(i,j,k+1,n) - mf(i,j,k-1,n))**2)
             enddo !- n
          ! compute 2nd derivative
             do n=1,nv
                d2va2(n) = sqrt((mf(i+1,j,k,n) - 2d0*mf(i,j,k,n) + mf(i-1,j,k,n))**2 +&
                                (mf(i,j+1,k,n) - 2d0*mf(i,j,k,n) + mf(i,j-1,k,n))**2 +&
                                (mf(i,j,k+1,n) - 2d0*mf(i,j,k,n) + mf(i,j,k-1,n))**2)
             enddo !- n
          ! compute ratio of derivatives
             do n=1,nv
                dva(1) = (f_weight * dva2(1) / (abs(mf(i,j,k,n)) + eps) + &
                         s_weight * d2va2(1) / (dva2(1) + thr_fl * (abs(mf(i,j,k,n)) + eps))) * lv_exp
             enddo !- n
          ! flag region
             tagbox(i,j,k) = any(dva > thr_va0)
          enddo !- i
       enddo !- j
    enddo !- k

  end subroutine tag_boxes_3d

end module tag_boxes_module
