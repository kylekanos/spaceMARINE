module write_plotfile_module

  use ml_layout_module
  use multifab_module
  use fabio_module

  implicit none

contains
  
  subroutine write_plotfile(mla,phi,istep,dx,nvars,time,prob_lo,prob_hi)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: phi(:)
    integer        , intent(in   ) :: istep,nvars
    real(dp_t)     , intent(in   ) :: dx(:),time
    real(dp_t)     , intent(in   ) :: prob_lo(mla%dim), prob_hi(mla%dim)

    ! local variables
    character(len=8)  :: plotfile_name
    character(len=20) :: variable_names(nvars)

    type(multifab), allocatable :: plotdata(:)

    ! dimensioned as an array with size dm for fabio_ml_multifab_write_d
    real(dp_t) :: dx_vec(mla%dim)

    integer n, nlevs, nc

    nlevs = mla%nlevel
    nc = ncomp(phi(1))
    dx_vec(:) = dx(1)

    variable_names(1) = "rho"
    variable_names(2) = "vx"
    variable_names(3) = "vy"
    variable_names(4) = "vz"
    variable_names(5) = "E"
    if(nvars > 5) then
       do n=5,nvars
          write(variable_names(n),'(a,i2.2)') "var",n
       enddo
    endif

    allocate(plotdata(nlevs))
    
    
    
    do n=1,nlevs
       ! build plotdata with 1 component and 0 ghost cells
       call multifab_build(plotdata(n),mla%la(n),nc,0)

       ! copy the state into plotdata
       call multifab_copy_c(plotdata(n),1,phi(n),1,nc,0)
    end do

    ! define the name of the plotfile that will be written
    write(unit=plotfile_name,fmt='("plt",i5.5)') istep
    if ( parallel_IOProcessor() ) print*, 'Writing ' // plotfile_name // ' ...'

    ! write the plotfile
    call fabio_ml_multifab_write_d(plotdata, mla%mba%rr(:,1), plotfile_name, variable_names, &
                                   mla%mba%pd(1), prob_lo, prob_hi, time, dx_vec)

    ! make sure to destroy the multifab or you'll leak memory
    do n=1,nlevs
       call multifab_destroy(plotdata(n))
    end do

  end subroutine write_plotfile

end module write_plotfile_module
