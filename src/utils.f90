! utils.f90
!
! Utilitary functions for ZC sequence in OTFS for sensing
!
! luizfelipe.coelho@smt.ufrj.br / luiz-felipe.da-silveira-coelho@lecnam.net
! Mar 27, 2025

module utils
    use iso_fortran_env, only: sp=>real32
    implicit none
    real(sp), parameter :: pi = 4 * atan(1.0_sp)

contains

! DFT Matrix
function dft(N) result(dft_mat)
    integer, intent(in) :: N
    integer :: nn, kk
    real(sp)  :: im
    complex(sp), dimension(N, N) :: dft_mat

    do nn = 1, N
        do kk = 1, N
            im = 2*pi*(nn-1)*(kk-1)/N
            dft_mat(nn, kk) = exp(cmplx(0, -im))
        end do
    end do

end function dft

! ZC Sequence
function zc(K, Q) result(zc_seq)
    integer, intent(in) :: K, Q
    integer :: i
    integer, dimension(K) :: K_axis
    complex(sp), dimension(K) :: zc_seq

    K_axis = [(i, i = 1, K)]
    if (mod(K, 2) .eq. 0) then
        zc_seq = exp(cmplx(0, Q*pi*K_axis*K_axis / K))
    else
        zc_sec = exp(cmplx(0, Q*pi*K_axis*(K_axis+1) / K))
    end if

end function zc

! Chirp Sequence
function chirp(T, Ts, BW, chirp_factor, K) result(chirp_seq)
    integer, intent(in) :: K
    real(sp), intent(in) :: T, Ts, BW, chirp_factor
    integer :: i
    real(sp) :: Tch
    real(sp), allocatable :: time_axis(:)
    complex(sp), allocatable :: chirp_pulse(:)
    complex(sp), dimension(K) :: chirp_seq

    Tch = T*chirp_factor

    time_axis = [(i, i = 0, Tch, Ts)]
    chirp_pulse = exp(cmplx(0, -pi*))

end function chirp

! Barker Sequence
function barker(K) result(barker_seq)
    integer, intent(in) :: K
    integer, dimension(K) :: barker_seq
    select case (K)
    case (2)
        barker_seq = [1, -1]
    case (3)
        barker_seq = [1, 1, -1]
    case (4)
        barker_seq = [1, -1, 1, 1]
    case (5)
        barker_seq = [1, 1, 1, -1, 1]
    case (7)
        barker_seq = [1, 1, 1, -1, -1, 1, -1]
    case (11)
        barker_seq = [1, 1, 1, -1, -1, -1, 1, -1, -1, 1, -1]
    case (13)
        barker_seq = [1, 1, 1, 1, 1, -1, -1, 1, 1, -1, 1, -1, 1]
    end select
    
end function barker

end module utils