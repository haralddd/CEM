program SolveReyleighEquation
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f03'

  ! Declare global parameters
  ! N_max is the max power and end of sum for series expansion in the I integral
  integer, parameter                                :: N_x = 400, N_q = 200, N_max = 10
  real(C_DOUBLE)                                    :: WaveLength, c, omega, pi = 3.14159266, theta_0, L, deltax
  real, dimension(N_q+1)                            :: pivot
  real(C_DOUBLE), dimension(N_q+1)                  :: realq, theta_array, NormRDiffTheta

  ! Declare and allocate global arrays and such
  integer                                           :: iterator, mIterator, nIterator, k_index, rc
  integer, dimension(N_q+1)                         :: n_values, m_values
  integer, dimension(N_max)                         :: FactorialArray
  
  complex(C_DOUBLE_COMPLEX)                         :: deltaq, deltap, deltak, q_value, p_value, mu, epsilon, alpha_value,&
   alpha_0_value, kappa_value, k_value, Q_big
  complex(C_DOUBLE_COMPLEX), dimension(N_q+1)       :: R_Vector, M_Minus_Array, q_values_vector, p_values_vector
  complex(C_DOUBLE_COMPLEX), dimension(N_x)         :: Surface, SurfacePower, Q_Vector
  complex(C_DOUBLE_COMPLEX), dimension(N_q+1,N_q+1) :: M_Plus_Array
  complex(C_DOUBLE_COMPLEX), dimension(N_x,N_max)   :: Q_ArrayShifted

  type(C_PTR)                                       :: plan

  ! ----------------------------------------------

  ! Define basic parameters
  WaveLength = 600E-9
  c          = 3.0E8
  omega      = 2*pi*c/WaveLength
  Q_big      = 4*(omega/c)
  
  ! Define material properties and polarization
  epsilon     = 0.001
  mu          = 0.001
  kappa_value = kappa(0,epsilon,mu)

  ! Define deltas
  deltaq = Q_big/N_q
  deltap = Q_big/N_q
  deltak = deltap
  deltax = pi/Q_big%re

  ! Define k-related values
  k_index = 10
  k_value = k_index*deltak
  theta_0 = asin(k_value%re*c/omega)

  ! Define length of sample
  L = deltax*N_x

  print*, "-----------"
  print*, "Variables displayed: "
  print*, "omega:               ", omega
  print*, "Q:      ", Q_big
  print*, "N_q:          ", N_q
  print*, "N_x:           ", N_x
  print*, "N_max:       ", N_max
  print*, "deltaq: ", deltaq
  print*, "deltap: ", deltap
  print*, "deltax:             ", deltax
  print*, "L:                  ", L
  print*, "Incident angle:     ", theta_0*180/pi
  print*, "k_value:", k_value
  print*, "-----------"

  ! Define the q and p values and m and n indexes
  do iterator = 1, N_q+1
    n_values(iterator)        = iterator-N_q/2-1
    m_values(iterator)        = iterator-N_q/2-1
    p_values_vector(iterator) = (iterator-N_q/2-1)*deltap
    q_values_vector(iterator) = (iterator-N_q/2-1)*deltaq
  end do


  ! Calculate surface
  call CreateSurface(Surface,N_x)

  ! Calculate factorials for later use
  call FactorialArrayCalculation(FactorialArray,N_max)
  

  ! Loop here?
  ! Split using OpenMP ---------------------------


  ! Calculate FFT of surface
  do iterator = 1, N_max
    !Taking the surface to the power of iterator
    SurfacePower = Surface**iterator
    ! Making an FFTW plan
    plan = fftw_plan_dft_1d(N_x, SurfacePower, Q_Vector, FFTW_FORWARD, FFTW_ESTIMATE)
    ! Executing the plan
    call fftw_execute_dft(plan, SurfacePower, Q_Vector)
    ! Using the fftshift to move the elements around
    call fftshift(Q_Vector,N_x)
    Q_ArrayShifted(:,iterator) = Q_Vector
    ! Destroying plan
    call fftw_destroy_plan(plan)
  end do


  ! Start loop over M⁺
  ! First loop over p values
  do mIterator = 1, N_q+1
    ! Then loop over q values
    do nIterator = 1, N_q+1
      ! Assign current p and q value
      q_value = q_values_vector(nIterator)
      p_value = p_values_vector(mIterator)
      ! Calculating alpha(p) and alpha_0(q)
      alpha_value = alpha(p_value,epsilon,mu,omega,c)
      alpha_0_value = alpha_0(q_value,omega,c)
      ! Call on the M_plus matrix element subroutine
      call M_plus_element(M_Plus_Array, nIterator, mIterator, q_value, p_value, n_values, m_values, alpha_value,&
       alpha_0_value, kappa_value, N_max, N_x, N_q, FactorialArray, Q_ArrayShifted)
    end do
  end do

  ! CUDA kernel?

  ! Endloop over M⁺


  ! Start loop over M⁻

  ! Loop over p values
  do mIterator = 1, N_q+1
    ! Assign current p value
    p_value = p_values_vector(mIterator)
    ! Calculate alpha(p) and alpha_0(k)
    alpha_value = alpha(p_value,epsilon,mu,omega,c)
    alpha_0_value = alpha_0(k_value,omega,c)
    ! Call on the M_minus vector element subroutine
    call M_minus_element(M_Minus_Array, mIterator, p_value, k_value, m_values, k_index, alpha_value, alpha_0_value,&
     kappa_value, N_max, N_x, N_q, FactorialArray, Q_ArrayShifted)
  end do

  ! CUDA kernel?

  ! End loop over M⁻


  ! Solve Ax=b or M⁺R = M⁻
  R_Vector = M_Minus_Array
  call zgesv(N_q+1, 1, M_Plus_Array, N_q+1, pivot, R_Vector, N_q+1, rc)

  ! Extract the real-part of q and calculate its corresponding angle
  realq = q_values_vector%re
  theta_array = asin(c*realq/omega)

  ! Calculate R_DifferentialMean for a SINGLE realization
  call CalculateNormDiffR(NormRDiffTheta, R_Vector, theta_0, theta_array, N_q, omega, pi, c, L)

  ! Remerge OpenMP -------------------------------

  print*, "Saving data!"





  ! Save data and reloop?
  call writeHDF5("R_DIff_Single_Plane_angle.h5",NormRDiffTheta, theta_array,N_q+1)






  ! Definine the functions to be used
  contains

  function alpha(q,epsilon,mu,omega,c) result(res)
    use, intrinsic :: iso_c_binding
    ! Input
    complex(C_DOUBLE_COMPLEX), intent (in) :: epsilon, mu, q
    real(C_DOUBLE), intent(in)             :: omega, c
    ! Output
    complex(C_DOUBLE_COMPLEX)              :: res
  
    res = sqrt(epsilon*mu*((omega*omega)/(c*c)) - q*q)
  
  end function

  function alpha_0(q,omega,c) result(res)
    use, intrinsic :: iso_c_binding
    ! Input
    complex(C_DOUBLE_COMPLEX), intent (in) :: q
    real(C_DOUBLE), intent(in)             :: omega, c
    ! Output
    complex(C_DOUBLE_COMPLEX)              :: res
  
    res = sqrt((omega*omega)/(c*c) - q*q)
  
  end function

  function kappa(nu, epsilon, mu) result(res)
    use, intrinsic :: iso_c_binding
    ! Input
    integer, intent(in)                    :: nu
    complex(C_DOUBLE_COMPLEX), intent(in)  :: epsilon, mu
    ! Output
    complex(C_DOUBLE_COMPLEX)              :: res

    if ( nu == 0 ) then
      ! p polarization
      res = epsilon
    else if ( nu == 1 ) then
      ! s polarization
      res = mu
    else
      ! Error as nu is neither 0 or 1!
      print*, "Neither p or s polarized light!"
      stop 1

    end if
  end function


end program SolveReyleighEquation


subroutine M_plus_element(M_Plus_Array, n, m, q_value, p_value, q_indexes, p_indexes, alpha_value, alpha_0_value,&
   kappa_value, N_max, N_x, N_q, FactorialArray, Q_ArrayShifted)
  use, intrinsic :: iso_c_binding
  ! Input
  integer, intent(in)                                             :: n, m, N_x, N_max, N_q
  integer, dimension(N_q+1), intent(in)                           :: p_indexes, q_indexes
  integer, dimension(N_max), intent(in)                           :: FactorialArray
  complex(C_DOUBLE_COMPLEX), intent(in)                           :: p_value, q_value, alpha_value, alpha_0_value, kappa_value
  complex(C_DOUBLE_COMPLEX), dimension(N_x,N_max)                 :: Q_ArrayShifted
  ! Output
  complex(C_DOUBLE_COMPLEX), dimension(N_q+1, N_q+1), intent(out) :: M_Plus_Array
  ! Allocate
  integer                                                         :: iterator, Q_index
  complex(C_DOUBLE_COMPLEX)                                       :: I_Integral, i = complex(0,1), PrefactorM

  ! q and p_indexes \in [\pm N_q], Q_indexes \in [\pm N_x] = [\pm 2N_q]
  Q_index = modulo(p_indexes(m) - q_indexes(n) + N_q, 2*N_q) + 1
  
  ! Iterate over sum stemming from Taylor series of FT-argument
  I_Integral = 0
  do iterator = 1, N_max
    I_Integral = I_Integral + (((i*(alpha_value-alpha_0_value))**iterator)/(FactorialArray(iterator)))*&
    Q_ArrayShifted(Q_index,iterator)
  end do

  
  PrefactorM = ((p_value+kappa_value*q_value)*(p_value-q_value))/(alpha_value-alpha_0_value) + alpha_value +&
   kappa_value*alpha_0_value

  M_Plus_Array(n,m) = PrefactorM*I_Integral
  

end subroutine

subroutine M_minus_element(M_Minus_Array, m, p_value, k_value, p_indexes, k_index, alpha_value, alpha_0_value,&
   kappa_value, N_max, N_x, N_q, FactorialArray, Q_ArrayShifted)
  use, intrinsic :: iso_c_binding
  ! Input
  integer, intent(in)                                           :: m, N_max, N_q, k_index
  integer, dimension(N_q+1), intent(in)                         :: p_indexes
  integer, dimension(N_max), intent(in)                         :: FactorialArray
  complex(C_DOUBLE_COMPLEX), intent(in)                         :: p_value, k_value, alpha_value, alpha_0_value, kappa_value
  complex(C_DOUBLE_COMPLEX), dimension(N_x,N_max), intent(in)   :: Q_ArrayShifted
  ! Output
  complex(C_DOUBLE_COMPLEX), dimension(N_q+1), intent(out)        :: M_Minus_Array
  ! Allocate
  integer                                                       :: iterator, Q_index
  complex(C_DOUBLE_COMPLEX)                                     :: I_Integral, i = complex(0,1), PrefactorM


  ! k and p_indexes \in [\pm N_q], Q_indexes \in [\pm N_x] = [\pm 2N_q]
  Q_index = modulo(p_indexes(m) - k_index + N_q, 2*N_q) + 1

  ! Iterate over sum stemming from Taylor series of FT-argument
  I_Integral = 0
  do iterator = 1, N_max
    I_Integral = I_Integral + (((i*(alpha_value+alpha_0_value))**iterator)/(FactorialArray(iterator)))*&
    Q_ArrayShifted(Q_index,iterator)
  end do

  
  PrefactorM = -(((p_value+kappa_value*k_value)*(p_value-k_value))/(alpha_value+alpha_0_value) + alpha_value -&
   kappa_value*alpha_0_value)

  M_Minus_Array(m) = PrefactorM*I_Integral
  

end subroutine

subroutine CalculateNormDiffR(NormRDiffTheta, R_Vector, theta_0, theta_array, N_q, omega, pi, c, L)
  use, intrinsic :: iso_c_binding
  ! Input
  integer, intent(in)                                     :: N_q
  real(C_DOUBLE), intent(in)                              :: theta_0, omega, pi, c, L
  real(C_DOUBLE), dimension(N_q+1), intent(in)            :: theta_array
  complex(C_DOUBLE_COMPLEX), dimension(N_q+1), intent(in) :: R_Vector
  ! Output
  real(C_DOUBLE), dimension(N_q+1), intent(out)           :: NormRDiffTheta
  ! Allocate
  integer                                                 :: Iterator
  real(C_DOUBLE)                                          :: AbsoluteSquare, costheta_0

  costheta_0 = cos(theta_0)

  do iterator = 1, N_q+1
    AbsoluteSquare = R_Vector(iterator)%re**2 + R_Vector(iterator)%im**2
    NormRDiffTheta(iterator) = ((1)/(L)) * ((omega)/(2*pi*c)) * ((cos(theta_array(iterator))**2)/(costheta_0)) &
     * AbsoluteSquare
  end do

end subroutine


subroutine CreateSurface(Surface, N_x)
  use, intrinsic :: iso_c_binding
  ! Input
  integer, intent(in)                                    :: N_x
  ! Output
  complex(C_DOUBLE_COMPLEX), dimension(N_x), intent(out) :: Surface
  ! Allocate
  integer                                                :: iterator

  do iterator = 1, N_x
    ! Flat surface
    Surface(iterator) = 1
  end do

end subroutine


subroutine fftshift(array,N)
  use, intrinsic :: iso_c_binding
  ! Input
  integer, intent(in)                                    :: N
  ! Input/Output
  complex(C_DOUBLE_COMPLEX), dimension(N), intent(inout) :: array
  ! Allocate
  integer                                                :: halfway
  complex(C_DOUBLE_COMPLEX), dimension(N)                :: arraycopy

  arraycopy = array

  halfway = ceiling((0.0+N)/2)

  do iterator = 1, N
    array(iterator) = arraycopy(modulo(iterator+halfway-1,N)+1)
  end do
  
end subroutine

subroutine FactorialArrayCalculation(FactorialArray,N_max)
use, intrinsic ::iso_c_binding
! Input
integer, intent(in)                    :: N_max
! Output
integer, dimension(N_max), intent(out) :: FactorialArray
! Allocate
integer                                :: iterator, res = 1

!Calculate the factorial of the numbers 1 through N_max once to save computation time
do iterator = 1, N_max
  res = res*iterator
  FactorialArray(iterator) = res
end do

end subroutine

subroutine writeHDF5(fileName,xarray,yarray,N)
  use hdf5
  ! Input
  integer, intent(in)              :: N
  real*8, dimension(N), intent(in) :: xarray, yarray
  character(len=*), intent(in)     :: fileName
  ! Allocate
  integer(4)                       :: error
  integer                          :: space_rank
  integer(HSIZE_T)                 :: data_dims(1)
  integer(HID_T)                   :: file_id, dspace_id, dset_id1, dset_id2, dset_id3, dset_id4

  !Interface
  call h5open_f(error)
  !Open file
  call h5fcreate_f(fileName,H5F_ACC_TRUNC_F,file_id,error)

  !Set sizes
  space_rank = 1
  data_dims(1) = N

  !Open dataspace
  call h5screate_simple_f(space_rank,data_dims,dspace_id,error)


  !Create dataset
  call h5dcreate_f(file_id,"rdiff",H5T_NATIVE_DOUBLE,dspace_id,dset_id1,error)
  call h5dcreate_f(file_id,"theta",H5T_NATIVE_DOUBLE,dspace_id,dset_id2,error)


  !Write dataset
  call h5dwrite_f(dset_id1,H5T_NATIVE_DOUBLE,xarray,data_dims,error)
  call h5dwrite_f(dset_id2,H5T_NATIVE_DOUBLE,yarray,data_dims,error)


  !Close dataset
  call h5dclose_f(dset_id1,error)
  call h5dclose_f(dset_id2,error)


  !Close dataspace
  call h5sclose_f(dspace_id,error)

  !Close file
  call h5fclose_f(file_id,error)
  !Close interface
  call h5close_f(error)

end subroutine
