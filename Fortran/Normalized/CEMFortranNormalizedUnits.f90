
module params
  use, intrinsic :: iso_c_binding 
  implicit none
  type :: param_t
    integer                   :: k_index

    real(C_DOUBLE)            :: WaveLength, c, pi, omega, theta_0, L, deltax    
    
    complex(C_DOUBLE_COMPLEX) :: deltaq, deltap, deltak, mu, epsilon, kappa_value, k_value, Q_big
    
  end type param_t

  
contains


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

!Initialize
subroutine Initialize(param_type, N_x, N_q)
  implicit none
  type(param_t), intent(inout) :: param_type
  integer, intent(in)          :: N_x, N_q

  ! Define basic parameters
  param_type%pi         = 3.14159265            !1
  param_type%WaveLength = 632.8E-9                !m
  param_type%c          = 3.0E8                 !m/s
  param_type%omega      = 2*param_type%pi*param_type%c/param_type%WaveLength     !1/s
  param_type%Q_big      = 4*(param_type%omega/param_type%c)           !1/m
  
  ! Define material properties and polarization
  param_type%epsilon     = 2.2487878787878787878787 !-17.5+0.48*complex(0,1)                      !1
  param_type%mu          = 1                                          !1
  param_type%kappa_value = kappa(1,param_type%epsilon,param_type%mu)  !1 0=p 1=s

  ! Define deltas
  param_type%deltaq = param_type%Q_big/N_q                            !1/m
  param_type%deltap = param_type%Q_big/N_q                            !1/m
  param_type%deltak = param_type%deltap                               !1/m
  param_type%deltax = param_type%pi/param_type%Q_big%re               !m

  ! Define k-related values
  param_type%k_index = 0                                            !1
  param_type%k_value = param_type%k_index*param_type%deltak           !1/m
  param_type%theta_0 = asin(param_type%k_value%re*param_type%c/param_type%omega) !rad

  ! Define length of sample
  param_type%L = param_type%deltax*N_x                                !m

end subroutine Initialize

subroutine Scale_Length( param_type, scale_const )
  implicit none
  type(param_t), intent(inout)        :: param_type
  real(C_DOUBLE_COMPLEX), intent(in)  :: scale_const

  param_type%Q_big    = param_type%Q_big* (1/scale_const) !1/m

  param_type%deltaq   = param_type%deltaq* (1/scale_const) !1/m
  param_type%deltap   = param_type%deltap* (1/scale_const) !1/m
  param_type%deltak   = param_type%deltak* (1/scale_const) !1/m
  param_type%deltax   = param_type%deltax*scale_const   !m

  param_type%k_value  = param_type%k_value*(1/scale_const) !1/m

  param_type%L        = param_type%L*scale_const   !m

end subroutine Scale_Length

end module params

! ---------------------------------------------------------------------------------------------------------------------

program SolveReyleighEquation
  use, intrinsic :: iso_c_binding
  use params
  implicit none
  include 'fftw3.f03'

  ! Declare global parameters
  ! N_max is the max power and end of sum for series expansion in the I integral
  integer, parameter                                :: N_x = 800  , N_q = N_x/2, N_max = 20, endpoint = N_q/4+1
  type(param_t)                                     :: parameters
  real, dimension(N_q+1)                            :: pivot
  real(C_DOUBLE), dimension(endpoint)                  :: Reflection
  real(C_DOUBLE), dimension(N_q+1)                  :: realq, theta_array, NormRDiffTheta

  ! Declare and allocate global arrays and such
  integer                                           :: iterator, Iterator2, mIterator, nIterator, rc
  integer, dimension(N_q+1)                         :: n_values, m_values
  integer*16, dimension(N_max)                      :: FactorialArray
  
  complex(C_DOUBLE_COMPLEX)                         :: q_value, p_value, alpha_value, alpha_0_value
  complex(C_DOUBLE_COMPLEX), dimension(N_q+1)       :: R_Vector, M_Minus_Array, q_values_vector, p_values_vector
  complex(C_DOUBLE_COMPLEX), dimension(N_x)         :: Surface, SurfacePower, Q_Vector
  complex(C_DOUBLE_COMPLEX), dimension(N_q+1,N_q+1) :: M_Plus_Array
  complex(C_DOUBLE_COMPLEX), dimension(N_x,N_max)   :: Q_ArrayShifted

  type(C_PTR)                                       :: plan

  ! ----------------------------------------------

  call Initialize(parameters, N_x, N_q)

  call Scale_Length(parameters, parameters%omega/parameters%c)

  ! Calculate factorials for later use
  call FactorialArrayCalculation(FactorialArray,N_max)

  

  do Iterator2 = 1, endpoint
    

    parameters%k_index = Iterator2-1
  
    parameters%k_value = parameters%k_index*parameters%deltak           !1/m
  
    parameters%theta_0 = asin(parameters%k_value%re) !rad

    print*, "-----------"
    print*, " Iteration: ", Iterator2
    print*, "Variables scaled displayed: "
    print*, "omega:               ", parameters%omega
    print*, "Q:      ", parameters%Q_big
    print*, "N_q:          ", N_q
    print*, "N_x:           ", N_x
    print*, "N_max:       ", N_max
    print*, "deltaq: ", parameters%deltaq
    print*, "deltap: ", parameters%deltap
    print*, "deltax:             ", parameters%deltax
    print*, "L:                  ", parameters%L
    print*, "Incident angle:     ", parameters%theta_0*180/parameters%pi
    print*, "k_value:", parameters%k_value
    print*, "epsilon", parameters%epsilon
    print*, "-----------"

    ! Define the q and p values and m and n indexes
    do iterator = 1, N_q+1
      n_values(iterator)        =  iterator-N_q/2-1
      m_values(iterator)        =  iterator-N_q/2-1
      p_values_vector(iterator) = (iterator-N_q/2-1)*parameters%deltap
      q_values_vector(iterator) = (iterator-N_q/2-1)*parameters%deltaq
    end do



    ! Calculate surface
    call CreateSurface(Surface,N_x)

    !print*, FactorialArray

    ! Loop here?
    ! Split using OpenMP ---------------------------



    ! Making an FFTW plan
    plan = fftw_plan_dft_1d(N_x, SurfacePower, Q_Vector, FFTW_FORWARD, FFTW_ESTIMATE)
    
    ! Calculate FFT of surface
    do iterator = 1, N_max
      
      !Taking the surface to the power of iterator
      SurfacePower = Surface**iterator
      
      ! Executing the plan
      call fftw_execute_dft(plan, SurfacePower, Q_Vector)
      
      ! Using the fftshift to move the elements around
      call fftshift(Q_Vector,N_x)
      ! Saving and multiplying by deltax
      Q_ArrayShifted(:,iterator) = Q_Vector*parameters%deltax

    end do

    ! Destroying plan
    call fftw_destroy_plan(plan)

    
    ! Start loop over M⁺
    ! First loop over p values
    do mIterator = 1, N_q+1
      ! Then loop over q values
      do nIterator = 1, N_q+1
        ! Assign current p and q value
        q_value = q_values_vector(nIterator)
        p_value = p_values_vector(mIterator)
        ! Calculating alpha(p) and alpha_0(q)
        alpha_value   = alpha(p_value,parameters)
        alpha_0_value = alpha_0(q_value)
        ! print*, alpha_0_value, alpha_value
        ! Call on the M_plus matrix element subroutine
        call M_plus_element(parameters, M_Plus_Array, nIterator, mIterator, q_value, p_value, n_values, m_values, alpha_value,&
        alpha_0_value, N_max, N_x, N_q, FactorialArray, Q_ArrayShifted)
        !print*, M_Plus_Array(mIterator,nIterator), alpha_0_value, alpha_value, p_value, q_value, nIterator, mIterator
      end do
    end do

    print*, "here"

    ! print*, M_Plus_Array

    ! CUDA kernel?

    ! Endloop over M⁺

    ! Start loop over M⁻

    ! Loop over p values
    do mIterator = 1, N_q+1
      ! Assign current p value
      p_value = p_values_vector(mIterator)
      ! Calculate alpha(p) and alpha_0(k)
      alpha_value = alpha(p_value,parameters)
      alpha_0_value = alpha_0(parameters%k_value)
      ! Call on the M_minus vector element subroutine
      call M_minus_element(parameters, M_Minus_Array, mIterator, p_value, m_values, alpha_value, alpha_0_value,&
      N_max, N_x, N_q, FactorialArray, Q_ArrayShifted)
    end do

    print*, "here2"

    ! CUDA kernel?

    ! End loop over M⁻

    ! Solve Ax=b or M⁺R = M⁻
    R_Vector = M_Minus_Array

    print*, "here33"

    call zgesv(N_q+1, 1, M_Plus_Array, N_q+1, pivot, R_Vector, N_q+1, rc)

    print*, "here3"

    ! Extract the real-part of q and calculate its corresponding angle
    realq = q_values_vector%re
    theta_array = asin(realq)

  !  print*, M_Plus_Array

    ! Calculate R_DifferentialMean for a SINGLE realization
    call CalculateNormDiffR(NormRDiffTheta, R_Vector, theta_array, N_q, parameters)

    ! Remerge OpenMP -------------------------------

  !  print*, sqrt(R_Vector(k_index+N_q/2+1)%re**2 + R_Vector(k_index+N_q/2+1)%im**2),&
  !  (parameters%omega/parameters%c), parameters%deltax, parameters%L

    Reflection(Iterator2) = sqrt(R_Vector(parameters%k_index+N_q/2+1)%re**2 + R_Vector(parameters%k_index+N_q/2+1)%im**2)

  end do

  print*, Reflection
  print*, "-----"
  print*, theta_array(N_q/2+1:N_q/2+endpoint)

  print*, "Saving data!"

  ! Save data and reloop?
  call writeHDF5("R_DIff_Single_Plane_s.h5",Reflection, theta_array(N_q/2+1:N_q*3/4),endpoint)


  ! Definine the functions to be used
  contains

  function alpha(q,parameters) result(res)
    use, intrinsic :: iso_c_binding
    ! Input
    type(param_t)                          :: parameters
    complex(C_DOUBLE_COMPLEX), intent (in) :: q
    ! Output
    complex(C_DOUBLE_COMPLEX)              :: res
  
    res = zsqrt(parameters%epsilon*parameters%mu - q*q)
  
  end function

  function alpha_0(q) result(res)
    use, intrinsic :: iso_c_binding
    ! Input
    complex(C_DOUBLE_COMPLEX), intent (in) :: q
    ! Output
    complex(C_DOUBLE_COMPLEX)              :: res
  
    res = zsqrt(1 - q*q)

  end function


end program SolveReyleighEquation

!----------------------------------------------------------------------------------------------------------------------

subroutine M_plus_element(parameters, M_Plus_Array, n, m, q_value, p_value, q_indexes, p_indexes, alpha_value, &
   alpha_0_value, N_max, N_x, N_q, FactorialArray, Q_ArrayShifted)
  use params
  use, intrinsic :: iso_c_binding
  ! Input
  type(param_t), intent(in)                                       :: parameters
  integer, intent(in)                                             :: n, m, N_x, N_max, N_q
  integer, dimension(N_q+1), intent(in)                           :: p_indexes, q_indexes
  integer*16, dimension(N_max), intent(in)                           :: FactorialArray
  complex(C_DOUBLE_COMPLEX), intent(in)                           :: p_value, q_value, alpha_value, alpha_0_value
  complex(C_DOUBLE_COMPLEX), dimension(N_x,N_max)                 :: Q_ArrayShifted
  ! Output
  complex(C_DOUBLE_COMPLEX), dimension(N_q+1, N_q+1), intent(out) :: M_Plus_Array
  ! Allocate
  integer                                                         :: Q_index
  complex(C_DOUBLE_COMPLEX)                                       :: I_Integral, PrefactorM

  ! q and p_indexes \in [\pm N_q], Q_indexes \in [\pm N_x] = [\pm 2N_q]
  Q_index = modulo(p_indexes(m) - q_indexes(n) + N_q, 2*N_q) + 1
  
  ! Iterate over sum stemming from Taylor series of FT-argument

  call Integral_plus(parameters, Q_index, Q_ArrayShifted, alpha_value, alpha_0_value, FactorialArray, I_Integral, N_max)

  ! Calculate the prefactor  
  PrefactorM = ((p_value+parameters%kappa_value*q_value)*(p_value-q_value))/(alpha_value-alpha_0_value) + &
   alpha_value + parameters%kappa_value*alpha_0_value
  !print*, PrefactorM

  M_Plus_Array(m,n) = PrefactorM*I_Integral
  

end subroutine

subroutine Integral_plus(parameters, Q_index, Q_ArrayShifted, alpha_value, alpha_0_value, FactorialArray, I_Integral,N_max)
  use params
  use, intrinsic :: iso_c_binding
  !Input
  type(param_t), intent(in)                                     :: parameters
  integer, intent(in) :: Q_index, N_max
  integer*16, dimension(N_max), intent(in)                      :: FactorialArray
  complex(C_DOUBLE_COMPLEX), dimension(N_x,N_max), intent(in)   :: Q_ArrayShifted
  complex(C_DOUBLE_COMPLEX), intent(in)                         :: alpha_value, alpha_0_value
  !Output
  complex(C_DOUBLE_COMPLEX), intent(inout)                                 :: I_Integral
  !Allocate
  integer                                                       :: iterator
  complex(C_DOUBLE_COMPLEX)                                     :: i = complex(0,1)

  I_Integral = 0
  do iterator = 1, N_max
    I_Integral = I_Integral + (((i*(alpha_value-alpha_0_value))**iterator)/(FactorialArray(iterator)))*&
    Q_ArrayShifted(Q_index,iterator)
  end do
  

end subroutine Integral_plus

subroutine M_minus_element(parameters, M_Minus_Array, m, p_value, p_indexes, alpha_value, &
   alpha_0_value, N_max, N_x, N_q, FactorialArray, Q_ArrayShifted)
  use params
  use, intrinsic :: iso_c_binding
  ! Input
  type(param_t), intent(in)                                     :: parameters
  integer, intent(in)                                           :: m, N_max, N_q
  integer, dimension(N_q+1), intent(in)                         :: p_indexes
  integer*16, dimension(N_max), intent(in)                         :: FactorialArray
  complex(C_DOUBLE_COMPLEX), intent(in)                         :: p_value, alpha_value, alpha_0_value
  complex(C_DOUBLE_COMPLEX), dimension(N_x,N_max), intent(in)   :: Q_ArrayShifted
  ! Output
  complex(C_DOUBLE_COMPLEX), dimension(N_q+1), intent(out)      :: M_Minus_Array
  ! Allocate
  integer                                                       :: Q_index
  complex(C_DOUBLE_COMPLEX)                                     :: I_Integral, PrefactorM


  ! k and p_indexes \in [\pm N_q], Q_indexes \in [\pm N_x] = [\pm 2N_q]
  Q_index = modulo(p_indexes(m) - parameters%k_index + N_q, 2*N_q) + 1

  ! Iterate over sum stemming from Taylor series of FT-argument
  call Integral_minus(parameters, Q_index, Q_ArrayShifted, alpha_value, alpha_0_value, FactorialArray, I_Integral, N_max)

  !print*, I_Integral, alpha_0_value, alpha_value

  ! Calculate the prefactor
  PrefactorM = -(((p_value+parameters%kappa_value*parameters%k_value)*(p_value-parameters%k_value))/(alpha_value+alpha_0_value) + &
   alpha_value - parameters%kappa_value*alpha_0_value)

  ! Finally calculate the element
  M_Minus_Array(m) = PrefactorM*I_Integral
  
end subroutine

subroutine Integral_minus(parameters, Q_index, Q_ArrayShifted, alpha_value, alpha_0_value, FactorialArray, I_Integral, N_max)
  use params
  use, intrinsic :: iso_c_binding
  !Input
  type(param_t), intent(in)                                     :: parameters
  integer, intent(in) :: Q_index, N_max
  integer*16, dimension(N_max), intent(in)                      :: FactorialArray
  complex(C_DOUBLE_COMPLEX), dimension(N_x,N_max), intent(in)   :: Q_ArrayShifted
  complex(C_DOUBLE_COMPLEX), intent(in)                         :: alpha_value, alpha_0_value
  !Output
  complex(C_DOUBLE_COMPLEX), intent(inout)                                 :: I_Integral
  !Allocate
  integer                                                       :: iterator
  complex(C_DOUBLE_COMPLEX)                                     :: i = complex(0,1)

  I_Integral = 0
  do iterator = 1, N_max
    I_Integral = I_Integral + (((i*(alpha_value+alpha_0_value))**iterator)/(FactorialArray(iterator)))*&
    Q_ArrayShifted(Q_index,iterator)
  end do

end subroutine Integral_minus

subroutine CalculateNormDiffR(NormRDiffTheta, R_Vector, theta_array, N_q, parameters)
  use params
  use, intrinsic :: iso_c_binding
  ! Input
  type(param_t)                                           :: parameters
  integer, intent(in)                                     :: N_q
  real(C_DOUBLE), dimension(N_q+1), intent(in)            :: theta_array
  complex(C_DOUBLE_COMPLEX), dimension(N_q+1), intent(in) :: R_Vector
  ! Output
  real(C_DOUBLE), dimension(N_q+1), intent(out)           :: NormRDiffTheta
  ! Allocate
  integer                                                 :: Iterator
  real(C_DOUBLE)                                          :: AbsoluteSquare, costheta_0

  costheta_0 = cos(parameters%theta_0)

  do iterator = 1, N_q+1
    AbsoluteSquare = R_Vector(iterator)%re**2 + R_Vector(iterator)%im**2
    NormRDiffTheta(iterator) = ((1)/(parameters%L)) * (1/(2*parameters%pi)) * ((cos(theta_array(iterator))**2)/(costheta_0)) &
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
    Surface(iterator) = 10
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
integer*16, dimension(N_max), intent(out) :: FactorialArray
! Allocate
integer                                :: iterator
integer*16                             :: res = 1

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