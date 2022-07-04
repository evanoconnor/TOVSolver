!-*-f90-*-
subroutine mass_polytrope(central_density,gamma,K,out_mass,out_bmass,out_radius,N,dr)

  implicit none

  !inputs
  real*8 central_density
  real*8 out_mass, out_radius, out_pradius,out_bmass
  real*8 dr !step size
  integer N !maximum number of steps
  real*8 gamma, K
  real*8 temp_eps,temp_rho

  !outputables
  real*8,allocatable :: TOVpressure(:)
  real*8,allocatable :: TOVrho(:)
  real*8,allocatable :: TOVphi(:)
  real*8,allocatable :: TOVrad(:)
  real*8,allocatable :: TOVprad(:)
  real*8,allocatable :: TOVmass(:)
  real*8,allocatable :: TOVbmass(:)
  real*8,allocatable :: TOVeps(:)

  !internal
  integer i
  real*8 h
  real*8 :: pi = 3.14159265358979323846d0
  real*8 :: c = 29979245800.0d0
  real*8 :: G = 6.6742d-8

  !RK variables
  real*8 k1p,k2p,k3p,k4p
  real*8 k1m,k2m,k3m,k4m
  real*8 k1ph,k2ph,k3ph,k4ph
  real*8 k1mb,k2mb,k3mb,k4mb

  allocate(TOVpressure(N))
  allocate(TOVrho(N))
  allocate(TOVphi(N))
  allocate(TOVrad(N))
  allocate(TOVprad(N))
  allocate(TOVmass(N))
  allocate(TOVbmass(N))
  allocate(TOVeps(N))

  TOVpressure(:) = 0.0d0
  TOVrho(:) = 0.0d0
  TOVphi(:) = 0.0d0
  TOVrad(:) = 0.0d0
  TOVprad(:) = 0.0d0
  TOVmass(:) = 0.0d0
  TOVbmass(:) = 0.0d0
  TOVeps(:) = 0.0d0

  !first get center cell
  !find ye given rho and temp, this is call to ott_eos with keytemp=2
  TOVrho(1) = central_density
  TOVrad(1) = 0.0d0
  TOVprad(1) = 0.0d0
  TOVmass(1) = 0.0d0
  TOVbmass(1) = 0.0d0
  TOVphi(1) = 0.0d0


  call press(TOVpressure(1),TOVrho(1),gamma,K)
  call eps(TOVeps(1),TOVpressure(1),gamma,K)

!  write(*,*) TOVpressure(1),TOVeps(1)
!  stop
  
  !now do RK4 to get new values of press, mass and phi
  do i=2,N
!     write(*,*) i-1,TOVmass(i-1),TOVrad(i-1),0.0d0,TOVrho(i-1),0.0d0,TOVpressure(i-1),TOVeps(i-1)
 !    write(*,*) i-1,TOVmass(i-1),TOVrad(i-1),0.0d0,TOVrho(i-1),0.0d0,0.0d0,0.0d0
     h = 1.0d0+TOVeps(i-1)/c**2+TOVpressure(i-1)/TOVrho(i-1)/c**2
     call dpdr(TOVmass(i-1),TOVpressure(i-1),TOVrho(i-1)*h,TOVrad(i-1),k1p)
     call dmdr(TOVpressure(i-1),TOVrho(i-1)*h,TOVrad(i-1),k1m)
     call dphidr(TOVmass(i-1),TOVpressure(i-1),TOVrho(i-1)*h,TOVrad(i-1),k1ph)
     call dmbdr(TOVmass(i-1),TOVrho(i-1),TOVrad(i-1),k1mb)
    
     if (TOVpressure(i-1)+dr*0.5d0*k1p.lt.0.0d0) then
        out_mass = TOVmass(i-1)
        out_bmass = TOVbmass(i-1)
        out_radius = TOVrad(i-1)
        return
     endif
    
     call rho(temp_rho,TOVpressure(i-1)+0.5d0*dr*k1p,gamma,K)
     call eps(temp_eps,TOVpressure(i-1)+0.5d0*dr*k1p,gamma,K)

     h = 1.0d0+temp_eps/c**2+(TOVpressure(i-1)+dr*0.5d0*k2p)/temp_rho/c**2
     call dpdr(TOVmass(i-1)+0.5d0*dr*k1m,TOVpressure(i-1)+0.5d0*dr*k1p, &
          temp_rho*h,TOVrad(i-1)+0.5d0*dr,k2p)
     call dmdr(TOVpressure(i-1)+0.5d0*dr*k1p,temp_rho*h,TOVrad(i-1)+0.5d0*dr,k2m)
     call dphidr(TOVmass(i-1)+0.5d0*dr*k1m,TOVpressure(i-1)+0.5d0*dr*k1p, &
          temp_rho*h,TOVrad(i-1)+0.5d0*dr,k2ph)
     call dmbdr(TOVmass(i-1)+0.5d0*dr*k1m,temp_rho,TOVrad(i-1)+0.5d0*dr,k2mb)

     if (TOVpressure(i-1)+dr*0.5d0*k2p.lt.0.0d0) then
        out_mass = TOVmass(i-1)
        out_bmass = TOVbmass(i-1)
        out_radius = TOVrad(i-1)
        return
     endif

     call rho(temp_rho,TOVpressure(i-1)+0.5d0*dr*k2p,gamma,K)
     call eps(temp_eps,TOVpressure(i-1)+0.5d0*dr*k2p,gamma,K)

     h = 1.0d0+temp_eps/c**2+(TOVpressure(i-1)+dr*0.5d0*k2p)/temp_rho/c**2
     call dpdr(TOVmass(i-1)+dr*0.5d0*k2m,TOVpressure(i-1)+dr*0.5d0*k2p, &
          temp_rho*h,TOVrad(i-1)+0.5d0*dr,k3p)
     call dmdr(TOVpressure(i-1)+dr*0.5d0*k2p,temp_rho*h,TOVrad(i-1)+0.5d0*dr,k3m)
     call dphidr(TOVmass(i-1)+dr*0.5d0*k2m,TOVpressure(i-1)+dr*0.5d0*k2p, &
          temp_rho*h,TOVrad(i-1)+0.5d0*dr,k3ph)
     call dmbdr(TOVmass(i-1)+0.5d0*dr*k2m,temp_rho,TOVrad(i-1)+0.5d0*dr,k3mb)

     if (TOVpressure(i-1)+dr*k3p.lt.0.0d0) then
        out_mass = TOVmass(i-1)
        out_bmass = TOVbmass(i-1)
        out_radius = TOVrad(i-1)
        return
     endif


     call rho(temp_rho,TOVpressure(i-1)+dr*k3p,gamma,K)
     call eps(temp_eps,TOVpressure(i-1)+dr*k3p,gamma,K)

     h = 1.0d0+temp_eps/c**2+(TOVpressure(i-1)+dr*k3p)/temp_rho/c**2
     call dpdr(TOVmass(i-1)+dr*k3m,TOVpressure(i-1)+dr*k3p, &
          temp_rho*h,TOVrad(i-1)+dr,k4p)
     call dmdr(TOVpressure(i-1)+dr*k3p,temp_rho*h,TOVrad(i-1)+dr,k4m)
     call dphidr(TOVmass(i-1)+dr*k3m,TOVpressure(i-1)+dr*k3p, &
          temp_rho*h,TOVrad(i-1)+dr,k4ph)
     call dmbdr(TOVmass(i-1)+dr*k3m,temp_rho,TOVrad(i-1)+dr,k4mb)

     TOVpressure(i) = TOVpressure(i-1) + dr/6.0d0*(k1p + 2.0d0*k2p + &
          2.0d0*k3p + k4p)
     TOVmass(i) = TOVmass(i-1) + dr/6.0d0*(k1m + 2.0d0*k2m + &
          2.0d0*k3m + k4m)
     TOVrad(i) = real(i-1)*dr
     TOVprad(i) = TOVprad(i-1) + dr*sqrt(1.0d0-2.0d0*G*TOVmass(i)/c**2/TOVrad(i))
     TOVphi(i) = TOVphi(i-1) + dr/6.0d0*(k1ph + 2.0d0*k2ph + &
          2.0d0*k3ph + k4ph)
     TOVbmass(i) = TOVbmass(i-1) + dr/6.0d0*(k1mb + 2.0d0*k2mb + &
          2.0d0*k3mb + k4mb)

!     write(*,*) TOVpressure(i),temp_rho
     if (TOVpressure(i).ne.TOVpressure(i)) then
        stop "NaN"
     else

     endif

     call rho(TOVrho(i),TOVpressure(i),gamma,K)
     call eps(TOVeps(i),TOVpressure(i),gamma,K)

     if (TOVrho(i).lt.1.0d2) then
        out_mass = TOVmass(i)
        out_bmass = TOVbmass(i)
        out_radius = TOVrad(i)
        out_pradius = TOVprad(i)
        return
     endif

     if (TOVpressure(i).lt.0.0d0) then
        out_mass = TOVmass(i-1)
        out_bmass = TOVbmass(i-1)
        out_radius = TOVrad(i-1)
        out_pradius = TOVprad(i-1)
        return
     endif

  enddo

  out_mass = TOVmass(N-1)
  out_bmass = TOVbmass(N-1)
  out_radius = N*dr

end subroutine mass_polytrope

subroutine rho(outrho,inpress,gamma,K)

  implicit none
  
  !inputs
  real*8 inpress, gamma, K
  
  !output
  real*8 outrho

  !internal
  real*8,parameter :: press_cgs_to_gc = 1.80171810d-39
  real*8,parameter :: rho_cgs_to_gc = 1.61930347d-18
  
  outrho = (inpress*press_cgs_to_gc/K)**(1.0d0/gamma)/rho_cgs_to_gc

end subroutine rho

subroutine eps(outeps,inpress,gamma,K)

  implicit none
  
  !inputs
  real*8 inpress, gamma, K
  
  !output
  real*8 outeps
  
  !internal
  real*8 inrho
  real*8,parameter :: press_cgs_to_gc = 1.80171810d-39
  real*8,parameter :: rho_cgs_to_gc = 1.61930347d-18
  real*8,parameter :: eps_cgs_to_gc = 1.11265006d-21

  call rho(inrho,inpress,gamma,K)

  outeps = inpress*press_cgs_to_gc/(inrho*rho_cgs_to_gc)/(gamma-1.0d0)/eps_cgs_to_gc


end subroutine eps

subroutine press(outpress,inrho,gamma,K)

  implicit none
  
  !inputs
  real*8 inrho, gamma, K
  
  !output
  real*8 outpress

  !internal
  real*8,parameter :: press_cgs_to_gc = 1.80171810d-39
  real*8,parameter :: rho_cgs_to_gc = 1.61930347d-18
  
  outpress = K*(inrho*rho_cgs_to_gc)**(gamma)/press_cgs_to_gc

!  write(*,*) K,inrho,outpress,K*rho_cgs_to_gc**2.0d0/press_cgs_to_gc
!  stop
  
end subroutine press
