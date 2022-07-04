!-*-f90-*-
module hybrid_eos_module

  real*8 gamma1
  real*8 gamma2
  real*8 gammath
  real*8 K1,K2,E1,E2,E3
  real*8 rhonuc

end module hybrid_eos_module

subroutine mass_hybrid(central_density,out_mass,out_bmass,out_radius,N,dr)

  implicit none

  !inputs
  real*8 central_density
  real*8 out_mass, out_radius, out_pradius,out_bmass
  real*8 dr !step size
  integer N !maximum number of steps
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

  call press_hy(TOVpressure(1),TOVrho(1))
  call eps_hy(TOVeps(1),TOVpressure(1))
  
  !now do RK4 to get new values of press, mass and phi
  do i=2,N
!     write(*,*) i-1,TOVmass(i-1),TOVrad(i-1),0.0d0,TOVrho(i-1),0.0d0,0.0d0, &
!          pi/(1.0d0+(TOVrad(i-1)/1.0d8)**2)
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
    
     call rho_hy(temp_rho,TOVpressure(i-1)+0.5d0*dr*k1p)
     call eps_hy(temp_eps,TOVpressure(i-1)+0.5d0*dr*k1p)

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

     call rho_hy(temp_rho,TOVpressure(i-1)+0.5d0*dr*k2p)
     call eps_hy(temp_eps,TOVpressure(i-1)+0.5d0*dr*k2p)

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


     call rho_hy(temp_rho,TOVpressure(i-1)+dr*k3p)
     call eps_hy(temp_eps,TOVpressure(i-1)+dr*k3p)

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

     call rho_hy(TOVrho(i),TOVpressure(i))
     call eps_hy(TOVeps(i),TOVpressure(i))

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

end subroutine mass_hybrid

subroutine rho_hy(outrho,inpress)

  use hybrid_eos_module
  implicit none
  
  !inputs
  real*8 inpress
  
  !output
  real*8 outrho

  !internal
  
  real*8 Kx,Ex,Gx,Ex3

  if (inpress.lt.K1*rhonuc**gamma1) then
     Kx=K1
     Ex=E1
     Gx=gamma1
     Ex3=0.d0
  else
     Kx=K2
     Ex=E2
     Gx=gamma2
     Ex3=E3
  endif

  outrho = (inpress/Kx)**(1.0d0/Gx)

end subroutine rho_hy

subroutine eps_hy(outeps,inpress)

  use hybrid_eos_module
  implicit none
  
  !inputs
  real*8 inpress
  
  !output
  real*8 outeps
  
  !internal
  real*8 inrho
  real*8 Kx,Ex,Gx,Ex3

  if(inrho .lt. rhonuc) then
     Kx=K1
     Ex=E1
     Gx=gamma1
     Ex3=0.d0
  else
     Kx=K2
     Ex=E2
     Gx=gamma2
     Ex3=E3
  endif

  call rho_hy(inrho,inpress)

  outeps = inpress/inrho/(Gx-1.0d0)

end subroutine eps_hy

subroutine press_hy(outpress,inrho)

  use hybrid_eos_module
  implicit none
  
  !inputs
  real*8 inrho
  
  !output
  real*8 outpress

  !internal
  real*8 Kx,Ex,Gx,Ex3

  if(inrho .lt. rhonuc) then
     Kx=K1
     Ex=E1
     Gx=gamma1
     Ex3=0.d0
  else
     Kx=K2
     Ex=E2
     Gx=gamma2
     Ex3=E3
  endif
  
  outpress = Kx*inrho**Gx

end subroutine press_hy

subroutine init_hybrid_eos

  use hybrid_eos_module
  implicit none

  E1 = K1/(gamma1-1.d0)
  E2 = (gamma1 - 1.d0)/(gamma2-1.d0)*E1*rhonuc**(gamma1-gamma2)
  K2 = (gamma2 - 1.d0)*E2
  write(*,*) K2
  stop
  E3 = (gamma2 - gamma1)/(gamma2-1.d0)*E1*rhonuc**(gamma1-1.d0)

end subroutine init_hybrid_eos
