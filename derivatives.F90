!-*-f90-*-
subroutine dpdr(mass,press,rho,rad,xx)

  implicit none

  !inputs
  real*8 mass,press,rho,rad

  !outputs
  real*8 xx

  !internal
  real*8 :: G = 6.6742d-8
  real*8 :: c = 29979245800.0d0
  real*8 :: pi = 3.14159265d0

  if (mass.eq.0.0d0) then
     xx = -G*rho*4.0d0*pi*press*rad/c**2
  else
     xx = -G*rho*mass/rad**2*(1.0d0+4.0d0*pi*press/c/c*rad**3/mass)* &
          (1.0d0-2.0d0*G*mass/c/c/rad)**(-1)
  endif

end subroutine dpdr

subroutine dmdr(press,rho,rad,xx)

  implicit none

  !inputs
  real*8 mass,press,rho,rad

  !Output
  real*8 xx

  !internal
  real*8 :: pi = 3.14159265d0
  real*8 :: c = 29979245800.0d0

  xx = 4.0d0*pi*rad**2*(rho-press/c**2)

end subroutine dmdr

subroutine dmbdr(mass,rho,rad,xx)

  implicit none

  !inputs
  real*8 mass,press,rho,rad

  !Output
  real*8 xx

  !internal
  real*8 :: pi = 3.14159265d0
  real*8 :: c = 29979245800.0d0
  real*8 :: G = 6.6742d-8

  if (mass.eq.0.0d0) then
     xx = 4.0d0*pi*rad**2*rho
  else
     xx = 4.0d0*pi*rad**2*rho*(1.0d0-2.0d0*G*mass/(rad*c**2))**(-0.5d0)
  endif

end subroutine dmbdr

subroutine dphidr(mass,press,rho,rad,xx)

  implicit none

  !inputs
  real*8 mass,press,rho,rad

  !outputs
  real*8 xx
  
  !internal
  real*8 :: pi = 3.14159265d0
  real*8 :: G = 6.6742d-8
  real*8 :: c = 29979245800.0d0
  
  if (mass.eq.0.0d0) then
     xx = G*4.0d0*pi*press*rad/c**4
  else
     xx = (G/c**2)*mass/rad**2*(1.0d0+4.0d0*pi*press/c/c*rad**3/mass)* &
          (1.0d0-2.0d0*G*mass/c/c/rad)**(-1)
  endif


end subroutine dphidr

