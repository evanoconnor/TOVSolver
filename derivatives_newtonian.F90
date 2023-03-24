!-*-f90-*-
subroutine dpdr_newtonian(mass,rho,rad,xx)

  implicit none

  !inputs
  real*8 mass,rho,rad

  !outputs
  real*8 xx

  !internal
  real*8 :: G = 6.6742d-8

  if (mass.eq.0.0d0) then
     xx = 0.0d0
  else
     xx = -G*rho*mass/rad**2
  endif

end subroutine dpdr_newtonian

subroutine dmdr_newtonian(rho,rad,xx)

  implicit none

  !inputs
  real*8 mass,rho,rad

  !Output
  real*8 xx

  !internal
  real*8 :: pi = 3.14159265d0

  xx = 4.0d0*pi*rad**2*rho

end subroutine dmdr_newtonian

subroutine dphidr_newtonian(mass,rad,xx)

  implicit none

  !inputs
  real*8 mass,rad

  !outputs
  real*8 xx
  
  !internal
  real*8 :: G = 6.6742d-8
  
  xx = G*mass/rad**2


end subroutine dphidr_newtonian

