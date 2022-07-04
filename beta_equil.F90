!-*-f90-*-
module betaequil 

  implicit none

  real*8, allocatable :: profileeps(:)
  real*8, allocatable :: profilerho(:)
  real*8, allocatable :: profilepressure(:)
  
end module betaequil

subroutine mass_betaequilEOS(central_density,out_mass,out_bmass,out_radius,out_pradius,N,dr)

  implicit none

  !inputs
  real*8 central_density
  real*8,intent(out) :: out_mass, out_radius, out_pradius,out_bmass
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

!  write(*,*) TOVrho(1)
  write(*,*) ""
  call press0(TOVpressure(1),TOVrho(1))
!     write(*,*) "two",TOVpressure(1)
  call eps0(TOVeps(1),TOVpressure(1))
!     write(*,*) "three"
!  write(*,*) TOVrho(1),TOVpressure(1),TOVeps(1)
!  stop
  
  !now do RK4 to get new values of press, mass and phi
  do i=2,N
!write(*,*) i
  !        write(*,"(1P10E18.9)") i-1,TOVmass(i-1),TOVrad(i-1),0.0d0,TOVrho(i-1),0.0d0,TOVpressure(i-1),TOVeps(i-1)
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
    
     call rho0(temp_rho,TOVpressure(i-1)+0.5d0*dr*k1p)
     call eps0(temp_eps,TOVpressure(i-1)+0.5d0*dr*k1p)

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

     call rho0(temp_rho,TOVpressure(i-1)+0.5d0*dr*k2p)
     call eps0(temp_eps,TOVpressure(i-1)+0.5d0*dr*k2p)

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

     call rho0(temp_rho,TOVpressure(i-1)+dr*k3p)
     call eps0(temp_eps,TOVpressure(i-1)+dr*k3p)

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

!     write(*,*) TOVpressure(i)
     call rho0(TOVrho(i),TOVpressure(i))
!     write(*,*) TOVpressure(i)
!     stop

     call eps0(TOVeps(i),TOVpressure(i))

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

end subroutine mass_betaequilEOS

subroutine rho0(outrho,inpress)

  use betaequil
  implicit none
  
  !inputs
  real*8 inpress
  
  !output
  real*8 outrho

  !internal
  real*8,parameter :: press_cgs_to_gc = 1.80171810d-39
  real*8,parameter :: rho_cgs_to_gc = 1.61930347d-18
  integer i
  logical :: cont = .true.

  !first find points that bracket press
  cont = .true.
!  write(*,*) "incoming" ,inpress
  i = 1
  do while(cont)
     if (profilepressure(i+1).gt.inpress) then
!        write(*,*) i,profilepressure(i)
        outrho = 10.0d0**(log10(profilerho(i)) + log10(profilerho(i+1)/profilerho(i))/ &
             log10(profilepressure(i+1)/profilepressure(i))*log10(inpress/profilepressure(i)))
!        write(*,*) i,profilepressure(i),outrho
        cont = .false.
     else
        i=i+1
     endif
!     write(*,*) i
  enddo

!  write(*,*) i

end subroutine rho0

subroutine eps0(outeps,inpress)

  use betaequil
  implicit none
  
  !inputs
  real*8 inpress
  
  !output
  real*8 outeps
  
  !internal
  real*8 inrho
  real*8,parameter :: press_cgs_to_gc = 1.80171810d-39
  real*8,parameter :: rho_cgs_to_gc = 1.61930347d-18
  real*8,parameter :: eps_cgs_to_gc = 1.11265006d-21

  integer i
  logical :: cont = .true.

  !first find points that bracket press
  i = 1
  cont = .true.
 ! profileeps = profileeps + 1.0d19
  do while(cont)
     if (profilepressure(i+1).gt.inpress) then
        outeps = 10.0d0**(log10(profileeps(i)) + log10(profileeps(i+1)/profileeps(i))/ &
             log10(profilepressure(i+1)/profilepressure(i))*log10(inpress/profilepressure(i)))
        cont = .false.
     else
        i=i+1
     endif
     
  enddo
!  profileeps = profileeps - 1.0d19
!  outeps = outeps - 1.0d19

end subroutine eps0

subroutine press0(outpress,inrho)

  use betaequil
  implicit none
  
  !inputs
  real*8 inrho
  
  !output
  real*8 outpress

  !internal
  real*8,parameter :: press_cgs_to_gc = 1.80171810d-39
  real*8,parameter :: rho_cgs_to_gc = 1.61930347d-18
  

  integer i
  logical :: cont = .true.

  !first find points that bracket rho
  i = 1
  cont = .true.
  do while(cont)
     if (profilerho(i+1).gt.inrho) then
        outpress = 10.0d0**(log10(profilepressure(i)) + log10(profilepressure(i+1)/profilepressure(i))/ &
             log10(profilerho(i+1)/profilerho(i))*log10(inrho/profilerho(i)))
        cont = .false.
     else
        i=i+1
!        write(*,*) i
     endif
  enddo

!  write(*,*) K,inrho,outpress,K*rho_cgs_to_gc**2.0d0/press_cgs_to_gc
!  stop
  
end subroutine press0

subroutine read_betaequilEOS(filename)

  use betaequil
  implicit none

  character*256 filename
  integer :: plen,iindex,i
  real*8 temp

  open(unit=473,file=trim(adjustl(filename)))

  read(473,*) plen

  allocate(profilerho(plen))
  allocate(profilepressure(plen))
  allocate(profileeps(plen))

  profilerho = 0.0d0
  profileeps = 0.0d0
  profilepressure = 0.0d0
  write(*,*) plen
!  do i=1,plen
!     read(473,*) iindex, temp, profilerho(i), temp,profilepressure(i), profileeps(i)
!  enddo

!  do i=plen,1,-1
  do i=1,plen
     read(473,*) profilerho(i), profileeps(i),profilepressure(i) 
  enddo

  profileeps = profileeps*9.5714d17
  profilepressure = profilepressure*1.60217646d33
  profilerho = profilerho/5.97400341d-16

  

  close(473)

end subroutine read_betaequilEOS
