!-*-f90-*-
subroutine mass_constant_entropy(central_density,entropy,out_mass, &
     out_radius,out_pradius,N,dr,fixed_ye)

  use eosmodule
  implicit none

  !inputs
  real*8 central_density
  real*8 fixed_ye
  real*8 entropy !to be made constant throughout star
  real*8 out_mass(2), out_radius, out_pradius
  real*8 dr !step size
  integer N !maximum number of steps

  !outputables
  real*8,allocatable :: TOVentropy(:)
  real*8,allocatable :: TOVye(:)
  real*8,allocatable :: TOVpressure(:)
  real*8,allocatable :: TOVrho(:)
  real*8,allocatable :: TOVphi(:)
  real*8,allocatable :: TOVrad(:)
  real*8,allocatable :: TOVprad(:)
  real*8,allocatable :: TOVmass(:)
  real*8,allocatable :: TOVbaryonmass(:)
  real*8,allocatable :: TOVcs2(:)
  real*8,allocatable :: TOVeps(:)
  real*8,allocatable :: TOVtemperature(:)

  !internal
  integer i
  real*8 dummy(17)
  real*8 ye_guess,ye_lower,ye_upper
  real*8 TOVmue,TOVmun,TOVmup,TOVmunu
  real*8 temp_rho,temp_ye,temp_eps
  integer keytemp,keyerr
  logical cont
  real*8 h
  real*8 :: c = 29979245800.0d0
  real*8 :: G = 6.6742d-8

  !RK variables
  real*8 k1p,k2p,k3p,k4p
  real*8 k1m,k2m,k3m,k4m
  real*8 k1ph,k2ph,k3ph,k4ph
  real*8 k1mb,k2mb,k3mb,k4mb
  real*8 rhocgs2nbfm,rho_nuc_cgs,nb,my_e


  allocate(TOVentropy(N))
  allocate(TOVye(N))
  allocate(TOVpressure(N))
  allocate(TOVrho(N))
  allocate(TOVphi(N))
  allocate(TOVrad(N))
  allocate(TOVprad(N))
  allocate(TOVmass(N))
  allocate(TOVbaryonmass(N))
  allocate(TOVcs2(N))
  allocate(TOVeps(N))
  allocate(TOVtemperature(N))

  TOVentropy(:) = 0.0d0
  TOVye(:) = 0.0d0
  TOVpressure(:) = 0.0d0
  TOVrho(:) = 0.0d0
  TOVphi(:) = 0.0d0
  TOVrad(:) = 0.0d0
  TOVprad(:) = 0.0d0
  TOVmass(:) = 0.0d0
  TOVcs2(:) = 0.0d0
  TOVeps(:) = 0.0d0
  TOVtemperature(:) = 0.8d0

  keytemp = 2 !constant entropy
  keyerr = 0


  !first get center cell
  !find ye given rho and temp, this is call to ott_eos with keytemp=2
  call get_ye_entropy(central_density,entropy,TOVye(1),fixed_ye)
  TOVrho(1) = central_density
  TOVrad(1) = 0.0d0
  TOVprad(1) = 0.0d0
  TOVmass(1) = 0.0d0
  TOVbaryonmass(1) = 0.0d0
  TOVphi(1) = 0.0d0
  TOVentropy(1) = entropy
  call nuc_eos_full(TOVrho(1),TOVtemperature(1),TOVye(1),TOVeps(1), &
       TOVpressure(1),entropy,TOVcs2(1),dummy(5),dummy(6), &
       dummy(7),dummy(8),dummy(9),dummy(10),dummy(11),dummy(12), &
       dummy(13),dummy(14),dummy(15),dummy(16),dummy(17),keytemp, &
       keyerr,precision)

  i=1
  write(*,"(i6,1P10E18.9)") i,TOVmass(i),TOVrad(i), &
       TOVtemperature(i)*1.1604447522806d10,TOVrho(i),0.0d0,TOVye(i),0.0d0  
  if(keyerr.ne.0) then
     write(*,*) "EOS error"
     stop "EOS error"
  endif
  
  !next cell
  TOVye(2) = TOVye(1)
  TOVpressure(2) = TOVpressure(1)
  rhocgs2nbfm = 1.0d-39/1.66d-24
  rho_nuc_cgs = 1.66d14
  nb = rhocgs2nbfm*central_density
  my_e = central_density*(1.0d0+TOVeps(1)/c**2)
  h = log10((my_e+TOVpressure(1)/c**2)/(10.0d0*nb*rho_nuc_cgs))

  TOVbaryonmass(2) = 0.0d0
  TOVphi(2) = 0.0d0
  TOVeps(2) = TOVeps(1)
  TOVentropy(2) = entropy
  TOVtemperature(2) = TOVtemperature(1)
  TOVrho(2) = TOVrho(1)
  TOVrad(2) = dr
  TOVmass(2) = TOVrho(1)*h*4.0d0*pi/3.0d0*TOVrad(2)**3
  TOVprad(2) = dr
  TOVcs2(2) = TOVcs2(1)
  i=2
  write(*,"(i6,1P10E18.9)") i,TOVmass(i),TOVrad(i), &
       TOVtemperature(i)*1.1604447522806d10,TOVrho(i),0.0d0,TOVye(i),0.0d0  
  !now do RK4 to get new values of press, mass and phi
  do i=3,N
     h = 1.0d0+TOVeps(i-1)/c**2+TOVpressure(i-1)/TOVrho(i-1)/c**2
     call dpdr(TOVmass(i-1),TOVpressure(i-1),TOVrho(i-1)*h,TOVrad(i-1),k1p)
     call dmdr(TOVpressure(i-1),TOVrho(i-1)*h,TOVrad(i-1),k1m)
     call dphidr(TOVmass(i-1),TOVpressure(i-1),TOVrho(i-1)*h,TOVrad(i-1),k1ph)
     call dmbdr(TOVmass(i-1),TOVrho(i-1),TOVrad(i-1),k1mb)
    
     if (TOVpressure(i-1)+dr*k1p.lt.0.0d0) then
        out_mass(1) = TOVmass(i-1)
        out_mass(2) = TOVbaryonmass(i-1)
        out_radius = TOVrad(i-1)
        return
     endif
    
     call get_rho_and_ye_entropy(TOVpressure(i-1),0.5d0*dr*k1p,entropy,temp_rho, &
          temp_ye,temp_eps,TOVrho(i-1),TOVye(i-1),TOVcs2(i-1),fixed_ye)

     h = 1.0d0+temp_eps/c**2+(TOVpressure(i-1)+0.5d0*dr*k1p)/temp_rho/c**2
     call dpdr(TOVmass(i-1)+0.5d0*dr*k1m,TOVpressure(i-1)+0.5d0*dr*k1p, &
          temp_rho*h,TOVrad(i-1)+0.5d0*dr,k2p)
     call dmdr(TOVpressure(i-1)+0.5d0*dr*k1p,temp_rho*h,TOVrad(i-1)+0.5d0*dr,k2m)
     call dphidr(TOVmass(i-1)+0.5d0*dr*k1m,TOVpressure(i-1)+0.5d0*dr*k1p, &
          temp_rho*h,TOVrad(i-1)+0.5d0*dr,k2ph)
     call dmbdr(TOVmass(i-1)+0.5d0*dr*k1m,temp_rho,TOVrad(i-1)+0.5d0*dr,k2mb)

     if (TOVpressure(i-1)+dr*0.5d0*k2p.lt.0.0d0) then
        out_mass(1) = TOVmass(i-1)
        out_mass(2) = TOVbaryonmass(i-1)
        out_radius = TOVrad(i-1)
        return
     endif

     call get_rho_and_ye_entropy(TOVpressure(i-1),dr*0.5d0*k2p,entropy,temp_rho, &
          temp_ye,temp_eps,TOVrho(i-1),TOVye(i-1),TOVcs2(i-1),fixed_ye)

     h = 1.0d0+temp_eps/c**2+(TOVpressure(i-1)+dr*0.5d0*k2p)/temp_rho/c**2
     call dpdr(TOVmass(i-1)+dr*0.5d0*k2m,TOVpressure(i-1)+dr*0.5d0*k2p, &
          temp_rho*h,TOVrad(i-1)+0.5d0*dr,k3p)
     call dmdr(TOVpressure(i-1)+dr*0.5d0*k2p,temp_rho*h,TOVrad(i-1)+0.5d0*dr,k3m)
     call dphidr(TOVmass(i-1)+dr*0.5d0*k2m,TOVpressure(i-1)+dr*0.5d0*k2p, &
          temp_rho*h,TOVrad(i-1)+0.5d0*dr,k3ph)
     call dmbdr(TOVmass(i-1)+0.5d0*dr*k2m,temp_rho,TOVrad(i-1)+0.5d0*dr,k3mb)

     if (TOVpressure(i-1)+dr*0.5d0*k3p.lt.0.0d0) then
        out_mass(1) = TOVmass(i-1)
        out_mass(2) = TOVbaryonmass(i-1)
        out_radius = TOVrad(i-1)
        return
     endif

     call get_rho_and_ye_entropy(TOVpressure(i-1),dr*k3p,entropy,temp_rho,temp_ye, &
          temp_eps,TOVrho(i-1),TOVye(i-1),TOVcs2(i-1),fixed_ye)

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
     TOVentropy(i) = entropy
     TOVbaryonmass(i) = TOVbaryonmass(i-1) + dr/6.0d0*(k1mb + 2.0d0*k2mb + &
          2.0d0*k3mb + k4mb)

     call get_rho_and_ye_entropy(TOVpressure(i-1), &
          dr/6.0d0*(k1p + 2.0d0*k2p + 2.0d0*k3p + k4p), &
          entropy,TOVrho(i),TOVye(i),TOVeps(i), &
          TOVrho(i-1),TOVye(i-1),TOVcs2(i-1),fixed_ye)

     if (TOVrho(i).lt.1.0d2) then
        out_mass(1) = TOVmass(i)
        out_mass(2) = TOVbaryonmass(i)
        out_radius = TOVrad(i)
        out_pradius = TOVprad(i)
        return
     endif

     call nuc_eos_full(TOVrho(i),TOVtemperature(i),TOVye(i),TOVeps(i), &
          dummy(4),entropy,TOVcs2(i),dummy(5),dummy(6), &
          dummy(7),dummy(8),dummy(9),dummy(10),dummy(11),dummy(12), &
          dummy(13),dummy(14),dummy(15),dummy(16),dummy(17),keytemp, &
          keyerr,precision)
     if(keyerr.ne.0) then
        write(*,*) "EOS error"
        stop "EOS error"
     endif

     if (abs(dummy(4)/TOVpressure(i)-1.0d0).gt.1.0d-7) then
        write(*,*) "pressure wrong", dummy(4),TOVpressure(i), &
             dummy(4)/TOVpressure(i), i,TOVmass(i),TOVrad(i), TOVrho(i), &
             TOVye(i),TOVtemperature(i)
        stop
     endif

!     if (TOVpressure(i).lt.1.0d-6*TOVpressure(1)) then
     if (TOVpressure(i).lt.0.0d0) then
        out_mass(1) = TOVmass(i-1)
        out_mass(2) = TOVbaryonmass(i-1)
        out_radius = TOVrad(i-1)
        out_pradius = TOVprad(i-1)
        return
     endif

!     if (mod(i,50).eq.0) then
!     write(*,"(i6,1P10E18.9)") i,TOVmass(i)/1.9889e33,TOVrad(i), &
!          TOVtemperature(i)*1.1604447522806d10,TOVrho(i),0.0d0,TOVye(i),0.0d0
!     endif

  enddo

  stop
  
  out_mass(1) = TOVmass(N-1)
  out_mass(2) = TOVbaryonmass(N-1)
  out_radius = N*dr

end subroutine mass_constant_entropy

subroutine get_ye_entropy(density,entropy,xye,type)

  use eosmodule
  implicit none

  !inputs
  real*8 density,entropy,type

  !outputs
  real*8 xye
  
  !internal
  integer keytemp,keyerr
  real*8 ye_guess, ye_upper, ye_lower
  real*8 dummy(14)
  real*8 TOVmue,TOVmun,TOVmup,TOVmunu
  logical cont
  real*8 tol
  integer counter

  keytemp = 2
  keyerr = 0

  ye_upper = 0.5d0
  ye_lower = 0.01d0
  tol = 1.0d-9

  if (type.gt.0.0d0) then
     xye = type
     return
  endif

  counter = 0
  cont = .true.
  do while(cont)
     counter = counter + 1
     ye_guess = (ye_upper + ye_lower)/2.0d0
     call nuc_eos_full(density,dummy(1),ye_guess,dummy(2),dummy(3), &
          entropy,dummy(4),dummy(5),dummy(6),dummy(7),dummy(8), &
          dummy(9),dummy(10),dummy(11),dummy(12),dummy(13),TOVmue,TOVmun, &
          TOVmup,dummy(14),keytemp,keyerr,precision)
     if(keyerr.ne.0) then
        write(*,*) "EOS error"
        stop "EOS error"
     endif
     TOVmunu = TOVmue-TOVmun+TOVmup

     if (TOVmunu.gt.0.0d0) then
        ye_upper = ye_guess
     else
        ye_lower = ye_guess
     endif

     if (ye_upper-ye_lower.lt.tol) then
        cont = .false.
     endif
     
     if (ye_guess.lt.0.001d0) then
        cont = .false.
     endif

  enddo

  xye = ye_guess

end subroutine get_ye_entropy

subroutine get_rho_and_ye_entropy(xpressure,dpress,entropy,newrho, &
     newye,neweps,oldrho,oldye,oldcs2,type)

  use eosmodule

  implicit none

  !inputs
  real*8 xpressure,entropy,dpress
  real*8 oldrho,oldye,oldcs2
  real*8 type

  !outputs
  real*8 newrho,newye,neweps

  !internal
  real*8 ye_guess,rho_guess,press_guess,eps_guess
  real*8 dummy(17)
  integer keytemp,keyerr
  logical cont
  integer counter
  real*8 newcs2
  real*8 TOVdpdrhoe
  real*8 press_compare

  keytemp = 2
  keyerr = 0

  press_compare = xpressure + dpress

  ye_guess = oldye
  rho_guess = oldrho + dpress/oldcs2

  counter = 0
  cont = .true.
  do while(cont)
     counter = counter + 1
     !get ye based on rho and entropy (if needing beta-equilibrium)
     call get_ye_entropy(rho_guess,entropy,ye_guess,type)
     !need to invert by hand since need entropy and pressure constant.
     call nuc_eos_full(rho_guess,dummy(1),ye_guess,eps_guess,press_guess, &
          entropy,newcs2,dummy(5),dummy(6),TOVdpdrhoe,dummy(8), &
          dummy(9),dummy(10),dummy(11),dummy(12),dummy(13),dummy(14), &
          dummy(15),dummy(16),dummy(17),keytemp,keyerr,precision)
     if(keyerr.ne.0) then
        write(*,*) "EOS error"
        stop "EOS error"
     endif

     if (abs(1.0d0-press_guess/press_compare).lt.1.0d-8) then
        cont = .false.
     else
        rho_guess = rho_guess + 0.5d0*(press_compare-press_guess)/newcs2
     endif

  enddo

  newye = ye_guess
  newrho = rho_guess
  neweps = eps_guess

end subroutine get_rho_and_ye_entropy
