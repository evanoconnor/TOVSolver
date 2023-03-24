!-*-f90-*-
subroutine mass_constant_temperature_newtonian(central_density,temperature, &
     out_mass,out_phi,out_radius,out_pradius,N,dr,fixed_ye)

  use eosmodule
  implicit none

  !inputs
  real*8 central_density
  real*8 temperature !to be made constant throughout star
  real*8 fixed_ye
  real*8 out_mass(2), out_radius, out_pradius
  real*8 dr !step size
  integer N !maximum number of steps
  real*8 out_phi(N)

  !outputables
  real*8,allocatable :: TOVentropy(:)
  real*8,allocatable :: TOVye(:)
  real*8,allocatable :: TOVpressure(:)
  real*8,allocatable :: TOVrho(:)
  real*8,allocatable :: TOVphi(:)
  real*8,allocatable :: TOVrad(:)
  real*8,allocatable :: TOVmass(:)
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
  integer yeflag
  real*8 :: G = 6.6742d-8

  !RK variables
  real*8 k1p,k2p,k3p,k4p
  real*8 k1m,k2m,k3m,k4m
  real*8 k1ph,k2ph,k3ph,k4ph
  
  real*8 x
  
  integer flag_return

  allocate(TOVentropy(N))
  allocate(TOVye(N))
  allocate(TOVpressure(N))
  allocate(TOVrho(N))
  allocate(TOVphi(N))
  allocate(TOVrad(N))
  allocate(TOVmass(N))
  allocate(TOVcs2(N))
  allocate(TOVeps(N))
  allocate(TOVtemperature(N))

  TOVentropy(:) = 0.0d0
  TOVye(:) = 0.0d0
  TOVpressure(:) = 0.0d0
  TOVrho(:) = 0.0d0
  TOVphi(:) = 0.0d0
  TOVrad(:) = 0.0d0
  TOVmass(:) = 0.0d0
  TOVcs2(:) = 0.0d0
  TOVeps(:) = 0.0d0
  TOVtemperature(:) = 0.0d0

  keytemp = 1
  keyerr = 0
  


  !first get center cell
  !find ye given rho and temp, this is call to ott_eos with keytemp=2
  yeflag=0
  call get_ye_temperature(central_density,temperature,TOVye(1),fixed_ye,yeflag)
  TOVrho(1) = central_density
  TOVrad(1) = 0.0d0
  TOVmass(1) = 0.0d0
  TOVphi(1) = 0.0d0
  TOVtemperature(1) = temperature
  i=1
  call nuc_eos_full(TOVrho(1),temperature,TOVye(1),TOVeps(1), &
       TOVpressure(1),TOVentropy(1),TOVcs2(1),dummy(5),dummy(6), &
       dummy(7),dummy(8),dummy(9),dummy(10),dummy(11),dummy(12), &
       dummy(13),dummy(14),dummy(15),dummy(16),dummy(17),keytemp, &
       keyerr,precision)
  if(keyerr.ne.0) then
     write(*,*) "EOS error"
     stop "EOS error"
  endif

  open(unit=473,file="neutron_star_structure.dat")

  write(473,"(i8,1P10E18.9)") i,TOVmass(i),TOVrad(i),TOVtemperature(i)*temp_mev_to_kelvin,TOVrho(i),0.0d0,TOVye(i),0.0d0
 
  !now do RK4 to get new values of press, mass and phi
  flag_return = 0
  do i=2,N
     call dpdr_newtonian(TOVmass(i-1),TOVpressure(i-1),TOVrho(i-1),TOVrad(i-1),k1p)
     call dmdr_newtonian(TOVrho(i-1),TOVrad(i-1),k1m)
     call dphidr_newtonian(TOVmass(i-1),TOVpressure(i-1),TOVrho(i-1),TOVrad(i-1),k1ph)

     if (TOVpressure(i-1)+dr*0.5d0*k1p.lt.0.0d0) then
        out_mass(1) = TOVmass(i-1)
        out_phi(:) = TOVphi(:)
        out_radius = TOVrad(i-1)
        return
     endif
    
     call get_rho_and_ye_temperature(TOVpressure(i-1),0.5d0*dr*k1p,temperature,temp_rho, &
          temp_ye,temp_eps,TOVrho(i-1),TOVye(i-1),TOVcs2(i-1),fixed_ye,flag_return,yeflag)
     if (flag_return.eq.1) then
        write(*,*) "here1",i
        out_mass(1) = TOVmass(i-1)
        out_phi(:) = TOVphi(:)
        out_radius = TOVrad(i-1)
        return
     endif

     if (temp_rho.lt.2.0d6) then
        out_mass(1) = TOVmass(i-1)
        out_radius = TOVrad(i-1)
        out_phi(:) = TOVphi(:)
        return
     endif

     call dpdr_newtonian(TOVmass(i-1)+0.5d0*dr*k1m,TOVpressure(i-1)+0.5d0*dr*k1p, &
          temp_rho,TOVrad(i-1)+0.5d0*dr,k2p)
     call dmdr_newtonian(temp_rho,TOVrad(i-1)+0.5d0*dr,k2m)
     call dphidr_newtonian(TOVmass(i-1)+0.5d0*dr*k1m,TOVpressure(i-1)+0.5d0*dr*k1p, &
          temp_rho,TOVrad(i-1)+0.5d0*dr,k2ph)

     if (TOVpressure(i-1)+dr*0.5d0*k2p.lt.0.0d0) then
        out_mass(1) = TOVmass(i-1)
        out_phi(:) = TOVphi(:)
        out_radius = TOVrad(i-1)
        return
     endif

     call get_rho_and_ye_temperature(TOVpressure(i-1),dr*0.5d0*k2p, &
          temperature,temp_rho,temp_ye,temp_eps,TOVrho(i-1), &
          TOVye(i-1),TOVcs2(i-1),fixed_ye,flag_return,yeflag)
     if (flag_return.eq.1) then
        write(*,*) "here2", i
        out_mass(1) = TOVmass(i-1)
        out_radius = TOVrad(i-1)
        out_phi(:) = TOVphi(:)
        return
     endif

     if (temp_rho.lt.2.0d6) then
        out_mass(1) = TOVmass(i-1)
        out_radius = TOVrad(i-1)
        out_phi(:) = TOVphi(:)
        return
     endif

     call dpdr_newtonian(TOVmass(i-1)+dr*0.5d0*k2m,TOVpressure(i-1)+dr*0.5d0*k2p, &
          temp_rho,TOVrad(i-1)+0.5d0*dr,k3p)
     call dmdr_newtonian(temp_rho,TOVrad(i-1)+0.5d0*dr,k3m)
     call dphidr_newtonian(TOVmass(i-1)+dr*0.5d0*k2m,TOVpressure(i-1)+dr*0.5d0*k2p, &
          temp_rho,TOVrad(i-1)+0.5d0*dr,k3ph)

     if (TOVpressure(i-1)+dr*k3p.lt.0.0d0) then
        out_mass(1) = TOVmass(i-1)
        out_radius = TOVrad(i-1)
        out_phi(:) = TOVphi(:)
        return
     endif
     
     call get_rho_and_ye_temperature(TOVpressure(i-1),dr*k3p,temperature, &
          temp_rho,temp_ye,temp_eps,TOVrho(i-1),TOVye(i-1),TOVcs2(i-1), &
          fixed_ye,flag_return,yeflag)
     if (flag_return.eq.1) then
        write(*,*) "here3", i
        out_mass(1) = TOVmass(i-1)
        out_radius = TOVrad(i-1)
        out_phi(:) = TOVphi(:)
        return
     endif

     if (temp_rho.lt.2.0d6) then
        out_mass(1) = TOVmass(i-1)
        out_radius = TOVrad(i-1)
        out_phi(:) = TOVphi(:)
        return
     endif

     call dpdr_newtonian(TOVmass(i-1)+dr*k3m,TOVpressure(i-1)+dr*k3p, &
          temp_rho,TOVrad(i-1)+dr,k4p)
     call dmdr_newtonian(temp_rho,TOVrad(i-1)+dr,k4m)
     call dphidr_newtonian(TOVmass(i-1)+dr*k3m,TOVpressure(i-1)+dr*k3p, &
          temp_rho,TOVrad(i-1)+dr,k4ph)

     TOVpressure(i) = TOVpressure(i-1) + dr/6.0d0*(k1p + 2.0d0*k2p + &
          2.0d0*k3p + k4p)
     TOVmass(i) = TOVmass(i-1) + dr/6.0d0*(k1m + 2.0d0*k2m + &
          2.0d0*k3m + k4m)
     TOVrad(i) = TOVrad(i-1)+dr
     TOVphi(i) = TOVphi(i-1) + dr/6.0d0*(k1ph + 2.0d0*k2ph + &
          2.0d0*k3ph + k4ph)
     TOVtemperature(i) = temperature

     call get_rho_and_ye_temperature(TOVpressure(i-1), &
          dr/6.0d0*(k1p + 2.0d0*k2p + 2.0d0*k3p + k4p), &
          temperature,TOVrho(i),TOVye(i),TOVeps(i), &
          TOVrho(i-1),TOVye(i-1),TOVcs2(i-1),fixed_ye,flag_return,yeflag)
     if (flag_return.eq.1) then
        write(*,*) "here4", i
        out_mass(1) = TOVmass(i-1)
         out_radius = TOVrad(i-1)
        out_phi(:) = TOVphi(:)
        return
     endif

      if (TOVrho(i).lt.1.0d6) then
        out_mass(1) = TOVmass(i)
        out_radius = TOVrad(i)
        out_phi(:) = TOVphi(:)
        return
     endif

     call nuc_eos_full(TOVrho(i),temperature,TOVye(i),TOVeps(i), &
          dummy(4),TOVentropy(i),TOVcs2(i),dummy(5),dummy(6), &
          dummy(7),dummy(8),dummy(9),dummy(10),dummy(11),dummy(12), &
          dummy(13),dummy(14),dummy(15),dummy(16),dummy(17),keytemp, &
          keyerr,precision)
     if(keyerr.ne.0) then
        write(*,*) "EOS error"
        stop "EOS error"
     endif

     if (abs(dummy(4)/TOVpressure(i)-1.0d0).gt.1.0d-7) then
        write(*,*) "pressure wrong", dummy(4),TOVpressure(i), &
             dummy(4)/TOVpressure(i)
        stop
     endif

     if (TOVpressure(i).lt.0.0d0) then
        out_mass(1) = TOVmass(i-1)
        out_radius = TOVrad(i-1)
        out_phi = TOVphi
        return
     endif

     if (TOVpressure(i).lt.1.6021765d25) then
        out_mass(1) = TOVmass(i-1)
        out_radius = TOVrad(i-1)
        out_phi = TOVphi
        return
     endif

     write(473,"(i8,1P10E18.9)") i,TOVmass(i),TOVrad(i),TOVtemperature(i)*temp_mev_to_kelvin,TOVrho(i),0.0d0,TOVye(i),0.0d0


  enddo

  close(473)

  out_mass(1) = TOVmass(N-1)
  out_radius = real(N-1)*dr
  out_phi = TOVphi
!  stop

end subroutine mass_constant_temperature_newtonian

subroutine get_ye_temperature(density,temperature,xye,type,yeflag)

  use eosmodule
  implicit none
  
  !inputs
  real*8 density,temperature,type
  
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
  integer yeflag
  
  keytemp = 1
  keyerr = 0
  
  ye_upper = 0.5d0
  ye_lower = eos_yemin
  tol = 1.0d-9
  if (type.gt.0.0d0) then
     xye=type
     return
  endif
  
  counter = 0
  cont = .true.
  do while(cont)
     counter = counter + 1
     ye_guess = (ye_upper + ye_lower)/2.0d0
     call nuc_eos_full(density,temperature,ye_guess,dummy(2),dummy(3), &
          dummy(1),dummy(4),dummy(5),dummy(6),dummy(7),dummy(8), &
          dummy(9),dummy(10),dummy(11),dummy(12),dummy(13),TOVmue,TOVmun, &
          TOVmup,dummy(14),keytemp,keyerr,precision)
     if(keyerr.ne.0) then
        write(*,*) "EOS error"
        stop "EOS error"
     endif
     !mue includes rest mass, mun,mup do not, yes they do!
     TOVmunu = (TOVmue)-TOVmun+TOVmup-1.295d0
     
     if (TOVmunu.gt.0.0d0) then
        ye_upper = ye_guess
     else
        ye_lower = ye_guess
     endif
     
     if (ye_upper-ye_lower.lt.tol) then
        cont = .false.
     endif
   
  enddo

xye = ye_guess

end subroutine get_ye_temperature

subroutine get_rho_and_ye_temperature(xpressure,dpress,temperature, &
     newrho,newye,neweps,oldrho,oldye,oldcs2,type,flag_return,yeflag)

  use eosmodule

  implicit none

  !inputs
  real*8 xpressure,temperature,dpress
  real*8 oldrho,oldye,oldcs2
  real*8 type
  integer yeflag

  !outputs
  real*8 newrho,newye,neweps

  integer flag_return

  !internal
  real*8 ye_guess,rho_guess,press_guess,eps_guess
  real*8 rho_guess2,press_guess2
  real*8 dummy(17)
  integer keytemp,keyerr
  logical cont
  integer counter
  real*8 newcs2
  real*8 TOVdpdrhoe
  real*8 press_compare
  real*8 fac
  real*8 mydpdrho

  keytemp = 3 !fix pressure, find density
  keyerr = 0

  !desired pressure
  press_compare = xpressure + dpress

  ye_guess = oldye
  rho_guess = oldrho

  counter = 0
  cont = .true.
  fac = 1.0d0
  
  !first get ye (if need beta-equilibrium), don't iterate this with press/rho
  call get_ye_temperature(rho_guess,temperature,ye_guess,type,yeflag)

  !call EOS to invert press to density.
  call nuc_eos_full(rho_guess,temperature,ye_guess,eps_guess,press_compare, &
       dummy(1),newcs2,dummy(5),dummy(6),TOVdpdrhoe,dummy(8), &
       dummy(9),dummy(10),dummy(11),dummy(12),dummy(13),dummy(14), &
       dummy(15),dummy(16),dummy(17),keytemp,keyerr,precision)
  if(keyerr.ne.0) then
     if (keyerr.eq.473) then
        write(*,*) "warning, dpdrho very small.  Hopefully close to edge of star"
        flag_return = 1
        return        
     endif
     if (rho_guess.lt.1.0d5) then
        write(*,*) "warning, rho<rho_min.  Hopefully close to edge of star"
        flag_return = 1
        return
     endif
     write(*,*) "EOS error"
     stop "EOS error"
  endif


  if (rho_guess.lt.eos_rhomin) then
     write(*,*) "warning, rho<rho_min.  Hopefully close to edge of star"
     flag_return = 1
     return
  endif

  newye = ye_guess
  newrho = rho_guess
  neweps = eps_guess

end subroutine get_rho_and_ye_temperature
