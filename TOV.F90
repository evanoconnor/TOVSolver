!-*-f90-*-
 program TOV

   use eosmodule
   implicit none

   real*8 out_mass_bary(2) !in grams
   real*8 out_mass,out_bmass !in grams
   real*8 out_radius !in cm
   real*8 out_pradius !in cm
   real*8 dr !step size in cm
   integer N !number of steps
   real*8 central_density !in g/cmm
   real*8 entropy ! in units of k_b, per baryon
   real*8 temperature ! in units of MeV
   real*8,allocatable :: phi(:)
   integer i
   real*8 maxlogrho,minlogrho
   real*8 fixed_ye
   real*8 TOVgamma
   real*8 K1
   logical temp, poly,hybrid,betaequilEOS,newt
   real*8 phi_bound, alpha_c
   real*8,save :: c = 29979245800.0d0
   real*8,save :: G = 6.67384d-8

   character (256) :: arg  
   real*8 density

   real*8 dummy(14),TOVmue,TOVmup,TOVmun,outye
   integer :: keytemp = 1
   integer :: keyerr = 0
   
   call readtable("./path/to/EOS.h5",1)

   !individual TOV star parameters
   dr =500.0d0
   N = 20000
   entropy = 1.0d0
   temperature = 0.1d0
   newt = .false. !for newtonian TOV, only available for constant temp
   temp = .true. !set to false to fix entropy (be warned of EOS errors)
   betaequilEOS = .false. !i.e. P(rho) file
   fixed_ye = -0.05d0 !make negative for beta equilibrium

   !for polytrope
   TOVgamma = 2.0d0
   K1 = 1.0036d13/(2.0d0)**(5.0d0/3.0d0)/1.23924525e9!0.467614341)
   K1 = 30.0d0

   if (betaequilEOS) then
      call read_betaequilEOS("./new_FSU_rho_e_P.dat")
  endif
   
   allocate(phi(N))

  call getarg(1,arg)
  read(arg,"(E18.9)") density

  minlogrho = log10(2.13450d14)
  maxlogrho = log10(3.13450d15)

! if using command line arguement, then change loop  to 1 below
!  minlogrho = log10(density)
!  maxlogrho = log10(density)

  do i = 1,30
     central_density = 10.0d0**(minlogrho)*10.0d0**(real(i)/30.0d0 * &
          (maxlogrho-minlogrho))
     if (betaequilEOS) then
        call mass_betaequilEOS(central_density,out_mass,out_bmass,out_radius,out_pradius,N,dr)
        write(*,"(1P20E18.8)") central_density/1.0d15,out_mass/1.98892d33,out_bmass/1.98892d33,out_radius
     endif
     if (temp) then
        if (newt) then
           call mass_constant_temperature_newtonian(central_density,temperature, &
                out_mass_bary,phi,out_radius,out_pradius,N,dr,fixed_ye)
        else
           call mass_constant_temperature(central_density,temperature, &
                out_mass_bary,phi,out_radius,out_pradius,N,dr,fixed_ye)
        endif
     else
        call mass_constant_entropy(central_density,entropy, &
            out_mass_bary,out_radius,out_pradius,N,dr,fixed_ye)
     endif
     write(*,"(1P20E18.8)") central_density/1.0d15,out_mass_bary(1)/1.98892d33,out_mass_bary(2)/1.98892d33,out_radius/1.0e5
  enddo


end program TOV
