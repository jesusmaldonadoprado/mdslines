!<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! MSDlines subroutine for determining masses, gravities and radius
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
! last change:  mar mar  3 11:57:42 CET 2015
! jmaldonado at inaf-oapa
!><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! --- Definition of variables
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine msdlgmr
  implicit none

  real,parameter    :: null_value            =  -99.
  real,parameter    :: stephan_boltzmann_cte =  5.67051e-5  ! (erg/cm^2*s*K^4) From Allen
  real,parameter    :: solar_radius          =  6.95508e10  ! (cm) From Allen
  real,parameter    :: solar_luminosity      =  3.845e33    ! (erg/s)  From Allen
  real,parameter    :: pi_number             =  3.141592654
  integer,parameter :: n_coefficients        =  5 
  integer,parameter :: n_outputs             =  14
  real,parameter    :: sigma_mass            =  0.0206811
  real,parameter    :: sigma_radius          =  0.0180707
  real,parameter    :: sigma_logg            =  0.0181498
  character(50),parameter  :: teff_file      =  'derived_temperatures.dat'
  character(50),parameter  :: metal_file     =  'derived_metallicities.dat'
  character(50),parameter  :: sptype_file    =  'derived_sptypes.dat'
  character(50),parameter  :: ofile          =  'msd_summary_stellar_parameters.txt'
  character(12),dimension(:),allocatable    :: star_identifier
  character(12)                             :: borra  
  real,dimension(:),allocatable             :: teff, eteff, metal, emetal, sptype
  real,dimension(:),allocatable             :: mass, emass, radius, eradius, logg, elogg
  real,dimension(:),allocatable             :: luminosity, eluminosity   
  real,dimension(n_coefficients)            :: c_mass, c_radius, c_logg
  integer                                   :: n_stars     
  character(13),dimension(n_outputs)        :: headers
  integer                                   :: ii, jj, kk, ios, jos

! ------------------------------------------------------
!  Open and read necessary files 
! ------------------------------------------------------

   ! >>> Computing n_stars
   n_stars = 0
   open(1,file=teff_file,status='old',action='read')
   read(1,*) ! >>> Headers
   do 
   read(1,*,iostat=ios)
       if (ios/=0) exit
           n_stars = n_stars + 1
   end do 
   close(1)

   ! >>> Allocate vectors
   allocate(star_identifier(n_stars))
   allocate(teff(n_stars))
   allocate(eteff(n_stars))
   allocate(metal(n_stars))
   allocate(emetal(n_stars))
   allocate(sptype(n_stars))
   allocate(mass(n_stars))
   allocate(emass(n_stars))
   allocate(radius(n_stars))
   allocate(eradius(n_stars))
   allocate(logg(n_stars))
   allocate(elogg(n_stars))
   allocate(luminosity(n_stars))
   allocate(eluminosity(n_stars))      
 
   ! >>> Reading file containing teff
   open(1,file=teff_file,status='old',action='read')
   read(1,*) ! >>> Headers
   do ii = 1, n_stars
      read(1,*) star_identifier(ii), teff(ii), eteff(ii)
   end do
   close(1)

   ! >> Reading file containing metallicities
   open(2,file=metal_file,status='old',action='read')
   read(2,*) ! >>> Headers
   do jj = 1, n_stars
      read(2,*) borra, metal(jj), emetal(jj)
   end do 
   close(2)

   ! >> Reading file containing spectral types
   open(3,file=sptype_file,status='old',action='read')
   read(3,*) ! >>> Headers
   do kk = 1, n_stars
      read(3,*) borra, sptype(kk)
   end do 
   close(3)
   
! ------------------------------------------------------
!  Initialization of variables
! ------------------------------------------------------

  ! *** Output file headers
  headers = (/'# Star       ', 'Teff         ', 'eTeff        ', 'SpType       ',  &
              '[Fe/H]       ', 'e[Fe/H]      ', 'Mass         ', 'eMass        ',  &
              'Radius       ', 'eRadius      ', 'logg         ', 'elogg        ',  &
              'log(L*/Lsun) ', 'elog(L*/Lsun)'/)   

  ! *** Stellar Mass:
   c_mass   = (/-171.616, 0.139436, -3.776e-05, 3.41915e-09, 0.382057/)
  ! *** Stellar Radii:
   c_radius = (/-159.857, 0.130206, -3.53437e-05, 3.20786e-09, 0.346746/)
  ! *** Stellar logg:
   c_logg   = (/174.462, -0.137614, 3.72836e-05, -3.37556e-09, -0.332336/)

! ------------------------------------------------------
!  Computing Mass, Radius, logg, and Luminosity
! ------------------------------------------------------

  do ii = 1, n_stars ! For each star 

    ! *** Null value case 
     if (teff(ii).eq.null_value.OR.metal(ii).eq.null_value) then 
         mass(ii)        = null_value
         emass(ii)       = null_value
         radius(ii)      = null_value   
         eradius(ii)     = null_value
         logg(ii)        = null_value
         elogg(ii)       = null_value
         luminosity(ii)  = null_value
         eluminosity(ii) = null_value  
     else

      ! *** Mass
      mass(ii)    = c_mass(1) + c_mass(2)*teff(ii) + c_mass(3)*(teff(ii)**2) + c_mass(4)*(teff(ii)**3)  &
                  +  c_mass(5)*(metal(ii))

      ! *** Uncertainty in mass
      emass(ii)   = sqrt(((c_mass(2) + 2*c_mass(3)*(teff(ii))  + 3*c_mass(4)*(teff(ii)**2))*eteff(ii))**2 &
                  + (c_mass(5)*emetal(ii))**2 + (sigma_mass**2))
 
      ! *** Radii
      radius(ii)  = c_radius(1) + c_radius(2)*teff(ii) + c_radius(3)*(teff(ii)**2) + c_radius(4)*(teff(ii)**3)  &
                  + c_radius(5)*(metal(ii))          

      ! *** Uncertainty in radii
      eradius(ii) = sqrt(((c_radius(2) + 2*c_radius(3)*(teff(ii))  + 3*c_radius(4)*(teff(ii)**2))*eteff(ii))**2 &
                  + (c_radius(5)*emetal(ii))**2 + (sigma_radius**2))
 
      
      ! *** Gravity   
      logg(ii)    = c_logg(1) + c_logg(2)*teff(ii) + c_logg(3)*(teff(ii)**2) + c_logg(4)*(teff(ii)**3)  &
                  + c_logg(5)*(metal(ii))


      ! *** Uncertainty in gravity 
      elogg(ii)   = sqrt(((c_logg(2) + 2*c_logg(3)*(teff(ii)) + 3*c_logg(4)*(teff(ii)**2))*eteff(ii))**2 &
                  + (c_logg(5)*emetal(ii))**2 + (sigma_logg**2))


      ! *** Luminosity in solar units
      luminosity(ii) = 4*pi_number*((radius(ii)*solar_radius)**2)*stephan_boltzmann_cte*(teff(ii)**4)
      luminosity(ii) = luminosity(ii)/solar_luminosity

    
      ! *** Uncertainty in luminosity   
      eluminosity(ii) = sqrt(((2*luminosity(ii)*eradius(ii))/radius(ii))**2 + &
                        ((4*luminosity(ii)*eteff(ii))/teff(ii))**2) 

    
      ! *** Better provide logarithms
      ! >>> Note that errors should be computed before
      eluminosity(ii) = abs((1/(luminosity(ii)*log(10.0)))*eluminosity(ii))
      luminosity(ii)  = log10(luminosity(ii))
 

  
  end if ! End of 'null value' case
  end do ! End of 'for each star ...'

! ------------------------------------------------------
!   Writing results in output file 
! ------------------------------------------------------

  open(4,file=ofile,status='unknown')
  write(4,10) headers(:) ! headers in ofile
  do jj = 1, n_stars
      write(4,20) star_identifier(jj), teff(jj), eteff(jj), sptype(jj), metal(jj), emetal(jj), &
                  mass(jj), emass(jj), radius(jj), eradius(jj), logg(jj), elogg(jj),  &
                  luminosity(jj), eluminosity(jj)
  end do 
  close(4)


 ! >>> Deallocate vectors
   deallocate(star_identifier)
   deallocate(teff)
   deallocate(eteff)
   deallocate(metal)
   deallocate(emetal)
   deallocate(sptype)
   deallocate(mass)
   deallocate(emass)
   deallocate(radius)  
   deallocate(eradius) 
   deallocate(logg) 
   deallocate(elogg) 
   deallocate(luminosity) 
   deallocate(eluminosity) 


! --- Formats and end of program -----------------------------------
10   format(14(A17,1x))
20   format(A12,1x,2(f12.0,1x),1x,f6.1,1x,8(f12.2,1x),1x,2(f12.3,1x))

     end subroutine 

