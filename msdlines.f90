!<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! MAIN DRIVER FOR MSDLINES: It reads the spectra, computes the
! pseudo EWs, and calls the subroutines for stellar parameters
! computation  
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
! last change: dom mar  1 11:38:34 CET 2015 
! jmaldonado at inaf-oapa
!><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  program msdlines
  implicit none

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! --- Definition of variables
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  character(80),parameter                   :: input_file = 'stars_list.txt'
  character(5)                              :: input_type 
  character(12)                             :: star_identifier 
  character(80)                             :: stellar_spectra
  double precision,dimension(:),allocatable :: wavelength, flux
  integer                                   :: n_wavelengths
  integer                                   :: im, ier, readonly, blocksize
  character(80)                             :: comment, keyword
  integer                                   :: naxis, naxis1
  double precision                          :: crval
  real                                      :: cdelt
  integer                                   :: group, fpixel
  real                                      :: nulv
  logical                                   :: anyf
  integer                                   :: nbuffer
  real,dimension(:),allocatable             :: x
  integer                                   :: ii, jj, kk, ios, jos
  integer                                   :: longitud


!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! Reading input spectra
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  ! >>> Show help text
 
  write(*,*) '----------------------------------------------'
  write(*,*) ' Pseudo equivalent widths computation         '
  write(*,*) '----------------------------------------------'
 
  ! >>> Master file containing all the spectra to be analyzed
  open(3,file=input_file,status='old',action='read') 
  read(3,*) ! >> Line containing the headers
  read(3,*) input_type   !>> Input type should be either ascii or fits 
  read(3,*) ! >> More headers
  do
  read(3,*,iostat=ios) star_identifier, stellar_spectra
  if (ios/=0) exit

  ! *** Input spectra in ASCII format
  !----------------------------------------
  if (input_type.eq.'ascii') then 

       ! >>> Reading input spectra once to compute number of wavelengths
       n_wavelengths = 0
       open(1,file=stellar_spectra,status='old',action='read')
       do
       read(1,*,iostat=jos)
           if (jos/=0) exit
               n_wavelengths = n_wavelengths + 1
       end do
       close(1)

       ! >>> Allocate vectors
       allocate(wavelength(n_wavelengths))
       allocate(flux(n_wavelengths))

       ! >>> Reading wavelenghts and fluxes
       open(1,file=stellar_spectra,status='old',action='read')
       do ii = 1, n_wavelengths
          read(1,*) wavelength(ii), flux(ii)
       end do
       close(1) 

  ! *** Input spectra in FITS format
  !---------------------------------------- 

   else if (input_type.eq.'fits') then

       !*** Initialization of variables / Reading the fits file 
       readonly = 0
       call ftgiou(im,ier)
       call ftopen (im,stellar_spectra,readonly,blocksize,ier)

       !*** Checking if fits is one-dimensional
        call ftgkyj(im,'NAXIS',naxis,comment,ier)
             if (naxis.eq.0.OR.naxis.gt.1) then
                 write(*,*) '>> ERROR!: Multidimensionnal spectra not yet implemented ...'
                 write(*,*) 'Exiting now ...'
                stop
             end if
 
       !*** Reading headers data to reconstruct the wavelengths
       call ftgkyj(im,'NAXIS1',naxis1,comment,ier)
       call ftgkys(im,'CRVAL1',keyword,comment,ier)
       read( keyword, '(f24.12)' ) crval
       call ftgkys(im,'CDELT1',keyword,comment,ier)
       read( keyword, '(f12.8)' ) cdelt

       !*** Reading the fluxes
        nulv    = -1.0
        fpixel  = 1

       allocate(x(naxis1))
       call ftgpve(im,group,fpixel,naxis1,nulv,x,anyf,ier)
   
       !*** Writing fluxes and wavelengths   
       allocate(wavelength(naxis1))
       allocate(flux(naxis1))
        
       n_wavelengths = naxis1 

       do ii = 1, naxis1
          flux(ii) = x(ii)
       end do

       wavelength(1) = crval
       do ii = 2, naxis1
          wavelength(ii) = wavelength(ii-1) + cdelt
       end do

      !*** closing files 
      deallocate(x) 
      call ftclos(im,ier)
      call ftfiou(im,ier)
  
  ! *** 'Unknown' input spectra format
  !----------------------------------------
  else
     write(*,*) '>> ERROR!: Format of input spectra not recognized. Exiting now ...'
                 stop 
  end if 

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! Calling the subroutine which computes the pseudo EWs
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  call pseudoEW (star_identifier,n_wavelengths,wavelength,flux)

  ! >>> Deallocate wavelength and flux
  deallocate(wavelength)
  deallocate(flux)

  end do ! >>> End of 'For each star in the input file ...'
  close(3)  

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! 2nd part of the program: Stellar parameters computation
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  ! >>> Show help text

  write(*,*) '----------------------------------------------'
  write(*,*) ' Stellar parameters computation               '
  write(*,*) '----------------------------------------------'

 ! >>> Calling the subroutine for effective temperature computation
  call  msdteff 
 
 ! >>> Calling the subroutine for metallicity computation
  call  msdmetal

 ! >>> Calling the subroutine for spectral type determination
  call  msdstypes

 ! >>> Calling the subroutine for luminisoty, mass, gravity
  call msdlgmr

  write(*,*) '-----------------------------------------------'
  write(*,*) ' DONE! See msd_summary_stellar_parameters file '
  write(*,*) '-----------------------------------------------'


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   end program

!====================================================================
! Functions and subroutines
!====================================================================

! -------------------------------------------------------------------
! Subroutine for computing the pseudo EWs
! -------------------------------------------------------------------

  subroutine pseudoEW(star_id,n_wave,wave,flx)
  implicit none

  real,parameter                          :: tol             =  5e-3
  real,parameter                          :: nwave_step      =  0.01
  integer,parameter                       :: n_features      =  4424
  character(80),parameter                 :: features_file   =  'harpn_all_4424_identified_features.dat'
  real,parameter                          :: null_ew         =  -99.0
  real,parameter                          :: min_ew          =    8.0
  real,parameter                          :: max_ew          =  120.0
  character(50)                           :: output_file
  character(12)                           :: star_id
  integer                                 :: n_wave
  double precision,dimension(n_wave)      :: wave, flx 
  integer                                 :: length
  integer,dimension(3)                    :: today, now
  character(5),dimension(n_features)      :: feature_id
  real,dimension(n_features)              :: feature_center, feature_width
  real,dimension(n_features)              :: feature_wave_low, feature_wave_upp
  real,dimension(n_features)              :: feature_flux_low, feature_flux_upp
  real,dimension(n_features)              :: ew, e_ew, slope, ordinate
  real,dimension(n_features)              :: n_points, mean_flux, sigma_flux
  real,dimension(n_features)              :: mean_ptp, mean_wave_step
  real                                    :: wave_step 
  integer                                 :: ii, jj

  ! >>> Name and headers in the output file
  length       = len_trim(star_id)
  output_file  = star_id(1:length)//'_psdEWs.dat'

  open(2,file=output_file,status='unknown')
  call idate(today)
  call itime(now)
  write(2,30) today(2), today(1), today(3), now
  write(2,35) '# Pseudo EW of selected features (mA) for star: ', star_id
  write(2,45) '# Feature ', 'Wavelength', 'pseudoEW  ', 'error    '

  ! >>> Show some info in the screen 
  write(*,50) ' Computing pseudo EWs for star : ', star_id

  ! >>> Reading the file containing the features information

   open(5,file=features_file,status='old',action='read')
   read(5,*)
      do ii = 1, n_features
         read(5,*) feature_id(ii), feature_center(ii), feature_width(ii)
         feature_wave_low(ii) = feature_center(ii) - (feature_width(ii)/2.0)
         feature_wave_upp(ii) = feature_center(ii) + (feature_width(ii)/2.0)
      end do
   close(5)


  !>>> For each line, computes the pseudo equivalent width 
  do ii = 1, n_features
     ew(ii)  = 0.0  ! >> EWs initialized to zero 
     n_points(ii)       = 0.0  ! *** Variables for error computation
     mean_wave_step(ii) = 0.0
     mean_flux(ii)      = 0.0
     sigma_flux(ii)     = 0.0
     mean_ptp(ii)       = 0.0

     do jj = 1, n_wave
        ! >> Fluxes at the limits of the integration (wave_low, wave_upp)
        if (abs(feature_wave_low(ii) - wave(jj)).lt.tol) then
             feature_flux_low(ii) = flx(jj)

        end if
        if (abs(feature_wave_upp(ii) - wave(jj)).lt.tol) then
             feature_flux_upp(ii) = flx(jj)
        end if
     end do


     ! >> Slope and ordinate of the line connecting the two integration limits
     slope(ii) = (feature_flux_low(ii) - feature_flux_upp(ii))/(feature_wave_low(ii) - feature_wave_upp(ii))
     ordinate(ii) = feature_flux_low(ii) - (slope(ii)*feature_wave_low(ii))

     ! >> Read again the spectra once more to compute the pseudo equivalent width 
     ! Flux peak to peak = "  slope(i)*wavelength(j) + ordinate(i) "
     do jj = 2, n_wave
        if (wave(jj).ge.feature_wave_low(ii).AND.wave(jj).le.feature_wave_upp(ii)) then

        wave_step = wave(jj) -  wave(jj-1)
        ew(ii) = ew(ii) +  ((((slope(ii)*wave(jj) + ordinate(ii)) - flx(jj))/ &
                           (slope(ii)*wave(jj) + ordinate(ii)))*wave_step)

        mean_flux(ii)      = mean_flux(ii) + flx(jj)
        mean_wave_step(ii) = mean_wave_step(ii) + wave_step
        mean_ptp(ii)       = mean_ptp(ii) + (slope(ii)*wave(jj) + ordinate(ii))
        n_points(ii)       = n_points(ii) + 1

       end if
     end do

     mean_flux(ii)      = mean_flux(ii)/n_points(ii)
     mean_wave_step(ii) = mean_wave_step(ii)/(n_points(ii))
     mean_ptp(ii)       = mean_ptp(ii)/(n_points(ii))

     ! >> Read a third time to compute the flux variance
     do jj = 2, n_wave
        if (wave(jj).ge.feature_wave_low(ii).AND.wave(jj).le.feature_wave_upp(ii)) then
            sigma_flux(ii) = sigma_flux(ii) + (flx(jj) - mean_flux(ii))**2
       end if
     end do

     sigma_flux(ii) = sqrt(sigma_flux(ii))/(n_points(ii) - 1.0)

    ! >> Sigma in the measured EWs
     e_ew(ii) = (mean_wave_step(ii)/mean_ptp(ii))*sigma_flux(ii)

    ! >> Fom A to mA
     ew(ii)   = ew(ii)*1000.
     e_ew(ii) = e_ew(ii)*1000.

    ! >> Set null_ew value for too weak/too strong lines 
     if (ew(ii).gt.max_ew.OR.ew(ii).lt.min_ew) then
         ew(ii)   = null_ew
         e_ew(ii) = null_ew
     end if

    ! >> Write results in output file (pseudo EWs in mA)
    write(2,55) feature_id(ii), feature_center(ii), ew(ii), e_ew(ii)

  end do  ! >> End of 'For each feature ...'

  close(2)
  
30 format('# Date ', i2.2, '/', i2.2, '/', i4.4, '; time ', i2.2, ':', i2.2, ':', i2.2 )
35 format(A48,1x,A12)
45 format(4(A10,1x))
50 format(A33,1x,A12)
55 format(A5,1x,3(f12.2,1x))

  end subroutine pseudoEW

! -------------------------------------------------------------------
! Fitsio subroutine printerror
! -------------------------------------------------------------------

  subroutine printerror(status)
  implicit none

  ! Print out the FITSIO error messages to the user

  integer  :: status
  character(30) :: errtext
  character(80) :: errmessage

  ! Check if status is OK (no error); if so, simply return
  if (status .le. 0)return

 ! Get the text string which describes the error
  call ftgerr(status,errtext)
  print *,'FITSIO Error Status =',status,': ',errtext

  ! Read and print out all the error messages on the FITSIO stack
  call ftgmsg(errmessage)
  do while (errmessage .ne. ' ')
     print *,errmessage
     call ftgmsg(errmessage)
  end do

  end subroutine printerror
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

