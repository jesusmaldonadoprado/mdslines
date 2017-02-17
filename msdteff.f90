!<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! MSDlines subroutine for computing the effective temperature
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
! last change: vie mar 20 15:04:17 CET 2015
! jmaldonado at inaf-oapa
!><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! --- Definition of variables
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine msdteff
   implicit none
   
   integer,parameter                  :: n_features      =  4424 
   integer,parameter                  :: n_ratios        =   112 
   character(80),parameter            :: input_file      = 'stars_list.txt'
   character(60),parameter            :: coeff_file      = 'coefficients_Mstars_ratios_teff_ver02.dat' 
   character(80),parameter            :: ofile           = 'derived_temperatures.dat'
   real,parameter                     :: null_value      = -99.
   character(30)                      :: ews_file
   character(15)                      :: star_identifier
   character(35)                      :: stellar_spectra 
   character(8),dimension(n_ratios)   :: numerador, denominador
   character(2),dimension(n_ratios)   :: equation
   real,dimension(n_ratios)           :: variance, aH, bH, cH, sH, aMH, bMH, cMH, sMH
   real,dimension(n_ratios)           :: aPL, bPL, sPL, aE, bE, sE, aL, bL, sL, bestS
   real,dimension(n_ratios)           :: min_ratio_value, max_ratio_value   
   real,dimension(n_features)         :: wavelength, width
   character(5),dimension(n_features) :: feature_id 
   real,dimension(n_features)         :: ew, e_ew 
   real,dimension(n_ratios)           :: ratio_numerador, ratio_denominador
   real,dimension(n_ratios)           :: ratio_enumerador, ratio_edenominador
   real,dimension(n_ratios)           :: ratio, e_ratio     
   real,dimension(n_ratios)           :: teff_ratio, eteff_ratio
   integer                            :: length, n_good
   real                               :: mean_teff, emean_teff
   integer                            :: ii, jj, kk, ios, jos

! ------------------------------------------------------
!  Open and read necessary files 
! ------------------------------------------------------

   ! *** File containing the fit to the selected ratios
   open(1,file=coeff_file,status='old',action='read')
   read(1,*)
   do ii = 1, n_ratios
   read(1,*) numerador(ii), denominador(ii), aH(ii), bH(ii), cH(ii), sH(ii),       &
    aMH(ii), bMH(ii), cMH(ii), sMH(ii), aPL(ii), bPL(ii), sPL(ii), aE(ii), bE(ii), &
    sE(ii), aL(ii), bL(ii), sL(ii), bestS(ii), equation(ii), min_ratio_value(ii), max_ratio_value(ii)
   end do
   close(1)

   ! *** Headers in the output file
   open(3,file=ofile,status='unknown')
   write(3,10)  '# Star_identifier, Teff(K),  eTeff(K)'

! ------------------------------------------------------
!  Open the files containing the EWs and compute Teff
! ------------------------------------------------------

   open(5,file=input_file,status='old',action='read')
   read(5,*) ! First three lines are useless 
   read(5,*)
   read(5,*)
   do
   read(5,*,iostat=ios) star_identifier, stellar_spectra
   if (ios/=0) exit

      length    = len_trim(star_identifier)
      ews_file  = star_identifier(1:length)//'_psdEWs.dat'

      open(2,file=ews_file,status='old',action='read')
      read(2,*)  ! Threee first lines are useless
      read(2,*)
      read(2,*)
           do ii = 1, n_features 
              read(2,*)  feature_id(ii), wavelength(ii), ew(ii), e_ew(ii)
           end do
      close(2) 

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  Compute effective temperature for each ratio
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      do ii = 1, n_ratios   ! *** For each ratio 
    
         do jj = 1, n_features  ! *** Searching for the numerator
              if (feature_id(jj).eq.numerador(ii)) then
                  ratio_numerador(ii)  = ew(jj)
                  ratio_enumerador(ii) = e_ew(jj)
              end if
         end do
    
           
         do jj = 1, n_features  ! *** Searching for the denominator
              if (feature_id(jj).eq.denominador(ii)) then 
                  ratio_denominador(ii) = ew(jj)
                  ratio_edenominador(ii) = e_ew(jj)
              end if
         end do   


         ! *** Computing the ratio and its uncertainty
         ! *** Checking if pseudoEWs have good values

         if (ratio_numerador(ii).eq.null_value.OR.ratio_denominador(ii).eq.null_value) then
             ratio(ii)       = null_value
             e_ratio(ii)     = null_value
             teff_ratio(ii)  = null_value  
             eteff_ratio(ii) = null_value    
         else


         ratio(ii)   = ratio_numerador(ii)/ratio_denominador(ii)  
         e_ratio(ii) = ratio(ii)*sqrt((ratio_enumerador(ii)/ratio_numerador(ii))**2 + &
                     (ratio_edenominador(ii)/ratio_denominador(ii))**2)

        ! *** Computing the effective temperature for each ratio
        ! *** Checking if ratio is within the calibrations limits

          
        if (ratio(ii).gt.max_ratio_value(ii).OR.ratio(ii).lt.min_ratio_value(ii)) then 
               ratio(ii)      = null_value
               e_ratio(ii)    = null_value
               teff_ratio(ii) = null_value
              eteff_ratio(ii) = null_value
         else
             
           ! *** Checking the functional form of the calibration
           if (equation(ii).eq.'HO') then  ! *** Hoerl function

               teff_ratio(ii)  = aH(ii)*(bH(ii)**(ratio(ii)))*(ratio(ii)**cH(ii))
                
               eteff_ratio(ii) = ((teff_ratio(ii)*log(bH(ii)))  + &
                                 (teff_ratio(ii)*cH(ii))/ratio(ii))*e_ratio(ii) 

           else if (equation(ii).eq.'MH') then ! *** Modified Hoerl function
                       
               teff_ratio(ii)  = aMH(ii)*(bMH(ii)**(1.0/ratio(ii)))*(ratio(ii)**cMH(ii))
              
               eteff_ratio(ii) = ((teff_ratio(ii)*log(bH(ii))/(ratio(ii)**2))  + & 
                                 (teff_ratio(ii)*cH(ii))/ratio(ii))*e_ratio(ii)

           else if (equation(ii).eq.'PL') then ! *** Power-law function
               
               teff_ratio(ii)  = aPL(ii)*(ratio(ii)**bPL(ii))
               eteff_ratio(ii) = (aPL(ii)*bPL(ii)*(ratio(ii)**(bPL(ii)-1.0)))*e_ratio(ii)

           else if (equation(ii).eq.'EX') then  ! *** Teff = a * b^(r)

               teff_ratio(ii)  = aE(ii)*(bE(ii)**ratio(ii))
               eteff_ratio(ii) = teff_ratio(ii)*log(bE(ii))*e_ratio(ii) 

           else if (equation(ii).eq.'LO') then  ! *** Teff = a * bln(r)

               teff_ratio(ii)  =  aL(ii) + bL(ii)*log(ratio(ii))    
               eteff_ratio(ii) =  bL(ii)*(e_ratio(ii)/ratio(ii))


           end if 
       
               eteff_ratio(ii) = sqrt(eteff_ratio(ii)**2 + bestS(ii)**2) 

           end if ! *** End of 'check if ratio is within the calibration's limits'

           end if ! *** End of 'check for null_value case' 


        end do  ! *** End of "for each ratiio ...'


    ! *** "Final" teff and error
       mean_teff  = 0.0
       emean_teff = 0.0
       n_good = 0   


    do ii = 1, n_ratios
        ! *** Just in case some fit gives a ridiculus value for some star
        if (teff_ratio(ii).gt.1000.AND.teff_ratio(ii).lt.6000.AND.teff_ratio(ii).ne.null_value) then
            mean_teff   = mean_teff  + teff_ratio(ii)
            emean_teff  = emean_teff + eteff_ratio(ii)
            n_good = n_good + 1
        end if
    end do


  ! *** At least one good teff is needed to give a median value
   if (n_good.ge.1) then
   mean_teff  = mean_teff/(n_good*1.0)
   emean_teff = emean_teff/(n_good*1.0)
   else
   mean_teff  = null_value
   emean_teff = null_value
   end if

   ! *** Writing results in output file
     write(3,15) star_identifier, mean_teff,  emean_teff

   end do ! *** End of 'For each star in input file ...'
   close(3)
   close(5)

10 format(A40) 
15 format(A15,1x,2(f12.2,1x)) 

! ---- end of program ----------------------------------
   end subroutine

