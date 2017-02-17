!<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! MSDlines subroutine for computing the stellar metallicity
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
! last change: dom mar  1 11:46:09 CET 2015
! jmaldonado at inaf-oapa
!><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! --- Definition of variables
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine msdmetal
   implicit none

   integer,parameter                  :: n_features     =  4424
   integer,parameter                  :: n_ratios       =   696
   character(80),parameter            :: input_file     = 'stars_list.txt'
   character(60),parameter            :: coeff_file     = 'coefficients_Mstars_ratios_metallicity_ver02.dat'
   character(80),parameter            :: ofile          = 'derived_metallicities.dat'
   real,parameter                     :: null_value     =  -99.
   character(30)                      :: ews_file
   character(15)                      :: star_identifier
   character(35)                      :: stellar_spectra
   real,dimension(n_ratios)           :: numerador, denominador, feature
   real,dimension(n_ratios)           :: variance, sigma, p0, p1, p2, rap
   real,dimension(n_features)         :: wavelength, width
   character(5),dimension(n_features) :: feature_id
   integer,dimension(n_ratios)        :: n1, n2, n3 
   real,dimension(n_features)         :: ew, e_ew
   real,dimension(n_ratios)           :: ratio_numerador, ratio_denominador
   real,dimension(n_ratios)           :: ratio_enumerador, ratio_edenominador
   real,dimension(n_ratios)           :: ratio_feature, ratio_efeature
   real,dimension(n_ratios)           :: ratio, e_ratio     
   real,dimension(n_ratios)           :: metal_ratio, emetal_ratio
   real,dimension(n_ratios)           :: min_ratio, max_ratio, min_feature, max_feature
   integer                            :: length 
   integer                            :: n_good
   real                               :: mean_metal, emean_metal
   integer                            :: ii, jj, kk, ios, jos

! ------------------------------------------------------
!  Open and read necessary files 
! ------------------------------------------------------

   ! *** File containing the fit to the selected ratios
   open(1,file=coeff_file,status='old',action='read')
   read(1,*)
   do ii = 1, n_ratios
   read(1,*) n1(ii), n2(ii), n3(ii), numerador(ii), denominador(ii), feature(ii), & 
             p0(ii), p1(ii), p2(ii), sigma(ii), rap(ii), min_ratio(ii), max_ratio(ii), &
             min_feature(ii), max_feature(ii)
   end do
   close(1)

   ! *** Headers in the output file
   open(3,file=ofile,status='unknown')
   write(3,10)  '# Star_identifier, Metallicity, eMetallicity'

! ------------------------------------------------------
!  Open the files containing the EWs and compute Teff
! ------------------------------------------------------

   open(5,file=input_file,status='old',action='read')
   read(5,*)  ! Three first lines
   read(5,*)
   read(5,*)
   do
   read(5,*,iostat=ios) star_identifier, stellar_spectra
   if (ios/=0) exit

      length    = len_trim(star_identifier)
      ews_file  = star_identifier(1:length)//'_psdEWs.dat'

      open(2,file=ews_file,status='old',action='read')
      read(2,*)  ! Three first lines are useless
      read(2,*)
      read(2,*)
           do ii = 1, n_features
              read(2,*)  feature_id(ii), wavelength(ii), ew(ii), e_ew(ii)
           end do
      close(2)


        do ii = 1, n_ratios   ! *** For each ratio 
    
           do jj = 1, n_features  ! *** Searching for the numerator
              if (wavelength(jj).eq.numerador(ii)) then
                  ratio_numerador(ii)  = ew(jj)
                  ratio_enumerador(ii) = e_ew(jj)
              end if
           end do
    
           
           do jj = 1, n_features  ! *** Searching for the denominator
              if (wavelength(jj).eq.denominador(ii)) then 
                  ratio_denominador(ii)   = ew(jj)
                  ratio_edenominador(ii)  = e_ew(jj)
              end if
           end do   


           do jj = 1, n_features ! *** Searching for the feature
              if (wavelength(jj).eq.feature(ii)) then
                  ratio_feature(ii)  = ew(jj)
                  ratio_efeature(ii) = e_ew(jj)
              end if
           end do 

                         
           ! *** Computing the ratio and its corresponding metallicity
           if (ratio_numerador(ii).ne.null_value.AND.ratio_denominador(ii).ne.null_value &
              .AND.ratio_feature(ii).ne.null_value) then

           ratio(ii)   = ratio_numerador(ii)/ratio_denominador(ii)  
           e_ratio(ii) = ratio(ii)*sqrt((ratio_enumerador(ii)/ratio_numerador(ii))**2 + &
                                 (ratio_edenominador(ii)/ratio_denominador(ii))**2)


           
           ! *** Cheking the calibration's limits
             if (ratio(ii).gt.min_ratio(ii).AND.ratio(ii).lt.max_ratio(ii).AND.  &
                 ratio_feature(ii).gt.min_feature(ii).AND.ratio_feature(ii).lt.max_feature(ii)) then  

           metal_ratio(ii)  = p0(ii)*ratio_feature(ii) + p1(ii)*ratio(ii) + p2(ii) 
           emetal_ratio(ii) = sqrt((p0(ii)*ratio_efeature(ii))**2 + (p1(ii)*e_ratio(ii))**2 &
                                  + sigma(ii)**2)  

              else
                metal_ratio(ii)  = null_value
                emetal_ratio(ii) = null_value
             end if  ! End of 'checking the calibration's limits  

            else
                ratio(ii)        = null_value
                e_ratio(ii)      = null_value  
                metal_ratio(ii)  = null_value
                emetal_ratio(ii) = null_value
            end if 


        end do  ! *** End of "for each ratiio ...'


    ! *** "Final" Metallicity and error
     mean_metal  = 0.0
     emean_metal = 0.0
     n_good  = 0  

 
     do ii = 1, n_ratios
       ! *** Just in case some fit gives a ridiculus value for some star
       if (metal_ratio(ii).ge.-3.0.AND.metal_ratio(ii).le.3.0.AND.metal_ratio(ii).ne. &
           null_value) then 
           mean_metal  = mean_metal  + metal_ratio(ii)
           emean_metal = emean_metal + emetal_ratio(ii)  
           n_good = n_good + 1
       end if
    end do  


  ! *** At least one good metallicity value is needed to give a median value
   if (n_good.ge.1) then
   mean_metal  = mean_metal/(n_good*1.0)
   emean_metal = emean_metal/(n_good*1.0)
   else
   mean_metal  = null_value
   emean_metal = null_value
   end if

  ! *** Writing results in output file
  write(3,15) star_identifier, mean_metal,  emean_metal

   end do ! *** End of 'For each star in input file ...'
   close(3)
   close(5)



10 format(A50) 
15 format(A35,1x,2(f12.2,1x)) 
20 format(A35,1x,13(f12.2,1x))

! ---- end of program ----------------------------------
   end subroutine

