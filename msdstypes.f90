!<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! MSDlines subroutine for determining the spectral type
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
! last change: dom mar  1 11:44:09 CET 2015
! jmaldonado at inaf-oapa
!><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! --- Definition of variables
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine msdstypes
   implicit none

   integer,parameter                  :: n_features     =  4424
   integer,parameter                  :: n_ratios       =    82
   character(80),parameter            :: input_file     = 'stars_list.txt'
   character(60),parameter            :: coeff_file     = 'coefficients_Mstars_ratios_sptype_ver02.dat'
   character(80),parameter            :: ofile          = 'derived_sptypes.dat'
   real,parameter                     :: null_value     =  -99.0 
   character(30)                      :: ews_file
   character(15)                      :: star_identifier
   character(35)                      :: stellar_spectra
   character(5),dimension(n_ratios)   :: numerador, denominador
   character(2),dimension(n_ratios)   :: equation
   real,dimension(n_ratios)           :: p0, p1, p2, p3, sigma, rap
   real,dimension(n_features)         :: wavelength, width
   character(8),dimension(n_features) :: feature_id 
   real,dimension(n_features)         :: ew, e_ew
   real,dimension(n_ratios)           :: ratio_numerador, ratio_denominador
   real,dimension(n_ratios)           :: ratio      
   real,dimension(n_ratios)           :: sptype_ratio
   real,dimension(n_ratios)           :: min_ratio, max_ratio
   real                               :: mean_sptype
   integer                            :: length
   integer                            :: n_good, entero
   real                               :: deci
   integer                            :: ii, jj, kk, ios, jos

! ------------------------------------------------------+
!  Open and read necessary files 
! ------------------------------------------------------+

   ! *** File containing the fit to the selected ratios
   open(1,file=coeff_file,status='old',action='read')
   read(1,*)
   do ii = 1, n_ratios
   read(1,*) numerador(ii), denominador(ii), p0(ii), p1(ii), p2(ii), p3(ii), sigma(ii), rap(ii), &
             min_ratio(ii), max_ratio(ii) 
   end do
   close(1)

   ! *** Headers in the output file
   open(3,file=ofile,status='unknown')
   write(3,10)  '# Star_identifier, sptype'


! ------------------------------------------------------
!  Open the files containing the EWs and compute SpType
! ------------------------------------------------------

   open(5,file=input_file,status='old',action='read')
   read(5,*) ! Three first lines
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
            if (feature_id(jj).eq.numerador(ii)) then
                ratio_numerador(ii)  = ew(jj)
            end if
         end do
    
           
         do jj = 1, n_features  ! *** Searching for the denominator
            if (feature_id(jj).eq.denominador(ii)) then 
                ratio_denominador(ii) = ew(jj)
            end if
         end do   


        ! *** Computing the ratio and its corresponding SpType
        if (ratio_numerador(ii).ne.null_value.AND.ratio_denominador(ii).ne.null_value) then   

        ratio(ii)   = ratio_numerador(ii)/ratio_denominador(ii) 

        if (ratio(ii).gt.min_ratio(ii).AND.ratio(ii).lt.max_ratio(ii)) then  

        sptype_ratio(ii) =  p0(ii) + p1(ii)*ratio(ii) + p2(ii)*(ratio(ii)**2)  &
                            +  p3(ii)*(ratio(ii)**3)    

        else

        ! *** If ratio is outside calibration's limits  
        sptype_ratio(ii) = null_value

        end if 

        else

        ! *** Null value case        
        ratio(ii)        = null_value
        sptype_ratio(ii) = null_value  
        
        end if

      end do  ! *** End of "for each ratiio ...'

    ! *** Mean spectral type
     mean_sptype  = 0.0
     n_good     = 0
   
     do ii = 1, n_ratios
        ! **** Just in case
        if (sptype_ratio(ii).gt.-3.0.AND.sptype_ratio(ii).lt.7.0.AND.sptype_ratio(ii) &
            .ne.null_value) then
        mean_sptype   = mean_sptype  + sptype_ratio(ii) 
        n_good        = n_good +1  
        end if
     end do 


  ! *** At least one good sptype is needed to give a median value
   if (n_good.ge.1) then
   mean_sptype  = mean_sptype/(n_good*1.0)
   else
   mean_sptype  = null_value
   end if


     if (mean_sptype.ne.null_value) then      
     ! *** Redondeamos
     entero = int(mean_sptype)
     deci   = mean_sptype - entero

     if (mean_sptype.lt.-1) then
         write(*,*) 'Oh no! You have a late-K star, check: ', star_identifier
     else if (mean_sptype.lt.0.AND.deci.le.-0.75) then
         mean_sptype = -1
     else if (mean_sptype.lt.0.AND.deci.gt.-0.75.AND.deci.le.-0.25) then
         mean_sptype = -0.5
     else if (mean_sptype.lt.0.AND.deci.gt.-0.25) then
         mean_sptype =  0
     else if (deci.le.0.25) then
         mean_sptype = mean_sptype - deci
     else if (deci.gt.0.25.AND.deci.le.0.75) then
         mean_sptype = entero + 0.5
     else if (deci.gt.0.75) then
         mean_sptype = entero + 1
     end if

     end if  ! End of null value case 


   ! *** Writing results in output file

     write(3,15) star_identifier, mean_sptype

   end do  ! End of "For each star in the input file ...'
   close(5)
   close(3)


10 format(A50) 
15 format(A15,1x,f12.1) 

   end subroutine

