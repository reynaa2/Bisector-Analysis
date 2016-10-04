;@/Volumes/CoveyData/APOGEE_Spectra/idl_routines/critical_apogeereduce/apload.pro
;---------------------------
FUNCTION measure_r, comb, velocity_axis, peak_loc

;peak location -> Needs to be defined according to the index that the maximum height occurs. But! You do not want the height of the CCF itself.
;measure r statistic for this cross correlation routine
;comb = the cross correlation functions 
;find the length of the cc function
n_cc = N_ELEMENTS(comb)
ccf_max = MAX(comb, peak_loc) 

;split the cross-correlation array on its
;peak
IF n_cc GT peak_loc*2 THEN BEGIN
     low_comb = comb[0:peak_loc]
     high_comb = comb[peak_loc:2*peak_loc]
ENDIF ELSE BEGIN
     low_comb = comb[2*peak_loc-n_cc+1:peak_loc]
     high_comb = comb[peak_loc:n_cc-1]
ENDELSE

;flip the low side
low_comb = REVERSE(low_comb)

IF N_ELEMENTS(low_comb) EQ N_ELEMENTS(high_comb) THEN BEGIN

    one_over_2n = 1./(2*N_ELEMENTS(low_comb))
    sigma_a = SQRT( one_over_2n * TOTAL( (high_comb-low_comb)^2.  ) )
    height = comb[peak_loc]
    r = height/(SQRT(2)*sigma_a)    


RETURN, r

ENDIF ELSE PRINT, 'problem with my r calculation!!!'

END
;---------------------------
FUNCTION SIMPLE_LORENTZIAN, X, P

bottom = ( (x - P[0])^2. + P[1]^2.)
RETURN,  P[2] * ( P[1]^2. / (bottom) )

END
;---------------------------
FUNCTION DUAL_LORENTZIANS, X, P

  RETURN,  simple_lorentzian(X, P[0:2]) + simple_lorentzian(X, P[3:5])

END
;---------------------------
FUNCTION LINKED_LORENTZIANS, X, P

  ;link the heights of the lorentzians by taking the relative height of the
  ;second peak from the input array, and converting to an absolute height for
  ;use by the simple lorentzian function
  secondary_array = [P[3],P[4],P[5]*P[2]]
  
  RETURN,  simple_lorentzian(X, P[0:2]) + simple_lorentzian(X, secondary_array)

END
;---------------------------
FUNCTION DOWNHILL_MASK, comb, lag, peak_lag

;LOADCT, 13, /SILENT
;OPLOT, [lag[primary_lag], lag[primary_lag]],[comb[primary_lag],comb[primary_lag]], PSYM = 5, SYMSIZE = 2, COLOR = color1

;find the subset of the cc function centered on the peak
still_looking_lower = 1
peak_low = peak_lag-5
WHILE still_looking_lower GT 0 DO BEGIN

   IF comb[peak_low-1] - comb[peak_low] LT 0 AND comb[peak_low-1] GT 0*comb[peak_lag] THEN peak_low = peak_low-1 ELSE still_looking_lower = 0
     
ENDWHILE

still_looking_upper = 1
peak_high = peak_lag+5
WHILE still_looking_upper GT 0 DO BEGIN

   IF comb[peak_high+1] - comb[peak_high] LT 0 AND comb[peak_high+1] GT 0*comb[peak_lag] THEN peak_high = peak_high+1 ELSE still_looking_upper = 0
     
ENDWHILE

;OPLOT, lag[primary_low:primary_lag], comb[primary_low:primary_lag], COLOR = color1 - 20
;OPLOT, lag[primary_lag:primary_high], comb[primary_lag:primary_high], COLOR = color1 + 20

n_lag = n_elements(lag)
mask = INTARR(n_lag)
mask[*] = 0
mask[peak_low:peak_high] = 1

RETURN, mask
END
;---------------------------
FUNCTION FIND_APOGEE_CCF_PEAKS, comb, lag

;find the primary peak
primary_peak_height = MAX(comb, primary_lag)
primary_peak = DOWNHILL_MASK(comb, lag, primary_lag)

search_for_secondary = WHERE(primary_peak EQ 0, n_secondary)
secondary_peak_height = MAX(comb[search_for_secondary], secondary_lag)

;make an array to return the results
RETURN, [lag[primary_lag], primary_peak_height, lag[search_for_secondary[secondary_lag]], secondary_peak_height]

END
;---------------------------
PRO FIT_SEPARATED_APOGEE_CCF_PEAKS, comb, lag, resolution, positions, widths, heights, primary_peak_fit, secondary_peak_fit, final_residuals, LORENTZIAN = lorentzian

weights = FLTARR(N_ELEMENTS(comb))

;find the primary peak and set initial guess accordingly
primary_peak_height = MAX(comb, primary_lag)
initial_primary_guess = [lag[primary_lag], resolution, primary_peak_height]

;make the mask for the primary peak
primary_mask = DOWNHILL_MASK(comb, lag, primary_lag)
in_peak = WHERE(primary_mask EQ 1, n_in_peak)

weights[*] = 1000.
weights[in_peak] = 1.

IF KEYWORD_SET(lorentzian) THEN BEGIN
  ; primary_fit = MPFITFUN('SIMPLE_LORENTZIAN', lag, comb, weights, initial_primary_guess, /SILENT)
   primary_peak_fit = SIMPLE_LORENTZIAN( lag , primary_fit)
ENDIF ELSE BEGIN
 ;  primary_fit = MPFITFUN('GAUSS1', lag, comb, weights, initial_primary_guess, /SILENT)
   primary_peak_fit = GAUSS1( lag , primary_fit)
ENDELSE

primary_residuals = comb-primary_peak_fit

;find the secondary peak in the residuals and set initial guess accordingly
secondary_peak_height = MAX(primary_residuals, secondary_lag)
initial_secondary_guess = [lag[secondary_lag], resolution, secondary_peak_height]

;make the mask for the secondary peak in the residuals
secondary_mask = DOWNHILL_MASK(primary_residuals, lag, secondary_lag)
in_sec_peak = WHERE(secondary_mask EQ 1, n_in_sec_peak)

weights[*] = 1000.
weights[in_sec_peak] = 1.
; Filler from last summer
IF KEYWORD_SET(lorentzian) THEN BEGIN
  ; secondary_fit = MPFITFUN('SIMPLE_LORENTZIAN', lag, primary_residuals, weights, initial_secondary_guess, /SILENT)
   secondary_peak_fit = SIMPLE_LORENTZIAN( lag , secondary_fit)
ENDIF ELSE BEGIN
 ;  secondary_fit = MPFITFUN('GAUSS1', lag, primary_residuals, weights, initial_secondary_guess, /SILENT)
   secondary_peak_fit = GAUSS1( lag , secondary_fit)
ENDELSE

final_residuals = primary_residuals-secondary_peak_fit
combined_peak_fit = primary_peak_fit + secondary_peak_fit

positions = FLTARR(2,2)
positions[1,*] = [0,0]
widths = FLTARR(2,2)
widths[1,*] = [0,0]
heights = FLTARR(2,2)
heights[1,*] = [0,0]

IF primary_fit[2] GT secondary_fit[2] THEN BEGIN
   positions[0,*] = [primary_fit[0],secondary_fit[0]]
   widths[0,*] = ABS([primary_fit[1],secondary_fit[1]])
   heights[0,*] = [primary_fit[2],secondary_fit[2]]
ENDIF ELSE BEGIN
   positions[0,*] = [secondary_fit[0],primary_fit[0]]
   widths[0,*] = ABS([secondary_fit[1],primary_fit[1]])
   heights[0,*] = [secondary_fit[2],primary_fit[2]]
ENDELSE

;now plot the results
LOADCT, 0, /SILENT
PLOT, lag, comb, XTITLE = 'Lag (pixel space)', CHARSIZE = 1.5
LOADCT, 13, /SILENT
OPLOT, [lag[primary_lag], lag[primary_lag]], [primary_peak_height, primary_peak_height], PSYM = 6, SYMSIZE = 2, COLOR = primary_color
OPLOT, [lag[secondary_lag], lag[secondary_lag]], [secondary_peak_height, secondary_peak_height], PSYM = 5, SYMSIZE = 2, COLOR = secondary_color

primary_color = 100
secondary_color = 200

OPLOT, lag, primary_peak_fit, COLOR = primary_color, LINESTYLE = 4
OPLOT, lag, secondary_peak_fit, COLOR = secondary_color, LINESTYLE = 4
OPLOT, lag, combined_peak_fit, COLOR = 255, LINESTYLE = 2
OPLOT, lag, final_residuals, COLOR = (primary_color+secondary_color)/2., LINESTYLE = 3

;wait = GET_KBRD(1)

END
;-------------------------------
PRO FIT_LINKED_APOGEE_CCF_PEAKS, comb, lag, resolution, max_sep, positions, widths, heights, primary_peak_fit, secondary_peak_fit, final_residuals, CONSTRAINTS = constraints

initial_guess = FLTARR(6)
peak_parameters = FIND_APOGEE_CCF_PEAKS(comb, lag)

;set weights
weights = FLTARR(N_ELEMENTS(comb))
weights[*] = 1.

we_care = WHERE(comb GT 0.4*peak_parameters[1], COMPLEMENT = dont_care)
weights[dont_care] = 1000.

;test if there is any info from previous fits
IF KEYWORD_SET(constraints) THEN BEGIN
   
   pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},6)
   secondary_relative_height = peak_parameters[3]/peak_parameters[1]

      initial_guess[0] = peak_parameters[0]
      initial_guess[3] = peak_parameters[2]

      PRINT, 'Key numbers: ',ABS(initial_guess[0] - initial_guess[3]), max_sep
      
      IF ABS(initial_guess[0] - initial_guess[3]) GT max_sep*1.1 THEN initial_guess[3] = initial_guess[0]

      IF positions[1,0] EQ 1 THEN BEGIN 
         pi[0].limited[0] = 1
         pi[0].limits[0] = initial_guess[0]-resolution*0.5
         pi[0].limited[1] = 1
         pi[0].limits[1] = initial_guess[0]+resolution*0.5
      ENDIF

      IF positions[1,1] EQ 1 THEN BEGIN 
         pi[3].limited[0] = 1
         pi[3].limits[0] = initial_guess[0]-max_sep*1.1
         pi[3].limited[1] = 1
         pi[3].limits[1] = initial_guess[0]+max_sep*1.1
         PRINT, 'secondary position limits', pi[3].limits[0], pi[3].limits[1] 
      ENDIF

      IF positions[1,0] EQ 2 THEN BEGIN 
         pi[0].fixed[0] = 1
      ENDIF

      IF positions[1,1] EQ 2 THEN BEGIN 
         pi[3].fixed[0] = 1
      ENDIF

      initial_guess[1] = widths[0,0]
      initial_guess[4] = widths[0,1]

      IF widths[1,0] EQ 1 THEN BEGIN 
         pi[1].limited[0] = 1
         pi[1].limits[0] = widths[0,0]*0.75
         pi[1].limited[1] = 1
         pi[1].limits[1] = widths[0,0]*1.25
      ENDIF

      IF widths[1,1] EQ 1 THEN BEGIN 
         pi[4].limited[0] = 1
         pi[4].limits[0] = widths[0,1]*0.75
         pi[4].limited[1] = 1
         pi[4].limits[1] = widths[0,1]*1.25
      ENDIF

      IF widths[1,0] EQ 2 THEN BEGIN 
         pi[1].fixed[0] = 1
      ENDIF

      IF widths[1,1] EQ 2 THEN BEGIN 
         pi[4].fixed[0] = 1
      ENDIF

      initial_guess[2] = heights[0,0]
      initial_guess[5] = heights[0,1]

      IF heights[1,0] EQ 1 THEN BEGIN 
         pi[2].limited[0] = 1
         pi[2].limits[0] = initial_guess[2]*0.75
         pi[2].limited[1] = 1
         pi[2].limits[1] = initial_guess[2]*1.25
      ENDIF

      IF heights[1,1] EQ 1 THEN BEGIN 
         pi[5].limited[0] = 1
         pi[5].limits[0] = initial_guess[5]*0.9
         pi[5].limited[1] = 1
         pi[5].limits[1] = initial_guess[5]*1.1
      ENDIF

      IF heights[1,0] EQ 2 THEN BEGIN 
         pi[2].fixed[0] = 1
      ENDIF

      IF heights[1,1] EQ 2 THEN BEGIN 
         pi[5].fixed[0] = 1
      ENDIF

  ; linked_lorentzian_fit = MPFITFUN('LINKED_LORENTZIANS', lag, comb, weights, initial_guess, PARINFO = pi, /SILENT)

ENDIF ELSE BEGIN

   initial_guess[0] = peak_parameters[0]
   initial_guess[3] = peak_parameters[2]
   initial_guess[1] = resolution
   initial_guess[4] = resolution
   initial_guess[2] = peak_parameters[1]
   initial_guess[5] = peak_parameters[3]/peak_parameters[1]

;   linked_lorentzian_fit = MPFITFUN('LINKED_LORENTZIANS', lag, comb, weights, initial_guess, /SILENT)

ENDELSE

;save the linked_lorentzian_fit info into the position-width-height
;arrays, after checking to make sure the primary is bigger (but
;because AREA is conserved, checking on height might not make sense...

IF linked_lorentzian_fit[5] LT 1 THEN BEGIN
   positions[0,*] = [linked_lorentzian_fit[0], linked_lorentzian_fit[3]]
   widths[0,*] = ABS([linked_lorentzian_fit[1], linked_lorentzian_fit[4]])
   heights[0,*] = [linked_lorentzian_fit[2], linked_lorentzian_fit[5]]
ENDIF ELSE BEGIN
   positions[0,*] = [linked_lorentzian_fit[3], linked_lorentzian_fit[0]]
   widths[0,*] = ABS([linked_lorentzian_fit[4], linked_lorentzian_fit[1]])
   heights[0,*] = [linked_lorentzian_fit[5]*linked_lorentzian_fit[2], linked_lorentzian_fit[2]/(linked_lorentzian_fit[5]*linked_lorentzian_fit[2])];
ENDELSE


primary_peak_fit = SIMPLE_LORENTZIAN(lag, linked_lorentzian_fit[0:2])
secondary_peak_fit = SIMPLE_LORENTZIAN(lag, [linked_lorentzian_fit[3],linked_lorentzian_fit[4],linked_lorentzian_fit[5]*linked_lorentzian_fit[2]] )
combined_peak_fit = primary_peak_fit+secondary_peak_fit
final_residuals = comb-combined_peak_fit

;now plot the results
LOADCT, 0, /SILENT
PLOT, lag, comb, XTITLE = 'Lag (pixel space)', CHARSIZE = 1.5, XRANGE = [-40,40], /XSTYLE
LOADCT, 13, /SILENT
;OPLOT, [lag[primary_lag], lag[primary_lag]], [primary_peak_height, primary_peak_height], PSYM = 6, SYMSIZE = 2, COLOR = primary_color
;OPLOT, [lag[secondary_lag], lag[secondary_lag]], [secondary_peak_height, secondary_peak_height], PSYM = 5, SYMSIZE = 2, COLOR = secondary_color

primary_color = 100
secondary_color = 200

OPLOT, lag, primary_peak_fit, COLOR = primary_color, LINESTYLE = 4
OPLOT, lag, secondary_peak_fit, COLOR = secondary_color, LINESTYLE = 4
OPLOT, lag, combined_peak_fit, COLOR = 255, LINESTYLE = 2
OPLOT, lag, final_residuals, COLOR = (primary_color+secondary_color)/2., LINESTYLE = 3
END
;-------------------------------
FUNCTION CCF_WIDTH, ccf, ccflag
max_lag = MAX(ccflag)
min_lag = MIN(ccflag)
ccf_max = MAX(ccf)
ccf_min = MIN(ccf)
height = ccf_max - ccf_min

above_half = WHERE(ccf GE 0.5*ccf_max, n_above_half)
IF n_above_half GT 1 THEN BEGIN
  low_above_half = MIN(ccflag[above_half])
  high_above_half = MAX(ccflag[above_half])
  fwhm = high_above_half - low_above_half
ENDIF ELSE fwhm = 0.0

RETURN, FWHM

END
;---------------------------
FUNCTION CCF_HEIGHT , ccf, ccflag
  max_lag = MAX(ccflag)
  min_lag = MIN(ccflag)
  ccf_max = MAX(ccf)
  ccf_min = MIN(ccf)
  height = ccf_max - ccf_min
  
 RETURN, height
 END
;---------------------
FUNCTION CCF_THRESHOLD, ccf, ccflag
  max_lag = MAX(ccflag)
  min_lag = MIN(ccflag)
  ccf_max = MAX(ccf)
  ccf_min = MIN(ccf)
  height = ccf_max - ccf_min
  ;Normalize the CCF in the y-direction so peak=1 and floor=0
  norm_ccf = (ccf - ccf_min)/(ccf_max-ccf_min)

  ;make finer grid via interpolation. Create grid size
  finer_grid = FINDGEN((max_lag-min_lag)*20.)/20. + min_lag
  finer_ccf = INTERPOL(norm_ccf, ccflag,finer_grid)
  
  ;Find the ccf threshold
  ccf_median = MEDIAN(finer_ccf)
  threshold = 2.0 * ccf_median

RETURN, threshold
END
;----------------------
FUNCTION MEASURE_BISECTOR, ccflag, ccf

;find the height of the CCF
ccf_max = MAX(ccf)
ccf_min = MIN(ccf)
height = ccf_max-ccf_min
;print, 'Real CCF Max and Min: ', ccf_max, ' ' ,ccf_min
;find the boundaries that split the CCF into n_zones
n_boundaries = 6

zone_boundaries = ccf_min+(height*FINDGEN(n_boundaries)/(n_boundaries - 1)) ;This generates the horizontal slices to find the bisector pt for each region
;print, 'Zone Boundaries: ', zone_boundaries ; This is calulating the proper incriments between 1 and 0.8

;make arrays to store the vertical middle of each zone, and the
;resultant horizontal bisector [note that we can only measure
;n_zones-1 bisectors, as there are more boundaries]
bisector_positions = FLTARR(2,n_boundaries-1)

;loop through the zones and measure the bisector of each 
FOR i=0, n_boundaries - 2 DO BEGIN
  
  ;save y position for this part of the bisector
  bisector_positions[1,i] = (zone_boundaries[i] + zone_boundaries[i+1]) / 2.   ;(This seems like this is the range for the x-value bisector not the y)
  
  ;identify portion of CCF in between boundaries
  in_boundaries = WHERE(ccf GE zone_boundaries[i] AND ccf LE zone_boundaries[i+1], n_in_boundaries)
  ;PRINT, 'In-boundaries value: ', in_boundaries

 
                                ;find the mean x position between
                                ;these boundaries, if it exists; if
                                ;not, save -9999. to indicate that
                                ;measurement is not real. 
  IF n_in_boundaries GT 1 THEN BEGIN
      bisector_positions[0,i] = MEAN(ccflag[in_boundaries])
   ENDIF ELSE BEGIN
      bisector_positions[0,i] = -9999.
      ;Giving false data to signify the error in finding a position in a boundary
   ENDELSE  
ENDFOR
RETURN, bisector_positions
END
  
;---------------------------- 
FUNCTION measure_APOGEE_CCF_bisector, CCF, ccflag

  LOADCT, 0, /SILENT
  ;PLOT, ccflag, ccf, PSYM =10
;Mid point of ccf lies on (0,0)

;print, ccf where ccf is the data that lies within the bottom cut off pt in bright red and the threshold in dark red;
;print, 'The difference in ccf: ', MAX(ccf) -ABS(MIN(ccf)), 'Max: ', MAX(ccf), '          ', 'Min: ', MIN(ccf)

  ;find the horizontal range of the CCF
  max_lag = MAX(ccflag)
  min_lag = MIN(ccflag)
  
  ;find the vertical range of the CCF
  ccf_max = MAX(ccf)
 ; print, 'ccf max for vertical range: ', max_lag ; This looks like the horizontal range of the CCF not the vertical!
  ccf_min = MIN(ccf)
  
;Actual height before normalizing ccf
  height = ccf_max - ccf_min
  
  ;normalize the CCF vertically so that peak = 1 and floor = 0
  norm_ccf = (ccf - ccf_min)/(ccf_max - ccf_min)

  ;make a grid with 10x sampling (so
  ;that we are less susceptible to
  ;sampling effects.
  finer_grid = FINDGEN((max_lag-min_lag)*20)/20.+min_lag

  ;interpolate the CCF onto this finer grid
  finer_ccf = INTERPOL(norm_ccf, ccflag, finer_grid)
 
  ;find the CCF's base level
  ccf_median = MEDIAN(finer_ccf)
  truncate = 0.0

  ;Test the bisector calculation without the threshold:
  to_analyze = WHERE(finer_ccf GT truncate, n_analyze)
  IF n_analyze GT 3 THEN BEGIN
    
     ; bisector = MEASURE_BISECTOR(finer_grid[to_analyze], finer_ccf[to_analyze])
  ;  ;Array that extracts the endpoints to match the total number of visits provided
  
      bisector = MEASURE_BISECTOR(finer_grid[to_analyze-2], finer_ccf[to_analyze-2])
      
  ;print, 'this is the bisector:' , bisector
  ;  x_positions = bisector[0,*]
  ;  y_positions = bisector[1,*]
  ;  print, 'X- Positions: ' ,x_positions , ' ' , 'Y-Positions: ', y_positions
  ;  x_median = MEDIAN(x_positions)
  ;  y_median = MEDIAN(y_positions)
  ;  Print, 'X Median: ', x_median ;, '     ', 'Y Median: ', y_median
    ;Extracting the 2 endpoints of the array to have values that match the # of visits
    ;bisector_edited_mean = TOTAL(bisector)/N_ELEMENTS(bisector)
    ;print, 'Edited mean bisector x-values:' , bisector_edited_mean
  ;wait = GET_KBRD(1)

    ;zoom in to the part of the CCF we are analyzing; find the plot
    ;range we need to see that clearly.
    min_lag_to_analyze = MIN(finer_grid[to_analyze])
    max_lag_to_analyze = MAX(finer_grid[to_analyze])
    zoom_range = max_lag_to_analyze - min_lag_to_analyze
  
  ;Double check bisector y values to makes sure they lie in the appropriate range for each zone (why isn't the math working accordingly?)
 ; OPLOT, bisector[0,*], bisector[1,*], COLOR = 160 ; <- Plotting the Bisector
 ;wait = GET_KBRD(1)
  ; ENDIF ELSE BEGIN
;  bisector = FLTARR(2,9)
;  bisector[*,*] = -9999.
;
;ENDELSE
                                ;set the threshold above the base
                                ;level that we will measure the
  ;                              ;bisector for
  ;threshold = 2.0
  ;bthreshold =threshold*ccf_median
  ;print, 'Threshold: ', bthreshold
  ;show finer CCF
 ;LOADCT, 0, /SILENT
; PLOT, finer_grid, finer_ccf, PSYM=10   ;, COLOR = 255
;LOADCT, 13, /SILENT
 ;OPLOT, [-250,250], [ccf_median, ccf_median], COLOR = 255
; OPLOT, [-250,250], [threshold*ccf_median, threshold*ccf_median], COLOR = 80
  ;only analyze sources whose CCF
  ;exceeds 2.5x the median level (ie,
  ;ignore the broadest, smoothest sources)
  ;to_analyze = WHERE(finer_ccf GT threshold*ccf_median, n_to_analyze)  
 ; IF n_to_analyze GT 3 THEN BEGIN
    ; bisector = MEASURE_BISECTOR(finer_grid[to_analyze], finer_ccf[to_analyze])
     ;Array that extracts the endpoints to match the total number of visits provided
    ; bisectors = MEASURE_BISECTOR(finer_grid[to_analyze-2], finer_ccf[to_analyze-2])
     ;print, 'this is the bisector:' , bisector
    ; x_positions = bisector[0,*]
    ; y_positions = bisector[1,*]
   ;  print, 'X- Positions: ' ,x_positions , ' ' , 'Y-Positions: ', y_positions
    ; x_median = MEDIAN(x_positions)
    ; y_median = MEDIAN(y_positions)   
    ; Print, 'X Median: ', x_median ;, '     ', 'Y Median: ', y_median 
     ;Extracting the 2 endpoints of the array to have values that match the # of visits
    ; bisector_edited_mean = TOTAL(bisectors)/N_ELEMENTS(bisectors)
     ; print, 'Edited mean bisector x-values:' , bisector_edited_mean
     ; 
  ;wait = GET_KBRD(1)
  ;zoom in to the part of the CCF we are analyzing; find the plot 
  ;range we need to see that clearly.
 ; min_lag_to_analyze = MIN(finer_grid[to_analyze])
 ; max_lag_to_analyze = MAX(finer_grid[to_analyze])
 ; zoom_range = max_lag_to_analyze - min_lag_to_analyze

  LOADCT, 0, /SILENT
  ; -> ->PLOT, finer_grid, finer_ccf, XRANGE = [min_lag_to_analyze-0.1*zoom_range, max_lag_to_analyze+0.1*zoom_range], /XSTY , PSYM=10  ;, COLOR = 255
  ;Trying to plot on white backgrounds and enable saving procedures
  plot1 = PLOT(finer_grid, finer_ccf, XRANGE = [min_lag_to_analyze-0.1*zoom_range, max_lag_to_analyze+0.1*zoom_range],LINESTYLE=0,/buffer)
  
  LOADCT, 13, /SILENT
  ;-> -> OPLOT, [-250,250], [ccf_median, ccf_median], COLOR = 255
  ;-> -> OPLOT, [-250,250], [ccf_median, ccf_median], COLOR = 80
  ;plot2 = PLOT([-250,250], [ccf_median, ccf_median], /overplot,LINESTYLE=2,/buffer)
  ;plot3 = PLOT([-250,250], [ccf_median, ccf_median],/overplot, COLOR = 80,/buffer )
  plot5 = PLOT([-250,250],[0.1,0.1], /overplot, LINESTYLE=2,/buffer)
  plot6 = PLOT([-250,250],[0.3,0.3], /overplot, LINESTYLE=2, /buffer)
  plot7 = PLOT([-250,250], [0.5,0.5], /overplot, LINESTYLE=2,/buffer)
  plot8 = PLOT([-250,250], [0.7,0.7], /overplot, LINESTYLE=2,/buffer)
  plot9 = PLOT([-250,250], [0.9,0.9], /overplot, LINESTYLE=2,/buffer)

 ;Double check bisector y values to makes sure they lie in the appropriate range for each zone (why isn't the math working accordingly?)
  ;OPLOT, bisector[0,*], bisector[1,*], COLOR = 160 ; <- Plotting the Bisector
  plot_bisector = PLOT(bisector[0,*], bisector[1,*],/overplot, 'b', XTITLE = 'CCF Lag', YTITLE='$\wedge{CCF Units}$',/buffer)
  
 ; wait = GET_KBRD(1)
  ENDIF ELSE BEGIN

     bisector = FLTARR(2,5)
     bisector[*,*] = -9999.
;
  ENDELSE
  RETURN, bisector
END
;----------------------------
PRO APOGEE_spec_binary_classify, apStarFile, apogeeID, field
;read in the apstar file
apload, apStarFile, binary

;help,binary
;
;help, binary.head.head , /st
;print, binary.head

;make a name for plotting the spectrum and the table
short_name = STRMID(apStarFile,15,23)
;print, short_name

;make my life a **lot** easier by stripping off the combined CCFs
all_ccfs = N_ELEMENTS(binary.rv.ccf[0,*])
;print, 'all ccfs: ', all_ccfs

IF all_ccfs GE 3 THEN BEGIN
   only_visit_ccfs = binary.rv.ccf[*,2:all_ccfs-1]
   n_visits = all_ccfs-2
ENDIF ELSE BEGIN  ;If the source has only 1 visit it hsouldn't be read into the next for loop
   only_visit_ccfs = binary.rv.ccf[*,*]
   n_visits = all_ccfs
ENDELSE

PRINT, 'Number of Visits:',n_visits

;set guess for APOGEE_resolution
resolution = 3.
cspeed = 2.99792458d5
;!p.multi=[0,1,1]
SET_PLOT, 'x'
;!P.FONT=3
!P.THICK=3
;!P.CHARSIZE=1.25
;!P.CHARTHICK=2.
;!X.THICK=5
;!Y.THICK=5
;LOADCT, 0, /SILENT
;PLOTSYM, 0, /FILL

bisectors = FLTARR(n_visits,2,5)
;Zone = FLTARR(n_visits,2,5)
Rs = FLTARR(n_visits)
bx_range = FLTARR(n_visits)
CCFWidth = FLTARR(n_visits)

FOR i=0L, n_visits - 1 DO BEGIN
  
  ;Finding the Visit IDs with the plate-date-fiber
  paramTofind = 'SFILE'+STRTRIM(STRING(i+1),2)
  visitID = SXPAR(binary.head, paramTofind)
  visit_ID = STRMID(visitID,11,14)

  ;print, 'Visit ID: ', visit_ID
  ;wait = get_kbrd(1)
  ;
; Location of the stats to be calculated such as:
  ;Range of x values in the bisector (DONE)
  ;Mean, median and dfference of the x values compared to the two in the bisector (DONE)
  ;Median of Y values (DONE)
  ;All the x and y positions and the number of data points used for each bisector pt (DONE)
  ;(Number of points to make a bisector pt. -> 2
  ;(all x and y positions of bisector is in bisector[0,*], (x-values) and  bisector[1,*] (y-values) 
  ;Height of the peak (DONE)
  ; R statistic (DONE)
  ; Threshold stat (median for the y values) (DONE)
  ;Width of the peak

  ;Things to add later:
  ; Slope and correaltion coeffienct to linear bisector -> Heliocentric RV to make Helio vs. RV slope 
  ; Divide known binaries into well separated and merged peaks or single stars
;----------

   ;send this CCF to our bisector measuring code
   bisectors[i,*,*] = measure_APOGEE_CCF_bisector(only_visit_ccfs[*,i], binary.rv.ccflag)
   
   
   x_positions = REFORM(bisectors[i,0,*])
   y_positions = bisectors[i,1,*]
 
    x_median = MEDIAN(x_positions)
    y_median = MEDIAN(y_positions)
 
   ;Try moving the calculation of the bisector ranges to cover each bisector:
   max_x = MAX(bisectors[i,0,*])
   min_x = MIN(bisectors[i,0,*])

    ;Mean of x-values in the bisectors, this includes the beginning and end elements of the array.
    bisector_xmean = MEAN(bisectors[i,0,*])
  
    ;Difference between x mean and x median
    diff_in_median_mean = x_median - bisector_xmean

    ;Calling the R function
    R_value =  MEASURE_R(only_visit_ccfs[*,i], binary.rv.ccflag, ccf_height)

    ;Finding the x-range for the bisectors
    bisector_x_range = max_x - min_x
    
    ;Find minimum R values out of a visits for a star
    Rs[i] = MEASURE_R(only_visit_ccfs[*,i], binary.rv.ccflag, ccf_height)
    Rmin = [MIN(Rs),STRTRIM(i,1)]
   ; print, 'R min and visit: ', Rmin
    
    bx_range[i] = max_x - min_x
    maxXR = [MAX(bx_range), STRTRIM(i,1)] 
    ;print, 'max x-range and visit: ', maxXR

    XR_average = MEAN(bx_range)
    FWHM = CCF_WIDTH(only_visit_ccfs[*,i], binary.rv.ccflag)
    height_ccf = CCF_HEIGHT(only_visit_ccfs[*,i],binary.rv.ccflag)
    threshold = CCF_THRESHOLD(only_visit_ccfs[*,i],binary.rv.ccflag)

    CCFWidth[i] = CCF_WIDTH(only_visit_ccfs[*,i], binary.rv.ccflag)
    MaxWidth = MAX(CCFWIDTH)
    
  ;Write to file
   PRINTF,1,FORMAT = '(A20, A20,3X,F6.2, 3X,F6.2,3X, F10.2,3X,F6.2,3X,F6.2,3X,F10.2 ,3X,F10.2, 3X,F10.2, 3X, F10.2,3X,F10.2,3X,F8.2,3X,F8.2,3X,F10.2)', STRTRIM(apogeeID), STRTRIM(visit_ID),R_value, $;, STRTRIM(threshold), $
   bisector_x_range,x_median,bisector_xmean,diff_in_median_mean,x_positions[0],x_positions[1],x_positions[2] ,x_positions[3] ,x_positions[4], FWHM, height_ccf, threshold ;, STRTRIM(), STRTRIM(), STRTRIM(),STRTRIM(),STRTRIM()

   PRINTF,3, FORMAT = '(A10,3X,A20, 3X,F10.2 ,3X,F10.2, 3X,F10.2, 3X, F10.2,3X,F10.2,3X,F10.2)',STRTRIM(Field),STRTRIM(apogeeID), x_positions[0],x_positions[1],x_positions[2] ,x_positions[3] ,x_positions[4], bisector_x_range

ENDFOR

PRINTF,2, FORMAT = '(A10,3X,A20,3X,F10.2,2X,I5 ,3X , F10.2, 2X,I5, 3X ,F10.2,3X,F10.2)' ,STRTRIM(field), STRTRIM(apogeeID), STRTRIM(Rmin), STRTRIM(maxXR), MaxWidth, XR_average


LOADCT, 0, /SILENT
PLOT, [-9999.,-9999.],[0,1], XRANGE = [min_x, max_x], YRANGE = [0,1], /XSTY, /YSTY

FOR i=0, n_visits-1 DO BEGIN
   ;-> -> OPLOT, bisectors[i,0,*], bisectors[i,1,*]
   ; Trying to plot with white bkacgrounds anf saving features:
   plot4 = PLOT(bisectors[i,0,*], bisectors[i,1,*], 'r', /overplot, XTITLE='CCF Lag', YTITLE='$\wedge{CCF Units}$', TITLE='APOGEE ID:' + STRTRIM(string(field)) +'/'+ STRTRIM(string(apogeeID)) $
     +'   ' + 'Visit:'+STRTRIM(string(i)),/buffer)

ENDFOR
plot4.save, '/Volumes/CoveyData/APOGEE_Spectra/APOGEE2_DR13/Bisector/DR13_CCF_Bisector_Plots/' + apogeeID+'.png'
;wait = GET_KBRD(1)
END
;-----
PRO classify_all_APOGEE_spec_binaries, spec_list
  ;feed a list of stars and record the values into a single file (writing the .txt should be located here)

output = 'DR13_'
READCOL, spec_list, FORMAT = '(X,X,X,A,I,X,A)', twoMASSID,field, path ;(For all Lists)

;  output = 'DR12_SB9'
 ; READCOL, spec_list, FORMAT = '(X,X,X,X,A,I,X,A)', twoMASSID,field, path
 ;ENDIF  ELSE BEGIN
;
 ;  ;read formatted for Drew's List and APOGEE control sample
 ;  READCOL, spec_list,FORMAT = '(A,I,A)',twoMASSID,field, path
 ;  
; ENDELSE
;
;read formatted for all DR12 list
;READCOL, spec_list, FORMAT = '(X,X,X,I,X,A)', field, path

  ;0   0                VESTA     1                 apStar-r5-VESTA.fits                    1/apStar-r5-VESTA.fits

n_spectra = N_ELEMENTS(path)
;print, 'n_spectra: ', n_spectra

CLOSE, 1
OPENW, 1,'Bisector/' + output+ 'Bisector_Stats.lis'
; Missing piece -> 'Threshold Value',
PRINTF,1,'2MASS_ID', 'Visit_ID','R',' x_range ', ' x_median', ' x_mean ','  Diff_in_x', $
      '  x@y=0.1', '  x@y=0.3', '  x@y=0.5', '  x@y=0.7', '  x@y=0.9','FWHM','Height','  Thresh ' ,FORMAT = '(A20,A20,A10,A10,A10,A10,A11,A13,A13,A13,A13,A13,A13, A10,A10,A10,A10)';,  FORMAT = '(A30, A20, A20, A20, A20,A20,A20,A20,A20,A20,A20,A20,A20,A20,A20)'

CLOSE,2
 OPENW,2, 'Bisector/'+output+'Rmin_XR_Max.lis'
 PRINTF,2,'Field','2MASS_ID','R_min', 'XR_Max','Max_Width','XR_average', FORMAT = '(A10,A20,A20,A20,A20,A13)'
 CLOSE,3
 OPENW,3, 'Bisector/'+ output+ 'Bisector_pts.lis'
 PRINTF,3,'Location ID','2MASS_ID', '  x@y=0.1', '  x@y=0.3', '  x@y=0.5', '  x@y=0.7', '  x@y=0.9','x_range', FORMAT ='(A12,A20,A13,A13,A13,A13,A13,A13)'
    
FOR i=0L, n_spectra-1 DO BEGIN
   IF field[i] NE 1 THEN APOGEE_SPEC_BINARY_CLASSIFY, path[i], twoMASSID[i], field[i]
   
ENDFOR

CLOSE,1
CLOSE,2
CLOSE,3
END

