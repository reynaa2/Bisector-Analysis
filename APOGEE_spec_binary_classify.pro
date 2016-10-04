
FUNCTION measure_r, comb, velocity_axis, peak_loc

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

;This generates the horizontal slices to find the bisector pt for each region
zone_boundaries = ccf_min+(height*FINDGEN(n_boundaries)/(n_boundaries - 1)) 
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

  ;find the horizontal range of the CCF
  max_lag = MAX(ccflag)
  min_lag = MIN(ccflag)
  
  ;find the vertical range of the CCF
  ccf_max = MAX(ccf)
 ; print, 'ccf max for vertical range: ', max_lag
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
  ;wait = GET_KBRD(1)

    ;zoom in to the part of the CCF we are analyzing; find the plot
    ;range we need to see that clearly.
    min_lag_to_analyze = MIN(finer_grid[to_analyze])
    max_lag_to_analyze = MAX(finer_grid[to_analyze])
    zoom_range = max_lag_to_analyze - min_lag_to_analyze

  LOADCT, 0, /SILENT
  ;Trying to plot on white backgrounds and enable saving procedures
  plot1 = PLOT(finer_grid, finer_ccf, XRANGE = [min_lag_to_analyze-0.1*zoom_range, max_lag_to_analyze+0.1*zoom_range],LINESTYLE=0,/buffer)
  
  LOADCT, 13, /SILENT
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
;make a name for plotting the spectrum and the table
short_name = STRMID(apStarFile,15,23)

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
   PRINTF,1,FORMAT = '(A20, A20,3X,F6.2, 3X,F6.2,3X, F10.2,3X,F6.2,3X,F6.2,3X,F10.2 ,3X,F10.2, 3X,F10.2, 3X, $
   F10.2,3X,F10.2,3X,F8.2,3X,F8.2,3X,F10.2)', STRTRIM(apogeeID), STRTRIM(visit_ID),R_value, $;, STRTRIM(threshold), $
   bisector_x_range,x_median,bisector_xmean,diff_in_median_mean,x_positions[0],x_positions[1],x_positions[2] ,x_positions[3], $
   x_positions[4], FWHM, height_ccf, threshold

   PRINTF,3, FORMAT = '(A10,3X,A20, 3X,F10.2 ,3X,F10.2, 3X,F10.2, 3X, F10.2,3X,F10.2,3X,F10.2)',STRTRIM(Field),STRTRIM(apogeeID), $
    x_positions[0],x_positions[1],x_positions[2] ,x_positions[3] ,x_positions[4], bisector_x_range

ENDFOR

PRINTF,2, FORMAT = '(A10,3X,A20,3X,F10.2,2X,I5 ,3X , F10.2, 2X,I5, 3X ,F10.2,3X,F10.2)' ,STRTRIM(field), STRTRIM(apogeeID), $
STRTRIM(Rmin), STRTRIM(maxXR), MaxWidth, XR_average


LOADCT, 0, /SILENT
PLOT, [-9999.,-9999.],[0,1], XRANGE = [min_x, max_x], YRANGE = [0,1], /XSTY, /YSTY

FOR i=0, n_visits-1 DO BEGIN
   ;-> -> OPLOT, bisectors[i,0,*], bisectors[i,1,*]
   ; Trying to plot with white bkacgrounds anf saving features:
   plot4 = PLOT(bisectors[i,0,*], bisectors[i,1,*], 'r', /overplot, XTITLE='CCF Lag', YTITLE='$\wedge{CCF Units}$', $
   TITLE='APOGEE ID:' + STRTRIM(string(field)) +'/'+ STRTRIM(string(apogeeID)) $
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

n_spectra = N_ELEMENTS(path)
;print, 'n_spectra: ', n_spectra

CLOSE, 1
OPENW, 1,'Bisector/' + output+ 'Bisector_Stats.lis'
; Missing piece -> 'Threshold Value',
PRINTF,1,'2MASS_ID', 'Visit_ID','R',' x_range ', ' x_median', ' x_mean ','  Diff_in_x', $
      '  x@y=0.1', '  x@y=0.3', '  x@y=0.5', '  x@y=0.7', '  x@y=0.9','FWHM','Height','  Thresh ' , $
      FORMAT = '(A20,A20,A10,A10,A10,A10,A11,A13,A13,A13,A13,A13,A13, A10,A10,A10,A10)'

CLOSE,2
 OPENW,2, 'Bisector/'+output+'Rmin_XR_Max.lis'
 PRINTF,2,'Field','2MASS_ID','R_min', 'XR_Max','Max_Width','XR_average', FORMAT = '(A10,A20,A20,A20,A20,A13)'
 CLOSE,3
 OPENW,3, 'Bisector/'+ output+ 'Bisector_pts.lis'
 PRINTF,3,'Location ID','2MASS_ID', '  x@y=0.1', '  x@y=0.3', '  x@y=0.5', '  x@y=0.7', '  x@y=0.9','x_range', $
 FORMAT ='(A12,A20,A13,A13,A13,A13,A13,A13)'
    
FOR i=0L, n_spectra-1 DO BEGIN
   IF field[i] NE 1 THEN APOGEE_SPEC_BINARY_CLASSIFY, path[i], twoMASSID[i], field[i]
   
ENDFOR

CLOSE,1
CLOSE,2
CLOSE,3
END

