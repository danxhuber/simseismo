; various smoothing functions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; boxcar mean
function mean_smooth,y,step
n = n_elements(y)
s = dblarr(n)
for j=0.,n-1 do s[j] = mean(y[(j-step/2)>0:(j+step/2)<n-1])
return,s
end

; boxcar mean with independent datapoints
function mean_smooth_ind,x,y,width
step = width-1
n = n_elements(y)
fs = dblarr(n)
s = dblarr(n)
sig = dblarr(n)
j=0.
while (j+step lt n-1) do begin
    fs[j] = mean(x[j:j+step])
    s[j] = mean(y[j:j+step])
    sig[j] = stddev(y[j:j+step])/sqrt(width)
    j += step+1
endwhile
fsn=fs[where(fs ne 0)]
sn=s[where(s ne 0)]
sign=sig[where(s ne 0)]
use = where(sign eq 0)
if (use[0] ne -1) then sign[use] = median(sign)
return,[[fsn],[sn],[sign]]
end


; boxcar mean with independent datapoints
function mean_smooth_ind_log,x,y,width
x = alog10(x)
step = width-1
n = n_elements(y)
fs = dblarr(n)
s = dblarr(n)
sig = dblarr(n)
j=0.
while (j+step lt n-1) do begin
    fs[j] = mean(x[j:j+step])
    s[j] = mean(y[j:j+step])
    sig[j] = stddev(y[j:j+step])/sqrt(width)
    j += step+1
endwhile
fs=fs[where(fs ne 0)]
s=s[where(s ne 0)]
sig=sig[where(sig ne 0)]
x = 10^(x)
fs = 10^(fs)
return,[[fs],[s],[sig]]
end



; boxcar median
function median_smooth,y,step
n = n_elements(y)
s = dblarr(n)
for j=0.,n-1 do s[j] = median(y[(j-step/2)>0:(j+step/2)<n-1])
return,s
end

; mean or median boxcar in log
function log_smooth,x,y,step,median=median
n=n_elements(x)
if (x[0] eq 0) then begin
	xs = x[1:n-1]
	ys = y[1:n-1]
endif else begin
	xs = x
	ys = y
endelse
xs = alog10(xs)
xs = xs-min(xs)
n=n_elements(xs)
s = dblarr(n)
if (keyword_set(median)) then for j=0.,n-1 do s[j] = median(ys[where((xs[j]-xs)>0 le step/2 and $
(xs-xs[j])>0 le step/2)]) $
else for j=0.,n-1 do s[j] = mean(ys[where((xs[j]-xs)>0 le step/2 and (xs-xs[j])>0 le step/2)])
return,s
end

; FFT gauss smooth by Dennis
function gauss_smooth, arr, fwhm, FFT=fft, SILENT=silent
;PURPOSE: Smooth an array of equidistant elements with a moving gauss
;         (instead of a boxcar). Note this program's FFT approach avoids
;         edge effects by gluing a "mirror" of arr on both ends of arr
;         before taking fft. This is however, not the Nurep apporach to
;         deal with edge effects in fft convolution. 
;INPUT  : 1: name of array to be smoothed.
;         2: fwhm of gaussian filter (in units of #elements in array).
;         3: set keyword FFT if smoothing should be done in fouier space
;            (which is much faster).
;OUTPUT : 1: smoothed array.
;CALL BY: e.g. gauss_smooth, mag, 10, magsmooth, /FFT
;
;
;14/11-07 by Dennis Stello 
;****************************************************************************************
ON_ERROR,2

n    = 2*N_ELEMENTS(arr)
sig  = fwhm/SQRT(8*alog(2))

IF keyword_set( FFT ) THEN BEGIN
 IF keyword_set( SILENT ) EQ 0 THEN print,'     ->FFT approach'
 w   = n*GAUSSIAN( FINDGEN(n), [1,n/2.,sig])/TOTAL(GAUSSIAN( FINDGEN(n), [1,n/2.,sig]))
 rev = REVERSE(arr)
 out = SHIFT(FFT(FFT([rev(n/4:n/2-1),arr,rev(0:n/4-1)],-1)*FFT(w,-1),1),LONG(n/2.))
 out = out(n/4:n*3/4)
 if (n/2. mod 2. eq 0.) then out = out[0.:n_elements(out)-2.] else out = out[1.:n_elements(out)-1.]
ENDIF ELSE BEGIN
 IF keyword_set( SILENT ) EQ 0 THEN print,'     ->Weighted mean approach'
 n   = N_ELEMENTS(arr)
 out = arr
 FOR i=0L,n-1 DO BEGIN
  w      = GAUSSIAN( FINDGEN(n), [1,i,sig])
  out(i) = TOTAL(arr*w)/TOTAL(w)
 ENDFOR
ENDELSE

return,out

END


function lorentzian,n,f,eta
return,1./(1.+4.*(n-f)^2/eta^2)
end

function lor_smooth, arr, fwhm, out, FFT=fft, SILENT=silent

ON_ERROR,2

n    = 2*N_ELEMENTS(arr)

IF keyword_set( FFT ) THEN BEGIN
 IF keyword_set( SILENT ) EQ 0 THEN print,'     ->FFT approach'
 w = n*lorentzian(findgen(n),n/2.,fwhm)/total(lorentzian(findgen(n),n/2.,fwhm))
 rev = REVERSE(arr)
 out = SHIFT(FFT(FFT([rev(n/4:n/2-1),arr,rev(0:n/4-1)],-1)*FFT(w,-1),1),LONG(n/2.))
 out = out(n/4:n*3/4)
 if (n/2. mod 2. eq 0.) then out = out[0.:n_elements(out)-2.] else out = out[1.:n_elements(out)-1.]
ENDIF ELSE BEGIN
 IF keyword_set( SILENT ) EQ 0 THEN print,'     ->Weighted mean approach'
 n   = N_ELEMENTS(arr)
 out = arr
 FOR i=0L,n-1 DO BEGIN
  w = lorentzian(findgen(n),i,fwhm)
  out(i) = TOTAL(arr*w)/TOTAL(w)
 ENDFOR
ENDELSE

return,real_part(out)

END

; boxcar smooth over equidistant time
function smooth_time,x,y,int
smarr = dblarr(n_elements(x))
for j=0.,n_elements(x)-1 do smarr[j] = mean(y[where(x gt x[j]-int and x le x[j]+int)])
return,smarr
end






