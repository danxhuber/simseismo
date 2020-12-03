;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; simulate solar-like dipole modes with a given rotational splitting and inclination
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

@makeprofile.pro
@smooth_funct.pro

pro main

!EXCEPT = 0

COMMON FREQDATA,f,a,w
COMMON TDDATA,x,y,ysynth
COMMON FDDATA,freq,amp,phase,lorentzian
COMMON SETUP,col
COMMON RANDOM,seed  	    ; IMPORTANT FOR INDEPENDENT REALIZATIONS!

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;; Start to edit parameters here ;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; star parameters
rad=1.0    	; Rsun
teff=5800.  	; K
rotperiod=10.	; days
inclination=90. ; stellar inclination in degrees

; optional. mass if known (gives more accurate amplitude estimate)
mass=1.0

; timeseries parameters
length=100.	; timeseries length in days
hrs=24.		; duty cycle (hours per day)
spl=1.	    	; sampling in minutes
rmsnoise=100.0    ; white noise in ppm

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;; End edit parameters here ;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


; estimate amplitude from stellar pars
s=0.886
t=1.89
r=2.0
lum=rad^2*(teff/5777.)^4
amp=(lum^s)/(mass^t * (teff/5777D)^(r-1D) * (teff/5934D)^0.8 )
amp=lum/mass
amp=3.6*amp
rvamp=((amp*(600./550.)*((teff)/5777.)^2D)/20.1)

; estimate numax 
;logg=4.44
;numax=3090.*(10^logg/10^4.44 * (teff/5777.)^(-0.5))
numax=3100.*rad^(-1.85)*(teff/5777.)^0.92

; linewidth from Lund+ 2017
lwidth=1.02*exp((teff-5777.)/436.)
lifetime=1./(lwidth*!DPI)/0.0864    ; mode lifetime\

rvampin=amp
numaxin=numax
teffin=teff

nyq=1./(2.*spl/60./24)/0.0864

print,'------------------------------------------------------'
print,'teff (K):',teffin
print,'mass (solar):',mass
print,'lum (solar):',lum
print,'numax (muHz):',numaxin
print,'length (nights):',length
print,'sampling (minutes):',spl
print,'nyquist (muHz):',nyq
print,'mode liftime (days):',lifetime
print,'stellar inclination (degrees):',inclination
print,'duty cycle (hrs/day):',hrs
print,'amplitude (ppm):',rvampin
print,'precision per cadence (ppm):',rmsnoise
print,'------------------------------------------------------'

ndata=length*hrs*60./spl
sn=rvampin/(rmsnoise/sqrt(ndata))
print,'expected S/N:',sn

rvampin=rvampin+rvampin*1.9

rdfloat,'VIRGO_low_profile.dat',xvirg,yvirg,/double
yvirg = sqrt(yvirg)
yvirg = yvirg-min(yvirg)
yvirg /= max(yvirg)

splt=spl/60./24.
x=dindgen(length/splt)*splt
x=x*86400.*1e-6

; fixed simulation parameters. these should be the same as in the real application
burn = 0.1D		; burn in time in days
lag = [10.,300.]    	; lag range
noisefloor = 5.     	; noisefloor of simulations in power
step = 1.*86400D*1e-6   ; stepsize of co-adding
box = 5.*86400D*1e-6    ; boxsize of co-adding
nopeak = 5. 	    	; number of peaks to consider for delta_nu
noisedet = 8100.    	; frequency above which noise is determined
acrange = 10.	    	; times*delnu_exp to calculate AC

xstart = x[0]
x = (x-x[0]+burn)	; zero normalize in time

n=n_elements(x)
samp = mean(x[1:n-1]-x[0:n-2])
z = (burn/samp)
xburn = findgen(z)*samp
xnew = dblarr(n+n_elements(xburn))
xnew[0:n_elements(xburn)-1]=xburn
xnew[n_elements(xburn):n_elements(xnew)-1]=x
x=xnew
length = x[n_elements(x)-1]-x[0]
fres = 1D/length

xo = x

noise = dblarr(n_elements(x))
nsims = 100.
delnu_m = dblarr(nsims)
sns = dblarr(nsims)

os = 1D
fres = 1D/((length-burn)*os)
print,fres

bison = makebisonprofile('sun_pmodes_use.dat',rotperiod,inclination*!DPI/180.)
freqs_sun = bison[0,*]
amps_sun = bison[1,*]
m = max(amps_sun,pos)
diff = 3050.- freqs_sun[pos]	; difference between numax and max amp!
freqs_sun = freqs_sun - freqs_sun[pos]

maxv = xvirg[where(yvirg eq max(yvirg))]
xvirg = xvirg - maxv[0]



dnu = 134.9 * (numax/3050.)^0.8

; scale the frequencies at zeropoint
f = freqs_sun/(134.9D/dnu)
; same for the difference between numax and max am
diff2 = diff/(134.9D/dnu)
; shift frequencies
f = (f+numax-diff2)

xv = xvirg/(134.9D/dnu)
xv = xv+numax

f = f*2D*!DPI



a = rvampin*amps_sun


plot,f,a,psym=4
l = dblarr(n_elements(f))
l[*] = lifetime*0.0864
w = 1D/(l*!DPI)	

	


pps=rmsnoise
noise=randomn(seed,n_elements(x))*pps

; generate synthetic time series
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
td,length,0,0,50D,0,burn    	    	    	    	
use = where(x ge burn)
x=x[use]
x=x-burn+xstart
yo=ysynth
y = ysynth+noise

x=(x*1e6)/60./60./24.

p = x mod 1.
u=where(p lt hrs/24.)
x=x[u]
y=y[u]
yo=yo[u]

!p.charsize=2
!p.multi=[0,1,3]

plot,x,y,psym=-4,/xs
oplot,x,yo,psym=-4,color=1

openw,1,'ts_noisefree.txt'
w,2,x,yo
close,1

openw,1,'ts_noise.txt'
w,2,x,y
close,1

d=x
fc=yo

os=5

powspec = lnp_test(d,fc,wk1=freq,wk2=amp,ofac=os)
freq = 1000.*freq/86.4
corfac = nyq/max(freq)    
powspec = lnp_test(d,fc,wk1=freq,wk2=amp,ofac=os,hifac=corfac)
freq = 1000.*freq/86.4
bin = freq[1]-freq[0]
amp = 2.*amp*variance(fc)/(total(amp)*bin)
pssm = real_part(gauss_smooth(amp,4.*dnu/(freq[1]-freq[0]),/fft,/silent))	
plot,freq,amp,xrange=[0,8000],/xs
;oplot,freq,pssm,color=1,thick=2

freq_sim=freq
amp_sim=amp

openw,1,'psd_noisefree.txt'
w,2,freq,amp
close,1

print,'actual input amplitude:',sqrt(max(pssm)*dnu/4.09)



d=x
fc=y

powspec = lnp_test(d,fc,wk1=freq,wk2=amp,ofac=os)
freq = 1000.*freq/86.4
corfac = nyq/max(freq)    
powspec = lnp_test(d,fc,wk1=freq,wk2=amp,ofac=os,hifac=corfac)
freq = 1000.*freq/86.4
bin = freq[1]-freq[0]
amp = 2.*amp*variance(fc)/(total(amp)*bin)


openw,1,'psd_noise.txt'
w,2,freq,amp
close,1


!p.charsize=3
!p.multi=[0,1,3]

plot,x,y,xtitle='days',ytitle='ppm',title='timeseries'
oplot,x,yo,color=1;,psym=-4,color=1

plot,freq_sim,amp_sim,xrange=[0,8000],/xs,title='noise-free',xtitle='muhz'

plot,freq,amp,xrange=[0,8000],/xs,/ys,title='simulated',xtitle='muhz'
;oplot,freq,pssm2,color=1,thick=2

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; generate data in time domain using existing dataset

pro tdex,tsfile,freqfile,t_ex,lfreq,ufreq,whitenoise,limita

COMMON FREQDATA,f,a,w
COMMON TDDATA,x,y,ysynth
COMMON FDDATA,freq,amp,phase,lorentzian
COMMON SETUP,col

setup

burn = 100D						; burn in time in days

print,' '
print,'reading dataset and frequencies ...'

rdfloat,tsfile,x,y,/double

xstart = x[0]
x = (x-x[0]+burn)	; zero normalize in time

; sampling for burn-in period is determined through mean of original sampling
; histogram would be better, implement once available for GDL

n=n_elements(x)
samp = mean(x[1:n-1]-x[0:n-2])
z = (burn/samp)
xburn = findgen(z)*samp
xnew = dblarr(n+n_elements(xburn))
xnew[0:n_elements(xburn)-1]=xburn
xnew[n_elements(xburn):n_elements(xnew)-1]=x

x=xnew
length = x[n_elements(x)-1]-x[0]
fres = 1D/length

readfreqs,freqfile,fres,0

print,' '
print,'length of dataset:',length-burn
print,'mean sampling used for burn-in period:',samp
print,'frequency resolution:',1D/(length-burn)

print,' '
print,'calculating synthetic dataset ...'

td,length,0,0,t_ex,limita,burn

use = where(x ge burn)
x=x[use]
x=x-burn+xstart

if (whitenoise ne 0D) then begin
	print,' '
	print,'adding white noise ...'
		noise = dblarr(n_elements(y))
		for j=0.,n_elements(y)-1 do begin
			noise[j] = randomn(seed)*whitenoise
		endfor
		ysynth=ysynth+noise
endif

wset,0
plot,x,ysynth,xs=1,title='synthetic data'
oplot,x,ysynth,color=col

y = y+ysynth

print,' '
print,'calculating DFT ...'
dft,lfreq,ufreq,5D,0

wset,1
plot,freq,amp,title='synthetic spectrum'
if (limita eq 1) then oplot,freq,lorentzian,color=col

print,' '
openw,1,'synthdata_tdex.dat'
for j=0.,n_elements(x)-1 do begin
	printf,1,format='(d15.8,A,d15.8)',x[j],'    ',y[j]
endfor
c

openw,1,'synthspectrum_tdex.dat'
for j=0.,n_elements(freq)-1 do begin;
	printf,1,format='(d15.8,A,d15.8,A,d15.8)',freq[j],'    ',amp[j],'    ',lorentzian[j]
endfor
c

print,' '
print,'done'

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; generate data in time-domain (Chaplin et al. 1997)

pro td,length,ini_d,ini_v,t_ex,limita,burn

COMMON FREQDATA,f,a,w
COMMON TDDATA,x,y,ysynth
COMMON FDDATA,freq,amp,phase,lorentzian
COMMON SETUP,col
COMMON RANDOM,seed  	    	; IMPORTANT FOR INDEPENDENT REALIZATIONS!!!!!!!!!!!

tau =  1./(!DPI*w)			; mode lifetime
damp = 1./tau				; damping constant eta
t_ex = t_ex*1e-6;/(60D*60D*24D)		; kick rate in days

n = n_elements(x)
samp = mean(x[1:n-1]-x[0:n-2])		; mean sampling
nel = FIX(t_ex/samp)+1			; minimum size of re-excitation time array

use = where(x ge burn)

ytmp = dblarr(n_elements(x))		; calculated displacement
ytmpp = dblarr(n_elements(x))		; calculated velocity
ysynth = dblarr(n_elements(use))	; displacement of all frequencies (without burn-in)

xstep = x

for k=0.,n_elements(f)-1 do begin

	;print,'         calculating frequency',strtrim(FIX(k),0)
	sigmaA = sqrt((a[k]^2*(length-burn))/(4D*tau[k]))	; rms of amplitude needed to generate lorentzian of height a
	amp = a[k]
	sq = sqrt(f[k]^2 - damp[k]^2)
	n = 0.
	init = x[0]

	count = 0.

	for j=0.,(n_elements(x))-2 do begin
		;if (x[j]-init gt 3D*t_ex) then goto,skipexc
		if ((j gt 1) && (x[j]-init gt t_ex)) then begin	; re-excitation
			ini_d = ytmp[j-1]
			ini_v = ytmpp[j-1]
			amp = randomn(seed)*a[k];*sigmaA + a[k]
			xstep = x[j:(j+nel)]-x[j-1]
			n = 0.
			init = x[j]
		endif 

		si = Sin(sq*xstep[n])
		co = Cos(sq*xstep[n])
		e = exp(-damp[k]*xstep[n])

		C = (amp+ini_v+ini_d*damp[k])/sq

		ytmp[j] = C * e * si + ini_d * e * co
		ytmpp[j] = (C*e*sq - ini_d*damp[k]*e) * co - (ini_d*e*sq + C*damp[k]*e) * si

		n = n + 1
		
	endfor

if (limita eq 1) then begin
	ysynth += (ytmp[use]/stddev(ytmp[use]))*(sigmaA)	; set stddev of spectrum so the input amplitudes correspond to resuling height of LP
endif else begin
	ysynth += (ytmp[use]/stddev(ytmp[use]))*a[k]
endelse

endfor

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; read frequencies

pro readfreqs,file,fres,type

COMMON FREQDATA,f,a,w
COMMON TDDATA,x,y,ysynth
COMMON FDDATA,freq,amp,phase
	
rdfloat,file,f,a,l,/double
w = 1D/(l*!DPI)			; convert lifetime to width of lorentzian profile

if (type eq 1) then begin

	u = where((f - (FIX(f/fres)*fres)) ne 0)
	if (u[0] ne -1) then begin
		print,' '
		print,'one or more frequencies are not multiples of frequency resolution stepsize ...' 
		print,'readjusting frequencies:'
		print,' '
		for i=0.,n_elements(f)-1 do begin
			fold = f[i]
			f[i] = double(fres*round(f[i]/fres))
			print,fold,' replaced by ',strtrim(f[i],1)
		endfor

	openw,1,file+'_adjusted'
	for j=0.,n_elements(f)-1 do begin;
		printf,1,format='(d15.8,A,d15.8,A,d15.8)',f[j],'    ',a[j],'    ',l[j]
	endfor
	c

endif

endif

end
