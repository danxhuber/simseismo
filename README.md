# Simulate asteroseismic oscillations

This code simulates solar-like dipole oscillations for any given input set of stellar parameters, rotation period and inclination using a stochastically excited and damped harmonic oscillator model and scaling observed properties for the Sun. The core methodology follows Chaplin et al. 1997 (https://ui.adsabs.harvard.edu/abs/1997MNRAS.287...51C/abstract) - please cite this paper if using the code.

### Input (to be edited directly in sim.pro):

; star parameters  
rad=1.0     	  ; stellar radius in Rsun  <br/>
teff=5800.  	  ; effective temperature in K. <br/>
rotperiod=2.	  ; rotation period days. <br/>
inclination=45. ; stellar inclination in degrees. <br/>
mass=1.0        ; stellar mass in Msun (optional, gives more accurate amplitude estimate). <br/>

; timeseries parameters <br/>
length=100.	    ; timeseries length in days <br/>
hrs=24.		      ; duty cycle (hours per day) <br/>
spl=1.	    	  ; sampling in minutes <br/>
rmsnoise=100.0  ; white noise in ppm <br/>


### Output:

ts_noisefree.txt: timeseries without noise added <br/>
ts_noise.txt: timeseries with noise added <br/>
psd_noisefree.txt: power density spectrum without noise added <br/>
psd_noise.txt: power density spectrum with noise added <br/>

### Usage:
idl <br/>
.compile sim <br/>
main <br/>
