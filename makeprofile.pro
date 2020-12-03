; generates a normalized solar p-mode envelope with degrees ell=0-2 with mean 
; mode rations and threshold amplitude

function makevirgoprofile,noorders,thresh

rdfloat,'VIRGO_low_profile.dat',x,y,/double

y = sqrt(y)	; convert profile to amplitudes
y = y-min(y)
y = y/max(y)

usemax = where(y eq max(y))
dev = 3100.-x[usemax[0]]
x = x+dev

delnu = 135.
dnu = 11.

sp1 = (delnu/2.)-dnu/2.

nofreqs = (noorders*3.)
freqs = dblarr(nofreqs)

freqs[fix(nofreqs/2)] = 3100.

start = fix(nofreqs/2.)+1
k = start-2
for j=start,nofreqs-1,3 do begin
	freqs[j] = freqs[j-1]+sp1
	if (j+1 lt nofreqs) then freqs[j+1] = freqs[j] + dnu
	if (j+2 lt nofreqs) then freqs[j+2] = freqs[j+1] + sp1
	freqs[k] = freqs[k+1] - sp1
	if (k-1 ge 0.) then freqs[k-1] = freqs[k] - dnu
	if (k-2 ge 0.) then freqs[k-2] = freqs[k-1] - sp1
	k = k-3
endfor

amps = dblarr(nofreqs)

for j=0.,nofreqs-1 do begin
	use = where(abs(freqs[j]-x) eq min(abs(freqs[j]-x)))
	amps[j] = y[use[0]]
endfor

ell1 = [fix(nofreqs/2) + findgen(fix(nofreqs/2)/3)*3,fix(nofreqs/2) - 3.-findgen(fix(nofreqs/2)/3)*3]
ell2 = [fix(nofreqs/2)+1 + findgen(fix(nofreqs/2)/3)*3,fix(nofreqs/2)+1 - 3.-findgen(fix(nofreqs/2)/3)*3]
ell0 = [fix(nofreqs/2)+2 + findgen(fix(nofreqs/2)/3)*3,fix(nofreqs/2)+2 - 3.-findgen(fix(nofreqs/2)/3)*3]

amps[ell0] /= 1.25
amps[ell2] /= 1.25/0.75
print,mean(amps[ell1])/mean(amps[ell0])
print,mean(amps[ell2])/mean(amps[ell0])

usel0 = where(amps[ell0] gt thresh)
usel1 = where(amps[ell1] gt thresh)
usel2 = where(amps[ell2] gt thresh)

plot,x,y
oplot,freqs[ell0[usel0]],amps[ell0[usel0]],psym=6,color=255
oplot,freqs[ell1[usel1]],amps[ell1[usel1]],psym=5,color=60000
oplot,freqs[ell2[usel2]],amps[ell2[usel2]],psym=4,color=90000000

freqs = [freqs[ell0[usel0]],freqs[ell1[usel1]],freqs[ell2[usel2]]]
amps = [amps[ell0[usel0]],amps[ell1[usel1]],amps[ell2[usel2]]]

sorted = sort(freqs)

res = dblarr(2,n_elements(freqs))
res[0,*] = freqs[sorted]
res[1,*] = amps[sorted]

return,res
end


function makebisonprofile,file,rotperiod,inclination

rdfloat,'VIRGO_low_profile.dat',x,y,/double

y = sqrt(y)	; convert profile to amplitudes
y = y-min(y)
y = y/max(y)

usemax = where(y eq max(y))
dev = 3050.-x[usemax[0]]
x = x+dev

readcol,file,ell,n,freqs,errs;,skipline=2

um=where(ell eq 1)
freqs=freqs[um]

nf=n_elements(freqs)
freqs=[freqs,freqs+(1./rotperiod/0.0864),freqs-(1./rotperiod/0.0864)]
em=[indgen(nf)*0,indgen(nf)*0+1,indgen(nf)*0-1]

nofreqs = n_elements(freqs)
amps = replicate(1.,n_elements(freqs))

for j=0.,nofreqs-1 do begin
	use = where(abs(freqs[j]-x) eq min(abs(freqs[j]-x)))
	amps[j] = y[use[0]]
	if (ell[j] eq 0) then amps[j] /= 1.25
	if (ell[j] eq 2) then amps[j] /= 1.25/0.75
	if (em[j] eq 0) then amps[j]=amps[j]*cos(inclination)*cos(inclination)
	if (em[j] ne 0) then amps[j]=amps[j]*0.5*sin(inclination)*sin(inclination)
endfor

ell0 = where(ell eq 0)
ell1 = where(ell eq 1)
ell2 = where(ell eq 2)

plot,x,y
oplot,freqs[ell0],amps[ell0],psym=6,color=255
oplot,freqs[ell1],amps[ell1],psym=5,color=60000
oplot,freqs[ell2],amps[ell2],psym=4,color=90000000

res = dblarr(2,n_elements(freqs))
res[0,*] = freqs
res[1,*] = amps

return,res
end



