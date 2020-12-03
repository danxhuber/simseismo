pro w,no,a,b,c,d,e

if (no eq 2) then begin
    for j=0.,n_elements(a)-1 do begin 
    	printf,1,format='(d15.8,d15.8)',a[j],b[j]
    endfor
endif

if (no eq 3) then begin
    for j=0.,n_elements(a)-1 do begin 
    	printf,1,format='(d15.8,d15.8,d15.8)',a[j],b[j],c[j]
    endfor
endif

if (no eq 4) then begin
    for j=0.,n_elements(a)-1 do begin 
    	printf,1,format='(d15.8,d15.8,d15.8,d15.8)',a[j],b[j],c[j],d[j]
    endfor
endif

if (no eq 5) then begin
    for j=0.,n_elements(a)-1 do begin 
    	printf,1,format='(d15.8,d15.8,d15.8,d15.8,d15.8)',a[j],b[j],c[j],d[j],e[j]
    endfor
endif

end
