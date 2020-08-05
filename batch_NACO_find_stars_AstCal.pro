pro batch_NACO_find_stars_AstCal, pathr, fwhm

file = file_search(pathr+'AstrometricField*.fits', count=nfiles)
roundlim = [-0.35,0.35]
fact = 0.5

for xx=0,nfiles-1 do begin

  im = mrdfits(file[xx],0,/silent)
  imo = im
  pos = strpos(file[xx], '/', /reverse_search)
  out = strmid(file[xx],pos+1, strlen(file[xx])-pos-6)

  ;mmm, im, skymod, sigma, skew
  im = im-median(im);skymod

  fn = pathr+out+'_list.dat'
;   find, im, dx, dy, flux, sharp, roundness, 0.5*stddev(im), 2.*fwhm, [-1.,1.], [0.2,1.0], print=fn
  find, im, dx, dy, flux, sharp, roundness, fact*stddev(im), 2.*fwhm, roundlim, [0.2,1.0], print=fn

  if (xx eq 0) then begin

    sz = n_elements(im[*,0])

    if (sz gt 800) then wfact = 1. else wfact = 1.5
    window, 0, xsize=wfact*sz, ysize=wfact*sz
    device, cursor_standard=2

    scale = [-20,50]
    cgimage, im, minvalue=scale[0], maxvalue=scale[1], stretch=1, /axis

    quest = ''
    repeat begin

      read, 'Image scale OK? (y/n): ', quest

      if (quest eq 'n') then begin

	scale = dblarr(2)
	read, 'Enter new low/high values: ', scale
	cgimage, im, minvalue=scale[0], maxvalue=scale[1], stretch=1, /axis

      endif

    endrep until (quest eq 'y')

  endif


  window, 2, xs=1000, ys=1000	;this is only for comparison
  cgimage, im, /axes, minvalue=scale[0], maxvalue=scale[1]
  window, 0, xs=1000, ys=1000
  cgimage, im, /axes, minvalue=scale[0], maxvalue=scale[1]
  oplot, dx, dy, psym=sym(11), color=cgcolor('green')


  if (xx eq 0) then begin

    quest = ''
    repeat begin

      print, ''
      read, 'Detections OK? y/n: ', quest
      if (quest eq 'n') then begin

	read, 'Enter new value (usually between 0.5 and 1.2): ', fact

	find, im, dx, dy, flux, sharp, roundness, fact*stddev(im), 2.*fwhm, roundlim, [0.2,1.0], print=fn

	window, 0, xs=1000, ys=1000
	cgimage, im, /axes, minvalue=scale[0], maxvalue=scale[1]
	oplot, dx, dy, psym=sym(11), color=cgcolor('green')

      endif

    endrep until (quest eq 'y')

  endif

;   oplot, dx, dy, psym=sym(11), color=cgcolor('red')
; 
;   !mouse.button = 0
;   print, ''
;   print, 'Select bad detections'
; 
;   while !mouse.button ne 4 do begin
; 
;     cursor, x, y, 3, /data
; 
;     if (!mouse.button eq 1) then begin
; 
;       !mouse.button = 0
; 
;       dist = sqrt((dx-x)^2+(dy-y)^2)
;       dum = min(dist, idxmin)
; 
;       dx = arrdelete(dx, at=idxmin, length=1, /overwrite)
;       dy = arrdelete(dy, at=idxmin, length=1, /overwrite)
;       flux = arrdelete(flux, at=idxmin, length=1, /overwrite)
;       sharp = arrdelete(sharp, at=idxmin, length=1, /overwrite)
;       roundness = arrdelete(roundness, at=idxmin, length=1, /overwrite)
; 
;       window, 0, xsize=1000, ysize=1000 ;xsize=wfact*sz, ysize=wfact*sz
;       cgimage, im, minvalue=scale[0], maxvalue=scale[1], /axis
;       oplot, dx, dy, psym=sym(11), color=cgcolor('red')
; 
;     endif
; 
;   endwhile


  ;=======================================================================

  ;refine position of found stars
  im = imo
  nstars = n_elements(dx)
  dim = n_elements(im[*,0])-1	;-1 because of the index 
  radiuso = 7
  amp = dblarr(nstars)
  fwhmx = dblarr(nstars)
  fwhmy = dblarr(nstars)
  txc = dblarr(nstars)
  tyc = dblarr(nstars)
  xpos = dblarr(nstars)
  ypos = dblarr(nstars)

  for i=0,nstars-1 do begin

    ;in case of detections close to the image border adjust radius
    tmpmin = min([dx[i]-radiuso, dy[i]-radiuso])
    if (tmpmin lt 0.) then tmpmin = floor(tmpmin) else tmpmin = 0.
    tmpmax = max([dx[i]+radiuso, dy[i]+radiuso])
    if (tmpmax gt dim) then tmpmax = floor(dim-tmpmax) else tmpmax = 0.

    if (tmpmin ne 0.) then rmin = radiuso+tmpmin else rmin = radiuso
    if (tmpmax ne 0.) then rmax = radiuso+tmpmax else rmax = radiuso

    radius = min([rmin,rmax])

    if (radius gt 1) then begin
      cutim = im[round(dx[i])-radius:round(dx[i])+radius, round(dy[i])-radius:round(dy[i])+radius]
      xa = dindgen(2*radius+1)+1. & ya = xa

  ;      weights = cutim
  ;       idx1 = where(cutim eq 0.)	;e.g. beta Pic close to the center which is masked out
  ;       if (idx1[0] ne -1) then begin
  ; 
  ; 	idx2 = array_indices(weights, idx1)
  ; 	for k=0,n_elements(idx2[0,*])-1 do weights[idx2[0,k],idx2[1,k]] = 1.d5
  ; 
  ;       endif
  ; 
  ;       ;work around for having not signum because of IDL < v8.3
  ;       signdum = fltarr(n_elements(cutim[*,0]),n_elements(cutim[0,*]))
  ;       for isign=0,n_elements(cutim[*,0])-1 do begin
  ; 
  ; 	for jsign=0,n_elements(cutim[0,*])-1 do begin
  ; 
  ; 	  signdum[isign,jsign] = sign((reform(cutim[isign,jsign]))[0])
  ; 
  ; 	endfor
  ; 
  ;       endfor
  ;       sign = signdum
  ;       ;sign = signum(cutim)
  ; 
  ;       idx1 = where(sign eq -1.)
  ;       if (idx1[0] ne -1) then begin
  ; 
  ; 	idx2 = array_indices(weights, idx1)
  ; 	for k=0,n_elements(idx2[0,*])-1 do weights[idx2[0,k],idx2[1,k]] = 1.d5
  ; 
  ;       endif
  ;       ;weights = 1./sqrt(weights)

      estimates = [median(cutim), max(cutim), 3.5, 3.5, radius, radius, 0, 0]
      yfit = mpfit2dpeak(cutim, A, xa, ya, /gaussian, estimates=estimates, dof=dof, chisq=chisq, perror=perror, sigma=sigma, /tilt, /quiet, weights=weights)

      amp[i] = A[1]
      fwhmx[i] = A[2]
      fwhmy[i] = A[3]
      txc[i] = A[4]
      tyc[i] = A[5]
      xpos[i] = txc[i]-radius+round(dx[i])
      ypos[i] = tyc[i]-radius+round(dy[i])

    endif else begin

      amp[i] = -1
      fwhmx[i] = -1
      fwhmy[i] = -1
      txc[i] = -1
      tyc[i] = -1

    endelse


;     print, i
;     window, 0
;     cgimage, cutim
;     window, 1
;     cgimage, yfit

  endfor

  ;remove bad fits/outliers

  xa = indgen(nstars)
  resistant_mean, fwhmx, 3, t1, t2, nbad, /double, goodvec=idxg, badvec=idxb
  amp = amp[idxg]
  fwhmx = fwhmx[idxg]
  fwhmy = fwhmy[idxg]
  txc = txc[idxg]
  tyc = tyc[idxg]
  dx = dx[idxg]
  dy = dy[idxg]
  xpos = xpos[idxg]
  ypos = ypos[idxg]
  flux = flux[idxg]
  sharp = sharp[idxg]
  roundness = roundness[idxg]

  xa = indgen(nstars)
  resistant_mean, fwhmy, 3, t1, t2, nbad, /double, goodvec=idxg, badvec=idxb
  amp = amp[idxg]
  fwhmx = fwhmx[idxg]
  fwhmy = fwhmy[idxg]
  txc = txc[idxg]
  tyc = tyc[idxg]
  dx = dx[idxg]
  dy = dy[idxg]
  xpos = xpos[idxg]
  ypos = ypos[idxg]
  flux = flux[idxg]
  sharp = sharp[idxg]
  roundness = roundness[idxg]

  xa = indgen(nstars)
  resistant_mean, txc, 3, t1, t2, nbad, /double, goodvec=idxg, badvec=idxb
  amp = amp[idxg]
  fwhmx = fwhmx[idxg]
  fwhmy = fwhmy[idxg]
  txc = txc[idxg]
  tyc = tyc[idxg]
  dx = dx[idxg]
  dy = dy[idxg]
  xpos = xpos[idxg]
  ypos = ypos[idxg]
  flux = flux[idxg]
  sharp = sharp[idxg]
  roundness = roundness[idxg]

  xa = indgen(nstars)
  resistant_mean, tyc, 3, t1, t2, nbad, /double, goodvec=idxg, badvec=idxb
  amp = amp[idxg]
  fwhmx = fwhmx[idxg]
  fwhmy = fwhmy[idxg]
  txc = txc[idxg]
  tyc = tyc[idxg]
  dx = dx[idxg]
  dy = dy[idxg]
  xpos = xpos[idxg]
  ypos = ypos[idxg]
  flux = flux[idxg]
  sharp = sharp[idxg]
  roundness = roundness[idxg]

  nstars = n_elements(dx)

  ;=======================================================================

  openw, lun, fn, width=2000, /get_lun

  for i=0,nstars-1 do printf, lun, dx[i], dy[i], flux[i], sharp[i], roundness[i], format='(2f14.7,f15.3,3f14.9)'

  close, lun
  free_lun, lun

endfor

end
