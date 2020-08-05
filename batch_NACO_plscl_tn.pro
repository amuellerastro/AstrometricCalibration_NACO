pro batch_NACO_plscl_tn, pathr, field, ptflag, agpmflag

if (field eq '47Tuc') then catepoch = 2002.26d0


;===============================================================================================

if (field eq 'Trapezium') then begin	;no proper motion for this catalog

  offxy = dblarr(2)

  print, 'Reading Catalog...'

  ;catalog with reference coordinayes
  cat = '/home/amueller/work/IDLlibs/AO/AstCal_tet01OriC_Close2012/table.dat'
  readcol, cat, ra5o, dec5o, format='d,d', /silent

  ;-------------------------------------------------------------------------------

  file = file_search(pathr+'AstrometricField*.fits', count=nfiles)

  for xx=0,nfiles-1 do begin

    ra5 = ra5o
    dec5 = dec5o

    pos = strpos(file[xx], '/', /reverse_search)
    out = strmid(file[xx],pos+1, strlen(file[xx])-pos-6)

    ;read in list of detected stars from starfinder
    readcol, pathr+out+'_list.dat', xf, yf, ff, roundf, sharpf, format='d,d,d,d,d', /silent

    ;-------------------------------------------------------------------------------

    nstars = n_elements(xf)

    im = mrdfits(file[xx], 0, hdr, /silent)
    im = im-median(im)
    npix = double(get_eso_keyword(hdr, 'NAXIS1'))
    extast, hdr, astro

    if (agpmflag eq 0) then astro.cd = [[-7.5527800e-06,0.],[0.,7.5527800e-06]]

    jd  = double(get_eso_keyword(hdr, 'MJD-OBS'))+2400000.5

    ;-------------------------------------------------------------------------------

    ;x-y coordinates of cat stars
    xref = dblarr(n_elements(ra5))
    yref = dblarr(n_elements(ra5))

    for i=0,n_elements(ra5)-1 do begin

      ad2xy, ra5[i]*15., dec5[i], astro, t1, t2
      xref[i] = t1
      yref[i] = t2

    endfor

    ;-------------------------------------------------------------------------------

    ;ra, dec coordinates for stars in image
    raf = dblarr(nstars)
    decf = dblarr(nstars)

    for i=0,nstars-1 do begin

      xy2ad, xf[i], yf[i], astro, t1, t2
      raf[i] = t1/15.
      decf[i] = t2

    endfor
    

    ;-------------------------------------------------------------------------------

    qoff = ''
  ;   repeat begin

  ;   window, 0, xs=npix, ys=npix
  ;   plot, xf, yf, psym=2, symsize=2, xr=[0,npix], yr=[0,npix], xst=1, yst=1, xtitle='X [px]', ytitle='Y [px]', charsize=2
  ;   oplot, xref, yref, psym=1, symsize=2, color=cgcolor('red')

    if (xx eq 0) then begin

      window, 0, xs=npix, ys=npix
      scale = [-20,50]
      cgimage, im, /axis, minvalue=scale[0], maxvalue=scale[1], stretch=1;, xtitle='X [px]', ytitle='Y [px]'

      if (xx eq 0) then begin


	quest = ''
	repeat begin

	  read, 'Image scale OK? (y/n): ', quest

	  if (quest eq 'n') then begin

	    print, ''
	    print, 'Current values: ', scale
	    scale = dblarr(2)
	    read, 'Enter new low/high values: ', scale
	    cgimage, im, minvalue=scale[0], maxvalue=scale[1], /axis

	  endif

	endrep until (quest eq 'y')

      endif

    endif

  ; 
  ; 
  ;     read, 'Apply offsets? (y/n): ', qoff
  ;     if (qoff eq 'y') then begin
  ; 
  ;       read, 'X offset: ', xoff
  ;       read, 'Y offset: ', yoff
  ; 
  ;       xref = xref+xoff
  ;       yref = yref+yoff
  ; 
  ;       window, 0, xs=npix, ys=npix
  ;       plot, xf, yf, psym=2, symsize=2,xr=[0,npix], yr=[0,npix], xst=1, yst=1, xtitle='X [px]', ytitle='Y [px]', charsize=2
  ;       oplot, xref, yref, psym=1, symsize=2, color=cgcolor('red')
  ; 
  ;     endif
  ; 
  ; 
  ;   endrep until (qoff eq 'n')


    xref = xref+offxy[0]
    yref = yref+offxy[1]

    repeat begin

      window, 0, xs=npix, ys=npix
      cgimage, im, /axis, minvalue=scale[0], maxvalue=scale[1], xtitle='X [px]', ytitle='Y [px]'
      oplot, xf, yf, psym=sym(11), symsize=2, color=cgcolor('green')
      oplot, xref, yref, psym=sym(4), symsize=2, color=cgcolor('red')

      read, 'Apply offsets? (y/n): ', qoff
      if (qoff eq 'y') then begin

	read, 'X offset: ', xoff
	read, 'Y offset: ', yoff

	offxy[0] = offxy[0]+xoff
	offxy[1] = offxy[1]+yoff

	xref = xref+xoff
	yref = yref+yoff

	window, 0, xs=npix, ys=npix
	oplot, xf, yf, psym=sym(11), symsize=2, color=cgcolor('green')
	oplot, xref, yref, psym=1, symsize=2, color=cgcolor('red')

      endif

    endrep until (qoff eq 'n')

    window, 0, xs=npix, ys=npix
    cgimage, im, /axis, minvalue=scale[0], maxvalue=scale[1], xtitle='X [px]', ytitle='Y [px]'
    oplot, xf, yf, psym=sym(11), symsize=2, color=cgcolor('green')
    oplot, xref, yref, psym=sym(4), symsize=2, color=cgcolor('red')


    ;find match between the two list (measured and catalog)

    idxsort = intarr(nstars)
    dmin = dblarr(nstars)
    for i=0,nstars-1 do begin

      tmp = sqrt((xf[i]-xref)^2.+(yf[i]-yref)^2.)

      dum = min(tmp, idxmin)
      idxsort[i] = idxmin
      dmin[i] = dum

    endfor

  ;   xf = xf[idxsort]
  ;   yf = yf[idxsort]
    xref = xref[idxsort]
    yref = yref[idxsort]
    ra5 = ra5[idxsort]
    dec5 = dec5[idxsort]
  ;   raf = raf[idxsort]
  ;   decf = decf[idxsort]

    ;check if stars really match using distance
    idx = where(dmin gt median(dmin)+stddev(dmin))
    if (idx[0] ne -1) then begin

      xf[idx] = -99.
      yf[idx] = -99.
      xref[idx] = -99.
      yref[idx] = -99.
      ra5[idx] = -99.
      dec5[idx] = -99.
      raf[idx] = -99.
      decf[idx] = -99.

    endif

    idx0 = where(xf gt -1.)
    if (idx0[0] ne -1) then begin

      xf = xf[idx0]
      yf = yf[idx0]
      xref = xref[idx0]
      yref = yref[idx0]
      ra5 = ra5[idx0]
      dec5 = dec5[idx0]
      raf = raf[idx0]
      decf = decf[idx0]

    endif

;     window, 1, xs=npix, ys=npix, title='Selected Stars'
;     cgimage, im, /axis, minvalue=scale[0], maxvalue=scale[1], xtitle='X [px]', ytitle='Y [px]'
;     oplot, xf, yf, psym=sym(11), symsize=2, color=cgcolor('green')
;     oplot, xref, yref, psym=sym(4), symsize=2, color=cgcolor('red')
; 
;     print, ''
;     print, '***************************************'
;     print, 'Check plot if selection of stars is OK.'
;     print, '***************************************'
;     hak

    nuse = n_elements(xf)

    d5_deg = dblarr(nuse, nuse)
    ;d5_px = dblarr(nuse, nuse)
    df_px = dblarr(nuse, nuse)
    df_deg = dblarr(nuse, nuse)
    pa5 = dblarr(nuse, nuse)
    paf = dblarr(nuse, nuse)

    for i=0,nuse-1 do begin

      ;distance

      gcirc, 1, ra5[i], dec5[i], ra5, dec5, d
      d5_deg[i,*] = d	;arcsec

      gcirc, 1, raf[i], decf[i], raf, decf, d
      df_deg[i,*] = d

      ;d5_px[i,*] = sqrt((xref[i]-xref)^2.+(yref[i]-yref)^2.)

      df_px[i,*] = sqrt((xf[i]-xf)^2.+(yf[i]-yf)^2.)


      ;position angle
      ;computing it with star pairs - gives higher number than PA from reference coordinate, e.g. 5 stars available -> 10 PAs
      for j=0,nuse-1 do begin

  ;       tmpx = xref[i]-xref[j]
  ;       tmpy = yref[i]-yref[j]
  ;       tmp = cv_coord(from_rect=[tmpx, tmpy], /to_polar, /degrees)
  ;       if (tmp[0] lt 0.) then tmp[0] = tmp[0]+360.
  ;       pa5[i,j] = tmp[0]
  ; 
  ;       tmpx = xf[i]-xf[j]
  ;       tmpy = yf[i]-yf[j]
  ;       tmp = cv_coord(from_rect=[tmpx, tmpy], /to_polar, /degrees)
  ;       if (tmp[0] lt 0.) then tmp[0] = tmp[0]+360.
  ;       paf[i,j] = tmp[0]

	posang, 1, ra5[i], dec5[i], ra5[j], dec5[j], angle5
	if (angle5 lt 0.) then angle5 = angle5+360.
	pa5[i,j] = angle5

	posang, 1, raf[i], decf[i], raf[j], decf[j], anglef
	if (anglef lt 0.) then anglef = anglef+360.
	paf[i,j] = anglef

      endfor

    endfor

    ;convert from catalog dist and pixel into a platescale

    window, 1, xs=1000, ys=800, title='Computed Values'
    !p.multi=[0,1,2]
  ; 
  ;   ;plate scale cat
  ;   tmp = (d5_deg/d5_px)[*]
  ;   idxnan = where(finite(tmp) eq 1)
  ;   tmp2 = tmp[idxnan]
  ;   tmp3 = tmp2[uniq(tmp2, sort(tmp2))]
  ;   resistant_mean, tmp3*1.d3, 5., plscl5, eplscl5, nbad, /double, GoodVec=good
  ;   print, 'Number of bad plscl catalog values: ', nbad
  ;   ;plscl5 = median(tmp3*1.d3)
  ;   eplscl5 = stddev(tmp3[good]*1.d3)
  ;   plot, tmp3*1.d3, /yn, psym=2, charsize=2, ytitle='PLSCL CAT'
  ;   oplot, tmp3[good]*1.d3, color=cgcolor('red')

    ;plate scale
    tmp = (d5_deg/df_px)[*]
    idxun = uniq(tmp, sort(tmp))
    tmp2 = tmp[idxun]
    tmp2 = tmp2[where(tmp2 gt 0.)]
    resistant_mean, tmp2*1.d3, 5., plscl, eplscl, nbad, /double, GoodVec=good
    print, 'Number of bad plscl values        : ', nbad
    ;plscl = median(tmp3*1.d3)
    eplscl = stddev(tmp2[good]*1.d3)
    xaxis = indgen(n_elements(tmp2))
    plot, xaxis, tmp2*1.d3, /yn, psym=2, charsize=2, ytitle='PLSCL'
    oplot, xaxis[good], tmp2[good]*1.d3, color=cgcolor('red'), psym=2

    ;true north
    tmp = (paf-pa5)[*]
    tmp2 = tmp[idxun]
    tmp2 = tmp2[where(tmp2 ne 0.)]
    if (median(tmp2) lt 0.) then begin
      if (max(tmp2) gt 360.) then begin

	idxmax = where(tmp2 gt 360.)
	tmp2[idxmax] = tmp2[idxmax]-360.

      endif
    endif
    resistant_mean, tmp2, 5., tn, etn, nbad, /double, GoodVec=good
    print, 'Number of bad TN values           : ', nbad
    ;tn = median(tmp)
    etn = stddev(tmp2[good])
    xaxis = indgen(n_elements(tmp2))
    plot, xaxis, tmp2, /yn, psym=2, charsize=2, ytitle='TN'
    oplot, xaxis[good], tmp2[good], color=cgcolor('red'), psym=2

    !p.multi=[0,1,0]

    print, ''
    print, 'plate scale [mas]: ', sigfig(plscl, 5), ' +/- ', sigfig(eplscl, 2)
    print, 'true north [deg]: ', sigfig(tn, 3), ' +/- ', sigfig(etn, 2)
    print, ''

    if (eplscl gt 0.15) then begin

      print, '!!!!!!!!!!!!!!!!!!!!!'
      print, '!!! Check values! !!!'
      print, '!!!!!!!!!!!!!!!!!!!!!'

    endif

    mode = agpmflag
    fn = pathr+out+'_AstrometricCorrection.sav'
    save, plscl, eplscl, tn, etn, jd, mode, ptflag, field, nuse, offxy, filename=fn

  endfor

endif

;===============================================================================================
;===============================================================================================

if (field eq '47Tuc') then begin

  offxy = dblarr(2)

  print, 'Reading Catalog...'

  cat = '/home/amueller/work/IDLlibs/AO/AstCal_47Tuc_McLaughlin2006/table4.dat'

  readcol, cat, ra_c, dec_c, hstmag, ra4, dec4, ra_h, ra_m, ra_s, dec_d, dec_m, dec_s, magmin, flag, comp, frame4, xm, ym, xr, yr, format='d,d,d,d,d,d,d,d,d,d,d,d,d,d,a,d,d,d,d', /silent

  ;   cat = '/home/amueller/work/IDLlibs/AO/AstCal_47Tuc_McLaughlin2006/table5a.dat'
  ; 
  ;   readcol, cat, frame5,ra_c5,dec_c5,vmag5,umag5,hstmag5,ra5,dec5,pmra,epmra,chi2ra,pchi2ra,pmdec,epmdec,chi2dec,pchi2dec, format='a,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d', /silent
  ; 
  ; ra4 = ra5
  ; dec4 = dec5
  ; hstmag = hstmag5

  idx = where(flag eq 1)
  hstmag = hstmag[idx]
  ra4 = ra4[idx]
  dec4 = dec4[idx]
  comp = comp[idx]

  idx = where(hstmag lt 16.)
  hstmag = hstmag[idx]
  ra4 = ra4[idx]
  dec4 = dec4[idx]
  comp = comp[idx]

  ra5o = ra4
  dec5o = dec4
  hstmago = hstmag

  ;=======================================================================================

  file = file_search(pathr+'AstrometricField*.fits', count=nfiles)

  for xx=0,nfiles-1 do begin

    ra5 = ra5o
    dec5 = dec5o
    hstmag = hstmago

    pos = strpos(file[xx], '/', /reverse_search)
    out = strmid(file[xx],pos+1, strlen(file[xx])-pos-6)

    ;read in list of detected stars from starfinder
    readcol, pathr+out+'_list.dat', xf, yf, ff, roundf, sharpf, format='d,d,d,d,d', /silent

    nstars = n_elements(xf)

    ;-------------------------------------------------------------------------------

    im = mrdfits(file[xx], 0, hdr, /silent)
    im = im-median(im)
    npix = double(get_eso_keyword(hdr, 'NAXIS1'))
    extast, hdr, astro

;===================================================================================
; astro.cd = [[-7.5527800e-06,0.],[0.,7.5527800e-06]]
;===================================================================================

    jd  = double(get_eso_keyword(hdr, 'MJD-OBS'))+2400000.5

    ;-------------------------------------------------------------------------------

    ;x-y coordinates of cat stars
    xref = dblarr(n_elements(ra5))
    yref = dblarr(n_elements(ra5))

    for i=0,n_elements(ra5)-1 do begin

      ad2xy, ra5[i]*15., dec5[i], astro, t1, t2
      xref[i] = t1
      yref[i] = t2

    endfor

    ;-------------------------------------------------------------------------------

    ;ra, dec coordinates for stars in image
    raf = dblarr(nstars)
    decf = dblarr(nstars)

    for i=0,nstars-1 do begin

      xy2ad, xf[i], yf[i], astro, t1, t2
      raf[i] = t1/15.
      decf[i] = t2

    endfor  

    ;-------------------------------------------------------------------------------

    magbin = [10,11,12,13,14,15,16,17,18]
    plotsz = linspace(0.5,2.5,n_elements(magbin)-1)
    qoff = ''

    window, 0, xs=npix, ys=npix
    scale = [-20,50]
    cgimage, im, /axis, minvalue=scale[0], maxvalue=scale[1], xtitle='X [px]', ytitle='Y [px]'

    if (xx eq 0) then begin

      quest = ''
      repeat begin

	read, 'Image scale OK? (y/n): ', quest

	if (quest eq 'n') then begin

	  scale = dblarr(2)
	  read, 'Enter new low/high values: ', scale
	  cgimage, im, minvalue=scale[0], maxvalue=scale[1], /axis

	endif

      endrep until (quest eq 'y')

    endif

    ;-------------------------------------------------------------------------------

;     xaref = abs(floor(min(xref)))+abs(ceil(max(xref)))
;     yaref = abs(floor(min(yref)))+abs(ceil(max(yref)))
;     xaf = abs(floor(min(xf)))+abs(ceil(max(xf)))
;     yaf = abs(floor(min(yf)))+abs(ceil(max(yf)))
;     imref = fltarr(xaref,yaref)
;     ; imf = imref
;     xcoordref = xref+abs(floor(min(xref)))
;     ycoordref = yref+abs(floor(min(yref)))
;     szpsf = 25.
;     rszpsf = (szpsf-1.)/2.
;     psf = psf_Gaussian(npixel=szpsf, fwhm=10., /normalize)
;     for i=0L,n_elements(xref)-1 do begin
; 
;       t1 = round(xcoordref[i])
;       t2 = round(ycoordref[i])
; 
;       if (t1-rszpsf gt 0 and t2-rszpsf gt 0 and t1+rszpsf lt xaref and t2+rszpsf lt yaref) then imref[t1-rszpsf:t1+rszpsf,t2-rszpsf:t2+rszpsf] = psf
; 
;     endfor
; 
;     ; for i=0L,n_elements(xf)-1 do begin
;     ; 
;     ;   t1 = round(xf[i])+abs(floor(min(xref))) 
;     ;   t2 = round(yf[i])+abs(floor(min(yref)))
;     ;   if (t1-7 gt 0 and t2-7 gt 0 and t1+7 lt xaf+abs(floor(min(xref))) and t2+7 lt yaf+abs(floor(min(yref)))) then imf[t1-7:t1+7,t2-7:t2+7] = psf
;     ; 
;     ; endfor
; 
;     imf = fltarr(xaf,yaf)
;     for i=0L,n_elements(xf)-1 do begin
; 
;       t1 = round(xf[i])
;       t2 = round(yf[i])
;       if (t1-rszpsf gt 0 and t2-rszpsf gt 0 and t1+rszpsf lt xaf and t2+rszpsf lt yaf) then imf[t1-rszpsf:t1+rszpsf,t2-rszpsf:t2+rszpsf] = psf
; 
;     endfor
; 
;     szref = double(size(imref))
;     szf = double(size(imf))
;     reduc = 7.
;     imref = congrid(imref, szref[1]/reduc, szref[2]/reduc, cubic=-0.5)
;     imf = congrid(imf, szf[1]/reduc, szf[2]/reduc, cubic=-0.5)
; 
;     nx = n_elements(imref[*,0])-n_elements(imf[*,0])+1
;     ny = n_elements(imref[0,*])-n_elements(imf[0,*])+1
; 
;     mean_imref = mean(imref)
;     im_ov_B = imf - mean(imf)
;     totBB = total( im_ov_B * im_ov_B )
; 
;     correl_mat = dblarr(nx,ny)
; 
;     for xdir=0L,nx-1,2 do begin
; 
;       for ydir=0L,ny-1,2 do begin
; 
; 	tmp = imref[xdir:n_elements(imf[*,0])-1+xdir, ydir:n_elements(imf[0,*])-1+ydir]
; 
; 	im_ov_A = tmp - mean_imref
; 	totAA = total( im_ov_A * im_ov_A )
; 	if (totAA EQ 0) or (totBB EQ 0) then $
; 	  correl_mat[xdir, ydir] = 0.0 else $
; 	  correl_mat[xdir, ydir] = total( im_ov_A * im_ov_B ) / sqrt( totAA * totBB )
; 
;       endfor
; 
;       proceeding_text,loop=nx, i=xdir, prompt='> 2D Correlation   '+string(xdir+1,form='(I4)')
; 
;     endfor

    ;-------------------------------------------------------------------------------

    xref = xref+offxy[0]
    yref = yref+offxy[1]

    repeat begin

      window, 0, xs=npix, ys=npix
      cgimage, im, /axis, minvalue=scale[0], maxvalue=scale[1], xtitle='X [px]', ytitle='Y [px]'
      oplot, xf, yf, psym=sym(11), symsize=2, color=cgcolor('green')

      for i=0,n_elements(magbin)-2 do begin

	idx = where(hstmag ge magbin[i] and hstmag lt magbin[i+1])
	if (idx[0] ne -1) then oplot, xref[idx], yref[idx], psym=sym(4), symsize=plotsz[n_elements(magbin)-1-i], color=cgcolor('red')

      endfor

      read, 'Apply offsets? (y/n): ', qoff
      if (qoff eq 'y') then begin

	read, 'X offset: ', xoff
	read, 'Y offset: ', yoff

	offxy[0] = offxy[0]+xoff
	offxy[1] = offxy[1]+yoff

	xref = xref+xoff
	yref = yref+yoff

	window, 0, xs=npix, ys=npix
	oplot, xf, yf, psym=sym(11), symsize=2, color=cgcolor('green')
	for i=0,n_elements(magbin)-2 do begin

	  idx = where(hstmag ge magbin[i] and hstmag lt magbin[i+1])
	  if (idx[0] ne -1) then oplot, xref[idx], yref[idx], psym=sym(4), symsize=plotsz[n_elements(magbin)-1-i], color=cgcolor('red')

	endfor

      endif

    endrep until (qoff eq 'n')

    window, 0, xs=npix, ys=npix
    cgimage, im, /axis, minvalue=scale[0], maxvalue=scale[1], xtitle='X [px]', ytitle='Y [px]'

    oplot, xf, yf, psym=sym(11), symsize=2, color=cgcolor('green')
    for i=0,n_elements(magbin)-2 do begin

      idx = where(hstmag ge magbin[i] and hstmag lt magbin[i+1])
      if (idx[0] ne -1) then oplot, xref[idx], yref[idx], psym=sym(4), symsize=plotsz[n_elements(magbin)-1-i], color=cgcolor('red')

    endfor

    idx = where(xref le npix and xref ge 0.)
    xref = xref[idx]
    yref = yref[idx]
    hstmag = hstmag[idx]
    ra5 = ra5[idx]
    dec5 = dec5[idx]
    idx = where(yref le npix and yref ge 0.)
    xref = xref[idx]
    yref = yref[idx]
    hstmag = hstmag[idx]
    ra5 = ra5[idx]
    dec5 = dec5[idx]

  ;   fn = path+'temp_stars.txt'
  ;   openw, lun, fn, width=1400, /get_lun
  ; 
  ;     for i=0,n_elements(xref)-1 do printf, lun, xref[i], yref[i], hstmag[i], format='(2f12.6,f8.2)'
  ; 
  ;   close, lun
  ;   free_lun, lun


    ;find match between the two list (measured and catalog)

    idxsort = intarr(nstars)
    dmin = dblarr(nstars)
    for i=0,nstars-1 do begin

      tmp = sqrt((xf[i]-xref)^2.+(yf[i]-yref)^2.)

      dum = min(tmp, idxmin)
      idxsort[i] = idxmin
      dmin[i] = dum

    endfor

    xref = xref[idxsort]
    yref = yref[idxsort]
    ra5 = ra5[idxsort]
    dec5 = dec5[idxsort]
    hstmag = hstmag[idxsort]

    ;check if stars really match using distance
    idx = where(dmin gt median(dmin)+stddev(dmin))
    if (idx[0] ne -1) then begin

      xf[idx] = -99.
      yf[idx] = -99.
      xref[idx] = -99.
      yref[idx] = -99.
      ra5[idx] = -99.
      dec5[idx] = -99.
      raf[idx] = -99.
      decf[idx] = -99.
      hstmag[idx] = -99.

    endif

    idx0 = where(xf gt -1.)
    if (idx0[0] ne -1) then begin

      xf = xf[idx0]
      yf = yf[idx0]
      xref = xref[idx0]
      yref = yref[idx0]
      ra5 = ra5[idx0]
      dec5 = dec5[idx0]
      raf = raf[idx0]
      decf = decf[idx0]
      hstmag = hstmag[idx0]

    endif


;     window, 1, xs=npix, ys=npix, title='Selected Stars'
;     cgimage, im, /axis, minvalue=scale[0], maxvalue=scale[1], xtitle='X [px]', ytitle='Y [px]'
;     oplot, xf, yf, psym=sym(11), symsize=2, color=cgcolor('green')
;     for i=0,n_elements(magbin)-2 do begin
; 
;       idx = where(hstmag ge magbin[i] and hstmag lt magbin[i+1])
;       if (idx[0] ne -1) then oplot, xref[idx], yref[idx], psym=sym(4), symsize=plotsz[n_elements(magbin)-1-i], color=cgcolor('red')
; 
;     endfor
; 
;     print, ''
;     print, '***************************************'
;     print, 'Check plot if selection of stars is OK.'
;     print, '***************************************'
;     hak

    nuse = n_elements(xf)

    d5_deg = dblarr(nuse, nuse)
    df_px = dblarr(nuse, nuse)
    df_deg = dblarr(nuse, nuse)
    pa5 = dblarr(nuse, nuse)
    paf = dblarr(nuse, nuse)

    for i=0,nuse-1 do begin

      ;distance

      gcirc, 1, ra5[i], dec5[i], ra5, dec5, d
      d5_deg[i,*] = d	;arcsec

      gcirc, 1, raf[i], decf[i], raf, decf, d
      df_deg[i,*] = d

      df_px[i,*] = sqrt((xf[i]-xf)^2.+(yf[i]-yf)^2.)


      ;position angle
      ;computing it with star pairs - gives higher number than PA from reference coordinate, e.g. 5 stars available -> 10 PAs
      for j=0,nuse-1 do begin

	posang, 1, ra5[i], dec5[i], ra5[j], dec5[j], angle5
	if (angle5 lt 0.) then angle5 = angle5+360.
	pa5[i,j] = angle5

	posang, 1, raf[i], decf[i], raf[j], decf[j], anglef
	if (anglef lt 0.) then anglef = anglef+360.
	paf[i,j] = anglef

      endfor

    endfor

    ;convert from catalog dist and pixel into a platescale

    window, 1, xs=1000, ys=800, title='Computed Values'
    !p.multi=[0,1,2]
  ; 

    ;plate scale
    tmp = (d5_deg/df_px)[*]
    idxun = uniq(tmp, sort(tmp))
    tmp2 = tmp[idxun]
    tmp2 = tmp2[where(tmp2 gt 0.)]

    resistant_mean, tmp2*1.d3, 2., plscl, eplscl, nbad, /double, GoodVec=good, badvec=bad
    print, 'Number of bad plscl values        : ', nbad
    ;plscl = median(tmp3*1.d3)
    eplscl = stddev(tmp2[good]*1.d3)
    xaxis = indgen(n_elements(tmp2))
    plot, xaxis, tmp2*1.d3, /yn, psym=2, charsize=2, ytitle='PLSCL'
    oplot, xaxis[good], tmp2[good]*1.d3, color=cgcolor('red'), psym=2

    ;true north
    tmp = (paf-pa5)[*]
    tmp2 = tmp[idxun]
    tmp2 = tmp2[where(tmp2 ne 0.)]
    if (median(tmp2) lt 0.) then begin
      if (max(tmp2) gt 360.) then begin

	idxmax = where(tmp2 gt 360.)
	tmp2[idxmax] = tmp2[idxmax]-360.

      endif
    endif
    resistant_mean, tmp2, 2., tn, etn, nbad, /double, GoodVec=good
    print, 'Number of bad TN values           : ', nbad
    ;tn = median(tmp)
    etn = stddev(tmp2[good])
    xaxis = indgen(n_elements(tmp2))
    plot, xaxis, tmp2, /yn, psym=2, charsize=2, ytitle='TN'
    oplot, xaxis[good], tmp2[good], color=cgcolor('red'), psym=2

    !p.multi=[0,1,0]

    print, ''
    print, 'plate scale [mas]: ', sigfig(plscl, 5), ' +/- ', sigfig(eplscl, 2)
    print, 'true north [deg]: ', sigfig(tn, 3), ' +/- ', sigfig(etn, 2)
    print, ''

    if (eplscl gt 0.15) then begin

      print, '!!!!!!!!!!!!!!!!!!!!!'
      print, '!!! Check values! !!!'
      print, '!!!!!!!!!!!!!!!!!!!!!'

    endif

    mode = agpmflag
    fn = pathr+out+'_AstrometricCorrection.sav'
    save, plscl, eplscl, tn, etn, jd, mode, ptflag, field, nuse, offxy, filename=fn

  endfor

endif

end