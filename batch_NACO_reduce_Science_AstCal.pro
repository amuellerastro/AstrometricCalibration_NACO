function get_parangle, jd, epoch, ra, dec, pma, pmd, radvel, plx, lat, lon, altitude, temp, pres

  caldat, jd, mo, day, yr, hh, mm, ss

  ;get date of the form yyyy-mm-ddThh:mm:ss
  date = dblarr(6)
  date[0] = yr
  date[1] = mo
  date[2] = day
  date[3] = hh
  date[4] = mm
  date[5] = ss

  ;get year of observations for new RA, DEC, i.e. epoch
  hour = hh+mm/60.d0+ss/3600.d0
  if ((double(yr) mod 4.) eq 0.) then date2 = yr+mo/12.d0+day/366.d0+hour/8766.0d0 $
    else date2 = yr+mo/12.d0+day/365.d0+hour/8766.0d0	;takes leap year into account

  ;compute current cordinates and parallactic angle at time of observation
  dt = date2-epoch
  eq2hor_MOD, ra, dec, pma, pmd, radvel, plx, dt, jd, alt, az, hatmp, lat=lat, lon=lon, altitude=altitude, pres=pres, temp=temp, outra=outra, outdec=outdec;, /verbose
  ha = hatmp
  ha = ha/15.

  parang = parangle(ha, outdec, lat)

  if (outdec gt lat and parang lt 0.) then parang = parang+360.d0

  return, {parang:parang, ha:ha}

end

pro batch_NACO_reduce_Science_AstCal, path, field, ptflag

if (ptflag eq 1) then begin

  filestar = file_search('/home/amueller/work/IDLlibs/AO/TargetProperties/Targets/'+field+'.sav', count=nstars)
  restore, filestar, /v

  ra_st = ra
  dec_st = dec

endif

fileo = file_search(path+'NACO*fits', count=nfiles)

; out = strarr(nfiles)
object = strarr(nfiles)
type = strarr(nfiles)
catg = strarr(nfiles)
opti1 = strarr(nfiles)
opti3 = strarr(nfiles)
opti6 = strarr(nfiles)
opti7 = strarr(nfiles)
dit = strarr(nfiles)
naxis1 = strarr(nfiles)
naxis2 = strarr(nfiles)
ndit = strarr(nfiles)

st = strarr(nfiles)

for i=0,nfiles-1 do begin

  pos1 = strpos(fileo[i], '/', /reverse_search)
  pos2 = strpos(fileo[i], '.fits', /reverse_search)
;   out[i] = strmid(fileo[i], pos1+1, pos2-pos1-1)

  hdr = headfits(fileo[i], exten=0, /silent)

  object[i] = strcompress(get_eso_keyword(hdr,'OBJECT'),/rem)
  naxis1[i] = strcompress(get_eso_keyword(hdr,'NAXIS1'),/rem)
  naxis2[i] = strcompress(get_eso_keyword(hdr,'NAXIS2'),/rem)
  type[i] = strcompress(get_eso_keyword(hdr,'HIERARCH ESO DPR TYPE'),/rem)
  catg[i] = strcompress(get_eso_keyword(hdr,'HIERARCH ESO DPR CATG'),/rem)
  dit[i] = strcompress(get_eso_keyword(hdr,'HIERARCH ESO DET DIT'),/rem)
  dit[i] = sigfig(dit[i],2)
  opti1[i] = strcompress(get_eso_keyword(hdr,'HIERARCH ESO INS OPTI1 NAME'),/rem)	;mask
  opti3[i] = strcompress(get_eso_keyword(hdr,'HIERARCH ESO INS OPTI3 NAME'),/rem)	;pupil stop
  opti6[i] = strcompress(get_eso_keyword(hdr,'HIERARCH ESO INS OPTI6 NAME'),/rem)	;filter
  opti7[i] = strcompress(get_eso_keyword(hdr,'HIERARCH ESO INS OPTI7 NAME'),/rem)	;Objective
  ndit[i] = strcompress(get_eso_keyword(hdr, 'HIERARCH ESO DET NDIT'),/rem)

  st[i] = dit[i]+'s_'+naxis1[i]+'px'

endfor

;===============================================================================================

;NO AGPM used, sky from science frames itself

; if (agpmflag eq 0) then begin

idx = where(catg eq 'SCIENCE' and type eq 'OBJECT')
file = fileo[idx]
object = object[idx]
type = type[idx]
catg = catg[idx]
opti1 = opti1[idx]
opti3 = opti3[idx]
opti6 = opti6[idx]
opti7 = opti7[idx]
dit = dit[idx]
naxis1 = naxis1[idx]
naxis2 = naxis2[idx]
st = st[idx]

idx = uniq(st, sort(st))
ust = st[idx]
nu = n_elements(idx)
if (nu gt 1) then begin

  print, ''
  print, 'More than 1 set of different object files. There should be only 1 set of sky and 1 set of science! Stop.'
  stop

endif

;bad pixel
bpf = path+'../Reduced/'+'static_badpixels_'+ust+'.fits'
;dark
dkf = path+'../Reduced/'+'master_dark_'+ust+'.fits'
;flat field
fff = path+'../Reduced/'+'instrument_flat.fits'
;sky
skf = path+'../Reduced/'+'sky_background_'+ust+'.fits'

ff = mrdfits(fff, 0, hdrff, /silent)
bp = mrdfits(bpf, 0, hdrbg, /silent)
dk = mrdfits(dkf, 0, hdrdk, /silent)
sk = mrdfits(skf, 0, hdrsk, /silent)

sx = uint(get_eso_keyword(hdrdk, 'HIERARCH ESO DET WIN STARTX'))
sy = uint(get_eso_keyword(hdrdk, 'HIERARCH ESO DET WIN STARTY'))

nx = uint(naxis1[0]) & ny = nx

ff = ff[sx-1:sx+nx-2,sy-1:sy+ny-2,*]

for i=0,n_elements(file)-1 do begin

  pos1 = strpos(file[i], '/', /reverse_search)
  pos2 = strpos(file[i], '.fits', /reverse_search)
  outname = strmid(file[i], pos1+1, pos2-pos1-1)

  cube = mrdfits(file[i], 0, hdrraw, /silent)
  if (n_elements(cube[0,0,*]) gt 1) then cube = cube[*,0:nx-1,30:n_elements(cube[0,0,*])-2]
  nframes = n_elements(cube[0,0,*])

  exptime = double(get_eso_keyword(hdrraw, 'EXPTIME'))
  ndit = double(get_eso_keyword(hdrraw, 'HIERARCH ESO DET NDIT'))
  dit = double(get_eso_keyword(hdrraw,'HIERARCH ESO DET DIT'))
  jd = double(get_eso_keyword(hdrraw, 'MJD-OBS'))+2400000.5d0
  pres = double(get_eso_keyword(hdrraw, 'HIERARCH ESO TEL AMBI PRES START'))
  temp = double(get_eso_keyword(hdrraw, 'HIERARCH ESO TEL AMBI TEMP'))+274.15d0
  epoch = double(get_eso_keyword(hdrraw, 'HIERARCH ESO TEL TARG EPOCH'))
  pupilpos = double(get_eso_keyword(hdrraw, 'HIERARCH ESO ADA PUPILPOS'))

  tel = get_eso_keyword(hdrraw, 'TELESCOP')
  ;Observatory parameters
  ;from http://www.eso.org/sci/facilities/paranal/astroclimate/site.html

  if (strmatch(tel, '*U4*') eq 1) then begin

    lat = ten(-24.d0, 37.d0, 31.00d0)	;for UT4
    lon = ten(-70.d0, 24.d0, 8.00d0)	;for UT4

  endif

  if (strmatch(tel, '*U1*') eq 1) then begin

    lat = ten(-24.d0, 37.d0, 33.117d0)	;for UT1
    lon = ten(-70.d0, 24.d0, 11.642d0)	;for UT1

  endif

  altitude = 2635.43d0

  corim = dblarr(nx,nx,nframes)

  for j=0,nframes-1 do begin

    corim[*,*,j] = (cube[*,*,j]-sk)/ff
    tmp = corim[*,*,j]
    fixpix_mod, tmp, bp, outim, npix=24, /weight, /silent
    corim[*,*,j] = outim

    proceeding_text,loop=nframes, i=j, prompt='> Reducing Frame   '+string(j+1,form='(I4)')

  endfor

  if (ptflag eq 1) then begin

    ;compute parallactic angle for exposure start and end of a cube

    ra = [double(strmid(ra_st,0,2)),double(strmid(ra_st,2,2)),double(strmid(ra_st,4,6))]
    sign1 = strmid(dec_st,0,1)
    if (sign1 eq '-') then sign1 = -1.d0 else sign1 = 1.d0
    dec = [sign1*double(strmid(dec_st,1,2)), double(strmid(dec_st,3,2)), double(strmid(dec_st,5,6))]
    ra = ten(ra)*15.
    dec = ten(dec)

    if (finite(radvel[0] ne 1)) then radvel = 0.
    if (finite(plx[0] ne 1)) then plx = 1./100.

    parang = dblarr(nframes)

    for j=0,nframes-1 do begin

      ;+1 because I remove the first frame anyway
      st = get_parangle(jd+(dit*(j+1.))/86400.d0, epoch, ra, dec, pma, pmd, radvel, plx, lat, lon, altitude, temp, pres)
      parang[j] = st.parang

    endfor


    rotptoff = 90.d0+(89.44-pupilpos)	;user manual page 71

    PA_onsky = parang-rotptoff

    dim = n_elements(corim[*,0,0])
    rotcorim = corim
    for j=0,nframes-1 do rotcorim[*,*,j]=rot(corim[*,*,j],-PA_onsky[j],1.0,dim/2,dim/2,cubic=-0.5,/pivot)

    medim = median(rotcorim, dim=3, /even)

  endif else begin

    if (n_elements(corim[0,0,*]) gt 1) then medim = median(corim, dim=3, /even) else medim = corim

  endelse

  writefits, path+'../Reduced/'+'AstrometricField_'+outname+'.fits', medim, hdrraw

endfor


end
