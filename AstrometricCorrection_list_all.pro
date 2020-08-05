@readcol.pro
@remchar.pro
@gettok.pro
@strnumber.pro
@cgcolor.pro
@cggetcolorstate.pro
@oploterror.pro
@setdefaultvalue.pro
@cgquery.pro
@cgplot.pro
@cgsetcolorstate.pro
@cgcheckforsymbols.pro
@cgdefaultcolor.pro
@cgcolor24.pro
@colorsareidentical.pro
@cgdefcharsize.pro
@cgsymcat.pro
@sigfig.pro
@legend.pro
@showsym.pro

pro AstrometricCorrection_list_all

base = '/data/beegfs/astro-storage/groups/launhardt/amueller/NACO/data/'
resdir = base+'AstrometricCorrection/'
spawn, 'mkdir -p '+resdir

; file = file_search(base+'*/*/AstCal*/AstrometricCorrection.sav', count=n)
file = file_search(base+'*/*/AstCal*/Reduced/*AstrometricCorrection.sav', count=n)

jd_all = dblarr(n)
plscl_all = dblarr(n)
eplscl_all = dblarr(n)
tn_all = dblarr(n)
etn_all = dblarr(n)
mode_all = strarr(n)
field_all = strarr(n)
nuse_all = intarr(n)
ptflag_all = intarr(n)

for i=0,n-1 do begin

  restore, file[i]
  jd_all[i] = jd
  plscl_all[i] = plscl
  eplscl_all[i] = eplscl
  tn_all[i] = -1.*tn
  etn_all[i] = etn
  mode_all[i] = mode
  field_all[i] = field
  nuse_all[i] = nuse
  ptflag_all[i] = ptflag

endfor

idx = sort(jd_all)
jd_all = jd_all[idx]
field_all = field_all[idx]
mode_all = mode_all[idx]
plscl_all = plscl_all[idx]
eplscl_all = eplscl_all[idx]
tn_all = tn_all[idx]
etn_all = etn_all[idx]
nuse_all = nuse_all[idx]
ptflag_all = ptflag_all[idx]

ptflag_allo = ptflag_all
ptflag_all = strarr(n_elements(ptflag_allo))
idx = where(ptflag_allo eq 1)
ptflag_all[idx] = 'Pupil'
idx = where(ptflag_allo eq 0)
ptflag_all[idx] = 'Field'

mode_allo = mode_all
mode_all = strarr(n_elements(mode_allo))
idx = where(mode_allo eq 1)
mode_all[idx] = 'AGPM'
idx = where(mode_allo eq 0)
mode_all[idx] = 'satPSF'


fn = resdir+'AstrometricCorrection.txt'
openw, lun, fn, width=2000, /get_lun

  printf, lun, '       Field                JD    Mode   Plate Scale [mas]  Error   True North [deg]  Error  #Stars Tracking'
  for i=0,n-1 do printf, lun, field_all[i], jd_all[i], mode_all[i], plscl_all[i], eplscl_all[i], tn_all[i], etn_all[i], nuse_all[i], ptflag_all[i], format='(a12,f18.8,a8,f8.3,f19.3,f9.3,f17.3,i5,a7)'


close, lun
free_lun, lun

;resdir = '/home/amueller/Downloads/'

readcol, resdir+'AstrometricCorrection.txt', field, jd, mode, plscl, eplscl, tn, etn, nuse, ptflag, format='(a,d,a,d,d,d,d,i,a)', /silent

jd = jd-2450000.d0

ofield = field
ojd = jd
omode = mode
oplscl = plscl
oeplscl = eplscl
otn = tn
oetn = etn
onuse = nuse
optflag = ptflag

idx = where(finite(plscl) eq 1)
field = field[idx]
jd = jd[idx]
mode = mode[idx]
plscl = plscl[idx]
eplscl = eplscl[idx]
tn = tn[idx]
etn = etn[idx]
nuse = nuse[idx]
ptflag = ptflag[idx]

idx = where(finite(eplscl) eq 1)
field = field[idx]
jd = jd[idx]
mode = mode[idx]
plscl = plscl[idx]
eplscl = eplscl[idx]
tn = tn[idx]
etn = etn[idx]
nuse = nuse[idx]
ptflag = ptflag[idx]

idx = where(finite(tn) eq 1)
field = field[idx]
jd = jd[idx]
mode = mode[idx]
plscl = plscl[idx]
eplscl = eplscl[idx]
tn = tn[idx]
etn = etn[idx]
nuse = nuse[idx]
ptflag = ptflag[idx]

idx = where(finite(etn) eq 1)
field = field[idx]
jd = jd[idx]
mode = mode[idx]
plscl = plscl[idx]
eplscl = eplscl[idx]
tn = tn[idx]
etn = etn[idx]
nuse = nuse[idx]
ptflag = ptflag[idx]

resistant_mean, plscl, 4, t1, t2, nbad, /double, goodvec=idx, badvec=idxb
field = field[idx]
jd = jd[idx]
mode = mode[idx]
plscl = plscl[idx]
eplscl = eplscl[idx]
tn = tn[idx]
etn = etn[idx]
nuse = nuse[idx]
ptflag = ptflag[idx]

resistant_mean, eplscl, 4, t1, t2, nbad, /double, goodvec=idx, badvec=idxb
field = field[idx]
jd = jd[idx]
mode = mode[idx]
plscl = plscl[idx]
eplscl = eplscl[idx]
tn = tn[idx]
etn = etn[idx]
nuse = nuse[idx]
ptflag = ptflag[idx]

resistant_mean, tn, 4, t1, t2, nbad, /double, goodvec=idx, badvec=idxb
field = field[idx]
jd = jd[idx]
mode = mode[idx]
plscl = plscl[idx]
eplscl = eplscl[idx]
tn = tn[idx]
etn = etn[idx]
nuse = nuse[idx]
ptflag = ptflag[idx]

resistant_mean, etn, 4, t1, t2, nbad, /double, goodvec=idx, badvec=idxb
field = field[idx]
jd = jd[idx]
mode = mode[idx]
plscl = plscl[idx]
eplscl = eplscl[idx]
tn = tn[idx]
etn = etn[idx]
nuse = nuse[idx]
ptflag = ptflag[idx]

;=========================================================================

;idxgood provides the selected data sets used for the analysis

flag = intarr(n_elements(ojd))
for i=0,n_elements(ojd)-1 do begin

  idx = where(ojd[i] eq jd)
  if (idx[0] ne -1) then flag[i] = 1

endfor
idxgood = where(flag eq 1)
print, ''
print, 'Selected Measurements (IDL index):'
print, idxgood

print, ''
print, 'Bad Measurements (IDL index):'
print, where(flag eq 0)


;=========================================================================

idxs = where(mode eq 'satPSF')
idxa = where(mode eq 'AGPM')

; avens = median(nuse[idxs],/even)
; avena = median(nuse[idxa],/even)
avens = mean(nuse[idxs])
avena = mean(nuse[idxa])
ratio = avens/avena

print, ''
print, 'Average number of stars satPSF: ', round(avens)
print, 'Average number of stars w/ AGPM: ', round(avena)
print, 'Ratio: ', ratio

;=========================================================================

;weighted mean AGPM 
idxa = where(mode eq 'AGPM')
awmplscl = (total((1./eplscl[idxa]^2.)*plscl[idxa]))/total(1./eplscl[idxa]^2.)
awmplsclerr = sqrt(n_elements(plscl[idxa])/total(1./eplscl[idxa]^2.))
awmtn = (total((1./etn[idxa]^2.)*tn[idxa]))/total(1./etn[idxa]^2.)
awmtnerr = sqrt(n_elements(tn[idxa])/total(1./etn[idxa]^2.))

;weighted mean AGPM - field tracking
idxa_ft = where(mode eq 'AGPM' and ptflag eq 'Field')
awmplscl_ft = (total((1./eplscl[idxa_ft]^2.)*plscl[idxa_ft]))/total(1./eplscl[idxa_ft]^2.)
awmplsclerr_ft = sqrt(n_elements(plscl[idxa_ft])/total(1./eplscl[idxa_ft]^2.))
awmtn_ft = (total((1./etn[idxa_ft]^2.)*tn[idxa_ft]))/total(1./etn[idxa_ft]^2.)
awmtnerr_ft = sqrt(n_elements(tn[idxa_ft])/total(1./etn[idxa_ft]^2.))

;weighted mean AGPM - pupil tracking
idxa_pt = where(mode eq 'AGPM' and ptflag eq 'Pupil')
awmplscl_pt = (total((1./eplscl[idxa_pt]^2.)*plscl[idxa_pt]))/total(1./eplscl[idxa_pt]^2.)
awmplsclerr_pt = sqrt(n_elements(plscl[idxa_pt])/total(1./eplscl[idxa_pt]^2.))
awmtn_pt = (total((1./etn[idxa_pt]^2.)*tn[idxa_pt]))/total(1./etn[idxa_pt]^2.)
awmtnerr_pt = sqrt(n_elements(tn[idxa_pt])/total(1./etn[idxa_pt]^2.))

;weighted mean satPSF
idxs = where(mode eq 'satPSF')
swmplscl = (total((1./eplscl[idxs]^2.)*plscl[idxs]))/total(1./eplscl[idxs]^2.)
swmplsclerr = sqrt(n_elements(plscl[idxs])/total(1./eplscl[idxs]^2.))
swmtn = (total((1./etn[idxs]^2.)*tn[idxs]))/total(1./etn[idxs]^2.)
swmtnerr = sqrt(n_elements(tn[idxs])/total(1./etn[idxs]^2.))

;weighted mean satPSF - field tracking
idxs_ft = where(mode eq 'satPSF' and ptflag eq 'Field')
swmplscl_ft = (total((1./eplscl[idxs_ft]^2.)*plscl[idxs_ft]))/total(1./eplscl[idxs_ft]^2.)
swmplsclerr_ft = sqrt(n_elements(plscl[idxs_ft])/total(1./eplscl[idxs_ft]^2.))
swmtn_ft = (total((1./etn[idxs_ft]^2.)*tn[idxs_ft]))/total(1./etn[idxs_ft]^2.)
swmtnerr_ft = sqrt(n_elements(tn[idxs_ft])/total(1./etn[idxs_ft]^2.))

;weighted mean satPSF - pupil tracking
idxs_pt = where(mode eq 'satPSF' and ptflag eq 'Pupil')
swmplscl_pt = (total((1./eplscl[idxs_pt]^2.)*plscl[idxs_pt]))/total(1./eplscl[idxs_pt]^2.)
swmplsclerr_pt = sqrt(n_elements(plscl[idxs_pt])/total(1./eplscl[idxs_pt]^2.))
swmtn_pt = (total((1./etn[idxs_pt]^2.)*tn[idxs_pt]))/total(1./etn[idxs_pt]^2.)
swmtnerr_pt = sqrt(n_elements(tn[idxs_pt])/total(1./etn[idxs_pt]^2.))


;output
print, ''
print, 'Weighted means - UPDATE Astrometry_Photometry'
print, 'AGPM'
print, 'Plate scale', awmplscl, awmplsclerr
print, 'True North', awmtn, awmtnerr
print, 'AGPM - field tracking'
print, 'Plate scale', awmplscl_ft, awmplsclerr_ft
print, 'True North', awmtn_ft, awmtnerr_ft
print, 'AGPM - pupil tracking'
print, 'Plate scale', awmplscl_pt, awmplsclerr_pt
print, 'True North', awmtn_pt, awmtnerr_pt
print, ''
print, 'satPSF'
print, 'Plate scale', swmplscl, swmplsclerr
print, 'True North', swmtn, swmtnerr
print, 'satPSF - field tracking'
print, 'Plate scale', swmplscl_ft, swmplsclerr_ft
print, 'True North', swmtn_ft, swmtnerr_ft
print, 'satPSF - pupil tracking'
print, 'Plate scale', swmplscl_pt, swmplsclerr_pt
print, 'True North', swmtn_pt, swmtnerr_pt

;symsz = [2.0,1.5,1.0,0.7,0.5,0.3]
xr = [7300, max(jd)+50.]
yrtn = [27.0,27.4]

set_plot, 'ps'
fn = resdir+'AstrometricCorrectionMeasurements_Selected'
device, filename=fn+'.ps',/color, XSIZE=30, YSIZE=30
  !p.font=0
  !p.thick=3
  !x.thick=3
  !y.thick=3
  !p.multi=[0,1,2]

  plot, jd, plscl, /nodata, xtitle='JD-2450000 [days]', ytitle='Plate scale [mas]', /yn, yr=[yrtn[0],yrtn[1]], yst=1, charsize=2, xr=[xr[0],xr[1]], xst=1, title='Cleaned Data Set - PT and FT'
  oplot, jd[idxa], plscl[idxa], color=cgcolor('green')
  oploterror, jd[idxa], plscl[idxa], eplscl[idxa], /nohat, color=cgcolor('green'), errcolor=cgcolor('green'), psym=sym(1), symsize=1.2
  oplot, jd[idxs], plscl[idxs], color=cgcolor('blue')
  oploterror, jd[idxs], plscl[idxs], eplscl[idxs], /nohat, color=cgcolor('blue'), errcolor=cgcolor('blue'), psym=sym(1), symsize=1.2
  legend, ['AGPM '+sigfig(awmplscl,5)+'+/-'+sigfig(awmplsclerr,3)+' mas', 'satPSF '+sigfig(swmplscl,5)+'+/-'+sigfig(swmplsclerr,3)+' mas'], box=0, margin=0, /top, /left, textcolor=[cgcolor('green'),cgcolor('blue')], charsize=1.5

  plot, jd, tn, /nodata, xtitle='JD-2450000 [days]', ytitle='True North [deg]', /yn, yr=[0.,1.], yst=1, charsize=2, xr=[xr[0],xr[1]], xst=1
  oplot, jd[idxa], tn[idxa], color=cgcolor('green')
  oploterror, jd[idxa], tn[idxa], etn[idxa], /nohat, color=cgcolor('green'), errcolor=cgcolor('green'), psym=sym(1), symsize=1.2
  oplot, jd[idxs], tn[idxs], color=cgcolor('blue')
  oploterror, jd[idxs], tn[idxs], eplscl[idxs], /nohat, color=cgcolor('blue'), errcolor=cgcolor('blue'), psym=sym(1), symsize=1.2
  legend, ['AGPM '+sigfig(awmtn,3)+'+/-'+sigfig(awmtnerr,3)+' deg', 'satPSF '+sigfig(swmtn,3)+'+/-'+sigfig(swmtnerr,3)+' deg'], box=0, margin=0, /top, /left, textcolor=[cgcolor('green'),cgcolor('blue')], charsize=1.5
  !p.multi=[0,1,0]
  !p.thick=1
  !x.thick=1
  !y.thick=1

device,/close
set_plot,'x'
spawn, '/home/amueller/work/IDLlibs/epstopdf.pl '+fn+'.ps'
spawn, 'rm '+fn+'.ps'

;=========================================================================

set_plot, 'ps'
fn = resdir+'AstrometricCorrectionMeasurements_FieldTrack'
device, filename=fn+'.ps',/color, XSIZE=30, YSIZE=30
  !p.font=0
  !p.thick=3
  !x.thick=3
  !y.thick=3
  !p.multi=[0,1,2]

  plot, jd, plscl, /nodata, xtitle='JD-2450000 [days]', ytitle='Plate scale [mas]', /yn, yr=[yrtn[0],yrtn[1]], yst=1, charsize=2, xr=[xr[0],xr[1]], xst=1, title='Field Track'
  oplot, jd[idxa_ft], plscl[idxa_ft], color=cgcolor('green')
  oploterror, jd[idxa_ft], plscl[idxa_ft], eplscl[idxa_ft], /nohat, color=cgcolor('green'), errcolor=cgcolor('green'), psym=sym(1), symsize=1.2
  oplot, jd[idxs_ft], plscl[idxs_ft], color=cgcolor('blue')
  oploterror, jd[idxs_ft], plscl[idxs_ft], eplscl[idxs_ft], /nohat, color=cgcolor('blue'), errcolor=cgcolor('blue'), psym=sym(1), symsize=1.2
  legend, ['AGPM '+sigfig(awmplscl_ft,5)+'+/-'+sigfig(awmplsclerr_ft,3)+' mas', 'satPSF '+sigfig(swmplscl_ft,5)+'+/-'+sigfig(swmplsclerr_ft,3)+' mas'], box=0, margin=0, /top, /left, textcolor=[cgcolor('green'),cgcolor('blue')], charsize=1.5

  plot, jd, tn, /nodata, xtitle='JD-2450000 [days]', ytitle='True North [deg]', /yn, yr=[0.,1.], yst=1, charsize=2, xr=[xr[0],xr[1]], xst=1
  oplot, jd[idxa_ft], tn[idxa_ft], color=cgcolor('green')
  oploterror, jd[idxa_ft], tn[idxa_ft], etn[idxa_ft], /nohat, color=cgcolor('green'), errcolor=cgcolor('green'), psym=sym(1), symsize=1.2
  oplot, jd[idxs_ft], tn[idxs_ft], color=cgcolor('blue')
  oploterror, jd[idxs_ft], tn[idxs_ft], eplscl[idxs_ft], /nohat, color=cgcolor('blue'), errcolor=cgcolor('blue'), psym=sym(1), symsize=1.2
  legend, ['AGPM '+sigfig(awmtn_ft,3)+'+/-'+sigfig(awmtnerr_ft,3)+' deg', 'satPSF '+sigfig(swmtn_ft,3)+'+/-'+sigfig(swmtnerr_ft,3)+' deg'], box=0, margin=0, /top, /left, textcolor=[cgcolor('green'),cgcolor('blue')], charsize=1.5
  !p.multi=[0,1,0]
  !p.thick=1
  !x.thick=1
  !y.thick=1

device,/close
set_plot,'x'
spawn, '/home/amueller/work/IDLlibs/epstopdf.pl '+fn+'.ps'
spawn, 'rm '+fn+'.ps'

;=========================================================================

set_plot, 'ps'
fn = resdir+'AstrometricCorrectionMeasurements_PupilTrack'
device, filename=fn+'.ps',/color, XSIZE=30, YSIZE=30
  !p.font=0
  !p.thick=3
  !x.thick=3
  !y.thick=3
  !p.multi=[0,1,2]

  plot, jd, plscl, /nodata, xtitle='JD-2450000 [days]', ytitle='Plate scale [mas]', /yn, yr=[yrtn[0],yrtn[1]], yst=1, charsize=2, xr=[xr[0],xr[1]], xst=1, title='Pupil Track'
  oplot, jd[idxa_pt], plscl[idxa_pt], color=cgcolor('green')
  oploterror, jd[idxa_pt], plscl[idxa_pt], eplscl[idxa_pt], /nohat, color=cgcolor('green'), errcolor=cgcolor('green'), psym=sym(1), symsize=1.2
  oplot, jd[idxs_pt], plscl[idxs_pt], color=cgcolor('blue')
  oploterror, jd[idxs_pt], plscl[idxs_pt], eplscl[idxs_pt], /nohat, color=cgcolor('blue'), errcolor=cgcolor('blue'), psym=sym(1), symsize=1.2
  legend, ['AGPM '+sigfig(awmplscl_pt,5)+'+/-'+sigfig(awmplsclerr_pt,3)+' mas', 'satPSF '+sigfig(swmplscl_pt,5)+'+/-'+sigfig(swmplsclerr_pt,3)+' mas'], box=0, margin=0, /top, /left, textcolor=[cgcolor('green'),cgcolor('blue')], charsize=1.5

  plot, jd, tn, /nodata, xtitle='JD-2450000 [days]', ytitle='True North [deg]', /yn, yr=[0.,1.], yst=1, charsize=2, xr=[xr[0],xr[1]], xst=1
  oplot, jd[idxa_pt], tn[idxa_pt], color=cgcolor('green')
  oploterror, jd[idxa_pt], tn[idxa_pt], etn[idxa_pt], /nohat, color=cgcolor('green'), errcolor=cgcolor('green'), psym=sym(1), symsize=1.2
  oplot, jd[idxs_pt], tn[idxs_pt], color=cgcolor('blue')
  oploterror, jd[idxs_pt], tn[idxs_pt], eplscl[idxs_pt], /nohat, color=cgcolor('blue'), errcolor=cgcolor('blue'), psym=sym(1), symsize=1.2
  legend, ['AGPM '+sigfig(awmtn_pt,3)+'+/-'+sigfig(awmtnerr_pt,3)+' deg', 'satPSF '+sigfig(swmtn_pt,3)+'+/-'+sigfig(swmtnerr_pt,3)+' deg'], box=0, margin=0, /top, /left, textcolor=[cgcolor('green'),cgcolor('blue')], charsize=1.5
  !p.multi=[0,1,0]
  !p.thick=1
  !x.thick=1
  !y.thick=1

device,/close
set_plot,'x'
spawn, '/home/amueller/work/IDLlibs/epstopdf.pl '+fn+'.ps'
spawn, 'rm '+fn+'.ps'

;=========================================================================

set_plot, 'ps'
fn = resdir+'AstrometricCorrectionMeasurements_All'
device, filename=fn+'.ps',/color, XSIZE=30, YSIZE=30
  !p.font=0
  !p.thick=3
  !x.thick=3
  !y.thick=3
  !p.multi=[0,1,2]

  idxa = where(omode eq 'AGPM')
  idxs = where(omode eq 'satPSF')

  plot, ojd, oplscl, /nodata, xtitle='JD-2450000 [days]', ytitle='Plate scale [mas]', /yn, yr=[yrtn[0],yrtn[1]], yst=1, charsize=2, xr=[xr[0],xr[1]], xst=1, title='All Data'
  oplot, ojd[idxa], oplscl[idxa], color=cgcolor('green')
  oploterror, ojd[idxa], oplscl[idxa], oeplscl[idxa], /nohat, color=cgcolor('green'), errcolor=cgcolor('green'), psym=sym(1), symsize=1.2
  oplot, ojd[idxs], oplscl[idxs], color=cgcolor('blue')
  oploterror, ojd[idxs], oplscl[idxs], oeplscl[idxs], /nohat, color=cgcolor('blue'), errcolor=cgcolor('blue'), psym=sym(1), symsize=1.2
  legend, ['AGPM', 'satPSF'], box=0, margin=0, /top, /left, textcolor=[cgcolor('green'),cgcolor('blue')], charsize=1.5

  plot, ojd, otn, /nodata, xtitle='JD-2450000 [days]', ytitle='True North [deg]', /yn, yr=[-0.5,1.5], yst=1, charsize=2, xr=[xr[0],xr[1]], xst=1
  oplot, ojd[idxa], otn[idxa], color=cgcolor('green')
  oploterror, ojd[idxa], otn[idxa], oetn[idxa], /nohat, color=cgcolor('green'), errcolor=cgcolor('green'), psym=sym(1), symsize=1.2
  oplot, ojd[idxs], otn[idxs], color=cgcolor('blue')
  oploterror, ojd[idxs], otn[idxs], oeplscl[idxs], /nohat, color=cgcolor('blue'), errcolor=cgcolor('blue'), psym=sym(1), symsize=1.2
  legend, ['AGPM', 'satPSF'], box=0, margin=0, /top, /left, textcolor=[cgcolor('green'),cgcolor('blue')], charsize=1.5
  !p.multi=[0,1,0]
  !p.thick=1
  !x.thick=1
  !y.thick=1

device,/close
set_plot,'x'
spawn, '/home/amueller/work/IDLlibs/epstopdf.pl '+fn+'.ps'
spawn, 'rm '+fn+'.ps'

;=========================================================================



stop
end
