pro batch_NACO_SkyBG_AstCal, path;, agpmflag

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

idx = where(catg eq 'SCIENCE' and type eq 'SKY')
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
  print, 'More than 1 set of different sky files. There should be only 1 set of sky and 1 set of science! Stop.'
  stop

endif

hdr = headfits(file[0], exten=0, /silent)
nx = get_eso_keyword(hdr, 'NAXIS1')
ny = get_eso_keyword(hdr, 'NAXIS2')

sky = dblarr(nx,nx,n_elements(file))
for i=0,n_elements(file)-1 do begin

  im = mrdfits(file[i],0,/silent)
  im = im[*,0:nx-1,30:n_elements(im[0,0,*])-2]
  sky[*,*,i] = median(im,dim=3,/even)

endfor

skyname = path+'../Reduced/'+'sky_background_'+ust+'.fits'
if (n_elements(sky[0,0,*]) gt 1) then writefits, skyname, median(sky,dim=3,/even), hdr else writefits, skyname, sky, hdr



end
