@batch_NACO_Dark.pro
@batch_NACO_BPmap.pro
@batch_NACO_SkyFlat.pro
@batch_NACO_SkyBG_AstCal.pro
; @batch_NACO_SkyBG_v2.pro
@batch_NACO_reduce_Science_AstCal.pro
@batch_NACO_QC.pro
@circint_MOD.pro
@eq2hor_MOD.pro
@sign.pro
@readcol.pro
@remchar.pro
@gettok.pro
@headfits.pro
@fxposit.pro
@mrd_hread.pro
@get_eso_keyword.pro
@sigfig.pro
@mrdfits.pro
@fxpar.pro
@valid_num.pro
@mrd_skip.pro
@arrdelete.pro
@writefits.pro
@check_fits.pro
@fxaddpar.pro
@sxdelpar.pro
@sxpar.pro
@fits_add_checksum.pro
@checksum32.pro
@n_bytes.pro
@is_ieee_big.pro
@host_to_ieee.pro
@get_date.pro
@daycnv.pro
@sxaddpar.pro
@fits_ascii_encode.pro
@mpfitfun.pro
@mpfit.pro
@cgcolor.pro
@cggetcolorstate.pro
@cgsnapshot.pro
@cgcolor24.pro
@array_indices.pro
@robust_mean.pro
@avg.pro
@mkhdr.pro
@cgimage.pro
@image_dimensions.pro
@cgdefcharsize.pro
@setdefaultvalue.pro
@cgdefaultcolor.pro
@cgerase.pro
@cgsetcolorstate.pro
@clipscl.pro
@cgresizeimage.pro
@sixlin.pro
@proceeding_text.pro
@cgscalevector.pro
@fpufix.pro
@cghistoplot.pro                                                                                                                                                                     
@convert_to_type.pro                                                                                                                                                                 
@cgcheckforsymbols.pro     
@dist.pro 
@ten.pro
@closest.pro
@gaussscl.pro
@cgplot.pro
@cgbitget.pro
@colorsareidentical.pro
@fixpix_mod.pro
@dist_circle.pro
@caldat.pro
@precess.pro
@premat.pro
@co_nutate.pro
@nutate.pro
@poly.pro
@cirrange.pro
@isarray.pro
@co_aberration.pro
@sunpos.pro
@addpm.pro
@ct2lst.pro
@hadec2altaz.pro
@co_refract.pro
@parangle.pro
@mpfit2dpeak.pro
@mpfit2dfun.pro
@fftshift.pro
@interpol.pro
@resistant_mean.pro
@scale_image_am.pro
@sky.pro
@mmm.pro
@asinh.pro
@remove.pro
@ts_diff.pro
@showsym.pro
@batch_NACO_find_stars_AstCal.pro
@batch_NACO_plscl_tn.pro
@strnumber.pro
@find.pro

pro batch_NACO_reduce_AstCal, display=display

;===============================================================================================

readcol, '/home/amueller/work/IDLlibs/AO/NACO/datapaths.txt', tmp, format='a', /silent
print, ''
for i=0,n_elements(tmp)-1 do print, strcompress(i+1, /rem), ' ', tmp[i]

print, ''
read, 'Select Path: ', selp
; selp = 1
path = tmp[selp-1]+'RAW/'

pathr = path+'../Reduced/'
spawn, 'mkdir '+path+'../Reduced/'

if (strmatch(path, '*47Tuc*') eq 1) then field = '47Tuc'
if (strmatch(path, '*Trapezium*') eq 1) then field = 'Trapezium'


;===============================================================================================

if (strmatch(path, '*_PT*') eq 1) then ptflag = 1 else ptflag=0
agpmflag = strmatch(path, '*noAGPM*')
agpmflag = abs(agpmflag-1)	;no agpm = 0, w/ agpm = 1

;===============================================================================================

plsc = 27.19d-3
diam = 8.2
lambda = 3.8d-6	;L'
fwhm = lambda/diam*206265.d0/plsc

;===============================================================================================

;DARK
batch_NACO_Dark, path

;===============================================================================================

;Bad Pixel map
batch_NACO_BPmap, path, display

;===============================================================================================

;Sky Flat
batch_NACO_SkyFlat, path, display

;===============================================================================================

;Sky extraction
batch_NACO_SkyBG_AstCal, path

;===============================================================================================

;reduce cubes
batch_NACO_reduce_Science_AstCal, path, field, ptflag

;===============================================================================================

;find stars
batch_NACO_find_stars_AstCal, pathr, fwhm

;===============================================================================================

;compute TN and place scale
batch_NACO_plscl_tn, pathr, field, ptflag, agpmflag

;===============================================================================================

print, ''
print, 'Reduction finished!'
print, ''
stop
end
