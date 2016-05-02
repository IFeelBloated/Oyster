import vapoursynth as vs
import mvmulti

### Global Settings ###
fmtc_args              = dict (fulls=True, fulld=True)
msuper_args            = dict (hpad=32, vpad=32, sharp=2, levels=0)
manalyze_args          = dict (search=3, truemotion=True, trymany=True, levels=0, badrange=-24)
nnedi_args             = dict (field=1, dh=True, nns=4, qual=2, etype=1, nsize=0)

### No F**king Individual Cores Within Each Function ###
class vaporcore:
      core             = vs.get_core ()
      NNEDI            = core.nnedi3.nnedi3
      KNLMeansCL       = core.knlm.KNLMeansCL
      RGB2OPP          = core.bm3d.RGB2OPP
      OPP2RGB          = core.bm3d.OPP2RGB
      BM3D             = core.bm3d.VFinal
      Aggregate        = core.bm3d.VAggregate
      Resample         = core.fmtc.resample
      BlankClip        = core.std.BlankClip
      ShufflePlanes    = core.std.ShufflePlanes
      AddBorders       = core.std.AddBorders
      Crop             = core.std.CropRel
      CropAbs          = core.std.CropAbs
      Transpose        = core.std.Transpose
      StackHorizontal  = core.std.StackHorizontal
      StackVertical    = core.std.StackVertical
      Expr             = core.std.Expr
      MaskedMerge      = core.std.MaskedMerge
      MakeDiff         = core.std.MakeDiff
      MergeDiff        = core.std.MergeDiff
      MSuper           = core.mvsf.Super
      MAnalyze         = mvmulti.Analyze
      MDegrainN        = mvmulti.DegrainN

### Helpers ###
class helpers:
      def freq_merge (low, hi, p=8):
          def gauss (src):
              upsmp    = vaporcore.Resample (src, src.width * 2, src.height * 2, kernel="gauss", a1=100, **fmtc_args)
              clip     = vaporcore.Resample (upsmp, src.width, src.height, kernel="gauss", a1=p, **fmtc_args)
              return clip
          hif          = vaporcore.MakeDiff (hi, gauss (hi))
          clip         = vaporcore.MergeDiff (gauss (low), hif)
          return clip
      def padding (src, left=0, right=0, top=0, bottom=0):
          w            = src.width
          h            = src.height
          clip         = vaporcore.Resample (src, w+left+right, h+top+bottom, -left, -top, w+left+right, h+top+bottom, kernel="point", **fmtc_args)
          return clip
      def thr_merge (flt, src, ref=None, thr=0.0009765625, elast=None):
          ref          = src if ref is None else ref
          elast        = thr / 2 if elast is None else elast
          BExp         = ["x {thr} {elast} + z - 2 {elast} * / * y {elast} z + {thr} - 2 {elast} * / * +".format (thr=thr, elast=elast)]
          BDif         = vaporcore.Expr (src, "0.0")
          PDif         = vaporcore.Expr ([flt, src], "x y - 0.0 max")
          PRef         = vaporcore.Expr ([flt, ref], "x y - 0.0 max")
          PBLD         = vaporcore.Expr ([PDif, BDif, PRef], BExp)
          NDif         = vaporcore.Expr ([flt, src], "y x - 0.0 max")
          NRef         = vaporcore.Expr ([flt, ref], "y x - 0.0 max")
          NBLD         = vaporcore.Expr ([NDif, BDif, NRef], BExp)
          BLDD         = vaporcore.MakeDiff (PBLD, NBLD)
          BLD          = vaporcore.MergeDiff (src, BLDD)
          UDN          = vaporcore.Expr ([flt, ref, BLD], ["x y - abs {thr} {elast} - > z x ?".format (thr=thr, elast=elast)])
          clip         = vaporcore.Expr ([flt, ref, UDN, src], ["x y - abs {thr} {elast} + < z a ?".format (thr=thr, elast=elast)])
          return clip
      def genblockmask (src):
          clip         = vaporcore.BlankClip (src, 24, 24, color=0.0)
          clip         = vaporcore.AddBorders (clip, 4, 4, 4, 4, color=1.0)
          clip         = vaporcore.StackHorizontal ([clip, clip, clip, clip])
          clip         = vaporcore.StackVertical ([clip, clip, clip, clip])
          clip         = vaporcore.Resample (clip, 32, 32, kernel="point", **fmtc_args)
          clip         = vaporcore.Expr (clip, ["x 0.0 > 1.0 0.0 ?"])
          clip         = vaporcore.StackHorizontal ([clip, clip, clip, clip, clip, clip, clip, clip])
          clip         = vaporcore.StackVertical ([clip, clip, clip, clip, clip, clip])
          clip         = vaporcore.StackHorizontal ([clip, clip, clip, clip, clip, clip])
          clip         = vaporcore.StackVertical ([clip, clip, clip, clip, clip])
          clip         = vaporcore.CropAbs (clip, src.width, src.height, 0, 0)
          return clip
      def NLMeans (src, a, h, wref, rclip, color):
          pad          = helpers.padding (src, a+1, a+1, a+1, a+1)
          rclip        = helpers.padding (rclip, a+1, a+1, a+1, a+1) if rclip is not None else None
          nlm          = vaporcore.KNLMeansCL (pad, d=0, a=a, s=1, h=h, cmode=color, wref=wref, rclip=rclip)
          clip         = vaporcore.Crop (nlm, a+1, a+1, a+1, a+1)
          return clip
      def genpelclip (src, src2=None, pel=4):
          u2x          = vaporcore.Transpose (vaporcore.NNEDI (vaporcore.Transpose (vaporcore.NNEDI (src, **nnedi_args)), **nnedi_args))
          u2x2         = 0 if src2 is None else vaporcore.Transpose (vaporcore.NNEDI (vaporcore.Transpose (vaporcore.NNEDI (src2, **nnedi_args)), **nnedi_args))
          u4x          = vaporcore.Transpose (vaporcore.NNEDI (vaporcore.Transpose (vaporcore.NNEDI (u2x, **nnedi_args)), **nnedi_args))
          u4x2         = 0 if src2 is None else vaporcore.Transpose (vaporcore.NNEDI (vaporcore.Transpose (vaporcore.NNEDI (u2x2, **nnedi_args)), **nnedi_args))
          dif          = 0 if src2 is None else (vaporcore.MakeDiff (u2x, u2x2) if pel == 2 else vaporcore.MakeDiff (u4x, u4x2))
          clip         = (u2x if pel == 2 else u4x) if src2 is None else dif
          return clip

### Basic Estimation ###
def Basic (src, level=1, \
           radius=6, h=6.4, pel=4, pel_precise=True, thscd1=10000, thscd2=255, \
           deblock=True, deblock_thr=0.03125, deblock_elast=0.015625, \
           lowpass=8):
    _rgb               = False
    _color             = True
    _mdg_plane         = 4
    _colorspace        = src.format.color_family
    if src.format.bits_per_sample < 32:
       raise TypeError ("Oyster.Basic: 32bits floating point precision input required!")
    if src.format.subsampling_w > 0 or src.format.subsampling_h > 0:
       raise TypeError ("Oyster.Basic: subsampled stuff not supported!")
    if _colorspace == vs.RGB:
       _rgb            = True
       src             = vaporcore.RGB2OPP (src, 1)
    if _colorspace == vs.GRAY:
       _color          = False
       _mdg_plane      = 0
       blank           = vaporcore.Expr (src, "0.5")
    else:
       blank           = vaporcore.Expr (src, ["0.5", "0.0", "0.0"])
    coarse             = helpers.NLMeans (src, 32, h, 1.0, None, _color)
    coarse             = helpers.freq_merge (src, coarse, lowpass)
    dif                = vaporcore.MakeDiff (src, coarse)
    supsrh             = vaporcore.MSuper (coarse, pelclip=helpers.genpelclip (coarse, pel=pel) if pel_precise else None, rfilter=4, pel=pel, chroma=_color, **msuper_args)
    suprdr             = vaporcore.MSuper (dif, pelclip=helpers.genpelclip (dif, pel=pel) if pel_precise else None, rfilter=2, pel=pel, chroma=_color, **msuper_args)
    vmulti             = vaporcore.MAnalyze (supsrh, overlap=2, blksize=4, divide=0, tr=radius, dct=5, chroma=_color, **manalyze_args)
    dif                = vaporcore.MDegrainN (blank, suprdr, vmulti, tr=radius, thsad=10000, thscd1=thscd1, thscd2=thscd2, plane=_mdg_plane)
    _refine_coarse     = vaporcore.MergeDiff (coarse, dif)
    if level == 1:
       def _nlm_loop (flt, n):
           window      = 32 // pow (2, n)
           flt         = _refine_coarse if n == 4 else flt
           dif         = vaporcore.MakeDiff (src, flt)
           dif         = helpers.NLMeans (dif, window, h, 1.0, flt, _color)
           fnl         = vaporcore.MergeDiff (flt, dif)
           n           = n - 1
           return fnl if n == -1 else _nlm_loop (fnl, n)
       _refine_fine    = _nlm_loop (None, 4)
       if deblock == True:
          _coarse      = helpers.thr_merge (_refine_fine, _refine_coarse, thr=deblock_thr, elast=deblock_elast)
          _mask        = helpers.genblockmask (vaporcore.ShufflePlanes (src, 0, vs.GRAY))
          _refine_fine = vaporcore.MaskedMerge (_refine_fine, _coarse, _mask, first_plane=True)
    clip               = _refine_fine if level == 1 else _refine_coarse
    clip               = vaporcore.OPP2RGB (clip, 1) if _rgb else clip
    return clip

### Final Estimation ###
def Final (src, ref, level=1, \
           radius=6, h=6.4, sigma=12.0, pel=4, pel_precise=True, thsad=2000, thscd1=10000, thscd2=255, \
           block_size=8, block_step=1, group_size=32, bm_range=24, bm_step=1, ps_num=2, ps_range=8, ps_step=1, \
           deblock=True, deblock_thr=0.03125, deblock_elast=0.015625, \
           lowpass=8):
    _rgb               = False
    _color             = True
    _mdg_plane         = 4
    _hfine             = pow (1.1988568728336214663622280225868, h)
    _mse               = sigma * 160 + thsad / 1.7097673100105595991973500117303
    _matrix            = None
    _colorspace        = src.format.color_family
    if src.format.bits_per_sample < 32:
       raise TypeError ("Oyster.Final: 32bits floating point precision input required!")
    if src.format.subsampling_w > 0 or src.format.subsampling_h > 0:
       raise TypeError ("Oyster.Final: subsampled stuff not supported!")
    if _colorspace == vs.RGB:
       _rgb            = True
       _matrix         = 100
       src             = vaporcore.RGB2OPP (src, 1)
       ref             = vaporcore.RGB2OPP (ref, 1)
    if _colorspace == vs.GRAY:
       _color          = False
       _mdg_plane      = 0
    def _nlm_loop (flt, init, src, n):
        c1             = 1.0707892518365290738330599429051
        c2             = 0.4798695862246764421520306169363
        str            = n * h / 4 + _hfine * (1 - n / 4)
        weight         = pow (c1, str * (4 - n)) - c2
        window         = 32 // pow (2, n)
        flt            = init if n == 4 else flt
        dif            = vaporcore.MakeDiff (src, flt)
        dif            = helpers.NLMeans (dif, window, str, weight, flt, _color)
        fnl            = vaporcore.MergeDiff (flt, dif)
        n              = n - 1
        return fnl if n == -1 else _nlm_loop (fnl, init, src, n)
    supsrh             = vaporcore.MSuper (ref, pelclip=helpers.genpelclip (ref, pel=pel) if pel_precise else None, rfilter=4, pel=pel, chroma=_color, **msuper_args)
    suprdr             = vaporcore.MSuper (src, pelclip=helpers.genpelclip (src, pel=pel) if pel_precise else None, rfilter=2, pel=pel, chroma=_color, **msuper_args)
    vmulti             = vaporcore.MAnalyze (supsrh, overlap=2, blksize=4, divide=0, tr=radius, dct=5, chroma=_color, **manalyze_args)
    temporal_bm        = vaporcore.MDegrainN (src, suprdr, vmulti, tr=radius, thsad=thsad, thscd1=thscd1, thscd2=thscd2, plane=_mdg_plane)
    temporal_bm        = helpers.freq_merge (src, temporal_bm, lowpass)
    if level == 1:
       tmp_fine        = _nlm_loop (None, temporal_bm, src, 4)
       if deblock == True:
          _coarse      = helpers.thr_merge (tmp_fine, temporal_bm, thr=deblock_thr, elast=deblock_elast)
          _mask        = helpers.genblockmask (vaporcore.ShufflePlanes (src, 0, vs.GRAY))
          tmp_fine     = vaporcore.MaskedMerge (tmp_fine, _coarse, _mask, first_plane=True)
       temporal_bm     = tmp_fine
    bm3d               = vaporcore.BM3D (temporal_bm, ref, radius=radius, th_mse=_mse, sigma=sigma, \
                                         block_size=block_size, block_step=block_step, group_size=group_size, bm_range=bm_range, bm_step=bm_step, \
                                         ps_num=ps_num, ps_range=ps_range, ps_step=ps_step, matrix=_matrix)
    bm3d               = vaporcore.Aggregate (bm3d, radius, 1)
    clip               = _nlm_loop (None, bm3d, temporal_bm, 4)
    clip               = vaporcore.OPP2RGB (clip, 1) if _rgb else clip
    return clip
