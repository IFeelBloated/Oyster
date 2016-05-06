import vapoursynth as vs
import mvmulti

### Global Settings ###
fmtc_args                 = dict (fulls=True, fulld=True)
msuper_args               = dict (hpad=32, vpad=32, sharp=2, levels=0)
manalyze_args             = dict (search=3, truemotion=True, trymany=True, levels=0, badrange=-24)
nnedi_args                = dict (field=1, dh=True, nns=4, qual=2, etype=1, nsize=0)

### Helpers ###
class helpers:
      def freq_merge (low, hi, p=8):
          core            = vs.get_core ()
          Resample        = core.fmtc.resample
          MakeDiff        = core.std.MakeDiff
          MergeDiff       = core.std.MergeDiff
          def gauss (src):
              upsmp       = Resample (src, src.width * 2, src.height * 2, kernel="gauss", a1=100, **fmtc_args)
              clip        = Resample (upsmp, src.width, src.height, kernel="gauss", a1=p, **fmtc_args)
              return clip
          hif             = MakeDiff (hi, gauss (hi))
          clip            = MergeDiff (gauss (low), hif)
          return clip
      def padding (src, left=0, right=0, top=0, bottom=0):
          core            = vs.get_core ()
          Resample        = core.fmtc.resample
          w               = src.width
          h               = src.height
          clip            = Resample (src, w+left+right, h+top+bottom, -left, -top, w+left+right, h+top+bottom, kernel="point", **fmtc_args)
          return clip
      def thr_merge (flt, src, ref=None, thr=0.0009765625, elast=None):
          core            = vs.get_core ()
          Expr            = core.std.Expr
          MakeDiff        = core.std.MakeDiff
          MergeDiff       = core.std.MergeDiff
          ref             = src if ref is None else ref
          elast           = thr / 2 if elast is None else elast
          BExp            = ["x {thr} {elast} + z - 2 {elast} * / * y {elast} z + {thr} - 2 {elast} * / * +".format (thr=thr, elast=elast)]
          BDif            = Expr (src, "0.0")
          PDif            = Expr ([flt, src], "x y - 0.0 max")
          PRef            = Expr ([flt, ref], "x y - 0.0 max")
          PBLD            = Expr ([PDif, BDif, PRef], BExp)
          NDif            = Expr ([flt, src], "y x - 0.0 max")
          NRef            = Expr ([flt, ref], "y x - 0.0 max")
          NBLD            = Expr ([NDif, BDif, NRef], BExp)
          BLDD            = MakeDiff (PBLD, NBLD)
          BLD             = MergeDiff (src, BLDD)
          UDN             = Expr ([flt, ref, BLD], ["x y - abs {thr} {elast} - > z x ?".format (thr=thr, elast=elast)])
          clip            = Expr ([flt, ref, UDN, src], ["x y - abs {thr} {elast} + < z a ?".format (thr=thr, elast=elast)])
          return clip
      def genblockmask (src):
          core            = vs.get_core ()
          Resample        = core.fmtc.resample
          BlankClip       = core.std.BlankClip
          AddBorders      = core.std.AddBorders
          StackHorizontal = core.std.StackHorizontal
          StackVertical   = core.std.StackVertical
          Expr            = core.std.Expr
          CropAbs         = core.std.CropAbs
          clip            = BlankClip (src, 24, 24, color=0.0)
          clip            = AddBorders (clip, 4, 4, 4, 4, color=1.0)
          clip            = StackHorizontal ([clip, clip, clip, clip])
          clip            = StackVertical ([clip, clip, clip, clip])
          clip            = Resample (clip, 32, 32, kernel="point", **fmtc_args)
          clip            = Expr (clip, ["x 0.0 > 1.0 0.0 ?"])
          clip            = StackHorizontal ([clip, clip, clip, clip, clip, clip, clip, clip])
          clip            = StackVertical ([clip, clip, clip, clip, clip, clip])
          clip            = StackHorizontal ([clip, clip, clip, clip, clip, clip])
          clip            = StackVertical ([clip, clip, clip, clip, clip])
          clip            = CropAbs (clip, src.width, src.height, 0, 0)
          return clip
      def NLMeans (src, a, h, wref, rclip, color):
          core            = vs.get_core ()
          Crop            = core.std.CropRel
          KNLMeansCL      = core.knlm.KNLMeansCL
          pad             = helpers.padding (src, a+1, a+1, a+1, a+1)
          rclip           = helpers.padding (rclip, a+1, a+1, a+1, a+1) if rclip is not None else None
          nlm             = KNLMeansCL (pad, d=0, a=a, s=1, h=h, cmode=color, wref=wref, rclip=rclip)
          clip            = Crop (nlm, a+1, a+1, a+1, a+1)
          return clip
      def genpelclip (src, src2=None, pel=4):
          core            = vs.get_core ()
          NNEDI           = core.nnedi3.nnedi3
          Transpose       = core.std.Transpose
          MakeDiff        = core.std.MakeDiff
          MergeDiff       = core.std.MergeDiff
          u2x             = Transpose (NNEDI (Transpose (NNEDI (src, **nnedi_args)), **nnedi_args))
          u2x2            = 0 if src2 is None else Transpose (NNEDI (Transpose (NNEDI (src2, **nnedi_args)), **nnedi_args))
          u4x             = Transpose (NNEDI (Transpose (NNEDI (u2x, **nnedi_args)), **nnedi_args))
          u4x2            = 0 if src2 is None else Transpose (NNEDI (Transpose (NNEDI (u2x2, **nnedi_args)), **nnedi_args))
          dif             = 0 if src2 is None else (MakeDiff (u2x, u2x2) if pel == 2 else MakeDiff (u4x, u4x2))
          clip            = (u2x if pel == 2 else u4x) if src2 is None else dif
          return clip

### Motion Estimation ###
def Search (src, radius=6, pel=4, pel_precise=True, satd=True):
    core                  = vs.get_core ()
    RGB2OPP               = core.bm3d.RGB2OPP
    MSuper                = core.mvsf.Super
    MAnalyze              = mvmulti.Analyze
    _color                = True
    _colorspace           = src.format.color_family
    if src.format.bits_per_sample < 32:
       raise TypeError ("Oyster.Search: 32bits floating point precision input required!")
    if src.format.subsampling_w > 0 or src.format.subsampling_h > 0:
       raise TypeError ("Oyster.Search: subsampled stuff not supported!")
    if _colorspace == vs.RGB:
       src                = RGB2OPP (src, 1)
    if _colorspace == vs.GRAY:
       _color             = False
    supsrh                = MSuper (src, pelclip=helpers.genpelclip (src, pel=pel) if pel_precise else None, rfilter=4, pel=pel, chroma=_color, **msuper_args)
    vmulti                = MAnalyze (supsrh, overlap=2, blksize=4, divide=0, tr=radius, dct=5 if satd else 0, chroma=_color, **manalyze_args)
    return vmulti

### Basic Estimation ###
def Basic (src, vec, level=1, \
           radius=6, h=6.4, sigma=12.0, pel=4, pel_precise=True, thscd1=10000, thscd2=255, \
           block_size=8, block_step=1, group_size=32, bm_range=24, bm_step=1, ps_num=2, ps_range=8, ps_step=1, \
           deblock=True, deblock_thr=0.03125, deblock_elast=0.015625, \
           lowpass=8):
    core                  = vs.get_core ()
    BM3D                  = core.bm3d.VBasic
    Aggregate             = core.bm3d.VAggregate
    RGB2OPP               = core.bm3d.RGB2OPP
    OPP2RGB               = core.bm3d.OPP2RGB
    Expr                  = core.std.Expr
    MakeDiff              = core.std.MakeDiff
    MergeDiff             = core.std.MergeDiff
    MaskedMerge           = core.std.MaskedMerge
    ShufflePlanes         = core.std.ShufflePlanes
    MSuper                = core.mvsf.Super
    MDegrainN             = mvmulti.DegrainN
    _rgb                  = False
    _color                = True
    _mdg_plane            = 4
    _mse                  = sigma * 160 + 1200
    _matrix               = None
    _colorspace           = src.format.color_family
    if src.format.bits_per_sample < 32:
       raise TypeError ("Oyster.Basic: 32bits floating point precision input required!")
    if src.format.subsampling_w > 0 or src.format.subsampling_h > 0:
       raise TypeError ("Oyster.Basic: subsampled stuff not supported!")
    if _colorspace == vs.RGB:
       _rgb               = True
       _matrix            = 100
       src                = RGB2OPP (src, 1)
    if _colorspace == vs.GRAY:
       _color             = False
       _mdg_plane         = 0
       blank              = Expr (src, "0.5")
    else:
       blank              = Expr (src, ["0.5", "0.0", "0.0"])
    def _nlm_loop (flt, init, src, n):
        window            = 32 // pow (2, n)
        flt               = init if n == 4 else flt
        dif               = MakeDiff (src, flt)
        dif               = helpers.NLMeans (dif, window, h, 1.0, flt, _color)
        fnl               = MergeDiff (flt, dif)
        n                 = n - 1
        return fnl if n == -1 else _nlm_loop (fnl, init, src, n)
    coarse                = helpers.NLMeans (src, 32, h, 1.0, None, _color)
    coarse                = helpers.freq_merge (src, coarse, lowpass) if lowpass != 0 else coarse
    dif                   = MakeDiff (src, coarse)
    suprdr                = MSuper (dif, pelclip=helpers.genpelclip (dif, pel=pel) if pel_precise else None, rfilter=2, pel=pel, chroma=_color, **msuper_args)
    dif                   = MDegrainN (blank, suprdr, vec, tr=radius, thsad=10000, thscd1=thscd1, thscd2=thscd2, plane=_mdg_plane)
    temporal_bm           = MergeDiff (coarse, dif)
    if level == 1:
       tmp_fine           = _nlm_loop (None, temporal_bm, src, 4)
       if deblock == True:
          _coarse         = helpers.thr_merge (tmp_fine, temporal_bm, thr=deblock_thr, elast=deblock_elast)
          _mask           = helpers.genblockmask (ShufflePlanes (src, 0, vs.GRAY))
          tmp_fine        = MaskedMerge (tmp_fine, _coarse, _mask, first_plane=True)
       temporal_bm        = tmp_fine
    bm3d                  = BM3D (temporal_bm, src, radius=radius, th_mse=_mse, sigma=sigma, \
                                  block_size=block_size, block_step=block_step, group_size=group_size, bm_range=bm_range, bm_step=bm_step, \
                                  ps_num=ps_num, ps_range=ps_range, ps_step=ps_step, matrix=_matrix)
    bm3d                  = Aggregate (bm3d, radius, 1)
    clip                  = _nlm_loop (None, bm3d, temporal_bm, 4)
    clip                  = OPP2RGB (clip, 1) if _rgb else clip
    return clip

### Final Estimation ###
def Final (src, ref, vec, level=1, \
           radius=6, h=6.4, sigma=12.0, pel=4, pel_precise=True, thsad=2000, thscd1=10000, thscd2=255, \
           block_size=8, block_step=1, group_size=32, bm_range=24, bm_step=1, ps_num=2, ps_range=8, ps_step=1, \
           deblock=True, deblock_thr=0.03125, deblock_elast=0.015625, \
           lowpass=8):
    core                  = vs.get_core ()
    BM3D                  = core.bm3d.VFinal
    Aggregate             = core.bm3d.VAggregate
    RGB2OPP               = core.bm3d.RGB2OPP
    OPP2RGB               = core.bm3d.OPP2RGB
    Expr                  = core.std.Expr
    MakeDiff              = core.std.MakeDiff
    MergeDiff             = core.std.MergeDiff
    MaskedMerge           = core.std.MaskedMerge
    ShufflePlanes         = core.std.ShufflePlanes
    MSuper                = core.mvsf.Super
    MDegrainN             = mvmulti.DegrainN
    _rgb                  = False
    _color                = True
    _mdg_plane            = 4
    _hfine                = pow (1.1988568728336214663622280225868, h)
    _mse                  = sigma * 120 + thsad / 1.7097673100105595991973500117303
    _matrix               = None
    _colorspace           = src.format.color_family
    if src.format.bits_per_sample < 32:
       raise TypeError ("Oyster.Final: 32bits floating point precision input required!")
    if src.format.subsampling_w > 0 or src.format.subsampling_h > 0:
       raise TypeError ("Oyster.Final: subsampled stuff not supported!")
    if _colorspace == vs.RGB:
       _rgb               = True
       _matrix            = 100
       src                = RGB2OPP (src, 1)
       ref                = RGB2OPP (ref, 1)
    if _colorspace == vs.GRAY:
       _color             = False
       _mdg_plane         = 0
    def _nlm_loop (flt, init, src, n):
        c1                = 1.0707892518365290738330599429051
        c2                = 0.4798695862246764421520306169363
        str               = n * h / 4 + _hfine * (1 - n / 4)
        weight            = pow (c1, str * (4 - n)) - c2
        window            = 32 // pow (2, n)
        flt               = init if n == 4 else flt
        dif               = MakeDiff (src, flt)
        dif               = helpers.NLMeans (dif, window, str, weight, flt, _color)
        fnl               = MergeDiff (flt, dif)
        n                 = n - 1
        return fnl if n == -1 else _nlm_loop (fnl, init, src, n)
    suprdr                = MSuper (src, pelclip=helpers.genpelclip (src, pel=pel) if pel_precise else None, rfilter=2, pel=pel, chroma=_color, **msuper_args)
    temporal_bm           = MDegrainN (src, suprdr, vec, tr=radius, thsad=thsad, thscd1=thscd1, thscd2=thscd2, plane=_mdg_plane)
    temporal_bm           = helpers.freq_merge (src, temporal_bm, lowpass) if lowpass != 0 else temporal_bm
    if level == 1:
       tmp_fine           = _nlm_loop (None, temporal_bm, src, 4)
       if deblock == True:
          _coarse         = helpers.thr_merge (tmp_fine, temporal_bm, thr=deblock_thr, elast=deblock_elast)
          _mask           = helpers.genblockmask (ShufflePlanes (src, 0, vs.GRAY))
          tmp_fine        = MaskedMerge (tmp_fine, _coarse, _mask, first_plane=True)
       temporal_bm        = tmp_fine
    bm3d                  = BM3D (temporal_bm, ref, radius=radius, th_mse=_mse, sigma=sigma, \
                                  block_size=block_size, block_step=block_step, group_size=group_size, bm_range=bm_range, bm_step=bm_step, \
                                  ps_num=ps_num, ps_range=ps_range, ps_step=ps_step, matrix=_matrix)
    bm3d                  = Aggregate (bm3d, radius, 1)
    clip                  = _nlm_loop (None, bm3d, temporal_bm, 4)
    clip                  = OPP2RGB (clip, 1) if _rgb else clip
    return clip
