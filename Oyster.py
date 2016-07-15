import vapoursynth as vs
import math
import mvmulti

fmtc_args                 = dict (fulls=True, fulld=True)
msuper_args               = dict (hpad=32, vpad=32, sharp=2, levels=0)
manalyze_args             = dict (search=3, truemotion=False, trymany=True, levels=0, badrange=-24, divide=0, dct=0)
mrecalculate_args         = dict (truemotion=False, search=3, smooth=1, divide=0, dct=0)
nnedi_args                = dict (field=1, dh=True, nns=4, qual=2, etype=1, nsize=0)
dfttest_args              = dict (smode=0, sosize=0, tbsize=1, tosize=0, tmode=0)

class helpers:
      def freq_merge (low, hi, sbsize, sstring):
          core            = vs.get_core ()
          DFTTest         = core.dfttest.DFTTest
          MakeDiff        = core.std.MakeDiff
          MergeDiff       = core.std.MergeDiff
          hif             = MakeDiff (hi, DFTTest (hi, sbsize=sbsize, sstring=sstring, **dfttest_args))
          clip            = MergeDiff (DFTTest (low, sbsize=sbsize, sstring=sstring, **dfttest_args), hif)
          return clip
      def padding (src, left=0, right=0, top=0, bottom=0):
          core            = vs.get_core ()
          Resample        = core.fmtc.resample
          w               = src.width
          h               = src.height
          clip            = Resample (src, w+left+right, h+top+bottom, -left, -top, w+left+right, h+top+bottom, kernel="point", **fmtc_args)
          return clip
      def NLMeans (src, d, a, s, h, rclip, color):
          core            = vs.get_core ()
          Crop            = core.std.CropRel
          KNLMeansCL      = core.knlm.KNLMeansCL
          pad             = helpers.padding (src, a+s, a+s, a+s, a+s)
          rclip           = helpers.padding (rclip, a+s, a+s, a+s, a+s) if rclip is not None else None
          nlm             = KNLMeansCL (pad, d=d, a=a, s=s, h=h, cmode=color, wref=1.0, rclip=rclip)
          clip            = Crop (nlm, a+s, a+s, a+s, a+s)
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
          clip            = StackHorizontal ([clip, clip, clip, clip, clip, clip])
          clip            = StackVertical ([clip, clip, clip, clip, clip])
          clip            = CropAbs (clip, src.width, src.height, 0, 0)
          return clip

def Super (src, pel=4):
    core            = vs.get_core ()
    NNEDI           = core.nnedi3.nnedi3
    Transpose       = core.std.Transpose
    if src.format.bits_per_sample < 32:
       raise TypeError ("Oyster.Super: 32bits floating point precision input required!")
    if src.format.subsampling_w > 0 or src.format.subsampling_h > 0:
       raise TypeError ("Oyster.Super: subsampled stuff not supported!")
    u2x             = Transpose (NNEDI (Transpose (NNEDI (src, **nnedi_args)), **nnedi_args))
    u4x             = Transpose (NNEDI (Transpose (NNEDI (u2x, **nnedi_args)), **nnedi_args))
    clip            = u2x if pel == 2 else u4x
    return clip

def Basic (src, super=None, radius=6, pel=4, sad_me=200.0, sad_mc=2000.0, scd1=10000.0, scd2=255.0):
    core                  = vs.get_core ()
    RGB2OPP               = core.bm3d.RGB2OPP
    OPP2RGB               = core.bm3d.OPP2RGB
    MSuper                = core.mvsf.Super
    MAnalyze              = mvmulti.Analyze
    MRecalculate          = mvmulti.Recalculate
    MDegrainN             = mvmulti.DegrainN
    _mdg_plane            = 4
    _color                = True
    _rgb                  = False
    _colorspace           = src.format.color_family
    if src.format.bits_per_sample < 32:
       raise TypeError ("Oyster.Basic: 32bits floating point precision input required!")
    if src.format.subsampling_w > 0 or src.format.subsampling_h > 0:
       raise TypeError ("Oyster.Basic: subsampled stuff not supported!")
    if _colorspace == vs.RGB:
       src                = RGB2OPP (src, 1)
       super              = RGB2OPP (super, 1) if super is not None else None
       _rgb               = True
    if _colorspace == vs.GRAY:
       _color             = False
       _mdg_plane         = 0
    supersoft             = MSuper (src, pelclip=super, rfilter=4, pel=pel, chroma=_color, **msuper_args)
    supersharp            = MSuper (src, pelclip=super, rfilter=2, pel=pel, chroma=_color, **msuper_args)
    vmulti                = MAnalyze (supersoft, tr=radius, chroma=_color, overlap=16, blksize=32, **manalyze_args)
    vmulti                = MRecalculate (supersoft, vmulti, tr=radius, chroma=_color, overlap=8, blksize=16, thsad=sad_me, **mrecalculate_args)
    vmulti                = MRecalculate (supersharp, vmulti, tr=radius, chroma=_color, overlap=4, blksize=8, thsad=sad_me, **mrecalculate_args)
    vmulti                = MRecalculate (supersharp, vmulti, tr=radius, chroma=_color, overlap=2, blksize=4, thsad=sad_me, **mrecalculate_args)
    clip                  = MDegrainN (src, supersharp, vmulti, tr=radius, thsad=sad_mc, thscd1=scd1, thscd2=scd2, plane=_mdg_plane)
    clip                  = OPP2RGB (clip, 1) if _rgb else clip
    return clip

def Deringing (src, ref, radius=6, h=6.4, sigma=16.0, \
               mse=[None, None], hard_thr=3.2, block_size=8, block_step=1, group_size=32, bm_range=24, bm_step=1, ps_num=2, ps_range=8, ps_step=1, \
               lowpass=None):
    core                  = vs.get_core ()
    BMBasic               = core.bm3d.VBasic
    BMFinal               = core.bm3d.VFinal
    Aggregate             = core.bm3d.VAggregate
    RGB2OPP               = core.bm3d.RGB2OPP
    OPP2RGB               = core.bm3d.OPP2RGB
    MakeDiff              = core.std.MakeDiff
    MergeDiff             = core.std.MergeDiff
    c1                    = 1.1396386205122096184557327136584
    c2                    = 4.8995241035176996103733445761166
    _hfine                = ((math.exp (c1 * h) - 1.0) / (math.pow (h, h) / math.gamma (h + 1.0))) / c2
    _rgb                  = False
    _color                = True
    mse[0]                = sigma * 160.0 + 1200.0 if mse[0] is None else mse[0]
    mse[1]                = sigma * 120.0 + 800.0 if mse[1] is None else mse[1]
    lowpass               = "0.0:{sigma} 0.48:1024.0 1.0:1024.0".format (sigma=sigma) if lowpass is None else lowpass
    _matrix               = None
    _colorspace           = src.format.color_family
    if src.format.bits_per_sample < 32:
       raise TypeError ("Oyster.Deringing: 32bits floating point precision input required!")
    if src.format.subsampling_w > 0 or src.format.subsampling_h > 0:
       raise TypeError ("Oyster.Deringing: subsampled stuff not supported!")
    if _colorspace == vs.RGB:
       _rgb               = True
       _matrix            = 100
       src                = RGB2OPP (src, 1)
       ref                = RGB2OPP (ref, 1)
    if _colorspace == vs.GRAY:
       _color             = False
    def _nlm_loop (flt, init, src, n): 
        str               = n * h / 4 + _hfine * (1 - n / 4)
        window            = 32 // pow (2, n)
        long_term_mean    = [2, 3, 4, 3, 4]
        short_term_mean   = [1, 1, 1, 2, 2]
        flt               = init if n == 4 else flt
        dif               = MakeDiff (src, flt)
        long_term_dif     = helpers.NLMeans (dif, 0, window, long_term_mean[n], str, flt, _color)
        short_term_dif    = helpers.NLMeans (dif, 0, window, short_term_mean[n], str, flt, _color)
        dif               = helpers.freq_merge (long_term_dif, short_term_dif, block_size // 2 * 2 + 1, lowpass)
        fnl               = MergeDiff (flt, dif)
        n                 = n - 1
        return fnl if n == -1 else _nlm_loop (fnl, init, src, n)
    ref                   = helpers.freq_merge (src, ref, block_size // 2 * 2 + 1, lowpass)
    dif                   = MakeDiff (src, ref)
    dif                   = BMBasic (dif, ref, radius=radius, th_mse=mse[0], hard_thr=hard_thr, sigma=sigma, \
                                     block_size=block_size, block_step=block_step, group_size=group_size, bm_range=bm_range, bm_step=bm_step, \
                                     ps_num=ps_num, ps_range=ps_range, ps_step=ps_step, matrix=_matrix)
    dif                   = Aggregate (dif, radius, 1)
    ref                   = MergeDiff (ref, dif)
    refined               = _nlm_loop (None, ref, src, 4)
    bm3d                  = BMFinal (refined, ref, radius=radius, th_mse=mse[1], sigma=sigma, \
                                     block_size=block_size, block_step=block_step, group_size=group_size, bm_range=bm_range, bm_step=bm_step, \
                                     ps_num=ps_num, ps_range=ps_range, ps_step=ps_step, matrix=_matrix)
    bm3d                  = Aggregate (bm3d, radius, 1)
    bm3d                  = helpers.freq_merge (refined, bm3d, block_size // 2 * 2 + 1, lowpass)
    clip                  = _nlm_loop (None, bm3d, refined, 4)
    clip                  = OPP2RGB (clip, 1) if _rgb else clip
    return clip

def Destaircase (src, ref, radius=6, sigma=16.0, \
                 mse=[None, None], hard_thr=3.2, block_size=8, block_step=1, group_size=32, bm_range=24, bm_step=1, ps_num=2, ps_range=8, ps_step=1, \
                 thr=0.03125, elast=0.015625, lowpass=None):
    core                  = vs.get_core ()
    BMBasic               = core.bm3d.VBasic
    BMFinal               = core.bm3d.VFinal
    Aggregate             = core.bm3d.VAggregate
    RGB2OPP               = core.bm3d.RGB2OPP
    OPP2RGB               = core.bm3d.OPP2RGB
    MakeDiff              = core.std.MakeDiff
    MergeDiff             = core.std.MergeDiff
    MaskedMerge           = core.std.MaskedMerge
    ShufflePlanes         = core.std.ShufflePlanes
    _rgb                  = False
    _color                = True
    mse[0]                = sigma * 160.0 + 1200.0 if mse[0] is None else mse[0]
    mse[1]                = sigma * 120.0 + 800.0 if mse[1] is None else mse[1]
    lowpass               = "0.0:{sigma} 0.48:1024.0 1.0:1024.0".format (sigma=sigma) if lowpass is None else lowpass
    _matrix               = None
    _colorspace           = src.format.color_family
    if src.format.bits_per_sample < 32:
       raise TypeError ("Oyster.Destaircase: 32bits floating point precision input required!")
    if src.format.subsampling_w > 0 or src.format.subsampling_h > 0:
       raise TypeError ("Oyster.Destaircase: subsampled stuff not supported!")
    if _colorspace == vs.RGB:
       _rgb               = True
       _matrix            = 100
       src                = RGB2OPP (src, 1)
       ref                = RGB2OPP (ref, 1)
    if _colorspace == vs.GRAY:
       _color             = False
    mask                  = helpers.genblockmask (ShufflePlanes (src, 0, vs.GRAY))
    ref                   = helpers.freq_merge (src, ref, block_size // 2 * 2 + 1, lowpass)
    ref                   = helpers.thr_merge (src, ref, thr=thr, elast=elast)
    dif                   = MakeDiff (src, ref)
    dif                   = BMBasic (dif, ref, radius=radius, th_mse=mse[0], hard_thr=hard_thr, sigma=sigma, \
                                     block_size=block_size, block_step=block_step, group_size=group_size, bm_range=bm_range, bm_step=bm_step, \
                                     ps_num=ps_num, ps_range=ps_range, ps_step=ps_step, matrix=_matrix)
    dif                   = Aggregate (dif, radius, 1)
    ref                   = MergeDiff (ref, dif)
    dif                   = MakeDiff (src, ref)
    dif                   = BMFinal (dif, ref, radius=radius, th_mse=mse[1], sigma=sigma, \
                                     block_size=block_size, block_step=block_step, group_size=group_size, bm_range=bm_range, bm_step=bm_step, \
                                     ps_num=ps_num, ps_range=ps_range, ps_step=ps_step, matrix=_matrix)
    dif                   = Aggregate (dif, radius, 1)
    ref                   = MergeDiff (ref, dif)
    clip                  = MaskedMerge (src, ref, mask, first_plane=True)
    clip                  = OPP2RGB (clip, 1) if _rgb else clip
    return clip

def Deblocking (src, ref, radius=6, h=6.4, sigma=16.0, \
                mse=[None, None], hard_thr=3.2, block_size=8, block_step=1, group_size=32, bm_range=24, bm_step=1, ps_num=2, ps_range=8, ps_step=1, \
                lowpass="0.0:0.0 0.12:1024.0 1.0:1024.0"):
    core                  = vs.get_core ()
    BMBasic               = core.bm3d.VBasic
    BMFinal               = core.bm3d.VFinal
    Aggregate             = core.bm3d.VAggregate
    RGB2OPP               = core.bm3d.RGB2OPP
    OPP2RGB               = core.bm3d.OPP2RGB
    MakeDiff              = core.std.MakeDiff
    MergeDiff             = core.std.MergeDiff
    MaskedMerge           = core.std.MaskedMerge
    ShufflePlanes         = core.std.ShufflePlanes
    _rgb                  = False
    _color                = True
    mse[0]                = sigma * 160.0 + 1200.0 if mse[0] is None else mse[0]
    mse[1]                = sigma * 120.0 + 800.0 if mse[1] is None else mse[1]
    _matrix               = None
    _colorspace           = src.format.color_family
    if src.format.bits_per_sample < 32:
       raise TypeError ("Oyster.Deblocking: 32bits floating point precision input required!")
    if src.format.subsampling_w > 0 or src.format.subsampling_h > 0:
       raise TypeError ("Oyster.Deblocking: subsampled stuff not supported!")
    if _colorspace == vs.RGB:
       _rgb               = True
       _matrix            = 100
       src                = RGB2OPP (src, 1)
       ref                = RGB2OPP (ref, 1)
    if _colorspace == vs.GRAY:
       _color             = False
    mask                  = helpers.genblockmask (ShufflePlanes (src, 0, vs.GRAY))
    cleansed              = helpers.NLMeans (ref, radius, block_size, 4, h, ref, _color)
    dif                   = MakeDiff (ref, cleansed)
    dif                   = BMBasic (dif, cleansed, radius=radius, th_mse=mse[0], hard_thr=hard_thr, sigma=sigma, \
                                     block_size=block_size, block_step=block_step, group_size=group_size, bm_range=bm_range, bm_step=bm_step, \
                                     ps_num=ps_num, ps_range=ps_range, ps_step=ps_step, matrix=_matrix)
    dif                   = Aggregate (dif, radius, 1)
    cleansed              = MergeDiff (cleansed, dif)
    dif                   = MakeDiff (ref, cleansed)
    dif                   = BMFinal (dif, cleansed, radius=radius, th_mse=mse[1], sigma=sigma, \
                                     block_size=block_size, block_step=block_step, group_size=group_size, bm_range=bm_range, bm_step=bm_step, \
                                     ps_num=ps_num, ps_range=ps_range, ps_step=ps_step, matrix=_matrix)
    dif                   = Aggregate (dif, radius, 1)
    cleansed              = MergeDiff (cleansed, dif)
    ref                   = helpers.freq_merge (cleansed, ref, block_size // 2 * 2 + 1, lowpass)
    src                   = helpers.freq_merge (cleansed, src, block_size // 2 * 2 + 1, lowpass)
    clip                  = MaskedMerge (src, ref, mask, first_plane=True)
    clip                  = OPP2RGB (clip, 1) if _rgb else clip
    return clip
