import vapoursynth as vs
import math
import mvmulti

fmtc_args                 = dict(fulls=True, fulld=True)
msuper_args               = dict(hpad=32, vpad=32, sharp=2, levels=0)
manalyze_args             = dict(search=3, truemotion=False, trymany=True, levels=0, badrange=-24, divide=0, dct=0)
mrecalculate_args         = dict(truemotion=False, search=3, smooth=1, divide=0, dct=0)
nnedi_args                = dict(field=1, dh=True, nns=4, qual=2, etype=1, nsize=0)
dfttest_args              = dict(smode=0, sosize=0, tbsize=1, tosize=0, tmode=0)

class helpers:
      def freq_merge(low, hi, sbsize, sstring):
          core            = vs.get_core()
          DFTTest         = core.dfttest.DFTTest
          MakeDiff        = core.std.MakeDiff
          MergeDiff       = core.std.MergeDiff
          hif             = MakeDiff(hi, DFTTest(hi, sbsize=sbsize, sstring=sstring, **dfttest_args))
          clip            = MergeDiff(DFTTest(low, sbsize=sbsize, sstring=sstring, **dfttest_args), hif)
          return clip
      def padding(src, left, right, top, bottom):
          core            = vs.get_core()
          Resample        = core.fmtc.resample
          w               = src.width
          h               = src.height
          clip            = Resample(src, w+left+right, h+top+bottom, -left, -top, w+left+right, h+top+bottom, kernel="point", **fmtc_args)
          return clip
      def nlmeans(src, d, a, s, h, rclip, color):
          core            = vs.get_core()
          Crop            = core.std.CropRel
          KNLMeansCL      = core.knlm.KNLMeansCL
          def duplicate(src):
              if d > 0:
                 head     = src[0] * d
                 tail     = src[src.num_frames - 1] * d
                 clip     = head + src + tail
              else:
                 clip     = src
              return clip
          pad             = helpers.padding(src, a+s, a+s, a+s, a+s)
          pad             = duplicate(pad)
          if rclip is not None:
             rclip        = helpers.padding(rclip, a+s, a+s, a+s, a+s)
             rclip        = duplicate(rclip)
          nlm             = KNLMeansCL(pad, d=d, a=a, s=s, h=h, cmode=color, wref=1.0, rclip=rclip)
          clip            = Crop(nlm, a+s, a+s, a+s, a+s)
          return clip[d:clip.num_frames - d]
      def thr_merge(flt, src, ref=None, thr=0.0009765625, elast=None):
          core            = vs.get_core()
          Expr            = core.std.Expr
          MakeDiff        = core.std.MakeDiff
          MergeDiff       = core.std.MergeDiff
          ref             = src if ref is None else ref
          elast           = thr / 2 if elast is None else elast
          BExp            = ["x {thr} {elast} + z - 2 {elast} * / * y {elast} z + {thr} - 2 {elast} * / * +".format(thr=thr, elast=elast)]
          BDif            = Expr(src, "0.0")
          PDif            = Expr([flt, src], "x y - 0.0 max")
          PRef            = Expr([flt, ref], "x y - 0.0 max")
          PBLD            = Expr([PDif, BDif, PRef], BExp)
          NDif            = Expr([flt, src], "y x - 0.0 max")
          NRef            = Expr([flt, ref], "y x - 0.0 max")
          NBLD            = Expr([NDif, BDif, NRef], BExp)
          BLDD            = MakeDiff(PBLD, NBLD)
          BLD             = MergeDiff(src, BLDD)
          UDN             = Expr([flt, ref, BLD], ["x y - abs {thr} {elast} - > z x ?".format(thr=thr, elast=elast)])
          clip            = Expr([flt, ref, UDN, src], ["x y - abs {thr} {elast} + < z a ?".format(thr=thr, elast=elast)])
          return clip
      def genblockmask(src):
          core            = vs.get_core()
          Resample        = core.fmtc.resample
          BlankClip       = core.std.BlankClip
          AddBorders      = core.std.AddBorders
          StackHorizontal = core.std.StackHorizontal
          StackVertical   = core.std.StackVertical
          Expr            = core.std.Expr
          CropAbs         = core.std.CropAbs
          clip            = BlankClip(src, 24, 24, color=0.0)
          clip            = AddBorders(clip, 4, 4, 4, 4, color=1.0)
          clip            = StackHorizontal([clip, clip, clip, clip])
          clip            = StackVertical([clip, clip, clip, clip])
          clip            = Resample(clip, 32, 32, kernel="point", **fmtc_args)
          clip            = Expr(clip, ["x 0.0 > 1.0 0.0 ?"])
          clip            = StackHorizontal([clip, clip, clip, clip, clip, clip, clip, clip])
          clip            = StackVertical([clip, clip, clip, clip, clip, clip])
          clip            = StackHorizontal([clip, clip, clip, clip, clip, clip])
          clip            = StackVertical([clip, clip, clip, clip, clip])
          clip            = StackHorizontal([clip, clip, clip, clip, clip, clip])
          clip            = StackVertical([clip, clip, clip, clip, clip])
          clip            = CropAbs(clip, src.width, src.height, 0, 0)
          return clip

class internal:
      def super(src, pel):
          core            = vs.get_core()
          NNEDI           = core.nnedi3.nnedi3
          Transpose       = core.std.Transpose
          clip            = Transpose(NNEDI(Transpose(NNEDI(src, **nnedi_args)), **nnedi_args))
          if pel == 4:
             clip         = Transpose(NNEDI(Transpose(NNEDI(clip, **nnedi_args)), **nnedi_args))
          return clip
      def basic(src, super, radius, pel, sad_me, sad_mc, color):
          core            = vs.get_core()
          MSuper          = core.mvsf.Super
          MAnalyze        = mvmulti.Analyze
          MRecalculate    = mvmulti.Recalculate
          MDegrainN       = mvmulti.DegrainN
          plane           = 4 if color else 0
          supersoft       = MSuper(src, pelclip=super, rfilter=4, pel=pel, chroma=color, **msuper_args)
          supersharp      = MSuper(src, pelclip=super, rfilter=2, pel=pel, chroma=color, **msuper_args)
          vmulti          = MAnalyze(supersoft, tr=radius, chroma=color, overlap=16, blksize=32, **manalyze_args)
          vmulti          = MRecalculate(supersoft, vmulti, tr=radius, chroma=color, overlap=8, blksize=16, thsad=sad_me, **mrecalculate_args)
          vmulti          = MRecalculate(supersharp, vmulti, tr=radius, chroma=color, overlap=4, blksize=8, thsad=sad_me, **mrecalculate_args)
          vmulti          = MRecalculate(supersharp, vmulti, tr=radius, chroma=color, overlap=2, blksize=4, thsad=sad_me, **mrecalculate_args)
          clip            = MDegrainN(src, supersharp, vmulti, tr=radius, thsad=sad_mc, thscd1=10000.0, thscd2=255.0, plane=plane)
          return clip
      def deringing(src, ref, radius, h, sigma, \
                    mse, hard_thr, block_size, block_step, group_size, bm_range, bm_step, ps_num, ps_range, ps_step, \
                    lowpass, color, matrix):
          core            = vs.get_core()
          BMBasic         = core.bm3d.VBasic
          BMFinal         = core.bm3d.VFinal
          Aggregate       = core.bm3d.VAggregate
          MakeDiff        = core.std.MakeDiff
          MergeDiff       = core.std.MergeDiff
          c1              = 1.1396386205122096184557327136584
          c2              = 4.8995241035176996103733445761166
          strength        = [h, None, None]
          strength[1]     =((math.exp(c1 * h) - 1.0) /(math.pow(h, h) / math.gamma(h + 1.0))) / c2
          def loop(flt, init, src, n):
              strength[2] = n * strength[0] / 4 + strength[1] *(1 - n / 4)
              window      = 32 // pow(2, n)
              flt         = init if n == 4 else flt
              dif         = MakeDiff(src, flt)
              dif         = helpers.nlmeans(dif, 0, window, 1, strength[2], flt, color)
              fnl         = MergeDiff(flt, dif)
              n          -= 1
              return fnl if n == -1 else loop(fnl, init, src, n)
          ref             = helpers.freq_merge(src, ref, block_size // 2 * 2 + 1, lowpass)
          dif             = MakeDiff(src, ref)
          dif             = BMBasic(dif, ref, radius=radius, th_mse=mse[0], hard_thr=hard_thr, sigma=sigma, \
                                    block_size=block_size, block_step=block_step, group_size=group_size, bm_range=bm_range, bm_step=bm_step, \
                                    ps_num=ps_num, ps_range=ps_range, ps_step=ps_step, matrix=matrix)
          dif             = Aggregate(dif, radius, 1)
          ref             = MergeDiff(ref, dif)
          refined         = loop(None, ref, src, 4)
          bm3d            = BMFinal(refined, ref, radius=radius, th_mse=mse[1], sigma=sigma, \
                                    block_size=block_size, block_step=block_step, group_size=group_size, bm_range=bm_range, bm_step=bm_step, \
                                    ps_num=ps_num, ps_range=ps_range, ps_step=ps_step, matrix=matrix)
          bm3d            = Aggregate(bm3d, radius, 1)
          bm3d            = helpers.freq_merge(refined, bm3d, block_size // 2 * 2 + 1, lowpass)
          clip            = loop(None, bm3d, refined, 4)
          return clip
      def destaircase(src, ref, radius, sigma, \
                      mse, hard_thr, block_size, block_step, group_size, bm_range, bm_step, ps_num, ps_range, ps_step, \
                      thr, elast, lowpass, matrix):
          core            = vs.get_core()
          BMBasic         = core.bm3d.VBasic
          BMFinal         = core.bm3d.VFinal
          Aggregate       = core.bm3d.VAggregate
          MakeDiff        = core.std.MakeDiff
          MergeDiff       = core.std.MergeDiff
          MaskedMerge     = core.std.MaskedMerge
          ShufflePlanes   = core.std.ShufflePlanes
          mask            = helpers.genblockmask(ShufflePlanes(src, 0, vs.GRAY))
          ref             = helpers.freq_merge(src, ref, block_size // 2 * 2 + 1, lowpass)
          ref             = helpers.thr_merge(src, ref, thr=thr, elast=elast)
          dif             = MakeDiff(src, ref)
          dif             = BMBasic(dif, ref, radius=radius, th_mse=mse[0], hard_thr=hard_thr, sigma=sigma, \
                                    block_size=block_size, block_step=block_step, group_size=group_size, bm_range=bm_range, bm_step=bm_step, \
                                    ps_num=ps_num, ps_range=ps_range, ps_step=ps_step, matrix=matrix)
          dif             = Aggregate(dif, radius, 1)
          ref             = MergeDiff(ref, dif)
          dif             = MakeDiff(src, ref)
          dif             = BMFinal(dif, ref, radius=radius, th_mse=mse[1], sigma=sigma, \
                                    block_size=block_size, block_step=block_step, group_size=group_size, bm_range=bm_range, bm_step=bm_step, \
                                    ps_num=ps_num, ps_range=ps_range, ps_step=ps_step, matrix=matrix)
          dif             = Aggregate(dif, radius, 1)
          ref             = MergeDiff(ref, dif)
          clip            = MaskedMerge(src, ref, mask, first_plane=True)
          return clip
      def deblocking(src, ref, radius, h, sigma, \
                     mse, hard_thr, block_size, block_step, group_size, bm_range, bm_step, ps_num, ps_range, ps_step, \
                     lowpass, color, matrix):
          core            = vs.get_core()
          BMBasic         = core.bm3d.VBasic
          BMFinal         = core.bm3d.VFinal
          Aggregate       = core.bm3d.VAggregate
          MakeDiff        = core.std.MakeDiff
          MergeDiff       = core.std.MergeDiff
          MaskedMerge     = core.std.MaskedMerge
          ShufflePlanes   = core.std.ShufflePlanes
          mask            = helpers.genblockmask(ShufflePlanes(src, 0, vs.GRAY))
          cleansed        = helpers.nlmeans(ref, radius, block_size, 4, h, ref, color)
          dif             = MakeDiff(ref, cleansed)
          dif             = BMBasic(dif, cleansed, radius=radius, th_mse=mse[0], hard_thr=hard_thr, sigma=sigma, \
                                    block_size=block_size, block_step=block_step, group_size=group_size, bm_range=bm_range, bm_step=bm_step, \
                                    ps_num=ps_num, ps_range=ps_range, ps_step=ps_step, matrix=matrix)
          dif             = Aggregate(dif, radius, 1)
          cleansed        = MergeDiff(cleansed, dif)
          dif             = MakeDiff(ref, cleansed)
          dif             = BMFinal(dif, cleansed, radius=radius, th_mse=mse[1], sigma=sigma, \
                                    block_size=block_size, block_step=block_step, group_size=group_size, bm_range=bm_range, bm_step=bm_step, \
                                    ps_num=ps_num, ps_range=ps_range, ps_step=ps_step, matrix=matrix)
          dif             = Aggregate(dif, radius, 1)
          cleansed        = MergeDiff(cleansed, dif)
          ref             = helpers.freq_merge(cleansed, ref, block_size // 2 * 2 + 1, lowpass)
          src             = helpers.freq_merge(cleansed, src, block_size // 2 * 2 + 1, lowpass)
          clip            = MaskedMerge(src, ref, mask, first_plane=True)
          return clip

def Super(src, pel=4):
    core                  = vs.get_core()
    RGB2OPP               = core.bm3d.RGB2OPP
    SetFieldBased         = core.std.SetFieldBased
    if not isinstance(src, vs.VideoNode):
       raise TypeError("Oyster.Super: src has to be a video clip!")
    elif src.format.sample_type != vs.FLOAT or src.format.bits_per_sample < 32:
       raise TypeError("Oyster.Super: the sample type of src has to be single precision!")
    elif src.format.subsampling_w > 0 or src.format.subsampling_h > 0:
       raise RuntimeError("Oyster.Super: subsampled stuff not supported!")
    if not isinstance(pel, int):
       raise TypeError("Oyster.Super: pel has to be an integer!")
    elif pel != 2 and pel != 4:
       raise RuntimeError("Oyster.Super: pel has to be 2 or 4!")
    src                   = SetFieldBased(src, 0)
    colorspace            = src.format.color_family
    if colorspace == vs.RGB:
       src                = RGB2OPP(src, 1)
    clip                  = internal.super(src, pel)
    return clip

def Basic(src, super=None, radius=6, pel=4, sad_me=200.0, sad_mc=2000.0):
    core                  = vs.get_core()
    RGB2OPP               = core.bm3d.RGB2OPP
    OPP2RGB               = core.bm3d.OPP2RGB
    SetFieldBased         = core.std.SetFieldBased
    if not isinstance(src, vs.VideoNode):
       raise TypeError("Oyster.Basic: src has to be a video clip!")
    elif src.format.sample_type != vs.FLOAT or src.format.bits_per_sample < 32:
       raise TypeError("Oyster.Basic: the sample type of src has to be single precision!")
    elif src.format.subsampling_w > 0 or src.format.subsampling_h > 0:
       raise RuntimeError("Oyster.Basic: subsampled stuff not supported!")
    if not isinstance(super, vs.VideoNode) and super is not None:
       raise TypeError("Oyster.Basic: super has to be a video clip or None!")
    elif super is not None:
       if super.format.sample_type != vs.FLOAT or super.format.bits_per_sample < 32 or super.format.subsampling_w > 0 or super.format.subsampling_h > 0:
          raise RuntimeError("Oyster.Basic: corrupted super clip!")
    if not isinstance(radius, int):
       raise TypeError("Oyster.Basic: radius has to be an integer!")
    elif radius < 1:
       raise RuntimeError("Oyster.Basic: radius has to be greater than 0!")
    if not isinstance(pel, int):
       raise TypeError("Oyster.Basic: pel has to be an integer!")
    elif pel != 1 and pel != 2 and pel != 4:
       raise RuntimeError("Oyster.Basic: pel has to be 1, 2 or 4!")
    if not isinstance(sad_me, float) and not isinstance(sad_me, int):
       raise TypeError("Oyster.Basic: sad_me has to be a real number!")
    elif sad_me <= 0:
       raise RuntimeError("Oyster.Basic: sad_me has to be greater than 0!")
    if not isinstance(sad_mc, float) and not isinstance(sad_mc, int):
       raise TypeError("Oyster.Basic: sad_mc has to be a real number!")
    elif sad_mc <= 0:
       raise RuntimeError("Oyster.Basic: sad_mc has to be greater than 0!")
    color                 = True
    rgb                   = False
    colorspace            = src.format.color_family
    if colorspace == vs.RGB:
       src                = RGB2OPP(src, 1)
       rgb                = True
    if colorspace == vs.GRAY:
       color              = False
    src                   = SetFieldBased(src, 0)
    super                 = SetFieldBased(super, 0) if super is not None else None
    clip                  = internal.basic(src, super, radius, pel, sad_me, sad_mc, color)
    clip                  = OPP2RGB(clip, 1) if rgb else clip
    return clip

def Deringing(src, ref, radius=6, h=6.4, sigma=16.0, \
              mse=[None, None], hard_thr=3.2, block_size=8, block_step=1, group_size=32, bm_range=24, bm_step=1, ps_num=2, ps_range=8, ps_step=1, \
              lowpass=None):
    core                  = vs.get_core()
    RGB2OPP               = core.bm3d.RGB2OPP
    OPP2RGB               = core.bm3d.OPP2RGB
    SetFieldBased         = core.std.SetFieldBased
    if not isinstance(src, vs.VideoNode):
       raise TypeError("Oyster.Deringing: src has to be a video clip!")
    elif src.format.sample_type != vs.FLOAT or src.format.bits_per_sample < 32:
       raise TypeError("Oyster.Deringing: the sample type of src has to be single precision!")
    elif src.format.subsampling_w > 0 or src.format.subsampling_h > 0:
       raise RuntimeError("Oyster.Deringing: subsampled stuff not supported!")
    if not isinstance(ref, vs.VideoNode):
       raise TypeError("Oyster.Deringing: ref has to be a video clip!")
    elif ref.format.sample_type != vs.FLOAT or ref.format.bits_per_sample < 32:
       raise TypeError("Oyster.Deringing: the sample type of ref has to be single precision!")
    elif ref.format.subsampling_w > 0 or ref.format.subsampling_h > 0:
       raise RuntimeError("Oyster.Deringing: subsampled stuff not supported!")
    if not isinstance(radius, int):
       raise TypeError("Oyster.Deringing: radius has to be an integer!")
    elif radius < 1:
       raise RuntimeError("Oyster.Deringing: radius has to be greater than 0!")
    if not isinstance(h, float) and not isinstance(h, int):
       raise TypeError("Oyster.Deringing: h has to be a real number!")
    elif h <= 0:
       raise RuntimeError("Oyster.Deringing: h has to be greater than 0!")
    if not isinstance(mse, list):
       raise TypeError("Oyster.Deringing: mse parameter has to be an array!")
    elif len(mse) != 2:
       raise RuntimeError("Oyster.Deringing: mse parameter has to contain 2 elements exactly!")
    for i in range(2):
        if not isinstance(mse[i], float) and not isinstance(mse[i], int) and mse[i] is not None:
           raise TypeError("Oyster.Deringing: elements in mse must be real numbers or None!")
    if not isinstance(lowpass, str) and lowpass is not None:
       raise TypeError("Oyster.Deringing: lowpass has to be a string or None!")
    rgb                   = False
    color                 = True
    mse[0]                = sigma * 160.0 + 1200.0 if mse[0] is None else mse[0]
    mse[1]                = sigma * 120.0 + 800.0 if mse[1] is None else mse[1]
    lowpass               = "0.0:{sigma} 0.48:1024.0 1.0:1024.0".format(sigma=sigma) if lowpass is None else lowpass
    matrix                = None
    colorspace            = src.format.color_family
    if colorspace == vs.RGB:
       rgb                = True
       matrix             = 100
       src                = RGB2OPP(src, 1)
       ref                = RGB2OPP(ref, 1)
    if colorspace == vs.GRAY:
       color              = False
    src                   = SetFieldBased(src, 0)
    ref                   = SetFieldBased(ref, 0)
    clip                  = internal.deringing(src, ref, radius, h, sigma, \
                                               mse, hard_thr, block_size, block_step, group_size, bm_range, bm_step, ps_num, ps_range, ps_step, \
                                               lowpass, color, matrix)
    clip                  = OPP2RGB(clip, 1) if rgb else clip
    return clip

def Destaircase(src, ref, radius=6, sigma=16.0, \
                mse=[None, None], hard_thr=3.2, block_size=8, block_step=1, group_size=32, bm_range=24, bm_step=1, ps_num=2, ps_range=8, ps_step=1, \
                thr=0.03125, elast=0.015625, lowpass=None):
    core                  = vs.get_core()
    RGB2OPP               = core.bm3d.RGB2OPP
    OPP2RGB               = core.bm3d.OPP2RGB
    SetFieldBased         = core.std.SetFieldBased
    if not isinstance(src, vs.VideoNode):
       raise TypeError("Oyster.Destaircase: src has to be a video clip!")
    elif src.format.sample_type != vs.FLOAT or src.format.bits_per_sample < 32:
       raise TypeError("Oyster.Destaircase: the sample type of src has to be single precision!")
    elif src.format.subsampling_w > 0 or src.format.subsampling_h > 0:
       raise RuntimeError("Oyster.Destaircase: subsampled stuff not supported!")
    if not isinstance(ref, vs.VideoNode):
       raise TypeError("Oyster.Destaircase: ref has to be a video clip!")
    elif ref.format.sample_type != vs.FLOAT or ref.format.bits_per_sample < 32:
       raise TypeError("Oyster.Destaircase: the sample type of ref has to be single precision!")
    elif ref.format.subsampling_w > 0 or ref.format.subsampling_h > 0:
       raise RuntimeError("Oyster.Destaircase: subsampled stuff not supported!")
    if not isinstance(radius, int):
       raise TypeError("Oyster.Destaircase: radius has to be an integer!")
    elif radius < 1:
       raise RuntimeError("Oyster.Destaircase: radius has to be greater than 0!")
    if not isinstance(mse, list):
       raise TypeError("Oyster.Destaircase: mse parameter has to be an array!")
    elif len(mse) != 2:
       raise RuntimeError("Oyster.Destaircase: mse parameter has to contain 2 elements exactly!")
    for i in range(2):
        if not isinstance(mse[i], float) and not isinstance(mse[i], int) and mse[i] is not None:
           raise TypeError("Oyster.Destaircase: elements in mse must be real numbers or None!")
    if not isinstance(thr, float) and not isinstance(thr, int):
       raise TypeError("Oyster.Destaircase: thr has to be a real number!")
    elif thr < 0 or thr > 1:
       raise RuntimeError("Oyster.Destaircase: thr has to fall in [0, 1]!")
    if not isinstance(elast, float) and not isinstance(elast, int):
       raise TypeError("Oyster.Destaircase: elast has to be a real number!")
    elif elast < 0 or elast > thr:
       raise RuntimeError("Oyster.Destaircase: elast has to fall in [0, thr]!")
    if not isinstance(lowpass, str) and lowpass is not None:
       raise TypeError("Oyster.Destaircase: lowpass has to be a string or None!")
    rgb                   = False
    mse[0]                = sigma * 160.0 + 1200.0 if mse[0] is None else mse[0]
    mse[1]                = sigma * 120.0 + 800.0 if mse[1] is None else mse[1]
    lowpass               = "0.0:{sigma} 0.48:1024.0 1.0:1024.0".format(sigma=sigma) if lowpass is None else lowpass
    matrix                = None
    colorspace            = src.format.color_family
    if colorspace == vs.RGB:
       rgb                = True
       matrix             = 100
       src                = RGB2OPP(src, 1)
       ref                = RGB2OPP(ref, 1)
    src                   = SetFieldBased(src, 0)
    ref                   = SetFieldBased(ref, 0)
    clip                  = internal.destaircase(src, ref, radius, sigma, \
                                                 mse, hard_thr, block_size, block_step, group_size, bm_range, bm_step, ps_num, ps_range, ps_step, \
                                                 thr, elast, lowpass, matrix)
    clip                  = OPP2RGB(clip, 1) if rgb else clip
    return clip

def Deblocking(src, ref, radius=6, h=6.4, sigma=16.0, \
               mse=[None, None], hard_thr=3.2, block_size=8, block_step=1, group_size=32, bm_range=24, bm_step=1, ps_num=2, ps_range=8, ps_step=1, \
               lowpass="0.0:0.0 0.12:1024.0 1.0:1024.0"):
    core                  = vs.get_core()
    RGB2OPP               = core.bm3d.RGB2OPP
    OPP2RGB               = core.bm3d.OPP2RGB
    SetFieldBased         = core.std.SetFieldBased
    if not isinstance(src, vs.VideoNode):
       raise TypeError("Oyster.Deblocking: src has to be a video clip!")
    elif src.format.sample_type != vs.FLOAT or src.format.bits_per_sample < 32:
       raise TypeError("Oyster.Deblocking: the sample type of src has to be single precision!")
    elif src.format.subsampling_w > 0 or src.format.subsampling_h > 0:
       raise RuntimeError("Oyster.Deblocking: subsampled stuff not supported!")
    if not isinstance(ref, vs.VideoNode):
       raise TypeError("Oyster.Deblocking: ref has to be a video clip!")
    elif ref.format.sample_type != vs.FLOAT or ref.format.bits_per_sample < 32:
       raise TypeError("Oyster.Deblocking: the sample type of ref has to be single precision!")
    elif ref.format.subsampling_w > 0 or ref.format.subsampling_h > 0:
       raise RuntimeError("Oyster.Deblocking: subsampled stuff not supported!")
    if not isinstance(radius, int):
       raise TypeError("Oyster.Deblocking: radius has to be an integer!")
    elif radius < 1:
       raise RuntimeError("Oyster.Deblocking: radius has to be greater than 0!")
    if not isinstance(h, float) and not isinstance(h, int):
       raise TypeError("Oyster.Deblocking: h has to be a real number!")
    elif h <= 0:
       raise RuntimeError("Oyster.Deblocking: h has to be greater than 0!")
    if not isinstance(mse, list):
       raise TypeError("Oyster.Deblocking: mse parameter has to be an array!")
    elif len(mse) != 2:
       raise RuntimeError("Oyster.Deblocking: mse parameter has to contain 2 elements exactly!")
    for i in range(2):
        if not isinstance(mse[i], float) and not isinstance(mse[i], int) and mse[i] is not None:
           raise TypeError("Oyster.Deblocking: elements in mse must be real numbers or None!")
    if not isinstance(lowpass, str):
       raise TypeError("Oyster.Deblocking: lowpass has to be a string!")
    rgb                   = False
    color                 = True
    mse[0]                = sigma * 160.0 + 1200.0 if mse[0] is None else mse[0]
    mse[1]                = sigma * 120.0 + 800.0 if mse[1] is None else mse[1]
    matrix                = None
    colorspace            = src.format.color_family
    if colorspace == vs.RGB:
       rgb                = True
       matrix             = 100
       src                = RGB2OPP(src, 1)
       ref                = RGB2OPP(ref, 1)
    if colorspace == vs.GRAY:
       color              = False
    src                   = SetFieldBased(src, 0)
    ref                   = SetFieldBased(ref, 0)
    clip                  = internal.deblocking(src, ref, radius, h, sigma, \
                                                mse, hard_thr, block_size, block_step, group_size, bm_range, bm_step, ps_num, ps_range, ps_step, \
                                                lowpass, color, matrix)
    clip                  = OPP2RGB(clip, 1) if rgb else clip
    return clip
