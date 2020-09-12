from vapoursynth import *
from VaporMagik import *
import math

fmtc_args = {'fulls': True, 'fulld': True}
msuper_args = {'hpad': 0, 'vpad': 0, 'sharp': 2, 'levels': 0}
manalyze_args = {'search': 3, 'truemotion': False, 'trymany': True, 'levels': 0, 'badrange': -24, 'divide': 0, 'dct': 0}
mrecalculate_args = {'truemotion': False, 'search': 3, 'smooth': 1, 'divide': 0, 'dct': 0}
mdegrain_args = {'thscd1': 16320.0, 'thscd2': 255.0}
nnedi_args = {'field': 1, 'dh': True, 'nns': 4, 'qual': 2, 'etype': 1, 'nsize': 0}
dfttest_args = {'smode': 0, 'sosize': 0, 'tbsize': 1, 'tosize': 0, 'tmode': 0}

def CosineInterpolate(Begin, End, Step):
    def CurveStepDomain(x):
        x *= 0.5 * math.pi / (Step + 1)
        return 1 - math.cos(x)
    return [Begin + CurveStepDomain(x) * (End - Begin) for x in range(Step + 2)]

def ConvexAttenuate(Value, Scale, Multiplier):
    """
    lim x * log(1 + 1 / x) as x->inf = 1
    lim x * log(1 + 1 / x) as x->0 = 0
    """
    ScaledValue = math.pow(Scale * Value, Multiplier)
    Coefficient = ScaledValue * math.log(1 + 1 / ScaledValue)
    return Coefficient * Value

def ConcaveAttenuate(Value, Alpha, Beta):
    """
    lim log(1 + x) / x as x->inf = 0
    lim log(1 + x) / x as x->0 = 1
    """
    ScaledValues = [Alpha * Value, Beta * Value]
    Coefficients = [math.log(1 + x) / x for x in ScaledValues]
    return (Coefficients[0] + Coefficients[1]) * Value / 2

def _init():
    RegisterPlugin(core.std)
    RegisterPlugin(core.fmtc)
    RegisterPlugin(core.nnedi3)
    RegisterPlugin(core.knlm)
    RegisterPlugin(core.dfttest)
    RegisterPlugin(core.mvsf, lambda _, FilterName: 'M' + FilterName)
    RegisterPlugin(core.bm3d, lambda _, FilterName: 'BM3D' + FilterName + 'Native')

@Inject
def Mirror(self: VideoNode, left, right, top, bottom):
    clip = self
    vertical_filler = clip.FlipVertical()
    if top > 0:
        top_filler = vertical_filler.Crop(0, 0, vertical_filler.height - top - 1, 1)
        clip = [top_filler, clip].StackVertical()
    if bottom > 0:
        bottom_filler = vertical_filler.Crop(0, 0, 1, vertical_filler.height - bottom - 1)
        clip = [clip, bottom_filler].StackVertical()
    horizontal_filler = clip.FlipHorizontal()
    if left > 0:
        left_filler = horizontal_filler.Crop(horizontal_filler.width - left - 1, 1, 0, 0)
        clip = [left_filler, clip].StackHorizontal()
    if right > 0:
        right_filler = horizontal_filler.Crop(1, horizontal_filler.width - right - 1, 0, 0)
        clip = [clip, right_filler].StackHorizontal()
    return clip

@Inject
def TemporalMirror(self: VideoNode, radius):
    if radius > 0:
        head = self[1: 1 + radius].Reverse()
        tail = self[self.num_frames - 1 - radius: self.num_frames - 1].Reverse()
        return head + self + tail
    else:
        return self

@Inject
def NLMeans(self: VideoNode, d, a, s, h, rclip = None):
    clip = self.Mirror(a + s, a + s, a + s, a + s).TemporalMirror(d)
    rclip = rclip.Mirror(a + s, a + s, a + s, a + s).TemporalMirror(d) if rclip is not None else None
    clip = clip.KNLMeansCL(d = d, a = a, s = s, h = h, channels = 'Y', wref = 1.0, rclip = rclip)
    clip = clip.Crop(a + s, a + s, a + s, a + s)
    return clip[d: clip.num_frames - d]

@Inject
def BM3DBasic(self: VideoNode, ref, **kw):
    block_size = kw['block_size']
    radius = kw['radius']
    clip = self.Mirror(block_size, block_size, block_size, block_size).TemporalMirror(radius)
    ref = ref.Mirror(block_size, block_size, block_size, block_size).TemporalMirror(radius) if ref is not None else None
    if radius > 0:
        clip = clip.BM3DVBasicNative(ref, **kw).BM3DVAggregateNative(radius, 1)
    else:
        del kw['radius']
        clip = clip.BM3DBasicNative(ref, **kw)
    clip = clip.Crop(block_size, block_size, block_size, block_size)
    return clip[radius: clip.num_frames - radius]

@Inject
def BM3DFinal(self: VideoNode, ref, wref, **kw):
    block_size = kw['block_size']
    radius = kw['radius']
    clip = self.Mirror(block_size, block_size, block_size, block_size).TemporalMirror(radius)
    ref = ref.Mirror(block_size, block_size, block_size, block_size).TemporalMirror(radius)
    wref = wref.Mirror(block_size, block_size, block_size, block_size).TemporalMirror(radius) if wref is not None else None
    if radius > 0:
        clip = clip.BM3DVFinalNative(ref, wref, **kw).BM3DVAggregateNative(radius, 1)
    else:
        del kw['radius']
        clip = clip.BM3DFinalNative(ref, wref, **kw)
    clip = clip.Crop(block_size, block_size, block_size, block_size)
    return clip[radius: clip.num_frames - radius]

@Inject
def DrawMacroblockMask(self: VideoNode, left = 0, top = 0):
    ref_width = self.width + left
    ref_height = self.height + top
    clip = self.BlankClip(24, 24, color = 0.0)
    clip = clip.AddBorders(4, 4, 4, 4, color = 1.0)
    clip = [clip, clip, clip, clip].StackHorizontal()
    clip = [clip, clip, clip, clip].StackVertical()
    clip = clip.resample(32, 32, kernel = "point", **fmtc_args)
    clip = clip.Expr("x 0 > 1 0 ?")
    h_extend = [clip] * (ref_width // 32 + 1)
    clip = h_extend.StackHorizontal()
    v_extend = [clip] * (ref_height // 32 + 1)
    clip = v_extend.StackVertical()
    clip = clip.CropAbs(ref_width, ref_height, 0, 0)
    return clip.Crop(left, 0, top, 0)

@Inject
def ReplaceHighFrequencyComponent(self: VideoNode, ref, sbsize, slocation):
    clip = self.Mirror(sbsize, sbsize, sbsize, sbsize)
    ref = ref.Mirror(sbsize, sbsize, sbsize, sbsize)
    low_freq = clip.DFTTest(sbsize = sbsize, slocation = slocation, **dfttest_args)
    high_freq = ref.MakeDiff(ref.DFTTest(sbsize = sbsize, slocation = slocation, **dfttest_args))
    clip = low_freq.MergeDiff(high_freq)
    return clip.Crop(sbsize, sbsize, sbsize, sbsize)

@Inject
def ScanMotionVectors(self: VideoNode, superclip, radius, pel, me_sad_upperbound, me_sad_lowerbound, searchparam):
    me_sad = CosineInterpolate(me_sad_lowerbound, me_sad_upperbound, 3)
    clip = self.MSuper(pelclip = superclip, rfilter = 4, pel = pel, **msuper_args)
    vec = clip.MAnalyze(radius = radius, overlap = 32, blksize = 64, searchparam = searchparam, **manalyze_args)
    vec = clip.MRecalculate(vec, overlap = 16, blksize = 32, searchparam = searchparam * 2, thsad = me_sad[0], **mrecalculate_args)
    vec = clip.MRecalculate(vec, overlap = 8, blksize = 16, searchparam = searchparam * 4, thsad = me_sad[1], **mrecalculate_args)
    vec = clip.MRecalculate(vec, overlap = 4, blksize = 8, searchparam = searchparam * 8, thsad = me_sad[2], **mrecalculate_args)
    vec = clip.MRecalculate(vec, overlap = 2, blksize = 4, searchparam = searchparam * 16, thsad = me_sad[3], **mrecalculate_args)
    vec = clip.MRecalculate(vec, overlap = 1, blksize = 2, searchparam = searchparam * 32, thsad = me_sad[4], **mrecalculate_args)
    return vec

@Inject
def MSupersample(self: VideoNode, pel):
    if pel == 1:
        return None
    clip = self.nnedi3(**nnedi_args).Transpose().nnedi3(**nnedi_args).Transpose()
    if pel == 4:
        clip = clip.nnedi3(**nnedi_args).Transpose().nnedi3(**nnedi_args).Transpose()
    return clip

@Inject
def OysterBasicSpatial(self: VideoNode, sigma, sigma2, mse, mse2, **kw):
    ref = self.BM3DBasic(self, radius = 0, sigma = sigma, th_mse = mse, **kw)
    if 'hard_thr' in kw:
        del kw['hard_thr']
    return self.BM3DFinal(ref, ref, radius = 0, sigma = sigma2, th_mse = mse2, **kw)

@Inject
def OysterBasicTemporal(self: VideoNode, src, superclip, supersrc, radius, pel, sad, **kw):
    degrain_super = src.MSuper(pelclip = supersrc, rfilter = 2, pel = pel, **msuper_args)
    vec = self.ScanMotionVectors(superclip, radius, pel, **kw)
    return self.MDegrain(degrain_super, vec, thsad = sad, **mdegrain_args)

@Inject
def OysterBasic(self: VideoNode, radius, sigma, sigma2, \
                mse, mse2, hard_thr, block_size, block_step, group_size, bm_range, bm_step, \
                pel, sad, me_sad_upperbound, me_sad_lowerbound, searchparam):
    spatial_args = {'hard_thr': hard_thr, 'block_size': block_size, \
                    'block_step': block_step, 'group_size': group_size, \
                    'bm_range': bm_range, 'bm_step': bm_step}
    temporal_args = {'me_sad_upperbound': me_sad_upperbound, 'me_sad_lowerbound': me_sad_lowerbound, \
                     'searchparam': searchparam}
    clip = self.OysterBasicSpatial(sigma, sigma2, mse, mse2, **spatial_args).Mirror(64, 64, 64, 64).TemporalMirror(radius)
    src = self.Mirror(64, 64, 64, 64).TemporalMirror(radius)
    superclip = clip.MSupersample(pel)
    supersrc = src.MSupersample(pel)
    clip = clip.OysterBasicTemporal(src, superclip, supersrc, radius, pel, sad, **temporal_args).Crop(64, 64, 64, 64)
    return clip[radius: clip.num_frames - radius]

@Inject
def OysterDeblock(self: VideoNode, ref, radius, a, s, h, sigma, sigma2, \
                  mse, mse2, hard_thr, block_size, block_step, group_size, bm_range, bm_step, ps_num, ps_range, ps_step, \
                  sbsize, slocation, crop_left, crop_top):
    bm3d_args = {'block_size': block_size, 'block_step': block_step, \
                 'group_size': group_size, 'bm_range': bm_range, \
                 'bm_step': bm_step, 'ps_num': ps_num, 'ps_range': ps_range, 'ps_step': ps_step}
    mask = self.DrawMacroblockMask(crop_left, crop_top)  
    cleansed = ref.NLMeans(radius, a, s, h)
    dif = self.MakeDiff(cleansed)
    ref_dif = dif.BM3DBasic(cleansed, radius = radius, sigma = sigma, th_mse = mse, hard_thr = hard_thr, **bm3d_args)
    ref = cleansed.MergeDiff(ref_dif)
    dif = dif.BM3DFinal(ref, ref_dif, radius = radius, sigma = sigma2, th_mse = mse2, **bm3d_args)
    clean = cleansed.MergeDiff(dif)
    clip = clean.ReplaceHighFrequencyComponent(self, sbsize, slocation)
    return clip.MaskedMerge(clean, mask)