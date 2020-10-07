from vapoursynth import *
from VaporMagik import *
import math
import os

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
    RegisterPlugin(core.raws, lambda _, FilterName: 'Raw' + FilterName)

@Inject
def Materialize(self: VideoNode, filename, message = None):
    cache_dir = os.path.dirname(__script__) + '/Cache'
    cache_path = cache_dir + '/' + filename
    width = self.get_frame(0).width
    height = self.get_frame(0).height
    fmt = 'GRAYS' if self.format.bits_per_sample == 32 else 'Y8'
    def monitor(x, y):
        if x != y:
            print(f'Frame: {x}/{y}', end = '\r')
        else:
            print('Evaluation completed')
    if not os.path.exists(cache_path):
        if message is not None:
            print(message)
        if not os.path.exists(cache_dir):
            os.mkdir(cache_dir)
        with open(cache_path, 'w+') as f:
            self.CropAbs(width = width, height = height).output(f.buffer, progress_update = monitor)
    clip = RawSource(cache_path, width, height, self.fps_num, self.fps_den, src_fmt = fmt)
    if self.width == 0 or self.height == 0:
        placeholder = clip.BlankClip(width = 1, length = 1)
        clip = [placeholder, clip].Splice(mismatch = True)
        clip = clip[1:]
    return clip

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
def OysterBasic(self: VideoNode, radius, sigma, sigma2, \
                mse, mse2, hard_thr, block_size, block_step, group_size, bm_range, bm_step, \
                pel, sad, me_sad_upperbound, me_sad_lowerbound, searchparam, \
                breakdown, checkpoints):
    spatial_args = {'block_size': block_size, \
                    'block_step': block_step, 'group_size': group_size, \
                    'bm_range': bm_range, 'bm_step': bm_step}
    temporal_args = {'me_sad_upperbound': me_sad_upperbound, 'me_sad_lowerbound': me_sad_lowerbound, \
                     'searchparam': searchparam}
    ref = self.BM3DBasic(self, radius = 0, sigma = sigma, th_mse = mse, hard_thr = hard_thr, **spatial_args)
    if breakdown.lower() == 'spatialbasic':
        return ref
    if 'spatialbasic' in checkpoints:
        ref = ref.Materialize('ob0.bin', 'Evaluating basic spatial estimation')
    clip = self.BM3DFinal(ref, ref, radius = 0, sigma = sigma2, th_mse = mse2, **spatial_args)
    if breakdown.lower() == 'spatialfinal':
        return clip
    if 'spatialfinal' in checkpoints:
        clip = clip.Materialize('ob1.bin', 'Evaluating final spatial estimation')
    clip = clip.Mirror(64, 64, 64, 64).TemporalMirror(radius)
    src = self.Mirror(64, 64, 64, 64).TemporalMirror(radius)
    superclip = clip.MSupersample(pel)
    supersrc = src.MSupersample(pel)
    if 'supersample' in checkpoints and pel != 1:
        superclip = superclip.Materialize('ob2.bin', 'Supersampling')
        supersrc = supersrc.Materialize('ob3.bin', 'Supersampling')
    degrain_super = src.MSuper(pelclip = supersrc, rfilter = 2, pel = pel, **msuper_args)
    vec = clip.ScanMotionVectors(superclip, radius, pel, **temporal_args)
    if 'motionanalysis' in checkpoints:
        vec = vec.Materialize('ob4.bin', 'Analyzing motion vector field')
    clip = clip.MDegrain(degrain_super, vec, thsad = sad, **mdegrain_args).Crop(64, 64, 64, 64)
    return clip[radius: clip.num_frames - radius]

@Inject
def OysterDeblocking(self: VideoNode, ref, radius, a, s, h, sigma, sigma2, \
                  mse, mse2, hard_thr, block_size, block_step, group_size, bm_range, bm_step, ps_num, ps_range, ps_step, \
                  sbsize, slocation, crop_left, crop_top, \
                  breakdown, checkpoints):
    bm3d_args = {'block_size': block_size, 'block_step': block_step, \
                 'group_size': group_size, 'bm_range': bm_range, \
                 'bm_step': bm_step, 'ps_num': ps_num, 'ps_range': ps_range, 'ps_step': ps_step}
    mask = self.DrawMacroblockMask(crop_left, crop_top)  
    cleansed = ref.NLMeans(radius, a, s, h)
    if breakdown.lower() == 'cleanse':
        return cleansed
    if 'cleanse' in checkpoints:
        cleansed = cleansed.Materialize('odbk0.bin', 'Cleansing')
    dif = self.MakeDiff(cleansed)
    ref_dif = dif.BM3DBasic(cleansed, radius = radius, sigma = sigma, th_mse = mse, hard_thr = hard_thr, **bm3d_args)
    if 'refinebasic' in checkpoints:
        ref_dif = ref_dif.Materialize('odbk1.bin', 'Refining (basic estimation)')
    ref = cleansed.MergeDiff(ref_dif)
    if breakdown.lower() == 'refinebasic':
        return ref
    dif = dif.BM3DFinal(ref, ref_dif, radius = radius, sigma = sigma2, th_mse = mse2, **bm3d_args)
    if 'refinefinal' in checkpoints:
        dif = dif.Materialize('odbk2.bin', 'Refining (final estimation)')
    clean = cleansed.MergeDiff(dif)
    if breakdown.lower() == 'refinefinal':
        return clean
    clip = clean.ReplaceHighFrequencyComponent(self, sbsize, slocation)
    return clip.MaskedMerge(clean, mask)

@Inject
def RestoreHighFrequencyComponent(self: VideoNode, src, h_upperbound, h_lowerbound):
    def Refine(ref, a, h):
        dif = src.MakeDiff(ref).NLMeans(0, a, 1, h, ref)
        return ref.MergeDiff(dif)
    h = CosineInterpolate(h_upperbound, h_lowerbound, 3)
    ref = Refine(self, 2, h[0])
    ref = Refine(ref, 4, h[1])
    ref = Refine(ref, 8, h[2])
    ref = Refine(ref, 16, h[3])
    return Refine(ref, 32, h[4])

@Inject
def OysterDeringing(self: VideoNode, ref, radius, sigma, sigma2, mse, mse2, \
                    hard_thr, block_size, block_step, group_size, bm_range, bm_step, ps_num, ps_range, ps_step, \
                    h_upperbound, h_lowerbound, \
                    refine_sigma, refine_sigma2, refine_mse, refine_mse2, refine_hard_thr, \
                    sbsize, slocation, breakdown, checkpoints):
    bm3d_args = {'block_size': block_size, 'block_step': block_step, \
                 'group_size': group_size, 'bm_range': bm_range, \
                 'bm_step': bm_step, 'ps_num': ps_num, 'ps_range': ps_range, 'ps_step': ps_step}
    cleansed_basic = ref.BM3DBasic(ref, radius = radius, sigma = sigma, th_mse = mse, hard_thr = hard_thr, **bm3d_args)
    if breakdown.lower() == 'cleansebasic':
        return cleansed_basic
    if 'cleansebasic' in checkpoints:
        cleansed_basic = cleansed_basic.Materialize('dr0.bin', 'Cleansing (basic estimation)')
    cleansed = ref.BM3DFinal(cleansed_basic, cleansed_basic, radius = radius, sigma = sigma2, th_mse = mse2, **bm3d_args)
    if breakdown.lower() == 'cleansefinal':
        return cleansed
    if 'cleansefinal' in checkpoints:
        cleansed = cleansed.Materialize('dr1.bin', 'Cleansing (final estimation)')
    cleansed = cleansed.RestoreHighFrequencyComponent(self, h_upperbound, h_lowerbound)
    if breakdown.lower() == 'restore':
        return cleansed
    if 'restore' in checkpoints:
        cleansed = cleansed.Materialize('dr2.bin', 'Restoring high frequency component')
    dif = self.MakeDiff(cleansed)
    ref_dif = dif.BM3DBasic(cleansed, radius = radius, sigma = refine_sigma, th_mse = refine_mse, hard_thr = refine_hard_thr, **bm3d_args)
    if 'refinebasic' in checkpoints:
        ref_dif = ref_dif.Materialize('dr3.bin', 'Refining (basic estimation)')
    ref = cleansed.MergeDiff(ref_dif)
    if breakdown.lower() == 'refinebasic':
        return ref
    dif = dif.BM3DFinal(ref, ref_dif, radius = radius, sigma = refine_sigma2, th_mse = refine_mse2, **bm3d_args)
    if 'refinefinal' in checkpoints:
        dif = dif.Materialize('dr4.bin', 'Refining (final estimation)')
    clean = cleansed.MergeDiff(dif)
    if breakdown.lower() == 'refinefinal':
        return clean
    return self.ReplaceHighFrequencyComponent(clean, sbsize, slocation)

@Inject
def OysterDecombingSpatial(self: VideoNode, p):
    clip = self
    for _ in range(p):
        clip = clip.nnedi3(**nnedi_args).resample(self.width, self.height, sx = 0, sy = -0.5, kernel = "gauss", a1 = 100, **fmtc_args)
    return clip

@Inject
def OysterDecombing(self: VideoNode, radius, p, p2, \
                    pel, sad, me_sad_upperbound, me_sad_lowerbound, searchparam, \
                    breakdown, checkpoints):
    prefiltered = self.OysterDecombingSpatial(p)
    if breakdown.lower() == 'prefilter':
        return prefiltered
    if 'prefilter' in checkpoints:
        prefiltered = prefiltered.Materialize('odcb0.bin', 'Prefiltering')
    prefiltered = prefiltered.Mirror(64, 64, 64, 64).TemporalMirror(radius)
    clip = self.Mirror(64, 64, 64, 64).TemporalMirror(radius)
    superprefiltered = prefiltered.MSupersample(pel)
    superclip = clip.MSupersample(pel)
    if 'supersample' in checkpoints and pel != 1:
        superprefiltered = superprefiltered.Materialize('odcb1.bin', 'Supersampling')
        superclip = superclip.Materialize('odcb2.bin', 'Supersampling')
    degrain_super = clip.MSuper(pelclip = superclip, rfilter = 2, pel = pel, **msuper_args)
    vec = prefiltered.ScanMotionVectors(superprefiltered, radius, pel, me_sad_upperbound, me_sad_lowerbound, searchparam)
    if 'motionanalysis' in checkpoints:
        vec = vec.Materialize('odcb3.bin', 'Analyzing motion vector field')
    clip = clip.MDegrain(degrain_super, vec, thsad = sad, **mdegrain_args).Crop(64, 64, 64, 64)
    clip = clip[radius: clip.num_frames - radius]
    if breakdown.lower() == 'compensate':
        return clip
    if 'compensate' in checkpoints:
        clip = clip.Materialize('odcb4.bin', 'Applying motion compensation')
    dif = self.MakeDiff(clip).OysterDecombingSpatial(p2)
    return clip.MergeDiff(dif)