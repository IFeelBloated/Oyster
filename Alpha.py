from vapoursynth import *
from VaporMagik import *
import math

fmtc_args = dict(fulls=True, fulld=True)
msuper_args = dict(hpad=0, vpad=0, sharp=2, levels=0)
manalyze_args = dict(search=3, truemotion=False, trymany=True, levels=0, badrange=-24, divide=0, dct=0)
mrecalculate_args = dict(truemotion=False, search=3, smooth=1, divide=0, dct=0)
mdegrain_args = dict(thscd1=16320.0, thscd2=255.0)
nnedi_args = dict(field=1, dh=True, nns=4, qual=2, etype=1, nsize=0)

def CosineInterpolate(Begin, End, Step):
    def CurveStepDomain(x):
        x *= 0.5 * math.pi / (Step + 1)
        return 1 - math.cos(x)
    return [Begin + CurveStepDomain(x) * (End - Begin) for x in range(Step + 2)]

def _init():
    RegisterPlugin(core.std)
    RegisterPlugin(core.fmtc)
    RegisterPlugin(core.nnedi3)
    RegisterPlugin(core.mvsf, lambda _, FilterName: 'M' + FilterName)

@Inject
def Pad(self: VideoNode, left, right, top, bottom):
    w = self.width
    h = self.height
    return self.resample(w + left + right, h + top + bottom, -left, -top, w + left + right, h + top + bottom, kernel = "point", **fmtc_args)

@Inject
def TemporalMirror(self: VideoNode, radius):
    head = self[1: 1 + radius].Reverse()
    tail = self[self.num_frames - 1 - radius: self.num_frames - 1].Reverse()
    return head + self + tail

def _super(clip, pel):
    clip = clip.Pad(64, 64, 64, 64)
    clip = clip.nnedi3(**nnedi_args).Transpose().nnedi3(**nnedi_args).Transpose()
    if pel == 4:
        clip = clip.nnedi3(**nnedi_args).Transpose().nnedi3(**nnedi_args).Transpose()
    return clip

def _basic(clip, superclip, radius, pel, sad, me_sad_lowerbound, me_sad_upperbound):
    clip = clip.Pad(64, 64, 64, 64).TemporalMirror(radius)
    superclip = superclip.TemporalMirror(radius) if superclip is not None else None
    me_sad = CosineInterpolate(me_sad_lowerbound, me_sad_upperbound, 3)
    me_super = clip.MSuper(pelclip=superclip, rfilter=4, pel=pel, **msuper_args)
    degrain_super = clip.MSuper(pelclip=superclip, rfilter=2, pel=pel, **msuper_args)
    vec = me_super.MAnalyze(radius=radius, overlap=32, blksize=64, searchparam=4, **manalyze_args)
    vec = me_super.MRecalculate(vec, overlap=16, blksize=32, searchparam=8, thsad=me_sad[0], **mrecalculate_args)
    vec = me_super.MRecalculate(vec, overlap=8, blksize=16, searchparam=16, thsad=me_sad[1], **mrecalculate_args)
    vec = me_super.MRecalculate(vec, overlap=4, blksize=8, searchparam=32, thsad=me_sad[2], **mrecalculate_args)
    vec = me_super.MRecalculate(vec, overlap=2, blksize=4, searchparam=64, thsad=me_sad[3], **mrecalculate_args)
    vec = me_super.MRecalculate(vec, overlap=1, blksize=2, searchparam=128, thsad=me_sad[4], **mrecalculate_args)
    clip = clip.MDegrain(degrain_super, vec, thsad=sad, **mdegrain_args)
    clip = clip.Crop(64, 64, 64, 64)
    return clip[radius: clip.num_frames - radius]

def Super(clip, pel = 4):
    if not isinstance(clip, VideoNode):
       raise TypeError("Oyster.Super: clip has to be a video clip!")
    elif clip.format.sample_type != FLOAT or clip.format.bits_per_sample < 32:
       raise TypeError("Oyster.Super: clip has to be of fp32 sample type!")  
    elif clip.format.color_family != GRAY:
       raise RuntimeError("Oyster.Super: clip has to be of GRAY format!")
    if not isinstance(pel, int):
       raise TypeError("Oyster.Super: pel has to be an integer!")
    elif pel != 2 and pel != 4:
       raise RuntimeError("Oyster.Super: pel has to be 2 or 4!")
    _init()
    return _super(clip.SetFieldBased(0), pel)



