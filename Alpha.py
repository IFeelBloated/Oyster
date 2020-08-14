from vapoursynth import *

fmtc_args = dict(fulls=True, fulld=True)
nnedi_args = dict(field=1, dh=True, nns=4, qual=2, etype=1, nsize=0)

class Utilities:
      def __init__(self):
          self.Resample = core.fmtc.resample
          self.NNEDI = core.nnedi3.nnedi3
          self.Transpose = core.std.Transpose
          self.Reverse = core.std.Reverse
          self.SetFieldBased = core.std.SetFieldBased
      
      def Pad(self, src, left, right, top, bottom):
          w = src.width
          h = src.height
          return self.Resample(src, w+left+right, h+top+bottom, -left, -top, w+left+right, h+top+bottom, kernel="point", **fmtc_args)
      
      def TemporalMirror(self, src, radius):
          head = src[1:1+radius]
          head = self.Reverse(head)
          tail = src[src.num_frames-1-radius:src.num_frames-1]
          tail = self.Reverse(tail)
          return head + src + tail
 
class Clip:
      def __init__(self, node):
          self.node = node
          
      def __getattr__(self, func):
          return lambda *args, **kw: Clip(getattr(Utilities(), func)(self.node, *args, **kw))
          
def _super_(src, pel):
    src = src.Pad(32, 32, 32, 32)
    clip = src.NNEDI(**nnedi_args).Transpose().NNEDI(**nnedi_args).Transpose()
    if pel == 4:
        clip = clip.NNEDI(**nnedi_args).Transpose().NNEDI(**nnedi_args).Transpose()
    return clip

def Super(src, pel=4):
    if not isinstance(src, VideoNode):
       raise TypeError("Oyster.Super: src has to be a video clip!")
    elif src.format.sample_type != FLOAT or src.format.bits_per_sample < 32:
       raise TypeError("Oyster.Super: src has to be single precision!")  
    elif src.format.color_family != GRAY:
       raise RuntimeError("Oyster.Super: src has to be of GRAY format!")
    if not isinstance(pel, int):
       raise TypeError("Oyster.Super: pel has to be an integer!")
    elif pel != 2 and pel != 4:
       raise RuntimeError("Oyster.Super: pel has to be 2 or 4!")
    clip = _super_(Clip(src).SetFieldBased(0), pel)
    return clip.node