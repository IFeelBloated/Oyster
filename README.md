# Oyster
©2016 IFeelBloated, Oyster Python Module for VapourSynth
## License
LGPL v2.1
## Description
Oyster is an experimental implement of the Blocking Matching concept, designed specifically for compression artifacts removal.<br />
according to [Wikipedia](https://en.wikipedia.org/wiki/Compression_artifact#Images), when performing block-based coding for quantization, several types of artifacts can appear.
- Ringing
- Contouring
- Posterizing
- Staircase noise
- Blockiness

and oyster handles 3 of them, ringing, staircase noise and blockiness

## Requirements
- [NNEDI3](https://github.com/dubhater/vapoursynth-nnedi3)
- [KNLMeansCL](https://github.com/Khanattila/KNLMeansCL)
- [BM3D](https://github.com/HomeOfVapourSynthEvolution/VapourSynth-BM3D)
- [DFTTest](https://github.com/HomeOfVapourSynthEvolution/VapourSynth-DFTTest)
- [FMTConv](https://github.com/EleonoreMizo/fmtconv)
- [MVTools (floating point ver)](https://github.com/IFeelBloated/vapoursynth-mvtools-sf/tree/master)

## Function List
- Super
- Basic
- Deringing
- Destaircase
- Deblocking

## Formats
- Bit Depth: 32bits floating point
- Color Space: Gray, RGB, YUV 4:4:4 (subsampled YUV formats are not supported)
- Scan Type: Progressive

## Notes
- **DO NOT** upsample your video to YUV 4:4:4 or RGB before processing if it's not natively full-sampled, just pass Y as a gray clip and merge the result with UV from the source clip, low-res chroma will jeopardize the correctness of weight calculation (fatal, especially to NLMeans).
- **DO NOT** crop your video before processing, it will destroy the macroblock boundary detecting.
- **QUALITY**: cutting edge
- **PERFORMANCE**: abysmal, like, literally..

## Details
### Super
Optional, it helps improve the precision of sub-pixel motion estimation and compensation, use it and get a quality boost or don't and get a performance boost
```python
Super(src, pel=4)
```
- src<br />
  clip to be processed
- pel<br />
  sub-pixel precision, could be 2 or 4, 2 = precision by half a pixel, 4 = precision by quarter a pixel.

### Basic
The basic estimation does a wild block matching based motion compensation, it removes all significant artifacts and serves as the reference to all later specific artifacts removing functions
```python
Basic(src, super=None, radius=6, pel=4, sad_me=200.0, sad_mc=2000.0)
```
- super<br />
  optional, clip generated by Oyster.Super
- radius<br />
  temporal radius, frames that fall in [current frame - radius, current frame + radius] will be referenced
- sad_me<br />
  SAD threshold of motion vector refining, refer to the MRecalculate chapter in MVTools doc for more details
- sad_mc<br />
  SAD threshold of temporal averaging, refer to the MDeGrain chapter in MVTools doc for more details

### Deringing
Deringing removes ringing (aka. mosquito noise) artifacts caused by lossy compression.<br />

workflow:
- low frequencies of the basic estimation will be replaced with those from the source clip, as ringing is obviously high frequency artifacts and we don't wanna screw irrelevant things up
- coarse refining on basic estimation by BM3D (VBasic)
- another more delicate refining by NLMeans
- BM3D (VFinal) filtering to remove trivial artifacts, basic estimation handles significant artifacts only
- replace low frequencies of the filtered with source clip low frequencies yet again
- NLMeans refining again

```python
Deringing(src, ref, radius=6, h=6.4, sigma=16.0, mse=[None, None], hard_thr=3.2, block_size=8, block_step=1, group_size=32, bm_range=24, bm_step=1, ps_num=2, ps_range=8, ps_step=1, lowpass=None)
```
- ref<br />
  clip generated by Oyster.Basic
- h<br />
  filtering strength of NLMeans refining, greater value = more relaxed refining (less possible to have residual artifacts but more detail loss)
- sigma, mse, hard_thr, block_size, block_step, group_size, bm_range, bm_step, ps_num, ps_range, ps_step<br />
  refer to BM3D doc for more details and, mse[0] is the mse value for VBasic, mse[1] for VFinal, default mse[0] = sigma * 160.0 + 1200.0, mse[1] = sigma * 120.0 + 800.0
- lowpass<br />
  controls how lowpass filter works, refer to the sstring section in DFTTest doc for more details, default = "0.0:sigma 0.48:1024.0 1.0:1024.0"

### Destaircase
block based quantization sometimes zeros out high frequency coefficients and leaves the low frequency part (almost or completely) unprocessed, which yields staircase noise, a kind of artifacts that resembles blocking, and sometimes may considered as blocking, it shows as discontinuities (aliasing) along curving edges, and Destaircase kills it

workflow:
- low frequencies of the basic estimation will be replaced with those from the source clip like Deringing
- a threshold based limiter eliminates all small differences, discontinuities are large differences apparently
- coarse BM3D refining (VBasic)
- more delicate BM3D refining (VFinal)
- replace macroblock boundaries in the source clip with the filtered result

```python
Destaircase(src, ref, radius=6, sigma=16.0, mse=[None, None], hard_thr=3.2, block_size=8, block_step=1, group_size=32, bm_range=24, bm_step=1, ps_num=2, ps_range=8, ps_step=1, thr=0.03125, elast=0.015625, lowpass=None)
```
- thr<br />
  threshold of the limiter, ranges from 0.0 (no limit) to 1.0 (no filtering), differences between basic estimation and source clip < thr will be discarded, otherwise remain unaffected.
- elast<br />
  elasticity of the threshold, ranges from 0.0 to thr.

### Deblocking
generally, Destaircase + Deringing combo is enough to typical blocking artifacts, this one works in extreme cases with severe blocking artifacts, you can actually tell this filter is fairly destructive from the previous "extreme cases" statement, **DON'T** use it unless you have to, and use it with caution

workflow:
- Make a 100% free-of-artifacts copy of the input by appending an NLMeans filtering to the basic estimation, regardless of detail loss
- coarse BM3D refining (VBasic)
- more delicate BM3D refining (VFinal)
- replace low frequency components of both source clip and basic estimation with low frequencies from the filtered result
- replace macroblock boundaries in the source clip with the basic estimation

```python
Deblocking(src, ref, radius=6, h=6.4, sigma=16.0, mse=[None, None], hard_thr=3.2, block_size=8, block_step=1, group_size=32, bm_range=24, bm_step=1, ps_num=2, ps_range=8, ps_step=1, lowpass="0.0:0.0 0.12:1024.0 1.0:1024.0")
```

## Demos
- Destaircase<br />
```python
y    = core.std.ShufflePlanes(clip, 0, vs.GRAY)
ref  = Oyster.Basic(y, Oyster.Super(y))
y    = Oyster.Destaircase(y, ref, block_step=2)
clip = core.std.ShufflePlanes([y, clip], [0, 1, 2], vs.YUV)
```
![](http://i.imgur.com/WurNke9.png)
![](http://i.imgur.com/MFA0sNY.png)
- Deringing<br />
```python
y    = core.std.ShufflePlanes(clip, 0, vs.GRAY)
ref  = Oyster.Basic(y, Oyster.Super(y))
y    = Oyster.Destaircase(y, ref, block_step=2)
y    = Oyster.Deringing(y, ref, block_step=2)
clip = core.std.ShufflePlanes([y, clip], [0, 1, 2], vs.YUV)
```
![](http://i.imgur.com/ixqv2fj.png)
![](http://i.imgur.com/u85ILWz.png)
- Deringing (severe mosquito noise)<br />
```python
y    = core.std.ShufflePlanes(clip, 0, vs.GRAY)
ref  = Oyster.Basic(y, Oyster.Super(y))
y    = Oyster.Destaircase(y, ref, block_step=2, lowpass="0.0:1024 1.0:1024")
y    = Oyster.Deringing(y, ref, sigma=24.0, h=12.8, block_step=2, lowpass="0.0:1024 1.0:1024")
clip = core.std.ShufflePlanes([y, clip], [0, 1, 2], vs.YUV)
```
![](http://i.imgur.com/jCDUuJa.png)
![](http://i.imgur.com/JqDZrtD.png)
- Deringing (H.264 compression artifacts)<br />
  *click the image and view at full size*
```python
ref = Oyster.Basic (clp, Oyster.Super (clp))
clp = Oyster.Destaircase (clp, ref, sigma=24.0, block_step=2)
clp = Oyster.Deringing (clp, ref, sigma=24.0, h=12.8, block_step=2)
```
![](http://i.imgur.com/Iw0wy79.png)
![](http://i.imgur.com/NX8ugUu.png)
- Deblocking<br />
```python
ref = Oyster.Basic (clp, Oyster.Super (clp))
clp = Oyster.Deblocking (clp, ref, block_step=2)
clp = Oyster.Deringing (clp, ref, sigma=24.0, h=10.8, block_step=2)
```
![](http://i.imgur.com/CZzS4Ci.png)
![](http://i.imgur.com/YmFQVCg.png)
![](http://i.imgur.com/kgksDfR.png)
![](http://i.imgur.com/hHbHxxM.png)
