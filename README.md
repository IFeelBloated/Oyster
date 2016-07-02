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

## Requirements
- NNEDI3
- KNLMeansCL
- BM3D
- DFTTest
- FMTConv
- MVTools (floating point ver)

## Function List
- Super
- Basic
- Deringing
- Destaircase
- Deblocking

## Formats
- Bitdepth: 32bits floating point
- Color Space: Gray, RGB, YUV 4:4:4 (subsampled YUV formats are not supported)
- Scan Type: Progressive

## Notes
- DO NOT upsample your video to YUV 4:4:4 or RGB before processing if it's not natively full-sampled, just pass the luminance plane as a gray clip and merge the processed luma with the source chroma, fake 4:4:4 is toxic as the low-res chroma will jeopardize the correctness of weight calculation (especially on Pixel-Matching), and then the quality degradation on luma sets in.
- DO NOT crop your video before processing, it will destroy the macroblock boundary detecting.
- You might wanna try waifu2x instead if your video is of CGI-like content, Oyster is times slower than waifu2x and designed specifically for photographic videos.

## Details
### Super
Optional, it helps improve the precision of sub-pixel motion estimation and compensation, use it and get a quality boost or don't and get a performance boost
```python
Super (src, pel=4)
```
- src<br />
  clip to be processed
- pel<br />
  sub-pixel precision, could be 2 or 4, 2 = precision by half a pixel, 4 = precision by quarter a pixel.
