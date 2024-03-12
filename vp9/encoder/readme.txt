This version adds pixel histogram stats, mainly two files added:
    vid_shd.h
    vid_shd.cpp

The vid_shd(histogram) class takes two frames (frame N-1, and frame N).
we call frame N-1 reference frame, frame N source frame.
It generates 3 difference histrograms for Luma and Chroma:
   1. source histogram: source frame pixel histogram, for 8-bit video, it has 256 bins;
   2. reference histogram: source frame pixel histogram, for 8-bit video, it has 256 bins;
   3. difference histogram, the histogram of "source frame - reference frame", since it could be negative, so it has 512 bins.
