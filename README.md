## Template matching 
ref1: https://medium.com/@bkwtwn/%E7%89%B9%E5%BE%B5%E5%8C%B9%E9%85%8D-template-matching-f2de49998dcc
ref2: https://scribblethink.org/Work/nvisionInterface/nip.html
ref3: https://github.com/opencv/opencv/blob/4.x/modules/imgproc/src/templmatch.cpp
ref4: https://github.com/ttang10/FFT_2D_CONVOLUTION/blob/master/source%20codes/fft.c


This project is to learn template matching using normalized cross-correlation(NCC)
1. src/NormedCrossCorrelation/ncc.c:     original version of calculating NCC 
2. src/NormedCrossCorrelation/ncc_sum.c: optimized version of using sum table
3. src/NormedCrossCorrelation/ncc_fft.c: optimized version of using Fast Fourier Transfrom and sum table

## external libraries list:
To read/write image file formats
1. libjpeg: https://libjpeg.sourceforge.net/
2. LibPng: https://github.com/pnggroup/libpng
3. LibTiff: http://www.libtiff.org/
