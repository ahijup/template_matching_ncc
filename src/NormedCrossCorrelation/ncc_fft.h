#ifndef __NCC_FFT_H__
#define __NCC_FFT_H__

#include "ImageData.h"

#ifdef __cplusplus
extern "C"
{
#endif
	void test_fft2d(struct ImageData* src);
	void test_image_sum();
	void normed_cross_correlation_fft(struct ImageData* pattern, struct ImageData* target, struct ImageData* result);
	void test_cross_correlation_fft(struct ImageData* pattern, struct ImageData* target);
#ifdef __cplusplus
};
#endif

#endif