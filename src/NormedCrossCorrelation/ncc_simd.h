#ifndef __NCC_SIMD_H__
#define __NCC_SIMD_H__

#include "ImageData.h"

#ifdef __cplusplus
extern "C"
{
#endif
	void normed_cross_correlation_simd(struct ImageData* pattern, struct ImageData* target, struct ImageData* result);
#ifdef __cplusplus
};
#endif
#endif

