#ifndef __NCC_H__
#define __NCC_H__

#include "ImageData.h"
#ifdef __cplusplus
extern "C"
{
#endif
	void normed_cross_correlation(struct ImageData* pattern, struct ImageData* target, struct ImageData* result);

#ifdef __cplusplus
};
#endif
#endif
