#ifndef __NCC_COMMON_H__
#define __NCC_COMMON_H__

#include "ImageData.h"

#ifdef __cplusplus
extern "C"
{
#endif
	void normalize_cross_correlation(struct ImageData* pattern, struct ImageData* target, struct ImageData* cross_correlation, struct ImageData* result);
#ifdef __cplusplus
};
#endif
#endif

