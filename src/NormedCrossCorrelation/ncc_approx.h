#ifndef __NCC_APPROX_H__
#define __NCC_APPROX_H__

#include "ImageData.h"

#if __cplusplus
extern "C" {
#endif
	void ncc_approx_offline(struct ImageData* pattern);
	void ncc_approx_online(struct ImageData* target, struct ImageData* result);
#if __cplusplus
}
#endif
#endif
