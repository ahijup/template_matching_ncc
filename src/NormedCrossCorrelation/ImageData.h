#pragma once

#include <stdlib.h>
#include "image_api.h"

#ifndef PI
#define PI 3.14159265358979323846
#endif



struct ImageData
{
	int w, h, ch;
	unsigned char* data;
};

static void* alloc_data(int size)
{
	return _aligned_malloc(size, 32);
}

static inline struct ImageData read_image(const char* filename)
{
	struct ImageData img;
	img.data = (unsigned char*)read_image_ascii(filename, &img.w, &img.h, &img.ch, alloc_data);
	return img;
}

static inline void write_image(const char* filename, struct ImageData* img)
{
	write_image_ascii(filename, img->w, img->h, img->ch, img->data);
}

static inline char is_image_empty(struct ImageData* img)
{
	return img->data == 0;
}

static inline void free_image(struct ImageData* img)
{
	if (!is_image_empty(img))
		_aligned_free(img->data);
}

// find the maximum value and its index
static void find_max_value_index(struct ImageData* result, int* x, int* y, float* max_value)
{
	*max_value = -1;
	*x = 0;
	*y = 0;
	for (int j = 0; j < result->h; j++)
	{
		float* p = (float*)(result->data + j * result->w * 4);
		for (int i = 0; i < result->w; i++)
		{
			if (p[i] > *max_value)
			{
				*max_value = p[i];
				*x = i;
				*y = j;
			}
		}
	}
}
