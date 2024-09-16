#include "ncc.h"
#include <omp.h>

// original version
void normed_cross_correlation(struct ImageData* pattern, struct ImageData* target, struct ImageData* result)
{
	// result size
	result->w = target->w - pattern->w + 1;
	result->h = target->h - pattern->h + 1;
	result->ch = 1;

	// allocate memory for result
	result->data = (unsigned char*)alloc_data(result->w * result->h * 4);

	

	unsigned char** t = (unsigned char**)malloc(sizeof(unsigned char*) * target->h);
	for (int i = 0; i < target->h; i++)
	{
		t[i] = (unsigned char*)(target->data + i * target->w);
	}
	unsigned char** p = (unsigned char**)malloc(sizeof(unsigned char*) * pattern->h);
	for (int i = 0; i < pattern->h; i++)
	{
		p[i] = (unsigned char*)(pattern->data + i * pattern->w);
	}
	float** r = (float**)malloc(sizeof(float*) * result->h);
	for (int i = 0; i < result->h; i++)
	{
		r[i] = (float*)(result->data + i * result->w * 4);
	}


	// pre-calculating mean and standard deviation of pattern 
	float* pattern_diff_map = (float*)alloc_data(sizeof(float) * pattern->w * pattern->h);
	float pattern_mean = 0;
	float pattern_std = 0;
	for (int j = 0; j < pattern->h; j++)
	{
		for (int i = 0; i < pattern->w; i++)
		{
			pattern_mean += p[j][i];
		}
	}

	pattern_mean /= (pattern->w * pattern->h);
	for (int j = 0; j < pattern->h; j++)
	{
		for (int i = 0; i < pattern->w; i++)
		{
			pattern_diff_map[j * pattern->w + i] = p[j][i] - pattern_mean;
			pattern_std += pattern_diff_map[j * pattern->w + i] * pattern_diff_map[j * pattern->w + i];
		}
	}

	float** pattern_diff_map_2d = (float**)malloc(sizeof(float*) * pattern->h);
	for (int i = 0; i < pattern->h; i++)
	{
		pattern_diff_map_2d[i] = pattern_diff_map + i * pattern->w;
	}

	int y;
	// calculate normed cross correlation
#pragma omp parallel for shared(t, p, r, pattern_diff_map)
	for (y = 0; y < result->h; y++)
	{
		float* rr = r[y];
		for (int x = 0; x < result->w; x++)
		{
			float sum_target = 0;
			for (int j = 0; j < pattern->h; j++)
			{
				for (int i = 0; i < pattern->w; i++)
				{
					sum_target += (t[y + j][x + i]);
				}
			}

			float target_mean = sum_target / (pattern->w * pattern->h);


			float pattern_std = 0;
			float target_std = 0;

			float sum_crossed = 0;

			for (int j = 0; j < pattern->h; ++j)
			{
				float* p_d_ptr = pattern_diff_map_2d[j];
				for (int i = 0; i < pattern->w; ++i)
				{
					float p_d = p_d_ptr[i];
					float t_d = t[y + j][x + i] - target_mean;

					pattern_std += (p_d) * (p_d);
					target_std += (t_d) * (t_d);
					sum_crossed += (p_d) * (t_d);
				}
			}

			rr[x] = sum_crossed / sqrtf(pattern_std * target_std);
		}
	}

	free(t);
	free(p);
	free(r);
	_aligned_free(pattern_diff_map);
	free(pattern_diff_map_2d);
}
