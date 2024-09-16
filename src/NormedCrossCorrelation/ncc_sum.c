#include "ncc_sum.h"
#include <stdio.h>

// calculate integral image
void integral_image(struct ImageData* src, struct ImageData* dst)
{
	dst->w = src->w + 1;
	dst->h = src->h + 1;
	dst->ch = 1;
	dst->data = (unsigned char*)alloc_data(dst->w * dst->h * 4);

	unsigned char** s = (unsigned char**)malloc(sizeof(unsigned char*) * src->h);
	for (int i = 0; i < src->h; i++)
	{
		s[i] = (unsigned char*)(src->data + i * src->w);
	}
	float** d = (float**)malloc(sizeof(float*) * dst->h);
	for (int i = 0; i < dst->h; i++)
	{
		d[i] = (float*)(dst->data + i * dst->w * 4);
	}

	// calculate integral image
	for (int j = 0; j < dst->h; j++)
	{
		float* dd = d[j];
		for (int i = 0; i < dst->w; i++)
		{
			if (i == 0 || j == 0)
			{
				dd[i] = 0;
			}
			else
			{
				dd[i] = s[j - 1][i - 1] + d[j - 1][i] + d[j][i - 1] - d[j - 1][i - 1];
			}
		}
	}

	free(s);
	free(d);
}


void normed_cross_correlation_integral(struct ImageData* pattern, struct ImageData* target, struct ImageData* result)
{
	// result size
	result->w = target->w - pattern->w + 1;
	result->h = target->h - pattern->h + 1;
	result->ch = 1;

	// allocate memory for result
	result->data = (unsigned char*)alloc_data(result->w * result->h * 4);

	// calculate integral image
	struct ImageData target_integral;
	integral_image(target, &target_integral);

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

	float** t = (float**)malloc(sizeof(float*) * target_integral.h);
	for (int i = 0; i < target_integral.h; i++)
	{
		t[i] = (float*)(target_integral.data + 4 * i * target_integral.w);
	}
#define GET_SUM(x, y, w, h) (t[(y + h)][(x + w)] + t[(y)][(x)] - t[(y + h)][(x)] - t[(y)][(x + w)])

	unsigned char** tp = (unsigned char**)malloc(sizeof(unsigned char*) * target->h);
	for (int i = 0; i < target->h; i++)
	{
		tp[i] = (unsigned char*)(target->data + i * target->w);
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
#pragma omp parallel for shared(t, p, r, pattern_diff_map)
	for (y = 0; y < result->h; y++)
	{
		float* rr = r[y];
		for (int x = 0; x < result->w; x++)
		{
			float sum_target = GET_SUM(x, y, pattern->w, pattern->h);
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
					float t_d = tp[y + j][x + i] - target_mean;
					pattern_std += (p_d) * (p_d);
					target_std += (t_d) * (t_d);
					sum_crossed += (p_d) * (t_d);
				}
			}
			rr[x] = sum_crossed / sqrtf(pattern_std * target_std);
		}
	}

	free_image(&target_integral);
	free(p);
	free(r);
	_aligned_free(pattern_diff_map);
	free(t);
	free(tp);
#undef GET_SUM
}


void test_integral_image()
{
	struct ImageData img;
	img.data = (unsigned char*)alloc_data(5 * 5);
	img.w = 5;
	img.h = 5;
	img.ch = 1;

	for (int i = 0; i < 5 * 5; i++)
	{
		img.data[i] = i;
	}

	struct ImageData integral;
	integral_image(&img, &integral);

	float** t = (float**)malloc(sizeof(float*) * integral.h);
	for (int i = 0; i < integral.h; i++)
	{
		t[i] = (float*)(integral.data + i * integral.w * 4);
	}

#define GET_SUM(x, y, w, h) (t[(y + h)][(x + w)] + t[(y)][(x)] - t[(y + h)][(x)] - t[(y)][(x + w)])


	for (int j = 0; j < integral.h; j++)
	{
		float* p = (float*)(integral.data + j * integral.w * 4);
		for (int i = 0; i < integral.w; i++)
		{
			printf("%f ", p[i]);
		}
		printf("\n");
	}

	printf("sum1 : %f\n", GET_SUM(1, 0, 3, 3));
	printf("sum2 : %f\n", GET_SUM(0, 0, 5, 5));


	free_image(&img);
	free_image(&integral);
#undef GET_SUM
}
