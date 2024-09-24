#include "ncc_simd.h"
#include <math.h>

#ifndef MAX
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#endif

void image_sum_table(struct  ImageData* src, struct ImageData* sum, struct ImageData* sqSum)
{
	sum->w = src->w + 1;
	sum->h = src->h + 1;
	sum->ch = 1;
	sum->data = (unsigned char*)alloc_data(sum->w * sum->h * sizeof(double));

	sqSum->w = src->w + 1;
	sqSum->h = src->h + 1;
	sqSum->ch = 1;
	sqSum->data = (unsigned char*)alloc_data(sqSum->w * sqSum->h * sizeof(double));

	unsigned char** s = (unsigned char**)malloc(sizeof(unsigned char*) * src->h);
	for (int i = 0; i < src->h; i++)
	{
		s[i] = (unsigned char*)(src->data + i * src->w);
	}

	double** sum2d = (double**)malloc(sizeof(double*) * sum->h);
	for (int i = 0; i < sum->h; i++)
	{
		sum2d[i] = (double*)(sum->data + i * sum->w * sizeof(double));
	}
	double** sqSum2d = (double**)malloc(sizeof(double*) * sqSum->h);
	for (int i = 0; i < sqSum->h; i++)
	{
		sqSum2d[i] = (double*)(sqSum->data + i * sqSum->w * sizeof(double));
	}

	// calculate integral image
	for (int j = 0; j < sum->h; j++)
	{
		double* sum_row = sum2d[j];
		double* sqSum_row = sqSum2d[j];
		for (int i = 0; i < sum->w; i++)
		{
			if (i == 0 || j == 0)
			{
				sum_row[i] = 0;
				sqSum_row[i] = 0;
			}
			else
			{
				sum_row[i] = s[j - 1][i - 1] + sum2d[j - 1][i] + sum2d[j][i - 1] - sum2d[j - 1][i - 1];
				sqSum_row[i] = s[j - 1][i - 1] * s[j - 1][i - 1] + sqSum2d[j - 1][i] + sqSum2d[j][i - 1] - sqSum2d[j - 1][i - 1];
			}
		}
	}

	free(s);
	free(sum2d);
	free(sqSum2d);
}

void image_mean(struct ImageData* src, double* mean, double* sdv)
{
	unsigned char** s = (unsigned char**)malloc(sizeof(unsigned char*) * src->h);
	for (int i = 0; i < src->h; i++)
	{
		s[i] = (unsigned char*)(src->data + i * src->w);
	}

	double sum = 0;
	double sum_sq = 0;
	for (int j = 0; j < src->h; j++)
	{
		unsigned char* srow = s[j];
		for (int i = 0; i < src->w; i++)
		{
			sum += srow[i];
			sum_sq += srow[i] * srow[i];
		}
	}

	*mean = sum / (src->w * src->h);
	*sdv = sqrt(sum_sq / (src->w * src->h) - (*mean) * (*mean));
}

void test_image_sum()
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

	struct ImageData integral, integral_sq;
	image_sum_table(&img, &integral, &integral_sq);

	double** t = (double**)malloc(sizeof(double*) * integral.h);
	for (int i = 0; i < integral.h; i++)
	{
		t[i] = (double*)(integral.data + i * integral.w * sizeof(double));
	}

	double** t_sq = (double**)malloc(sizeof(double*) * integral_sq.h);
	for (int i = 0; i < integral_sq.h; i++)
	{
		t_sq[i] = (double*)(integral_sq.data + i * integral_sq.w * sizeof(double));
	}

#define GET_SUM(x, y, w, h) (t[(y + h)][(x + w)] + t[(y)][(x)] - t[(y + h)][(x)] - t[(y)][(x + w)])
#define GET_SQ_SUM(x, y, w, h) (t_sq[(y + h)][(x + w)] + t_sq[(y)][(x)] - t_sq[(y + h)][(x)] - t_sq[(y)][(x + w)])

	printf("Integral image\n");
	for (int j = 0; j < integral.h; j++)
	{
		double* p = (double*)(integral.data + j * integral.w * sizeof(double));
		for (int i = 0; i < integral.w; i++)
		{
			printf("%lf ", p[i]);
		}
		printf("\n");
	}

	printf("Integral image square\n");
	for (int j = 0; j < integral_sq.h; j++)
	{
		double* p = (double*)(integral_sq.data + j * integral_sq.w * sizeof(double));
		for (int i = 0; i < integral_sq.w; i++)
		{
			printf("%lf ", p[i]);
		}
		printf("\n");
	}


	printf("sum1: %lf\n", GET_SUM(1, 0, 3, 3));
	printf("sum2: %lf\n", GET_SUM(0, 0, 5, 5));

	printf("sq sum1: %lf\n", GET_SQ_SUM(1, 0, 3, 3));
	printf("sq sum2: %lf\n", GET_SQ_SUM(0, 0, 5, 5));


	free(t);
	free(t_sq);

	free_image(&img);
	free_image(&integral);
	free_image(&integral_sq);
#undef GET_SUM
#undef GET_SQ_SUM
}


void normalize_cross_correlation(struct ImageData* pattern, struct ImageData* target, struct ImageData* cross_correlation, struct ImageData* result)
{
	struct ImageData target_integral, target_integral_sq;
	image_sum_table(target, &target_integral, &target_integral_sq);

	double pattern_mean, pattern_sdv;
	image_mean(pattern, &pattern_mean, &pattern_sdv);
	double pattern_sum = pattern_mean * pattern->w * pattern->h;
	// sdv = sqrt(sum(x - mean)^2 / n)
	// sdv^2 = sum(x - mean)^2 / n
	// sdv^2 * n = sum(x - mean)^2
	// sum(x - mean)^2 = sqrt(sdv^2 * n)
	double pattern_norm_2 = pattern_sdv / sqrt(1.0 / (pattern->w * pattern->h));
	double pattern_norm = sqrt(pattern_sdv * pattern_sdv * pattern->w * pattern->h);


	double** targetSum2d = (double**)malloc(sizeof(double*) * target_integral.h);
	double** targetSqSum2d = (double**)malloc(sizeof(double*) * target_integral.h);
	for (int i = 0; i < target_integral.h; i++)
	{
		targetSum2d[i] = (double*)(target_integral.data + i * target_integral.w * sizeof(double));
		targetSqSum2d[i] = (double*)(target_integral_sq.data + i * target_integral_sq.w * sizeof(double));
	}
#define GET_SUM(x, y, w, h) (targetSum2d[(y + h)][(x + w)] + targetSum2d[(y)][(x)] - targetSum2d[(y + h)][(x)] - targetSum2d[(y)][(x + w)])
#define GET_SQ_SUM(x, y, w, h) (targetSqSum2d[(y + h)][(x + w)] + targetSqSum2d[(y)][(x)] - targetSqSum2d[(y + h)][(x)] - targetSqSum2d[(y)][(x + w)])


	result->w = cross_correlation->w;
	result->h = cross_correlation->h;
	result->ch = 1;
	result->data = (unsigned char*)alloc_data(result->w * result->h * sizeof(float));

	float** r2d = (float**)malloc(sizeof(float*) * result->h);
	for (int i = 0; i < result->h; i++)
	{
		r2d[i] = (float*)(result->data + i * result->w * 4);
	}

	float** resultcc_2d = (float**)malloc(sizeof(float*) * cross_correlation->h);
	for (int i = 0; i < cross_correlation->h; i++)
	{
		resultcc_2d[i] = (float*)(cross_correlation->data + i * cross_correlation->w * sizeof(float));
	}


	if (fabs(pattern_norm) > 1e-8f)
	{
		int j;
#pragma omp parallel for shared(r2d, targetSum2d, targetSqSum2d) private(j)
		for (j = 0; j < result->h; j++)
		{
			float* rrow = r2d[j];
			float* ccrow = resultcc_2d[j];

			for (int i = 0; i < result->w; i++)
			{
				float num = ccrow[i];
				double sum = GET_SUM(i, j, pattern->w, pattern->h);
				double sqSum = GET_SQ_SUM(i, j, pattern->w, pattern->h);
				num -= sum * pattern_mean;

				double wndMean = sum * sum / (pattern->w * pattern->h);
				double wndSum2 = sqSum;

				double diff2 = MAX(wndSum2 - wndMean, 0.0);

				if (fabs(diff2) < 1e-8)
				{
					rrow[i] = 0;
					continue;
				}
				double den = (sqrt(diff2) * pattern_norm);
				double num2 = num / den;

				rrow[i] = (float)num2;
			}
		}
	}
	else
	{
		// if pattern_norm is too small, return
		// to avoid division by zero
		// set the result to zero
		for (int j = 0; j < result->h; j++)
		{
			float* rrow = r2d[j];
			for (int i = 0; i < result->w; i++)
			{
				rrow[i] = 0;
			}
		}
	}

	free(targetSum2d);
	free(targetSqSum2d);
	free(r2d);
	free(resultcc_2d);
	free_image(&target_integral);
	free_image(&target_integral_sq);
	
#undef GET_SUM
#undef GET_SQ_SUM
}


