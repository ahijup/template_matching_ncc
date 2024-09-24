// NormedCrossCorrelation.cpp: 定義應用程式的進入點。
//
#include <stdio.h>
#include <stdlib.h>
#include <Windows.h>

#include "ImageData.h"
#include "ncc.h"
#include "ncc_sum.h"
#include "ncc_fft.h"
#include "ncc_approx.h"
#include "ncc_simd.h"

void split_3channels(struct ImageData* src, struct ImageData* dst1, struct ImageData* dst2, struct ImageData* dst3)
{
	// create 3 images with the same size as src
	dst1->data = (unsigned char*)alloc_data(src->w * src->h * 1);
	dst2->data = (unsigned char*)alloc_data(src->w * src->h * 1);
	dst3->data = (unsigned char*)alloc_data(src->w * src->h * 1);

	dst1->w = dst2->w = dst3->w = src->w;
	dst1->h = dst2->h = dst3->h = src->h;
	dst1->ch = dst2->ch = dst3->ch = 1;

	// split 3 channels
	for (int i = 0; i < src->w * src->h; i++)
	{
		dst1->data[i] = src->data[i * 3 + 0];
		dst2->data[i] = src->data[i * 3 + 1];
		dst3->data[i] = src->data[i * 3 + 2];
	}
}


void print_max_value_neighbor(struct ImageData* result, int x, int y, int ext_size)
{
	float* src = (float*)(result->data);
#define GET_VALUE(x, y) (x >= 0 && x < result->w && y >= 0 && y < result->h) ? src[y * result->w + x] : 0

	
	printf("max value: %f\n", GET_VALUE(x, y));
	printf("neighbor values:\n");
	for (int j = y - ext_size; j <= y + ext_size; j++)
	{
		for (int i = x - ext_size; i <= x + ext_size; i++)
		{
			//if (i >= 0 && i < result->w && j >= 0 && j < result->h)
			{
				printf("%f ", GET_VALUE(i, j));
			}
		}
		printf("\n");
	}
#undef GET_VALUE
}

// subpixel using 2nd order taylor expansion
void get_subpixel_peak(struct ImageData* result, int x, int y, float* subpixel_x, float* subpixel_y)
{
	float* src = (float*)(result->data);
	float v0 = src[y * result->w + x];
	float vx1 = src[y * result->w + x + 1];
	float vx_1 = src[y * result->w + x - 1];
	float vy1 = src[(y + 1) * result->w + x];
	float vy_1 = src[(y - 1) * result->w + x];

	*subpixel_x = x + 0.5f * (vx1 - vx_1) / (2 * v0 - vx1 - vx_1);
	*subpixel_y = y + 0.5f * (vy1 - vy_1) / (2 * v0 - vy1 - vy_1);
}

int main()
{
	// test_integral_image();
	// test_image_sum();
	//

	struct ImageData pattern = read_image(TEST_IMAGE_DIR"/pattern.png");
	if (is_image_empty(&pattern))
	{
		printf("Failed to read pattern image.\n");
		return -1;
	}


	struct ImageData target = read_image(TEST_IMAGE_DIR"/source.png");
	if (is_image_empty(&target))
	{
		printf("Failed to read target image.\n");
		return -1;
	}

	

	struct ImageData pattern_r, pattern_g, pattern_b;
	split_3channels(&pattern, &pattern_r, &pattern_g, &pattern_b);
	struct ImageData target_r, target_g, target_b;
	split_3channels(&target, &target_r, &target_g, &target_b);


	// ncc_approx_offline(&pattern_r);
	// test_fft2d(&target_r);

	write_image("pattern_r.bmp", &pattern_r);
	write_image("target_r.bmp", &target_r);


	// test_cross_correlation_fft(&pattern_r, &target_r);
	// test_normed_cross_correlation_fft(&pattern_r, &target_r);

	LARGE_INTEGER frequency; // Stores the frequency of the high-resolution performance counter
	LARGE_INTEGER tic, toc; // Variables to store the start and end times
	int x, y;
	float max_value;
	float x_sub, y_sub;

	// Get the frequency of the performance counter
	QueryPerformanceFrequency(&frequency);

	// Get the starting time
	QueryPerformanceCounter(&tic);
	struct ImageData result;
	normed_cross_correlation(&pattern_r, &target_r, &result);
	// Get the ending time
	QueryPerformanceCounter(&toc);

	printf("Elapsed time[original]: %f ms\n", (1000.0f) * (toc.QuadPart - tic.QuadPart) / frequency.QuadPart);
	
	
	find_max_value_index(&result, &x, &y, &max_value);
	printf("Red channel: x= %d, y= %d, max_value= %f\n", x, y, max_value);
	print_max_value_neighbor(&result, x, y, 1);
	get_subpixel_peak(&result, x, y, &x_sub, &y_sub);
	printf("subpixel peak: x= %f, y= %f\n", x_sub, y_sub);
	printf("###############################################\n");

	QueryPerformanceCounter(&tic);
	struct ImageData result_i;
	normed_cross_correlation_integral(&pattern_r, &target_r, &result_i);
	QueryPerformanceCounter(&toc);


	printf("Elapsed time[integral]: %f ms\n", (1000.0f) * (toc.QuadPart - tic.QuadPart) / frequency.QuadPart);

	find_max_value_index(&result_i, &x, &y, &max_value);
	printf("Red channel: x= %d, y= %d, max_value= %f\n", x, y, max_value);
	print_max_value_neighbor(&result_i, x, y, 1);
	get_subpixel_peak(&result_i, x, y, &x_sub, &y_sub);
	printf("subpixel peak: x= %f, y= %f\n", x_sub, y_sub);
	printf("###############################################\n");

	QueryPerformanceCounter(&tic);
	struct ImageData result_fft;
	normed_cross_correlation_fft(&pattern_r, &target_r, &result_fft);
	QueryPerformanceCounter(&toc);

	printf("Elapsed time[fft]: %f ms\n", (1000.0f) * (toc.QuadPart - tic.QuadPart) / frequency.QuadPart);

	find_max_value_index(&result_fft, &x, &y, &max_value);
	printf("Red channel: x= %d, y= %d, max_value= %f\n", x, y, max_value);
	print_max_value_neighbor(&result_fft, x, y, 1);
	
	get_subpixel_peak(&result_fft, x, y, &x_sub, &y_sub);
	printf("subpixel peak: x= %f, y= %f\n", x_sub, y_sub);
	printf("###############################################\n");

	QueryPerformanceCounter(&tic);
	struct ImageData result_simd;
	normed_cross_correlation_simd(&pattern_r, &target_r, &result_simd);
	QueryPerformanceCounter(&toc);

	printf("Elapsed time[simd]: %f ms\n", (1000.0f) * (toc.QuadPart - tic.QuadPart) / frequency.QuadPart);

	find_max_value_index(&result_simd, &x, &y, &max_value);
	printf("Red channel: x= %d, y= %d, max_value= %f\n", x, y, max_value);
	print_max_value_neighbor(&result_simd, x, y, 1);

	get_subpixel_peak(&result_simd, x, y, &x_sub, &y_sub);
	printf("subpixel peak: x= %f, y= %f\n", x_sub, y_sub);
	printf("###############################################\n");

	free_image(&pattern);
	free_image(&target);
	
	free_image(&pattern_r);
	free_image(&pattern_g);
	free_image(&pattern_b);

	free_image(&target_r);
	free_image(&target_g);
	free_image(&target_b);

	/*free_image(&result);
	free_image(&result_i);*/
	free_image(&result_fft);
	free_image(&result_simd);
	system("pause");
	return 0;
}
