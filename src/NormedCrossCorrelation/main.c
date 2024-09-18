// NormedCrossCorrelation.cpp: 定義應用程式的進入點。
//
#include <stdio.h>
#include <stdlib.h>
#include <Windows.h>

#include "ImageData.h"
#include "ncc.h"
#include "ncc_sum.h"
#include "ncc_fft.h"

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

	// test_fft2d(&target_r);

	write_image("pattern_r.bmp", &pattern_r);
	write_image("target_r.bmp", &target_r);


	// test_cross_correlation_fft(&pattern_r, &target_r);
	// test_normed_cross_correlation_fft(&pattern_r, &target_r);

	LARGE_INTEGER frequency; // Stores the frequency of the high-resolution performance counter
	LARGE_INTEGER tic, toc; // Variables to store the start and end times
	int x, y;
	float max_value;

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

	QueryPerformanceCounter(&tic);
	struct ImageData result_i;
	normed_cross_correlation_integral(&pattern_r, &target_r, &result_i);
	QueryPerformanceCounter(&toc);


	printf("Elapsed time[integral]: %f ms\n", (1000.0f) * (toc.QuadPart - tic.QuadPart) / frequency.QuadPart);

	find_max_value_index(&result_i, &x, &y, &max_value);
	printf("Red channel: x= %d, y= %d, max_value= %f\n", x, y, max_value);


	QueryPerformanceCounter(&tic);
	struct ImageData result_fft;
	normed_cross_correlation_fft(&pattern_r, &target_r, &result_fft);
	QueryPerformanceCounter(&toc);

	printf("Elapsed time[fft]: %f ms\n", (1000.0f) * (toc.QuadPart - tic.QuadPart) / frequency.QuadPart);

	find_max_value_index(&result_fft, &x, &y, &max_value);
	printf("Red channel: x= %d, y= %d, max_value= %f\n", x, y, max_value);

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
	system("pause");
	return 0;
}
