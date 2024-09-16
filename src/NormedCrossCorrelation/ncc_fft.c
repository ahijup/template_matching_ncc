#include "ncc_fft.h"
#include <stdio.h>
#include <math.h>
#include <windows.h>

typedef struct {
	float real;
	float imag;
} Complex;

int get_opt_fft_size(int n)
{
	int m = 1;
	while (m < n)
	{
		m *= 2;
	}
	return m;
}

// from: https://github.com/ttang10/FFT_2D_CONVOLUTION/blob/master/source%20codes/fft.c

static Complex ctmp;
#define C_SWAP(a,b) {ctmp=(a);(a)=(b);(b)=ctmp;}

void c_fft1d(Complex* r, int      n, int      isign)
{
	int     m, i, i1, j, k, i2, l, l1, l2;
	float   c1, c2, z;
	Complex t, u;

	if (isign == 0) return;

	/* Do the bit reversal */
	i2 = n >> 1;
	j = 0;
	for (i = 0; i < n - 1; i++) {
		if (i < j)
			C_SWAP(r[i], r[j]);
		k = i2;
		while (k <= j) {
			j -= k;
			k >>= 1;
		}
		j += k;
	}

	/* m = (int) log2((double)n); */
	for (i = n, m = 0; i > 1; m++, i /= 2);

	/* Compute the FFT */
	c1 = -1.0;
	c2 = 0.0;
	l2 = 1;
	for (l = 0; l < m; l++) {
		l1 = l2;
		l2 <<= 1;
		u.real = 1.0;
		u.imag = 0.0;
		for (j = 0; j < l1; j++) {
			for (i = j; i < n; i += l2) {
				i1 = i + l1;

				/* t = u * r[i1] */
				t.real = u.real * r[i1].real - u.imag * r[i1].imag;
				t.imag = u.real * r[i1].imag + u.imag * r[i1].real;

				/* r[i1] = r[i] - t */
				r[i1].real = r[i].real - t.real;
				r[i1].imag = r[i].imag - t.imag;

				/* r[i] = r[i] + t */
				r[i].real += t.real;
				r[i].imag += t.imag;
			}
			z = u.real * c1 - u.imag * c2;

			u.imag = u.real * c2 + u.imag * c1;
			u.real = z;
		}
		c2 = sqrtf((1.0f - c1) / 2.0f);
		if (isign == -1) /* FWD FFT */
			c2 = -c2;
		c1 = sqrtf((1.0f + c1) / 2.0f);
	}

	/* Scaling for inverse transform */
	if (isign == 1) {       /* IFFT*/
		for (i = 0; i < n; i++) {
			r[i].real /= n;
			r[i].imag /= n;
		}
	}
}

// Function to compute 1D FFT (Cooley-Tukey)
void fft1D(Complex* X, int N) {
	c_fft1d(X, N, -1);
}

void transpose_complex_matrix(struct ImageData* src, struct ImageData* dst)
{
	int i, j;
	Complex** s2d = (Complex**)malloc(sizeof(Complex*) * src->h);
	for (i = 0; i < src->h; i++)
	{
		s2d[i] = (Complex*)(src->data + i * src->w * sizeof(Complex));
	}

	Complex** d2d = (Complex**)malloc(sizeof(Complex*) * dst->h);

	for (i = 0; i < dst->h; i++)
	{
		d2d[i] = (Complex*)(dst->data + i * dst->w * sizeof(Complex));

	}

//#pragma omp parallel for shared(s2d, d2d)
	for (i = 0; i < dst->h; i++)
	{
		Complex* drow = d2d[i];
		for (j = 0; j < dst->w; j++)
		{
			drow[j] = s2d[j][i];
		}
	}

	free(s2d);
	free(d2d);
}

// fft-based version
// fft forward function
void fft2d(struct ImageData* src, struct ImageData* dst)
{
	// allocate memory for dst
	dst->w = get_opt_fft_size(src->w);
	dst->h = get_opt_fft_size(src->h);
	dst->ch = src->ch;
	dst->data = (unsigned char*)alloc_data(dst->w * dst->h * sizeof(Complex));

	unsigned char** s2d = (unsigned char**)malloc(sizeof(unsigned char*) * src->h);
	for (int i = 0; i < src->h; i++)
	{
		s2d[i] = (unsigned char*)(src->data + i * src->w);
	}

	Complex** d2d = (Complex**)malloc(sizeof(Complex*) * dst->h);

	for (int i = 0; i < dst->h; i++)
	{
		d2d[i] = (Complex*)(dst->data + i * dst->w * sizeof(Complex));
	}
	
	int i;
#pragma omp parallel for shared(s2d, d2d)
		// fft forward without any library
	for (i = 0; i < src->h; i++)
	{
		Complex* row = d2d[i];
		memset(row, 0x00, sizeof(Complex) * dst->w);
		unsigned char* srow = s2d[i];
		for (int j = 0; j < src->w; j++)
		{
			row[j].real = srow[j];
			row[j].imag = 0;
		}
		fft1D(row, dst->w);
	}

	for (int i = src->h; i < dst->h; i++)
	{
		Complex* row = d2d[i];
		memset(row, 0x00, sizeof(Complex) * dst->w);
	}

	struct ImageData dst_transposed;
	dst_transposed.w = dst->h;
	dst_transposed.h = dst->w;
	dst_transposed.ch = 1;
	dst_transposed.data = (unsigned char*)alloc_data(dst->w * dst->h * sizeof(Complex));
	transpose_complex_matrix(dst, &dst_transposed);

	Complex** d2d_transposed = (Complex**)malloc(sizeof(Complex*) * dst_transposed.h);
	for (i = 0; i < dst_transposed.h; i++)
	{
		d2d_transposed[i] = (Complex*)(dst_transposed.data + i * dst_transposed.w * sizeof(Complex));
	}

#pragma omp parallel for shared(d2d_transposed)
	for (i = 0; i < dst_transposed.h; i++)
	{
		Complex* row = d2d_transposed[i];
		fft1D(row, dst_transposed.w);
	}

	// transpose back
	transpose_complex_matrix(&dst_transposed, dst);

	// fft forward without any library
	free(s2d);
	free(d2d);
	free(d2d_transposed);
	free_image(&dst_transposed);
}


// Function to compute 1D Inverse FFT
void ifft1D(Complex* X, int N) {
	c_fft1d(X, N, 1);
}

void ifft2d(struct ImageData* src, struct ImageData* dst)
{
	// allocate memory for dst
	dst->w = get_opt_fft_size(src->w);
	dst->h = get_opt_fft_size(src->h);
	dst->ch = src->ch;
	dst->data = (unsigned char*)alloc_data(dst->w * dst->h * sizeof(Complex));

	Complex** s2d = (Complex**)malloc(sizeof(Complex*) * src->h);
	for (int i = 0; i < src->h; i++)
	{
		s2d[i] = (Complex*)(src->data + i * src->w * sizeof(Complex));
	}

	Complex** d2d = (Complex**)malloc(sizeof(Complex*) * dst->h);

	for (int i = 0; i < dst->h; i++)
	{
		d2d[i] = (Complex*)(dst->data + i * dst->w * sizeof(Complex));
	}
	
	int i;
#pragma omp parallel for shared(s2d, d2d)
			// fft forward without any library
	for (i = 0; i < src->h; i++)
	{
		Complex* row = d2d[i];
		Complex* srow = s2d[i];
		memcpy(row, srow, sizeof(Complex) * src->w);
		ifft1D(row, src->w);
	}

	for (i = src->h; i < dst->h; i++)
	{
		Complex* row = d2d[i];
		memset(row, 0x00, sizeof(Complex) * dst->w);
	}

	struct ImageData dst_transposed;
	dst_transposed.w = dst->h;
	dst_transposed.h = dst->w;
	dst_transposed.ch = 1;
	dst_transposed.data = (unsigned char*)alloc_data(dst->w * dst->h * sizeof(Complex));
	transpose_complex_matrix(dst, &dst_transposed);

	Complex** d2d_transposed = (Complex**)malloc(sizeof(Complex*) * dst_transposed.h);
	for (i = 0; i < dst_transposed.h; i++)
	{
		d2d_transposed[i] = (Complex*)(dst_transposed.data + i * dst_transposed.w * sizeof(Complex));
	}

#pragma omp parallel for shared(d2d_transposed)
	for (i = 0; i < dst_transposed.h; i++)
	{
		Complex* row = d2d_transposed[i];
		ifft1D(row, dst_transposed.w);
	}

	// transpose back
	transpose_complex_matrix(&dst_transposed, dst);

	free(s2d);
	free(d2d);
	free_image(&dst_transposed);
	free(d2d_transposed);
}




// test fft forward and inverse
void test_fft2d(struct ImageData* src)
{
	struct ImageData idst_u8;

	idst_u8.w = src->w;
	idst_u8.h = src->h;
	idst_u8.ch = 1;
	idst_u8.data = (unsigned char*)alloc_data(idst_u8.w * idst_u8.h * idst_u8.ch);

	LARGE_INTEGER frequency; // Stores the frequency of the high-resolution performance counter
	LARGE_INTEGER tic, toc; // Variables to store the start and end times

	// Get the frequency of the performance counter
	QueryPerformanceFrequency(&frequency);

	// Get the starting time
	QueryPerformanceCounter(&tic);

	struct ImageData dst;
	fft2d(src, &dst);

	struct ImageData idst;
	ifft2d(&dst, &idst);

	// Get the ending time
	QueryPerformanceCounter(&toc);

	// Calculate the elapsed time
	printf("Elapsed time: %f ms\n", (1000.0f) * (toc.QuadPart - tic.QuadPart) / frequency.QuadPart);
	

	// get the real part of the result
	for (int i = 0; i < idst_u8.h; ++i)
	{
		Complex* row = (Complex*)(idst.data + i * idst.w * sizeof(Complex));
		unsigned char* row_u8 = (unsigned char*)(idst_u8.data + i * idst_u8.w);
		for (int j = 0; j < idst_u8.w; ++j)
		{
			if (row[j].real > 255)
			{
				row[j].real = 255;
			}
			if (row[j].real < 0)
			{
				row[j].real = 0;
			}
			row_u8[j] = (unsigned char)row[j].real;
		}
	}
	write_image("ifft.bmp", &idst_u8);
	free_image(&idst_u8);
	free_image(&dst);
	free_image(&idst);
}

void complex_mul_fft(struct ImageData* pattern_fft, struct ImageData* target_fft,
	int w, int h,
	struct ImageData* result)
{
	// cross correlation in frequency domain
	// get the pointer of the 2D array
	Complex** p2d = (Complex**)malloc(sizeof(Complex*) * pattern_fft->h);
	for (int i = 0; i < pattern_fft->h; i++)
	{
		p2d[i] = (Complex*)(pattern_fft->data + i * pattern_fft->w * sizeof(Complex));
	}

	Complex** t2d = (Complex**)malloc(sizeof(Complex*) * target_fft->h);
	for (int i = 0; i < target_fft->h; i++)
	{
		t2d[i] = (Complex*)(target_fft->data + i * target_fft->w * sizeof(Complex));
	}

	// allocate memory for result
	result->ch = 1;
	result->w = w;
	result->h = h;
	result->data = (unsigned char*)alloc_data(w * h * sizeof(Complex));

	Complex** r2d = (Complex**)malloc(sizeof(Complex*) * result->h);
	for (int i = 0; i < result->h; i++)
	{
		r2d[i] = (Complex*)(result->data + i * result->w * sizeof(Complex));
	}
	int i;
#pragma omp parallel for shared(p2d, t2d)
	for (i = 0; i < h; i++)
	{
		Complex* t_row = t2d[i];
		Complex* r_row = r2d[i];
		Complex* p_row = p2d[i];

		for (int j = 0; j < w; ++j)
		{
			// complex multiplication
			//r_row[j].real = p_row[j].real * t_row[j].real - p_row[j].imag * t_row[j].imag;
			//r_row[j].imag = p_row[j].real * t_row[j].imag + p_row[j].imag * t_row[j].real;

			// complex multiplication with conjugate
			r_row[j].real = p_row[j].real * t_row[j].real + p_row[j].imag * t_row[j].imag;
			r_row[j].imag = p_row[j].real * t_row[j].imag - p_row[j].imag * t_row[j].real;
		}
	}


	free(p2d);
	free(t2d);
	free(r2d);
}

void get_real_part_fft(struct ImageData* src, struct ImageData* dst)
{
	dst->ch = 1;
	dst->data = (unsigned char*)alloc_data(dst->w * dst->h * 4);

	Complex** s2d = (Complex**)malloc(sizeof(Complex*) * src->h);
	for (int i = 0; i < src->h; i++)
	{
		s2d[i] = (Complex*)(src->data + i * src->w * sizeof(Complex));
	}

	float** d2d = (float**)malloc(sizeof(float*) * dst->h);

	for (int i = 0; i < dst->h; i++)
	{
		d2d[i] = (float*)(dst->data + i * dst->w * 4);
	}

	int i;
#pragma omp parallel for shared(s2d, d2d)
	for (i = 0; i < dst->h; i++)
	{
		Complex* row = s2d[i];
		float* drow = d2d[i];
		for (int j = 0; j < dst->w; j++)
		{
			drow[j] = row[j].real;
		}
	}

	free(s2d);
	free(d2d);

}

void abs_fft(struct ImageData* src, struct ImageData* dst)
{
	dst->ch = 1;
	dst->data = (unsigned char*)alloc_data(dst->w * dst->h * 4);

	Complex** s2d = (Complex**)malloc(sizeof(Complex*) * src->h);
	for (int i = 0; i < src->h; i++)
	{
		s2d[i] = (Complex*)(src->data + i * src->w * sizeof(Complex));
	}

	float** d2d = (float**)malloc(sizeof(float*) * dst->h);

	for (int i = 0; i < dst->h; i++)
	{
		d2d[i] = (float*)(dst->data + i * dst->w * 4);
	}

	// #pragma omp parallel for shared(s2d, d2d)
	for (int i = 0; i < dst->h; i++)
	{
		Complex* row = s2d[i];
		float* drow = d2d[i];
		for (int j = 0; j < dst->w; j++)
		{
			drow[j] = sqrtf(row[j].real * row[j].real + row[j].imag * row[j].imag);
		}
	}

	free(s2d);
	free(d2d);
}


void cross_correlation_fft(struct ImageData* pattern, struct ImageData* target, struct ImageData* result)
{
	// result size
	result->w = target->w - pattern->w + 1;
	result->h = target->h - pattern->h + 1;
	result->ch = 1;

	// allocate memory for result
	result->data = (unsigned char*)alloc_data(result->w * result->h * 4);

	struct ImageData target_fft;
	fft2d(target, &target_fft);

	// zero padding for pattern image to the size of target image
	struct ImageData pattern_padded;
	pattern_padded.w = target->w;
	pattern_padded.h = target->h;
	pattern_padded.ch = 1;
	pattern_padded.data = (unsigned char*)alloc_data(pattern_padded.w * pattern_padded.h);

	memset(pattern_padded.data, 0x00, pattern_padded.w * pattern_padded.h);
	for (int i = 0; i < pattern->h; i++)
	{
		memcpy(pattern_padded.data + i * pattern_padded.w, pattern->data + i * pattern->w, pattern->w);
	}

	struct ImageData pattern_fft;
	fft2d(&pattern_padded, &pattern_fft);

	struct ImageData result_fft;
	complex_mul_fft(&pattern_fft, &target_fft, pattern_fft.w, pattern_fft.h, &result_fft);


	struct ImageData idst;
	ifft2d(&result_fft, &idst);

	result->w = target->w - pattern->w + 1;
	result->h = target->h - pattern->h + 1;

	// get the real part of the result
	get_real_part_fft(&idst, result);

	free_image(&target_fft);
	free_image(&pattern_padded);
	free_image(&pattern_fft);
	free_image(&result_fft);
	free_image(&idst);
}

void image_sum_table(struct  ImageData* src, struct ImageData* sum, struct ImageData* sqSum)
{
	sum->w = src->w + 1;
	sum->h = src->h + 1;
	sum->ch = 1;
	sum->data = (unsigned char*)alloc_data(sum->w * sum->h * 4);

	sqSum->w = src->w + 1;
	sqSum->h = src->h + 1;
	sqSum->ch = 1;
	sqSum->data = (unsigned char*)alloc_data(sqSum->w * sqSum->h * 4);

	unsigned char** s = (unsigned char**)malloc(sizeof(unsigned char*) * src->h);
	for (int i = 0; i < src->h; i++)
	{
		s[i] = (unsigned char*)(src->data + i * src->w);
	}
	float** sum2d = (float**)malloc(sizeof(float*) * sum->h);
	for (int i = 0; i < sum->h; i++)
	{
		sum2d[i] = (float*)(sum->data + i * sum->w * 4);
	}
	float** sqSum2d = (float**)malloc(sizeof(float*) * sqSum->h);
	for (int i = 0; i < sqSum->h; i++)
	{
		sqSum2d[i] = (float*)(sqSum->data + i * sqSum->w * 4);
	}

	// calculate integral image
	for (int j = 0; j < sum->h; j++)
	{
		float* sum_row = sum2d[j];
		float* sqSum_row = sqSum2d[j];
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

void image_mean(struct ImageData* src, float* mean, float* sdv)
{
	unsigned char** s = (unsigned char**)malloc(sizeof(unsigned char*) * src->h);
	for (int i = 0; i < src->h; i++)
	{
		s[i] = (unsigned char*)(src->data + i * src->w);
	}

	float sum = 0;
	float sum_sq = 0;
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
	*sdv = sqrtf(sum_sq / (src->w * src->h) - (*mean) * (*mean));
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

	float** t = (float**)malloc(sizeof(float*) * integral.h);
	for (int i = 0; i < integral.h; i++)
	{
		t[i] = (float*)(integral.data + i * integral.w * 4);
	}

	float** t_sq = (float**)malloc(sizeof(float*) * integral_sq.h);
	for (int i = 0; i < integral_sq.h; i++)
	{
		t_sq[i] = (float*)(integral_sq.data + i * integral_sq.w * 4);
	}

#define GET_SUM(x, y, w, h) (t[(y + h)][(x + w)] + t[(y)][(x)] - t[(y + h)][(x)] - t[(y)][(x + w)])
#define GET_SQ_SUM(x, y, w, h) (t_sq[(y + h)][(x + w)] + t_sq[(y)][(x)] - t_sq[(y + h)][(x)] - t_sq[(y)][(x + w)])

	printf("Integral image\n");
	for (int j = 0; j < integral.h; j++)
	{
		float* p = (float*)(integral.data + j * integral.w * 4);
		for (int i = 0; i < integral.w; i++)
		{
			printf("%f ", p[i]);
		}
		printf("\n");
	}

	printf("Integral image square\n");
	for (int j = 0; j < integral_sq.h; j++)
	{
		float* p = (float*)(integral_sq.data + j * integral_sq.w * 4);
		for (int i = 0; i < integral_sq.w; i++)
		{
			printf("%f ", p[i]);
		}
		printf("\n");
	}


	printf("sum1: %f\n", GET_SUM(1, 0, 3, 3));
	printf("sum2: %f\n", GET_SUM(0, 0, 5, 5));

	printf("sq sum1: %f\n", GET_SQ_SUM(1, 0, 3, 3));
	printf("sq sum2: %f\n", GET_SQ_SUM(0, 0, 5, 5));


	free(t);
	free(t_sq);

	free_image(&img);
	free_image(&integral);
	free_image(&integral_sq);
#undef GET_SUM
#undef GET_SQ_SUM
}

#ifndef MAX
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#endif

void normed_cross_correlation_fft(struct ImageData* pattern, struct ImageData* target, struct ImageData* result)
{
	cross_correlation_fft(pattern, target, result);

	struct ImageData target_integral, target_integral_sq;
	image_sum_table(target, &target_integral, &target_integral_sq);

	float pattern_mean, pattern_sdv;
	image_mean(pattern, &pattern_mean, &pattern_sdv);
	float pattern_sum = pattern_mean * pattern->w * pattern->h;
	float pattern_norm = pattern_sdv / sqrtf(1.0f / (pattern->w * pattern->h));


	float** targetSum2d = (float**)malloc(sizeof(float*) * target_integral.h);
	float** targetSqSum2d = (float**)malloc(sizeof(float*) * target_integral.h);
	for (int i = 0; i < target_integral.h; i++)
	{
		targetSum2d[i] = (float*)(target_integral.data + i * target_integral.w * 4);
		targetSqSum2d[i] = (float*)(target_integral_sq.data + i * target_integral_sq.w * 4);
	}
#define GET_SUM(x, y, w, h) (targetSum2d[(y + h)][(x + w)] + targetSum2d[(y)][(x)] - targetSum2d[(y + h)][(x)] - targetSum2d[(y)][(x + w)])
#define GET_SQ_SUM(x, y, w, h) (targetSqSum2d[(y + h)][(x + w)] + targetSqSum2d[(y)][(x)] - targetSqSum2d[(y + h)][(x)] - targetSqSum2d[(y)][(x + w)])


	float** r2d = (float**)malloc(sizeof(float*) * result->h);
	for (int i = 0; i < result->h; i++)
	{
		r2d[i] = (float*)(result->data + i * result->w * 4);
	}

	// #pragma omp parallel for shared(r2d, targetSum2d, targetSqSum2d)
	for (int j = 0; j < result->h; j++)
	{
		float* rrow = r2d[j];
		for (int i = 0; i < result->w; i++)
		{
			float num = rrow[i];
			float sum = GET_SUM(i, j, pattern->w, pattern->h);
			float sqSum = GET_SQ_SUM(i, j, pattern->w, pattern->h);
			num -= sum * pattern_mean;

			float wndMean = sum * sum / (pattern->w * pattern->h);
			float wndSum2 = sqSum;

			float diff2 = MAX(wndSum2 - wndMean, 0.0f);
			float den = sqrtf(diff2) * pattern_norm;
			rrow[i] = num / den;
		}
	}

	free(targetSum2d);
	free(targetSqSum2d);
	free(r2d);
	free_image(&target_integral);
	free_image(&target_integral_sq);
#undef GET_SUM
#undef GET_SQ_SUM
}


void test_cross_correlation_fft(struct ImageData* pattern, struct ImageData* target)
{
	struct ImageData result;
	// get tick count in C
	LARGE_INTEGER frequency; // Stores the frequency of the high-resolution performance counter
	LARGE_INTEGER tic, toc; // Variables to store the start and end times

	// Get the frequency of the performance counter
	QueryPerformanceFrequency(&frequency);

	// Get the starting time
	QueryPerformanceCounter(&tic);

	cross_correlation_fft(pattern, target, &result);
	
	// Get the ending time
	QueryPerformanceCounter(&toc);

	printf("Elapsed time[cross correlation fft]: %f ms\n", (1000.0f) * (toc.QuadPart - tic.QuadPart) / frequency.QuadPart);

	int x, y;
	float max_value;
	find_max_value_index(&result, &x, &y, &max_value);
	printf("x= %d, y= %d, max_value= %f\n", x, y, max_value);
	free_image(&result);
}


