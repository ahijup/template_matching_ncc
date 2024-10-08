#include "ncc_fft.h"
#include <stdio.h>
#include <math.h>
#include <windows.h>
#include "ncc_common.h"

//#define FFT_DBL_TYPE double
//#define SQRT(a) sqrt((a))

 #define FFT_DBL_TYPE float
 #define SQRT(a) sqrtf((a))

typedef struct {
	FFT_DBL_TYPE real;
	FFT_DBL_TYPE imag;
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
void c_fft1d(Complex* r, int      n, int      isign)
{
	int     m, i, i1, j, k, i2, l, l1, l2;
	FFT_DBL_TYPE   c1, c2, z;
	Complex t, u;

	Complex ctmp;
#define C_SWAP(a,b) {ctmp=(a);(a)=(b);(b)=ctmp;}
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
		c2 = SQRT((1.0 - c1) / 2.0);
		if (isign == -1) /* FWD FFT */
			c2 = -c2;
		c1 = SQRT((1.0 + c1) / 2.0);
	}

	/* Scaling for inverse transform */
	if (isign == 1) {       /* IFFT*/
		for (i = 0; i < n; i++) {
			r[i].real /= n;
			r[i].imag /= n;
		}
	}
#undef C_SWAP
}

// Function to compute 1D FFT (Cooley-Tukey)
void fft1D(Complex* X, int N) {
	c_fft1d(X, N, -1);
}

void transpose_complex_matrix(struct ImageData* src, struct ImageData* dst)
{
	int i;
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

#pragma omp parallel for shared(s2d, d2d)
	for (i = 0; i < dst->h; i++)
	{
		Complex* drow = d2d[i];
		for (int j = 0; j < dst->w; j++)
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
	
	int i = 0;
#pragma omp parallel for shared(s2d, d2d) firstprivate(i)
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

#pragma omp parallel for shared(d2d_transposed) private(i)
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
	
	int i = 0;
#pragma omp parallel for shared(s2d, d2d) firstprivate(i)
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

#pragma omp parallel for shared(d2d_transposed) private(i)
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
// #pragma omp parallel for shared(p2d, t2d)
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
	dst->data = (unsigned char*)alloc_data(dst->w * dst->h * sizeof(FFT_DBL_TYPE));

	Complex** s2d = (Complex**)malloc(sizeof(Complex*) * src->h);
	for (int i = 0; i < src->h; i++)
	{
		s2d[i] = (Complex*)(src->data + i * src->w * sizeof(Complex));
	}

	FFT_DBL_TYPE** d2d = (FFT_DBL_TYPE**)malloc(sizeof(FFT_DBL_TYPE*) * dst->h);

	for (int i = 0; i < dst->h; i++)
	{
		d2d[i] = (FFT_DBL_TYPE*)(dst->data + i * dst->w * sizeof(FFT_DBL_TYPE));
	}

	int i;
// #pragma omp parallel for shared(s2d, d2d)
	for (i = 0; i < dst->h; i++)
	{
		Complex* row = s2d[i];
		FFT_DBL_TYPE* drow = d2d[i];
		for (int j = 0; j < dst->w; j++)
		{
			drow[j] = row[j].real;
		}
	}

	free(s2d);
	free(d2d);

}

void fabs_fft(struct ImageData* src, struct ImageData* dst)
{
	dst->ch = 1;
	dst->data = (unsigned char*)alloc_data(dst->w * dst->h * sizeof(FFT_DBL_TYPE));

	Complex** s2d = (Complex**)malloc(sizeof(Complex*) * src->h);
	for (int i = 0; i < src->h; i++)
	{
		s2d[i] = (Complex*)(src->data + i * src->w * sizeof(Complex));
	}

	FFT_DBL_TYPE** d2d = (FFT_DBL_TYPE**)malloc(sizeof(FFT_DBL_TYPE*) * dst->h);

	for (int i = 0; i < dst->h; i++)
	{
		d2d[i] = (FFT_DBL_TYPE*)(dst->data + i * dst->w * sizeof(FFT_DBL_TYPE));
	}

	// #pragma omp parallel for shared(s2d, d2d)
	for (int i = 0; i < dst->h; i++)
	{
		Complex* row = s2d[i];
		double* drow = d2d[i];
		for (int j = 0; j < dst->w; j++)
		{
			drow[j] = sqrt(row[j].real * row[j].real + row[j].imag * row[j].imag);
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
	result->data = (unsigned char*)alloc_data(result->w * result->h * sizeof(double));

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




void normed_cross_correlation_fft(struct ImageData* pattern, struct ImageData* target, struct ImageData* result)
{
	struct ImageData result_cc;
	cross_correlation_fft(pattern, target, &result_cc);
	normalize_cross_correlation(pattern, target, &result_cc, result);
	free_image(&result_cc);
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


