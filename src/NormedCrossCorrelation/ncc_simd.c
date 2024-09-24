#include "ncc_simd.h"
#include "ncc_common.h"
#include <emmintrin.h>
#include <immintrin.h>
#include <stdio.h>
#include <math.h>

// calculate integral image
static void integral_image(struct ImageData* src, struct ImageData* dst)
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

//From ImageShop
// 4個有符號的32位的數據相加的和。
inline int _mm_hsum_epi32(__m128i V)      // V3 V2 V1 V0
{
	// 實測這個速度要快些，_mm_extract_epi32最慢。
	__m128i T = _mm_add_epi32(V, _mm_srli_si128(V, 8));  // V3+V1   V2+V0  V1  V0  
	T = _mm_add_epi32(T, _mm_srli_si128(T, 4));    // V3+V1+V2+V0  V2+V0+V1 V1+V0 V0 
	return _mm_cvtsi128_si32(T);       // 提取低位 
}
// 基於SSE的字節數據的乘法。
// <param name="Kernel">需要卷積的核矩陣。 </param>
// <param name="Conv">卷積矩陣。 </param>
// <param name="Length">矩陣所有元素的長度。 </param>
inline int IM_Conv_SIMD(unsigned char* pCharKernel, unsigned char* pCharConv, int iLength)
{
	const int iBlockSize = 16, Block = iLength / iBlockSize;
	__m128i SumV = _mm_setzero_si128();
	__m128i Zero = _mm_setzero_si128();
	for (int Y = 0; Y < Block * iBlockSize; Y += iBlockSize)
	{
		__m128i SrcK = _mm_loadu_si128((__m128i*)(pCharKernel + Y));
		__m128i SrcC = _mm_loadu_si128((__m128i*)(pCharConv + Y));
		__m128i SrcK_L = _mm_unpacklo_epi8(SrcK, Zero);
		__m128i SrcK_H = _mm_unpackhi_epi8(SrcK, Zero);
		__m128i SrcC_L = _mm_unpacklo_epi8(SrcC, Zero);
		__m128i SrcC_H = _mm_unpackhi_epi8(SrcC, Zero);
		__m128i SumT = _mm_add_epi32(_mm_madd_epi16(SrcK_L, SrcC_L), _mm_madd_epi16(SrcK_H, SrcC_H));
		SumV = _mm_add_epi32(SumV, SumT);
	}
	int Sum = _mm_hsum_epi32(SumV);
	for (int Y = Block * iBlockSize; Y < iLength; Y++)
	{
		Sum += pCharKernel[Y] * pCharConv[Y];
	}
	return Sum;
}

inline int _mm256_hsum_epi32(__m256i V)      // V7 V6 V5 V4 V3 V2 V1 V0
{
	__m256i T = _mm256_add_epi32(V, _mm256_permute2f128_si256(V, V, 1));  // V7+V3 V6+V2 V5+V1 V4+V0 V3 V2 V1 V0
	T = _mm256_add_epi32(T, _mm256_shuffle_epi32(T, 0x4E));    // V7+V3+V6+V2 V6+V2+V5+V1 V5+V1+V4+V0 V4+V0+V3+V7 V3 V2 V1 V0
	T = _mm256_add_epi32(T, _mm256_shuffle_epi32(T, 0xB1));    // V7+V3+V6+V2+V5+V1+V4+V0 V6+V2+V5+V1+V4+V0+V3+V7 V5+V1+V4+V0+V3+V7+V6+V2 V4+V0+V3+V7+V6+V2+V5+V1 V3 V2 V1 V0
	return _mm256_extract_epi32(T, 0);       // 提取低位 
}

inline int IM_Conv_AVX2(unsigned char* pCharKernel, unsigned char* pCharConv, int iLength)
{
	const int iBlockSize = 32, Block = iLength / iBlockSize;
	__m256i SumV = _mm256_setzero_si256();
	__m256i Zero = _mm256_setzero_si256();

	for (int Y = 0; Y < Block * iBlockSize; Y += iBlockSize)
	{
		__m256i SrcK = _mm256_loadu_si256((__m256i*)(pCharKernel + Y));
		__m256i SrcC = _mm256_loadu_si256((__m256i*)(pCharConv + Y));
		// spilt head 128bits and tail 128bits of SrcK and SrcC
		__m256i SrcK_L = _mm256_unpacklo_epi8(SrcK, Zero);
		__m256i SrcK_H = _mm256_unpackhi_epi8(SrcK, Zero);

		__m256i SrcC_L = _mm256_unpacklo_epi8(SrcC, Zero);
		__m256i SrcC_H = _mm256_unpackhi_epi8(SrcC, Zero);

		__m256i Sum = _mm256_add_epi32(_mm256_madd_epi16(SrcK_L, SrcC_L), _mm256_madd_epi16(SrcK_H, SrcC_H));
		SumV = _mm256_add_epi32(SumV, Sum);
	}
	/*int Sum = 
		_mm256_extract_epi32(SumV, 0) + 
		_mm256_extract_epi32(SumV, 1) + 
		_mm256_extract_epi32(SumV, 2) + 
		_mm256_extract_epi32(SumV, 3) +

		_mm256_extract_epi32(SumV, 4) +
		_mm256_extract_epi32(SumV, 5) +
		_mm256_extract_epi32(SumV, 6) +
		_mm256_extract_epi32(SumV, 7);*/
	int Sum = _mm256_hsum_epi32(SumV);
	for (int Y = Block * iBlockSize; Y < iLength; Y++)
	{
		Sum += pCharKernel[Y] * pCharConv[Y];
	}
	return Sum;
}

void normed_cross_correlation_simd(struct ImageData* pattern, struct ImageData* target, struct ImageData* result)
{
	unsigned char** p = (unsigned char**)malloc(sizeof(unsigned char*) * pattern->h);
	for (int i = 0; i < pattern->h; i++)
	{
		p[i] = (unsigned char*)(pattern->data + i * pattern->w);
	}

	unsigned char** tp = (unsigned char**)malloc(sizeof(unsigned char*) * target->h);
	for (int i = 0; i < target->h; i++)
	{
		tp[i] = (unsigned char*)(target->data + i * target->w);
	}

	struct ImageData result_cc;
	result_cc.w = target->w - pattern->w + 1;
	result_cc.h = target->h - pattern->h + 1;
	result_cc.ch = 1;
	result_cc.data = (unsigned char*)alloc_data(result_cc.w * result_cc.h * sizeof(float));

	float** cc = (float**)malloc(sizeof(float*) * result_cc.h);
	for (int i = 0; i < result_cc.h; i++)
	{
		cc[i] = (float*)(result_cc.data + i * result_cc.w * 4);
	}

	int y;
#pragma omp parallel for shared(p, cc, tp)
	for (y = 0; y < result_cc.h; y++)
	{
		float* rr = cc[y];
		for (int x = 0; x < result_cc.w; x++)
		{
			rr[x] = 0;
			for (int j = 0; j < pattern->h; ++j)
			{
				rr[x] = rr[x] + //IM_Conv_AVX2(p[j], tp[y + j] + x, pattern->w);
				              IM_Conv_SIMD(p[j], tp[y + j] + x, pattern->w);
			}
		}
	}

	normalize_cross_correlation(pattern, target, &result_cc, result);


	free_image(&result_cc);
	free(p);
	free(tp);
	free(cc);
}


