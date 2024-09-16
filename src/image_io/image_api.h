#ifndef __IMAGE_API_H__
#define __IMAGE_API_H__

typedef void *(MallocBuffer)(int);
#if __cplusplus
extern "C" {
#endif
	__declspec(dllexport) void* read_image_ascii(const char *fileName, int *w, int *h, int *ch, MallocBuffer malloc);

	__declspec(dllexport) void write_image_ascii(const char *fileName, const int w, const int h, const int ch, void *data);

	__declspec(dllexport) int read_image_file_header(const char* fn);

#if _WIN32
	__declspec(dllexport) void* read_image_wstring(const wchar_t* fileName, int* w, int* h, int* ch, MallocBuffer malloc);
	__declspec(dllexport) int read_image_file_header_wstring(const wchar_t* fn);
#endif

	__declspec(dllexport) void* encode_image_to_jpeg_byte_stream(void *buffer, const int w, const int h, const int ch, unsigned long *output_length);

	__declspec(dllexport) void* decode_jpeg_byte_stream_to_image(void *byte_stream, unsigned long len, int *w, int *h, int *ch, MallocBuffer _malloc);
#if __cplusplus
}
#endif

#endif
