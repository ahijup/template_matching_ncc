#include "image_api.h"

#include "LibBmp\libbmp.h"
#include "LibPng\png.h"
#include "LibPng\pngpriv.h"
#include "LibTiff\tiff.h"
#include "LibTiff\tiffiop.h"
#include "LibTiff\tiffio.h"
#include "LibJpeg\jpeglib.h"

#include <stdio.h>
#include <string.h>
#include <xmmintrin.h>

#include <Shlwapi.h>
#include <comdef.h>

struct image_obj
{
	void *data;
	int w, h, ch;
	int res;
};

/**
	reading file from unicode filename
*/
inline image_obj ReadBmpFile(const wchar_t *fn, MallocBuffer malloc_)
{
	image_obj obj;
	void *data = NULL;

	bmp_header header;
	bmp_error err = bmp_img_read_ext_2(fn, header, data, (MallocImageBuffer)malloc_);

	obj.ch = header.biBitCount / 8;
	obj.w = header.biWidth;
	obj.h = header.biHeight;
	obj.res = err;
	obj.data = data;

	return obj;
};

inline image_obj ReadPngFile(const wchar_t *file, MallocBuffer malloc_) 
{
	image_obj obj;
	png_image image;

	/* Initialize the 'png_image' structure. */
	memset(&image, 0, (sizeof image));
	image.version = PNG_IMAGE_VERSION;

	if (png_image_begin_read_from_file(&image, file) != 0) 
	{
		png_control* png_type = image.opaque;
		if (png_type->info_ptr->bit_depth != 8) 
		{
			obj.res = -3;
		}
		else 
		{
			int w = image.width, h = image.height, ch = image.opaque->info_ptr->channels;
			if (ch == 4) 
			{
				image.format = PNG_FORMAT_BGRA;
			}
			else if (ch == 3) 
			{
				image.format = PNG_FORMAT_BGR;
			}
			else if (ch == 1) 
			{
				image.format = PNG_FORMAT_GRAY;
			}

			png_bytep buffer = (png_bytep)malloc_(PNG_IMAGE_SIZE(image));

			if (buffer != NULL) 
			{
				int png_res = png_image_finish_read(&image, NULL, buffer,
					0 /*row strike*/, NULL /*color map*/);

				obj.res = png_res;

				if (png_res == 0) 
				{
					obj.res = -2;
				}
				else 
				{
					obj.w = w;
					obj.h = h;
					obj.ch = ch;
					obj.data = buffer;
					obj.res = 0;
				}
			}
			else 
			{
				obj.res = -1;
				png_image_free(&image);
			}
		}
	}

	return obj;
};

inline image_obj ReadTiffFile(const wchar_t *file, MallocBuffer malloc_) 
{
	image_obj obj;
	TIFF *tif = TIFFOpenW(file, "r");
	if (tif) 
	{
		char emsg[1024];
		if (tif->tif_dir.td_bitspersample != 8 ||
			(tif->tif_dir.td_imagedepth != 1 && tif->tif_dir.td_imagedepth != 3 && tif->tif_dir.td_imagedepth != 4))
		{
			obj.res = -1;
			TIFFClose(tif);
		}
		else 
		{
			if (TIFFRGBAImageOK(tif, emsg)) 
			{
				TIFFRGBAImage img;
				char emsg[1024];

				if (TIFFRGBAImageBegin(&img, tif, -1, emsg)) {
					TIFFDirectory* td = &tif->tif_dir;
					size_t npixels;
					uint32* raster;

					npixels = img.width;
					raster = (uint32*)_TIFFmalloc(npixels * sizeof(uint32));
					
					obj.w = img.width;
					obj.h = img.height;
					obj.ch = td->td_imagedepth;
					void *data = malloc_(td->td_imagedepth * obj.w * obj.h);

					if (raster != NULL) 
					{
						//if (CWINDOW_DLL_LAST_ERROR == 0) 
						{
							bool succ = true;
							for (int i = 0; i < img.height; ++i, ++img.row_offset) 
							{
								if (!TIFFRGBAImageGet(&img, raster, img.width, 1)) 
								{
									succ = false;
									break;
								}
								uint32 *dstPtr = raster;
								BYTE* dst = ((BYTE *)(data) + obj.ch * obj.w * i);
								switch (td->td_imagedepth)
								{
								case 4:
									memcpy(dst, dstPtr, sizeof(uint32) * npixels);
									break;
								case 3:
									for (int j = 0; j < img.width; ++j, ++dstPtr) 
									{
										dst[j * 3 + 0] = TIFFGetB(*dstPtr);
										dst[j * 3 + 1] = TIFFGetG(*dstPtr);
										dst[j * 3 + 2] = TIFFGetR(*dstPtr);
									}
									break;
								case 1:
									for (int j = 0; j < img.width; ++j, ++dstPtr)
									{
										dst[j] = TIFFGetB(*dstPtr);
									}
									break;
								}
							}
							obj.data = data;
							
							if (!succ) 
							{
								obj.res = -2;
							}
							else
							{
								obj.res = 0;
							}
						}
						
						_TIFFfree(raster);
					}
					else 
					{
						obj.res = -4;
					}
					TIFFRGBAImageEnd(&img);
					TIFFClose(tif);
				}
				else 
				{
					obj.res = -1;
					TIFFClose(tif);
				}
			}
			else 
			{
				obj.res = -5;
				TIFFClose(tif);
			}
		}

	}
	else 
	{
		obj.res = -6;
	}
	return obj;
};

inline image_obj ReadJpegFile(const wchar_t *file, MallocBuffer malloc_) 
{
	image_obj obj;
	/* This struct contains the JPEG decompression parameters and pointers to
	* working space (which is allocated as needed by the JPEG library).
	*/
	struct jpeg_decompress_struct cinfo;
	/* More stuff */
	FILE * infile;		/* source file */
	JSAMPARRAY buffer;		/* Output row buffer */
	int row_stride;		/* physical row width in output buffer */

						//if ((infile = fopen(file, "rb")) == NULL) {
	if ((infile = _wfopen(file, L"rb")) == NULL) {
		//fprintf(stderr, "can't open %s\n", file);
		obj.res = -1;
		return obj;
	}

	struct jpeg_error_mgr pub;	/* "public" fields */

								// is single-thread application needed exit call back?
								//pub.error_exit = my_error_exit;
								/* We set up the normal JPEG error routines, then override error_exit. */
	cinfo.err = jpeg_std_error(&pub);

	/* Now we can initialize the JPEG decompression object. */
	jpeg_create_decompress(&cinfo);

	/* Step 2: specify data source (eg, a file) */
	jpeg_stdio_src(&cinfo, infile);

	/* Step 3: read file parameters with jpeg_read_header() */
	if (jpeg_read_header(&cinfo, TRUE) != JPEG_HEADER_OK) {
		jpeg_destroy_decompress(&cinfo);
		fclose(infile);
		obj.res = -2;
		return obj;
	};
	/* We can ignore the return value from jpeg_read_header since
	*   (a) suspension is not possible with the stdio data source, and
	*   (b) we passed TRUE to reject a tables-only JPEG file as an error.
	* See libjpeg.txt for more info.
	*/

	/* Step 4: set parameters for decompression */

	/* In this example, we don't need to change any of the defaults set by
	* jpeg_read_header(), so we do nothing here.
	*/

	/* Step 5: Start decompressor */

	if (!jpeg_start_decompress(&cinfo)) 
	{
		jpeg_destroy_decompress(&cinfo);
		fclose(infile);
		obj.res = -3;
		return obj;
	}
	/* We can ignore the return value since suspension is not possible
	* with the stdio data source.
	*/

	/* We may need to do some setup of our own at this point before reading
	* the data.  After jpeg_start_decompress() we have the correct scaled
	* output image dimensions available, as well as the output colormap
	* if we asked for color quantization.
	* In this example, we need to make an output work buffer of the right size.
	*/
	/* JSAMPLEs per row in output buffer */
	row_stride = cinfo.output_width * cinfo.output_components;



	//HCLIBIMAGE ret = CreateImage(cinfo.output_width, cinfo.output_height, cinfo.output_components, "BYTE");
	void *data = malloc_(cinfo.output_components * cinfo.output_width * cinfo.output_height);
	obj.ch = cinfo.output_components;
	obj.w = cinfo.output_width;
	obj.h = cinfo.output_height;
	obj.data = data;
	
	/* Make a one-row-high sample array that will go away when done with image */
	buffer = (*cinfo.mem->alloc_sarray)
		((j_common_ptr)&cinfo, JPOOL_IMAGE, row_stride, 1);

	bool readingError = false;
	if (cinfo.output_components == 1 ||
		cinfo.output_components == 3)
	{
		/* Step 6: while (scan lines remain to be read) */
		/*           jpeg_read_scanlines(...); */

		BYTE* dataPtr = (BYTE *)data;
		/* Here we use the library's state variable cinfo.output_scanline as the
		* loop counter, so that we don't have to keep track ourselves.
		*/
		while (cinfo.output_scanline < cinfo.output_height) {
			/* jpeg_read_scanlines expects an array of pointers to scanlines.
			* Here the array is only one element long, but you could ask for
			* more than one scanline at a time if that's more convenient.
			*/
			int res = jpeg_read_scanlines(&cinfo, buffer, 1);
			JSAMPROW bb = buffer[0];
			if (res != 1) {
				
				readingError = true;
				break;
			}
			/* Assume put_scanline_someplace wants a pointer and sample count. */
			//put_scanline_someplace(buffer[0], row_stride);


			if (cinfo.output_components == 1) 
			{
				memcpy(dataPtr, bb, row_stride);
			}
			else if (cinfo.output_components == 3) 
			{
				int j = 0;
				BYTE *srcPtr = bb;
				while (j < obj.w) 
				{
					dataPtr[j * 3] = srcPtr[j * 3 + 2];
					dataPtr[j * 3 + 1] = srcPtr[j * 3 + 1];
					dataPtr[j * 3 + 2] = srcPtr[j * 3];
					++j;
				}
				//memcpy( dataPtr , bb , row_stride) ;
			}
			dataPtr += (obj.w * obj.ch);
		};
	}
	else 
	{
		readingError = true;
		obj.res = -5;
		//CWINDOW_DLL_LAST_ERROR = CIMAGE_JPEG_NOT_SUPPORT_CHANNEL;
	}
	/* Step 7: Finish decompression */

	jpeg_finish_decompress(&cinfo);
	/* We can ignore the return value since suspension is not possible
	* with the stdio data source.
	*/

	/* Step 8: Release JPEG decompression object */

	/* This is an important step since it will release a good deal of memory. */
	jpeg_destroy_decompress(&cinfo);

	/* After finish_decompress, we can close the input file.
	* Here we postpone it until after no more JPEG errors are possible,
	* so as to simplify the setjmp error logic above.  (Actually, I don't
	* think that jpeg_destroy can do an error exit, but why assume anything...)
	*/
	fclose(infile);

	if (readingError) 
	{
		obj.res = -4;
	}
	else 
	{
		obj.res = 0;
	}
	return obj;
};
//////////////////////////////////////////////////////////////////////////

inline image_obj ReadBmpFile_ASCII(const char *fn, MallocBuffer malloc_)
{
	image_obj obj;
	void *data = NULL;

	bmp_header header;
	bmp_error err = bmp_img_read_ext_3(fn, header, data, (MallocImageBuffer)malloc_);

	obj.ch = header.biBitCount / 8;
	obj.w = header.biWidth;
	obj.h = header.biHeight;
	obj.res = err;
	obj.data = data;

	return obj;
};

inline image_obj ReadPngFile_ASCII(const char *file, MallocBuffer malloc_)
{
	image_obj obj;
	png_image image;

	/* Initialize the 'png_image' structure. */
	memset(&image, 0, (sizeof image));
	image.version = PNG_IMAGE_VERSION;

	FILE *fp = fopen(file, "rb");

	if (png_image_begin_read_from_stdio(&image, fp) != 0)
	{
		png_control* png_type = image.opaque;
		if (png_type->info_ptr->bit_depth != 8)
		{
			obj.res = -3;
		}
		else
		{
			int w = image.width, h = image.height, ch = image.opaque->info_ptr->channels;
			if (ch == 4)
			{
				image.format = PNG_FORMAT_BGRA;
			}
			else if (ch == 3)
			{
				image.format = PNG_FORMAT_BGR;
			}
			else if (ch == 1)
			{
				image.format = PNG_FORMAT_GRAY;
			}

			png_bytep buffer = (png_bytep)malloc_(PNG_IMAGE_SIZE(image));

			if (buffer != NULL)
			{
				int png_res = png_image_finish_read(&image, NULL, buffer,
					0 /*row strike*/, NULL /*color map*/);

				obj.res = png_res;

				if (png_res == 0)
				{
					obj.res = -2;
				}
				else
				{
					obj.w = w;
					obj.h = h;
					obj.ch = ch;
					obj.data = buffer;
					obj.res = 0;
				}
			}
			else
			{
				obj.res = -1;
				png_image_free(&image);
			}
		}
	}
	fclose(fp);

	return obj;
};

inline image_obj ReadTiffFile_ASCII(const char *file, MallocBuffer malloc_)
{
	image_obj obj;
	TIFF *tif = TIFFOpen(file, "r");
	if (tif)
	{
		char emsg[1024];
		if (tif->tif_dir.td_bitspersample != 8 ||
			(tif->tif_dir.td_imagedepth != 1 && tif->tif_dir.td_imagedepth != 3 && tif->tif_dir.td_imagedepth != 4))
		{
			obj.res = -1;
			TIFFClose(tif);
		}
		else
		{
			if (TIFFRGBAImageOK(tif, emsg))
			{
				TIFFRGBAImage img;
				char emsg[1024];

				if (TIFFRGBAImageBegin(&img, tif, -1, emsg)) {
					TIFFDirectory* td = &tif->tif_dir;
					size_t npixels;
					uint32* raster;

					npixels = img.width;
					raster = (uint32*)_TIFFmalloc(npixels * sizeof(uint32));

					obj.w = img.width;
					obj.h = img.height;
					obj.ch = td->td_imagedepth;
					void *data = malloc_(td->td_imagedepth * obj.w * obj.h);

					if (raster != NULL)
					{
						//if (CWINDOW_DLL_LAST_ERROR == 0) 
						{
							bool succ = true;
							for (int i = 0; i < img.height; ++i, ++img.row_offset)
							{
								if (!TIFFRGBAImageGet(&img, raster, img.width, 1))
								{
									succ = false;
									break;
								}
								uint32 *dstPtr = raster;
								BYTE* dst = ((BYTE *)(data)+obj.ch * obj.w * i);
								switch (td->td_imagedepth)
								{
								case 4:
									memcpy(dst, dstPtr, sizeof(uint32) * npixels);
									break;
								case 3:
									for (int j = 0; j < img.width; ++j, ++dstPtr)
									{
										dst[j * 3 + 0] = TIFFGetB(*dstPtr);
										dst[j * 3 + 1] = TIFFGetG(*dstPtr);
										dst[j * 3 + 2] = TIFFGetR(*dstPtr);
									}
									break;
								case 1:
									for (int j = 0; j < img.width; ++j, ++dstPtr)
									{
										dst[j] = TIFFGetB(*dstPtr);
									}
									break;
								}
							}
							obj.data = data;

							if (!succ)
							{
								obj.res = -2;
							}
							else
							{
								obj.res = 0;
							}
						}

						_TIFFfree(raster);
					}
					else
					{
						obj.res = -4;
					}
					TIFFRGBAImageEnd(&img);
					TIFFClose(tif);
				}
				else
				{
					obj.res = -1;
					TIFFClose(tif);
				}
			}
			else
			{
				obj.res = -5;
				TIFFClose(tif);
			}
		}

	}
	else
	{
		obj.res = -6;
	}
	return obj;
};

inline image_obj ReadJpegFile_ASCII(const char *file, MallocBuffer malloc_)
{
	image_obj obj;
	/* This struct contains the JPEG decompression parameters and pointers to
	* working space (which is allocated as needed by the JPEG library).
	*/
	struct jpeg_decompress_struct cinfo;
	/* More stuff */
	FILE * infile;		/* source file */
	JSAMPARRAY buffer;		/* Output row buffer */
	int row_stride;		/* physical row width in output buffer */

						//if ((infile = fopen(file, "rb")) == NULL) {
	if ((infile = fopen(file, "rb")) == NULL) {
		//fprintf(stderr, "can't open %s\n", file);
		obj.res = -1;
		return obj;
	}

	struct jpeg_error_mgr pub;	/* "public" fields */

								// is single-thread application needed exit call back?
								//pub.error_exit = my_error_exit;
								/* We set up the normal JPEG error routines, then override error_exit. */
	cinfo.err = jpeg_std_error(&pub);

	/* Now we can initialize the JPEG decompression object. */
	jpeg_create_decompress(&cinfo);

	/* Step 2: specify data source (eg, a file) */
	jpeg_stdio_src(&cinfo, infile);

	/* Step 3: read file parameters with jpeg_read_header() */
	if (jpeg_read_header(&cinfo, TRUE) != JPEG_HEADER_OK) {
		jpeg_destroy_decompress(&cinfo);
		fclose(infile);
		obj.res = -2;
		return obj;
	};
	/* We can ignore the return value from jpeg_read_header since
	*   (a) suspension is not possible with the stdio data source, and
	*   (b) we passed TRUE to reject a tables-only JPEG file as an error.
	* See libjpeg.txt for more info.
	*/

	/* Step 4: set parameters for decompression */

	/* In this example, we don't need to change any of the defaults set by
	* jpeg_read_header(), so we do nothing here.
	*/

	/* Step 5: Start decompressor */

	if (!jpeg_start_decompress(&cinfo))
	{
		jpeg_destroy_decompress(&cinfo);
		fclose(infile);
		obj.res = -3;
		return obj;
	}
	/* We can ignore the return value since suspension is not possible
	* with the stdio data source.
	*/

	/* We may need to do some setup of our own at this point before reading
	* the data.  After jpeg_start_decompress() we have the correct scaled
	* output image dimensions available, as well as the output colormap
	* if we asked for color quantization.
	* In this example, we need to make an output work buffer of the right size.
	*/
	/* JSAMPLEs per row in output buffer */
	row_stride = cinfo.output_width * cinfo.output_components;



	//HCLIBIMAGE ret = CreateImage(cinfo.output_width, cinfo.output_height, cinfo.output_components, "BYTE");
	void *data = malloc_(cinfo.output_components * cinfo.output_width * cinfo.output_height);
	obj.ch = cinfo.output_components;
	obj.w = cinfo.output_width;
	obj.h = cinfo.output_height;
	obj.data = data;

	/* Make a one-row-high sample array that will go away when done with image */
	buffer = (*cinfo.mem->alloc_sarray)
		((j_common_ptr)&cinfo, JPOOL_IMAGE, row_stride, 1);

	bool readingError = false;
	if (cinfo.output_components == 1 ||
		cinfo.output_components == 3)
	{
		/* Step 6: while (scan lines remain to be read) */
		/*           jpeg_read_scanlines(...); */

		BYTE* dataPtr = (BYTE *)data;
		/* Here we use the library's state variable cinfo.output_scanline as the
		* loop counter, so that we don't have to keep track ourselves.
		*/
		while (cinfo.output_scanline < cinfo.output_height) {
			/* jpeg_read_scanlines expects an array of pointers to scanlines.
			* Here the array is only one element long, but you could ask for
			* more than one scanline at a time if that's more convenient.
			*/
			int res = jpeg_read_scanlines(&cinfo, buffer, 1);
			JSAMPROW bb = buffer[0];
			if (res != 1) {

				readingError = true;
				break;
			}
			/* Assume put_scanline_someplace wants a pointer and sample count. */
			//put_scanline_someplace(buffer[0], row_stride);


			if (cinfo.output_components == 1)
			{
				memcpy(dataPtr, bb, row_stride);
			}
			else if (cinfo.output_components == 3)
			{
				int j = 0;
				BYTE *srcPtr = bb;
				while (j < obj.w)
				{
					dataPtr[j * 3] = srcPtr[j * 3 + 2];
					dataPtr[j * 3 + 1] = srcPtr[j * 3 + 1];
					dataPtr[j * 3 + 2] = srcPtr[j * 3];
					++j;
				}
				//memcpy( dataPtr , bb , row_stride) ;
			}
			dataPtr += (obj.w * obj.ch);
		};
	}
	else
	{
		readingError = true;
		obj.res = -5;
		//CWINDOW_DLL_LAST_ERROR = CIMAGE_JPEG_NOT_SUPPORT_CHANNEL;
	}
	/* Step 7: Finish decompression */

	jpeg_finish_decompress(&cinfo);
	/* We can ignore the return value since suspension is not possible
	* with the stdio data source.
	*/

	/* Step 8: Release JPEG decompression object */

	/* This is an important step since it will release a good deal of memory. */
	jpeg_destroy_decompress(&cinfo);

	/* After finish_decompress, we can close the input file.
	* Here we postpone it until after no more JPEG errors are possible,
	* so as to simplify the setjmp error logic above.  (Actually, I don't
	* think that jpeg_destroy can do an error exit, but why assume anything...)
	*/
	fclose(infile);

	if (readingError)
	{
		obj.res = -4;
	}
	else
	{
		obj.res = 0;
	}
	return obj;
};


int read_image_type(const unsigned char *ptr)
{
	//jpeg
	if (ptr[0] == 0xff && ptr[1] == 0xd8)
	{
		return 0;
	}
	//png
	else if (
		ptr[0] == 0x89 || ptr[1] == 0x50 ||
		ptr[2] == 0x4E || ptr[3] == 0x47 ||
		ptr[4] == 0x0D || ptr[5] == 0x0A ||
		ptr[6] == 0x1A || ptr[7] == 0x0A)
	{
		return 1;
	}
	//bmp
	else if (ptr[0] == 0x42 && ptr[1] == 0x4D)
	{
		return 2;
	}
	//tif
	else if (
		(ptr[0] == 0x4D && ptr[1] == 0x4D) ||
		(ptr[0] == 0x49 && ptr[1] == 0x49))
	{
		return 3;
	}
	return -1;
}


int read_image_file_header(const char *fn)
{
	unsigned char hdr[12];
	FILE *fp = fopen(fn, "r");
	if (fp)
	{
		fread(hdr, 1, 12, fp);
		fclose(fp);

		return read_image_type(hdr);
	}
	else
	{
		return -2;
	}
}

#if _WIN32
int read_image_file_header_wstring(const wchar_t* fn)
{
	unsigned char hdr[12];
	FILE* fp = _wfopen(fn, L"r");
	if (fp)
	{
		fread(hdr, 1, 12, fp);
		fclose(fp);

		return read_image_type(hdr);
	}
	else
	{
		return -2;
	}
}
#endif
void* read_image_ascii(const char *fileName, int *w, int *h, int *ch, MallocBuffer malloc_)
{
	_bstr_t b(fileName);

	LPSTR ext = PathFindExtensionA(b);

	image_obj obj;

	if (ext != NULL) 
	{
		size_t sLen = strlen(ext);
		for (size_t i = 0; i < sLen; ++i) 
		{
			if (ext[i] >= 'a' && ext[i] <= 'z') 
			{
				ext[i] = ext[i] - 'a' + 'A';
			}
		}


		/*
		if (strcmp(ext, ".BMP") == 0) 
		{
			obj = ReadBmpFile_ASCII(fileName, malloc_);
		}
		else if (strcmp(ext, ".PNG") == 0) 
		{
			obj = ReadPngFile_ASCII(fileName, malloc_);
		}
		else if (strcmp(ext, ".TIF") == 0 || strcmp(ext, ".TIFF") == 0) 
		{
			obj = ReadTiffFile_ASCII(fileName, malloc_);
		}
		else if (strcmp(ext, ".JPEG") == 0 || strcmp(ext, ".JPG") == 0) 
		{
			obj = ReadJpegFile_ASCII(fileName , malloc_);
		}
		*/

		int type = read_image_file_header(fileName);
		switch (type)
		{
		case 0:
			obj = ReadJpegFile_ASCII(fileName, malloc_);
			break;
		case 1:
			obj = ReadPngFile_ASCII(fileName, malloc_);
			break;
		case 2:
			obj = ReadBmpFile_ASCII(fileName, malloc_);
			break;
		case 3:
			obj = ReadTiffFile_ASCII(fileName, malloc_);
			break;
		default:
			return 0;
		}
	}
	*w = obj.w;
	*h = obj.h;
	*ch = obj.ch;

	return obj.data;
};

// ############# write image file ##################

inline bool WriteBmpFile(void* img, int w, int h, int ch, const char *fn) {
	bmp_header wHeader;
	bmp_header_init_df_2(&wHeader, w, h, ch);

	bmp_error st = bmp_img_write_2(wHeader, img, fn);

	if (st == BMP_OK) 
	{
		return true;
	}
	return false;
};

inline bool WriteBmpFile2(void* img, int w, int h, int ch, const wchar_t *fn) {
	bmp_header wHeader;
	bmp_header_init_df_2(&wHeader, w, h, ch);

	bmp_error st = bmp_img_write_3(wHeader, img, fn);

	if (st == BMP_OK) 
	{
		return true;
	}
	return false;
};

inline bool WriteTiffFile(void* img, int w, int h, int ch, const char *fn) {
	TIFF* tifW = TIFFOpen(fn, "wb");
	if (tifW) {
		TIFFSetField(tifW, TIFFTAG_IMAGEWIDTH, w);
		TIFFSetField(tifW, TIFFTAG_IMAGELENGTH, h);

		TIFFSetField(tifW, TIFFTAG_SAMPLESPERPIXEL, ch == 4 ? 3 : ch);
		TIFFSetField(tifW, TIFFTAG_BITSPERSAMPLE, 8);
		TIFFSetField(tifW, TIFFTAG_ORIENTATION, 4);    // set the origin of the image.
													   //   Some other essential fields to set that you do not have to understand for now.
		TIFFSetField(tifW, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
		TIFFSetField(tifW, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
		uint32 ss = TIFFDefaultStripSize(tifW, w);
		// We set the strip size of the file to be size of one row of pixels
		TIFFSetField(tifW, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize(tifW, w));
		TIFFSetField(tifW, TIFFTAG_SOFTWARE, "");

		BYTE *tmpRow = (BYTE *)malloc(w * ch == 4 ? 3 : ch);

		for (int i = 0; i < h; ++i) 
		{
			BYTE* buffer = ((BYTE *)img) + i * w * ch;
			if (ch == 1 || ch == 3) 
			{
				memcpy(tmpRow, buffer, w * ch == 4 ? 3 : ch);
			}
			else 
			{
				for (int j = 0; j < w; ++j) 
				{
					tmpRow[j * 3] = buffer[j * 4];
					tmpRow[j * 3 + 1] = buffer[j * 4 + 1];
					tmpRow[j * 3 + 2] = buffer[j * 4 + 2];
				}
			}

			if (TIFFWriteScanline(tifW, tmpRow, i, 0) < 0) 
			{
				break;
			}
		}

		free(tmpRow);
		TIFFClose(tifW);
		return true;

	}

	return false;
};

inline bool WriteTiffFile2(void* img, int w, int h, int ch, const wchar_t *fn) {
	TIFF* tifW = TIFFOpenW(fn, "wb");
	if (tifW) {
		TIFFSetField(tifW, TIFFTAG_IMAGEWIDTH, w);
		TIFFSetField(tifW, TIFFTAG_IMAGELENGTH, h);

		TIFFSetField(tifW, TIFFTAG_SAMPLESPERPIXEL, ch == 4 ? 3 : ch);
		TIFFSetField(tifW, TIFFTAG_BITSPERSAMPLE, 8);
		TIFFSetField(tifW, TIFFTAG_ORIENTATION, 4);    // set the origin of the image.
													   //   Some other essential fields to set that you do not have to understand for now.
		TIFFSetField(tifW, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
		TIFFSetField(tifW, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
		uint32 ss = TIFFDefaultStripSize(tifW, w);
		// We set the strip size of the file to be size of one row of pixels
		TIFFSetField(tifW, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize(tifW, w));
		TIFFSetField(tifW, TIFFTAG_SOFTWARE, "");

		BYTE *tmpRow = (BYTE *)malloc(w * ch == 4 ? 3 : ch);

		for (int i = 0; i < h; ++i) {
			BYTE* buffer = ((BYTE *)img) + i * w * ch;
			if (ch == 1 || ch == 3) {
				memcpy(tmpRow, buffer, w * ch == 4 ? 3 : ch);
			}
			else {
				for (int j = 0; j < w; ++j) {
					tmpRow[j * 3] = buffer[j * 4];
					tmpRow[j * 3 + 1] = buffer[j * 4 + 1];
					tmpRow[j * 3 + 2] = buffer[j * 4 + 2];
				}
			}

			if (TIFFWriteScanline(tifW, tmpRow, i, 0) < 0) 
			{
				break;
			}
		}

		free(tmpRow);
		TIFFClose(tifW);
		return true;

	}

	return false;
};

inline bool WritePngFile(void* img, int w, int h, int channel, const char *fn, const char *params) {
	png_image imgToWrite;
	/* Initialize the 'png_image' structure. */
	memset(&imgToWrite, 0, (sizeof imgToWrite));
	imgToWrite.version = PNG_IMAGE_VERSION;
	int ch = channel;
	switch (ch)
	{
	case 1:
		imgToWrite.format = PNG_FORMAT_GRAY;
		break;
	case 3:
		imgToWrite.format = PNG_FORMAT_BGR;
		break;
	case 4:
		imgToWrite.format = PNG_FORMAT_BGRA;
		break;
	}

	imgToWrite.width = w;
	imgToWrite.height = h;

	int res = png_image_write_to_file(&imgToWrite,
		fn, 0, img, 0, NULL);

	if (res == 0) 
	{
		return false;
	}
	else
	{
		return true;
	}
};

inline bool WritePngFile2(void* img, int w, int h, int channel, const wchar_t *fn, const char *params) {
	png_image imgToWrite;
	/* Initialize the 'png_image' structure. */
	memset(&imgToWrite, 0, (sizeof imgToWrite));
	imgToWrite.version = PNG_IMAGE_VERSION;
	int ch = channel;
	switch (ch)
	{
	case 1:
		imgToWrite.format = PNG_FORMAT_GRAY;
		break;
	case 3:
		imgToWrite.format = PNG_FORMAT_BGR;
		break;
	case 4:
		imgToWrite.format = PNG_FORMAT_BGRA;
		break;
	}

	imgToWrite.width = w;
	imgToWrite.height = h;

	FILE *fp = _wfopen(fn, L"wb");
	int res = -1;
	if (fp)
	{
		res = png_image_write_to_stdio(&imgToWrite,
			fp, 0, img, 0, NULL);
		fclose(fp);
	}
	if (res == 0) 
	{
		return false;
	}
	else {
		return true;
	}
};

inline char* TrimStr(char *ch) {
	int len = strlen(ch);
	char *ret = (char *)malloc(sizeof(char) * (len + 1));
	int i = 0, j = 0;
	while (i < len) {
		if (ch[i] != ' ') {
			ret[j++] = ch[i];
		}
		++i;
	}
	ch[j] = '\0';
	return ret;
};

inline wchar_t* TrimStr2(wchar_t *ch) {
	int len = wcslen(ch);
	wchar_t *ret = (wchar_t *)malloc(sizeof(wchar_t) * (len + 1));
	int i = 0, j = 0;
	while (i < len) {
		if (ch[i] != ' ') {
			ret[j++] = ch[i];
		}
		++i;
	}
	ch[j] = '\0';
	return ret;
};

inline bool WriteJpegFile(void* img, int w, int h, int ch, const char *fn, const char *params) {
	int quality = 100;
	char msg[2000];
	sprintf(msg, "%s", params);

	char *ptr = strtok(msg, " ");
	while (ptr) {
		char *buf = TrimStr(ptr);
		char data[1000];
		//Q=10
		switch ((char)toupper(buf[0]))
		{
		case 'Q':
			if (buf[1] == '=') {
				memcpy(buf + 2, data, strlen(buf) - 2);
				quality = atoi(data);
			}
			break;
		default:
			break;
		}
		free(buf);
		ptr = strtok(NULL, " ");
	}



	/* This struct contains the JPEG decompression parameters and pointers to
	* working space (which is allocated as needed by the JPEG library).
	*/
	struct jpeg_compress_struct cinfo;
	/* More stuff */
	FILE * outfile;		/* source file */
	int row_stride;		/* physical row width in output buffer */

	if ((outfile = fopen(fn, "wb")) == NULL) 
	{
		fprintf(stderr, "can't open %s\n", fn);
		return false;
	}

	/* Now we can initialize the JPEG decompression object. */
	jpeg_create_compress(&cinfo);

	struct jpeg_error_mgr pub;	/* "public" fields */

								/* We set up the normal JPEG error routines, then override error_exit. */
	cinfo.err = jpeg_std_error(&pub);

	/* Step 2: specify data destination (eg, a file) */
	jpeg_stdio_dest(&cinfo, outfile);


	cinfo.image_height = h;
	cinfo.image_width = w;
	cinfo.input_components = ch == 4 ? 3 : ch;

	if (ch == 1) 
	{
		cinfo.in_color_space = J_COLOR_SPACE::JCS_GRAYSCALE;
	}
	else if (ch == 3 || ch == 4) 
	{
		cinfo.in_color_space = J_COLOR_SPACE::JCS_RGB;
	}

	jpeg_set_defaults(&cinfo);
	jpeg_set_quality(&cinfo, quality, TRUE);

	BYTE *rowData = (BYTE *)malloc(sizeof(BYTE) * w * cinfo.input_components);
	if (!rowData) 
	{
		fclose(outfile);
		jpeg_destroy_compress(&cinfo);
		return false;
	};

	jpeg_start_compress(&cinfo, TRUE);
	row_stride = w * ch;
	bool succ = true;


	while (cinfo.next_scanline < h) {
		BYTE* buffer = ((BYTE *)img) + row_stride * cinfo.next_scanline;
		if (ch == 4 || ch == 3) {
			//shrink to 3 channels
			for (int j = 0; j < w; ++j) {
				rowData[j * 3] = buffer[j * ch + 2];
				rowData[j * 3 + 1] = buffer[j * ch + 1];
				rowData[j * 3 + 2] = buffer[j * ch + 0];
			}
		}
		else {

			memcpy(rowData, buffer, sizeof(BYTE) * row_stride);
		}
		int res = jpeg_write_scanlines(&cinfo, &rowData, 1);
		if (res != 1) 
		{
			succ = false;
			break;
		};
	};
	//clean-up all data
	free(rowData);
	jpeg_finish_compress(&cinfo);
	fclose(outfile);
	jpeg_destroy_compress(&cinfo);

	return succ;
};

inline bool WriteJpegFile2(void* img, int w, int h, int ch, const wchar_t *fn, const wchar_t *params) {
	int quality = 100;
	wchar_t msg[2000];
	wsprintfW(msg, L"%s", params);

	wchar_t *ptr = wcstok(msg, L" ");
	while (ptr) {
		wchar_t *buf = TrimStr2(ptr);
		wchar_t data[1000];
		//Q=10
		switch ((char)toupper(buf[0]))
		{
		case 'Q':
			if (buf[1] == '=') {
				memcpy(buf + 2, data, wcslen(buf) - 2);
				quality = _wtoi(data);
			}
			break;
		default:
			break;
		}
		free(buf);
		ptr = wcstok(NULL, L" ");
	}



	/* This struct contains the JPEG decompression parameters and pointers to
	* working space (which is allocated as needed by the JPEG library).
	*/
	struct jpeg_compress_struct cinfo;
	/* More stuff */
	FILE * outfile;		/* source file */
	int row_stride;		/* physical row width in output buffer */

	if ((outfile = _wfopen(fn, L"wb")) == NULL) 
	{
		fprintf(stderr, "can't open %s\n", fn);
		return false;
	}

	/* Now we can initialize the JPEG decompression object. */
	jpeg_create_compress(&cinfo);

	struct jpeg_error_mgr pub;	/* "public" fields */

								/* We set up the normal JPEG error routines, then override error_exit. */
	cinfo.err = jpeg_std_error(&pub);

	/* Step 2: specify data destination (eg, a file) */
	jpeg_stdio_dest(&cinfo, outfile);


	cinfo.image_height = h;
	cinfo.image_width = w;
	cinfo.input_components = ch == 4 ? 3 : ch;

	if (ch == 1) 
	{
		cinfo.in_color_space = J_COLOR_SPACE::JCS_GRAYSCALE;
	}
	else if (ch == 3 || ch == 4)
	{
		cinfo.in_color_space = J_COLOR_SPACE::JCS_RGB;
	}

	jpeg_set_defaults(&cinfo);
	jpeg_set_quality(&cinfo, quality, TRUE);

	BYTE *rowData = (BYTE *)malloc(sizeof(BYTE) * w * cinfo.input_components);
	if (!rowData)
	{
		fclose(outfile);
		jpeg_destroy_compress(&cinfo);
		return false;
	};

	jpeg_start_compress(&cinfo, TRUE);
	row_stride = w * ch;
	bool succ = true;


	while (cinfo.next_scanline < h) {
		BYTE* buffer = ((BYTE *)img) + row_stride * cinfo.next_scanline;
		if (ch == 4 || ch == 3) {
			//shrink to 3 channels
			for (int j = 0; j < w; ++j) {
				rowData[j * 3] = buffer[j * ch + 2];
				rowData[j * 3 + 1] = buffer[j * ch + 1];
				rowData[j * 3 + 2] = buffer[j * ch + 0];
			}
		}
		else {

			memcpy(rowData, buffer, sizeof(BYTE) * row_stride);
		}
		int res = jpeg_write_scanlines(&cinfo, &rowData, 1);
		if (res != 1) {
			succ = false;
			break;
		};
	};
	//clean-up all data
	free(rowData);
	jpeg_finish_compress(&cinfo);
	fclose(outfile);
	jpeg_destroy_compress(&cinfo);

	return succ;
};

// ########################################################

void write_image_ascii(const char *fileName, const int w, const int h, const int ch, void *data)
{
	_bstr_t b(fileName);

	LPSTR ext = PathFindExtensionA(b);

	image_obj obj = image_obj();

	if (ext != NULL)
	{
		size_t sLen = strlen(ext);
		for (size_t i = 0; i < sLen; ++i)
		{
			if (ext[i] >= 'a' && ext[i] <= 'z')
			{
				ext[i] = ext[i] - 'a' + 'A';
			}
		}

		int type = -1;
		
		if (strcmp(ext, ".BMP") == 0)
		{
			type = 2;
		}
		else if (strcmp(ext, ".PNG") == 0)
		{
			type = 1;
		}
		else if (strcmp(ext, ".TIF") == 0 || strcmp(ext, ".TIFF") == 0)
		{
			type = 3;
		}
		else if (strcmp(ext, ".JPEG") == 0 || strcmp(ext, ".JPG") == 0)
		{
			type = 0;
		}
		
		switch (type)
		{
		case 0:
			WriteJpegFile(data, w, h, ch, fileName, 0);
			break;
		case 1:
			WritePngFile(data, w, h, ch, fileName, 0);
			break;
		case 2:
			WriteBmpFile(data, w, h, ch, fileName);
			break;
		case 3:
			WriteTiffFile(data, w, h, ch, fileName);
			break;
		default:
			break;
		}
	}
};
inline bool WriteJpegFile2(void *buffer, const int w, const int h, const int ch, const wchar_t *fn) {
	int quality = 100;
	
	/* This struct contains the JPEG decompression parameters and pointers to
	* working space (which is allocated as needed by the JPEG library).
	*/
	struct jpeg_compress_struct cinfo;
	/* More stuff */
	FILE * outfile;		/* source file */
	int row_stride;		/* physical row width in output buffer */

	if ((outfile = _wfopen(fn, L"wb")) == NULL) {
		fprintf(stderr, "can't open %s\n", fn);
		return false;
	}

	/* Now we can initialize the JPEG decompression object. */
	jpeg_create_compress(&cinfo);

	struct jpeg_error_mgr pub;	/* "public" fields */

								/* We set up the normal JPEG error routines, then override error_exit. */
	cinfo.err = jpeg_std_error(&pub);

	/* Step 2: specify data destination (eg, a file) */
	jpeg_stdio_dest(&cinfo, outfile);


	cinfo.image_height = h;
	cinfo.image_width = w;
	cinfo.input_components = ch == 4 ? 3 : ch;

	if (ch == 1) {
		cinfo.in_color_space = J_COLOR_SPACE::JCS_GRAYSCALE;
	}
	else if (ch == 3 || ch == 4) {
		cinfo.in_color_space = J_COLOR_SPACE::JCS_RGB;
	}

	jpeg_set_defaults(&cinfo);
	jpeg_set_quality(&cinfo, quality, TRUE);

	BYTE *rowData = (BYTE *)malloc(sizeof(BYTE) * w * cinfo.input_components);
	if (!rowData) {
		fclose(outfile);
		jpeg_destroy_compress(&cinfo);
		return false;
	};

	jpeg_start_compress(&cinfo, TRUE);
	row_stride = w * ch;
	bool succ = true;


	while (cinfo.next_scanline < h) {
		BYTE* buffer = ((BYTE *)buffer) + row_stride * cinfo.next_scanline;
		if (ch == 4 || ch == 3) {
			//shrink to 3 channels
			for (int j = 0; j < w; ++j) {
				rowData[j * 3] = buffer[j * ch + 2];
				rowData[j * 3 + 1] = buffer[j * ch + 1];
				rowData[j * 3 + 2] = buffer[j * ch + 0];
			}
		}
		else {

			memcpy(rowData, buffer, sizeof(BYTE) * row_stride);
		}
		int res = jpeg_write_scanlines(&cinfo, &rowData, 1);
		if (res != 1) {
			succ = false;
			break;
		};
	};
	//clean-up all data
	free(rowData);
	jpeg_finish_compress(&cinfo);
	fclose(outfile);
	jpeg_destroy_compress(&cinfo);

	return succ;
};


void* encode_image_to_jpeg_byte_stream(void *data, const int w, const int h, const int ch, unsigned long *output_length)
{
	int quality = 100;

	/* This struct contains the JPEG decompression parameters and pointers to
	* working space (which is allocated as needed by the JPEG library).
	*/
	struct jpeg_compress_struct cinfo;
	/* More stuff */
	int row_stride;		/* physical row width in output buffer */
						/* Now we can initialize the JPEG decompression object. */
	jpeg_create_compress(&cinfo);

	struct jpeg_error_mgr pub;	/* "public" fields */

								/* Step 2: specify data destination (eg, a file) */
	unsigned char *tmpBuffer = NULL;
	jpeg_mem_dest(&cinfo, &tmpBuffer, output_length);
	//jpeg_stdio_dest(&cinfo, outfile);
								/* We set up the normal JPEG error routines, then override error_exit. */
	cinfo.err = jpeg_std_error(&pub);
	cinfo.image_height = h;
	cinfo.image_width = w;
	cinfo.input_components = ch == 4 ? 3 : ch;

	if (ch == 1) {
		cinfo.in_color_space = J_COLOR_SPACE::JCS_GRAYSCALE;
	}
	else if (ch == 3 || ch == 4) {
		cinfo.in_color_space = J_COLOR_SPACE::JCS_RGB;
	}

	jpeg_set_defaults(&cinfo);
	jpeg_set_quality(&cinfo, 60, TRUE);
	BYTE *rowData = (BYTE *)malloc(sizeof(BYTE) * w * cinfo.input_components);

	jpeg_start_compress(&cinfo, TRUE);
	row_stride = w * ch;

	bool succ = true;


	while (cinfo.next_scanline < h) {
		BYTE* buffer = ((BYTE *)data) + row_stride * cinfo.next_scanline;
		if (ch == 4 || ch == 3) {
			//shrink to 3 channels
			for (int j = 0; j < w; ++j) {
				rowData[j * 3] = buffer[j * ch + 2];
				rowData[j * 3 + 1] = buffer[j * ch + 1];
				rowData[j * 3 + 2] = buffer[j * ch + 0];
			}
		}
		else {
			memcpy(rowData, buffer, sizeof(BYTE) * row_stride);
		}
		int res = jpeg_write_scanlines(&cinfo, &rowData, 1);
		if (res != 1) {
			succ = false;
			break;
		};
	};
	//clean-up all data
	free(rowData);
	jpeg_finish_compress(&cinfo);
	jpeg_destroy_compress(&cinfo);

	void *ret = malloc(sizeof(unsigned char) * *output_length);
	memcpy(ret, tmpBuffer, sizeof(unsigned char) * *output_length);

	free(tmpBuffer);

	return ret;
};

void* decode_jpeg_byte_stream_to_image(void *byte_stream, unsigned long len, int *w, int *h, int *ch, MallocBuffer malloc_)
{
	image_obj obj;
	/* This struct contains the JPEG decompression parameters and pointers to
	* working space (which is allocated as needed by the JPEG library).
	*/
	struct jpeg_decompress_struct cinfo;
	/* More stuff */
	JSAMPARRAY buffer;		/* Output row buffer */
	int row_stride;		/* physical row width in output buffer */

						//if ((infile = fopen(file, "rb")) == NULL) {


	struct jpeg_error_mgr pub;	/* "public" fields */

								// is single-thread application needed exit call back?
								//pub.error_exit = my_error_exit;
								/* We set up the normal JPEG error routines, then override error_exit. */
	cinfo.err = jpeg_std_error(&pub);

	/* Now we can initialize the JPEG decompression object. */
	jpeg_create_decompress(&cinfo);

	/* Step 2: specify data source (eg, a file) */
	//jpeg_stdio_src(&cinfo, infile);
	jpeg_mem_src(&cinfo, (const unsigned char *)byte_stream, len);

	/* Step 3: read file parameters with jpeg_read_header() */
	if (jpeg_read_header(&cinfo, TRUE) != JPEG_HEADER_OK) {
		jpeg_destroy_decompress(&cinfo);
		obj.res = -2;
		return NULL;
	};
	/* We can ignore the return value from jpeg_read_header since
	*   (a) suspension is not possible with the stdio data source, and
	*   (b) we passed TRUE to reject a tables-only JPEG file as an error.
	* See libjpeg.txt for more info.
	*/

	/* Step 4: set parameters for decompression */

	/* In this example, we don't need to change any of the defaults set by
	* jpeg_read_header(), so we do nothing here.
	*/

	/* Step 5: Start decompressor */

	if (!jpeg_start_decompress(&cinfo))
	{
		jpeg_destroy_decompress(&cinfo);
		obj.res = -3;
		return NULL;
	}
	/* We can ignore the return value since suspension is not possible
	* with the stdio data source.
	*/

	/* We may need to do some setup of our own at this point before reading
	* the data.  After jpeg_start_decompress() we have the correct scaled
	* output image dimensions available, as well as the output colormap
	* if we asked for color quantization.
	* In this example, we need to make an output work buffer of the right size.
	*/
	/* JSAMPLEs per row in output buffer */
	row_stride = cinfo.output_width * cinfo.output_components;



	//HCLIBIMAGE ret = CreateImage(cinfo.output_width, cinfo.output_height, cinfo.output_components, "BYTE");
	void *data = malloc_(cinfo.output_components * cinfo.output_width * cinfo.output_height);
	obj.ch = cinfo.output_components;
	obj.w = cinfo.output_width;
	obj.h = cinfo.output_height;
	obj.data = data;

	/* Make a one-row-high sample array that will go away when done with image */
	buffer = (*cinfo.mem->alloc_sarray)
		((j_common_ptr)&cinfo, JPOOL_IMAGE, row_stride, 1);

	bool readingError = false;
	if (cinfo.output_components == 1 ||
		cinfo.output_components == 3)
	{
		/* Step 6: while (scan lines remain to be read) */
		/*           jpeg_read_scanlines(...); */

		BYTE* dataPtr = (BYTE *)data;
		/* Here we use the library's state variable cinfo.output_scanline as the
		* loop counter, so that we don't have to keep track ourselves.
		*/
		while (cinfo.output_scanline < cinfo.output_height) {
			/* jpeg_read_scanlines expects an array of pointers to scanlines.
			* Here the array is only one element long, but you could ask for
			* more than one scanline at a time if that's more convenient.
			*/
			int res = jpeg_read_scanlines(&cinfo, buffer, 1);
			JSAMPROW bb = buffer[0];
			if (res != 1) {

				readingError = true;
				break;
			}
			/* Assume put_scanline_someplace wants a pointer and sample count. */
			//put_scanline_someplace(buffer[0], row_stride);


			if (cinfo.output_components == 1)
			{
				memcpy(dataPtr, bb, row_stride);
			}
			else if (cinfo.output_components == 3)
			{
				int j = 0;
				BYTE *srcPtr = bb;
				while (j < obj.w)
				{
					dataPtr[j * 3] = srcPtr[j * 3 + 2];
					dataPtr[j * 3 + 1] = srcPtr[j * 3 + 1];
					dataPtr[j * 3 + 2] = srcPtr[j * 3];
					++j;
				}
				//memcpy( dataPtr , bb , row_stride) ;
			}
			dataPtr += (obj.w * obj.ch);
		};
	}
	else
	{
		readingError = true;
		obj.res = -5;
		//CWINDOW_DLL_LAST_ERROR = CIMAGE_JPEG_NOT_SUPPORT_CHANNEL;
	}
	/* Step 7: Finish decompression */

	jpeg_finish_decompress(&cinfo);
	/* We can ignore the return value since suspension is not possible
	* with the stdio data source.
	*/

	/* Step 8: Release JPEG decompression object */

	/* This is an important step since it will release a good deal of memory. */
	jpeg_destroy_decompress(&cinfo);


	if (readingError)
	{
		obj.res = -4;
		free(data);
	}
	else
	{
		obj.res = 0;
	}
	return data;
}

void* read_image_wstring(const wchar_t* fileName, int* w, int* h, int* ch, MallocBuffer malloc_)
{
	image_obj obj = image_obj();

	{

		int type = read_image_file_header_wstring(fileName);
		switch (type)
		{
		case 0:
			obj = ReadJpegFile(fileName, malloc_);
			break;
		case 1:
			obj = ReadPngFile(fileName, malloc_);
			break;
		case 2:
			obj = ReadBmpFile(fileName, malloc_);
			break;
		case 3:
			obj = ReadTiffFile(fileName, malloc_);
			break;
		default:
			break;
		}
	}
	*w = obj.w;
	*h = obj.h;
	*ch = obj.ch;

	return obj.data;
};
