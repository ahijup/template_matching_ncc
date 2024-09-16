#include <string.h>
#include <stdlib.h>
#include "libbmp.h"
#define ABS(A) (A<0?-A:A)
// BMP_HEADER
void
bmp_header_init_df (bmp_header *header,
                    const int   width,
                    const int   height)
{
	header->bfSize = (sizeof (bmp_pixel) * width + BMP_GET_PADDING (width))
	                  * abs (height);
	header->bfReserved = 0;
	header->bfOffBits = 54;
	header->biSize = 40;
	header->biWidth = width;
	header->biHeight = height;
	header->biPlanes = 1;
	header->biBitCount = 24;
	header->biCompression = 0;
	header->biSizeImage = 0;
	header->biXPelsPerMeter = 0;
	header->biYPelsPerMeter = 0;
	header->biClrUsed = 0;
	header->biClrImportant = 0;
}

enum bmp_error
bmp_header_write (const bmp_header *header,
                  FILE             *img_file)
{
	if (header == NULL)
	{
		return BMP_HEADER_NOT_INITIALIZED; 
	}
	else if (img_file == NULL)
	{
		return BMP_FILE_NOT_OPENED;
	}
	
	// Since an adress must be passed to fwrite, create a variable!
	const unsigned short magic = BMP_MAGIC;
	fwrite (&magic, sizeof (magic), 1, img_file);
	
	// Use the type instead of the variable because its a pointer!
	fwrite (header, sizeof (bmp_header), 1, img_file);
	return BMP_OK;
}

enum bmp_error
bmp_header_read (bmp_header *header,
                 FILE       *img_file)
{
	if (img_file == NULL)
	{
		return BMP_FILE_NOT_OPENED;
	}
	
	// Since an adress must be passed to fread, create a variable!
	unsigned short magic;
	
	// Check if its an bmp file by comparing the magic nbr:
	if (fread (&magic, sizeof (magic), 1, img_file) != 1 ||
	    magic != BMP_MAGIC)
	{
		return BMP_INVALID_FILE;
	}
	
	if (fread (header, sizeof (bmp_header), 1, img_file) != 1)
	{
		return BMP_ERROR;
	}

	return BMP_OK;
}

// BMP_PIXEL

void
bmp_pixel_init (bmp_pixel           *pxl,
                const unsigned char  red,
                const unsigned char  green,
                const unsigned char  blue)
{
	pxl->red = red;
	pxl->green = green;
	pxl->blue = blue;
}

// BMP_IMG

void
bmp_img_alloc (bmp_img *img)
{
	const size_t h = abs (img->img_header.biHeight);
	
	// Allocate the required memory for the pixels:
	img->img_pixels = (bmp_pixel **) malloc (sizeof (bmp_pixel*) * h);
	
	for (size_t y = 0; y < h; y++)
	{
		img->img_pixels[y] = (bmp_pixel *) malloc (sizeof (bmp_pixel) * img->img_header.biWidth);
	}
}

void
bmp_img_init_df (bmp_img   *img,
                 const int  width,
                 const int  height)
{
	// INIT the header with default values:
	bmp_header_init_df (&img->img_header, width, height);
	bmp_img_alloc (img);
}

void
bmp_img_free (bmp_img *img)
{
	const size_t h = abs (img->img_header.biHeight);
	
	for (size_t y = 0; y < h; y++)
	{
		free (img->img_pixels[y]);
	}
	free (img->img_pixels);
}

enum bmp_error
bmp_img_write (const bmp_img *img,
               const char    *filename)
{
	FILE *img_file = fopen (filename, "wb");
	
	if (img_file == NULL)
	{
		return BMP_FILE_NOT_OPENED;
	}
	
	// NOTE: This way the correct error code could be returned.
	const enum bmp_error err = bmp_header_write (&img->img_header, img_file);
	
	if (err != BMP_OK)
	{
		// ERROR: Could'nt write the header!
		fclose (img_file);
		return err;
	}
	
	// Select the mode (bottom-up or top-down):
	const size_t h = abs (img->img_header.biHeight);
	const size_t offset = (img->img_header.biHeight > 0 ? h - 1 : 0);
	
	// Create the padding:
	const unsigned char padding[3] = {'\0', '\0', '\0'};
	
	// Write the content:
	for (size_t y = 0; y < h; y++)
	{
		// Write a whole row of pixels to the file:
		fwrite (img->img_pixels[ABS (offset - y)], sizeof (bmp_pixel), img->img_header.biWidth, img_file);
		
		// Write the padding for the row!
		fwrite (padding, sizeof (unsigned char), BMP_GET_PADDING (img->img_header.biWidth), img_file);
	}
	
	// NOTE: All good!
	fclose (img_file);
	return BMP_OK;
}

enum bmp_error
bmp_img_read (bmp_img    *img,
              const char *filename)
{
	FILE *img_file = fopen (filename, "rb");
	
	if (img_file == NULL)
	{
		return BMP_FILE_NOT_OPENED;
	}
	
	// NOTE: This way the correct error code can be returned.
	const enum bmp_error err = bmp_header_read (&img->img_header, img_file);
	
	if (err != BMP_OK)
	{
		// ERROR: Could'nt read the image header!
		fclose (img_file);
		return err;
	}
	
	bmp_img_alloc (img);
	
	// Select the mode (bottom-up or top-down):
	const size_t h = abs (img->img_header.biHeight);
	const size_t offset = (img->img_header.biHeight > 0 ? h - 1 : 0);
	const size_t padding = BMP_GET_PADDING (img->img_header.biWidth);
	
	// Needed to compare the return value of fread
	const size_t items = img->img_header.biWidth;
	
	// Read the content:
	for (size_t y = 0; y < h; y++)
	{
		// Read a whole row of pixels from the file:
		if (fread (img->img_pixels[ABS (offset - y)], sizeof (bmp_pixel), items, img_file) != items)
		{
			fclose (img_file);
			return BMP_ERROR;
		}
		
		// Skip the padding:
		fseek (img_file, (long)padding, SEEK_CUR);
	}
	
	// NOTE: All good!
	fclose (img_file);
	return BMP_OK;
}


bmp_error	bmp_img_read2 (bmp_img_2* img,
							const char *filename) 
{
	FILE *img_file = fopen (filename, "rb");
	
	if (img_file == NULL)
	{
		return BMP_FILE_NOT_OPENED;
	}

	bmp_error err = bmp_header_read(&img->img_header , img_file) ;
	if ( err != BMP_OK ) {
		fclose(img_file) ;
		return err;
	}

	img->pixels = malloc(sizeof(unsigned char) * ( img->img_header.biBitCount / 8 ) * img->img_header.biHeight * img->img_header.biWidth) ;
	// Select the mode (bottom-up or top-down):
	const size_t h = abs (img->img_header.biHeight);
	const size_t offset = (img->img_header.biHeight > 0 ? h - 1 : 0);
	const size_t padding = BMP_GET_PADDING (img->img_header.biWidth);

	// Needed to compare the return value of fread
	const size_t items = img->img_header.biWidth * ( img->img_header.biBitCount / 8 ) ;

	size_t readSize = 0 ;
	// Read the content:
	for (size_t y = 0; y < h; y++)
	{
		void *px = ((unsigned char *)img->pixels) + (offset - y) * items;
		readSize = fread (px, sizeof (unsigned char), items, img_file);

		if ( readSize != items ) {
			fclose (img_file);
			return BMP_ERROR;
		}
		// Skip the padding:
		fseek (img_file, (long)padding, SEEK_CUR);
	}

	
	// NOTE: All good!
	fclose (img_file);
	return BMP_OK ;
};

bmp_error bmp_img_read_ext(const wchar_t *filename, bmp_header &img_header , void *& payload , ImagePointerCreationCallBack payloadCreation) {
	FILE *img_file = _wfopen(filename, L"rb");
	
	if (img_file == NULL)
	{
		return BMP_FILE_NOT_OPENED;
	}

	bmp_error err = bmp_header_read(&img_header , img_file) ;
	if ( err != BMP_OK ) {
		fclose(img_file) ;
		return err;
	}


	unsigned char *pale = NULL ;
	int palettes = 0 ;
	//skipping offset
	if ( img_header.bfOffBits > 0 ) {
		palettes = img_header.biClrUsed - img_header.biClrImportant;
		if ( palettes > 0 ) {
			pale = (unsigned char *)malloc(sizeof(unsigned char) * palettes * 4) ;
			fread ( pale , sizeof(unsigned char) , palettes * 4, img_file ) ;
		} else {
			fseek (img_file, (long)(img_header.bfOffBits - 54), SEEK_CUR);
		}
	}
	//char *offset1 = new char[img_header.bfOffBits - 54] ;
	//fread (offset1, sizeof (unsigned char), img_header.bfOffBits - 54, img_file);
	//delete offset1 ;

	// Select the mode (bottom-up or top-down):
	const size_t h = abs (img_header.biHeight);
	const size_t offset = (img_header.biHeight > 0 ? h - 1 : 0);
	int fnChannel = img_header.biBitCount / 8;

	const size_t padding = BMP_GET_PADDING (img_header.biWidth * fnChannel);

	
	size_t bufferSz = sizeof(unsigned char) * fnChannel * h * (img_header.biWidth);
	
	unsigned char* buffer = (unsigned char *)malloc(bufferSz) ;
	int channels = img_header.biBitCount / 8 ;

	if ( pale ) {
		channels = 3 ;
	}

	size_t sz = sizeof(unsigned char) * channels * h * (img_header.biWidth);
	
	payloadCreation(payload , sz) ;

	// Needed to compare the return value of fread(each row)
	const size_t items = ((img_header.biWidth) * fnChannel )  ;

	unsigned char *ptr = ((unsigned char *)payload) ;
	size_t readSize;

	if ( img_header.biHeight > 0 ) {
		for ( int y = (int)h - 1 ; y >= 0 ; --y ) {
			readSize = fread(buffer , sizeof(unsigned char) , items , img_file) ;
			if ( readSize != items ) {
				fclose(img_file) ;
				free(buffer) ;
				return BMP_ERROR;
			}

			//skipping padding bytes
			fseek (img_file, (long)padding , SEEK_CUR);


			unsigned char *curPtr = ptr + ( channels * y * img_header.biWidth ) ;
			//readSize = fread (curPtr, sizeof (unsigned char), items, img_file) ;
			//trans pixel to palette color pixel
			if ( palettes ) {
				for ( int i = 0 ; i < img_header.biWidth ; ++i , curPtr += channels ) {
					unsigned char idx = buffer[i] ;
					unsigned char *indexed = pale + idx * 4;
					memcpy(curPtr , indexed , 3) ;
				}
			} else {
				memcpy(curPtr , buffer , items) ;
			}
		}
	} else {
		for ( int y = 0 ; y < (int)h ; ++y) {
			readSize = fread(buffer , sizeof(unsigned char) , items , img_file) ;
			if ( readSize != items ) {
				fclose(img_file) ;
				free(buffer) ;
				return BMP_ERROR;
			}

			//skipping padding bytes
			fseek (img_file, (long)padding , SEEK_CUR);


			unsigned char *curPtr = ptr + ( channels * y * img_header.biWidth ) ;
			//readSize = fread (curPtr, sizeof (unsigned char), items, img_file) ;
			//trans pixel to palette color pixel
			if ( palettes ) {
				for ( int i = 0 ; i < img_header.biWidth ; ++i , curPtr += channels ) {
					unsigned char idx = buffer[i] ;
					unsigned char *indexed = pale + idx * 4;
					memcpy(curPtr , indexed , 3) ;
				}
			} else {
				memcpy(curPtr , buffer , items) ;
			}
		}
	}

	free(buffer) ;
	if ( pale ) {
		free(pale) ;
	}

	img_header.biBitCount = channels * 8 ;
	// Read the content:
	/*for (size_t y = 0; y < h; y++)
	{
		void *px = ((unsigned char *)payload) + (offset - y) * items;
		readSize = fread (px, sizeof (unsigned char), items, img_file);

		if ( readSize != items ) {
			fclose (img_file);
			return BMP_ERROR;
		}
		// Skip the padding:
		fseek (img_file, (long)padding, SEEK_CUR);
	}*/

	// NOTE: All good!
	fclose (img_file);
	return BMP_OK ;
} ;

void bmp_img_2_free(bmp_img_2 *img) {
	if(img->pixels) {
		free(img->pixels) ;
		img->pixels = NULL ;
	}
} ;

enum bmp_error
bmp_img_write_2 (const bmp_header &header , void *ptr ,
               const char    *filename)
{
	FILE *img_file = fopen (filename, "wb");
	
	if (img_file == NULL)
	{
		return BMP_FILE_NOT_OPENED;
	}

	// NOTE: This way the correct error code could be returned.
	const enum bmp_error err = bmp_header_write (&header, img_file);
	
	if (err != BMP_OK)
	{
		// ERROR: Could'nt write the header!
		fclose (img_file);
		return err;
	}
	
	// Select the mode (bottom-up or top-down):
	const size_t h = abs (header.biHeight);
	const size_t offset = (header.biHeight > 0 ? h - 1 : 0);
	const int ch = header.biBitCount / 8 ;
	
	if ( (header.bfOffBits - 54) > 0 ) {
		//write down pelette
		for ( int i = 0 ; i < 256 ; ++i ) {
			unsigned char pel[4] = { (unsigned char)i , (unsigned char)i , (unsigned char)i , 0xFF } ;
			fwrite(pel , sizeof(unsigned char) , 4 , img_file) ;
		}
	}
	
	
	int paddings = BMP_GET_PADDING(header.biWidth * ch) ;



	// Create the padding:
	const unsigned char padding[9] = {'\0', '\0', '\0' , '\0' , '\0' , '\0' , 
	'\0' , '\0' , '\0'};
	
	

	if ( header.biHeight > 0 ) {
		// Write the content:
		for (int y = h - 1; y >= 0; --y)
		{
			unsigned char* rowPtr = (unsigned char *)(ptr) + y * ch * header.biWidth;
			fwrite( rowPtr , sizeof(unsigned char) , ch * header.biWidth , img_file) ;
			// Write a whole row of pixels to the file:
			//fwrite (img->img_pixels[ABS (offset - y)], sizeof (bmp_pixel), img->img_header.biWidth, img_file);
			// Write the padding for the row!
			fwrite (padding, sizeof (unsigned char), paddings , img_file);
		}
	} else {
		// Write the content:
		for (size_t y = 0; y < h; y++)
		{
			unsigned char* rowPtr = (unsigned char *)(ptr) + y * ch * header.biWidth;
			fwrite( rowPtr , sizeof(unsigned char) , ch * header.biWidth , img_file) ;
			// Write a whole row of pixels to the file:
			//fwrite (img->img_pixels[ABS (offset - y)], sizeof (bmp_pixel), img->img_header.biWidth, img_file);
			// Write the padding for the row!
			fwrite (padding, sizeof (unsigned char), paddings , img_file);
		}
	}
	// NOTE: All good!
	fclose (img_file);
	return BMP_OK;
};

void bmp_header_init_df_2 (bmp_header* header,
                           const int width,
                           const int height, 
						   const int channels)
{
	int pad = BMP_GET_PADDING (width * channels);
	
	header->bfSize = (sizeof (unsigned char) * (width  * channels + pad) * abs (height)) + 1078;
	header->bfReserved = 0;
	header->bfOffBits = 1078;
	header->biSize = 40;
	header->biWidth = width;
	header->biHeight = height;
	header->biPlanes = 1;
	header->biBitCount = channels * 8;
	header->biCompression = 0;
	header->biSizeImage = (width * channels + pad) * abs( height ) ;
	header->biXPelsPerMeter = 0;
	header->biYPelsPerMeter = 0;
	header->biClrUsed = 256;
	header->biClrImportant = 256;
} ;

enum bmp_error bmp_img_read_ext_2(const wchar_t *filename, bmp_header &img_header, void *& payload, MallocImageBuffer malloc_)
{
	FILE *img_file = _wfopen(filename, L"rb");

	if (img_file == NULL)
	{
		return BMP_FILE_NOT_OPENED;
	}

	bmp_error err = bmp_header_read(&img_header, img_file);
	if (err != BMP_OK) {
		fclose(img_file);
		return err;
	}


	unsigned char *pale = NULL;
	int palettes = 0;
	//skipping offset
	if (img_header.bfOffBits > 0) {
		palettes = img_header.biClrUsed - img_header.biClrImportant;
		if (palettes > 0) {
			pale = (unsigned char *)malloc(sizeof(unsigned char) * palettes * 4);
			fread(pale, sizeof(unsigned char), palettes * 4, img_file);
		}
		else {
			fseek(img_file, (long)(img_header.bfOffBits - 54), SEEK_CUR);
		}
	}
	//char *offset1 = new char[img_header.bfOffBits - 54] ;
	//fread (offset1, sizeof (unsigned char), img_header.bfOffBits - 54, img_file);
	//delete offset1 ;

	// Select the mode (bottom-up or top-down):
	const size_t h = abs(img_header.biHeight);
	const size_t offset = (img_header.biHeight > 0 ? h - 1 : 0);
	int fnChannel = img_header.biBitCount / 8;

	const size_t padding = BMP_GET_PADDING(img_header.biWidth * fnChannel);


	size_t bufferSz = sizeof(unsigned char) * fnChannel * h * (img_header.biWidth);

	unsigned char* buffer = (unsigned char *)malloc(bufferSz);
	int channels = img_header.biBitCount / 8;

	if (pale) {
		channels = 3;
	}

	size_t sz = sizeof(unsigned char) * channels * h * (img_header.biWidth);

	payload = malloc_(sz);
	//payloadCreation(payload, sz);

	// Needed to compare the return value of fread(each row)
	const size_t items = ((img_header.biWidth) * fnChannel);

	unsigned char *ptr = ((unsigned char *)payload);
	size_t readSize;

	if (img_header.biHeight > 0) {
		for (int y = (int)h - 1; y >= 0; --y) {
			readSize = fread(buffer, sizeof(unsigned char), items, img_file);
			if (readSize != items) {
				fclose(img_file);
				free(buffer);
				return BMP_ERROR;
			}

			//skipping padding bytes
			fseek(img_file, (long)padding, SEEK_CUR);


			unsigned char *curPtr = ptr + (channels * y * img_header.biWidth);
			//readSize = fread (curPtr, sizeof (unsigned char), items, img_file) ;
			//trans pixel to palette color pixel
			if (palettes) {
				for (int i = 0; i < img_header.biWidth; ++i, curPtr += channels) {
					unsigned char idx = buffer[i];
					unsigned char *indexed = pale + idx * 4;
					memcpy(curPtr, indexed, 3);
				}
			}
			else {
				memcpy(curPtr, buffer, items);
			}
		}
	}
	else {
		for (int y = 0; y < (int)h; ++y) {
			readSize = fread(buffer, sizeof(unsigned char), items, img_file);
			if (readSize != items) {
				fclose(img_file);
				free(buffer);
				return BMP_ERROR;
			}

			//skipping padding bytes
			fseek(img_file, (long)padding, SEEK_CUR);


			unsigned char *curPtr = ptr + (channels * y * img_header.biWidth);
			//readSize = fread (curPtr, sizeof (unsigned char), items, img_file) ;
			//trans pixel to palette color pixel
			if (palettes) {
				for (int i = 0; i < img_header.biWidth; ++i, curPtr += channels) {
					unsigned char idx = buffer[i];
					unsigned char *indexed = pale + idx * 4;
					memcpy(curPtr, indexed, 3);
				}
			}
			else {
				memcpy(curPtr, buffer, items);
			}
		}
	}

	free(buffer);
	if (pale) {
		free(pale);
	}

	img_header.biBitCount = channels * 8;
	

	// NOTE: All good!
	fclose(img_file);
	return BMP_OK;
}

enum bmp_error bmp_img_read_ext_3(const char *filename, bmp_header &img_header, void *& payload, MallocImageBuffer malloc_)
{
	FILE *img_file = fopen(filename, "rb");

	if (img_file == NULL)
	{
		return BMP_FILE_NOT_OPENED;
	}

	bmp_error err = bmp_header_read(&img_header, img_file);
	if (err != BMP_OK) {
		fclose(img_file);
		return err;
	}


	unsigned char *pale = NULL;
	int palettes = 0;
	//skipping offset
	if (img_header.bfOffBits > 0) {
		palettes = img_header.biClrUsed - img_header.biClrImportant;
		if (palettes > 0) {
			pale = (unsigned char *)malloc(sizeof(unsigned char) * palettes * 4);
			fread(pale, sizeof(unsigned char), palettes * 4, img_file);
		}
		else {
			fseek(img_file, (long)(img_header.bfOffBits - 54), SEEK_CUR);
		}
	}
	//char *offset1 = new char[img_header.bfOffBits - 54] ;
	//fread (offset1, sizeof (unsigned char), img_header.bfOffBits - 54, img_file);
	//delete offset1 ;

	// Select the mode (bottom-up or top-down):
	const size_t h = abs(img_header.biHeight);
	const size_t offset = (img_header.biHeight > 0 ? h - 1 : 0);
	int fnChannel = img_header.biBitCount / 8;

	const size_t padding = BMP_GET_PADDING(img_header.biWidth * fnChannel);


	size_t bufferSz = sizeof(unsigned char) * fnChannel * h * (img_header.biWidth);

	unsigned char* buffer = (unsigned char *)malloc(bufferSz);
	int channels = img_header.biBitCount / 8;

	if (pale) {
		channels = 3;
	}

	size_t sz = sizeof(unsigned char) * channels * h * (img_header.biWidth);

	payload = malloc_(sz);
	//payloadCreation(payload, sz);

	// Needed to compare the return value of fread(each row)
	const size_t items = ((img_header.biWidth) * fnChannel);

	unsigned char *ptr = ((unsigned char *)payload);
	size_t readSize;

	if (img_header.biHeight > 0) {
		for (int y = (int)h - 1; y >= 0; --y) {
			readSize = fread(buffer, sizeof(unsigned char), items, img_file);
			if (readSize != items) {
				fclose(img_file);
				free(buffer);
				return BMP_ERROR;
			}

			//skipping padding bytes
			fseek(img_file, (long)padding, SEEK_CUR);


			unsigned char *curPtr = ptr + (channels * y * img_header.biWidth);
			//readSize = fread (curPtr, sizeof (unsigned char), items, img_file) ;
			//trans pixel to palette color pixel
			if (palettes) {
				for (int i = 0; i < img_header.biWidth; ++i, curPtr += channels) {
					unsigned char idx = buffer[i];
					unsigned char *indexed = pale + idx * 4;
					memcpy(curPtr, indexed, 3);
				}
			}
			else {
				memcpy(curPtr, buffer, items);
			}
		}
	}
	else {
		for (int y = 0; y < (int)h; ++y) {
			readSize = fread(buffer, sizeof(unsigned char), items, img_file);
			if (readSize != items) {
				fclose(img_file);
				free(buffer);
				return BMP_ERROR;
			}

			//skipping padding bytes
			fseek(img_file, (long)padding, SEEK_CUR);


			unsigned char *curPtr = ptr + (channels * y * img_header.biWidth);
			//readSize = fread (curPtr, sizeof (unsigned char), items, img_file) ;
			//trans pixel to palette color pixel
			if (palettes) {
				for (int i = 0; i < img_header.biWidth; ++i, curPtr += channels) {
					unsigned char idx = buffer[i];
					unsigned char *indexed = pale + idx * 4;
					memcpy(curPtr, indexed, 3);
				}
			}
			else {
				memcpy(curPtr, buffer, items);
			}
		}
	}

	free(buffer);
	if (pale) {
		free(pale);
	}

	img_header.biBitCount = channels * 8;


	// NOTE: All good!
	fclose(img_file);
	return BMP_OK;
};

enum bmp_error
	bmp_img_write_3(const bmp_header &header, void *ptr,
		const wchar_t    *filename)
{
	FILE *img_file = _wfopen(filename, L"wb");

	if (img_file == NULL)
	{
		return BMP_FILE_NOT_OPENED;
	}

	// NOTE: This way the correct error code could be returned.
	const enum bmp_error err = bmp_header_write(&header, img_file);

	if (err != BMP_OK)
	{
		// ERROR: Could'nt write the header!
		fclose(img_file);
		return err;
	}

	// Select the mode (bottom-up or top-down):
	const size_t h = abs(header.biHeight);
	const size_t offset = (header.biHeight > 0 ? h - 1 : 0);
	const int ch = header.biBitCount / 8;

	if ((header.bfOffBits - 54) > 0) {
		//write down pelette
		for (int i = 0; i < 256; ++i) {
			unsigned char pel[4] = { (unsigned char)i , (unsigned char)i , (unsigned char)i , 0xFF };
			fwrite(pel, sizeof(unsigned char), 4, img_file);
		}
	}


	int paddings = BMP_GET_PADDING(header.biWidth * ch);



	// Create the padding:
	const unsigned char padding[9] = { '\0', '\0', '\0' , '\0' , '\0' , '\0' ,
		'\0' , '\0' , '\0' };



	if (header.biHeight > 0) {
		// Write the content:
		for (int y = h - 1; y >= 0; --y)
		{
			unsigned char* rowPtr = (unsigned char *)(ptr)+y * ch * header.biWidth;
			fwrite(rowPtr, sizeof(unsigned char), ch * header.biWidth, img_file);
			// Write a whole row of pixels to the file:
			//fwrite (img->img_pixels[ABS (offset - y)], sizeof (bmp_pixel), img->img_header.biWidth, img_file);
			// Write the padding for the row!
			fwrite(padding, sizeof(unsigned char), paddings, img_file);
		}
	}
	else {
		// Write the content:
		for (size_t y = 0; y < h; y++)
		{
			unsigned char* rowPtr = (unsigned char *)(ptr)+y * ch * header.biWidth;
			fwrite(rowPtr, sizeof(unsigned char), ch * header.biWidth, img_file);
			// Write a whole row of pixels to the file:
			//fwrite (img->img_pixels[ABS (offset - y)], sizeof (bmp_pixel), img->img_header.biWidth, img_file);
			// Write the padding for the row!
			fwrite(padding, sizeof(unsigned char), paddings, img_file);
		}
	}
	// NOTE: All good!
	fclose(img_file);
	return BMP_OK;
};
