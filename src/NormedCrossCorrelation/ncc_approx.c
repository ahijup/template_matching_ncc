#include "ncc_approx.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Define the maximum number of basis functions
#define MAX_BASIS_FUNCTIONS 1000

// Define the acceptable approximation error threshold
#define ERROR_THRESHOLD 1e-3

// Structure to represent a basis function (rectangular region with a weight)
typedef struct {
    int x_start;  // Starting x-coordinate
    int x_end;    // Ending x-coordinate
    int y_start;  // Starting y-coordinate
    int y_end;    // Ending y-coordinate
    double weight; // Coefficient ki
} BasisFunction;

// Structure to hold template data and related information
typedef struct {
    int width;          // Width of the template
    int height;         // Height of the template
    double** data;      // 2D array of template pixel values
    double mean;        // Mean of the template
    double** zero_mean; // Zero-mean template data
} Template;

BasisFunction basis_functions[MAX_BASIS_FUNCTIONS];

void computeZeroMeanTemplate(Template* template);
double computeApproximationError(Template* template, BasisFunction* basis_functions, int num_basis_functions);
void computeBasisFunctionCoefficient(Template* template, BasisFunction* bf);
void splitBasisFunction(Template* template, BasisFunction* bf_to_split, BasisFunction* new_bfs, int* num_basis_functions);
void recursiveBasisFunctionDetermination(Template* template, BasisFunction* basis_functions, int* num_basis_functions);
void freeTemplate(Template* template);

// Function to compute the zero-mean template
void computeZeroMeanTemplate(Template* template) {
    double sum = 0.0;
    for (int y = 0; y < template->height; y++) {
        for (int x = 0; x < template->width; x++) {
            sum += template->data[y][x];
        }
    }
    template->mean = sum / (template->width * template->height);

    // Compute zero-mean template
    for (int y = 0; y < template->height; y++) {
        for (int x = 0; x < template->width; x++) {
            template->zero_mean[y][x] = template->data[y][x] - template->mean;
        }
    }
}

// Function to compute the approximation error J
double computeApproximationError(Template* template, BasisFunction* basis_functions, int num_basis_functions) {
    double error = 0.0;

    for (int y = 0; y < template->height; y++) {
        for (int x = 0; x < template->width; x++) {
            // Compute the approximation at (x, y)
            double approx = 0.0;
            for (int i = 0; i < num_basis_functions; i++) {
                BasisFunction* bf = &basis_functions[i];
                if (x >= bf->x_start && x <= bf->x_end && y >= bf->y_start && y <= bf->y_end) {
                    approx += bf->weight;
                }
            }
            // Compute the error at (x, y)
            double diff = template->data[y][x] - approx;
            error += diff * diff;
        }
    }
    return error;
}

// Function to compute the coefficient ki for a given basis function
void computeBasisFunctionCoefficient(Template* template, BasisFunction* bf) {
    double numerator = 0.0;
    double denominator = 0.0;

    for (int y = bf->y_start; y <= bf->y_end; y++) {
        for (int x = bf->x_start; x <= bf->x_end; x++) {
            numerator += template->zero_mean[y][x];
            denominator += 1.0; // Since t_i(x, y) = 1 within the region
        }
    }

    bf->weight = numerator / denominator;
}

// Function to split a basis function into two new basis functions
void splitBasisFunction(Template* template, BasisFunction* bf_to_split, BasisFunction* new_bfs, int* num_basis_functions) {
    // Decide where to split: along the longer dimension
    int x_length = bf_to_split->x_end - bf_to_split->x_start + 1;
    int y_length = bf_to_split->y_end - bf_to_split->y_start + 1;

    BasisFunction bf1, bf2;

    if (x_length >= y_length) {
        // Split along x-axis
        int x_mid = (bf_to_split->x_start + bf_to_split->x_end) / 2;

        // First new basis function
        bf1.x_start = bf_to_split->x_start;
        bf1.x_end = x_mid;
        bf1.y_start = bf_to_split->y_start;
        bf1.y_end = bf_to_split->y_end;

        // Second new basis function
        bf2.x_start = x_mid + 1;
        bf2.x_end = bf_to_split->x_end;
        bf2.y_start = bf_to_split->y_start;
        bf2.y_end = bf_to_split->y_end;
    }
    else {
        // Split along y-axis
        int y_mid = (bf_to_split->y_start + bf_to_split->y_end) / 2;

        // First new basis function
        bf1.x_start = bf_to_split->x_start;
        bf1.x_end = bf_to_split->x_end;
        bf1.y_start = bf_to_split->y_start;
        bf1.y_end = y_mid;

        // Second new basis function
        bf2.x_start = bf_to_split->x_start;
        bf2.x_end = bf_to_split->x_end;
        bf2.y_start = y_mid + 1;
        bf2.y_end = bf_to_split->y_end;
    }

    // Compute coefficients for new basis functions
    computeBasisFunctionCoefficient(template, &bf1);
    computeBasisFunctionCoefficient(template, &bf2);

    // Replace the old basis function with the first new one
    *bf_to_split = bf1;

    // Add the second new basis function to the array
    if (*num_basis_functions < MAX_BASIS_FUNCTIONS) {
        basis_functions[*num_basis_functions] = bf2;
        (*num_basis_functions)++;
    }
    else {
        fprintf(stderr, "Maximum number of basis functions reached.\n");
    }
}

// Recursive function to determine basis functions
void recursiveBasisFunctionDetermination(Template* template, BasisFunction* basis_functions, int* num_basis_functions) {
	if (*num_basis_functions >= MAX_BASIS_FUNCTIONS) {
		fprintf(stderr, "Maximum number of basis functions reached.\n");
		return;
	}
    // Compute the overall approximation error
    double error = computeApproximationError(template, basis_functions, *num_basis_functions);
	printf("iter[%d] => Error: %f\n", *num_basis_functions, error);
    // Check if error is acceptable
    if (error <= ERROR_THRESHOLD) {
        return;
    }

    // Find the basis function with the highest error contribution
    int max_error_index = -1;
    double max_error = -1.0;

    for (int i = 0; i < *num_basis_functions; i++) {
        BasisFunction* bf = &basis_functions[i];

        // Compute the error contribution of this basis function
        double bf_error = 0.0;

        for (int y = bf->y_start; y <= bf->y_end; y++) {
            for (int x = bf->x_start; x <= bf->x_end; x++) {
                // Compute the approximation at (x, y) without this basis function
                double approx = 0.0;
                for (int j = 0; j < *num_basis_functions; j++) {
                    if (j == i) continue; // Skip current basis function
                    BasisFunction* bf_other = &basis_functions[j];
                    if (x >= bf_other->x_start && x <= bf_other->x_end &&
                        y >= bf_other->y_start && y <= bf_other->y_end) {
                        approx += bf_other->weight;
                    }
                }
                // Compute the error at (x, y)
                double diff = template->zero_mean[y][x] - approx - bf->weight;
                bf_error += diff * diff;
            }
        }

        if (bf_error > max_error) {
            max_error = bf_error;
            max_error_index = i;
        }
    }

    // If no basis function contributes to error, return
    if (max_error_index == -1) {
        return;
    }

    // Split the basis function with the highest error contribution
    BasisFunction bf_to_split = basis_functions[max_error_index];
    splitBasisFunction(template, &basis_functions[max_error_index], basis_functions, num_basis_functions);

    // Recurse
    recursiveBasisFunctionDetermination(template, basis_functions, num_basis_functions);
}

// Function to free allocated memory for the template
void freeTemplate(Template* template) {
    if (template != NULL) {
        for (int y = 0; y < template->height; y++) {
            free(template->data[y]);
            free(template->zero_mean[y]);
        }
        free(template->data);
        free(template->zero_mean);
        free(template);
    }
}

Template templ;

void ncc_approx_offline(struct ImageData* pattern)
{
	// Initialize the template
	templ.width = pattern->w;
	templ.height = pattern->h;
	templ.data = (double**)malloc(templ.height * sizeof(double*));
	templ.zero_mean = (double**)malloc(templ.height * sizeof(double*));
	for (int y = 0; y < templ.height; y++) {
		templ.data[y] = (double*)malloc(templ.width * sizeof(double));
		templ.zero_mean[y] = (double*)malloc(templ.width * sizeof(double));
		for (int x = 0; x < templ.width; x++) {
			templ.data[y][x] = pattern->data[y * pattern->w + x];
		}
	}

	// Compute the zero-mean template
	computeZeroMeanTemplate(&templ);

	// Initialize the basis functions
	BasisFunction bf;
	bf.x_start = 0;
	bf.x_end = templ.width - 1;
	bf.y_start = 0;
	bf.y_end = templ.height - 1;
	computeBasisFunctionCoefficient(&templ, &bf);
	basis_functions[0] = bf;

	int num_basis_functions = 1;

	// Determine basis functions recursively
	recursiveBasisFunctionDetermination(&templ, basis_functions, &num_basis_functions);

	// Free the template
	// freeTemplate(&templ);
}

void ncc_approx_online(struct ImageData* target, struct ImageData* result)
{

}
