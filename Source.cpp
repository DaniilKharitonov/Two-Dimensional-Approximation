#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double Real_Function(double, double);
int Reflection_Method(int, double*, double*, double*);
double N_Fun_Polinom(int, double*, double*, double);
void Scanf(FILE *, int*, double*, double*, double*, double, double, double, double, int);
void Estimation(int, int*, double*, double*, double*, double*, double*);
void Estimation_RSA(int, int*, double*, double*, double*, double*);
void Pade_Method(int, double *, double *, double *);
double Pade_Function(int, double *, double);


double Real_Function(double x, double y) {
	return exp(-pow(x - 0.5, 2)-pow(y - 0.5, 2))-exp(-pow(x, 2) - pow(y, 2));
	//return 1 / (1 + 25 * pow(x + y * y, 2));
	//return 4.;
}

int Reflection_Method(int N, double* matrix, double* y, double* a) {

	int i;
	int j;
	int k;
	double tmp1;
	double tmp2;

	for (i = 0; i < N; ++i) {
		tmp1 = 0.0;
		for (j = i + 1; j < N; ++j)
			tmp1 += matrix[j * N + i] * matrix[j * N + i];

		tmp2 = sqrt(tmp1 + matrix[i * N + i] * matrix[i * N + i]);

		if (tmp2 < 1e-100)
			return -1;

		matrix[i * N + i] -= tmp2;

		tmp1 = sqrt(tmp1 + matrix[i * N + i] * matrix[i * N + i]);

		if (tmp1 < 1e-100)
		{
			matrix[i * N + i] += tmp2;
			continue;
		}

		tmp1 = 1.0 / tmp1;
		for (j = i; j < N; ++j)
			matrix[j * N + i] *= tmp1;

		tmp1 = 0.0;
		for (j = i; j < N; ++j)
			tmp1 += matrix[j * N + i] * y[j];

		tmp1 *= 2.0;
		for (j = i; j < N; ++j)
			y[j] -= tmp1 * matrix[j * N + i];

		for (k = i + 1; k < N; ++k)
		{
			tmp1 = 0.0;
			for (j = i; j < N; ++j)
				tmp1 += matrix[j * N + i] * matrix[j * N + k];

			tmp1 *= 2.0;
			for (j = i; j < N; ++j)
				matrix[j * N + k] -= tmp1 * matrix[j * N + i];
		}

		matrix[i * N + i] = tmp2;
	}

	for (i = N - 1; i >= 0; --i)
	{
		tmp1 = y[i];
		for (j = i + 1; j < N; ++j)
			tmp1 -= matrix[i * N + j] * a[j];
		a[i] = tmp1 / matrix[i * N + i];
	}

	return 0;

}

double N_Fun_Polinom(int N, double* x, double* y, double arg) {
	double phi = 1.;
	double fun = 0.;

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++)
			if (j != i)
				phi *= (arg - x[j]) / (x[i] - x[j]);
		fun += y[i] * phi;
		phi = 1.;
	}

	return fun;
}

void Scanf(FILE *input, int *existence, double *exis_values, double *x, double *y, double x0, double xN, double y0, double yN, int N) {
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++) {
			fscanf_s(input, "%d", &existence[i * N + j]);
			x[i] = x0 + (xN - x0) / (N - 1) * i;
			y[j] = y0 + (yN - y0) / (N - 1) * j;
			if (existence[i * N + j] == 1)
				exis_values[i * N + j] = Real_Function(x[i], y[j]);
		}
}

void Pade_Method(int N, double *x, double *values, double *vector) {
	double *matrix, *right;
	int tmp;

	right = (double*)malloc(N * sizeof(double));
	matrix = (double*)malloc(N * N * sizeof(double));

	for (int i = 0; i < N; i++) {
		right[i] = values[i];
		for (int j = 0; j < N; j++) {
			if (j < int(N / 2) + 1) {
				matrix[i * N + j] = pow(x[i], j);
			}
			else {
				matrix[i * N + j] = - values[i]*pow(x[i], j - int(N / 2));
			}
			printf_s("%f ", matrix[i * N + j]);
		}
		printf("\n");
	}

	tmp = Reflection_Method(N, matrix, right, vector);
	
	printf("\n\n");

	for (int i = 0; i < N; i++)
		printf_s("%f ", vector[i]);

	printf_s("\n");

	free(matrix);
	free(right);
}

double Pade_Function(int N, double *b, double x) {
	double num, denum;

	num = 0;
	denum = 1;

	for (int i = 0; i < int(N / 2) + 1; i++) {
		num += b[i] * pow(x, i);
	}


	for (int i = 0; i < N - int(N / 2) - 1; i++)
		denum += b[i + int(N / 2) + 1] * pow(x, i + 1);

	if (N == 2) {
		printf_s("Values = %f %f %f\n", x, num, denum);
	}

	return num / denum;
}

void Estimation(int N, int* existence, double* exis_values, double* x, double* y, double* estim_values1, double* estim_values2) {
	int count;
	double *x_i, *values, *y_j, *b;

	count = 0;
	x_i = (double*)malloc(N * sizeof(double));
	values = (double*)malloc(N * sizeof(double));
	y_j = (double*)malloc(N * sizeof(double));
	b = (double*)malloc(N * sizeof(double)); // a_0, ... , a_(int(N/2) - 1); b_0 = 1, b_1, ... , b_(N - int(N/2))
	
	// First method
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (existence[i * N + j] == 1) {
				x_i[count] = x[j];
				values[count] = exis_values[i * N + j];
				count += 1;
			}
		}

		for (int j = 0; j < N; j++) {
			if (existence[i * N + j] == 0) {
				estim_values1[i * N + j] = N_Fun_Polinom(count-1, x_i, values, x[j]);
			}
		}

		count = 0;
	}

	//Second Method Pade
	for (int j = 0; j < N; j++) {
		for (int i = 0; i < N; i++) {
			printf_s("%f \n", y[i]);
			if (existence[i * N + j] == 1) {
				y_j[count] = y[i];
				values[count] = exis_values[i * N + j];\
				count += 1;
			}
		}
		Pade_Method(count, y_j, values, b);

		for (int i = 0; i < N; i++) {
			if (existence[i * N + j] == 0) {
				estim_values2[i * N + j] = Pade_Function(count, b, y[i]);
			}
		}

		count = 0;
	}

	free(x_i);
	free(y_j);
	free(b);
	free(values);
}

void Estimation_RSA(int N, int* existence, double* exis_values, double* x, double* y, double* estim_values_RSA) {
	double sum_sq, sum;
	int w = 2;
	sum = 0.;
	sum_sq = 0.;

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			
			if (existence[i * N + j] == 0) {
				for (int k = fmax(0, i - w); k < fmin(N, i + w); k++)
					for (int l = fmax(0, j - w); l < fmin(N, j + w); l++)
						if (existence[k * N + l] == 1) {
							sum += exis_values[k * N + l] / (pow(y[i] - y[k], 2) + pow(x[j] - x[l], 2));
							sum_sq += 1/(pow(y[i] - y[k], 2) + pow(x[j] - x[l], 2));
							//printf_s("%f %f\n", sum, sum_sq);
						}
				estim_values_RSA[i * N + j] = sum / sum_sq;
				sum = 0.;
				sum_sq = 0.;
			}

		}
	}
}


int main() {

	double *exis_values, *estim_values1, *estim_values2, *estim_values_RSA, *x, *y;
	double x0, xN, y0, yN, x_i, y_j;
	int N, *existence;
	FILE *input;
	FILE *output_2DA, *output_RSA, *output_Real;
	errno_t err; // Special for VS17

	N = 10;

	err = fopen_s(&input,"input.txt", "r"); // Special for VS17
	err = fopen_s(&output_2DA, "output_2DA.txt", "w"); // Special for VS17
	err = fopen_s(&output_RSA, "output_RSA.txt", "w"); // Special for VS17
	err = fopen_s(&output_Real, "output_Real.txt", "w"); // Special for VS17

	fscanf_s(input, "%d", &N);

	x0 = -1.;
	xN = 1.;
	y0 = -1.;
	yN = 1.;

	existence = (int*)malloc(N * N * sizeof(int)); // 1 - value exist, 0 - doesn't
	exis_values = (double*)malloc(N * N * sizeof(double));
	estim_values1 = (double*)malloc(N * N * sizeof(double)); // for non existing values, estimated by x
	estim_values2 = (double*)malloc(N * N * sizeof(double)); // for non existing values, estimated by y
	estim_values_RSA = (double*)malloc(N * N * sizeof(double));
	x = (double*)malloc(N * sizeof(double));
	y = (double*)malloc(N * sizeof(double));

	Scanf(input, existence, exis_values, x, y, x0, xN, y0, yN, N);

	Estimation(N, existence, exis_values, x, y, estim_values1, estim_values2);
	Estimation_RSA(N, existence, exis_values, x, y, estim_values_RSA);

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (existence[i * N + j] == 1) {
				fprintf(output_2DA, "%f ", exis_values[i * N + j]);
				fprintf(output_RSA, "%f ", exis_values[i * N + j]);
				fprintf(output_Real, "%f ", exis_values[i * N + j]);
			} else {
				fprintf(output_2DA, "%f ", estim_values2[i * N + j]);
				fprintf(output_RSA, "%f ", estim_values_RSA[i * N + j]);
				fprintf(output_Real, "%f ", Real_Function(x[i], y[j]));
			}
		}
		fprintf(output_2DA, "\n");
		fprintf(output_Real, "\n");
		fprintf(output_RSA, "\n");
	}

	system("pause");

	free(existence);
	free(exis_values);
	free(estim_values1);
	free(estim_values2);
	free(x);
	free(y);

	fclose(input);
	fclose(output_Real);
	fclose(output_RSA);
	fclose(output_2DA);

	return 0;
}
