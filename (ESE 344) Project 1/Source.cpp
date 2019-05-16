// Jonathan Cruzate
// jcruzate
// 11030802
// ESE 344
// Project 1

#ifndef MATRIX_H
#define MATRIX_H
#define _USE_MATH_DEFINES
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <complex>
#include <math.h>
#include <cassert>

using namespace std;

template <typename Object>
class matrix
{
public:
	matrix(int rows, int cols) : array(rows)
	{
		for (auto & thisRow : array)
			thisRow.resize(cols);
	}

	matrix(vector<vector<Object>> v) : array{ v }
	{ }
	matrix(vector<vector<Object>> && v) : array{ std::move(v) }
	{ }

	const vector<Object> & operator[](int row) const
	{
		return array[row];
	}
	vector<Object> & operator[](int row)
	{
		return array[row];
	}

	int numrows() const
	{
		return array.size();
	}
	int numcols() const
	{
		return numrows() ? array[0].size() : 0;
	}

	void matprt()
	{
		int nr = numrows();
		int nc = numcols();

		for (int i = 0; i < nr; ++i)
		{
			for (int j = 0; j < nc; ++j)
			{
				cout << array[i][j] << "   ";
			}
			cout << endl;
		}
		cout << endl;
	}

	void matprtReal()
	{
		int nr = numrows();
		int nc = numcols();

		for (int i = 0; i < nr; ++i)
		{
			for (int j = 0; j < nc; ++j)
			{
				cout << array[i][j].real() << "   ";
			}
			cout << endl;
		}
		cout << endl;
	}
private:
	vector<vector<Object>> array;
};
#endif

/**
* Standard matrix multiplication.
* Arrays start at 0.
*/
matrix<complex<double>> operator*(const matrix<complex<double>> & a, const matrix<complex<double>> & b)
{
	int nra = a.numrows();
	int nca = a.numcols();
	int nrb = b.numrows();
	int ncb = b.numcols();

	if (nca != nrb) {
		cerr << "nca != nrb in matrix multiplication. Exiting." << endl;
		cout << endl << "Enter any char to exit : " << endl;
		char c;
		cin >> c;
		exit(0);
	}
	matrix<complex<double>> c{ nra, ncb };

	for (int i = 0; i < nra; ++i)
	{
		for (int j = 0; j < ncb; ++j)
		{
			c[i][j] = 0;

			for (int k = 0; k < nca; ++k)
			{
				c[i][j] += a[i][k] * b[k][j];
			}
		}
	}

	return c;
}

/**
* Compute Discrete Fourier Transform of a square matrix
*/
matrix<complex<double>> computeDFT(const int & n)
{

	matrix<complex<double>> a(n, n);
	complex<double> r, s;

	for (int i = 0; i < n; ++i)
	{
		{
			for (int j = 0; j < n; ++j)
			{
				r = complex<double>(cos(-2 * M_PI * i * j / n), sin(-2 * M_PI * i * j / n));
				s = n;
				a[i][j] = r / s;
			}
		}
	}

	return a;
}

/**
* Compute Inverse Discrete Fourier Transform of a square matrix
*/
matrix<complex<double>> computeIDFT(const int & n)
{

	matrix<complex<double>> a(n, n);
	complex<double> r;

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			r = complex<double>(cos(2 * M_PI * i * j / n), sin(2 * M_PI * i * j / n));
			a[i][j] = r;
		}
	}

	return a;
}

void matprint(const matrix<double> & a)
{
	int nr = a.numrows();
	int nc = a.numcols();

	for (int i = 0; i < nr; ++i)
	{
		for (int j = 0; j < nc; ++j)
		{
			cout << a[i][j] << "   ";
		}
		cout << endl;
	}
	cout << endl;

}

int main()
{
	int M, N;

	cin >> M;

	cin >> N;

	if ((M < 1) || (N < 1))
	{
		cout << "M and N cannot be less than 1" << endl;

		assert((1 <= M) && (1 <= N));

		return 1;

		exit(1);
	}

	matrix<double> f(M, N), h(M, N), g_prime(M, N);

	for (int i = 0; i < M; ++i) // Initialization
	{
		for (int j = 0; j < N; ++j)
		{
			cin >> f[i][j];
		}
	}

	cout << "f : " << endl;
	matprint(f);

	for (int i = 0; i < M; ++i) // Initialization
	{
		for (int j = 0; j < N; ++j)
		{
			cin >> h[i][j];
		}
	}

	cout << "h : " << endl;
	matprint(h);

	matrix<complex <double>> f_complex(M, N), F(M, N), f_prime(M, N), G(M, N), g(M, N), h_complex(M, N), H(M, N), h_prime(M, N), P(M, M), P_prime(M, M), Q(N, N), Q_prime(N, N);

	for (int i = 0; i < M; ++i) // Initialization
	{
		for (int j = 0; j < N; ++j)
		{
			f_complex[i][j] = complex<double>(f[i][j], 0.0);
			h_complex[i][j] = complex<double>(h[i][j], 0.0);
		}
	}

	P = computeDFT(M);
	Q = computeDFT(N);

	P_prime = computeIDFT(M);
	Q_prime = computeIDFT(N);

	F = P * f_complex;
	F = F * Q;

	cout << "F : " << endl;
	F.matprt();

	f_prime = P_prime * F;
	f_prime = f_prime * Q_prime;

	cout << "f' : " << endl;
	f_prime.matprtReal();

	H = P * h_complex;
	H = H * Q;

	cout << "H : " << endl;
	H.matprt();

	h_prime = P_prime * H;
	h_prime = h_prime * Q_prime;

	cout << "h' : " << endl;
	h_prime.matprtReal();

	for (int i = 0; i < M; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			G[i][j] = H[i][j] * F[i][j];
		}
	}

	cout << "G : " << endl;
	G.matprt();

	g = P_prime * G;
	g = g * Q_prime;

	cout << "g : " << endl;
	g.matprtReal();

	for (int m = 0; m < M; ++m)
	{
		for (int n = 0; n < N; ++n)
		{
			for (int u = 0; u < M; ++u)
			{
				for (int v = 0; v < N; ++v)
				{
					g_prime[m][n] += (f[(m - u + M) % M][(n - v + N) % N] * h[u][v]) / (M * N);
				}
			}
		}
	}

	cout << "g' : " << endl;
	matprint(g_prime);

	double SSD_error = 0;

	for (int i = 0; i < M; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			SSD_error += pow((g[i][j].real() - g_prime[i][j]), 2);
		}
	}

	cout << "SSD error : " << SSD_error << endl;

	return 0;
}