#include <stdio.h>
#include <fstream>
#include <stdlib.h>

#pragma once

using namespace std;

#ifndef BEGIN_FOR_MATRIX_LOOP_ROW_RI_COLUMN_CI
#define BEGIN_FOR_MATRIX_LOOP_ROW_RI_COLUMN_CI(N, M) \
	for (unsigned int ri = 0; ri < N; ri++) { \
		for (unsigned int ci = 0; ci < M; ci++){
#endif 

#ifndef END_FOR_MATRIX_LOOP_ROW_RI_COLUMN_CI 
#define END_FOR_MATRIX_LOOP_ROW_RI_COLUMN_CI \
	} } 
#endif 

#ifndef BEGIN_FOR_MATRIX_LOOP_RI_INDEX
#define BEGIN_FOR_MATRIX_LOOP_RI_INDEX(N) \
	for (unsigned int ri = 0; ri < N; ri++) \
	{ 
#endif 

#ifndef END_FOR_MATRIX_LOOP_RI_INDEX
#define END_FOR_MATRIX_LOOP_RI_INDEX \
				}
#endif 

#ifndef CLONE_MATRIX
#define CLONE_MATRIX (m1, m2, N, M) \
	for (unsigned int ri = 0; ri < N; ri++) \
	{ \
		for (unsigned int ci = 0; ci < M; ci++) \
		{ \
			m1(ri, ci) = m2(ri, ci); \
		} \
	}
#endif

	template<typename T, unsigned int N, unsigned int M>
	class Matrix
	{
	public:
		struct Array { typedef T Type[N][M]; };
		typedef Matrix* MatrixPointer;
	public:
		typename Matrix::Array::Type data;
	public:
		Matrix()
		{
			*this = { {} };
		}

		Matrix(Matrix& value) noexcept
		{
			*this = value;
		}

		Matrix(Matrix&& value) noexcept
		{
			*this = value;
		}

		Matrix(typename Matrix::Array::Type&& value)
		{
			*this = value;
		}

		Matrix(typename Matrix::Array::Type& value)
		{
			*this = value;
		}

	private:
		
		T Det1x1() { _ASSERT(N == M && N == 1); return data[0][0]; }
		T Det2x2()
		{
			_ASSERT(N == M && N == 2);
			return data[0][0] * data[1][1] - data[1][0] * data[0][1];
		}

	public:
		T& operator()(unsigned int row, unsigned int column)
		{
			return this->Data(row, column);
		}

		Matrix operator +(const Matrix& value)
		{
			Matrix result(*this);
			result += value;
			return result;
		}

		Matrix& Add(const Matrix& value)
		{
			BEGIN_FOR_MATRIX_LOOP_ROW_RI_COLUMN_CI(N, M);
			this->data[ri][ci] += value.data[ri][ci];
			END_FOR_MATRIX_LOOP_ROW_RI_COLUMN_CI;

			return *this;
		}

		Matrix& operator +=(const Matrix& value)
		{
			this->Add(value);
			return *this;
		}

		Matrix operator -(const Matrix& value)
		{
			Matrix result(*this);
			result -= value;
			return result;
		}

		Matrix& operator -=(const Matrix& value)
		{
			BEGIN_FOR_MATRIX_LOOP_ROW_RI_COLUMN_CI(N, M);
			this->data[ri][ci] -= value.data[ri][ci];
			END_FOR_MATRIX_LOOP_ROW_RI_COLUMN_CI;

			return *this;
		}

		Matrix& operator /(const T value)
		{
			BEGIN_FOR_MATRIX_LOOP_ROW_RI_COLUMN_CI(N, M);
			this->data[ri][ci] = this->data[ri][ci] / value;
			END_FOR_MATRIX_LOOP_ROW_RI_COLUMN_CI;

			return *this;
		}

		template<unsigned int P>
		Matrix<T, N, P> operator *(Matrix<T, M, P>& value)
		{
			Matrix<T, N, P> result({});
			BEGIN_FOR_MATRIX_LOOP_ROW_RI_COLUMN_CI(N, P);
			for (unsigned int x = 0; x < M; x++)
			{
				result(ri, ci) += this->data[ri][x] * value(x, ci);
			}
			END_FOR_MATRIX_LOOP_ROW_RI_COLUMN_CI;

			return result;
		}

		Matrix& operator =(Matrix& value)
		{
			BEGIN_FOR_MATRIX_LOOP_ROW_RI_COLUMN_CI(N, M)
				this->data[ri][ci] = value(ri, ci);
			END_FOR_MATRIX_LOOP_ROW_RI_COLUMN_CI;

			return *this;
		}

		Matrix& operator =(typename Matrix::Array::Type& value)
		{
			BEGIN_FOR_MATRIX_LOOP_ROW_RI_COLUMN_CI(N, M)
				data[ri][ci] = value[ri][ci];
			END_FOR_MATRIX_LOOP_ROW_RI_COLUMN_CI;

			return *this;
		}

		Matrix& operator =(typename Matrix::Array::Type&& value)
		{
			BEGIN_FOR_MATRIX_LOOP_ROW_RI_COLUMN_CI(N, M)
				data[ri][ci] = value[ri][ci];
			END_FOR_MATRIX_LOOP_ROW_RI_COLUMN_CI;

			return *this;
		}

		Matrix<T, M, N> Traspose()
		{
			Matrix<T, M, N> result({});

			BEGIN_FOR_MATRIX_LOOP_ROW_RI_COLUMN_CI(N, M)
				result(ci, ri) = this->data[ri][ci];
			END_FOR_MATRIX_LOOP_ROW_RI_COLUMN_CI;

			return result;
		}

		T& Data(unsigned int row, unsigned int column)
		{
			_ASSERT((row < N&& row >= 0) && (column < M&& column >= 0));
			return data[row][column];
		}

		Matrix Inverse()
		{
			_ASSERT(N == M && N > 1);
			T det = Determinant();
			_ASSERT(det != 0);

			return Adjoint() / det;
		}

		Matrix Adjoint()
		{
			_ASSERT(N == M && N > 1);
			int sign = 0;
			T det = 0;
			unsigned int n = 0, m = 0, p = 0, k = 0;
			Matrix matrix({});
			Matrix<T, M - 1, M - 1>* mxs = new Matrix<T, M - 1, M - 1>();
			for (unsigned int fi = 0; fi < M; fi++)
			{
				for (unsigned int fj = 0; fj < M; fj++)
				{
					k = 0;
					n = 0;
					while (n < M)
					{
						if (n != fi)
						{
							m = 0;
							p = 0;
							while (m < M)
							{
								if (m != fj)
								{
									mxs->Data(k, p) = data[n][m];
									p++;
								}
								m++;
							}
							k++;
						}
						n++;
					}

					det = mxs->Determinant();
					if (det != 0)
					{
						sign = ((fi + fj) & 1) == 0 ? 1 : -1;
						matrix.Data(fi, fj) = sign * det;
					}
				}
			}
			delete mxs;
			
			return matrix.Traspose();
		}

		T Determinant()
		{
			_ASSERT(N == M && N > 0);
			if (N > 2)
			{
				T det = 0;
				T tdet = 0;
				unsigned int s = 0, k = 0, p = 0, m = 0;
				int sign = 0;
				Matrix<T, M - 1, M - 1>* mxs = new Matrix<T, M - 1, M - 1>();
				for (unsigned int j = 0; j < M; j++)
				{
					if (data[0][j] != 0)
					{
						k = 0;
						for(unsigned int n = 1; n < M; n++)
						{
							p = 0;
							m = 0;
							while (m < M)
							{
								if (j != m)
								{
									mxs->Data(k, p) = data[n][m];
									p++;
								}
								m++;
							}
							k++;
						}

						tdet = mxs->Determinant();
						if (tdet != 0)
						{
							sign = (j & 1) == 0 ? 1 : -1;
							det = det + (sign * data[0][j] * tdet);
						}
					}
				}
				delete mxs;
				return det;
			}
			else
			{
				if (N == 1) return Det1x1();
				if (N == 2) return Det2x2();
			}

			return 0;
		}

		void ToIdentity()
		{
			_ASSERT(N == M && N > 0);
			BEGIN_FOR_MATRIX_LOOP_ROW_RI_COLUMN_CI(N, M);
			this->data[ri][ci] = (T)(ri != ci ? 0 : 1);
			END_FOR_MATRIX_LOOP_ROW_RI_COLUMN_CI;
		}

		void Clear()
		{
			BEGIN_FOR_MATRIX_LOOP_ROW_RI_COLUMN_CI(N, M);
			this->data[ri][ci] = 0;
			END_FOR_MATRIX_LOOP_ROW_RI_COLUMN_CI;
		}

		unsigned int NSize()
		{
			return N;
		}

		unsigned int MSize()
		{
			return M;
		}

		static Matrix Identity()
		{
			_ASSERT(N == M && N > 0);
			Matrix result({});
			for (unsigned int ri = 0; ri < N; ri++)
			{
				result(ri, ri) = 1;
			}

			return result;
		}
	};

	template<typename T>
	class Matrix<T, 0, 0>
	{
	private:
		T data = 0;
	public:
		T& Data(unsigned int row, unsigned int column)
		{
			return data;
		}
		unsigned int Count() { return 0; }
		T Determinant()
		{
			return 0;
		}
		unsigned int NSize()
		{
			return 0;
		}

		unsigned int MSize()
		{
			return 0;
		}
	};

	typedef Matrix<float, 1, 2> Mat1x2f;
	typedef Matrix<float, 1, 3> Mat1x3f;
	typedef Matrix<float, 1, 4> Mat1x4f;
	typedef Matrix<float, 2, 1> Mat2x1f;
	typedef Matrix<float, 2, 2> Mat2x2f;
	typedef Matrix<float, 2, 3> Mat2x3f;
	typedef Matrix<float, 2, 4> Mat2x4f;
	typedef Matrix<float, 3, 1> Mat3x1f;
	typedef Matrix<float, 3, 2> Mat3x2f;
	typedef Matrix<float, 3, 3> Mat3x3f;
	typedef Matrix<float, 3, 4> Mat3x4f;
	typedef Matrix<float, 4, 1> Mat4x1f;
	typedef Matrix<float, 4, 2> Mat4x2f;
	typedef Matrix<float, 4, 3> Mat4x3f;
	typedef Matrix<float, 4, 4> Mat4x4f;

	typedef Matrix<int, 1, 2> Mat1x2i;
	typedef Matrix<int, 1, 3> Mat1x3i;
	typedef Matrix<int, 1, 4> Mat1x4i;
	typedef Matrix<int, 2, 1> Mat2x1i;
	typedef Matrix<int, 2, 2> Mat2x2i;
	typedef Matrix<int, 2, 3> Mat2x3i;
	typedef Matrix<int, 2, 4> Mat2x4i;
	typedef Matrix<int, 3, 1> Mat3x1i;
	typedef Matrix<int, 3, 2> Mat3x2i;
	typedef Matrix<int, 3, 3> Mat3x3i;
	typedef Matrix<int, 3, 4> Mat3x4i;
	typedef Matrix<int, 4, 1> Mat4x1i;
	typedef Matrix<int, 4, 2> Mat4x2i;
	typedef Matrix<int, 4, 3> Mat4x3i;
	typedef Matrix<int, 4, 4> Mat4x4i;

	typedef Matrix<double, 1, 2> Mat1x2d;
	typedef Matrix<double, 1, 3> Mat1x3d;
	typedef Matrix<double, 1, 4> Mat1x4d;
	typedef Matrix<double, 2, 1> Mat2x1d;
	typedef Matrix<double, 2, 2> Mat2x2d;
	typedef Matrix<double, 2, 3> Mat2x3d;
	typedef Matrix<double, 2, 4> Mat2x4d;
	typedef Matrix<double, 3, 1> Mat3x1d;
	typedef Matrix<double, 3, 2> Mat3x2d;
	typedef Matrix<double, 3, 3> Mat3x3d;
	typedef Matrix<double, 3, 4> Mat3x4d;
	typedef Matrix<double, 4, 1> Mat4x1d;
	typedef Matrix<double, 4, 2> Mat4x2d;
	typedef Matrix<double, 4, 3> Mat4x3d;
	typedef Matrix<double, 4, 4> Mat4x4d;

	// matrix pointers
	typedef Mat1x2f* Mat1x2fPtr;
	typedef Mat1x3f* Mat1x3fPtr;
	typedef Mat1x4f* Mat1x4fPtr;
	typedef Mat2x1f* Mat2x1fPtr;
	typedef Mat2x2f* Mat2x2fPtr;
	typedef Mat2x3f* Mat2x3fPtr;
	typedef Mat2x4f* Mat2x4fPtr;
	typedef Mat3x1f* Mat3x1fPtr;
	typedef Mat3x2f* Mat3x2fPtr;
	typedef Mat3x3f* Mat3x3fPtr;
	typedef Mat3x4f* Mat3x4fPtr;
	typedef Mat4x1f* Mat4x1fPtr;
	typedef Mat4x2f* Mat4x2fPtr;
	typedef Mat4x3f* Mat4x3fPtr;
	typedef Mat4x4f* Mat4x4fPtr;

	typedef Mat1x2i* Mat1x2iPtr;
	typedef Mat1x3i* Mat1x3iPtr;
	typedef Mat1x4i* Mat1x4iPtr;
	typedef Mat2x1i* Mat2x1iPtr;
	typedef Mat2x2i* Mat2x2iPtr;
	typedef Mat2x3i* Mat2x3iPtr;
	typedef Mat2x4i* Mat2x4iPtr;
	typedef Mat3x1i* Mat3x1iPtr;
	typedef Mat3x2i* Mat3x2iPtr;
	typedef Mat3x3i* Mat3x3iPtr;
	typedef Mat3x4i* Mat3x4iPtr;
	typedef Mat4x1i* Mat4x1iPtr;
	typedef Mat4x2i* Mat4x2iPtr;
	typedef Mat4x3i* Mat4x3iPtr;
	typedef Mat4x4i* Mat4x4iPtr;

	typedef Mat1x2d* Mat1x2dPtr;
	typedef Mat1x3d* Mat1x3dPtr;
	typedef Mat1x4d* Mat1x4dPtr;
	typedef Mat2x1d* Mat2x1dPtr;
	typedef Mat2x2d* Mat2x2dPtr;
	typedef Mat2x3d* Mat2x3dPtr;
	typedef Mat2x4d* Mat2x4dPtr;
	typedef Mat3x1d* Mat3x1dPtr;
	typedef Mat3x2d* Mat3x2dPtr;
	typedef Mat3x3d* Mat3x3dPtr;
	typedef Mat3x4d* Mat3x4dPtr;
	typedef Mat4x1d* Mat4x1dPtr;
	typedef Mat4x2d* Mat4x2dPtr;
	typedef Mat4x3d* Mat4x3dPtr;
	typedef Mat4x4d* Mat4x4dPtr;

