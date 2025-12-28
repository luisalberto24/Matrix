	#pragma once

	#include <stdio.h>
	#include <iostream>
	#include <string.h>
	#include <fstream>
	#include <stdlib.h>
	#include <type_traits>
	#include <cassert>
	#include <concepts>

	using namespace std;

	template<typename Row, typename Column, typename Function>
	requires std::invocable<Function, unsigned int, unsigned int>
	void Matrix_For_Loop_Row_Column(Row N, Column M, Function functionCall)
	{
		for (unsigned int ri = 0; ri < N; ri++)
		{
			for (unsigned int ci = 0; ci < M; ci++)
			{
				functionCall(ri, ci);
			}
		}
	}

	template<typename T, unsigned int N, unsigned int M>
	class Matrix
	{
	public:
		typedef Matrix* MatrixPointer;
		using ArrayType = T[N][M];
	protected:
		Matrix::ArrayType data;
	public:
		Matrix()
		{
			*this = { {} };
		}

		Matrix(const Matrix& matrix) noexcept
		{
			*this = matrix;
		}

		Matrix(Matrix&& matrix) noexcept
		{
			*this = matrix;
		}

		Matrix(Matrix::ArrayType&& value) noexcept
		{
			*this = value;
		}

		Matrix(Matrix::ArrayType& value) noexcept
		{
			*this = value;
		}

	private:
		
		T Det1x1() const { assert(N == M && N == 1); return data[0][0]; }
		T Det2x2() const
		{
			assert(N == M && N == 2);
			return data[0][0] * data[1][1] - data[1][0] * data[0][1];
		}

	public:
		T& operator()(unsigned int row, unsigned int column)
		{
			assert(row < N&& column < M && N > 0 && M > 0 && (N > 1 || M > 1));
			return data[row][column];
		}

		Matrix operator +(const Matrix& matrix) const
		{
			Matrix result(*this);
			result += matrix;
			return result;
		}

		Matrix& operator +=(const Matrix& matrix)
		{
			Add(matrix);
			return *this;
		}

		Matrix operator -(const Matrix& matrix) const
		{
			Matrix result(*this);
			result -= matrix;
			return result;
		}

		Matrix& operator -=(const Matrix& matrix)
		{
			Matrix_For_Loop_Row_Column(N, M, [&_data = data, &matrix](unsigned int ri, unsigned int ci) { _data[ri][ci] -= matrix.data[ri][ci]; });
			return *this;
		}

		Matrix operator -(const T value) const
		{
			Matrix result{};
			Matrix_For_Loop_Row_Column(N, M, [&_data = data, &result, value](unsigned int ri, unsigned int ci) { result(ri, ci) = _data[ri][ci] - value; });
			return result;
		}

		Matrix& operator -=(const T value)
		{
			Matrix_For_Loop_Row_Column(N, M, [&_data = data, value](unsigned int ri, unsigned int ci) { _data[ri][ci] = _data[ri][ci] - value; });
			return *this;
		}

		Matrix operator +(const T value) const
		{
			Matrix result{};
			Matrix_For_Loop_Row_Column(N, M, [&_data = data, &result, value](unsigned int ri, unsigned int ci) { result(ri, ci) = _data[ri][ci] + value; });
			return result;
		}

		Matrix& operator +=(const T value)
		{
			Matrix_For_Loop_Row_Column(N, M, [&_data = data, value](unsigned int ri, unsigned int ci) { _data[ri][ci] = _data[ri][ci] + value; });
			return *this;
		}

		Matrix operator /(const T value) const
		{
			assert(value != 0);

			Matrix result{};
			Matrix_For_Loop_Row_Column(N, M, [&_data = data, &result, value](unsigned int ri, unsigned int ci) { result(ri, ci) = _data[ri][ci] / value; });
			return result;
		}

		Matrix& operator /=(const T value)
		{
			assert(value != 0);

			Matrix_For_Loop_Row_Column(N, M, [&_data = data, value](unsigned int ri, unsigned int ci) { _data[ri][ci] = _data[ri][ci] / value; });
			return *this;
		}

		Matrix operator *(const T value) const
		{
			Matrix result{};
			Matrix_For_Loop_Row_Column(N, M, [&_data = data, &result, value](unsigned int ri, unsigned int ci) { result(ri, ci) = _data[ri][ci] * value; });
			return result;
		}

		Matrix& operator *=(const T value)
		{
			Matrix_For_Loop_Row_Column(N, M, [&_data = data, value](unsigned int ri, unsigned int ci) { _data[ri][ci] = _data[ri][ci] * value; });
			return *this;
		}

		Matrix& operator ^=(const T exponent)
		{
			T absExponent = abs(exponent);
			assert((absExponent - (int)absExponent) == 0);
			Matrix_For_Loop_Row_Column(N, M, [&_data = data, exponent](unsigned int ri, unsigned int ci) 
				{ 
					assert(_data[ri][ci] >= 0 || _data[ri][ci] < 0);
					_data[ri][ci] = pow(_data[ri][ci], exponent);
				}
			);

			return *this;
		}

		template<unsigned int P>
		Matrix<T, N, P> operator *(const Matrix<T, M, P>& matrix) const 
		{
			Matrix<T, N, P> result{};
			Matrix_For_Loop_Row_Column(N, M, [&_data = data, &result, &matrix](unsigned int ri, unsigned int ci) {
				for (unsigned int x = 0; x < M; x++)
				{
					result.data[ri][ci] += _data[ri][x] * matrix.data[x][ci];
				}
			});

			return result;
		}

		Matrix& operator =(const Matrix& matrix)
		{
			Matrix_For_Loop_Row_Column(N, M, [&_data = data, &matrix](unsigned int ri, unsigned int ci) { _data[ri][ci] = matrix.data[ri][ci]; });
			return *this;
		}

		Matrix& operator =(const Matrix::ArrayType& value)
		{
			Matrix_For_Loop_Row_Column(N, M, [&_data = data, &value](unsigned int ri, unsigned int ci) { _data[ri][ci] = value[ri][ci]; });
			return *this;
		}

		Matrix& operator =(const Matrix::ArrayType&& value)
		{
			Matrix_For_Loop_Row_Column(N, M, [&_data = data, &value](unsigned int ri, unsigned int ci) { _data[ri][ci] = value[ri][ci]; });
			return *this;
		}

		const Matrix::ArrayType& Data() const 
		{
			return data;
		}

		Matrix::ArrayType& Data()
		{
			return data;
		}

		Matrix<T, 1, M> GetRow(int row)
		{
			assert(row >= 0 && row < N);
			Matrix<T, 1, M> result;
			for (unsigned int m = 0; m < M; m++)
			{
				result(0, m) = data[row][m];
			}

			return result;
		}

		Matrix<T, N, 1> GetColumn(int column)
		{
			assert(column >= 0 && column < M);
			Matrix<T, N, 1> result;
			for (unsigned int n = 0; n < N; n++)
			{
				result.Data()[n][0] = data[n][column];
			}

			return result;
		}
		
		Matrix& Add(const Matrix& matrix)
		{
			Matrix_For_Loop_Row_Column(N, M, [&_data = data, &matrix](unsigned int ri, unsigned int ci) { _data[ri][ci] += matrix.data[ri][ci]; });
			return *this;
		}

		Matrix<T, M, N> Traspose() const 
		{
			Matrix<T, M, N> result{};
			Matrix_For_Loop_Row_Column(N, M, [&_data = data, &result](unsigned int ri, unsigned int ci) { result(ci, ri) = _data[ri][ci]; });
			return result;
		}

		Matrix Inverse() const
		{
			assert(N == M && N > 1);
			T det = Determinant();
			assert(det != 0);

			return Adjoint() / det;
		}

		Matrix Adjoint() const
		{
			assert(N == M && N > 1);
			int sign = 0;
			T det = 0;
			unsigned int n = 0, m = 0, p = 0, k = 0;
			Matrix matrix{};
			Matrix<T, M - 1, M - 1> mxs{};
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
									mxs(k, p) = data[n][m];
									p++;
								}
								m++;
							}
							k++;
						}
						n++;
					}

					det = mxs.Determinant();
					if (det != 0)
					{
						sign = ((fi + fj) & 1) == 0 ? 1 : -1;
						matrix(fi, fj) = sign * det;
					}
				}
			}

			return matrix.Traspose();
		}

		T Determinant() const
		{
			assert(N == M && N > 0);
			if (N > 2)
			{
				T det = 0;
				T tdet = 0;
				unsigned int s = 0, k = 0, p = 0, m = 0;
				int sign = 0;
				Matrix<T, M - 1, M - 1> mxs{};
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
									mxs(k, p) = data[n][m];
									p++;
								}
								m++;
							}
							k++;
						}

						tdet = mxs.Determinant();
						if (tdet != 0)
						{
							sign = (j & 1) == 0 ? 1 : -1;
							det = det + (sign * data[0][j] * tdet);
						}
					}
				}
				return det;
			}
			else
			{
				if (N == 1) return Det1x1();
				if (N == 2) return Det2x2();
			}

			return 0;
		}
		bool IsInvertible()
		{
			assert(N == M && N > 0);
			return (Determinant() != 0);
		}

		void ToIdentity()
		{
			assert(N == M && N > 1);
			Matrix_For_Loop_Row_Column(N, M, [&_data = data](unsigned int ri, unsigned int ci) { _data[ri][ci] = (T)(ri != ci ? 0 : 1); });
		}

		bool IsIdentity() const
		{
			assert(N == M && N > 1);

			for (unsigned int ri = 0; ri < N; ri++)
			{
				for (unsigned int ci = 0; ci < M; ci++)
				{
					if (ri == ci)
					{
						if (data[ri][ci] != 1)
						{
							return false;
						}
					}
					else if (data[ri][ci] != 0)
					{
						return false;
					}
				}
			}

			return true;
		}
	
		void Clear()
		{
			Matrix_For_Loop_Row_Column(N, M, [&_data = data](unsigned int ri, unsigned int ci) { _data[ri][ci] = 0; });
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
			assert(N == M && N > 0);
			Matrix result{};
			for (unsigned int ri = 0; ri < N; ri++)
			{
				result.data[ri][ri] = 1;
			}

			return result;
		}

		template<typename U, unsigned int P, unsigned int Q>
		friend void PrintMatrix(Matrix<U, P, Q> matrix);
	};

	template<typename T>
	class Matrix<T, 0, 0>
	{
		protected:
			T data = 0;
		public:
			
			T& operator()(unsigned int row, unsigned int column)
			{
				return data;
			}

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

			template<typename U, unsigned int P, unsigned int Q>
			friend void PrintMatrix(Matrix<U, P, Q> matrix);
	};

	template<typename U, unsigned int P, unsigned int Q>
	void PrintMatrix(Matrix<U, P, Q> matrix)
	{
		printf("\n\n");
		for (unsigned int n = 0; n < P; n++)
		{
			printf("\t|\t");
			constexpr const char* format = !std::is_same_v<U, int> ? "%.8f\t" : "%d\t";
			for (unsigned int m = 0; m < Q; m++)
			{
				printf(format, matrix(n, m));
			}
			printf("\t|\n");
		}
	}

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

