	#pragma once

	#include <stdio.h>
	#include <iostream>
	#include <string.h>
	#include <cassert>
	#include <concepts>

	using namespace std;

	template <unsigned int N>
	concept GreaterThanZero = (N > 0);

	template<unsigned int N, unsigned int M, typename Function>
	requires GreaterThanZero<N> && GreaterThanZero<M> && std::invocable<Function, unsigned int, unsigned int>
	void MatrixForLoopRowColumn(Function executeFunction)
	{
		for (unsigned int r = 0; r < N; r++)
		{
			for (unsigned int c = 0; c < M; c++)
			{
				executeFunction(r, c);
			}
		}
	}

	template <typename T, unsigned int N, unsigned M>
	requires (GreaterThanZero<N> && GreaterThanZero<M>)
	class BaseMatrix
	{
		public:		
			using ArrayType = T[N][M];
		protected:
			ArrayType data;
		public:
			BaseMatrix()
			{
				*this = { {} };
			}

			BaseMatrix(const BaseMatrix& matrix) noexcept
			{
				*this = matrix;
			}

			BaseMatrix(BaseMatrix&& matrix) noexcept
			{
				*this = matrix;
			}

			BaseMatrix(ArrayType&& array) noexcept
			{
				*this = array;
			}

			BaseMatrix(ArrayType& array) noexcept
			{
				*this = array;
			}

		public:
			T& operator()(unsigned int r, unsigned int c)
			{
				assert(r < N && c < M);
				return data[r][c];
			}

			BaseMatrix operator +(const BaseMatrix& matrix) const
			{
				BaseMatrix result(*this);
				result += matrix;
				return result;
			}

			BaseMatrix& operator +=(const BaseMatrix& matrix)
			{
				Add(matrix);
				return *this;
			}

			BaseMatrix operator -(const BaseMatrix& matrix) const
			{
				BaseMatrix result(*this);
				result -= matrix;
				return result;
			}

			BaseMatrix& operator -=(const BaseMatrix& matrix)
			{
				MatrixForLoopRowColumn<N, M>([&_data = data, &matrix](unsigned int r, unsigned int c) { _data[r][c] -= matrix.data[r][c]; });
				return *this;
			}

			BaseMatrix operator -(const T value) const
			{
				BaseMatrix result{};
				MatrixForLoopRowColumn<N, M>([&_data = data, &result, &value](unsigned int r, unsigned int c) { result(r, c) = _data[r][c] - value; });
				return result;
			}

			BaseMatrix& operator -=(const T value)
			{
				MatrixForLoopRowColumn<N, M>([&_data = data, &value](unsigned int r, unsigned int c) { _data[r][c] = _data[r][c] - value; });
				return *this;
			}

			BaseMatrix operator +(const T value) const
			{
				BaseMatrix result{};
				MatrixForLoopRowColumn<N, M>([&_data = data, &result, &value](unsigned int r, unsigned int c) { result(r, c) = _data[r][c] + value; });
				return result;
			}

			BaseMatrix& operator +=(const T value)
			{
				MatrixForLoopRowColumn<N, M>([&_data = data, &value](unsigned int r, unsigned int c) { _data[r][c] = _data[r][c] + value; });
				return *this;
			}

			BaseMatrix operator /(const T value) const
			{
				assert(value != 0);

				BaseMatrix result{};
				MatrixForLoopRowColumn<N, M>([&_data = data, &result, &value](unsigned int r, unsigned int c) { result(r, c) = _data[r][c] / value; });
				return result;
			}

			BaseMatrix& operator /=(const T value)
			{
				assert(value != 0);

				MatrixForLoopRowColumn<N, M>([&_data = data, &value](unsigned int r, unsigned int c) { _data[r][c] = _data[r][c] / value; });
				return *this;
			}

			BaseMatrix operator *(const T value) const
			{
				BaseMatrix result{};
				MatrixForLoopRowColumn<N, M>([&_data = data, &result, &value](unsigned int r, unsigned int c) { result(r, c) = _data[r][c] * value; });
				return result;
			}

			BaseMatrix& operator *=(const T value)
			{
				MatrixForLoopRowColumn<N, M>([&_data = data, &value](unsigned int r, unsigned int c) { _data[r][c] = _data[r][c] * value; });
				return *this;
			}

			BaseMatrix& operator ^=(const T exponent)
			{
				T absExponent = abs(exponent);
				bool exponentHasDecimalPlaces = ((absExponent - (int)absExponent) != 0);
				MatrixForLoopRowColumn<N, M>([&_data = data, &exponent, &exponentHasDecimalPlaces](unsigned int r, unsigned int c)
					{
						assert(_data[r][c] >= 0 || (_data[r][c] < 0 && !exponentHasDecimalPlaces));
						_data[r][c] = pow(_data[r][c], exponent);
					}
				);

				return *this;
			}

			template<unsigned int P>
			BaseMatrix<T, N, P> operator *(const BaseMatrix<T, M, P>& matrix) const
			{
				BaseMatrix<T, N, P> result{};
				MatrixForLoopRowColumn<N, M>([&_data = data, &result, &matrix](unsigned int r, unsigned int c) {
					for (unsigned int x = 0; x < M; x++)
					{
						result.data[r][c] += _data[r][x] * matrix.data[x][c];
					}
					});

				return result;
			}

			BaseMatrix& operator =(const BaseMatrix& matrix)
			{
				MatrixForLoopRowColumn<N, M>([&_data = data, &matrix](unsigned int r, unsigned int c) { _data[r][c] = matrix.data[r][c]; });
				return *this;
			}

			BaseMatrix& operator =(const ArrayType& array)
			{
				MatrixForLoopRowColumn<N, M>([&_data = data, &array](unsigned int r, unsigned int c) { _data[r][c] = array[r][c]; });
				return *this;
			}

			BaseMatrix& operator =(const ArrayType&& array)
			{
				MatrixForLoopRowColumn<N, M>([&_data = data, &array](unsigned int r, unsigned int c) { _data[r][c] = array[r][c]; });
				return *this;
			}

			const ArrayType& Data() const
			{
				return data;
			}

			ArrayType& Data()
			{
				return data;
			}

			BaseMatrix<T, 1, M> GetRow(int row) const
			{
				assert(row >= 0 && row < N);
				BaseMatrix<T, 1, M> result{};
				for (unsigned int m = 0; m < M; m++)
				{
					result(0, m) = data[row][m];
				}

				return result;
			}

			BaseMatrix<T, N, 1> GetColumn(int column) const
			{
				assert(column >= 0 && column < M);
				BaseMatrix<T, N, 1> result{};
				for (unsigned int n = 0; n < N; n++)
				{
					result.Data()[n][0] = data[n][column];
				}

				return result;
			}

			BaseMatrix& Add(const BaseMatrix& matrix)
			{
				MatrixForLoopRowColumn<N, M>([&_data = data, &matrix](unsigned int r, unsigned int c) { _data[r][c] += matrix.data[r][c]; });
				return *this;
			}

			BaseMatrix<T, M, N> Traspose() const
			{
				BaseMatrix<T, M, N> result{};
				MatrixForLoopRowColumn<N, M>([&_data = data, &result](unsigned int r, unsigned int c) { result(c, r) = _data[r][c]; });
				return result;
			}

			void ToIdentity()
			{
				assert(N == M && N > 1);
				MatrixForLoopRowColumn<N, M>([&_data = data](unsigned int r, unsigned int c) { _data[r][c] = (T)(r != c ? 0 : 1); });
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
				MatrixForLoopRowColumn<N, M>([&_data = data](unsigned int r, unsigned int c) { _data[r][c] = 0; });
			}

			unsigned int NSize()
			{
				return N;
			}

			unsigned int MSize()
			{
				return M;
			}

			template<typename U, unsigned int P, unsigned int Q>
			friend void Print(BaseMatrix<U, P, Q> matrix);
	};

	template<typename U, unsigned int P, unsigned int Q>
	void Print(BaseMatrix<U, P, Q> matrix)
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

	template<typename T, unsigned int N, unsigned int M> 
	class Matrix : public BaseMatrix<T, N, M>
	{
		public:
		
			using ArrayType = typename BaseMatrix<T, N, M>::ArrayType;
			using MatrixPointer = Matrix*;

		public:

			Matrix() : BaseMatrix<T, N, M>(){}
			Matrix(const BaseMatrix<T, N, M>& matrix) noexcept : BaseMatrix<T, N, M>(matrix) {}
			Matrix(BaseMatrix<T, N, M>&& matrix) noexcept : BaseMatrix<T, N, M>(std::move(matrix)) {}
			Matrix(ArrayType&& array) noexcept : BaseMatrix<T, N, M>(std::move(array)) {}
			Matrix(ArrayType& array) noexcept : BaseMatrix<T, N, M>(array) {}

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
										mxs(k, p) = this->data[n][m];
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
				assert(N == M);
				if (N > 2)
				{
					T det = 0;
					T tdet = 0;
					unsigned int s = 0, k = 0, p = 0, m = 0;
					int sign = 0;
					Matrix<T, M - 1, M - 1> mxs{};
					for (unsigned int j = 0; j < M; j++)
					{
						if (this->data[0][j] != 0)
						{
							k = 0;
							for (unsigned int n = 1; n < M; n++)
							{
								p = 0;
								m = 0;
								while (m < M)
								{
									if (j != m)
									{
										mxs(k, p) = this->data[n][m];
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
								det = det + (sign * this->data[0][j] * tdet);
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

			Matrix Inverse() const
			{
				assert(N == M && N > 1);
				T det = Determinant();
				assert(det != 0);

				return Adjoint() / det;
			}

			bool IsInvertible()
			{
				assert(N == M);
				return (Determinant() != 0);
			}

			static Matrix Identity()
			{
				assert(N == M);
				Matrix result{};
				for (unsigned int ri = 0; ri < N; ri++)
				{
					result.data[ri][ri] = 1;
				}

				return result;
			}
		
			template<typename U, unsigned int P, unsigned int Q>
			friend void Print(BaseMatrix<U, P, Q> matrix);

		private:

				T Det1x1() const { assert(N == M && N == 1); return this->data[0][0]; }
				T Det2x2() const
				{
					assert(N == M && N == 2);
					return this->data[0][0] * this->data[1][1] - this->data[1][0] * this->data[0][1];
				}
	};

	template<typename T>
	class Matrix<T, 1, 1> : public BaseMatrix<T, 1, 1>
	{
		public:

			using ArrayType = typename BaseMatrix<T, 1, 1>::ArrayType;
			using MatrixPointer = Matrix*;

		public:

			Matrix() : BaseMatrix<T, 1, 1>() {}
			Matrix(const BaseMatrix<T, 1, 1>& matrix) noexcept : BaseMatrix<T, 1, 1>(matrix) {}
			Matrix(BaseMatrix<T, 1, 1>&& matrix) noexcept : BaseMatrix<T, 1, 1>(std::move(matrix)) {}
			Matrix(ArrayType&& array) noexcept : BaseMatrix<T, 1, 1>(std::move(array)) {}
			Matrix(ArrayType& array) noexcept : BaseMatrix<T, 1, 1>(array) {}

			Matrix<T, 1, 1> Adjoint() const
			{
				return Matrix{ {1} };
			}

			T Determinant() const
			{
				return this->data[0][0];
			}

			Matrix<T, 1, 1> Inverse() const
			{
				assert(data[0][0] != (T)0);
				Matrix result{};
				result.data[0][0] = (T)1/data[0][0];
				return result;
			}

			bool IsInvertible()
			{
				return (Determinant() != 0);
			}

			static Matrix<T, 1, 1> Identity()
			{
				return Matrix<T, 1, 1>{ {1} };
			}

			template<typename U, unsigned int P, unsigned int Q>
			friend void Print(BaseMatrix<U, 1, 1> matrix);
	};

	// typedef Matrix<float, 1, 0> Mat1f;
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

