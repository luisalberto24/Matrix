	#pragma once

	#include <stdio.h>
	#include <iostream>
	#include <string.h>
	#include <cassert>
	#include <stdlib.h>
	#include <concepts>
	#include "Concepts.h"

	using namespace std;

	template<unsigned int N, unsigned int M, typename Function>
	requires 
		nsConcepts::GreaterThanZero<N> && 
		nsConcepts::GreaterThanZero<M> &&
		std::invocable<Function, unsigned int, unsigned int>
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

	enum class ViewType { Row, Column };
	template<typename T, ViewType viewType>
	class VectorView
	{
		public:
			VectorView() noexcept = delete;
			VectorView(T* array, unsigned int size, unsigned int stride) noexcept :
				viewData(array),
				viewSize(size),
				viewStride(stride)
			{
			}

			VectorView& operator=(const std::initializer_list<T>& array) noexcept
			{
				assert(array.size() == viewSize);

				unsigned int r = 0;
				for (const T& value : array) { viewData[r * viewStride] = value; r++; }

				return *this;
			}

			T& operator()(unsigned int i) noexcept
			{
				assert(i < viewSize);
				return viewData[i * viewStride];
			}

			const T& operator()(unsigned int i) const noexcept
			{
				assert(i < viewSize);
				return viewData[i * viewStride];
			}

			unsigned int size() const noexcept
			{
				return viewSize;
			}

			template<typename P, ViewType Q>
			friend void Print(const VectorView<P, Q>&);

		private:

			T* viewData;
			unsigned int viewSize;
			unsigned int viewStride;
	};

	template <typename T, unsigned int N, unsigned M>
	requires 
		std::is_arithmetic_v<T> && 
		nsConcepts::GreaterThanZero<N> && 
		nsConcepts::GreaterThanZero<M>
	class BaseMatrix
	{
		public:
			using ArrayType = T[N][M];
			using RowArrayType = T[N];
			using ColumnArrayType = T[M];

		protected:
			ArrayType data{};
		public:
			BaseMatrix() noexcept = default;
			BaseMatrix(const BaseMatrix& matrix) noexcept { *this = matrix; }
			BaseMatrix(BaseMatrix&& matrix) noexcept { *this = matrix; }
			explicit BaseMatrix(const ArrayType&& array) noexcept { *this = array; };
			explicit BaseMatrix(const ArrayType& array) noexcept{ *this = array; }
			BaseMatrix(const VectorView<T, ViewType::Row>& row) noexcept { *this = row; }
			BaseMatrix(const VectorView<T, ViewType::Column>& column) noexcept { *this = column;}
			BaseMatrix(const std::initializer_list<std::initializer_list<T>>& array) noexcept
			{ 
				const size_t rowSize = array.size();
				if (!(rowSize == 1  && array.begin()->size() == 0))
				{
					assert(rowSize == N);
					unsigned int r = 0;
					for (const std::initializer_list<T>* row = array.begin(); row != array.end() && r < N; ++row, ++r)
					{
						assert(row->size() == M);
						unsigned int c = 0;
						for (const T* item = row->begin(); item != row->end() && c < M; ++item, ++c)
						{
							data[r][c] = *item;
						}
					}
				}
			}

		public:

			VectorView<T, ViewType::Row> Row(unsigned int r)
			{
				assert(r < N);
				return VectorView<T, ViewType::Row>(&data[r][0], N, 1);
			}

			VectorView<T, ViewType::Column> Column(unsigned int c)
			{
				assert(c < M);
				return VectorView<T, ViewType::Column>(&data[0][c], N, M);
			}

			const VectorView<const T, ViewType::Row> Row(unsigned int r) const
			{
				assert(r < N);
				return VectorView<const T, ViewType::Row>(&data[r][0], N, 1);
			}

			const VectorView<const T, ViewType::Column> Column(unsigned int c) const
			{
				assert(c < M);
				return VectorView<const T, ViewType::Column>(&data[0][c], N, M);
			}

			T& operator()(unsigned int r, unsigned int c)
			{
				assert(r < N && c < M);
				return data[r][c];
			}

			const T& operator()(unsigned int r, unsigned int c) const
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
				MatrixForLoopRowColumn<N, M>([&_data = data, &matrix](unsigned int r, unsigned int c) { _data[r][c] -= matrix(r, c); });
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

				MatrixForLoopRowColumn<N, P>([&_data = data, &result, &matrix](unsigned int r, unsigned int c) 
					{
						for (unsigned int x = 0; x < M; x++)
						{
							result(r, c) += _data[r][x] * matrix(x, c);
						}
					}
				);

				return result;
			}

			BaseMatrix& operator =(const VectorView<T, ViewType::Row>& row)
			{
				assert(N == 1 && row.size() == M);
				for(unsigned int c = 0; c < M; c++)
				{
					data[0][c] = row(c);
				}

				return *this;
			}

			BaseMatrix& operator =(const VectorView<T, ViewType::Column>& column)
			{
				assert(M == 1 && column.size() == N);
				for (unsigned int r = 0; r < N; r++)
				{
					data[r][0] = column(r);
				}

				return *this;
			}

			BaseMatrix& operator =(const BaseMatrix& matrix)
			{
				MatrixForLoopRowColumn<N, M>([&_data = data, &matrix](unsigned int r, unsigned int c) { _data[r][c] = matrix(r, c); });
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

			BaseMatrix& Add(const BaseMatrix& matrix)
			{
				MatrixForLoopRowColumn<N, M>([&_data = data, &matrix](unsigned int r, unsigned int c) { _data[r][c] += matrix(r, c); });
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

	template<typename T, unsigned int N, unsigned int M> 
	class Matrix : public BaseMatrix<T, N, M>
	{
		public:
		
			using ArrayType = typename BaseMatrix<T, N, M>::ArrayType;
			using MatrixPointer = Matrix*;

		public:

			Matrix() noexcept : BaseMatrix<T, N, M>(){}
			Matrix(const BaseMatrix<T, N, M>& matrix) noexcept : BaseMatrix<T, N, M>(matrix) {}
			Matrix(BaseMatrix<T, N, M>&& matrix) noexcept : BaseMatrix<T, N, M>(std::move(matrix)) {}
			explicit Matrix(const ArrayType&& array) noexcept : BaseMatrix<T, N, M>(std::move(array)) {}
			explicit Matrix(const ArrayType& array) noexcept : BaseMatrix<T, N, M>(array) {}
			Matrix(const VectorView<T, ViewType::Row>& row) noexcept : BaseMatrix<T, N, M>(row) {}
			Matrix(const VectorView<T, ViewType::Column>& column) noexcept : BaseMatrix<T, N, M>(column) {}
			Matrix(const std::initializer_list<std::initializer_list<T>>& array) noexcept : BaseMatrix<T, N, M>(array) {}

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

			Matrix() noexcept : BaseMatrix<T, 1, 1>() {}
			Matrix(const BaseMatrix<T, 1, 1>& matrix) noexcept : BaseMatrix<T, 1, 1>(matrix) {}
			Matrix(BaseMatrix<T, 1, 1>&& matrix) noexcept : BaseMatrix<T, 1, 1>(std::move(matrix)) {}
			explicit Matrix(const ArrayType&& array) noexcept : BaseMatrix<T, 1, 1>(std::move(array)) {}
			explicit Matrix(const ArrayType& array) noexcept : BaseMatrix<T, 1, 1>(array) {}
			Matrix(const VectorView<T, ViewType::Row>& row) noexcept : BaseMatrix<T, 1, 1>(row) {}
			Matrix(const VectorView<T, ViewType::Column>& column) noexcept : BaseMatrix<T, 1, 1>(column) {}
			Matrix(const std::initializer_list<std::initializer_list<T>>& array) noexcept : BaseMatrix<T, 1, 1>(array) {}

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
				return (this->data[0][0] != 0);
			}

			static Matrix<T, 1, 1> Identity()
			{
				return Matrix<T, 1, 1>{ {1} };
			}

			template<typename U, unsigned int P, unsigned int Q>
			friend void Print(BaseMatrix<U, 1, 1> matrix);
	};

	template<typename U, unsigned int P, unsigned int Q>
	void Print(BaseMatrix<U, P, Q> matrix)
	{
		constexpr const char* format = !std::is_same_v<U, int> ? "%.8f\t" : "%d\t";
		printf("\n\n");
		for (unsigned int n = 0; n < P; n++)
		{
			printf("\t|\t");			
			for (unsigned int m = 0; m < Q; m++)
			{
				printf(format, matrix(n, m));
			}
			printf("\t|\n");
		}
	}

	template<typename P, ViewType viewType>
	void Print(const VectorView<P, viewType>& vectorView)
	{
		constexpr const char* format = !std::is_same_v<P, int> ? "%.8f\t" : "%d\t";
		printf("\n\n");
		printf("View Type: %s", viewType == ViewType::Row ? "Row" : "Column");
		printf("\t|\t");
		
		for (unsigned int n = 0; n < vectorView.size(); n++)
		{
			printf(format, vectorView(n));
		}
		
		printf("\t|\n");
	}

	using Mat1x2f = Matrix<float, 1, 2>;
	using Mat1x3f = Matrix<float, 1, 3>;
	using Mat1x4f = Matrix<float, 1, 4>;
	using Mat2x1f = Matrix<float, 2, 1>;
	using Mat2x2f = Matrix<float, 2, 2>;
	using Mat2x3f = Matrix<float, 2, 3>;
	using Mat2x4f = Matrix<float, 2, 4>;
	using Mat3x1f = Matrix<float, 3, 1>;
	using Mat3x2f = Matrix<float, 3, 2>;
	using Mat3x3f = Matrix<float, 3, 3>;
	using Mat3x4f = Matrix<float, 3, 4>;
	using Mat4x1f = Matrix<float, 4, 1>;
	using Mat4x2f = Matrix<float, 4, 2>;
	using Mat4x3f = Matrix<float, 4, 3>;
	using Mat4x4f = Matrix<float, 4, 4>;

	using Mat1x2i = Matrix<int, 1, 2>;
	using Mat1x3i = Matrix<int, 1, 3>;
	using Mat1x4i = Matrix<int, 1, 4>;
	using Mat2x1i = Matrix<int, 2, 1>;
	using Mat2x2i = Matrix<int, 2, 2>;
	using Mat2x3i = Matrix<int, 2, 3>;
	using Mat2x4i = Matrix<int, 2, 4>;
	using Mat3x1i = Matrix<int, 3, 1>;
	using Mat3x2i = Matrix<int, 3, 2>;
	using Mat3x3i = Matrix<int, 3, 3>;
	using Mat3x4i = Matrix<int, 3, 4>;
	using Mat4x1i = Matrix<int, 4, 1>;
	using Mat4x2i = Matrix<int, 4, 2>;
	using Mat4x3i = Matrix<int, 4, 3>;
	using Mat4x4i = Matrix<int, 4, 4>;

	using Mat1x2d = Matrix<double, 1, 2>;
	using Mat1x3d = Matrix<double, 1, 3>;
	using Mat1x4d = Matrix<double, 1, 4>;
	using Mat2x1d = Matrix<double, 2, 1>;
	using Mat2x2d = Matrix<double, 2, 2>;
	using Mat2x3d = Matrix<double, 2, 3>;
	using Mat2x4d = Matrix<double, 2, 4>;
	using Mat3x1d = Matrix<double, 3, 1>;
	using Mat3x2d = Matrix<double, 3, 2>;
	using Mat3x3d = Matrix<double, 3, 3>;
	using Mat3x4d = Matrix<double, 3, 4>;
	using Mat4x1d = Matrix<double, 4, 1>;
	using Mat4x2d = Matrix<double, 4, 2>;
	using Mat4x3d = Matrix<double, 4, 3>;
	using Mat4x4d = Matrix<double, 4, 4>;

	using Mat1x2fPtr = Mat1x2f*;
	using Mat1x3fPtr = Mat1x3f*;
	using Mat1x4fPtr = Mat1x4f*;
	using Mat2x1fPtr = Mat2x1f*;
	using Mat2x2fPtr = Mat2x2f*;
	using Mat2x3fPtr = Mat2x3f*;
	using Mat2x4fPtr = Mat2x4f*;
	using Mat3x1fPtr = Mat3x1f*;
	using Mat3x2fPtr = Mat3x2f*;
	using Mat3x3fPtr = Mat3x3f*;
	using Mat3x4fPtr = Mat3x4f*;
	using Mat4x1fPtr = Mat4x1f*;
	using Mat4x2fPtr = Mat4x2f*;
	using Mat4x3fPtr = Mat4x3f*;
	using Mat4x4fPtr = Mat4x4f*;

	using Mat1x2iPtr = Mat1x2i*;
	using Mat1x3iPtr = Mat1x3i*;
	using Mat1x4iPtr = Mat1x4i*;
	using Mat2x1iPtr = Mat2x1i*;
	using Mat2x2iPtr = Mat2x2i*;
	using Mat2x3iPtr = Mat2x3i*;
	using Mat2x4iPtr = Mat2x4i*;
	using Mat3x1iPtr = Mat3x1i*;
	using Mat3x2iPtr = Mat3x2i*;
	using Mat3x3iPtr = Mat3x3i*;
	using Mat3x4iPtr = Mat3x4i*;
	using Mat4x1iPtr = Mat4x1i*;
	using Mat4x2iPtr = Mat4x2i*;
	using Mat4x3iPtr = Mat4x3i*;
	using Mat4x4iPtr = Mat4x4i*;

	using Mat1x2dPtr = Mat1x2d*;
	using Mat1x3dPtr = Mat1x3d*;
	using Mat1x4dPtr = Mat1x4d*;
	using Mat2x1dPtr = Mat2x1d*;
	using Mat2x2dPtr = Mat2x2d*;
	using Mat2x3dPtr = Mat2x3d*;
	using Mat2x4dPtr = Mat2x4d*;
	using Mat3x1dPtr = Mat3x1d*;
	using Mat3x2dPtr = Mat3x2d*;
	using Mat3x3dPtr = Mat3x3d*;
	using Mat3x4dPtr = Mat3x4d*;
	using Mat4x1dPtr = Mat4x1d*;
	using Mat4x2dPtr = Mat4x2d*;
	using Mat4x3dPtr = Mat4x3d*;
	using Mat4x4dPtr = Mat4x4d*;

