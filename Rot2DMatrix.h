#pragma once 

#include "Matrix.h"

namespace nsMatrix
{
	template<typename T>
    requires std::is_floating_point_v<T>
	constexpr T MATH_PI = T(3.14159265358979323846);

    template<typename T, unsigned int N, unsigned M>
    requires
        EqualsValue<N, 2> &&
        EqualsValue<M, 2>
    class Rot2DMatrix : public BaseMatrix<T, N, M> {
    public:
        Rot2DMatrix() {};
        Rot2DMatrix(T x1, T y1, T x2, T y2) 
        {
			this->data[0][0] = x1;
			this->data[0][1] = y1;
            this->data[1][0] = x2;
            this->data[1][1] = y2;
        };

        Rot2DMatrix(const std::initializer_list<std::initializer_list<T>>& array) : BaseMatrix<T, N, M>(array) {}

        static Rot2DMatrix<T, N, M> fromRadians(T radians)
        {
            return  Rot2DMatrix<T, N, M>(cos(radians), -sin(radians), sin(radians), cos(radians));
        }

        static Rot2DMatrix<T, N, M> fromAngle(T degrees)
        {
            T radians = degrees * (MATH_PI<T> / T{180});
            return  Rot2DMatrix<T, N, M>(cos(radians), -sin(radians), sin(radians), cos(radians));
        }

        template<typename U, unsigned int P, unsigned int Q>
        friend void Print(BaseMatrix<U, P, Q> matrix);
    };
    
    typedef Rot2DMatrix<float, 2, 2> Rot2DMatrixf;
    typedef Rot2DMatrix<double, 2, 2> Rot2DMatrixd;
}
