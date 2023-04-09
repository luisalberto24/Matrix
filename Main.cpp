#include "stdio.h"
#include "Matrix.h"

template<typename T, unsigned int N, unsigned int M>
void PrintMatrix(Matrix<T, N, M> matrix)
{
    printf("\n\n");    
    for (unsigned int n = 0; n < N; n++)
    {
        printf("\t|\t");
        for (unsigned int m = 0; m < M; m++)
        {
            printf("%.8f\t", matrix(n, m));
        }
        printf("\t|\n");
    }
}

int main()
{
    Matrix<double, 2, 2> matrix( { { 2.0f, 3.0f }, { -1.0f, 5.0f } });
    matrix ^= 8;
    PrintMatrix(matrix);
    printf("\n\nExtracting GetColumn:");
    Mat2x1d mxt = matrix.GetColumn(0);
    PrintMatrix(mxt);

    printf("\n\nExtracting GetRow:");
    Mat1x2d mxt2 = matrix.GetRow(0);
    PrintMatrix(mxt2.Traspose());
    
    printf("\n\nFloat Matrix (2x2):");
    PrintMatrix(matrix);
    if (matrix.IsInvertible())
    {
        printf("\n\tDeterminant: %.2f", matrix.Determinant());
    }
    else
    {
        printf("Matrix is not invertible\n");
    }

    Mat4x4f mat4x4f({ {3,5,8,54}, {32, 92,-34,33}, {8,4,3,2}, {-24,-98,-43,-3} });
    printf("\n\nFloat Matrix (4x4):");
    PrintMatrix(mat4x4f);
    printf("\n\tDeterminant: %.2f", mat4x4f.Determinant());
    printf("\n\nMatrix Identity:\n");
    PrintMatrix(mat4x4f.Identity());
    printf("\n\nMatrix Traspose:\n");
    PrintMatrix(mat4x4f.Traspose());
    printf("\n\nMatrix Inverse:\n");
    Mat4x4f mat4f = mat4x4f.Inverse();
    Mat4x4f Mat4x4fIdent = Mat4x4f::Identity();

    printf("%s", Mat4x4fIdent.IsIdentity() ? "Is an identity." : "Is not an identity.");

    PrintMatrix(mat4f);
    printf("\n\nMatrix Inverse * Matrix = Identity:\n");
    Mat4x4f mat4tf = mat4f * mat4x4f;

    
    PrintMatrix(mat4tf);

    printf("\n\nPress any key.\n");
    getchar();

    return 0;
}