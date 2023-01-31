#include "stdio.h"
#include "Matrix.h"

template<typename T, unsigned int N, unsigned int M>
void PrintMatrix(Matrix<T, N, M> matrix)
{
    printf("\n\n");    
    for (unsigned int n = 0; n < matrix.NSize(); n++)
    {
        printf("\t|\t");
        for (unsigned int m = 0; m < matrix.MSize(); m++)
        {
            printf("%.8f\t", matrix(n, m));
        }
        printf("\t|\n");
    }
}

int main()
{
    Matrix<float, 2, 2> matrix( { { 2.0f, 3.0f }, { -1.0f, 5.0f } });
    printf("\n\nFloat Matrix (2x2):");
    PrintMatrix(matrix);
    printf("\n\tDeterminant: %.2f", matrix.Determinant());

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
    PrintMatrix(mat4f);
    printf("\n\nMatrix Inverse * Matrix = Identity:\n");
    PrintMatrix(mat4f * mat4x4f);

    printf("\n\nPress any key.\n");
    getchar();

    return 0;
}