#include "stdio.h"
#include "Matrix.h"
#include "TypeFactory.h"
#include <direct.h>
#include <limits.h>
#include <string.h>

int main()
{
    printf("\nStart\n");
    std::unique_ptr<Mat2x2i> mat2_2i_u = std::make_unique<Mat2x2i>(Mat2x2i::ArrayType{ {2,3}, {4,5} });
    Print(*mat2_2i_u.get());
    printf("Item: %d\n", mat2_2i_u->Data()[0][1]);
    printf("\nEnd\n");

    printf("\nxxxx==================================================================================================xxx\n");
    Mat4x4d::ArrayType array = { {3,5,8,54}, {32, 92,-34,33}, {800.23,4,3,2}, {-24,-98,-43,-3} };
    printf("\n==================================================================================================\n");
    printf("\nCreate unique ptr matrix with TypeFactory:\n");
    std::unique_ptr<Mat4x4d> mat4_4d_u = TypeFactory::Create_Unique_Ptr<Mat4x4d>(array);
    Print(*mat4_4d_u);
    printf("\n==================================================================================================\n");

    printf("\nCreate unique shared matrix with TypeFactory:\n");
    std::shared_ptr<Mat4x4d> mat4_4d_s = TypeFactory::Create_Shared_Ptr<Mat4x4d>(array);
    Print(*mat4_4d_s);
    printf("\n==================================================================================================\n");

    printf("\nCreate direct pointer for matrix with TypeFactory:\n");
    Mat4x4d* mat4_4d = TypeFactory::Create<Mat4x4d>(array);
    Print(*mat4_4d);
    delete mat4_4d;
    printf("\n==================================================================================================\n");

    Matrix<double, 2, 2> matrix({ { 2.0f, 3.0f }, { -1.0f, 5.0f } });
    matrix ^= 8;
    Print(matrix);
    printf("\n\nExtracting GetColumn:");
    Mat2x1d mxt = matrix.GetColumn(0);
    Print(mxt);
    printf("\n\nExtracting GetRow:");
    Mat1x2d mxt2 = matrix.GetRow(0);
    Print(mxt2.Traspose());

    printf("\n\nFloat Matrix (2x2):");
    Print(matrix);
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
    Print(mat4x4f);
    printf("\n\tDeterminant: %.2f", mat4x4f.Determinant());
    printf("\n\nMatrix Identity:\n");
    Print(mat4x4f.Identity());
    printf("\n\nMatrix Traspose:\n");
    Print(mat4x4f.Traspose());
    printf("\n\nMatrix Inverse:\n");
    Mat4x4f mat4f = mat4x4f.Inverse();
    Mat4x4f Mat4x4fIdent = Mat4x4f::Identity();

    printf("%s", Mat4x4fIdent.IsIdentity() ? "Is an identity." : "Is not an identity.");

    Print(mat4f);
    printf("\n\nMatrix Inverse * Matrix = Identity:\n");
    Mat4x4f mat4tf = mat4f * mat4x4f;


    Print(mat4tf);

    printf("\n\nPress any key.\n");
    int v = getchar();

    return 0;
}