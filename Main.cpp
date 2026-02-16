#include "stdio.h"
#include "Matrix.h"
#include "TypeFactory.h"
#include "Buffer.h"
#include "BufferView.h"
#include <direct.h>
#include <limits.h>
#include <string.h>

int main()
{
    Mat3x3f mat3_3f = { {-0.25f, -1.0f, 0.35f}, {1.0f, 10.32f, -98.32f}, {0.89f, 0.98f, -1.25f} };
	Print(mat3_3f);
    printf("\nIs Orthonormal: %s", mat3_3f.IsOrthonormal() ? "True" : "False");

    Mat1x3f mat1_3f = mat3_3f.Row(0);
    Print(mat1_3f);

	Buffer<int, 5> bv;
    Buffer<std::string, 2> info;
    info[0] = "this is a test\n";
    printf("\n%s", info[0].c_str());

    std::string b[7] = { "One","Two", "Three", "Four", "Five", "Six", "Seven"};
    
    BufferView<std::string, 7> bufferView1(b);
    BufferView<std::string, 7> bufferView2(&b[0], &b[6], 6);
    int xx = 0;
    printf("\nBufferView 1\n\n");
    for (const auto& item : bufferView1)
    {
        printf("integer %d: %s\n", ++xx, item.c_str());
    }
    
    printf("\nBufferView 2\n\n");
    for (const auto& item : bufferView2)
    {
        printf("integer %d: %s\n", ++xx, item.c_str());
    }
    
    printf("\nPrint Mat6f - Start Process...\n");
    Mat6f matTest(1.0f);
	Print(matTest);
    printf("\nPrint Mat6f - End Process...\n");

    Buffer<Mat2x2f, 2> mat2x2fBuffer;
    mat2x2fBuffer[0] = Mat2x2f::Identity();
    mat2x2fBuffer[1] = { {-1, 2}, {356, 4} };

    // Begin -- Accessing matrix elements using Row and Column views.
    mat2x2fBuffer[1].Row(0) = {-800.0f, -5072.3f};
    mat2x2fBuffer[1].Row(1) = { -99.1f, -9111.3f };
    mat2x2fBuffer[1].Column(0) = { -11.35f, -88.1f };
    // End -- Accessing matrix elements using Row and Column views.

    printf("\nBegin -- Accessing matrix elements using Row and Column views and printing them.\n");
    // Print(mat2x2fBufferView[1].Row(0));
    // Print(mat2x2fBufferView[1].Row(1));
    Print(mat2x2fBuffer[1].Column(0));
    // Print(mat2x2fBufferView[1].Column(1));
    printf("\nEnd -- Accessing matrix elements using Row and Column views and printing them.\n");

    printf("\nBegin -- Accessing matrix elements using Row and Column views to convert them to Matrixes.\n");
    printf("Mat 1x2 Float from Row:");
    Mat1x2f matRow = mat2x2fBuffer[1].Row(0);
    Print(matRow);
    printf("Mat 2x1 Float from Column:");
    Mat2x1f matColumn = mat2x2fBuffer[1].Column(0);
    Print(matColumn);
    printf("\nEnd -- Accessing matrix elements using Row and Column views to convert them to Matrixes.\n");

    auto x = 0;
    for(const auto & mat : mat2x2fBuffer)
    {
        printf("\nMatrix %d:", ++x);
        Print(mat);
	}

    printf("\nPointer (Begin): %zu", reinterpret_cast<uintptr_t>(bv.begin()));
    printf("\nPointer (End): %zu", reinterpret_cast<uintptr_t>(bv.end()));
    printf("\nSize with(Pointer): %d\n", static_cast<int>(bv.end() - bv.begin()));
    printf("\nSize: %d\n", bv.size());
    
    bv[0] = 10;
	bv[1] = 20;
    bv[2] = 25;
    bv[3] = 800;
    bv[4] = -40;

    unsigned int i = 0;
    for (const auto& p : bv)
    {
        i++;
        printf("%d) Item: %d\n", i, p);
    }

    i = 0;
    for (const int* ptr = bv.cbegin(); ptr != bv.cend(); ++ptr)
    {
        i++;
        printf("%d) Item: %d\n", i, *ptr);
    }

    Mat3x2f mat2_1f = { {-1.3f, -92.35f}, {2.3f, 7.5f}, {-25.2f, 56.1f} };
    mat2_1f.Column(1) = { -1.3f, -92.35f, 0.569f};
    Print(mat2_1f);

    printf("\nStart\n");
    std::unique_ptr<Mat2x2i> mat2_2i_u = std::make_unique<Mat2x2i>(Mat2x2i::ArrayType{ {2,3}, {4,5} });
    Print(*mat2_2i_u.get());
    printf("Item: %d\n", mat2_2i_u->Data()[0][1]);
    printf("\nEnd\n");
    printf("New test\n");
    Print((Mat2x2f{ {1.0f, 1.5f}, {1.0f, 1.5f} }) * (Mat2x2f{ {1.0f, 1.5f}, {1.0f, 1.5f} }));
    printf("End new test\n");

    printf("\nxxxx==================================================================================================xxx\n");
    Mat4x4d::ArrayType array = { {7, 1,4,70}, {32, 92,-34,33}, {800.23,4,3,2}, {-24,-98,-43,-3} };
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
    Print(mat4f * mat4x4f);

    Mat2x2f mat2x2f_1 = { {-5.1f, 2.3f }, {-8.1f, 4.3f } };
    Mat2x3f mat2x3f_2 = {{ -3.1f, 15.3f, 10.30f }, { -5.1f, 2.3f, 30.54f } };
    printf("\n\nMatrix Multiplication Sample (2 x 3):\n");
    Print(mat2x2f_1 * mat2x3f_2);

    Print(mat2x2f_1);
    Print(mat2x3f_2);

    printf("\n\nPress any key.\n");
    int v = getchar();

    return 0;
}