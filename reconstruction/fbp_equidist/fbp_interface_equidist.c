#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PI 3.14159265358979

typedef struct TestStruct
{
    double ScanR;
    double DecLength;
    int YL;
    int AngleNumber;
    double DistD;
    double Radius;
    int RecSize;
    int centerX;
    int centerY;
    int FOILength;
    int FOIWidth;
    double** GF;
    double** RecIm;
} TestStruct;


extern void fbp(TestStruct* t)
{
    double ScanR = NAN;
    double DecL = NAN;
    double DistD = NAN;
    double Radius = NAN;
    double DeltaY = NAN;
    double YCtr = NAN;
    double DeltaR = NAN;
    double DeltaL = NAN;
    double RCtr = NAN;
    double RadiusSquare = NAN;
    int YL = 0;
    int AngleNumber = 0;
    int RecSize = 0;
    int FOILength = 0;
    int FOIWidth = 0;

    ScanR = t->ScanR;
    DecL = t->DecLength;
    YL = t->YL;
    AngleNumber = t->AngleNumber;
    DistD = t->DistD;
    Radius = t->Radius;
    RecSize = t->RecSize;

    // centerX = t->centerX;
    // centerY = t->centerY;
    FOILength = t->FOILength;
    FOIWidth = t->FOIWidth;

    DeltaY = DecL / YL;
    YCtr = (YL - 1) * 0.5;
    DeltaR = 2 * Radius / RecSize;
    RCtr = (RecSize - 1) * 0.5;
    RadiusSquare = Radius * Radius;
    DeltaL = 2 * PI / AngleNumber;

    // xl = (int)(centerX - FOILength*0.5);
    // xr = (int)(centerX + FOILength*0.5);
    // yl = (int)(centerY - FOIWidth*0.5);
    // yr = (int)(centerY + FOIWidth*0.5);

    /* Compute the coordinates of scanning locus and its corresponding local coordinate */

    double* VectorS = NULL;
    double* VectorE = NULL;
    double* xCor = NULL;
    double* yCor = NULL;
    int loop = 0;
    double temp = NAN;

    VectorS = (double*)malloc(sizeof(double) * 2 * AngleNumber);
    VectorE = (double*)malloc(sizeof(double) * 2 * AngleNumber);
    xCor = (double*)malloc(sizeof(double) * RecSize);
    yCor = (double*)malloc(sizeof(double) * RecSize);

    for (loop = 0; loop < AngleNumber; loop++) {
        temp = (loop + 180) * DeltaL;
        VectorS[(ptrdiff_t)loop * 2] = ScanR * cos(temp);
        VectorS[loop * 2 + 1] = ScanR * sin(temp);
        VectorE[(ptrdiff_t)loop * 2] = cos(temp);
        VectorE[loop * 2 + 1] = sin(temp);
    }

    for (loop = 0; loop < RecSize; loop++) {
        xCor[loop] = (loop - RCtr) * DeltaR;
    }

    for (loop = 0; loop < RecSize; loop++) {
        yCor[loop] = (loop - RCtr) * DeltaR;
    }

    /* argument about object:  ObjR RecMX RecMY RecMZ*/

    int j = 0;
    int i = 0;
    int ProjIndex = 0;
    int UU = 0;
    int UL = 0;
    double dis = NAN;
    double Dlocal = NAN;
    double m1 = NAN;
    double m2 = NAN;
    double UCor = NAN;
    double alfa = NAN;
    double W2 = NAN;

    for (i = 0; i < FOILength; i++) {
        for (j = 0; j < FOIWidth; j++) {
            if ((xCor[i] * xCor[i] + yCor[j] * yCor[j]) < 2 * RadiusSquare) {
                t->RecIm[i][j] = 0;
                for (ProjIndex = 0; ProjIndex < AngleNumber; ProjIndex++) {
                    m1 = xCor[i] - VectorS[(ptrdiff_t)ProjIndex * 2];
                    m2 = yCor[j] - VectorS[ProjIndex * 2 + 1];
                    Dlocal = sqrt(m1 * m1 + m2 * m2);
                    dis = fabs(xCor[i] * VectorE[(ptrdiff_t)ProjIndex * 2] + yCor[j] * VectorE[ProjIndex * 2 + 1] -
                               ScanR);
                    UCor = (m1 * VectorE[ProjIndex * 2 + 1] - m2 * VectorE[(ptrdiff_t)ProjIndex * 2]) * DistD / dis;

                    W2 = sqrt(DistD * DistD + UCor * UCor) / Dlocal;
                    UCor = UCor / DeltaY + YCtr + 0;
                    UL = (int)UCor;
                    UU = UL + 1;
                    alfa = UU - UCor;

                    if ((UL > 0) && (UU < YL)) {
                        t->RecIm[i][j] =
                            (t->GF[UU][ProjIndex] * (1 - alfa) + t->GF[UL][ProjIndex] * alfa) * W2 + t->RecIm[i][j];
                    }
                }
                t->RecIm[i][j] = t->RecIm[i][j] / (2 * AngleNumber);
            }
        }
    }

    free(VectorS);
    free(VectorE);
    free(xCor);
    free(yCor);
}
