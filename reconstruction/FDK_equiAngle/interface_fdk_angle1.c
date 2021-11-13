#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PI 3.14159265358979

typedef struct TestStruct
{
    double ScanR;
    double DecFanAng;
    double DecHeigh;
    int YL;
    int ZL;
    double YOffSet;
    double ZOffSet;
    int AngleNumber;
    double DistD;
    double Radius;
    int RecSize;
    int centerX;
    int centerY;
    int centerZ;
    int FOILength;
    int FOIWidth;
    int FOIHeigh;
    double*** GF;
    double*** RecIm;
} TestStruct;

extern void fbp(TestStruct* t)
{

    double ScanR = NAN;
    double DecL = NAN;
    double DecHeigh = NAN;
    double DistD = NAN;
    double Radius = NAN;
    double DeltaY = NAN;
    double DeltaZ = NAN;
    double DeltaL = NAN;
    double YCtr = NAN;
    double ZCtr = NAN;
    double DeltaR = NAN;
    double RCtr = NAN;
    double RadiusSquare = NAN;
    double YOffSet = NAN;
    double ZOffSet = NAN;
    int YL = 0;
    int ZL = 0;
    int ProjNum = 0;
    int RecSize = 0;
    int centerX = 0;
    int centerY = 0;
    int centerZ = 0;
    int FOILength = 0;
    int FOIWidth = 0;
    int FOIHeigh = 0;
    int xl = 0;
    int yl = 0;
    int zl = 0;

    ScanR = t->ScanR;         /*source object distance*/
    DecL = t->DecFanAng;      /*Project Angle*/
    YL = t->YL;               /*Detection Number 672*/
    ProjNum = t->AngleNumber; /*Projection Number 1160*/
    Radius = t->Radius;       /*Radius of the field of view*/
    RecSize = t->RecSize;     /*reconstruction size x*/
    DistD = t->DistD;
    DecHeigh = t->DecHeigh;
    ZL = t->ZL;

    centerX = t->centerX;
    centerY = t->centerY;
    centerZ = t->centerZ;
    FOILength = t->FOILength;
    FOIWidth = t->FOIWidth;
    FOIHeigh = t->FOIHeigh;
    YOffSet = t->YOffSet;
    ZOffSet = t->ZOffSet;

    DeltaY = DecL / YL;     /*Each move Angle of projection*/
    DeltaZ = DecHeigh / ZL; /*Each move Angle of projection*/
    YCtr = (YL - 1) * 0.5 + YOffSet;
    ZCtr = (ZL - 1) * 0.5 + ZOffSet;
    DeltaR = 2 * Radius / RecSize;
    RCtr = (RecSize - 1) * 0.5;
    RadiusSquare = Radius * Radius;
    DeltaL = 2 * PI / ProjNum;

    xl = (int)(centerX - FOILength * 0.5);
    yl = (int)(centerY - FOIWidth * 0.5);
    zl = (int)(centerZ - FOIHeigh * 0.5);

    double* VectorS = NULL;
    double* VectorE = NULL;
    double* xCor = NULL;
    double* yCor = NULL;
    double* zCor = NULL;
    int loop = 0;
    double temp = NAN;

    VectorS = (double*)malloc(sizeof(double) * 2 * ProjNum);
    VectorE = (double*)malloc(sizeof(double) * 2 * ProjNum);
    xCor = (double*)malloc(sizeof(double) * RecSize);
    yCor = (double*)malloc(sizeof(double) * RecSize);
    zCor = (double*)malloc(sizeof(double) * RecSize);

    for (loop = 0; loop < ProjNum; loop++) {
        temp = (loop + 180) * DeltaL;
        VectorS[(ptrdiff_t)loop * 2] = ScanR * cos(temp);
        VectorS[loop * 2 + 1] = ScanR * sin(temp);
        VectorE[(ptrdiff_t)loop * 2] = cos(temp);
        VectorE[loop * 2 + 1] = sin(temp);
    }

    for (loop = 0; loop < RecSize; loop++) {
        xCor[loop] = (loop - RCtr + xl) * DeltaR;
    }

    for (loop = 0; loop < RecSize; loop++) {
        yCor[loop] = (loop - RCtr + yl) * DeltaR;
    }

    for (loop = 0; loop < RecSize; loop++) {
        zCor[loop] = (loop - RCtr + zl) * DeltaR;
    }

    /* argument about object:  ObjR RecMX RecMY RecMZ */

    int i = 0;
    int j = 0;
    int k = 0;
    int ProjIndex = 0;
    int UU = 0;
    int UL = 0;
    int VV = 0;
    int VL = 0;
    double dis = NAN;
    double Dlocal = NAN;
    double DSX[2];
    double UCor = NAN;
    double VCor = NAN;
    double alfa = NAN;
    double beta = NAN;

    for (i = 0; i < FOILength; i++) {
        for (j = 0; j < FOIWidth; j++) {
            for (k = 0; k < FOIHeigh; k++)
            //                for (k=78;k<82;k++)
            {

                if ((xCor[i] * xCor[i] + yCor[j] * yCor[j]) < 2 * RadiusSquare) {
                    t->RecIm[i][j][k] = 0;
                    for (ProjIndex = 0; ProjIndex < ProjNum; ProjIndex++) {
                        DSX[0] = xCor[i] - VectorS[ProjIndex * 2];
                        DSX[1] = yCor[j] - VectorS[ProjIndex * 2 + 1];
                        Dlocal = sqrt(DSX[0] * DSX[0] + DSX[1] * DSX[1]);
                        dis = -fabs(xCor[i] * VectorE[ProjIndex * 2] + yCor[j] * VectorE[ProjIndex * 2 + 1] - ScanR);
                        UCor = atan((DSX[0] * VectorE[ProjIndex * 2 + 1] - DSX[1] * VectorE[ProjIndex * 2]) / dis);
                        VCor = (zCor[k] * DistD) / dis;
                        UCor = UCor / DeltaY + YCtr + 0.5;
                        VCor = VCor / DeltaZ + ZCtr + 0.5;
                        UL = (int)UCor;
                        VL = (int)VCor;
                        UU = UL + 1;
                        VV = VL + 1;
                        alfa = UU - UCor;
                        beta = VV - VCor;

                        if ((UL > 0) && (UU < YL) && (VL > 0) && (VV < ZL)) {
                            t->RecIm[i][j][k] = (t->GF[UU][VL][ProjIndex] * (1 - alfa) * beta +
                                                 t->GF[UU][VV][ProjIndex] * (1 - alfa) * (1 - beta) +
                                                 t->GF[UL][VL][ProjIndex] * alfa * beta +
                                                 t->GF[UL][VV][ProjIndex] * alfa * (1 - beta)) /
                                                    Dlocal +
                                                t->RecIm[i][j][k];
                        }
                    } // for(projindex=0;Projindex<ProjNum;Projindex++)
                    t->RecIm[i][j][k] = -t->RecIm[i][j][k] / (2 * PI * PI * 72);
                } // if(tpdata<0)
            }     // for (k=0;k<RecSize;k++)
        }         // for (j=0;j<RecSize;j++)
    }             // for (i=0;i<RecSize;i++)

    free(VectorS);
    free(VectorE);
    free(xCor);
    free(yCor);
    free(zCor);
}
