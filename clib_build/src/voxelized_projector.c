// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

// 1. Enabled: materials can have different dimensions, to save memory and computing time.
// 2. Added: free phantom memory at the last projection.
// 3. Bug fixed: memory not cleaned up during view loop (xds, yds, zds).
// Mingye Wu, Nov 8 2020

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <stdarg.h>
#ifndef WIN32
#include <pthread.h>
#include <unistd.h>
#endif

#include "DD3_roi_notrans_mm.h"
#include "getMemorySize.h"
#include "voxelized_projector.h"

#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) > (b)) ? (b) : (a))
#define VERY_BIG 1e300

// #define DEBUG
#if defined(DEBUG)

#define DEBUG_00 // Prints "where am I" info
#define DEBUG_10 // Prints info in convert_modular_detector
// #define DEBUG_20 // Prints info in voxelized_projector

#endif

// Global Variable Definitions

static char TempString[10000];
static char OutputString[10000];
const int PrintReportOutput = 1;

static module_info modules = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, 0, 0 };
static material_info materials = { 0, 0, NULL };
const static int debug_flag = 0;
static phantom_info phantom = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL };

//////////////////////////////////////////////
static void dbug(int level, const char* s, ...);
static int* my_memcpyi(int* src, int* dest, int bytes);
static float* my_memcpyf(float* src, float* dest, int bytes);

static void Report(int Unused)
{

    if (PrintReportOutput) {
        sprintf(TempString, "in C> %s", OutputString);
        fprintf(stdout, "%s", TempString);
        fflush(stdout);
    }
}


static void dbug(int level, const char* s, ...)
{
    if (debug_flag >= level) {
        va_list ap = NULL;
        va_start(ap, s);
        vprintf(s, ap);
        va_end(ap);
        fflush(stdout);
    }
}


static int* my_memcpyi(int* src, int* dest, int bytes)
{
    if (dest != NULL) {
        free(dest);
        dest = NULL;
    }

    dest = (int*)malloc(bytes);
    memcpy(dest, src, bytes);

    return dest;
}


static float* my_memcpyf(float* src, float* dest, int bytes)
{
    if (dest != NULL) {
        free(dest);
        dest = NULL;
    }

    dest = (float*)malloc(bytes);
    memcpy(dest, src, bytes);

    return dest;
}


DLLEXPORT void set_src_info_vox(float* sourceWeights, int nSubSources)
{
#if defined(DEBUG_00)
    Report(sprintf(OutputString, "In set_src_info_vox\n"));
#endif

    modules.sourceWeights = my_memcpyf(sourceWeights, modules.sourceWeights, nSubSources * sizeof(float));
    modules.nSubSources = nSubSources;
}


DLLEXPORT void set_module_info_vox(float* Height,
                                   float* Width,
                                   int* Pix,
                                   float* Coords,
                                   int* Sub,
                                   float* Sampling,
                                   float* Weight,
                                   int nModuleTypes,
                                   int maxPix,
                                   int maxSubDets,
                                   int moduleOverlapType)
{
    int i = 0;

#if defined(DEBUG_00)
    Report(sprintf(OutputString, "In set_module_info_vox\n"));
#endif

    modules.Height = my_memcpyf(Height, modules.Height, nModuleTypes * (int)sizeof(float));
    modules.Width = my_memcpyf(Width, modules.Width, nModuleTypes * (int)sizeof(float));
    for (i = 0; i < nModuleTypes; i++) {
        modules.Height[i] = MAX(1e-7, modules.Height[i]);
        modules.Width[i] = MAX(1e-7, modules.Width[i]);
    }

    modules.Pix = my_memcpyi(Pix, modules.Pix, nModuleTypes * (int)sizeof(int));
    modules.Coords = my_memcpyf(Coords, modules.Coords, maxPix * 2 * nModuleTypes * (int)sizeof(float));
    modules.Sub = my_memcpyi(Sub, modules.Sub, nModuleTypes * (int)sizeof(int));
    modules.Sampling = my_memcpyf(Sampling, modules.Sampling, 2 * maxSubDets * nModuleTypes * (int)sizeof(float));
    modules.Weight = my_memcpyf(Weight, modules.Weight, maxSubDets * nModuleTypes * (int)sizeof(float));
    modules.maxPixPerModule = maxPix;
    modules.maxSubDets = maxSubDets;
    modules.moduleOverlapType = moduleOverlapType;
    modules.nModuleTypes = nModuleTypes;
}


DLLEXPORT void set_material_info_vox(int materialCount, int eBinCount, float* muTable)
{
#if defined(DEBUG_00)
    Report(sprintf(OutputString, "In set_material_info_vox\n"));
#endif

    materials.materialCount = materialCount;
    materials.eBinCount = eBinCount;
    materials.muTable = my_memcpyf(muTable, materials.muTable, eBinCount * materialCount * (int)sizeof(float));
}


DLLEXPORT void set_phantom_info_vox(int* Status,
                                    float* vol,
                                    int* dims,
                                    float xoff,
                                    float yoff,
                                    float zoff,
                                    float dxy,
                                    float dz,
                                    unsigned char* xy_mask,
                                    int MaterialIndex,
                                    int NumOfMaterials)
{
    static uint64_t previously_allocated_memory_size;
    static uint64_t system_memory_size;
    uint64_t phantom_num_voxels_this_material = (uint64_t)dims[0] * (uint64_t)dims[1] * (uint64_t)dims[2];
    uint64_t required_memory_size_this_material = phantom_num_voxels_this_material * (uint64_t)sizeof(float);
    const uint64_t reserved_memory_size =
        (uint64_t)2 * 1024 * 1024 * 1024; // reserve 2GB memory to store sinogram and other things;

    *Status = 0;

    if (phantom.volume == NULL) {
        previously_allocated_memory_size = 0;
        system_memory_size = getMemorySize();

        Report(sprintf(OutputString, "Preparing to allocate memory for material volume data...\n"));
        phantom.volume = (float**)malloc(sizeof(float*) * NumOfMaterials);
        phantom.dims = (int**)malloc(sizeof(int*) * NumOfMaterials);
        phantom.xy_mask = (unsigned char**)malloc(sizeof(unsigned char*) * NumOfMaterials);

        phantom.xoff = (float*)malloc(sizeof(float) * NumOfMaterials);
        phantom.yoff = (float*)malloc(sizeof(float) * NumOfMaterials);
        phantom.zoff = (float*)malloc(sizeof(float) * NumOfMaterials);
        phantom.dxy = (float*)malloc(sizeof(float) * NumOfMaterials);
        phantom.dz = (float*)malloc(sizeof(float) * NumOfMaterials);

        if (phantom.volume == NULL || phantom.dims == NULL || phantom.xy_mask == NULL) {
            Report(sprintf(OutputString,
                           "Memory allocation error - couldn't allocate memory for pointers to materials.\n"));
            *Status = -1;
            return;
        }
    }

    if (system_memory_size - reserved_memory_size <
        previously_allocated_memory_size + required_memory_size_this_material) {
        Report(sprintf(OutputString, "Insuffucient system memory available.\n"));

        *Status = -2;
        return;
    }
    else {
        (phantom.volume)[MaterialIndex - 1] = (float*)malloc(required_memory_size_this_material);
        (phantom.dims)[MaterialIndex - 1] = (int*)malloc(4 * sizeof(int));
        phantom.xy_mask[MaterialIndex - 1] =
            (unsigned char*)malloc((size_t)dims[0] * dims[1] * 2 * sizeof(unsigned char));

        if (phantom.volume[MaterialIndex - 1] == NULL || phantom.dims[MaterialIndex - 1] == NULL) {
            Report(sprintf(
                OutputString, "Memory allocation error - couldn't allocate memory for material %i.\n", MaterialIndex));
            *Status = -1;
            return;
        }

        Report(sprintf(OutputString, "Allocated memory for image volume for material %2i\n", MaterialIndex));
        previously_allocated_memory_size += required_memory_size_this_material;
    }

    Report(sprintf(OutputString, "Copying data for material %2d into C memory...", MaterialIndex));
    memcpy((phantom.volume)[MaterialIndex - 1], vol, dims[0] * dims[1] * dims[2] * sizeof(float));
    Report(sprintf(OutputString, "done.\n"));

    memcpy(phantom.dims[MaterialIndex - 1], dims, 4 * sizeof(int));
    phantom.xoff[MaterialIndex - 1] = xoff;
    phantom.yoff[MaterialIndex - 1] = yoff;
    phantom.zoff[MaterialIndex - 1] = zoff;
    phantom.dxy[MaterialIndex - 1] = dxy;
    phantom.dz[MaterialIndex - 1] = dz;

    // phantom.xy_mask[MaterialIndex-1] = new unsigned char[dims[0]*dims[1]*2];
    memcpy(phantom.xy_mask[MaterialIndex - 1], xy_mask, sizeof(unsigned char) * dims[0] * dims[1]);

#if defined(DEBUG_00)
    Report(sprintf(OutputString, "Copied mask into C memory\n"));
#endif

    unsigned char* transposeMask = phantom.xy_mask[MaterialIndex - 1] + (ptrdiff_t)dims[0] * dims[1];
    for (int i = 0; i < dims[0]; i++) {
        for (int j = 0; j < dims[1]; j++) {
            transposeMask[j + i * dims[1]] = xy_mask[i + j * dims[0]];
#if defined(DEBUG_00)
            Report(sprintf(OutputString, "Transposed mask\n"));
#endif
        }
    }

    if (MaterialIndex == NumOfMaterials) {
        Report(sprintf(OutputString,
                       "Allocated a total of %6llu MB.\n",
                       previously_allocated_memory_size / ((uint32_t)(1024 * 1024))));
    }
}


int convert_modular_detector(float** xds,
                             float** yds,
                             float** zds,
                             int* nrdetcols,
                             int* nrdetrows,
                             int nModulesIn,
                             int* modTypeInds,
                             float* Up,
                             float* Right,
                             float* Center)
{
    float tiny = 1e-8F;
    float firstU = NAN;
    float thisU = NAN;
    float tol = 1e-6F;
    float* coords = NULL;
    int ModTypeIndex = 0;
    int ModIndex = 0;
    int ColIndex = 0;
    int RowIndex = 0;
    int NumPixels = 0;
    int thisNumRows = 0;
    int thisNumCols = 0;
    float* up = NULL;
    float* center = NULL;
    float* right = NULL;
    float inPlane = NAN;
    float inZ = NAN;

#if defined(DEBUG_00)
    Report(sprintf(OutputString, "In convert_modular_detector\n"));
#endif

    // first we verify that each module:
    // 1) has columns oriented in z and rows oriented orth. to z
    // 2) has (u,v) coords that are separable into an xy component and a z component (i.e., a grid)
    // 3) has the same number of rows

    *nrdetrows = 0;
    *nrdetcols = 0;
    // 1) has columns oriented in z and rows oriented orth. to z
    for (int ModIndex = 0; ModIndex < nModulesIn; ModIndex++) {
#if defined(DEBUG_10)
        Report(sprintf(OutputString, "Evaluating module %i of %i for correct orientation\n", ModIndex + 1, nModulesIn));
#endif
        up = Up + (ptrdiff_t)3 * ModIndex;
        right = Right + (ptrdiff_t)3 * ModIndex;
        inPlane = up[0] * up[0] + up[1] * up[1];
        inZ = up[2] * up[2];
#if defined(DEBUG_10)
        Report(sprintf(OutputString, "up = %f, right = %f, inPlane = %f, inZ = %f\n", *up, *right, inPlane, inZ));
#endif
        if ((inPlane > inZ) || (inPlane / inZ > tiny)) {
            Report(
                sprintf(OutputString,
                        "ERROR: Modules columns must be parallel to the z direction for the voxelized projector.\n"));
            // dbug(-1,"\r\nERROR: Modules columns must be parallel to the z direction for the voxelized
            // projector.\r\n");
            return (-2);
        }

        inPlane = right[0] * right[0] + right[1] * right[1];
        inZ = right[2] * right[2];
#if defined(DEBUG_10)
        Report(sprintf(OutputString, "inPlane = %f, inZ = %f\n", inPlane, inZ));
#endif
        if ((inZ > inPlane) || (inZ / inPlane > tiny)) {
            Report(sprintf(OutputString,
                           "ERROR: Modules rows must be orthogonal to the z direction for the voxelized projector.\n"));
            // dbug(-1,"\r\nERROR: Modules rows must be orthogonal to the z direction for the voxelized
            // projector.\r\n");
            return (-2);
        }
    }

    // 2) has (u,v) coords that are separable into an xy component and a z component (i.e., a grid)
    for (ModTypeIndex = 0; ModTypeIndex < modules.nModuleTypes; ModTypeIndex++) {
        NumPixels = modules.Pix[ModTypeIndex];
        coords = modules.Coords + (ptrdiff_t)ModTypeIndex * modules.maxPixPerModule * 2;

#if defined(DEBUG_10)
        Report(sprintf(OutputString,
                       "Evaluating (u,v) coordinates of module type %i of %i\n",
                       ModTypeIndex + 1,
                       modules.nModuleTypes));
        Report(sprintf(OutputString, "Module type %i has %i pixels.\n", ModTypeIndex + 1, modules.Pix[ModTypeIndex]));
        Report(sprintf(OutputString, "coords[0 1 2 3]: %f %f %f %f \r\n", coords[0], coords[1], coords[2], coords[3]));
#endif
        // dbug(2,"coords[0 1 2 3]: %f %f %f %f \r\n",coords[0],coords[1],coords[2],coords[3]);

        thisNumRows = 1;
        firstU = coords[0];
        thisU =
            coords[(ptrdiff_t)2 * thisNumRows]; // this will be bogus if only one pixel in module but that's ok because
                                                // the loop below will terminate because it's the last pixel.
        // if we're on the same column (same U) and we haven't reached the last pixel, ...
        while ((fabsf(thisU - firstU) < tol) && (thisNumRows != NumPixels)) {
            // increment the total number of rows,
            thisNumRows++;
            // and next time through the loop, check the next pixel.
            thisU = coords[(ptrdiff_t)2 * thisNumRows];
        }

        thisNumCols = modules.Pix[ModTypeIndex] / thisNumRows;
#if defined(DEBUG_10)
        Report(sprintf(OutputString, "firstU = %f, thisU = %f\n", firstU, thisU));
#endif
// dbug(2,"thisU: %f firstU: %f\r\n",thisU,firstU);
#if defined(DEBUG_10)
        Report(sprintf(OutputString, "thisNumRows = %i, thisNumCols = %i\n", thisNumRows, thisNumCols));
#endif
        // dbug(2,"pix: %d this_rows: %d\r\n",modules.Pix[ModTypeIndex],thisNumRows);

        // 3) has the same number of rows
        if (ModTypeIndex == 0) {
            *nrdetrows = thisNumRows;
        }
        else if (*nrdetrows != thisNumRows) {
            Report(sprintf(
                OutputString,
                "ERROR: All module types must have the same number of rows for voxelized projector (%d, %d, %d).\n",
                *nrdetrows,
                thisNumRows,
                ModTypeIndex));
            // dbug(-1,"\r\nERROR: All module types must have the same number of rows for voxelized projector (%d, %d,
            // %d).\r\n",*nrdetrows,thisNumRows,ModTypeIndex);
            return (-2);
        }

        // dbug(2,"pix: %d thisNumRows: %d thisNumCols: %d\r\n",modules.Pix[ModTypeIndex],thisNumRows,thisNumCols);

        if ((thisNumCols * thisNumRows) != modules.Pix[ModTypeIndex]) {
            Report(
                sprintf(OutputString,
                        "ERROR: All module types must be rectilinear grid pixel sampling for voxelized projector.\n"));
            // dbug(-1,"\r\nERROR: All module types must be rectilinear grid pixel sampling for voxelized projector
            // (errorcode 1).\r\n");
            return (-2);
        }

        for (ColIndex = 0; ColIndex < thisNumCols; ColIndex++) {
            for (RowIndex = 0; RowIndex < thisNumRows; RowIndex++) {
                if ((fabsf(coords[2 * (RowIndex + ColIndex * thisNumRows) + 1] -
                           coords[2 * (RowIndex + 0 * thisNumRows) + 1]) > tol) ||
                    (fabsf(coords[(ptrdiff_t)2 * (RowIndex + ColIndex * thisNumRows)] -
                           coords[(ptrdiff_t)2 * (0 + ColIndex * thisNumRows)]) > tol)) {
                    Report(sprintf(OutputString, "tol: %f \n", tol));
                    Report(sprintf(OutputString, "ColIndex: %d RowIndex: %d \r\n", ColIndex, RowIndex));
                    Report(sprintf(
                        OutputString,
                        "coords[2*(RowIndex+ColIndex*thisNumRows)+1]: %f coords[2*(RowIndex+0*thisNumRows)+1]: %f \r\n",
                        coords[2 * (RowIndex + ColIndex * thisNumRows) + 1],
                        coords[2 * (RowIndex + 0 * thisNumRows) + 1]));
                    Report(sprintf(
                        OutputString,
                        "coords[2*(RowIndex+ColIndex*thisNumRows)]:   %f coords[2*(0+ColIndex*thisNumRows)]:   %f \r\n",
                        coords[(ptrdiff_t)2 * (RowIndex + ColIndex * thisNumRows)],
                        coords[(ptrdiff_t)2 * (0 + ColIndex * thisNumRows)]));
                    Report(sprintf(
                        OutputString,
                        "ERROR: All module types must be rectilinear grid pixel sampling for voxelized projector.\n"));
                    return (-2);
                }
            }
        }
    }

    for (ModIndex = 0; ModIndex < nModulesIn; ModIndex++) {
        thisNumCols = modules.Pix[modTypeInds[ModIndex]] / (*nrdetrows);
        *nrdetcols += thisNumCols;
    }
#if defined(DEBUG_10)
    Report(sprintf(OutputString, "*nrdetcols = %i, *nrdetrows = %i\n", *nrdetcols, *nrdetrows));
#endif

    *xds = (float*)malloc(sizeof(float) * (*nrdetcols));
    *yds = (float*)malloc(sizeof(float) * (*nrdetcols));
    *zds = (float*)malloc(sizeof(float) * (*nrdetrows));

    // initial verifications complete... carry on
    // we also verify here that all modules have the same z locations for the rows.
    *nrdetcols = 0;
    for (ModIndex = 0; ModIndex < nModulesIn; ModIndex++) {
#if defined(DEBUG_10)
        Report(sprintf(
            OutputString, "Assigning system coordinates to pixels of module %i of %i\n", ModIndex + 1, nModulesIn));
#endif
        up = Up + (ptrdiff_t)3 * ModIndex;
        center = Center + (ptrdiff_t)3 * ModIndex;
        right = Right + (ptrdiff_t)3 * ModIndex;
#if defined(DEBUG_10)
        Report(sprintf(OutputString, "up = %f, center = %f, right = %f\n", *up, *center, *right));
#endif

        thisNumCols = modules.Pix[modTypeInds[ModIndex]] / (*nrdetrows);
        coords = modules.Coords + (ptrdiff_t)modTypeInds[ModIndex] * modules.maxPixPerModule * 2;
        if (ModIndex == 0) {
            for (RowIndex = 0; RowIndex < (*nrdetrows); RowIndex++) {
                (*zds)[RowIndex] = center[2] + coords[1 + RowIndex * 2] * up[2];
            }
        }
        else {
            for (RowIndex = 0; RowIndex < (*nrdetrows); RowIndex++) {
                if (fabsf((*zds)[RowIndex] - (center[2] + coords[1 + RowIndex * 2] * up[2])) > tol) {
                    Report(sprintf(
                        OutputString,
                        "nzds[RowIndex] = %f, center[2] = %f, coords[1+RowIndex*2] = %f, up[2] = %f, tol = %f\n",
                        (*zds)[RowIndex],
                        center[2],
                        coords[1 + RowIndex * 2],
                        up[2],
                        tol));
                    Report(sprintf(OutputString,
                                   "ERROR: All modules must line up in z direction (ModIndex = %d, RowIndex = %d).\n",
                                   ModIndex,
                                   RowIndex));
                    return (-2);
                    // dbug(0,"\r\nzds[RowIndex]: %f   center[2]: %f   coords[1+RowIndex*2]: %f   up[2]: %f  tol:
                    // %f\r\n",(*zds)[RowIndex],center[2],coords[1+RowIndex*2],up[2],tol); dbug(-1,"\r\nERROR: All
                    // modules must line up in z direction (ModIndex: %d, RowIndex: %d).\r\n",ModIndex,RowIndex);
                }
            }
        }

        for (ColIndex = 0; ColIndex < thisNumCols; ColIndex++) {
            (*xds)[ColIndex + *nrdetcols] = center[0] + coords[(ptrdiff_t)ColIndex * 2 * thisNumRows] * right[0];
            (*yds)[ColIndex + *nrdetcols] = center[1] + coords[(ptrdiff_t)ColIndex * 2 * thisNumRows] * right[1];
        }

        *nrdetcols += thisNumCols;
    }

#if defined(DEBUG_10)
    Report(sprintf(OutputString, "*nrdetcols = %i, *nrdetrows = %i\n", *nrdetcols, *nrdetrows));
    Report(sprintf(OutputString, "Pixel(%3s,%3s) is at (%8s,%8s,%8s):\n", "col", "row", "x", "y", "z"));
    Report(sprintf(
        OutputString, "Pixel(%3i,%3i) is at (%8.3f,%8.3f,%8.3f)\n", 1, 1, *xds[1 - 1], *yds[1 - 1], *zds[1 - 1]));
    Report(sprintf(OutputString,
                   "Pixel(%3i,%3i) is at (%8.3f,%8.3f,%8.3f)\n",
                   *nrdetcols,
                   *nrdetrows,
                   *xds[*nrdetcols - 1],
                   *yds[*nrdetcols - 1],
                   *zds[*nrdetrows - 1]));
#endif

    // dbug(1,"xds[e]: %f\n\r\n",(*xds)[100]);

#if defined(DEBUG_00)
    Report(sprintf(OutputString, "Returning from convert_modular_detector\n"));
#endif

    free(xds);
    free(yds);
    free(zds);

    return (0);
}


DLLEXPORT
void voxelized_projector(int* Status,               // [out] scalar. 0=normal, 1=Detector definition error
                         float unused1,             // [in] scalar
                         float* thisView,           // [out] [nrdetcols*nrdetrows][materials.eBinCount]
                         const float* sourcePoints, // [in] [nSubSources][3]
                         int nSubSources,           // [in] scalar
                         float* unused2,            // [in] pointer
                         int unused3,               // [in] scalar
                         int* unused4,              // [in] pointer
                         int nModulesIn,            // [in] scalar
                         int* modTypeInds,          // [in] pointer
                         float* Up,                 // [in] pointer
                         float* Right,              // [in] pointer
                         float* Center,             // [in] pointer
                         int unused5,               // [in] scalar
                         int unused6,               // [in] scalar
                         int MaterialIndex,         // [in] scalar, Material Index to look up mu table
                         int MaterialIndexInMemory, // [in] scalar, Material Index to materials stored in phantom.volume
                         float unused7,             // [in] scalar
                         int freeTheMemory)         // [in] scalar, free and clean up memory if it's 1
{
    float* xds = NULL;
    float* yds = NULL;
    float* zds = NULL;
    float x0 = 0.F;
    float y0 = 0.F;
    float z0 = 0.F;
    float totalWt = 0.F;
    float viewangle = NAN;
    int nrdetcols = 0;
    int nrdetrows = 0;

    int EnergyBin = 0;

    *Status = 0;

#if defined(DEBUG_00)
    Report(sprintf(OutputString, "In voxelized_projector\n"));
#endif

    // condense source to one point at center of mass
    for (int i = 0; i < nSubSources; i++) {
        x0 += modules.sourceWeights[i] * sourcePoints[0];
        y0 += modules.sourceWeights[i] * sourcePoints[1];
        z0 += modules.sourceWeights[i] * sourcePoints[2];
        totalWt += modules.sourceWeights[i];
        sourcePoints += 3;
    }

    x0 /= totalWt;
    y0 /= totalWt;
    z0 /= totalWt;

    // viewangle = atan2f(-x0,y0);
    viewangle = 0;

    *Status =
        convert_modular_detector(&xds, &yds, &zds, &nrdetcols, &nrdetrows, nModulesIn, modTypeInds, Up, Right, Center);

    if (0 != (*Status)) {
        Report(sprintf(OutputString, "Error code %i returned by convert_modular_detector\n", *Status));
        return;
    }

    float* thisViewProj = (float*)malloc(sizeof(float) * nrdetcols * nrdetrows);

    dbug(1, "xds[e]: %f\n\r\n", xds[100]);
    dbug(1, "materials.muTable[0]: %f\r\n", materials.muTable[0]);

#if defined(DEBUG_00)
    Report(sprintf(OutputString, "Calling DD3Proj_roi_notrans_mm\n"));
#endif

    float zoffset = 0.0F;
    int m = MaterialIndexInMemory - 1;

    DD3Proj_roi_notrans_mm(x0,
                           y0,
                           z0,
                           nrdetcols,
                           nrdetrows,
                           xds,
                           yds,
                           zds,
                           phantom.xoff[m],
                           phantom.yoff[m],
                           phantom.zoff[m],
                           &viewangle,
                           &zoffset,
                           1,
                           thisViewProj,
                           phantom.dims[m][0],
                           phantom.dims[m][1],
                           phantom.dims[m][2],
                           (phantom.volume)[m],
                           phantom.dxy[m],
                           phantom.dz[m],
                           phantom.xy_mask[m]);

#if defined(DEBUG_20)
    Report(sprintf(OutputString, "materials.muTable[0] = %f\n", materials.muTable[0]));
    Report(sprintf(
        OutputString, "thisView[center=%d] = %f\n", nrdetcols * nrdetrows / 2, thisView[nrdetcols * nrdetrows / 2]));
#endif

    for (int detIndex = 0; detIndex < nrdetcols * nrdetrows; detIndex++) {
        for (EnergyBin = 0; EnergyBin < materials.eBinCount; EnergyBin++) {
            thisView[detIndex * materials.eBinCount + EnergyBin] =
                thisViewProj[detIndex] * materials.muTable[EnergyBin * materials.materialCount + (MaterialIndex - 1)];

#if defined(DEBUG_20)
            Report(sprintf(OutputString, "thisView(index1)=%f\n", thisView[detIndex]));
            Report(sprintf(OutputString,
                           "materials.muTable(index2) = %f\n",
                           materials.muTable[EnergyBin * materials.materialCount + (MaterialIndex - 1)]));
#endif
        }
    }

#if defined(DEBUG_00)
    Report(sprintf(OutputString, "Returning from voxelized_projector\n"));
#endif

    free(thisViewProj);
    free(xds); // xds, yds, zds are newed in calling convert_modular_detector
    free(yds);
    free(zds);

    // Free phantom memory
    if (freeTheMemory == 1) {
        for (int i = 0; i < MaterialIndexInMemory; i++) {
            free(phantom.volume[i]);
            free(phantom.dims[i]);
            free(phantom.xy_mask[i]);
        }

        free(phantom.xoff);
        free(phantom.yoff);
        free(phantom.zoff);
        free(phantom.dxy);
        free(phantom.dz);

        phantom.dims = NULL;
        phantom.volume = NULL;
        phantom.xoff = NULL;
        phantom.yoff = NULL;
        phantom.zoff = NULL;
        phantom.dxy = NULL;
        phantom.dz = NULL;
        phantom.xy_mask = NULL;
    }
}
