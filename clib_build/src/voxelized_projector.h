// Copyright 2020, General Electric Company. All rights reserved. See https://github.com/xcist/code/blob/master/LICENSE

// 1. Enabled: materials can have different dimensions, to save memory and computing time.
// 2. Added: free phantom memory at the last projection.
// 3. Bug fixed: memory not cleaned up during view loop (xds, yds, zds).
// Mingye Wu, Nov 8 2020

#ifndef VOXELIZED_PROJECTOR_H
#define VOXELIZED_PROJECTOR_H

#ifdef WIN32
#define DLLEXPORT __declspec(dllexport)
#else
#define DLLEXPORT
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    float* Height;
    float* Width;
    int* Pix;
    float* Coords;
    int* Sub;
    float* Sampling;
    float* Weight;
    float* sourceWeights;
    int nSubSources;
    int maxPixPerModule;
    int maxSubDets;
    int moduleOverlapType;
    int nModuleTypes;
} module_info;

typedef struct
{
    int materialCount;
    int eBinCount;
    float* muTable;
} material_info;

typedef struct
{
    int** dims;
    float** volume; // volume is a pointer to pointers because the voxelized phantom may have multiple materials
    float* xoff;
    float* yoff;
    float* zoff;
    float* dxy;
    float* dz;
    unsigned char** xy_mask;
} phantom_info;

//////////////////////////////////////////////

DLLEXPORT void set_src_info_vox(float* sourceWeights, int nSubSources);

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
                                   int moduleOverlapType);

DLLEXPORT void set_material_info_vox(int materialCount, int eBinCount, float* muTable);

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
                                    int NumOfMaterials);

int convert_modular_detector(float** xds,
                             float** yds,
                             float** zds,
                             int* nrdetcols,
                             int* nrdetrows,
                             int nModulesIn,
                             int* modTypeInds,
                             float* Up,
                             float* Right,
                             float* Center);

DLLEXPORT void voxelized_projector(
    int* Status,               // [out] scalar. 0=normal, 1=Detector definition error
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
    int freeTheMemory);        // [in] scalar, free and clean up memory if it's 1


#ifdef __cplusplus
}
#endif

#endif