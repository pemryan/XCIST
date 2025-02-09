import ctypes as ct
import numpy as np
import scipy.io as scio
import matplotlib.pyplot as plt

# Init ctypes types
DOUBLE = ct.c_double
PtrDOUBLE = ct.POINTER(DOUBLE)
PtrPtrDOUBLE = ct.POINTER(PtrDOUBLE)


class TestStruct(ct.Structure):
    _fields_ = [
        ("ScanR", ct.c_double),
        ("DecFanAng", ct.c_double),
        ("YL", ct.c_int),
        ("AngleNumber", ct.c_int),
        ("Radius", ct.c_double),
        ("RecSizeX", ct.c_int),
        ("RecSizeY", ct.c_int),
        ("centerX", ct.c_int),
        ("centerY", ct.c_int),
        ("FOILength", ct.c_int),
        ("FOIWidth", ct.c_int),
        ("GF", PtrPtrDOUBLE),
        ("RecIm", PtrPtrDOUBLE),
    ]


def double2darray2pointer(array):
    # Converts a 2D numpy into a array ctypes 2D array.
    arr_dimx = DOUBLE * array.shape[1]
    arr_dimy = PtrDOUBLE * array.shape[0]
    arr_ptr = arr_dimy()

    for i, row in enumerate(array):
        arr_ptr[i] = arr_dimx()
        for j, val in enumerate(row):
            arr_ptr[i][j] = val

    return arr_ptr


def double2dpointer2array(ptr, n, m):
    # Converts ctypes 2D array into a 2D numpy array.
    arr = np.zeros(shape=(n, m))
    for i in range(n):
        for j in range(m):
            arr[i][j] = ptr[i][j]
    return arr


# Load the compiled library
recon = ct.CDLL("./fbp_equiAngular.dll")
# Define arguments of the C function
recon.fbp.argtypes = [ct.POINTER(TestStruct)]
# Define the return type of the C function
recon.fbp.restype = None

# Load the data
dataFile = "./data/Res_Filtering_Angle.mat"
data = scio.loadmat(dataFile)

# init the struct
t = TestStruct()

t.ScanR = data["ScanR"]
t.DecFanAng = data["DecFanAng"]
t.YL = data["YL"]
t.AngleNumber = len(data["GF"])
t.Radius = data["Radius"]

# These are flexible parameters.
t.RecSizeX = 256
t.RecSizeY = 256
t.centerX = 128
t.centerY = 128
t.FOILength = 256
t.FOIWidth = 256

# Generate a 2D ctypes array from numpy array
GF = data["GF"]
GF = GF.T
GF_ptr = double2darray2pointer(GF)
t.GF = GF_ptr

RecIm = np.zeros(shape=(t.FOILength, t.FOIWidth))
RecIm_ptr = double2darray2pointer(RecIm)
t.RecIm = RecIm_ptr

# interface with C function
recon.fbp(ct.byref(t))

# Convert ctypes 2D arrays to numpy arrays
RecA = double2dpointer2array(RecIm_ptr, *RecIm.shape)
RecA = RecA.T

# save result
dataNew = "./data/Res_equiAngular.mat"
scio.savemat(dataNew, {"Rec": RecA})
plt.figure()
plt.imshow(RecA, cmap="gray")
plt.show()
