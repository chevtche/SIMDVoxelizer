
EXAMPLE=SIMDVoxelizer
CPP_SRC=SIMDVoxelizer.cpp
ISPC_SRC=SIMDVoxelizer.ispc
ISPC_IA_TARGETS=sse2-i32x4,sse4-i32x8,avx1-i32x16,avx2-i32x16,avx512knl-i32x16,avx512skx-i32x16
ISPC_ARM_TARGETS=neon

include common.mk
