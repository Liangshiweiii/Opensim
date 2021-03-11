#include <fenv.h>

#define OP_LOOP_INTERIOR 0

#define OP_LOOP_WHOLE 1

#define STORAGE_LOOP_BEGIN(i,j,k,T__,op_loop_bcells__) \
{\
    const int op_loop_lower_X__ = std::min((long int)(-(op_loop_bcells__)),(long int)0); \
    const int op_loop_lower_Y__ = std::min((long int)(-(op_loop_bcells__)),(long int)0); \
    const int op_loop_lower_Z__ = std::min((long int)(-(op_loop_bcells__)),(long int)0); \
    if ( op_loop_bcells__ > (long int)((T__).Bcells() )){ std::cerr << "BOUNDARY TOO SMALL! PLEASE ADJUST! BCELLS NEEDED " << op_loop_bcells__ << std::endl; throw std::invalid_argument("Bcells");}\
    const int op_loop_upper_X__ = std::max((long int)((T__).sizeX() + (op_loop_bcells__)),(long int)((T__).sizeX())); \
    const int op_loop_upper_Y__ = std::max((long int)((T__).sizeY() + (op_loop_bcells__)),(long int)((T__).sizeY())); \
    const int op_loop_upper_Z__ = std::max((long int)((T__).sizeZ() + (op_loop_bcells__)),(long int)((T__).sizeZ())); \
    for (long int i = op_loop_lower_X__; i < op_loop_upper_X__; ++i) \
    for (long int j = op_loop_lower_Y__; j < op_loop_upper_Y__; ++j) \
    for (long int k = op_loop_lower_Z__; k < op_loop_upper_Z__; ++k) \
    {

#define STORAGE_LOOP_END \
    } \
}          //使用宏初始化结构体

#define STRINGIFY(...) #__VA_ARGS__

#define OMP_PARALLEL_FOR_LOOP_BEGIN(i,begin__,end__,...) \
{\
    _Pragma(STRINGIFY(omp parallel for schedule(OMP_SCHEDULING_TYPE,OMP_DYNAMIC_CHUNKSIZE) __VA_ARGS__) ) \
    for (long int i = begin__; i < end__; ++i) \
    {

#define OMP_PARALLEL_FOR_LOOP_END \
    } \
}


#define OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,T__,op_loop_bcells__,...) \
{\
    const int op_loop_lower_X__ = std::min((long int)(-(op_loop_bcells__)),(long int)0); \
    const int op_loop_lower_Y__ = std::min((long int)(-(op_loop_bcells__)),(long int)0); \
    const int op_loop_lower_Z__ = std::min((long int)(-(op_loop_bcells__)),(long int)0); \
    if ( op_loop_bcells__ > (long int)((T__).Bcells() )){std::cerr << "BOUNDARY TOO SMALL! PLEASE ADJUST! BCELLS NEEDED " << op_loop_bcells__ << std::endl; throw std::invalid_argument("Bcells");}\
    const int op_loop_upper_X__ = std::max((long int)((T__).sizeX() + (op_loop_bcells__)),(long int)((T__).sizeX())); \
    const int op_loop_upper_Y__ = std::max((long int)((T__).sizeY() + (op_loop_bcells__)),(long int)((T__).sizeY())); \
    const int op_loop_upper_Z__ = std::max((long int)((T__).sizeZ() + (op_loop_bcells__)),(long int)((T__).sizeZ())); \
    _Pragma(STRINGIFY(omp parallel for collapse(OMP_COLLAPSE_LOOPS) schedule(OMP_SCHEDULING_TYPE,OMP_DYNAMIC_CHUNKSIZE) __VA_ARGS__) ) \
    for (long int i = op_loop_lower_X__; i < op_loop_upper_X__; ++i) \
    for (long int j = op_loop_lower_Y__; j < op_loop_upper_Y__; ++j) \
    for (long int k = op_loop_lower_Z__; k < op_loop_upper_Z__; ++k) \
    {


#define OMP_PARALLEL_STORAGE_LOOP_END \
    } \
}

#define OMP_PARALLEL_REDUCTION_STORAGE_LOOP_BEGIN(i,j,k,T__,op_loop_bcells__,...) \
{\
    const int op_loop_lower_X__ = std::min((long int)(-(op_loop_bcells__)),(long int)0); \
    const int op_loop_lower_Y__ = std::min((long int)(-(op_loop_bcells__)),(long int)0); \
    const int op_loop_lower_Z__ = std::min((long int)(-(op_loop_bcells__)),(long int)0); \
    if ( op_loop_bcells__ > (long int)((T__).Bcells() )){std::cerr << "BOUNDARY TOO SMALL! PLEASE ADJUST! BCELLS NEEDED " << op_loop_bcells__ << std::endl; throw std::invalid_argument("Bcells");}\
    const int op_loop_upper_X__ = std::max((long int)((T__).sizeX() + (op_loop_bcells__)),(long int)((T__).sizeX())); \
    const int op_loop_upper_Y__ = std::max((long int)((T__).sizeY() + (op_loop_bcells__)),(long int)((T__).sizeY())); \
    const int op_loop_upper_Z__ = std::max((long int)((T__).sizeZ() + (op_loop_bcells__)),(long int)((T__).sizeZ())); \
    _Pragma(STRINGIFY(omp for collapse(OMP_COLLAPSE_LOOPS) schedule(OMP_SCHEDULING_TYPE,OMP_DYNAMIC_CHUNKSIZE) __VA_ARGS__) ) \
    for (long int i = op_loop_lower_X__; i < op_loop_upper_X__; ++i) \
    for (long int j = op_loop_lower_Y__; j < op_loop_upper_Y__; ++j) \
    for (long int k = op_loop_lower_Z__; k < op_loop_upper_Z__; ++k) \
    {


#define OMP_PARALLEL_REDUCTION_STORAGE_LOOP_END \
    } \
}

#ifdef DEBUG
    #define exit(i); raise(SIGABRT);
#endif

#ifndef PARALLEL
    #ifdef _OPENMP
        #define myclock_t double
        #define mygettime() omp_get_wtime()
        #define OP_CLOCKS_PER_SEC 1.0
        #define OMP_COLLAPSE_LOOPS 2
        #define OMP_DYNAMIC_CHUNKSIZE 8
        #define OMP_CHUNKSIZE 8
        #define OMP_SCHEDULING_TYPE dynamic
    #else
        #define myclock_t clock_t
        #define mygettime() clock()
        #define OP_CLOCKS_PER_SEC CLOCKS_PER_SEC
    #endif

    #define STORAGE_LOOP_BEGIN_NEW(i,j,k,T__,op_loop_type__,op_stencil_size__) \
    {\
        const int op_loop_lower_X__ = std::min((long int)(((op_loop_type__==0)?0:-(long int)((T__).Bcells()))+op_stencil_size__),(long int)0); \
        const int op_loop_lower_Y__ = std::min((long int)(((op_loop_type__==0)?0:-(long int)((T__).Bcells()))+op_stencil_size__),(long int)0); \
        const int op_loop_lower_Z__ = std::min((long int)(((op_loop_type__==0)?0:-(long int)((T__).Bcells()))+op_stencil_size__),(long int)0); \
        const int op_loop_upper_X__ = std::max((long int)((T__).sizeX() - (((op_loop_type__==0)?0:-(long int)((T__).Bcells()))+op_stencil_size__)),(long int)((T__).sizeX())); \
        const int op_loop_upper_Y__ = std::max((long int)((T__).sizeY() - (((op_loop_type__==0)?0:-(long int)((T__).Bcells()))+op_stencil_size__)),(long int)((T__).sizeY())); \
        const int op_loop_upper_Z__ = std::max((long int)((T__).sizeZ() - (((op_loop_type__==0)?0:-(long int)((T__).Bcells()))+op_stencil_size__)),(long int)((T__).sizeZ())); \
        if (op_loop_type__) { \
            if ( op_stencil_size__ > (long int)((T__).Bcells() )){std::cerr << "BOUNDARY TOO SMALL! PLEASE ADJUST! BCELLS NEEDED " << op_stencil_size__ << std::endl; throw std::invalid_argument("Bcells");}\
        } \
        for (long int i = op_loop_lower_X__; i < op_loop_upper_X__; ++i) \
        for (long int j = op_loop_lower_Y__; j < op_loop_upper_Y__; ++j) \
        for (long int k = op_loop_lower_Z__; k < op_loop_upper_Z__; ++k) \
        {

    #define OMP_PARALLEL_STORAGE_LOOP_BEGIN_NEW(i,j,k,T__,op_loop_type__,op_stencil_size__,...) \
    {\
        const int op_loop_lower_X__ = std::min((long int)(((op_loop_type__==0)?0:-(long int)((T__).Bcells()))+op_stencil_size__),(long int)0); \
        const int op_loop_lower_Y__ = std::min((long int)(((op_loop_type__==0)?0:-(long int)((T__).Bcells()))+op_stencil_size__),(long int)0); \
        const int op_loop_lower_Z__ = std::min((long int)(((op_loop_type__==0)?0:-(long int)((T__).Bcells()))+op_stencil_size__),(long int)0); \
        const int op_loop_upper_X__ = std::max((long int)((T__).sizeX() - (((op_loop_type__==0)?0:-(long int)((T__).Bcells()))+op_stencil_size__)),(long int)((T__).sizeX())); \
        const int op_loop_upper_Y__ = std::max((long int)((T__).sizeY() - (((op_loop_type__==0)?0:-(long int)((T__).Bcells()))+op_stencil_size__)),(long int)((T__).sizeY())); \
        const int op_loop_upper_Z__ = std::max((long int)((T__).sizeZ() - (((op_loop_type__==0)?0:-(long int)((T__).Bcells()))+op_stencil_size__)),(long int)((T__).sizeZ())); \
        if (op_loop_type__) { \
            if ( op_stencil_size__ > (long int)((T__).Bcells() )){std::cerr << "BOUNDARY TOO SMALL! PLEASE ADJUST! BCELLS NEEDED " << op_stencil_size__ << std::endl; throw std::invalid_argument("Bcells");} \
        } \
        _Pragma(STRINGIFY(omp parallel for collapse(OMP_COLLAPSE_LOOPS) schedule(OMP_SCHEDULING_TYPE,OMP_DYNAMIC_CHUNKSIZE) __VA_ARGS__) ) \
        for (long int i = op_loop_lower_X__; i < op_loop_upper_X__; ++i) \
        for (long int j = op_loop_lower_Y__; j < op_loop_upper_Y__; ++j) \
        for (long int k = op_loop_lower_Z__; k < op_loop_upper_Z__; ++k) \
        {

#else
    #include "mpi.h"
    #ifdef _OPENMP
        #define myclock_t double
        #define mygettime() omp_get_wtime()
        #define OP_CLOCKS_PER_SEC 1.
        #define OMP_COLLAPSE_LOOPS 2
        #define OMP_DYNAMIC_CHUNKSIZE 8
        #define OMP_CHUNKSIZE 8
        #define OMP_SCHEDULING_TYPE dynamic
    #else
        #define myclock_t clock_t
        #define mygettime() clock()
        #define OP_CLOCKS_PER_SEC CLOCKS_PER_SEC
    #endif

    extern int op_loop_stencil_count;

    #define STORAGE_LOOP_BEGIN_NEW(i,j,k,T__,op_loop_type__,op_stencil_size__) \
    {\
        const int op_loop_lower_X__ = std::min((long int)(((op_loop_type__==0)?0:-(long int)((T__).Bcells()))+op_stencil_size__+op_loop_stencil_count),(long int)0); \
        const int op_loop_lower_Y__ = std::min((long int)(((op_loop_type__==0)?0:-(long int)((T__).Bcells()))+op_stencil_size__+op_loop_stencil_count),(long int)0); \
        const int op_loop_lower_Z__ = std::min((long int)(((op_loop_type__==0)?0:-(long int)((T__).Bcells()))+op_stencil_size__+op_loop_stencil_count),(long int)0); \
        const int op_loop_upper_X__ = std::max((long int)((T__).sizeX() - (((op_loop_type__==0)?0:-(long int)((T__).Bcells()))+op_stencil_size__+op_loop_stencil_count)),(long int)((T__).sizeX())); \
        const int op_loop_upper_Y__ = std::max((long int)((T__).sizeY() - (((op_loop_type__==0)?0:-(long int)((T__).Bcells()))+op_stencil_size__+op_loop_stencil_count)),(long int)((T__).sizeY())); \
        const int op_loop_upper_Z__ = std::max((long int)((T__).sizeZ() - (((op_loop_type__==0)?0:-(long int)((T__).Bcells()))+op_stencil_size__+op_loop_stencil_count)),(long int)((T__).sizeZ())); \
        if (op_loop_type__) { \
            if ( op_stencil_size__+op_loop_stencil_count > (long int)((T__).Bcells() )){std::cerr << "BOUNDARY TOO SMALL! PLEASE ADJUST! BCELLS NEEDED " << op_stencil_size__ + op_loop_stencil_count<< std::endl; throw std::invalid_argument("Bcells");}\
        } \
        for (long int i = op_loop_lower_X__; i < op_loop_upper_X__; ++i) \
        for (long int j = op_loop_lower_Y__; j < op_loop_upper_Y__; ++j) \
        for (long int k = op_loop_lower_Z__; k < op_loop_upper_Z__; ++k) \
        {

    #define OMP_PARALLEL_STORAGE_LOOP_BEGIN_NEW(i,j,k,T__,op_loop_type__,op_stencil_size__,...) \
    {\
        const int op_loop_lower_X__ = std::min((long int)(((op_loop_type__==0)?0:-(long int)((T__).Bcells()))+op_stencil_size__+op_loop_stencil_count),(long int)0); \
        const int op_loop_lower_Y__ = std::min((long int)(((op_loop_type__==0)?0:-(long int)((T__).Bcells()))+op_stencil_size__+op_loop_stencil_count),(long int)0); \
        const int op_loop_lower_Z__ = std::min((long int)(((op_loop_type__==0)?0:-(long int)((T__).Bcells()))+op_stencil_size__+op_loop_stencil_count),(long int)0); \
        const int op_loop_upper_X__ = std::max((long int)((T__).sizeX() - (((op_loop_type__==0)?0:-(long int)((T__).Bcells()))+op_stencil_size__+op_loop_stencil_count)),(long int)((T__).sizeX())); \
        const int op_loop_upper_Y__ = std::max((long int)((T__).sizeY() - (((op_loop_type__==0)?0:-(long int)((T__).Bcells()))+op_stencil_size__+op_loop_stencil_count)),(long int)((T__).sizeY())); \
        const int op_loop_upper_Z__ = std::max((long int)((T__).sizeZ() - (((op_loop_type__==0)?0:-(long int)((T__).Bcells()))+op_stencil_size__+op_loop_stencil_count)),(long int)((T__).sizeZ())); \
        if (op_loop_type__) { \
            if ( op_stencil_size__+op_loop_stencil_count > (long int)((T__).Bcells() )){std::cerr << "BOUNDARY TOO SMALL! PLEASE ADJUST! BCELLS NEEDED " << op_stencil_size__ + op_loop_stencil_count << std::endl; throw std::invalid_argument("Bcells");}\
        } \
        _Pragma(STRINGIFY(omp parallel for collapse(OMP_COLLAPSE_LOOPS) schedule(OMP_SCHEDULING_TYPE,OMP_DYNAMIC_CHUNKSIZE) __VA_ARGS__) ) \
        for (long int i = op_loop_lower_X__; i < op_loop_upper_X__; ++i) \
        for (long int j = op_loop_lower_Y__; j < op_loop_upper_Y__; ++j) \
        for (long int k = op_loop_lower_Z__; k < op_loop_upper_Z__; ++k) \
        {

    #ifndef MACROS_H
        #define MACROS_H
        extern int MPI_RANK;
        extern int MPI_SIZE;
        extern int OPP_BlockSize_X;
        extern int OPP_BlockSize_Y;
        extern int OPP_BlockSize_Z;
        extern int OPP_BlockNumber_X;
        extern int OPP_BlockNumber_Y;
        extern int OPP_BlockNumber_Z;
    #endif
#endif