#######################################################################
#  This is the makefile to create a library for the BLAS.
#  The files are grouped as follows:
#
#       SBLAS1 -- Single precision real BLAS routines
#       CBLAS1 -- Single precision complex BLAS routines
#       DBLAS1 -- Double precision real BLAS routines
#       ZBLAS1 -- Double precision complex BLAS routines
#
#       CB1AUX -- Real BLAS routines called by complex routines
#       ZB1AUX -- D.P. real BLAS routines called by d.p. complex
#                 routines
#
#      ALLBLAS -- Auxiliary routines for Level 2 and 3 BLAS
#
#       SBLAS2 -- Single precision real BLAS2 routines
#       CBLAS2 -- Single precision complex BLAS2 routines
#       DBLAS2 -- Double precision real BLAS2 routines
#       ZBLAS2 -- Double precision complex BLAS2 routines
#
#       SBLAS3 -- Single precision real BLAS3 routines
#       CBLAS3 -- Single precision complex BLAS3 routines
#       DBLAS3 -- Double precision real BLAS3 routines
#       ZBLAS3 -- Double precision complex BLAS3 routines
#
#######################################################################

#---------------------------------------------------------
#  Level 1 BLAS
#---------------------------------------------------------
set(SBLAS1 isamax0 sasum.f90 saxpy.f90 scopy.f90 sdot.f90 snrm2.f90
	srot0 srotg.f90 sscal.f90 sswap.f90 sdsdot.f90 srotmg.f90 srotm.f90)

set(CBLAS1 scabs10 scasum.f90 scnrm2.f90 icamax.f90 caxpy.f90 ccopy.f90
	cdotc0 cdotu.f90 csscal.f90 crotg.f90 cscal.f90 cswap.f90 csrot.f90)

set(DBLAS1 idamax0 dasum.f90 daxpy.f90 dcopy.f90 ddot.f90 dnrm2.f90
	drot0 drotg.f90 dscal.f90 dsdot.f90 dswap.f90 drotmg.f90 drotm.f90)

set(ZBLAS1 dcabs10 dzasum.f90 dznrm2.f90 izamax.f90 zaxpy.f90 zcopy.f90
	zdotc0 zdotu.f90 zdscal.f90 zrotg.f90 zscal.f90 zswap.f90 zdrot.f90)

set(CB1AUX isamax0 sasum.f90 saxpy.f90 scopy.f90 snrm2.f90 sscal.f90)

set(ZB1AUX idamax0 dasum.f90 daxpy.f90 dcopy.f90 dnrm2.f90 dscal.f90)

#---------------------------------------------------------------------
#  Auxiliary routines needed by both the Level 2 and Level 3 BLAS
#---------------------------------------------------------------------
set(ALLBLAS lsame0 xerbla.f90 xerbla_array.f90)

#---------------------------------------------------------
#  Level 2 BLAS
#---------------------------------------------------------
set(SBLAS2 sgemv0 sgbmv.f90 ssymv.f90 ssbmv.f90 sspmv.f90
	strmv0 stbmv.f90 stpmv.f90 strsv.f90 stbsv.f90 stpsv.f90
	sger0 ssyr.f90 sspr.f90 ssyr2.f90 sspr2.f90)

set(CBLAS2 cgemv0 cgbmv.f90 chemv.f90 chbmv.f90 chpmv.f90
	ctrmv0 ctbmv.f90 ctpmv.f90 ctrsv.f90 ctbsv.f90 ctpsv.f90
	cgerc0 cgeru.f90 cher.f90 chpr.f90 cher2.f90 chpr2.f90)

set(DBLAS2 dgemv0 dgbmv.f90 dsymv.f90 dsbmv.f90 dspmv.f90
	dtrmv0 dtbmv.f90 dtpmv.f90 dtrsv.f90 dtbsv.f90 dtpsv.f90
	dger0 dsyr.f90 dspr.f90 dsyr2.f90 dspr2.f90)

set(ZBLAS2 zgemv0 zgbmv.f90 zhemv.f90 zhbmv.f90 zhpmv.f90
	ztrmv0 ztbmv.f90 ztpmv.f90 ztrsv.f90 ztbsv.f90 ztpsv.f90
	zgerc0 zgeru.f90 zher.f90 zhpr.f90 zher2.f90 zhpr2.f90)

#---------------------------------------------------------
#  Level 3 BLAS
#---------------------------------------------------------
set(SBLAS3 sgemm0 ssymm.f90 ssyrk.f90 ssyr2k.f90 strmm.f90 strsm.f90)

set(CBLAS3 cgemm0 csymm.f90 csyrk.f90 csyr2k.f90 ctrmm.f90 ctrsm.f90
	chemm0 cherk.f90 cher2k.f90)

set(DBLAS3 dgemm0 dsymm.f90 dsyrk.f90 dsyr2k.f90 dtrmm.f90 dtrsm.f90)

set(ZBLAS3 zgemm0 zsymm.f90 zsyrk.f90 zsyr2k.f90 ztrmm.f90 ztrsm.f90
	zhemm0 zherk.f90 zher2k.f90)


set(SOURCES)
if(BUILD_SINGLE)
  list(APPEND SOURCES ${SBLAS1} ${ALLBLAS} ${SBLAS2} ${SBLAS3})
endif()
if(BUILD_DOUBLE)
  list(APPEND SOURCES ${DBLAS1} ${ALLBLAS} ${DBLAS2} ${DBLAS3})
endif()
if(BUILD_COMPLEX)
  list(APPEND SOURCES ${CBLAS1} ${CB1AUX} ${ALLBLAS} ${CBLAS2} ${CBLAS3})
endif()
if(BUILD_COMPLEX16)
  list(APPEND SOURCES ${ZBLAS1} ${ZB1AUX} ${ALLBLAS} ${ZBLAS2} ${ZBLAS3})
endif()
list(REMOVE_DUPLICATES SOURCES)

add_library(${BLASLIB} ${SOURCES})
set_target_properties(
  ${BLASLIB} PROPERTIES
  VERSION ${LAPACK_VERSION}
  SOVERSION ${LAPACK_MAJOR_VERSION}
  )
lapack_install_library(${BLASLIB})
