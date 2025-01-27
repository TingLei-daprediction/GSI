####################################################################
# FLAGS COMMON TO ALL BUILD TYPES
####################################################################

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -traceback -assume byterecl -convert big_endian -implicitnone")

####################################################################
# RELEASE FLAGS
####################################################################

#clt set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fp-model strict")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fp-model strict -fimf-precision=hight -assume protect_parens -fma -qopt-report=3")

####################################################################
# DEBUG FLAGS
####################################################################

set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -init=snan,arrays -fp-model source -debug -ftrapuv -warn all,nointerfaces -check all,noarg_temp_created -fp-stack-check -fstack-protector -fpe0")

####################################################################
# LINK FLAGS
####################################################################

set(CMAKE_Fortran_LINK_FLAGS "")

####################################################################
# FLAGS FOR AUTOPROFILING
####################################################################

set(Fortran_AUTOPROFILING_FLAGS "-finstrument-functions")

####################################################################

# Meaning of flags
# ----------------
# todo
