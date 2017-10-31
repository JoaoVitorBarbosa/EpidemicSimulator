#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux
CND_DLIB_EXT=so
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/EpidemicAnalysis.o \
	${OBJECTDIR}/EpidemicManager.o \
	${OBJECTDIR}/ManipulaGrafo.o \
	${OBJECTDIR}/RandomGenerator.o \
	${OBJECTDIR}/RandomWalk.o \
	${OBJECTDIR}/Simulator.o \
	${OBJECTDIR}/Utils.o \
	${OBJECTDIR}/Vertex.o \
	${OBJECTDIR}/main.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-pthread
CXXFLAGS=-pthread

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=-L/usr/lib/x86_64-linux-gnu -lboost_filesystem -lboost_iostreams -lboost_system

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/epidemic_simulator

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/epidemic_simulator: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/epidemic_simulator ${OBJECTFILES} ${LDLIBSOPTIONS}

${OBJECTDIR}/EpidemicAnalysis.o: EpidemicAnalysis.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I/usr/include -I../../../../Desenvolvimento/gnuplot-iostream -I../ThirdParties/JSON -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/EpidemicAnalysis.o EpidemicAnalysis.cpp

${OBJECTDIR}/EpidemicManager.o: EpidemicManager.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I/usr/include -I../../../../Desenvolvimento/gnuplot-iostream -I../ThirdParties/JSON -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/EpidemicManager.o EpidemicManager.cpp

${OBJECTDIR}/ManipulaGrafo.o: ManipulaGrafo.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I/usr/include -I../../../../Desenvolvimento/gnuplot-iostream -I../ThirdParties/JSON -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ManipulaGrafo.o ManipulaGrafo.cpp

${OBJECTDIR}/RandomGenerator.o: RandomGenerator.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I/usr/include -I../../../../Desenvolvimento/gnuplot-iostream -I../ThirdParties/JSON -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/RandomGenerator.o RandomGenerator.cpp

${OBJECTDIR}/RandomWalk.o: RandomWalk.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I/usr/include -I../../../../Desenvolvimento/gnuplot-iostream -I../ThirdParties/JSON -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/RandomWalk.o RandomWalk.cpp

${OBJECTDIR}/Simulator.o: Simulator.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I/usr/include -I../../../../Desenvolvimento/gnuplot-iostream -I../ThirdParties/JSON -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Simulator.o Simulator.cpp

${OBJECTDIR}/Utils.o: Utils.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I/usr/include -I../../../../Desenvolvimento/gnuplot-iostream -I../ThirdParties/JSON -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Utils.o Utils.cpp

${OBJECTDIR}/Vertex.o: Vertex.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I/usr/include -I../../../../Desenvolvimento/gnuplot-iostream -I../ThirdParties/JSON -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Vertex.o Vertex.cpp

${OBJECTDIR}/main.o: main.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I/usr/include -I../../../../Desenvolvimento/gnuplot-iostream -I../ThirdParties/JSON -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/main.o main.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
