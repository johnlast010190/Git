#!/usr/bin/make

CXX = /usr/bin/g++

CXX_FLAGS = -m64 -ftemplate-depth=1024 -Wall -Wextra -fuse-ld=bfd -Wnon-virtual-dtor -Wno-unused-parameter -Wno-overloaded-virtual -Wno-old-style-cast -O3 -DNDEBUG
CXX_FLAGS += -Dlinux -DHELYX_ARCH_OPTION=64 -DHELYX_DP -DHELYX_LABEL_SIZE=32 -DNoRepository -DHELYX_API=40400 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -DRUNTIME_DEBUG_FLAGS
# These are treated specially in CMake
CXX_FLAGS += -fPIC -std=c++17

INCLUDES = -I$(HELYX_PROJECT_DIR)/src/OpenFOAM -I$(HELYX_PROJECT_DIR)/src/OSspecific/POSIX -I$(MPI_ARCH_PATH)/include

# Options from the dynamicCode dictionary (evaluated at run-time)
DYNAMIC_CODE_OPTIONS = |DYNAMIC_CODE_OPTIONS|

LFLAGS = -shared -Xlinker --add-needed -Xlinker --no-as-needed  -L$(HELYX_PROJECT_DIR)/platforms/$(HELYX_OPTIONS)/lib -L$(HELYX_PROJECT_DIR)/platforms/$(HELYX_OPTIONS)/lib/$(HELYX_MPI_NAME)

# Define some default libs
LIBS = -lOpenFOAM -lm -ldl |LINK_LIBS|

# Evaluated source files at run-time
SRCS = |SOURCE_FILES|

OBJS = $(SRCS:.C=.o)

# define the output file
MAIN = |LIBRARY_NAME|


.PHONY: depend clean

all: $(MAIN)
	@echo  |LIBRARY_NAME| compiled

$(MAIN): $(OBJS)
	$(CXX) $(CXX_FLAGS) $(INCLUDES) $(DYNAMIC_CODE_OPTIONS) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

.C.o:
	$(CXX) $(CXX_FLAGS) $(INCLUDES) $(DYNAMIC_CODE_OPTIONS) -c $<  -o $@

clean:
	$(RM) *.o *~ $(MAIN)

depend: $(SRCS)
	makedepend $(INCLUDES) $^
