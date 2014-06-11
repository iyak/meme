CXX = g++
CXXFLAGS += -Wall -Wextra -MMD -std=c++11 -O3
LDFLAGS += -lm
TARGET = meme
SRCS = $(wildcard *.cc)
OBJS = $(SRCS:.cc=.o)
DEPS = $(OBJS:.o=.d)
RM = rm -fv
MODE = release
ifeq ($(MODE),debug)
	CXXFLAGS += -ggdb -DNDEBUG
	CXXFLAGS := $(filter-out -O3, $(CXXFLAGS))
endif
ifeq ($(MODE),profile)
	CXXFLAGS += -pg -fno-inline-functions-called-once -fno-optimize-sibling-calls
endif

.PHONY:all information
all: information $(TARGET) 
information:
ifneq ($(MODE),release)
ifneq ($(MODE),debug)
ifneq ($(MODE),profile)
	@echo "Invalid build mode." 
	@exit 1
endif
endif
endif
	@echo "Building $(TARGET) on "$(MODE)" MODE"
	@echo "............................."

$(TARGET): $(OBJS) Makefile
	@echo linking $(OBJS)
	@$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS)

-include $(DEPS)
%.o:%.cc Makefile
	@echo compiling $<
	@$(CXX) -c $(CXXFLAGS) -o $@ $<

.SILENT:debug profile clean d prof p c
debug d:
	$(MAKE) MODE=debug   --always-make --no-print-directory
profile prof p:
	$(MAKE) MODE=profile --always-make --no-print-directory
clean c:
	$(RM) $(TARGET) $(OBJS) $(DEPS) gmon.out .*.swp .*.swo *~
