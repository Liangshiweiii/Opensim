DEPTH = ../..

# Include makefile definitions
include $(DEPTH)/Makefile.defs

# Object SUBDIRS
OBJDIR := $(DEPTH)$(OBJDIR)
# Finds cpp files in src/ folder and remove path
SRC := $(wildcard *.cpp)
SRC := $(notdir $(SRC))
# Exlude entries from list above
SRC := $(filter-out $(EXCLUDE), $(SRC))

# List of object files to create
OBJS = $(SRC:%.cpp=$(OBJDIR)/%.o)
# List of directories in src/ 
SUBDIRS = $(sort $(dir $(wildcard */)))
# List of header files in include and its subfolders
HEADERS = $(wildcard $(DEPTH)/include/*.h) $(wildcard $(DEPTH)/include/**/*.h)
# Dependency files
DEPS = $(patsubst %.o,%.d,$(OBJS))

INCLUDES = -I $(DEPTH)/include/
LIB = $(DEPTH)/lib/$(LIBNAME)

.PHONY : all clean subdirs $(SUBDIRS)

all: subdirs $(LIB)
subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

$(LIB): $(OBJS)
	ar $(AROPTIONS) $(LIB) $(OBJS)
	ranlib $(LIB)

$(OBJDIR)/%.o: %.cpp Makefile
	$(CC) -MMD -MP $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:

# Include List of dependency files in /obj folder. Will not complain
# if not present at first compilation.
-include $(DEPS)