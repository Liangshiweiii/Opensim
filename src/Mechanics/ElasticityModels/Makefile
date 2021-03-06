# This file is part of OPENPHASE BASE.
# Copyright (c) 2020 OpenPhase Solutions GmbH, Bochum, Germany
# For more details visit https://www.openphase-solutions.com
# 
# OPENPHASE BASE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3.
#  
# OPENPHASE BASE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#  
# You should have received a copy of the GNU General Public License
# along with OPENPHASE BASE.  If not, see <http://www.gnu.org/licenses/>.

# Makefile for subfolder ElasticityModels

DEPTH = ../../..

# Include makefile definitions
include $(DEPTH)/Makefile.defs

# Object SUBDIRS
OBJDIR := $(DEPTH)$(OBJDIR)

# Finds cpp files in current folder and remove path
SRC := $(wildcard *.cpp)
SRC := $(notdir $(SRC))
# Exlude entries from EXCLUDE list in src/Makefile
SRC := $(filter-out $(EXCLUDE), $(SRC))

# List of object files to create
OBJS = $(SRC:%.cpp=$(OBJDIR)/%.o)
# List of header files in include and its subfolders
HEADERS = $(wildcard $(DEPTH)/include/*.h) $(wildcard $(DEPTH)/include/**/*.h)
# Dependency files
DEPS = $(patsubst %.o,%.d,$(OBJS))

INCLUDES = -I $(DEPTH)/include/
LIB = $(DEPTH)/lib/$(LIBNAME)

.PHONY : all clean

all: $(LIB)

$(LIB): $(OBJS)
	ar $(AROPTIONS) $(LIB) $(OBJS)
	ranlib $(LIB)

$(OBJDIR)/%.o: %.cpp Makefile
	$(CC) -MMD -MP $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:

# Include List of dependency files in /obj folder. Will not complain
# if not present at first compilation.
-include $(DEPS)
