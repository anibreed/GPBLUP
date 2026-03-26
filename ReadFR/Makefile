FC = gfortran
PROJECT_ROOT = $(HOME)/GPBLUP
SRC_DIR = $(PROJECT_ROOT)/source
BIN_DIR = $(PROJECT_ROOT)/bin
LIB_DIR = $(PROJECT_ROOT)/lib

FFLAGS = -O3 -g -fcheck=all -Wall -I$(SRC_DIR)
LDFLAGS = -L$(LIB_DIR) -lgpblup_common

TARGET = ReadFR
MAIN_SRC = ReadFR.f90
MAIN_OBJ = ReadFR.o

# Common Library Path
LIB_COMMON = $(LIB_DIR)/libgpblup_common.a

all: $(BIN_DIR)/$(TARGET)

# Library build trigger
$(LIB_COMMON):
	$(MAKE) -C $(SRC_DIR)

$(MAIN_OBJ): $(MAIN_SRC) $(LIB_COMMON)
	$(FC) $(FFLAGS) -c -o $@ $<

$(BIN_DIR)/$(TARGET): $(MAIN_OBJ)
	@mkdir -p $(BIN_DIR)
	$(FC) $(FFLAGS) -o $@ $(MAIN_OBJ) $(LDFLAGS)

clean:
	rm -f *.o *.mod
	rm -f $(BIN_DIR)/$(TARGET)
	$(MAKE) -C $(SRC_DIR) clean

.PHONY: all clean
