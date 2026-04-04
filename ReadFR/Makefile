TARGET := ReadFR
ROOT_DIR := ..
BUILD_DIR := $(ROOT_DIR)/build-readfr
BIN_DIR := $(ROOT_DIR)/bin

all: $(BIN_DIR)/$(TARGET)

$(BUILD_DIR)/Makefile:
	cmake -S $(ROOT_DIR) -B $(BUILD_DIR)

$(BIN_DIR)/$(TARGET): $(BUILD_DIR)/Makefile
	$(MAKE) -C $(BUILD_DIR) $(TARGET)

clean:
	@rm -f $(BIN_DIR)/$(TARGET)
	@rm -rf $(BUILD_DIR)/CMakeFiles/$(TARGET).dir

clean-all:
	$(MAKE) -C $(BUILD_DIR) clean || true

rebuild:
	$(MAKE) clean
	$(MAKE) all

.PHONY: all clean clean-all rebuild
