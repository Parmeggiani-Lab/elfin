EXE=util
CC=gcc
CFLAGS=-MMD -O3 -fmax-errors=1
DEFS=-D_TEST_JUTIL
COMPILE=$(CC) $(CFLAGS) $(DEFS)
LDFLAGS=
LDLIBS=

TEST_ARGS:=

SRC_DIR:=src
OBJ_DIR:=.obj
$(shell mkdir -p 	$(OBJ_DIR)/$(SRC_DIR))

C_SRC:=$(wildcard $(SRC_DIR)/*.c)
CC_SRC:=$(wildcard $(SRC_DIR)/*.cc)
	
OBJS:=$(C_SRC:%.c=$(OBJ_DIR)/%.o) $(CC_SRC:%.cc=$(OBJ_DIR)/%.o)
DEPS:=$(C_SRC:%.c=$(OBJ_DIR)/%.d) $(CC_SRC:%.cc=$(OBJ_DIR)/%.d)

all: $(EXE)

EXTS=c cc
define make_rule
$(OBJ_DIR)/%.o: %.$1
	$$(COMPILE) -o $$@ -c $$<
endef
$(foreach EXT,$(EXTS),$(eval $(call make_rule,$(EXT))))

$(EXE): $(OBJS)
	$(COMPILE) $(LDFLAGS) $^ $(LDLIBS) -o $@

-include $(DEPS)

fresh: clean all

test: $(EXE)
	./$(EXE) $(TEST_ARGS)

PHONY: clean

clean:
	rm -rf $(EXE) $(OBJ_DIR)/* *.dSYM *.gch