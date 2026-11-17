CARGO_TARGET_DIR ?= ../../target-stata
RUST_LIB := $(CARGO_TARGET_DIR)/release/libspatial_stata.a
PLUGIN_OUT ?= ../../../scpc-stata/src/scpc_rust.plugin

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	PLATFORM_CFLAGS := -bundle -DSYSTEM=APPLEMAC
	PLATFORM_LDFLAGS := -framework Accelerate
	ifeq ($(origin CC),default)
		CC := clang
	endif
else
	PLATFORM_CFLAGS := -shared -fPIC -DSYSTEM=OPUNIX
	PLATFORM_LDFLAGS := -ldl -lm -lpthread
endif

CFLAGS ?= -O3 -Wall -Wextra
ALL_CFLAGS := $(PLATFORM_CFLAGS) -Ic $(CFLAGS)
ALL_LDFLAGS := $(PLATFORM_LDFLAGS) $(LDFLAGS)

.PHONY: all rust-lib test clean

all: $(PLUGIN_OUT)

rust-lib:
	CARGO_TARGET_DIR=$(CARGO_TARGET_DIR) cargo build --release --offline --locked -p spatial-stata

$(RUST_LIB): src/lib.rs src/stata_spi.rs Cargo.toml ../spatial-core/Cargo.toml ../spatial-core/src/*.rs
	CARGO_TARGET_DIR=$(CARGO_TARGET_DIR) cargo build --release --offline --locked -p spatial-stata

$(PLUGIN_OUT): $(RUST_LIB) c/stplugin.c c/stata_bridge.c c/stplugin.h
	$(CC) $(ALL_CFLAGS) c/stplugin.c c/stata_bridge.c $(RUST_LIB) $(ALL_LDFLAGS) -o $@

test:
	CARGO_TARGET_DIR=$(CARGO_TARGET_DIR) cargo test --offline --locked -p spatial-stata

clean:
	rm -f $(PLUGIN_OUT)
