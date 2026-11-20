CARGO_TARGET_DIR ?= ../../target-stata
RUST_LIB := $(CARGO_TARGET_DIR)/release/libspatial_stata.a
PLUGIN_OUT ?= ../../../scpc-stata/src/scpc_rust.plugin
SPUR_PLUGIN_OUT ?= ../../../spur-stata/spur_rust.plugin

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	PLATFORM_CFLAGS := -bundle -DSYSTEM=APPLEMAC
	PLATFORM_LDFLAGS := -Wl,-exported_symbol,_pginit,-exported_symbol,_stata_call
	ifeq ($(origin CC),default)
		CC := clang
	endif
else ifeq ($(UNAME_S),Linux)
	PLATFORM_CFLAGS := -shared -fPIC -DSYSTEM=OPUNIX
	PLATFORM_EXPORTS := c/exports.map
	PLATFORM_LDFLAGS := -Wl,--version-script=$(PLATFORM_EXPORTS) -ldl -lm -lpthread
else
	$(error Building the Stata plugin is supported by this Makefile on macOS and Linux)
endif

CFLAGS ?= -O3 -Wall -Wextra
ALL_CFLAGS := $(PLATFORM_CFLAGS) -Ic $(CFLAGS)
ALL_LDFLAGS := $(PLATFORM_LDFLAGS) $(LDFLAGS)

.PHONY: all rust-lib test clean

all: $(PLUGIN_OUT) $(SPUR_PLUGIN_OUT)

rust-lib:
	CARGO_TARGET_DIR=$(CARGO_TARGET_DIR) cargo build --release --offline --locked -p spatial-stata

$(RUST_LIB): src/*.rs Cargo.toml ../../Cargo.toml ../../Cargo.lock ../spatial-core/Cargo.toml ../spatial-core/src/*.rs
	CARGO_TARGET_DIR=$(CARGO_TARGET_DIR) cargo build --release --offline --locked -p spatial-stata

# Publish only the entry points in the Stata plugin ABI.
$(PLUGIN_OUT) $(SPUR_PLUGIN_OUT): $(RUST_LIB) c/stplugin.c c/stata_bridge.c c/stplugin.h $(PLATFORM_EXPORTS) Makefile
	$(CC) $(ALL_CFLAGS) c/stplugin.c c/stata_bridge.c $(RUST_LIB) $(ALL_LDFLAGS) -o $@.tmp
	mv $@.tmp $@

test:
	CARGO_TARGET_DIR=$(CARGO_TARGET_DIR) cargo test --offline --locked -p spatial-stata

clean:
	rm -f $(PLUGIN_OUT) $(SPUR_PLUGIN_OUT)
