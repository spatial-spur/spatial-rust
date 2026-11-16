CARGO_TARGET_DIR ?= ../../target-stata
RUST_LIB := $(CARGO_TARGET_DIR)/release/libspatial_stata.a
PLUGIN_OUT ?= /Users/pauldavidboll/Github/spatial-spur/scpc-stata/src/scpc_rust.plugin
CC ?= clang
CFLAGS := -bundle -DSYSTEM=APPLEMAC -O3 -Wall -Wextra -Ic
LDFLAGS := -framework Accelerate

.PHONY: all rust-lib test clean

all: $(PLUGIN_OUT)

rust-lib:
	CARGO_TARGET_DIR=$(CARGO_TARGET_DIR) cargo build --release --offline --locked -p spatial-stata

$(RUST_LIB): src/lib.rs src/stata_spi.rs Cargo.toml ../spatial-core/src/*.rs
	CARGO_TARGET_DIR=$(CARGO_TARGET_DIR) cargo build --release --offline --locked -p spatial-stata

$(PLUGIN_OUT): $(RUST_LIB) c/stplugin.c c/stata_bridge.c c/stplugin.h
	$(CC) $(CFLAGS) c/stplugin.c c/stata_bridge.c $(RUST_LIB) $(LDFLAGS) -o $@

test:
	CARGO_TARGET_DIR=$(CARGO_TARGET_DIR) cargo test --offline --locked -p spatial-stata

clean:
	rm -f $(PLUGIN_OUT)
