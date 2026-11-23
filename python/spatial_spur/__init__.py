"""Spatial correlation-robust inference and spatial unit root diagnostics."""

from ._native import build_info
from .scpc import SCPCResult, scpc
from .spur import (
    HalfLifeResult,
    PipelineResult,
    Simulation,
    TestResult,
    load_chetty_data,
    spur,
    spurhalflife,
    spurtest,
    spurtest_i0,
    spurtest_i0resid,
    spurtest_i1,
    spurtest_i1resid,
    spurtransform,
    standardize,
)

__all__ = [
    "HalfLifeResult",
    "PipelineResult",
    "SCPCResult",
    "Simulation",
    "TestResult",
    "build_info",
    "load_chetty_data",
    "scpc",
    "spur",
    "spurhalflife",
    "spurtest",
    "spurtest_i0",
    "spurtest_i0resid",
    "spurtest_i1",
    "spurtest_i1resid",
    "spurtransform",
    "standardize",
]
