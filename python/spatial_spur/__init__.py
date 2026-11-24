"""Spatial correlation-robust inference and spatial unit root diagnostics."""

from importlib.metadata import PackageNotFoundError, version

from ._native import TiedEigenspaceError, build_info
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

try:
    __version__ = version("spatial-spur")
except PackageNotFoundError:
    __version__ = build_info().get("package_version", "0.6.0")

__all__ = [
    "HalfLifeResult",
    "PipelineResult",
    "SCPCResult",
    "Simulation",
    "TestResult",
    "TiedEigenspaceError",
    "__version__",
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
