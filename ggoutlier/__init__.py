"""
Python module to read a point cloud or raster file of depths, identify outliers,
write out points to a shape or laz file for further QC.
"""
from importlib.metadata import version
from .ggoutlier import main

# version is configured in pyproject.toml
__version__ = version("ggoutlier")

__all__ = [
    "main",
    "__version__",
]
