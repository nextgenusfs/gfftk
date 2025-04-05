import importlib.metadata

try:
    __version__ = importlib.metadata.version("gfftk")
except importlib.metadata.PackageNotFoundError:
    __version__ = "24.4.1"  # Default version if package is not installed
