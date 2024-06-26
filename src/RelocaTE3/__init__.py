"""Required actions when initiating the package."""

try:
    from importlib_metadata import entry_points, metadata
except ImportError:
    from importlib.metadata import entry_points, metadata

import logging


def _map_entry_point_module(project_name):
    """Get the dictionary that map the scripts to their entry point names."""
    d = {}
    for x in entry_points().select(group="console_scripts"):
        if x.dist.name == project_name:
            d[x.module] = x.name
    return d


_meta = metadata("RelocaTE3")
__version__ = _meta["Version"]
__author__ = _meta["Author-email"]
__entry_points__ = _map_entry_point_module(_meta["Name"])

logging.basicConfig(format="%(asctime)s %(name)s %(levelname)s %(message)s", level="INFO")
logger = logging.getLogger(__name__)
