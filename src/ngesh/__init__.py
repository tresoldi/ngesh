"""Standard-library phylogenetic benchmark simulation."""
from .tree import Node, Tree
from .rng import RandomStream, derive_seed

__version__ = "2.0.0a1"
__all__ = ["Node", "Tree", "RandomStream", "derive_seed", "simulate"]

# Import after __version__ to avoid manifest import cycles.
from .simulation import simulate
