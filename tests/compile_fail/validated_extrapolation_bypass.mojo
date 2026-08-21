from nagare import ExtrapolationPolicy
from std.collections import Optional


def main():
    # The former Optional[Bool] representation is no longer a public API.
    _ = ExtrapolationPolicy(Optional[Bool](None))
