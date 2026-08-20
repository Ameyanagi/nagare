from nagare import ExtrapolationPolicy


def main():
    # Integer discriminants are deliberately not part of the public API.
    _ = ExtrapolationPolicy(999)
