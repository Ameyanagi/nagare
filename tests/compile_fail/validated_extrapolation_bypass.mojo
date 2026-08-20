from nagare import ExtrapolationPolicy


def main():
    # This was the former invariant-breaking public constructor overload.
    _ = ExtrapolationPolicy(999, _validated=True)
