"""Precise one-dimensional interpolation primitives."""

from .cubic import CubicSplineInterpolator
from .extrapolation import ExtrapolationPolicy
from .linear import LinearInterpolator
from .search import locate_interval
from .step import StepInterpolator, StepMode
