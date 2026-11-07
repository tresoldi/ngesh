"""Versioned samplers using only Random.random(), whose stream Python guarantees.

Gamma: Marsaglia–Tsang (2000), shape boosting for shape < 1.
Poisson: exponential arrivals in chunks of mean <= 20 (additivity).
Normal: Box–Muller without a cached partner; beta: ratio of gammas.
Discrete gamma: equal-probability bins, conditional means via incomplete gamma.
"""
from __future__ import annotations

import hashlib
import math
import random
from typing import TypeAlias

Seed: TypeAlias = str | int | bytes


def derive_seed(seed: Seed, component: str = "simulation", replicate: int = 0) -> int:
    """Type-tagged, length-prefixed SHA-256 encoding; no Python hash()."""
    if isinstance(seed, bool) or not isinstance(seed, (str, int, bytes)):
        raise TypeError("seed must be str, int or bytes (not bool)")
    if not isinstance(replicate, int) or replicate < 0:
        raise ValueError("replicate must be a nonnegative integer")
    parts = [type(seed).__name__.encode(), seed if isinstance(seed, bytes) else str(seed).encode(),
             component.encode(), str(replicate).encode()]
    payload = b"ngesh-rng-v1" + b"".join(len(x).to_bytes(8, "big") + x for x in parts)
    return int.from_bytes(hashlib.sha256(payload).digest(), "big")


def positive(value: float, name: str, zero: bool = False) -> None:
    if not math.isfinite(value) or value < 0 or (not zero and value == 0):
        raise ValueError(f"{name} must be finite and {'nonnegative' if zero else 'positive'}")


class RandomStream:
    """One private Random instance; components can derive independent streams."""

    def __init__(self, seed: Seed) -> None:
        self.seed = seed
        self.integer_seed = derive_seed(seed)
        self.random = random.Random(self.integer_seed).random

    def index(self, n: int) -> int:
        if n <= 0:
            raise ValueError("empty choice")
        return int(self.random() * n)

    def exponential(self, rate: float) -> float:
        positive(rate, "rate")
        return -math.log1p(-self.random()) / rate

    def normal(self) -> float:
        return math.sqrt(-2 * math.log1p(-self.random())) * math.cos(2 * math.pi * self.random())

    def gamma(self, shape: float, scale: float = 1.0) -> float:
        positive(shape, "shape")
        positive(scale, "scale")
        if shape < 1:
            return self.gamma(shape + 1, scale) * (1 - self.random()) ** (1 / shape)
        d = shape - 1 / 3
        c = 1 / math.sqrt(9 * d)
        while True:
            x = self.normal()
            v = 1 + c * x
            if v <= 0:
                continue
            v **= 3
            u = 1 - self.random()
            if u < 1 - 0.0331 * x ** 4 or math.log(u) < 0.5 * x * x + d * (1 - v + math.log(v)):
                return scale * d * v

    def poisson(self, mean: float) -> int:
        positive(mean, "mean", zero=True)
        count = 0
        while mean > 0:
            chunk = min(mean, 20)
            t = self.exponential(1)
            while t <= chunk:
                count += 1
                t += self.exponential(1)
            mean -= chunk
        return count

    def beta(self, a: float, b: float) -> float:
        x, y = self.gamma(a), self.gamma(b)
        if x + y == 0:
            raise ValueError("beta shapes too small for floating-point precision")
        return x / (x + y)

    def lognormal(self, sigma: float) -> float:
        """Mean-one multiplier exp(N(-sigma²/2, sigma²))."""
        positive(sigma, "sigma", zero=True)
        return math.exp(sigma * self.normal() - sigma * sigma / 2) if sigma else 1.0


def gamma_cdf(a: float, x: float) -> float:
    """Regularized lower incomplete gamma, series / continued fraction."""
    if x <= 0:
        return 0.0
    factor = math.exp(a * math.log(x) - x - math.lgamma(a))
    if x < a + 1:
        term = total = 1 / a
        for n in range(1, 10000):
            term *= x / (a + n)
            total += term
            if abs(term) <= abs(total) * 1e-14:
                return min(1.0, factor * total)
    else:
        b = x + 1 - a
        c, d = 1e300, 1 / b
        h = d
        for n in range(1, 10000):
            an = -n * (n - a)
            b += 2
            d = an * d + b
            c = b + an / c
            d = 1 / (d if abs(d) > 1e-300 else 1e-300)
            c = c if abs(c) > 1e-300 else 1e-300
            delta = d * c
            h *= delta
            if abs(delta - 1) < 1e-14:
                return max(0.0, 1 - factor * h)
    raise ArithmeticError("incomplete gamma did not converge")


def gamma_categories(alpha: float, categories: int) -> list[float]:
    """Yang-style equal-probability conditional means, exactly mean one."""
    positive(alpha, "alpha")
    if not isinstance(categories, int) or categories < 1:
        raise ValueError("categories must be a positive integer")
    boundaries = [0.0]
    for i in range(1, categories):
        lo, hi = 0.0, max(1.0, alpha)
        while gamma_cdf(alpha, hi) < i / categories:
            hi *= 2
        for _ in range(100):
            mid = (lo + hi) / 2
            if gamma_cdf(alpha, mid) < i / categories:
                lo = mid
            else:
                hi = mid
        boundaries.append((lo + hi) / 2)
    cumulative = [gamma_cdf(alpha + 1, x) for x in boundaries] + [1.0]
    return [categories * (b - a) for a, b in zip(cumulative, cumulative[1:])]
