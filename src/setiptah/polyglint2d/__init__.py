"""polyglint2d: Computing two-dimensional bounded linear integrals."""

"""

MIT License

Copyright (c) 2025 Kyle Treleaven

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""
import itertools
import numpy as np
from numbers import Real
from typing import Collection, Tuple
from . import _geom

__all__ = [
    "integrate",
    "integrate_over_convexhull",
    "lint_trap",
    "enumerate_vertices_2d",
]

DEFAULT_SENSITIVITY = 10e-12  # vertex inclusion buffer


polyval = np.polyval
polyadd = np.polyadd
polysub = np.polysub
polymul = np.polymul
polyint = np.polyint


def integrate(c: np.array, d: Real, A: np.array, b: np.array) -> Real:
    r"""Compute the integral (2D) of a linear function over a convex polygon:

    $$\int_{ Ax \leq b } c'x \, {\rm d} x,$$

    i.e., the integral of the linear function $c'x$, over the 2-dimensional polygon described by $Ax \leq b$.
    (The method assumes that the polygon is closed.)

    Args:
        c: Integrand coefficients; ${\mathbb R}^2$
        d: Integrand bias; ${\mathbb R}$
        A: Constraint matrix; ${\mathbb R}^{n \times 2}$
        b: Constraint vector; ${\mathbb R}^n$

    Returns:
        value: The value of the integral.

    Example:

    Area of the triangle $x + y \leq 1$, $x, y \geq 0$:

        >>> Ab = np.array([[1, 1, 1], [-1, 0, 0], [0, -1, 0]])
        >>> A, b = Ab[:, :-1], Ab[:, -1]
        >>> integrate([0, 0], 1, A, b.T)
        np.float64(0.5)

    """
    vertices = enumerate_vertices_2d(A, b)
    return integrate_over_convexhull(c, d, vertices)


def integrate_over_convexhull(c: np.array, d: Real, vertices: Collection) -> Real:
    r"""Compute the 2D integral of a linear function over the convex hull of given points.

    Args:
        c: Integrand coefficients; ${\mathbb R}^2$
        d: Integrand bias; ${\mathbb R}$
        vertices: Collection of vertices; [${\mathbb R}^2$]

    Returns:
        value: The value of the integral.

    Example:

    Calculating the volume contained in $0 \leq x, y \leq 1$,
    and between the $z=0$ and $z=x$ planes:

        >>> corners = [(0, 0), (0, 1), (1, 0), (1, 1)]
        >>> integrate_over_convexhull([1, 0], 0, corners)
        np.float64(0.5)

    """
    assert len(c) == 2
    cx, cy = c

    upper = _geom.upperHull(vertices, strict=True)
    lower = _geom.lowerHull(vertices, strict=True)
    traps = _geom.trapezoids2d(upper, lower)

    total = 0.0
    for (a, b), (ln_upper, ln_lower) in traps.items():
        m1, b1 = ln_lower
        m2, b2 = ln_upper
        total += lint_trap(cx, cy, d, a, b, m1, b1, m2, b2)

    return total


def lint_trap(cx, cy, d, a, b, m1, b1, m2, b2) -> Real:
    r"""Compute the integral of a linear function over a horizontal trapezoid:

    $$\int_{x=a}^b \int_{y=m_1 x + b_1}^{m_2 x + b_2} ( c_x x + c_y y + d ) \, {\rm d}y \, {\rm d}x.$$

    The closed form given by Wolfram Alpha is surprisingly ugly, so
    this method uses a semi-symbolic approach:
    It uses `numpy.{polyint, polymul}` to evaluate intermediate polynomials.

    Returns:
        value: The value of the integral.

    Examples:
        >>> lint_trap(0, 0, 1, 0, 1, 0, 0, -1, 1)
        np.float64(0.5)

    """
    # polynomial representations (in variable x) of the inner integral bounds
    fx = np.array([m1, b1])
    gx = np.array([m2, b2])
    # representation of the part of f which is a polynomial in x [alone]
    px = np.array([cx, d])

    # result of inner integral of c2*y
    term1 = 0.5 * cy * polysub(polymul(gx, gx), polymul(fx, fx))
    term2 = polymul(px, gx - fx)
    Px = polyint(polyadd(term1, term2))

    # print fx, gx, px, term1, term2, term1+term2, Px
    return polyval(Px, b) - polyval(Px, a)


def enumerate_vertices_2d(A: np.array, b: np.array, **kwargs) -> Tuple:
    r"""Enumerate the vertices of $Ax \leq b$.

    This is currently an inefficient implementation of 2-dimensional vertex enumeration.
    It simply enumerates all bases (row pairs) and checks for containment/feasibility.

    (Assumes a bounded polygon, and doesn't bother to check.)

    Args:
        A: Constraint matrix; ${\mathbb R}^{n \times 2}$
        b: Constraint vector; ${\mathbb R}^n$

    Returns:
        vertices: The vertices of the convex region bounded.

    Example:

    Enumerate the vertices of the unit square.

        >>> Ab = np.array([[-1, 0, 0], [0, -1, 0], [1, 0, 1], [0, 1, 1]])
        >>> A, b = Ab[:, :-1], Ab[:, -1]
        >>> set(enumerate_vertices_2d(A, b)) == set((x, y) for x in [0, 1] for y in [0, 1])
        True

    """
    sensitivity = kwargs.get("sensitivity", DEFAULT_SENSITIVITY)

    rows, cols = A.shape
    assert cols == 2
    assert len(b) == rows

    vertices = set()  # hopefully takes care of degenerate cases (but probably doesn't)

    E = range(rows)
    for i, j in itertools.combinations(E, 2):
        AB = A[[i, j], :]
        bB = b[[i, j]]

        if np.linalg.det(AB) == 0:
            continue  # no finite intersection
        xB = np.linalg.solve(AB, bB)
        xB = tuple(xB)

        """
        This can cause very small trapezoids, which introduce numerical issues
        (The algorithm was rejecting some vertices due to numerical issues;
        trying to fix that with a small sensitivity buffer)
        """
        if np.all(np.dot(A, xB) <= b + sensitivity):
            vertices.add(xB)

    return vertices


def vertices_to_hull_inequality(vertices):
    hull = _geom.convexHull(vertices)  # obtains a clock-wise circulation

    # print len(vertices), len(hull)
    # assert len(hull) == len(vertices)

    pairs = list(zip(hull, hull[1:] + hull[:1]))
    n = len(pairs)

    A = np.zeros((n, 2))
    b = np.zeros(n)

    for k, (p, q) in enumerate(pairs):
        xp, yp = p
        xq, yq = q

        theta = np.arctan2(yq - yp, xq - xp)
        u = [-np.sin(theta), np.cos(theta)]  # rotate 90deg to point "away"
        A[k, :] = u

        bb = np.dot(p, u)
        # print np.dot( q, u ) - bb
        # assert abs( np.dot( q, u ) - bb ) < 10e-10
        b[k] = bb

    return A, b
