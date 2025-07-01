"""Geometry utilities modified from various sources."""

"""convexhull.py

https://www.oreilly.com/library/view/python-cookbook/0596001673/ch17s19.html

Calculate the convex hull of a set of n 2D-points in O(n log n) time.
Taken from Berg et al., Computational Geometry, Springer-Verlag, 1997.
Prints output as EPS file.

When run from the command line it generates a random set of points
inside a square of given length and finds the convex hull for those,
printing the result as an EPS file.

Usage:

    convexhull.py <numPoints> <squareLength> <outFile>

Dinu C. Gherman

"""

WIDTH_THRESHOLD = 0.0001


def _myDet(p, q, r):
    """Calc. determinant of a special matrix with three 2D points.

    The sign, "-" or "+", determines the side, right or left,
    respectivly, on which the point r lies, when measured against
    a directed vector from p to q.
    """
    # We use Sarrus' Rule to calculate the determinant.
    # (could also use the Numeric package...)
    sum1 = q[0] * r[1] + p[0] * q[1] + r[0] * p[1]
    sum2 = q[0] * p[1] + r[0] * q[1] + p[0] * r[1]
    return sum1 - sum2


def _isRightTurn(p, q, r):
    "Do the vectors pq:qr form a right turn, or not?"
    assert p != q and q != r and p != r

    if _myDet(p, q, r) < 0:
        return True
    else:
        return False


def upperHull(P, strict=False):
    # Get a local list copy of the points and sort them lexically.
    points = sorted(P)

    # Build upper half of the hull
    upper = points[:2]
    for p in points[2:]:
        upper.append(p)
        while len(upper) > 2 and not _isRightTurn(*upper[-3:]):
            del upper[-2]

    if not strict:
        res = upper

    # extra logic to compute "strict" upper hull,
    # which throws away vertical segments on the left or right boundaries
    else:
        res = []
        most_recent_q = None
        pairs = zip(upper[:-1], upper[1:])
        for p, q in pairs:
            xp, yp = p
            xq, yq = q
            if xq == xp:
                if yp < yq:  # on the left side
                    continue  # i.e., throw away p, pick up q on the next iteration
                if yp > yq:  # on the right side, I DON'T THINK THIS CAN EVEN HAPPEN
                    break  # don't want this segment, so stop, p is covered (from last iteration) by most_recent_q
                else:
                    raise "duplicate point detected"

            # normal operation
            res.append(p)
            most_recent_q = q
        if most_recent_q is not None:
            res.append(most_recent_q)

    return res


def lowerHull(P, strict=False, clockwise=False):  # still left-to-right
    """cheating a little bit here, rotating by 180deg, using upper hull, de-rotating"""
    Q = [(-x, -y) for x, y in P]
    upper = upperHull(Q, strict=strict)
    lower = [(-x, -y) for x, y in upper]
    if not clockwise:
        lower.reverse()
    return lower


def convexHull(P):
    upper = upperHull(P)
    lower = lowerHull(P, strict=False, clockwise=True)
    return upper + lower[1:-1]


def _getLineParams(p, q):
    """

    >>> _getLineParams((-1, -2), (1, 2))
    (2.0, 0.0)

    """
    xp, yp = p
    xq, yq = q
    assert xp != xq
    m = float(yq - yp) / (xq - xp)
    b = yp - m * xp
    return m, b


def _intervals_on_xaxis(P):
    points = sorted(P)

    res = {}

    intervals = zip(points[:-1], points[1:])
    for p, q in intervals:
        xp, yp = p
        xq, yq = q
        res[(xp, xq)] = _getLineParams(p, q)

    return res


def intersectIntervals(arrgt1, arrgt2, combine=None, **kwargs):
    """
    accepts two dictionaries, and a combining function
    computes the overlay of the two arrangements, issues the combine function on the values to produce a new arrangement

    TODO: This can definitely be reduced from O(n^2).

    """
    overlay = {}
    width_threshold = kwargs.get("width_threshold", WIDTH_THRESHOLD)

    # ( not super efficient )
    for (a1, b1), val1 in arrgt1.items():
        for (a2, b2), val2 in arrgt2.items():
            a = max(a1, a2)
            b = min(b1, b2)

            # ensure the trapezoid exists and has sufficient width
            if b > a + width_threshold:
                val = combine(val1, val2)
                overlay[(a, b)] = val

    return overlay


def trapezoids2d(upper, lower):
    upper_arrgt = _intervals_on_xaxis(upper)
    lower_arrgt = _intervals_on_xaxis(lower)
    return intersectIntervals(upper_arrgt, lower_arrgt, lambda x, y: (x, y))


def convex2d_to_traps(vertices):
    upper = upperHull(vertices, strict=True)
    lower = lowerHull(vertices, strict=True)
    return trapezoids2d(upper, lower)
