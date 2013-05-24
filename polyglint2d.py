
import itertools
import numpy as np
import scipy as sp


class test :
    def __init__(self) :
        A = np.zeros( (3,2) )
        A[0,:] = [-1,0]
        A[1,:] = [0,-1]
        A[2,:] = [1,1]
        self.A = A
        
        self.b = np.array([ 0,0,1] )
        
        # compute area!
        self.c = np.array([0,0])
        self.d = 1
        

""" first step is to enumerate vertices of the enclosing polygon """
def enumerate_vertices_2d( A, b ) :
    """
    an inefficient implementation of 2-dimensional vertex enumeration.
    simply enumerates all bases and checks feasibility.
    returns a set of points as tuples.
    """
    # assume bounded, don't bother to check
    rows, cols = A.shape
    assert cols == 2
    assert len( b ) == rows

    vertices = set()    # hopefully takes care of degenerate cases

    E = range( rows )
    for i, j in itertools.combinations( E, 2 ) :
        AB = A[ [i,j], : ]
        bB = b[ [i,j] ]

        if np.linalg.det( AB ) == 0 : continue    # no finite intersection
        xB = np.linalg.solve( AB, bB ) ; xB = tuple(xB)
        
        if np.all( np.dot( A, xB ) <= b ) :
            vertices.add( xB )
            
    return vertices




""" next step is to compute upper- and lower- hull traversals """

"""convexhull.py

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

def _myDet( p, q, r ) :
    """Calc. determinant of a special matrix with three 2D points.
    
    The sign, "-" or "+", determines the side, right or left,
    respectivly, on which the point r lies, when measured against
    a directed vector from p to q.
    """
    # We use Sarrus' Rule to calculate the determinant.
    # (could also use the Numeric package...)
    sum1 = q[0]*r[1] + p[0]*q[1] + r[0]*p[1]
    sum2 = q[0]*p[1] + r[0]*q[1] + p[0]*r[1]
    return sum1 - sum2

def _isRightTurn( p, q, r ) :
    "Do the vectors pq:qr form a right turn, or not?"
    assert p != q and q != r and p != r
    
    if _myDet(p, q, r) < 0 :
        return True
    else :
        return False
       
def upperHull( P, strict=False ) :
    # Get a local list copy of the points and sort them lexically.
    points = map(None, P)
    points.sort()
    
    # Build upper half of the hull
    upper = points[:2]
    for p in points[2:] :
        upper.append( p )
        while len(upper) > 2 and not _isRightTurn( *upper[-3:] ) :
            del upper[-2]
            
    if not strict :
        res = upper
        
    # extra logic to compute "strict" upper hull,
    # which throws away vertical segments on the left or right boundaries
    else :
        res = []
        most_recent_q = None
        pairs = zip( upper[:-1], upper[1:] )
        for p,q in pairs :
            xp, yp = p
            xq, yq = q
            if xq == xp :
                if yp < yq :    # on the left side
                    continue        # i.e., throw away p, pick up q on the next iteration
                if yp > yq :    # on the right side, I DON'T THINK THIS CAN EVEN HAPPEN
                    break           # don't want this segment, so stop, p is covered (from last iteration) by most_recent_q
                else :
                    raise 'duplicate point detected'
                
            # normal operation
            res.append( p )
            most_recent_q = q
        if most_recent_q is not None : res.append( most_recent_q )
    
    return res


def lowerHull( P, strict=False, clockwise=False ) :        # still left-to-right
    """ cheating a little bit here, rotating by 180deg, using upper hull, de-rotating """
    Q = [ (-x,-y) for x,y in P ]
    upper = upperHull( Q, strict=strict )
    lower = [ (-x,-y) for x,y in upper ]
    if not clockwise : lower.reverse()
    return lower


""" just for fun... we don't use it """
def convexHull( P ) :
    upper = upperHull( P )
    lower = lowerHull( P, strict=False, clockwise=True )
    return upper + lower[1:-1]
    
    


""" then we can break up the x axis into trapezoidal chunks """

def _getLineParams( p, q ) :
    xp, yp = p
    xq, yq = q
    assert xp != xq
    m = float( yq - yp ) / ( xq - xp )
    b = yp - m * xp
    return m, b


def _intervals_on_xaxis( P ) :
    points = map( None, P )
    points.sort()
    
    res= {}
    
    intervals = zip( points[:-1], points[1:] )
    for p, q in intervals :
        xp, yp = p
        xq, yq = q
        res[ (xp,xq) ] = _getLineParams( p, q )
        
    return res


def intersectIntervals( arrgt1, arrgt2, combine=None ) :
    """
        accepts two dictionaries, and a combining function
        computes the overlay of the two arrangements, issues the combine function on the values to produce a new arrangement
    """
    overlay = {}
    # ( not super efficient )
    for (a1,b1), val1 in arrgt1.iteritems() :
        for (a2,b2), val2 in arrgt2.iteritems() :
            a = max( a1, a2 )
            b = min( b1, b2 )
            if a < b :
                val = combine( val1, val2 ) 
                overlay[ (a,b) ] = val
                
    return overlay


def trapezoids2d( upper, lower ) :
    upper_arrgt = _intervals_on_xaxis( upper )
    lower_arrgt = _intervals_on_xaxis( lower )
    return intersectIntervals( upper_arrgt, lower_arrgt, lambda x, y : ( x, y ) )

def convex2d_to_traps( vertices ) :
    upper = upperHull( vertices, strict=True )
    lower = lowerHull( vertices, strict=True )
    return trapezoids2d( upper, lower )
    




""" and integrate! """

def lint_trap( cx, cy, d, a, b, m1, b1, m2, b2 ) :
    """
    computes the integral :
        \int_{x=a}^b
            int_{y=m1*x+b1}^{m2*x+b2}
                f(x,y) = ( c1*x + c2*y + d )
            dy
        dx
    using scipy polyint / polymul / polyval;
    the closed form given by Wolfram Alpha is surprisingly ugly, so I've broken the integral down into stages
    """
    # polynomial representations (in variable x) of the inner integral bounds
    fx = np.array( [ m1, b1 ] )
    gx = np.array( [ m2, b2 ] )
    # representation of the part of f which is a polynomial in x [alone]
    px = np.array( [ cx, d ] )
    
    # result of inner integral of c2*y
    term1 = .5 * cy * sp.polysub( sp.polymul( gx, gx ), sp.polymul( fx, fx ) )
    term2 = sp.polymul( px, gx - fx )
    Px = sp.polyint( sp.polyadd( term1, term2 ) )
    
    #print fx, gx, px, term1, term2, term1+term2, Px
    return sp.polyval( Px, b ) - sp.polyval( Px, a )




def integrate_over_convexhull( c, d, vertices ) :
    assert len( c ) == 2
    cx, cy = c
    
    upper = upperHull( vertices, strict=True )
    lower = lowerHull( vertices, strict=True )
    traps = trapezoids2d( upper, lower )
    
    total = 0.
    for (a,b), (ln_upper,ln_lower) in traps.iteritems() :
        m1, b1 = ln_lower
        m2, b2 = ln_upper
        total += lint_trap( cx, cy, d, a, b, m1, b1, m2, b2 )
        
    return total

    # do more stuff


def integrate( c, d, A, b ) :
    """
    this routine computes the integral:
        \int_{ Ax <= b } c'x dx,
    i.e., the integral of the linear function c'x, over the 2-dimensional polygon described by Ax <= b;
    (assumes that the polygon is closed.) 
    """
    vertices = enumerate_vertices_2d( A, b )
    return integrate_over_convexhull( c, d, vertices )















