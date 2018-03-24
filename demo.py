

import numpy as np
import setiptah.polyglint2d as pglint


import bintrees     # for trapezoid containment query DS



def trapezoid_tree( trap_dict ) :
    # for trapezoid containment query DS
    tree = bintrees.RBTree()
    
    for I in trap_dict :
        a, b = I
        L1, L2 = trap_dict[I]
        
        tree[a] = ( a, b, L1, L2 )
        
    return tree


def _bounding_box( vertices ) :
    X = [ x for x,y in vertices ]
    Y = [ y for x,y in vertices ]
    return [ min(X), max(X), min(Y), max(Y) ]



class RectGrid :
    
    def __init__(self, **kwargs ) :
        self.xstart = kwargs.get( 'xstart', 0. )
        self.ystart = kwargs.get( 'ystart', 0. )
        
        delta = kwargs.get( 'delta', 1. )
        xdelta = kwargs.get( 'xdelta', None )
        ydelta = kwargs.get( 'ydelta', None )
        
        if xdelta is None : xdelta = delta
        self.xdelta = xdelta
        
        if ydelta is None : ydelta = delta
        self.ydelta = ydelta
        
    def grid_cell(self, point ) :
        x, y = point
        i = float( x - self.xstart ) / self.xdelta
        i = int( np.floor(i) )
        j = float( y - self.ystart ) / self.ydelta
        j = int( np.floor(j) )
        return i, j
    
    def cell_bottomleft( self, cell ) :
        i, j = cell
        x = self.xstart + i * self.xdelta
        y = self.ystart + j * self.ydelta
        return x,y

    def cell_topright( self, cell ) :
        i, j = cell
        x = self.xstart + (i+1) * self.xdelta
        y = self.ystart + (j+1) * self.ydelta
        return x,y




def sample_inequality_interior( A, b, bbox=None, **kwargs ) :
    
    # get a bounding box
    if bbox is None :
        vertices = pglint.enumerate_vertices_2d( A, b )
        bbox = _bounding_box( vertices )
        
    xmin, xmax, ymin, ymax = bbox

    # round it to the grid
    grid = RectGrid(**kwargs)
    xmin, ymin = grid.cell_bottomleft( grid.grid_cell( (xmin,ymin) ) )
    xmax, ymax = grid.cell_topright( grid.grid_cell( (xmax,ymax) ) )
    
    xspace = np.arange(xmin,xmax, grid.xdelta )
    yspace = np.arange(ymin,ymax, grid.ydelta )
    
    res = []
    
    for x in xspace :
        for y in yspace :
            xx = x,y
            if np.all( np.dot( A, xx ) < b ) :      # strict?
                res.append( xx )
    
    return res


def sample_hull_interior( vertices, **kwargs ) :
    A, b = pglint.vertices_to_hull_inequality( vertices )
    bbox = _bounding_box( vertices )
    return sample_inequality_interior( A, b, bbox=bbox, **kwargs )


def sample_connecting_line( p, q, delta ) :
    p = np.array(p)
    q = np.array(q)
    width = np.linalg.norm( q - p )
    
    n = np.ceil( width / delta )
    n = max( n, 2 )
    tspace = np.linspace(0,1,n)
    
    res = [ p + t * ( q - p ) for t in tspace ]
    return res
    
    
class two_point_interpolator :
    def __init__(self, p, q ) :
        self.p = np.array(p)
        self.q = np.array(q)
        
        v = self.q - self.p
        self.length = np.linalg.norm( v )
        assert self.length > 0
        
        self._norm = v / self.length
        
    def __call__(self, x ) :
        return self.p + x * self._norm
    
def arange_minprec(a,b,eps) :
    n = int( np.ceil( float(b-a) / eps ) )
    n = max(n,2)
    return np.linspace(a,b,n)
    
    
    
        
        
    
def hull_perimeter( vertices ) :
    
    hull = pglint.convexHull( vertices )    # clock-wise circulation, whoops
    cycle = zip( hull, hull[1:] + hull[:1] )
    
    res = []
    for p,q in cycle :
        line = two_point_interpolator(p,q)
        res.append( line )
        #line = sample_connecting_line( p, q, delta )
        #res.append( line )
        
    return res




""" some drawing convenience """

import matplotlib.lines as lines
import matplotlib.patches as patches
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D



def draw_hull( vertices, ax=None, **kwargs ) :
    if ax is None : ax = plt.gca()
    
    hull = pglint.convexHull(vertices)
    path = zip( hull[:1], hull[1:] )
    
    hull_art = patches.Polygon( path, closed=True, **kwargs )
    ax.add_patch( hull_art ) 
    
    return ax
    



def draw_traps( traps, ax=None ) :
    
    if ax is None : ax = plt.gca()
    
    for I in traps :
        a,b = I
        L1, L2 = traps[I]
        m1, c1 = L1
        m2, c2 = L2
        #Y = [ m*x + b for m,b in traps[I] for x in I ]
        """
            3---4
            |   |
            1---2
        """
        y1 = m1 * a + c1
        y2 = m1 * b + c1
        y3 = m2 * a + c2
        y4 = m2 * b + c2

        path = [ (a,y1), (b,y2), (b,y4), (a,y3) ]
        trap_art = patches.Polygon( path, closed=True, fill=False )
        
        if False :
            # left vertical
            l1 = lines.Line2D( [a,a], [y1,y3] )
            # right vertical
            l2 = lines.Line2D( [b,b], [y2,y4] )
            # bottom
            l3 = lines.Line2D( [a,b], [y1,y2] )
            # top
            l4 = lines.Line2D( [a,b], [y3,y4] )
            
        ax.add_patch( trap_art )
        
    return ax













if __name__ == '__main__' :
        
    plt.close('all')
    
    
    """ instance """
    
    if False :
        # choose N points randomly in the plane
        N = 10
        X = np.random.rand(N)
        Y = np.random.rand(N)
        P = zip(X,Y)
        #points = np.random.rand(N,2)
        #points = [ tuple( points[k,:] for k in xrange(N) ]

    else :
        P = [
             (-1,0), (-.3,-1), (-.1,1.3), (1,.9), (2,-.5)
             ] 


    def mysurf(x,y) :
        return ( np.cos(2.2 * np.pi * x) + 1. ) * ( np.sin(1.7 * np.pi * y) + 1. ) + 1



    """ DEMO code """
    
    
    """ 0. Earlier geometry demo """
    # compute convex hull
    hull = pglint.convexHull( P )
    
    if False :
        # scatter plot of the hull points
        plt.figure()
        Xh = [ x for x,y in hull ]
        Yh = [ y for x,y in hull ]
        
        #plt.scatter(X,Y)        # all points
        plt.scatter(Xh,Yh, c='r', marker='x' )  # just hull
        
        for k,x in enumerate(Xh) :
            y = Yh[k]
            plt.annotate( '%d' % k, (x,y) )

    """ 0.1 But *do* get the hull inequality """            
    # obtain the hull inequality
    A, b = pglint.vertices_to_hull_inequality( P )
    
    if False :    
        # obtain the hull indirectly from the inequalities
        hull2 = pglint.enumerate_vertices_2d( A, b )
        Xh2 = [ x for x,y in hull2 ]
        Yh2 = [ y for x,y in hull2 ]
        
        # scatter plot the indirect hull
        plt.scatter(Xh2,Yh2, c='b', marker='o' )
    
            
        
    
    
    
    """ 1. show the surface, alone """
    fig = plt.figure()
    ax = Axes3D(fig)
    
    bbox = _bounding_box( P )
    
    def enlarge_bbox( bbox, alpha ) :
        xmin,xmax,ymin,ymax = bbox
        W,H = xmax-xmin, ymax-ymin
        xmid = .5 * ( xmin + xmax )
        ymid = .5 * ( ymin + ymax )
        
        xmin = xmid - .5 * alpha * W
        xmax = xmid + .5 * alpha * W
        ymin = ymid - .5 * alpha * H
        ymax = ymid + .5 * alpha * H
        return xmin,xmax,ymin,ymax
    
    xmin,xmax,ymin,ymax = enlarge_bbox( bbox, 1.15 )
    N = 50
    xs = np.linspace(xmin,xmax,N)
    ys = np.linspace(ymin,ymax,N)
    X,Y = np.meshgrid(xs,ys)
    
    func = np.vectorize( mysurf )
    Z = func(X,Y)
    
    ax.plot_surface(X,Y,Z,alpha=.2)
    ax.set_aspect('equal')
    
    if False :
        # THEN!!!
        # trisurf the interior
        X = [ x for x,y in interior_samples ]
        Y = [ y for x,y in interior_samples ]
        Z = [ mysurf(x,y) for x,y in interior_samples ]
        
        
        import matplotlib.cm as cm
        #ax.plot_trisurf( X, Y, Z, cmap=cm.jet, alpha=.2 )
        ax.plot_trisurf( X, Y, Z, alpha=.2 )
        ax.set_aspect('equal')


    
    
    # int_{x=a}^b \int_{y=m_1*x+c_1}^{m_2*x+c_2} c_x*x + c_y*y + c_0 dy dx
    
    
    
    
    """ 2. show the region of interest """
    plt.figure()
    ax = plt.gca()

    DELTA = .1
    interior_samples = sample_hull_interior( P, delta=DELTA )
    
    if True :
        # draw hull
        cycle = zip( hull, hull[1:] + hull[:1] )
        for p,q in cycle :
            xp,yp = p
            xq,yq = q
            ax.plot( [xp,xq], [yp,yq], c='k', linewidth=1.5 )
        
    if True :
        # show interior samples
        xs = [ x for x,y in interior_samples ]
        ys = [ y for x,y in interior_samples ]
        plt.scatter( xs, ys )
    
    if False :
        # draw bounding box
        xmin, xmax, ymin, ymax = _bounding_box( P )
        xy = xmin, ymin
        W, H = xmax-xmin, ymax-ymin
        bbox_art = patches.Rectangle( xy, W, H, fill=False )
        ax.add_patch( bbox_art )
        
    if False :
        # draw trapezoids
        traps = pglint.convex2d_to_traps( P )
        draw_traps( traps, ax )

    ax.set_aspect('equal')
    
    
    
    """ 3d stuff """
    
    if False :

    
        fig = plt.figure()
        ax = Axes3D(fig)
    
    
        if False :
            # wall-up the perimeter
            perim = hull_perimeter( P )
            for line in perim :
                tspace = arange_minprec( 0., line.length, DELTA )
                seq = [ line(t) for t in tspace ]
                #X = [ x for x,y in seq ]
                #Y = [ y for x,y in seq ]
                
                # just try to draw the top...
                xs = [ x for x,y in seq ]
                ys = [ y for x,y in seq ]
                zs = [ mysurf(x,y) for x,y in seq ]
                
                options = dict( c='b', linewidth=2. )
                myplot = lambda xs,ys,zs : ax.plot( xs,ys,zs, **options )
                
                myplot(xs,ys,zs)
                #ax.plot( xs,ys,zs, **options )
                x1,y1,z1 = xs[0],ys[0],zs[0]
                x2,y2,z2 = xs[-1],ys[-1],zs[-1]
                myplot([x1,x2],[y1,y2],[0,0])
                myplot( [x1,x1], [y1,y1], [0,z1] )
                myplot( [x2,x2], [y2,y2], [0,z2] )
                
                
                
                if False :
                    path = [ (0,0) ] + zip( tspace, Z ) + [ (line.length,0) ]
                    wall = patches.Polygon( path, closed=True, alpha=.2 )
                    ax.add_patch( wall )
                
                    # orientation --- functions operate on the patch directly
                    from mpl_patch_placement import my_pathpatch_2d_to_3d
                    #pathpatch_2d_to_3d, pathpatch_translate
                    
                    p, q = line(0.), line(line.length)
                    v = np.zeros(3) ; v[:2] = q - p     # parallel vector
                    UP = np.array([0,0,1])
                    normal_vec = np.cross( v, UP )      # drawing plane normal vector
                    
                    root = np.zeros(3) ; root[:2] = p   # root at p
                    my_pathpatch_2d_to_3d( wall, normal_vec, root )    # place wall in normal plane, rooted at root
    
            
        if True :
            verts = []
                    
            # wall-up the perimeter
            perim = hull_perimeter( P )
            for line in perim :
                tspace = arange_minprec( 0., line.length, DELTA )
                seq = [ line(t) for t in tspace ]
                coords = [ ( x,y, mysurf(x,y) ) for x,y in seq ]
    
                x1,y1 = seq[0]
                x2,y2 = seq[-1]    
                path = [ (x1,y1,0) ] + coords + [ (x2,y2,0) ]
                verts.append(path)
                
            from mpl_toolkits.mplot3d.art3d import Poly3DCollection
            polys = Poly3DCollection( verts, alpha=.1 )
            #X = [ x for x,y in seq ]
            #Y = [ y for x,y in seq ]
            #Z = [ mysurf(x,y) for x,y in seq ]
            #coords = zip(X,Y,Z)
            
            ax.add_collection3d( polys )
            
                

    
    
                
        if True :
            # THEN!!!
            # trisurf the interior
            X = [ x for x,y in interior_samples ]
            Y = [ y for x,y in interior_samples ]
            Z = [ mysurf(x,y) for x,y in interior_samples ]
            
            
            import matplotlib.cm as cm
            #ax.plot_trisurf( X, Y, Z, cmap=cm.jet, alpha=.2 )
            ax.plot_trisurf( X, Y, Z, alpha=.2 )
            ax.set_aspect('equal')
    
    
    
        
        
    
    
    
    