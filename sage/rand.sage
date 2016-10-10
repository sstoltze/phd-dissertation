
# The variables used in the computations. These are named after the coordinates they represent.
variables = 'g2,h1,e,g,h,a,d,i,f,y,C,D,w,z' 
# w is y_2
# g2 and h1 are special versions of g and h for X_011 and X_101 respectively. They are first because this saves some issues calculating signs when we use the boundary map.

A = FreeAlgebra(QQ,14,variables)

(g2,h1,e,g,h,a,d,i,f,y,C,D,w,z) = A.gens()

# Setup the skew-symmetric part of the algebra
dim_one = [g2,h1,g,h,a,d,i,f,y,C,D] 
relations = {}
for r in dim_one:
    for s in dim_one:
        if r < s: relations[s*r] = -r*s

# The basic ring we work in.
R = A.g_algebra(relations)
(g2,h1,e,g,h,a,d,i,f,y,C,D,w,z) = R.gens()
var = R.gens()
# Elements of dimension one
dim_one = [g2,h1,g,h,a,d,i,f,y,C,D] 

# A dictionary describing the first inclusion map on generators
# From the first column to the second:
first_inclusion = { a : (a, a-g2),
                    d : (d-h1,f),
                    i : (i,i),
                    h : (h1,h),
                    g : (g,g2),
                    D : (d+y-h1,d),
                    C : (f+a,a-g2+y),
                    z : ((a-d+f+h1)*y + w, w+(a-d+f-g2)*y) }

# From the second column to the third:
second_inclusion_one = { a : a-g,
                         d : d,
                         i : i,
                         g : g, 
                         f : y,
                         y : y,
                         w : w,
                         e : e,
                         h1 : g}
second_inclusion_two = { a : a,
                         d : d-g+y,
                         i : i,
                         h : g,
                         f : d-g,
                         y : y,
                         w : w-d*y+g*y,
                         e : e,
                         g2 : g}


def build_monomial(ex,boundary_one=False,boundary_two=False):
    """Given an exponent, builds a monomial. 
    The options are for removing terms using the boundary map."""
    res = R.one()
    n = 0
    # For each power
    for i in range(len(ex)):
        if ex[i] >= 2: # All powers greater than 1 are zero
            return R.zero()
        elif boundary_one and (var[i]==h1 or var[i]==g2):
            res = res*e^ex[i]
            n += ex[i]
        else:
            # Multiply the variable to the right power
            res = res * var[i]^ex[i] 
    # Special cases
    if ex[var.index(y)] + ex[var.index(w)] >= 2:
        return R.zero() # y*w = 0
    if boundary_one and n == 0:
        return R.zero() # Remove terms that do not contain e
    if boundary_two and (ex[var.index(e)] < 1 or ex[var.index(g)]<1):
        return R.zero() # Remove terms that do not contain e and g
    return res

def clean_element(x,boundary_one=False,boundary_two=False):
    """Removes multiple powers of the same element, since these are 0."""
    r = R.zero()
    for ex in x.exponents():
        # Build the monomial from the exponent and multiply by the correct coefficient
        r += x.dict()[ex]*build_monomial(ex,
                                         boundary_one,
                                         boundary_two)
    return r

def first_boundary(x,opt=True):
    """The first boundary map, H^*(U) -> H^*(V,U)"""
    x = clean_element(x)
    temp = (R.zero(),R.zero())
    # Loop through (exponents of) terms
    for k in x.dict().keys():
        (te,mp) = (R.one(),R.one())
        # Loop through factors in terms
        for t in range(len(k)):
            if k[t] >= 1:
                # Use the inclusion on each factor and multiply out
                te = te * first_inclusion[var[t]][0]^k[t]
                mp = mp * first_inclusion[var[t]][1]^k[t]
        v = x.dict()[k]
        temp = (temp[0] + v*te,temp[1]+v*mp)
    # Do cleanup and return
    return (clean_element(temp[0],boundary_one=opt),
            clean_element(temp[1],boundary_one=opt))

def second_boundary(x,opt=True):
    """The second boundary map, H^*(V,U) -> H^*(X,V)"""
    (r,s) = clean_element(x[0]),clean_element(x[1])
    res = R.zero()
    # Loop through terms and factors as above.
    for k in r.dict().keys():
        temp = R.one()
        for i in range(len(k)):
            if k[i] >= 1:
                temp = temp * second_inclusion_one[var[i]]^k[i]
        temp = r.dict()[k] * temp
        res += temp
    for k in s.dict().keys():
        temp = R.one() # +-?
        for i in range(len(k)):
            if k[i] >= 1:
                temp = temp * second_inclusion_two[var[i]]^k[i]
        temp = temp * s.dict()[k]    
        res += temp
    # Cleanup
    return clean_element(res,boundary_two=opt)
    
def boundary_one(x):
    return first_boundary(x)

def boundary_two(x):
    return second_boundary(x)

def combined(x):
    """The composition of the two boundary maps"""
    return second_boundary(first_boundary(x))

def check_boundary(n=8):
    """Takes all products of exactly n elements and runs the combined boundary map on them, reporting when it gets a non-zero result"""
    if n == 0:
        r = combined(R.one())
        if r != R.zero():
            print "Element 1 gives ", r
    tvar = [a,d,h,i,g,C,D,z]
    p = CartesianProduct(*[tvar for k in range(n)])
    dyn = {}
    for x in p:
        start = clean_element(mul(x))
        if not dyn.get(start,False):
            r = combined(start)
            dyn[start] = True
            if r != R.zero():
                print "Elements ", x, " gives ", r

def degree(x):
    """Returns the degree of x as an element of R"""
    d = 0
    if type(x) == type(0): return 0
    temp = clean_element(x).exponents()
    if not temp: return 0
    exp = temp[0]
    for i in range(len(exp)):
        if exp[i]>0:
            if var[i] in dim_one:
                d += 1
            else:
                d += 2
    return d

def bit_to_elem(b,l):
    """Builds an element from a bitstring, by multiplying together elements of l corresponding to the bits that are 1"""
    res = R.one()
    for i in range(len(b)):
        if b[i]: res *= l[i]
    if res == R.one(): res = R.zero()
    return clean_element(res)
 
# Build the cohomology groups and split them by degree
H1_elems = [g,h,a,d,i,C,D,z]
H1 = {0 : [R.one()],
      1 : [],
      2 : [],
      3 : [],
      4 : [],
      5 : [],
      6 : [],
      7 : [],
      8 : [],
      9 : []}
for s in range(1,256):
    r = clean_element(bit_to_elem(Integer(s).bits(),
                                  H1_elems))
    H1[degree(r)].append(r) # H1[d] is all elements of degree d in H^*(U)

# First coord is X_{2}, second is X_{1}
H2_elems_one = [g,a,d,i,f,y,w]
H2_elems_two = [h,a,d,i,f,y,w]
H2 = {1 : [],
      2 : [(e,R.zero()),(R.zero(),e)],
      3 : [],
      4 : [],
      5 : [],
      6 : [],
      7 : [],
      8 : [],
      9 : []}
for s in range(1,2^len(H2_elems_one)):
    b = Integer(s).bits()
    r = (clean_element(e*bit_to_elem(b,H2_elems_one)),
         R.zero())
    if degree(r[0]):
        H2[degree(r[0])].append(r)
    r = (R.zero(),
         clean_element(e*bit_to_elem(b,H2_elems_two)))
    if degree(r[1]):
        H2[degree(r[1])].append(r)

H3_elems = [a,d,i,y,w]
H3 = {1 : [],
      2 : [],
      3 : [],
      4 : [e*g],
      5 : [],
      6 : [],
      7 : [],
      8 : [],
      9 : []}
for s in range(1,2^len(H3_elems)):
    r = clean_element(e*g*bit_to_elem(Integer(s).bits(),
                                      H3_elems))
    if degree(r):
        H3[degree(r)+1].append(r)


def find_first_kernel(get_matrix=False):
    res_m = {}
    res_ker = {0 : [],
           1 : [],
           2 : [],
           3 : [],
           4 : [],
           5 : [],
           6 : [],
           7 : [],
           8 : [],
           9 : []}
    res_im = {0 : [],
           1 : [],
           2 : [],
           3 : [],
           4 : [],
           5 : [],
           6 : [],
           7 : [],
           8 : [],
           9 : []}
    F = build_free_module(H1)
    G = build_free_module(H2)
    for deg in H1.keys():
        if G.get(deg+1) == None:
            res_ker[deg].extend(F[deg].gens())
            continue
        m = matrix([convert_element_H2(
            G[deg+1],
            first_boundary(x)).to_vector()
                    for x in H1[deg]])
        if get_matrix:
           res_m[deg] = m
           continue
        for v in m.kernel().gens():
            res_ker[deg].append(F[deg].from_vector(v))
        for v in m.image().gens():
            res_im[deg+1].append(G[deg+1].from_vector(v))
    if get_matrix: return res_m
    return (res_ker,res_im)

def build_free_module(H):
    """Builds a dictionary of free modules from the graded group H"""
    res = {}
    for deg in H.keys():
        res[deg] = CombinatorialFreeModule(ZZ,H[deg])
    return res

def find_second_kernel(get_matrix=False):
    """Finds the kernel of the second boundary map H^*(V,U) -> H^*(X,V)"""
    res_m = {}
    res_ker = {0 : [],
           1 : [],
           2 : [],
           3 : [],
           4 : [],
           5 : [],
           6 : [],
           7 : [],
           8 : [],
           9 : []}
    res_im = {0 : [],
           1 : [],
           2 : [],
           3 : [],
           4 : [],
           5 : [],
           6 : [],
           7 : [],
           8 : [],
           9 : []}
    F = build_free_module(H2)
    G = build_free_module(H3)
    # For each degree
    for deg in H2.keys():
        if G.get(deg+1) == None:
            # If our codomain is None, it is zero. Hence everything is in the kernel
            res_ker[deg].extend(F[deg].gens())
            continue
        # Build a matrix representing the boundary map
        m = matrix([convert_element_H3(
            G[deg+1],
            second_boundary(x)).to_vector()
                    for x in H2[deg]])
        if get_matrix:
           res_m[deg] = m
           continue
        # List the generators of the kernel and image
        for v in m.kernel().gens():
            res_ker[deg].append(F[deg].from_vector(v))
        for v in m.image().gens():
            res_im[deg+1].append(G[deg+1].from_vector(v))
    if get_matrix: return res_m
    return (res_ker,res_im)

def convert_element_H3(F,x):
    """Helper function. Does the same as convert_element_H2"""
    res = F.zero()
    if type(x) == type(0): return res
    for e in x.exponents():
        res += F.term(clean_element(build_monomial(e)),
                      x.dict()[e])
    return res

def convert_element_H2(F,x):
    """Turns an element of R into an element of the free group F"""
    res = F.zero()
    if not type(x[0]) == type(0):
        for e in x[0].exponents():
            res += F.term((clean_element(build_monomial(e)),
                           R.zero()),
                          x[0].dict()[e])
    if not type(x[1]) == type(0):
        for e in x[1].exponents():
            res += F.term((R.zero(),
                           clean_element(build_monomial(e))),
                          x[1].dict()[e])
    return res

def setup_kernel(get_matrix=False):
    if get_matrix:
        return (find_first_kernel(True),
                find_second_kernel(True))
    (K1,I2) = find_first_kernel()
    (K2,I3) = find_second_kernel()
    return (K1,K2,I2,I3)

def calculate_quotients():
    """Calculates the quotients on the E2 page"""
    (K1,K2,I2,I3) = setup_kernel()
    (M1,M2) = setup_kernel(get_matrix=True)
    res_one = {}
    res_two = {}
    for deg in M1.keys():
        codomain = K2.get(deg+1)
        # We only look at the subspace given by the kernel of the second boundary, so we use the basis previously calculated
        V = span(ZZ,map(lambda x:x.to_vector(),codomain))
        res_one[deg+1] = V/(M1[deg].image())
    for deg in M2.keys():
        codomain = H3.get(deg+1)
        # Use ZZ^len since the codomain is all of H3
        res_two[deg+1] = ZZ^len(codomain)/(M2[deg].image())
    return (res_one,res_two)

def print_e1():
    print "E1:\t\tH^i(U)\tH^(i+1)(V,U)\tH^(i+2)(X,V)"
    for i in H1.keys()[::-1]:
        print "Row {0}:\t\t {1:>3}\t\t {2:>3}\t\t {3:>3}".format(
            i,
            len(H1[i]),
            len(H2.get(i+1,[])),
            len(H3.get(i+2,[])))

def pretty_print():
    """Prints the E2 page."""
    (K1,K2,I2,I3) = setup_kernel()
    (Q1,Q2) = calculate_quotients()
    print "E2:\t\t H^i(U)\t H^(i+1)(V,U)\t H^(i+2)(X,V)\n\t\t ker\tker / im\t ker / im"
    for i in K1.keys()[::-1]:
        print "Row {0}:\t\t {1:>3}\t {2:>3} / {3:>2}\t {4:>3} / {5:>2}".format(
            i,
            len(K1[i]),
            len(K2.get(i+1,[])),
            len(I2.get(i+1,[])),
            len(H3.get(i+2,[])),
            len(I3.get(i+2,[])))
        
    print "\nQuotients:"
    for i in Q1.keys()[::-1]:
        print "Row {0}:\t\t\t {1:^10}\t {2:^10}".format(
            i-1,
            Q1[i].invariants(),
            Q2.get(i+1,ZZ^1/ZZ^1).invariants())
    print "0 means Z = Z/0. () means the zero group.\n"
    
    print "Cohomology of X:"
    for i in K1.keys()[::-1]:
        print "Rank {0}:\t\t Z^{1}".format(
            i,
            len(K1[i]) + len(Q1.get(i,ZZ^1/ZZ^1).invariants()))

def print_all():
    """Prints the E2 page with generators of the groups."""
    pretty_print()
    print "\n"
    
    K1,K2,I2,I3 = setup_kernel()
    for i in K1.keys():
        print "Kernel of d1 in rank", i, ": ", K1[i]
        
    print "\n"

    for i in K2.keys():
        print "Kernel of d2 in rank", i, ": ", K2[i]
        print "\n"
        print "Image of d1 in rank", i, ": ", I2[i]
        print "\n"
        
    print "\n"
    
    for i in I3.keys():
        print "Image of d2 in rank", i, ": ", I3[i]

def print_mod():
    T = {0 : 1, 1: 3, 2: 3, 3: 1}
    (K1,K2,I2,I3) = setup_kernel()
    (Q1,Q2) = calculate_quotients()
    R1 = {}
    R2 = {}
    R3 = {}
    for i in H1.keys():
        R1[i] = len(H1[i])
        R2[i] = len(H2.get(i,[]))
        R3[i] = len(H3.get(i,[]))
        for j in range(i):
            R1[i] -= T.get(i-j,0)*R1[j]
            if R2[i] != 0:
                R2[i] -= T.get(i-j,0)*R2[j]
            if R3[i] != 0:
                R3[i] -= T.get(i-j,0)*R3[j]
    print "E1 for X/T:\t H^i(U/T)\t H^(i+1)(V/T,U/T)\t H^(i+2)(X/T,V/T)"
    for i in H1.keys()[::-1]:
        print "Row {0}:\t\t {1:>3}\t\t\t {2:>3}\t\t\t{3:>3}".format(
            i,
            R1[i],
            R2.get(i+1,0),
            R3.get(i+2,0))
        
    print "Cohomology of X/T:"
    R1 = {}
    R2 = {}
    for i in K1.keys():
        R1[i] = len(K1[i])
        R2[i] = len(Q1.get(i,ZZ^1/ZZ^1).invariants())
        for j in range(i):
            R1[i] -= T.get(i-j,0)*R1[j]
            if R2[i] != 0:
                R2[i] -= T.get(i-j,0)*R2[j]
    for i in R1.keys()[::-1]:
        print "Rank {0}:\t\t {1}{2}".format(
            i,
            "Z^{0}".format(R1[i]) if R1[i] != 0 else "",
            " + Z^{0}".format(R2[i]) if R2[i]!= 0 else "")
