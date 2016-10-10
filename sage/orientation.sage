
var('a','l','t','f','h','b','c','d', 'g')

small = 10^(-10000)

print
print "##########################"
print "#                        #"
print "#  # #       #   # #     #"
print "#   #     #####   #      #"
print "#  # #001    #   # #011  #"
print "#                        #"
print "##########################"
print "\n"

x = b*h/d
y = c*d/h

# (d,1-bc)
xd = x.subs(b=1,c=1-l)
yd = y.subs(b=1,c=1-l)
firstd = (xd.subs(d=e^(i*t),
                  l=-1).diff(t).subs(t=0).mul(-i),
          yd.subs(d=e^(i*t),
                  l=-1).diff(t).subs(t=0).mul(-i))
secondd = (xd.subs(d=1,
                   l=-e^(i*t)).diff(t).subs(t=0).mul(-i),
           yd.subs(d=1,
                   l=-e^(i*t)).diff(t).subs(t=0).mul(-i))
resd = Matrix((firstd,secondd)).determinant().expand()
print "For (d,1-bc):"
print resd
print "Preserves orientation:", resd.subs(h=small).is_positive()
print

# (1-bc,h)
xh = x.subs(d=1,b=1,c=1-l)
yh = y.subs(d=1,b=1,c=1-l)
firsth = (xh.subs(l=-e^(i*t)).diff(t).subs(t=0).mul(-i),
          yh.subs(l=-e^(i*t)).diff(t).subs(t=0).mul(-i))
secondh = (xh.subs(c=2,h=h*e^(i*t)).diff(t).subs(t=0).mul(-i),
           yh.subs(c=2,h=h*e^(i*t)).diff(t).subs(t=0).mul(-i))
resh = Matrix((firsth,secondh)).determinant().expand()
print "For (1-bc,h):"
print resh
print "Preserves orientation: ", resh.subs(h=small).is_positive()
print

# b=c=2
xbc = x.subs(d=1)
ybc = y.subs(d=1)
firstbc  = (xbc.subs(b=2*e^(i*t),
                     c=2).diff(t).subs(t=0).mul(-i),
            ybc.subs(b=2*e^(i*t),
                     c=2).diff(t).subs(t=0).mul(-i))
secondbc = (xbc.subs(c=2*e^(i*t),
                     b=2).diff(t).subs(t=0).mul(-i),
            ybc.subs(c=2*e^(i*t),
                     b=2).diff(t).subs(t=0).mul(-i))
resbc = Matrix((firstbc,secondbc)).determinant().expand()
print "For b=c=2:"
print resbc
print "Preserves orientation: ", resbc.subs(h=small).is_positive()
print

print 
print "##########################"
print "#                        #"
print "#  # #       #   # #     #"
print "#   #     #####   #      #"
print "#  # #001    #   # #111  #"
print "#                        #"
print "##########################"
print "\n"

x = -b*a*h/(g*d)
y = -c*d*g/(a*h*(1-b*c))+1/(1-b*c)

# b = c = 2
xbc = x.subs(a=1,d=1)
ybc = y.subs(a=1,d=1)
firstbc = (xbc.subs(b=-2*e^(i*t),
                    c=-2).diff(t).subs(t=0).mul(-i),
           ybc.subs(b=-2*e^(i*t),
                    c=-2).diff(t).subs(t=0).mul(-i))
secondbc = (xbc.subs(c=-2*e^(i*t),
                     b=-2).diff(t).subs(t=0).mul(-i),
            ybc.subs(c=-2*e^(i*t),
                     b=-2).diff(t).subs(t=0).mul(-i))
resbc = Matrix((firstbc,secondbc)).determinant().expand()
print "For b=c=2:"
print resbc
print "Preserves orientation:", resbc.subs(g=small,h=small).is_positive()
print

x = x.subs(b=1,c=1-l)
y = y.subs(b=1,c=1-l)

# a,l
xa = x.subs(d=1)
ya = y.subs(d=1)
firsta = (xa.subs(a=-e^(i*t),
                  l=-1).diff(t).subs(t=0).mul(-i),
          ya.subs(a=-e^(i*t),
                  l=-1).diff(t).subs(t=0).mul(-i))
seconda = (xa.subs(l=-e^(i*t),
                   a=-1).diff(t).subs(t=0).mul(-i),
           ya.subs(l=-e^(i*t),
                   a=-1).diff(t).subs(t=0).mul(-i))
resa = Matrix((firsta,seconda)).determinant().expand()
print "For (a,1-bc)"
print resa
print "Preserves orientation:", resa.subs(g=small,h=small).is_positive()
print

# d,l
xa = x.subs(a=1)
ya = y.subs(a=1)
firsta = (xa.subs(d=-e^(i*t),
                  l=-1).diff(t).subs(t=0).mul(-i),
          ya.subs(d=-e^(i*t),
                  l=-1).diff(t).subs(t=0).mul(-i))
seconda = (xa.subs(l=-e^(i*t),
                   d=-1).diff(t).subs(t=0).mul(-i),
           ya.subs(l=-e^(i*t),
                   d=-1).diff(t).subs(t=0).mul(-i))
resa = Matrix((firsta,seconda)).determinant().expand()
print "For (d,1-bc)"
print resa
print "Preserves orientation:", resa.subs(g=small,h=small).is_positive()
print

# g,l 
xa = x.subs(a=1,d=1)
ya = y.subs(d=1,a=1)
firsta = (xa.subs(g=-g*e^(i*t),
                  l=-1).diff(t).subs(t=0).mul(-i),
          ya.subs(g=-g*e^(i*t),
                  l=-1).diff(t).subs(t=0).mul(-i))
seconda = (xa.subs(l=-e^(i*t),
                   g=-g).diff(t).subs(t=0).mul(-i),
           ya.subs(l=-e^(i*t),
                   g=-g).diff(t).subs(t=0).mul(-i))
resa = Matrix((firsta,seconda)).determinant().expand()
print "For (g,1-bc)"
print resa
print "Preserves orientation:", resa.subs(g=small,h=small).is_positive()
print

# l,h 
xa = x.subs(a=1,d=1)
ya = y.subs(d=1,a=1)
seconda = (xa.subs(h=-g*e^(i*t),
                   l=-1).diff(t).subs(t=0).mul(-i),
           ya.subs(h=-g*e^(i*t),
                   l=-1).diff(t).subs(t=0).mul(-i))
firsta = (xa.subs(l=-e^(i*t),
                  h=-h).diff(t).subs(t=0).mul(-i),
          ya.subs(l=-e^(i*t),
                  h=-h).diff(t).subs(t=0).mul(-i))
resa = Matrix((firsta,seconda)).determinant().expand()
print "For (1-bc,h)"
print resa
print "Preserves orientation:", resa.subs(g=small,h=small).is_positive()


print 
print "##########################"
print "#                        #"
print "#  # #       #   # #     #"
print "#   #     #####   #      #"
print "#  # #101    #   # #111  #"
print "#                        #"
print "##########################"
print "\n"

C = -a*f-a^2*b*h/d+a^2*b*h*f/d+a^2*b^2*c*h/d
D = d*(1-b*c)/h
x = b*C/D
y = c*d/(h*C)-a/C

# (a,1-bc)
xa = x.subs(b=1,c=1-l,d=1,f=1)
ya = y.subs(b=1,c=1-l,d=1,f=1)
firsta = (xa.subs(a=e^(i*t),
                  l=-1).diff(t).subs(t=0).mul(-i),
          ya.subs(a=e^(i*t),
                  l=-1).diff(t).subs(t=0).mul(-i))
seconda = (xa.subs(a=1,
                   l=e^(i*t)).diff(t).subs(t=pi).mul(-i),
           ya.subs(a=1,
                   l=e^(i*t)).diff(t).subs(t=pi).mul(-i))
resa = Matrix((firsta,seconda)).determinant().expand()
print "For (a,1-bc):"
print resa
print "Preserves orientation: ", resa.subs(h=small).is_positive()
print

# (d, 1-bc)
xd = x.subs(b=1,c=1-l,a=1,f=1)
yd = y.subs(b=1,c=1-l,a=1,f=1)
firstd = (xd.subs(d=e^(i*t),
                  l=-1).diff(t).subs(t=0).mul(-i),
          yd.subs(d=e^(i*t),
                  l=-1).diff(t).subs(t=0).mul(-i))
secondd = (xd.subs(d=1,
                   l=-e^(i*t)).diff(t).subs(t=0).mul(-i),
           yd.subs(d=1,
                   l=-e^(i*t)).diff(t).subs(t=0).mul(-i))
resd = Matrix((firstd,secondd)).determinant().expand()
print "For (d,1-bc)"
print resd
print "Preserves orientation: ", resd.subs(h=small).is_positive()
print

# (f,1-bc)
xf = x.subs(b=1,c=1-l,a=1,d=1)
yf = y.subs(b=1,c=1-l,a=1,d=1)
firstf = (xf.subs(f=e^(i*t),
                  l=-1).diff(t).subs(t=0).mul(-i),
          yf.subs(f=e^(i*t),
                  l=-1).diff(t).subs(t=0).mul(-i))
secondf = (xf.subs(f=1,
                   l=-e^(i*t)).diff(t).subs(t=0).mul(-i),
           yf.subs(f=1,
                   l=-e^(i*t)).diff(t).subs(t=0).mul(-i))
resf = Matrix((firstf,secondf)).determinant().expand()
print "For (f,1-bc)"
print resf
print "Preserves orientation: ", resf.subs(h=small).is_positive()
print

# (1-bc,h)
xh = x.subs(b=1,c=1-l,a=1,d=1,f=1)
yh = y.subs(b=1,c=1-l,a=1,d=1,f=1)
var('epsilon')
firsth = (xh.subs(l=e^(i*t),
                  h=epsilon).diff(t).subs(t=pi).mul(-i),
          yh.subs(l=e^(i*t),
                  h=epsilon).diff(t).subs(t=pi).mul(-i))
secondh = (xh.subs(h=epsilon*e^(i*t),
                   l=-1).diff(t).subs(t=0).mul(-i),
           yh.subs(h=epsilon*e^(i*t),
                   l=-1).diff(t).subs(t=0).mul(-i))
resh = Matrix((firsth,secondh)).determinant().expand()
print "For (1-bc,h)"
print resh
print "Preserves orientation: ", resh.subs(epsilon=small).is_positive()
print

# b=c=2
xbc = x.subs(a=1,d=1,f=1)
ybc = y.subs(a=1,d=1,f=1)
firstbc = (xbc.subs(b=2*e^(i*t),
                    c=2).diff(t).subs(t=0).mul(-i),
           ybc.subs(b=2*e^(i*t),
                    c=2).diff(t).subs(t=0).mul(-i))
secondbc = (xbc.subs(c=2*e^(i*t),
                     b=2).diff(t).subs(t=0).mul(-i),
            ybc.subs(c=2*e^(i*t),
                     b=2).diff(t).subs(t=0).mul(-i))
resbc = Matrix((firstbc,secondbc)).determinant().expand()
print "For b=c=2"
print resbc
print "Preserves orientation: ", resbc.subs(h=small).is_positive()
print


print 
print "##########################"
print "#                        #"
print "#  # #       #   # #     #"
print "#   #     #####   #      #"
print "#  # #011    #   # #111  #"
print "#                        #"
print "##########################"
print "\n"

C = c-a/g-c*d/f-b*c^2+a*b*c/g
x = b*f*C/d
y = c/C-a/(g*C)

xt = x.subs(a=1,f=1,d=1)
yt = y.subs(a=1,f=1,d=1)
firstt = (xt.subs(b=2*e^(i*t),
                  c=2).diff(t).subs(t=0).mul(-i),
          yt.subs(b=2*e^(i*t),
                  c=2).diff(t).subs(t=0).mul(-i))
secondt = (xt.subs(b=2,
                   c=2*e^(i*t)).diff(t).subs(t=0).mul(-i),
           yt.subs(b=2,
                   c=2*e^(i*t)).diff(t).subs(t=0).mul(-i))
rest = Matrix((firstt,secondt)).determinant().expand()
print "For b=c=2:"
print rest
print "Preserves orientation: ", rest.subs(g=small).is_positive()

x = x.subs(c=1,b=1-l)
y = y.subs(c=1,b=1-l)

# (a,l)
xa = x.subs(d=1,f=1)
ya = y.subs(d=1,f=1)
firsta = (xa.subs(a=e^(i*t),
                  l=-1).diff(t).subs(t=0).mul(-i),
          ya.subs(a=e^(i*t),
                  l=-1).diff(t).subs(t=0).mul(-i))
seconda = (xa.subs(a=1,
                   l=-e^(i*t)).diff(t).subs(t=0).mul(-i),
           ya.subs(a=1,
                   l=-e^(i*t)).diff(t).subs(t=0).mul(-i))
resa = Matrix((firsta,seconda)).determinant().expand()
print "For (a,1-bc):"
print resa
print "Preserves orientation: ", resa.subs(g=small).is_positive()
print

xa = x.subs(a=1,f=1)
ya = y.subs(a=1,f=1)
firsta = (xa.subs(d=e^(i*t),
                  l=-1).diff(t).subs(t=0).mul(-i),
          ya.subs(d=e^(i*t),
                  l=-1).diff(t).subs(t=0).mul(-i))
seconda = (xa.subs(d=1,
                   l=-e^(i*t)).diff(t).subs(t=0).mul(-i),
           ya.subs(d=1,
                   l=-e^(i*t)).diff(t).subs(t=0).mul(-i))
resa = Matrix((firsta,seconda)).determinant().expand()
print "For (d,1-bc):"
print resa
print "Preserves orientation: ", resa.subs(g=small).is_positive()
print

xa = x.subs(d=1,a=1)
ya = y.subs(d=1,a=1)
firsta = (xa.subs(f=e^(i*t),
                  l=-1).diff(t).subs(t=0).mul(-i),
          ya.subs(f=e^(i*t),
                  l=-1).diff(t).subs(t=0).mul(-i))
seconda = (xa.subs(f=1,
                   l=-e^(i*t)).diff(t).subs(t=0).mul(-i),
           ya.subs(f=1,
                   l=-e^(i*t)).diff(t).subs(t=0).mul(-i))
resa = Matrix((firsta,seconda)).determinant().expand()
print "For (f,1-bc):"
print resa
print "Preserves orientation: ", resa.subs(g=small).is_positive()
print

xa = x.subs(d=1,f=1,a=1)
ya = y.subs(d=1,f=1,a=1)
seconda = (xa.subs(g=g*e^(i*t),
                   l=-1).diff(t).subs(t=0).mul(-i).expand(),
           ya.subs(g=g*e^(i*t),
                   l=-1).diff(t).subs(t=0).mul(-i).expand())
firsta = (xa.subs(l=-e^(i*t)).diff(t).subs(t=0).mul(-i),
          ya.subs(l=-e^(i*t)).diff(t).subs(t=0).mul(-i))
resa = Matrix((firsta,seconda)).determinant().expand()
print "For (1-bc,g):"
print resa
print "Preserves orientation: ", resa.subs(g=small).is_positive()
print
