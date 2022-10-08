import sys
from sage.all_cmdline import *   # import sage library
cimport cython
from libc.stdlib cimport malloc, free
from sage.misc.search import search

cdef extern from "renumber.c":
    int renumber(int vp, int vr, int* R, int L)

#_sage_const_3 = Integer(3); _sage_const_2 = Integer(2); _sage_const_1 = Integer(1); _sage_const_2016590912589408 = Integer(2016590912589408); _sage_const_6 = Integer(6); _sage_const_4 = Integer(4); _sage_const_1268 = Integer(1268); _sage_const_1574794321679839 = Integer(1574794321679839); _sage_const_3514 = Integer(3514); _sage_const_76 = Integer(76); _sage_const_2235 = Integer(2235); _sage_const_1015 = Integer(1015)
#K = NumberField(x**_sage_const_2 +_sage_const_3 *x+_sage_const_1 , names=('y',)); (y,) = K._first_ngens(1)
#R = K['x']; (x,) = R._first_ngens(1)
#L1 = K.extension(_sage_const_2016590912589408 *y*x**_sage_const_2  + _sage_const_1574794321679839 *x + _sage_const_2016590912589408 *y, names=('x1',)); (x1,) = L1._first_ngens(1)
#L2 = K.extension((_sage_const_3 *y + _sage_const_1 )*x**_sage_const_4  + (_sage_const_6 *y + _sage_const_76 )*x**_sage_const_2  + (_sage_const_3 *y + _sage_const_1 ), names=('x2',)); (x2,) = L2._first_ngens(1)
#A=(_sage_const_1015 *y - _sage_const_1268 )*x1 + (-_sage_const_2235 *y + _sage_const_3514 )
#J1=L1.ideal(_sage_const_1 ,x1)**-_sage_const_1 
#J2=L2.ideal(_sage_const_1 ,x2)**-_sage_const_1 
#A=A*J1
#A.is_integral()
#A.absolute_norm();
#print A.factor()

_sage_const_2 = Integer(2); _sage_const_1 = Integer(1); _sage_const_356 = Integer(356); _sage_const_257 = Integer(257); _sage_const_4 = Integer(4); _sage_const_5 = Integer(5)
K = NumberField(x**_sage_const_2 -x+_sage_const_1 ,'y', names=('y',)); (y,) = K._first_ngens(1)
R = K['x']; (x,) = R._first_ngens(1)
L1 = K.extension(_sage_const_356 *y*x**_sage_const_2  - _sage_const_257 *x + _sage_const_356 *y, 'x1', names=('x1',)); (x1,) = L1._first_ngens(1)
L2 = K.extension((y - _sage_const_1 )*x**_sage_const_4  + (_sage_const_2 *y - _sage_const_5 )*x**_sage_const_2  + (y - _sage_const_1 ), 'x2', names=('x2',)); (x2,) = L2._first_ngens(1)
J1 = L1.ideal(_sage_const_1 ,x1)**-_sage_const_1 
J2 = L2.ideal(_sage_const_1 ,x2)**-_sage_const_1 

if len(sys.argv) != 3:
    print("You must supply an input filename and an output filename")
    sys.exit()

inputfilename = sys.argv[1]
outputfilename = sys.argv[2]

hh = 18446744073709551616
list0 = []
g = 0
t = 0
with open(inputfilename,'r') as inputfile:
    for line in inputfile:
        [Astr,Bstr] = line.split(":")[0].split(",")
        [A,B]=[int(Astr,16),int(Bstr,16)]
        [a,b,c,d]=[(A>>24)-2**23,(A & (2**24-1))-2**23,(B>>24)-2**23,(B & (2**24-1))-2**23]
        abcd = str(a)+","+str(b)+","+str(c)+","+str(d)
        #abcd = line.split(":")[0]
        relstr = line.split(":")[1]
        #[a,b,c,d]=[Integer(abcd.split(",")[i]) for i in (ellipsis_range(0,Ellipsis,3))]
        A1 = ((a+b*y)+(c+d*y)*x1)*J1
        A2 = ((a+b*y)+(c+d*y)*x2)*J2
        N1 = A1.absolute_norm()
        N2 = A2.absolute_norm()
        if gcd(N1,N2) == 1:
            h = hash((a+b*y)/(c+d*y)) % hh
            list0.append([h,abs(N1*N2),abcd,relstr])
        else:
            g += 1
        t += 1
        if t % 10000 == 0:
            print str(t) + " relations read..."

print str(t) + " relations read.\n"
print str(g) + " relations with gcd(N1,N2) > 1."
print "\nfinished reading input.\n"
list0.sort()
h1 = 0
t = 0
with open(outputfilename,'w') as outputfile:
    for item in list0:
        h2 = item[0]
        if h2 != h1:
            h1 = h2
            abcd = item[2]
            relstr = item[3]
            outputfile.write(abcd+":"+relstr)
            t += 1
            if t % 10000 == 0:
                print str(t) + " relations written..."

print str(t) + " relations written.\n"

