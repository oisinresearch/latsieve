import sys
from sage.all_cmdline import *   # import sage library
cimport cython
from libc.stdlib cimport malloc, free

cdef extern from "renumber.c":
    int renumber(int vp, int vr, int* R, int L)

_sage_const_3 = Integer(3); _sage_const_2 = Integer(2); _sage_const_1 = Integer(1); _sage_const_2016590912589408 = Integer(2016590912589408); _sage_const_6 = Integer(6); _sage_const_4 = Integer(4); _sage_const_1268 = Integer(1268); _sage_const_1574794321679839 = Integer(1574794321679839); _sage_const_3514 = Integer(3514); _sage_const_76 = Integer(76); _sage_const_2235 = Integer(2235); _sage_const_1015 = Integer(1015)
K = NumberField(x**_sage_const_2 +_sage_const_3 *x+_sage_const_1 , names=('y',)); (y,) = K._first_ngens(1)
R = K['x']; (x,) = R._first_ngens(1)
L1 = K.extension(_sage_const_2016590912589408 *y*x**_sage_const_2  + _sage_const_1574794321679839 *x + _sage_const_2016590912589408 *y, names=('x1',)); (x1,) = L1._first_ngens(1)
L2 = K.extension((_sage_const_3 *y + _sage_const_1 )*x**_sage_const_4  + (_sage_const_6 *y + _sage_const_76 )*x**_sage_const_2  + (_sage_const_3 *y + _sage_const_1 ), names=('x2',)); (x2,) = L2._first_ngens(1)
A=(_sage_const_1015 *y - _sage_const_1268 )*x1 + (-_sage_const_2235 *y + _sage_const_3514 )
J1=L1.ideal(_sage_const_1 ,x1)**-_sage_const_1 
J2=L2.ideal(_sage_const_1 ,x2)**-_sage_const_1 
A=A*J1
#A.is_integral()
#A.absolute_norm();
#print A.factor()

off0 = 2
off1 = 2

if len(sys.argv) != 5:
    print("You must supply the renumber filename, the rerenumber filename, an input " +
    "filename and an output filename")
    sys.exit()

rnifilename = sys.argv[1]
rrnifilename = sys.argv[2]
inputfilename = sys.argv[3]
outputfilename = sys.argv[4]

with open(rnifilename) as f:
    List0 = []
    List1 = []
    i = 0
    for line in f:
        if i == 0:
            v = [int(x) for x in line.split(' ')] # read first line
            p1max = v[8]
            i1max = v[9]
            off1 = i1max - 1
        if i >= 3 and i < i1max:
            List0.append(int(line,16))
        elif i >= i1max:
            List1.append(int(line,16))
        i += 1

with open(rrnifilename) as f:
    List2 = []
    i = 0
    for line in f:
        if i == 0:
            imax = int(line) # read first line
            off2 = 1
        if i >= 1:
            List2.append(int(line,16))
        i += 1

cdef int *list0
cdef int *list1
cdef int *list2

list0 = <int *>malloc(len(List0)*cython.sizeof(int))
list1 = <int *>malloc(len(List1)*cython.sizeof(int))
list2 = <int *>malloc(len(List2)*cython.sizeof(int))

# populate c arrays
len0 = len(List0)
len1 = len(List1)
len2 = len(List2)
for i in range(0,len0):
    list0[i] = int(List0[i])
for i in range(0,len1):
    list1[i] = int(List1[i])
for i in range(0,len2):
    list2[i] = int(List2[i])

t = 0
with open(inputfilename,'r') as inputfile, open(outputfilename,'w') as outputfile:
    for line in inputfile:
        abcd = line.split(":")[0]
        [a,b,c,d]=[Integer(abcd.split(",")[i]) for i in (ellipsis_range(0,Ellipsis,3))]
        A = ((a+b*y)+(c+d*y)*x1)*J1
        #print [a,b,c,d]
        AF = A.factor()
        V = [0,1]
        for k in range(0, len(AF)):
            [Ip,e] = AF[k]
            [p,f] = Ip.gens_two()
            p = Integer(p)
            f1 = f.list()
            [r0,r1,r2,r3] = [f1[0][0]%p,f1[0][1]%p,f1[1][0]%p,f1[1][1]%p]
            if p <= p1max:
                r = r0+r1*p+r2*p**2+r3*p**3
                [vp,vr] = [2*p**4+1,r]
                vn = off0+renumber(vp,vr,list0,len0)
            else:
                r = r0+r1+r2+r3
                [vp,vr] = [8*p+1,r]
                vn = off1+renumber(vp,vr,list1,len1)
            V.extend([vn]*e)
        A = ((a+b*y)+(c+d*y)*x2)*J2
        AF = A.factor()
        for k in range(0, len(AF)):
            [Ip,e] = AF[k]
            [p,f] = Ip.gens_two()
            p = Integer(p)
            f1 = f.list()
            [r0,r1,r2,r3] = [f1[0][0]%p,f1[0][1]%p,f1[1][0]%p,f1[1][1]%p]
            if p <= p1max:
                r = r0+r1*p+r2*p**2+r3*p**3
                [vp,vr] = [2*p**4+1,p**4+r]
                vn = off0+renumber(vp,vr,list0,len0)
            else:
                r = r0+r1+r2+r3
                [vp,vr] = [8*p+1,4*p+r]
                vn = off1+renumber(vp,vr,list1,len1)
            V.extend([vn]*e)
        V.sort()
        m = len(V)-1
        V = [format(V[i],'x') for i in (ellipsis_range(0,Ellipsis,m))]
        str1 = abcd + ":"+ ",".join(V)
        #print str1
        outputfile.write(str1+"\n")

free(list1)
free(list0)

