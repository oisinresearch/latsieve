
# This file was *autogenerated* from the file renumber_file.sage
from sage.all_cmdline import *   # import sage library

_sage_const_3 = Integer(3); _sage_const_2 = Integer(2); _sage_const_1 = Integer(1); _sage_const_2016590912589408 = Integer(2016590912589408); _sage_const_6 = Integer(6); _sage_const_4 = Integer(4); _sage_const_1268 = Integer(1268); _sage_const_8 = Integer(8); _sage_const_1574794321679839 = Integer(1574794321679839); _sage_const_3514 = Integer(3514); _sage_const_76 = Integer(76); _sage_const_10 = Integer(10); _sage_const_2235 = Integer(2235); _sage_const_1015 = Integer(1015); _sage_const_0 = Integer(0); _sage_const_50 = Integer(50)
import datetime
import sys

K = NumberField(x**_sage_const_2 +_sage_const_3 *x+_sage_const_1 ,'y', names=('y',)); (y,) = K._first_ngens(1)
R = K['x']; (x,) = R._first_ngens(1)
L1 = K.extension(_sage_const_2016590912589408 *y*x**_sage_const_2  + _sage_const_1574794321679839 *x + _sage_const_2016590912589408 *y,'x1', names=('x1',)); (x1,) = L1._first_ngens(1)
L2 = K.extension((_sage_const_3 *y + _sage_const_1 )*x**_sage_const_4  + (_sage_const_6 *y + _sage_const_76 )*x**_sage_const_2  + (_sage_const_3 *y + _sage_const_1 ),'x2', names=('x2',)); (x2,) = L2._first_ngens(1)
A = (_sage_const_1015 *y - _sage_const_1268 )*x1 + (-_sage_const_2235 *y + _sage_const_3514 )
J1 = L1.ideal(_sage_const_1 ,x1)**-_sage_const_1 
A = A*J1

if len(sys.argv) != _sage_const_4 :
    print "You must supply a prime upper bound, total number t of cpus, and cpu number between 0 and t-1"
    sys.exit()

t = Integer(sys.argv[_sage_const_2 ])
cpu = Integer(sys.argv[_sage_const_3 ])
B = Integer(sys.argv[_sage_const_1 ])
start = floor(B*cpu/t)
if start == _sage_const_0 :
    start = _sage_const_2 
finish = ceil(B*(cpu+_sage_const_1 )/t)

filename = 'renfp4.' + str(cpu)
output = open(filename,'a')
pc = _sage_const_0 
tenpc = floor(B/t/_sage_const_10 )
for p in range(start,finish):
    if p % tenpc == _sage_const_0 :
        pc += _sage_const_10 
        print pc,'%'
    if is_prime(p):
        if p < _sage_const_50 :
            vp = _sage_const_2 *p**_sage_const_4 +_sage_const_1 
            output.write(vp.hex()+'\n')
            pi = L2.ideal(p).factor()
            n = len(pi)
            vri = []
            for i in range(_sage_const_0 ,n):
                [Ip,e] = pi[i]
                f = Ip.gens_two()[_sage_const_1 ]
                f1 = f.list()
                [r0,r1,r2,r3] = [f1[_sage_const_0 ][_sage_const_0 ]%p,f1[_sage_const_0 ][_sage_const_1 ]%p,f1[_sage_const_1 ][_sage_const_0 ]%p,f1[_sage_const_1 ][_sage_const_1 ]%p]
                vr = p**_sage_const_4 +r0+r1*p+r2*p**_sage_const_2 +r3*p**_sage_const_3 
                vri.append(vr)
            pi = L1.ideal(p).factor()
            n = len(pi)
            for i in range(_sage_const_0 ,n):
                [Ip,e] = pi[i]
                f = Ip.gens_two()[_sage_const_1 ]
                f1 = f.list()
                [r0,r1,r2,r3] = [f1[_sage_const_0 ][_sage_const_0 ]%p,f1[_sage_const_0 ][_sage_const_1 ]%p,f1[_sage_const_1 ][_sage_const_0 ]%p,f1[_sage_const_1 ][_sage_const_1 ]%p]
                vr = r0+r1*p+r2*p**_sage_const_2 +r3*p**_sage_const_3 
                vri.append(vr)
            vri.sort(reverse=True)
            for i in range(_sage_const_1 ,len(vri)):
                output.write(vri[i].hex()+'\n')
        else:
            vp = _sage_const_8 *p+_sage_const_1 
            output.write(vp.hex()+'\n')
            pi = L2.ideal(p).factor()
            n = len(pi)
            vri = []
            for i in range(_sage_const_0 ,n):
                [Ip,e] = pi[i]
                f = Ip.gens_two()[_sage_const_1 ]
                f1 = f.list()
                [r0,r1,r2,r3] = [f1[_sage_const_0 ][_sage_const_0 ]%p,f1[_sage_const_0 ][_sage_const_1 ]%p,f1[_sage_const_1 ][_sage_const_0 ]%p,f1[_sage_const_1 ][_sage_const_1 ]%p]
                vr = _sage_const_4 *p+r0+r1+r2+r3
                vri.append(vr)
            pi = L1.ideal(p).factor()
            n = len(pi)
            for i in range(_sage_const_0 ,n):
                [Ip,e] = pi[i]
                f = Ip.gens_two()[_sage_const_1 ]
                f1 = f.list()
                [r0,r1,r2,r3] = [f1[_sage_const_0 ][_sage_const_0 ]%p,f1[_sage_const_0 ][_sage_const_1 ]%p,f1[_sage_const_1 ][_sage_const_0 ]%p,f1[_sage_const_1 ][_sage_const_1 ]%p]
                vr = r0+r1+r2+r3
                vri.append(vr)
            vri.sort(reverse=True)
            for i in range(_sage_const_1 ,len(vri)):
                output.write(vri[i].hex()+'\n')

