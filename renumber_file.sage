import datetime
import sys

K.<y> = NumberField(x^2+3*x+1,'y')
R.<x> = K[]
L1.<x1> = K.extension(2016590912589408*y*x^2 + 1574794321679839*x + 2016590912589408*y,'x1')
L2.<x2> = K.extension((3*y + 1)*x^4 + (6*y + 76)*x^2 + (3*y + 1),'x2')
A = (1015*y - 1268)*x1 + (-2235*y + 3514)
J1 = L1.ideal(1,x1)^-1
A = A*J1

if len(sys.argv) != 4:
    print "You must supply a prime upper bound, total number t of cpus, and cpu number between 0 and t-1"
    sys.exit()

t = Integer(sys.argv[2])
cpu = Integer(sys.argv[3])
B = Integer(sys.argv[1])
start = floor(B*cpu/t)
if start == 0:
    start = 2
finish = ceil(B*(cpu+1)/t)

filename = 'renfp4.' + str(cpu)
output = open(filename,'a')
pc = 0
tenpc = floor(B/t/10)
for p in range(start,finish):
    if p % tenpc == 0:
        pc += 10
        print pc,'%'
    if is_prime(p):
        if p < 50:
            vp = 2*p^4+1
            output.write(vp.hex()+'\n')
            pi = L2.ideal(p).factor()
            n = len(pi)
            vri = []
            for i in range(0,n):
                [Ip,e] = pi[i]
                f = Ip.gens_two()[1]
                f1 = f.list()
                [r0,r1,r2,r3] = [f1[0][0]%p,f1[0][1]%p,f1[1][0]%p,f1[1][1]%p]
                vr = p^4+r0+r1*p+r2*p^2+r3*p^3
                vri.append(vr)
            pi = L1.ideal(p).factor()
            n = len(pi)
            for i in range(0,n):
                [Ip,e] = pi[i]
                f = Ip.gens_two()[1]
                f1 = f.list()
                [r0,r1,r2,r3] = [f1[0][0]%p,f1[0][1]%p,f1[1][0]%p,f1[1][1]%p]
                vr = r0+r1*p+r2*p^2+r3*p^3
                vri.append(vr)
            vri.sort(reverse=True)
            for i in range(1,len(vri)):
                output.write(vri[i].hex()+'\n')
        else:
            vp = 8*p+1
            output.write(vp.hex()+'\n')
            pi = L2.ideal(p).factor()
            n = len(pi)
            vri = []
            for i in range(0,n):
                [Ip,e] = pi[i]
                f = Ip.gens_two()[1]
                f1 = f.list()
                [r0,r1,r2,r3] = [f1[0][0]%p,f1[0][1]%p,f1[1][0]%p,f1[1][1]%p]
                vr = 4*p+r0+r1+r2+r3
                vri.append(vr)
            pi = L1.ideal(p).factor()
            n = len(pi)
            for i in range(0,n):
                [Ip,e] = pi[i]
                f = Ip.gens_two()[1]
                f1 = f.list()
                [r0,r1,r2,r3] = [f1[0][0]%p,f1[0][1]%p,f1[1][0]%p,f1[1][1]%p]
                vr = r0+r1+r2+r3
                vri.append(vr)
            vri.sort(reverse=True)
            for i in range(1,len(vri)):
                output.write(vri[i].hex()+'\n')
