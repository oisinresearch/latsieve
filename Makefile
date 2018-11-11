all:
	g++ intpoly.cc factorsmall.cc -DPOLROOTSMAIN -o polroots -O3 -std=c++11
	g++ L2lu64.cc factorsmall.cc intpoly.cc latsieve3d.cc -o latsieve3d -O3 -lgmp -lgmpxx -std=c++11
