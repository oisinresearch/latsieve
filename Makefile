all:
	g++ intpoly.cc factorsmall.cc -DPOLROOTSMAIN -o polroots -O3 -std=c++11
