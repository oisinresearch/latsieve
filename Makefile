all:
	g++ intpoly.cc factorsmall.cc -DPOLROOTSMAIN -o polroots -O3 -std=c++11
	g++ mpz_poly.cpp factorsmall.cc intpoly.cc presort.cc -o presort -lgmp -lgmpxx -std=c++11 -g
	g++ encode3d.cc -o encode3d -O3
	g++ -o descend -lgmp -lgmpxx -lm -lpari -lpthread -O0 descend.cc -std=c++11 -g -I/usr/local/include -L/usr/local/lib
	g++ mpz_poly.cpp factorsmall.cc intpoly.cc galoisexpand.cc -o galoisexpand -lgmp -lgmpxx -std=c++11 -g
	g++ mpz_poly.cpp L2lu64.cc factorsmall.cc intpoly.cc latsieve3d.cc -o latsieve3d -lgmp -lgmpxx -std=c++11 -fopenmp -O3
	g++ -o makefb -lgmp -lgmpxx makefb.cpp intpoly.cc mpz_poly.cpp factorsmall.cc -std=c++11 -fopenmp -O3
	g++ rootsofunity.cc -o rootsofunity -O3
	g++ -o deduce -lgmp -lgmpxx -lm -lpari -lpthread deduce.cc -std=c++11 -g -I/usr/local/include -L/usr/local/lib -O0 -g
	g++ -o deduce_one_side -lgmp -lgmpxx -lm -lpari -lpthread deduce_one_side.cc -std=c++11 -g -I/usr/local/include -L/usr/local/lib -O0 -g
	g++ -o deduce_full -lgmp -lgmpxx -lm -lpari -lpthread deduce_full.cc -std=c++11 -g -I/usr/local/include -L/usr/local/lib -O0 -g
