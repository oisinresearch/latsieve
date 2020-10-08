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
	g++ mpz_poly.cpp L2lu64.cc L2lu128.cc factorsmall.cc intpoly.cc latsieve3dQ.cc -o latsieve3dQ -lgmp -lgmpxx -std=c++11 -fopenmp -O3
	g++ -o rotate_full -lgmp -lgmpxx -lm -lpari -lpthread rotate_full.cc -std=c++11 -g -I/ichec/work/ndmat031c/local/gmp-6.1.2/include -I/ichec/work/ndmat031c/local/pari-2.11.0/include -L/ichec/work/ndmat031c/local/gmp-6.1.2/lib -L/ichec/work/ndmat031c/local/pari-2.11.0/lib -O0 -g
	g++ -o deduce_full -lgmp -lgmpxx -lm -lpari -lpthread deduce_full.cc -std=c++11 -g -I/ichec/work/ndmat031c/local/gmp-6.1.2/include -I/ichec/work/ndmat031c/local/pari-2.11.0/include -L/ichec/work/ndmat031c/local/gmp-6.1.2/lib -L/ichec/work/ndmat031c/local/pari-2.11.0/lib -O0 -g
	g++ mpz_poly.cpp L2lu64.cc factorsmall.cc intpoly.cc latsieve2d.cc -o latsieve2d -lgmp -lgmpxx -std=c++11 -fopenmp -O3
	g++ mpz_poly.cpp L2lu64.cc factorsmall.cc intpoly.cc cflat.cc -o cflat -lgmp -lgmpxx -std=c++11 -fopenmp -O3
	g++ mpz_poly.cpp mpz_poly_bivariate.c L2lu64.cc factorsmall.cc intpoly.cc latsieve4d.cc -o latsieve4d -lgmp -lgmpxx -std=c++11 -fopenmp -O3
	g++ -o makefbtnfs -lgmp -lgmpxx makefbtnfs.cc mpz_poly_bivariate.c intpoly.cc mpz_poly.cpp factorsmall.cc -std=c++11 -fopenmp -O3
	g++ galoisexpand4d.cc mpz_poly_bivariate.c mpz_poly.cpp factorsmall.cc intpoly.cc -o galoisexpand4d -lgmp -lgmpxx -std=c++11 -O3

clean:
	rm -f cflat
	rm -f deduce
	rm -f deduce_full
	rm -f deduce_one_side
	rm -f descend
	rm -f encode3d
	rm -f galoisexpand
	rm -f latsieve2d
	rm -f latsieve3d
	rm -f latsieve3dQ
	rm -f latsieve4d
	rm -f makefb
	rm -f makefbtnfs
	rm -f polroots
	rm -f presort
	rm -f rareoutliers
	rm -f rootsofunity
	rm -f rotate_full
	rm -f galoisexpand4d
