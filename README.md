# Utilities for higher-dimensional sieving for exTNFS
PhD work of Oisin Robinson, see for example [arXiv:2212.04999 [math.cs.CR]](https://arxiv.org/abs/2212.04999) (2024)

### Overview
This is the raw git project maintained during my PhD work.  I used this to set a
computational number theory record in fully breaking the discrete log in a finite
field of shape F_{p^4}, of size 512 bits (see https://dldb.loria.fr/?filter=all&sort=date).

This code needs to be tidied heavily, but in theory provides production-ready 3d/4d
orthotope sieves and an experimental 6d 'pruned enumeration' sieve (working and
tested but 'unproven').

There are many other utilities I wrote, for various routines required in completing
a record size discrete log computation using exTNFS, that are included here, but
documentation on these at present is severely lacking, with barely more than a program
parameter list available in many cases.

Other software that this repository uses/depends on include GNU-MP, Pari-gp and Sage, and
last but not least, CADO-NFS.

I plan to work on this repository and fine-tune it so that it is more presentable, but
for the moment, this is the raw original.

### License
&copy; 2024, Oisin Robinson.
This work is published under the GNU General Pupblic License (GPL) v3.0.
SEE the LICENSE file for complete statement.

