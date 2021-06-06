#include <stdio.h>

int renumber(int vp, int vr, int* R, int L)
{
	int min = 0;
	int max = L;
	int i = max;

	if (vp == R[0]) i = min;
	else {
		i =  min + (max - min)/2;
		
		while (1) {
			int i_bak = i;
			while (i > 0 && R[i-1] > R[i]) i--;

			if (i != min) {
				if (vp == R[i])
					break;
				else if (vp < R[i])
					max = i;
				else
					min = i;
			}
			else {
				i = i_bak + 1;
				while (i < max && R[i-1] > R[i]) i++;

				if (i != max) {
					if (vp == R[i])
						break;
					else if (vp < R[i])
						max = i;
					else
						min = i;
				}
				else { // shouldn't get here
					printf("ERROR:  DATA MISSING IN RENUMBER TABLE!\n");
					exit(1);
				}
			}

			i = min + (max - min)/2;
		}
	}

	// now we have the first element vp of a decreasing list
	int j = i;
	while (R[i+1] >= vr && R[i+1] < R[j]) i++;

	return i;
}

/*
{
	int a = 0; int b = L;
	int i = L/2;
	while (R[i] != vp || R[i-1] > R[i]) {
		i = (a+b)/2;
		if (R[i-1] > R[i]) {
			if (R[i] < vp) i--;
			else {
				b = i;
			}
		}
		else {
			if (R[i] < vp) a = i;
			else b = i;
		}
		printf("%d\n",i);
		if (i == 0 || a == b) break;
	}
	int j = i;
	while (R[i+1] >= vr && R[i+1] < R[j]) i++;

	return i;
}
*/
