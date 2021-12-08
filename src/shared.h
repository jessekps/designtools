#ifndef DESIGNTOOLS_SHARED_
#define DESIGNTOOLS_SHARED_
#include <RcppArmadillo.h>

inline void swp(int& a, int& b)
{
	const int c=a;
	a=b;
	b=c;
};
void swp_pth(int* pth, int& s1, int& s2);
void sample2(int& s1, int& s2, const int n);


#endif