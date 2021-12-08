#include "shared.h"

void swp_pth(int* pth, int& s1, int& s2)
{
	if(s1>s2) swp(s1,s2);
	const int n = s1 + (s2-s1 + 1)/2;
	for(int i=s1, i2=s2; i<n; i++, i2--)
	{
		int a = *(pth+i);
		*(pth+i) = *(pth+i2);
		*(pth+i2) = a;		
	}	
}

void sample2(int& s1, int& s2, const int n)
{
	s1 = (int)R::runif(0,n);
	s2 = (int)R::runif(0,n-1);
	s2 += (s2 >= s1);
}