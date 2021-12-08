#include <RcppArmadillo.h>
#include <omp.h>
#include "shared.h"

using namespace arma;

using Rcpp::List;
using Rcpp::Named;

int omp_thread_count() 
{
    int n = 0;
    #pragma omp parallel reduction(+:n)
    n += 1;
    return n;
}

/*
complete is O(n!)

optimized version is usable for n <= 13

n=13
m=matrix(sample(1:100,n**2,TRUE),n,n)
m[lower.tri(m)] = t(m)[lower.tri(m)] 

system.time({a=designtools:::tsp_complete(m)})
# 6.64
system.time({a=designtools:::tsp_complete2(m)})
# 1.28


# linear time / reasonably good / non deterministic

sa = designtools:::tsp_sa(m)
*/


// basic complete version

//[[Rcpp::export]]
List tsp_complete(const arma::imat& m)
{
	const int n = m.n_rows;
	int cost, min_cost = accu(m);;
	std::vector<int> pth(n), shortest_path(n);
	std::iota(pth.begin(), pth.end(), 0);
	
	do {
		cost = m.at(pth[0],pth[n-1]);
		for(int i=1; i<n; i++)
			cost += m.at(pth[i],pth[i-1]);
		if(cost < min_cost)
		{
			min_cost = cost;
			shortest_path = pth;
		}
    } while(std::next_permutation(pth.begin()+1, pth.end()));
	
	return List::create(Named("path") = shortest_path, Named("cost") = min_cost);
}


// optimized, same O as complete
// n needs to be >= 3

//[[Rcpp::export]]
List tsp_complete2(const arma::imat& m)
{
	const int n=m.n_rows;
	const int n2 = n-2, nth = omp_thread_count();
	
	ivec min_cost(nth, fill::value(accu(m)));
	imat paths(n, nth, fill::zeros);
	
#pragma omp parallel
	{
		ivec pth(n2);
		const int th = omp_get_thread_num();

#pragma omp for
		for(int nd2=1; nd2<n-1; nd2++)
		{
			int j=0;
			for(int i=1; i<n; i++) if(i!=nd2) pth[j++] = i;
			
			do {
				if(pth[n2-1] < nd2) continue;
				int cost = m.at(0,pth[n2-1]) + m.at(0,nd2) + m.at(nd2,pth[0]);
				for(int i=1; i<n2; i++)
					cost += m.at(pth[i],pth[i-1]);
				if(cost < min_cost[th])
				{
					min_cost[th] = cost;
					paths.at(1,th) = nd2;
					for(int i=0;i<n2;i++)
						paths.at(i+2,th) = pth[i];						
				} 			
			} while(std::next_permutation(pth.begin(), pth.end()));				
		}
	}
	
	return List::create(Named("path") = paths.col(min_cost.index_min()), Named("cost") = min(min_cost));
}


// SA, the temperature parameters were obtained caballistically

//[[Rcpp::export]]
List tsp_sa(const arma::imat& m, const int max_iter=25000, const double t_amp=4000, const double t_center=0, const double t_width=3000)
{
	const int n = m.n_rows;
	const int n1 = n-1;
	ivec min_pth(n1), cost(2);
	imat pth(n1,2);

	int s1,s2, j=0;

	std::iota(pth.col(j).begin(), pth.col(j).end(),1);
		
	int min_cost =	m.at(0, pth.at(0,j)) + m.at(0, pth.at(n1-1,j));
	for(int i=1; i<n1; i++) 
		min_cost += m.at(pth.at(i-1,j),pth.at(i,j));
	
	cost[j] = min_cost;
	cost[1-j] = min_cost+1;
	min_pth = pth.col(j);
	
	for(int iter=0; iter<max_iter; iter++)
	{
		double temp = t_amp/(1+std::exp((iter-t_center)/t_width));
		pth.col(1-j) = pth.col(j);
		cost[1-j] = cost[j];

		sample2(s1,s2,n1);
		swp_pth(pth.colptr(j),s1,s2);

		cost[j] = m.at(0, pth.at(0,j)) + m.at(0, pth.at(n1-1,j));
		for(int i=1; i<n1; i++) 
			cost[j] += m.at(pth.at(i-1,j),pth.at(i,j));
		
		
		if(cost[j] < cost[1-j])
		{
			if(cost[j]<min_cost)
			{
				min_cost = cost[j];
				min_pth = pth.col(j);
			}			
		} else if(R::runif(0,1) > std::exp((cost[1-j] - cost[j]) / temp)) // _not_ accepted probabilistically
		{
			j = 1-j;
		} 		
	}
	
	ivec out(n);
	out[0] = 1;
	out.tail(n1) = min_pth+1;
	
	return List::create(Named("path") = out, Named("cost") = min_cost);
}


