#include <RcppArmadillo.h>
#include "shared.h"

using namespace arma;
using Rcpp::List;
using Rcpp::Named;



int bin_loss(const ivec& bin_prop, const ivec& bin_min, const ivec& bin_max)
{
	return accu((bin_prop < bin_min) % (bin_min - bin_prop)) + \
			accu((bin_prop > bin_max) % (bin_prop - bin_max));
}

// the simplest possible greedy matching with a tie breaker;

//[[Rcpp::export]]
List bdivide(const arma::imat& properties, const arma::ivec& bin_min, const arma::ivec& bin_max, const int nbins)
{
	// number of items should be larger than number of bins
	
	const int nit = properties.n_cols, m = properties.n_rows;
	ivec test(m);
	imat cur(m, nbins, fill::zeros);
	ivec ibin(nit, fill::value(-1));
	
	// absolute match if possible
	for(int i=0; i<nit; i++)
	{
		for(int j=0; j<nbins; j++) 
		{
			test = cur.col(j) + properties.col(i);
			if(all(test <= bin_max) && all(test >= bin_min - bin_max[0] + test[0])) // like betting an inside strait
			{
				ibin[i] = j;
				cur.col(j) = test;				
				break;
			}
		}
	}
	
	arma::ivec bin_test(nbins, fill::zeros);
	for(int j=0;j<nbins;j++)
		bin_test[j] = all(cur.col(j) >= bin_min);
	
	// assign somewhere least damaging
	for(int i=0; i< nit; i++) if(ibin[i] < 0)
	{
		int min_loss = 10000;
		int jbest = -1;
		for(int j=0; j<nbins; j++) if(bin_test[j] == 0)
		{
			test = cur.col(j) + properties.col(i);
			int loss = accu((test > bin_max) % (test - bin_max)) + accu((test < bin_min) % (bin_min - test)) + 10 * (test[0] > bin_max[0]);
			loss -= accu((cur.col(j) < bin_min) % (bin_min - cur.col(j)));
			if(loss < min_loss)
			{
				min_loss = loss;
				jbest = j;
			}
		}
		if(jbest >= 0)
		{
			ibin[i] = jbest;	
			cur.col(jbest) += properties.col(i);
		}
	}
	
	return List::create(Named("item_bin") = ibin+1, Named("bin_properties") = cur, Named("bin_test") = bin_test);
}


int bin_items(const arma::imat& properties, const arma::ivec& perm, const arma::ivec& bin_min, const arma::ivec& bin_max, const int nbins)
{
	// number of items should be larger than number of bins
	
	const int nit = properties.n_cols, m = properties.n_rows;
	ivec test(m);
	imat cur(m, nbins, fill::zeros);
	int missed=nit;
	
	// absolute match if possible
	for(int i=0; i<nit; i++)
	{
		for(int j=0; j<nbins; j++) 
		{
			test = cur.col(j) + properties.col(perm[i]);
			if(all(test <= bin_max) && all(test >= bin_min - bin_max[0] + test[0])) // like betting an inside strait
			{
				missed--;
				cur.col(j) = test;				
				break;
			}
		}
	}
	
	int res = 0;
	for(int j=0;j<nbins;j++)
		res += all(cur.col(j) >= bin_min);
	
	return missed - res;
}

//number of violations, all items allocated
int bin_items_nv(const arma::imat& properties, const arma::ivec& perm, const arma::ivec& bin_min, const arma::ivec& bin_max, const int nbins)
{
	// number of items should be larger than number of bins
	
	const int nit = properties.n_cols, m = properties.n_rows;
	ivec test(m);
	imat cur(m, nbins, fill::zeros);
	int missed=nit;
	
	// absolute match if possible
	for(int i=0; i<nit; i++)
	{
		for(int j=0; j<nbins; j++) 
		{
			test = cur.col(j) + properties.col(perm[i]);
			if(all(test <= bin_max) && all(test >= bin_min - bin_max[0] + test[0])) // like betting an inside strait
			{
				missed--;
				cur.col(j) = test;				
				break;
			}
		}
	}
	
	int res = 0;
	for(int j=0;j<nbins;j++)
		res += all(cur.col(j) >= bin_min);
	
	return missed - res;
}



// works like a charm, can also try number of violations as cost
//[[Rcpp::export]]
arma::ivec perm_bdivide(const arma::imat& properties, const arma::ivec& bin_min, const arma::ivec& bin_max, const int nbins,
				const int max_iter=25000, const double t_amp=4000, const double t_center=0, const double t_width=3000)
{
	const int nit = properties.n_cols;
	int s1,s2,j=0;
	imat perm(nit,2);
	std::iota(perm.col(j).begin(), perm.col(j).end(), 0);
	ivec min_perm = perm.col(j), cost(2);
	int min_cost = bin_items(properties, perm.col(j), bin_min, bin_max, nbins);
	cost[j] = min_cost;
	cost[1-j] = min_cost+1;
	
	for(int iter=0; iter<max_iter; iter++)
	{
		double temp = t_amp/(1+std::exp((iter-t_center)/t_width));
		perm.col(1-j) = perm.col(j);
		cost[1-j] = cost[j];

		sample2(s1,s2,nit);
		swp_pth(perm.colptr(j),s1,s2);
		//swp(perm[s1],perm[s2]);
		cost[j] = bin_items(properties, perm.col(j), bin_min, bin_max, nbins);

		if(cost[j] <= cost[1-j])
		{
			if(cost[j]<min_cost)
			{
				min_cost = cost[j];
				min_perm = perm.col(j);
			}			
		} else if(R::runif(0,1) > std::exp((cost[1-j] - cost[j]) / temp)) // _not_ accepted probabilistically
		{
			j = 1-j;
		} 		
	}
	
	return min_perm+1;
}



double bin_items2(const arma::ivec& nit, const arma::mat& properties, const arma::ivec& perm, const arma::vec& bin_opt, const int nbins)
{
	// number of items should be larger than number of bins
	
	const int n = properties.n_cols, m = properties.n_rows;
	mat cur(m,nbins,fill::zeros);
	ivec nit_bin(nbins, fill::zeros);
	
	for(int ii=0; ii<n; ii++)
	{
		int b = nit_bin.index_min(), i=perm[ii];
		nit_bin[b] += nit[i];
		cur.col(b) += properties.col(i);
	}
	
	return accu(abs(cur.each_col() - bin_opt));
}

//[[Rcpp::export]]
arma::ivec perm_bdivide2(const arma::ivec& n_nit, const arma::mat& properties, const arma::vec& bin_opt, const int nbins, 
						const int max_iter=25000, const double t_amp=4000, const double t_center=0, const double t_width=3000)
{
	const int nit = properties.n_cols;
	int s1,s2,j=0;
	imat perm(nit,2);
	std::iota(perm.col(j).begin(), perm.col(j).end(), 0);
	ivec min_perm = perm.col(j);
	vec cost(2);
	double min_cost = bin_items2(n_nit, properties, perm.col(j), bin_opt, nbins);
	cost[j] = min_cost;
	cost[1-j] = min_cost+1;
	
	for(int iter=0; iter<max_iter; iter++)
	{
		double temp = t_amp/(1+std::exp((iter-t_center)/t_width));
		perm.col(1-j) = perm.col(j);
		cost[1-j] = cost[j];

		sample2(s1,s2,nit);
		swp_pth(perm.colptr(j),s1,s2);
		//swp(perm[s1],perm[s2]);
		cost[j] = bin_items2(n_nit, properties, perm.col(j), bin_opt, nbins);

		if(cost[j] <= cost[1-j])
		{
			if(cost[j]<min_cost)
			{
				min_cost = cost[j];
				min_perm = perm.col(j);
			}			
		} else if(R::runif(0,1) > std::exp((cost[1-j] - cost[j]) / temp)) // _not_ accepted probabilistically
		{
			j = 1-j;
		} 		
	}
	
	return min_perm+1;
}

//[[Rcpp::export]]
List bdivide2(const arma::ivec& nit, const arma::mat& properties,  const arma::vec& bin_opt, const int nbins)
{
	// number of items should be larger than number of bins
	
	const int n = properties.n_cols, m = properties.n_rows;
	mat cur(m,nbins,fill::zeros);
	ivec nit_bin(nbins, fill::zeros), ibin(n, fill::value(-1));
	
	for(int i=0; i<n; i++)
	{
		int b = nit_bin.index_min();
		nit_bin[b] += nit[i];
		cur.col(b) += properties.col(i);
		ibin[i] = b;
	}
	
	return List::create(Named("item_bin") = ibin+1, Named("bin_properties") = cur, Named("cost") = accu(abs(cur.each_col() - bin_opt)));
}


