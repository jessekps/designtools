

tsalesman = function(m)
{
  n=nrow(m)
  if(n<=3)
    return(list(path=1:m, cost=sum(m)/2))
  
  if(n<13)
  {
    res = tsp_complete2(m)  
  } else
  {
    res = tsp_sa(m)
    for(i in 1:5)
    {
      r = tsp_sa(m)
      if(r$cost < res$cost) res = r
    }
  }
  res$path = drop(res$path)
  res
}

#' Simple pretest designs
#'
#' @param items data.frame with columns item_id and variables listed under "balance" or "friends"
#' @param max_nit maximum number of items in a booklet.
#' @param balance categorical variables in items that will be balanced across booklets
#' @param friends name of an id column in items that specifies friend groups (e.g. reading texts)
#' 
#' @return data.frame
#' 
#' @details 
#' This function returns a simple chained block design which minimizes the number of booklet versions
#'
pretest_design = function(items, max_nit, balance=NULL, friends=NULL)
{
  overlap = 2L
  
  items = arrange(items, .data$item_id)
  items = mutate_if(items, is.numeric,as.character)
  nit = nrow(items)
  nblocks = ceiling(2*nit/max_nit)
  
  res = bin_items2(items,nblocks,balance,friends=friends)
  
  items = rename(res$items,block='bin') 

  if(!all(apply(res$bin_properties,1,n_distinct)==1))
  {
    # now make a cost matrix for block association
    # optimal
    bk_opt = rowMeans(res$bin_properties) * 2
    m = matrix(0L,nblocks,nblocks)
    for(i in 1:(nblocks-1))
    {
      for(j in (i+1):nblocks)
      {
        m[i,j] = m[i,j] + sum(abs(rowSums(res$bin_properties[,c(i,j),drop=FALSE]) - bk_opt))
        if(sum(res$bin_properties[1,c(i,j)]) > max_nit) m[i,j] = m[i,j] + 1000
      }
      
    }
    diag(m)=0
    m[lower.tri(m)] = t(m)[lower.tri(m)]
    
    #alleen nog maar voor integers geprogrammeerd
    
    m = m*1000
    mode(m)='integer'
    
    res2 = tsalesman(m)
    
  
    r = rep(res2$path, each=2)
  } else
  {
    r = rep(1:nblocks,each=2) - 1L
  }
  r = c(r[-1],r[1])
  tibble( booklet_id = rep(1:nblocks, each=2),
          block = (1:nblocks)[r+1L],
          bn = c(1L,rep(2:nblocks,each=2),1L)
  ) %>%
    inner_join(items,by='block') %>%
    select(-.data$block) %>%
    rename(block='bn')
  
}

