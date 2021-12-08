
# bizar genoeg was de oude beter?????

pretest_design2 = function(items, max_nit, balance=NULL, friends=NULL)
{
  overlap = 2L
  
  items = arrange(items, .data$item_id)
  items = mutate_if(items, is.numeric,as.character)
  nit = nrow(items)
  nblocks = ceiling(2*nit/max_nit)
  
  items = bin_items(items,nblocks,balance,friends=friends)$items %>%
    rename(block='bin')
  
  bd = count(items,.data$block,name='n_items')
  
  constraints = list(n_items = c(2*min(bd$n_items), max_nit))

  for(cn in balance)
  {
    bd = items %>%
      group_by_at(c('block',cn)) %>%
      tally() %>%
      pivot_wider(id_cols='block', names_from=cn, values_from='n', 
                  values_fill=0L, names_prefix = paste0(cn,'_')) %>%
      inner_join(bd, by='block')
    
    p = items[[cn]]
    for(catg in unique(p))
    {
      constraints[[paste(cn, catg, sep='_')]] = mean(p == catg)
    }
  }
  
  res = block_design(bd,constraints)
  
  bd$block = as.integer(bd$block)
  items$block = as.integer(items$block)
  r = rep(res$solution, each=2)
  r = c(r[-1],r[1])
  
  tibble( booklet_id = rep(1:nblocks, each=2),
          block = bd$block[r],
          bn = c(1L,rep(2:nblocks,each=2),1L)
  ) %>%
    inner_join(items,by='block') %>%
    select(-.data$block) %>%
    rename(block='bn')

}


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
      r= tsp_sa(m)
      if(r$cost < res$cost) res = r
    }
  }
  res$path = drop(res$path)
  res
}


pretest_design3 = function(items, max_nit, balance=NULL, friends=NULL)
{
  overlap = 2L
  
  items = arrange(items, .data$item_id)
  items = mutate_if(items, is.numeric,as.character)
  nit = nrow(items)
  nblocks = ceiling(2*nit/max_nit)
  
  res = bin_items2(items,nblocks,balance,friends=friends)
  
  items = rename(res$items,block='bin') 

  # now make a cost matrix for block association
  # optimal
  bk_opt = rowMeans(res$bin_properties) * 2
  m = matrix(0L,nblocks,nblocks)
  for(i in 1:(nblocks-1))
  {
    for(j in (i+1):nblocks)
    {
      m[i,j] = m[i,j] + sum(abs(rowSums(res$bin_properties[,c(i,j)]) - bk_opt))
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
  r = c(r[-1],r[1])
  
  tibble( booklet_id = rep(1:nblocks, each=2),
          block = (1:nblocks)[r],
          bn = c(1L,rep(2:nblocks,each=2),1L)
  ) %>%
    inner_join(items,by='block') %>%
    select(-.data$block) %>%
    rename(block='bn')
  
}

