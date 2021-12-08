utils::globalVariables(c("."))

#' Create a grouping variable by size
#'
#' divide objects into bins of equal size. Uses a sort of greedy method
#'
#' @param size typically an integer (e.g. number of persons within schools)
#' @param k number of groups (groupsize is variable)
#' @param n maximum size of groups (number of groups is variable)
#'
#' @details
#' You must specify only one of k or n. Since bin packing and the multi-way partition problem are NP-hard,
#' the solution is not guaranteed to be optimal. 
#' 
#' Bin packing (i.e. maximum size groups), puts items of larger size together, 
#' the number of items in each bin thus tends to be very unequal
#'
size_groups = function(size, k=NULL, n=NULL)
{
  if(!xor(is.null(k),is.null(n)))
    stop("you must specify one of `k` or `n`")

  if(!is.null(k))
  {
    if(length(size)<k)
      stop('size may not be smaller than k')
    if(length(size)==k)
      return(1:length(size))
    
    prt = kmk_partition(size, k) %>% arrange(size)
    
    s = tibble(size=size, i=1:length(size)) %>% arrange(size)
    
    s$group = prt$group
    
    s %>% arrange(i) %>% pull(group)
    
  } else
  {
    if(n<max(size))
      stop("n is smaller than maximum size, grouping is not possible")
    # dit kan veel beter
    add = 0
    while(TRUE)
    {
      s = as.list(sort(size, decreasing=TRUE))
      names(s) = as.character(seq_along(s))

      k = ceiling(sum(size)/n) + add

      g = vector(mode='list',k)
      gn = integer(k)

      for(ns in names(s))
      {
        for(i in seq_along(g))
        {
          if(gn[i] + s[[ns]] <= n)
          {
            gn[i] = gn[i] + s[[ns]]
            s[[ns]] = NULL
            g[[i]] = c(g[[i]],ns)
            break
          }
        }
      }
      if(length(s)==0) break
      add=add+1
    }
    a = tibble(sz = size, indx=1:length(size)) %>%
      arrange(desc(.data$sz))

    b = tibble(i = as.integer(unlist(g)), group = rep(1:length(g), sapply(g,length))) %>%
      arrange(.data$i)

    a$group = b$group
    a %>%
      arrange(.data$indx) %>%
      pull(.data$group) %>%
      as.integer()
  }
}

divide_items = function(iprop, bmin, bmax, nbins)
{
  res = bdivide(iprop, bmin, bmax, as.integer(nbins))
  for(iter in 1:8)
  {
    perm = perm_bdivide(iprop, bmin, bmax, as.integer(nbins))
    res2 = bdivide(iprop[,perm], bmin, bmax, as.integer(nbins))
    if(sum(res2$bin_test) > sum(res$bin_test))
    {
      res = res2
      res$item_bin = res$item_bin[order(perm)]
    }
  }
  # possibly iterate to not have max items larger than desired
  res
  
  
  
}

# dit kan beter
# zorg dat min-nit + max_nit iig <= max_nit_bk
# de 0-1 verdeling lijkt me niet ok (zou alleen ok zijn als niet wordt geschonden, maar zo niet
# dan kunnen de 2 'grootste' items bij elkaar komen)
# dus wellicht -> dit proberen, als resultaat niet perfect dan
# bmin en bmax echte getallen maken iig waar het wordt geschonden
# we kunnen ook nog iets anders bedenken als we een echte cost/quality kunnen bedenken van een bin, maar indelen wordt
# dan lastiger
bin_items = function(items, nbins, balance=NULL, friends=NULL)
{
  
  if(!is.null(balance))
  {
    np = 1 + sum(sapply(balance, function(n) n_distinct(items[[n]])))
    nit = nrow(items)
    if(is.null(friends))
    {
      iprop = matrix(0L,np,nit)
      iprop[1,] = 1L
      rownames(iprop) = c('n', rep('',nrow(iprop)-1))
      j=2
      for(cn in balance)
      {
        for(v in unique(items[[cn]]))
        {
          rownames(iprop)[j] = paste0(cn,'_',v)
          iprop[j,] = as.integer(items[[cn]] == v)
          j = j+1
        }
      }
      bmin = as.integer(floor(rowSums(iprop)/nbins))
      bmax = as.integer(ceiling(rowSums(iprop)/nbins))
      
      res = divide_items(iprop, bmin, bmax, nbins)
      
      items$bin = res$item_bin
      
      list(items=items, bin_prop = res$bin_properties)
      
    } else
    {
      bd = count(items, .data[[friends]]) %>% arrange(desc(n))
      for(cn in balance)
      {
        bd = items %>%
          group_by_at(c(friends,cn)) %>%
          tally() %>%
          pivot_wider(id_cols=friends, names_from=cn, values_from='n', 
                      values_fill=0L, names_prefix = paste0(cn,'_')) %>%
          inner_join(bd,., by=friends)
      }
      iprop = t(as.matrix(bd[,-1]))
      
      bmin = integer(nrow(iprop))
      bmax = rep(1L,nrow(iprop))
      
      for(j in 1:nrow(iprop))
      {
        if(sum(iprop[j,]>0) <= nbins)
        {
          iprop[j,] = as.integer(iprop[j,] > 0)
        } else
        {
          sz = iprop[j,iprop[j,]>0]
          s = rowsum(sz,size_groups(sz,k=nbins))
          bmin[j] = min(s)
          bmax[j] = max(s)
        }
      }
      res = divide_items(iprop, bmin, bmax, nbins)

        
        # tibble(bin=1:nbins,
        #         too_few = apply(res$bin_properties,2, function(x) sum((bmin-x)*(x<bmin))),
        #         too_many = apply(res$bin_properties,2, function(x) sum((x-bmax)*(x>bmax))))
      
      bins = select(bd,.data[[friends]])
      bins$bin=res$item_bin
      
      list(items=inner_join(bins,items, by=friends), 
        bin_prop = res$bin_properties) 
    }
  } else
  {
    if(is.null(friends))
    {
      r = nrow(items) %% nbins
      m = rep(floor(nrow(items)/nbins), nbins)
      m[1:r] = m[1:r] + 1L
      items$bin = rep(1:nbins, m)
      list(items=items,bin_prop=matrix(m,nrow=1))
    } else
    {
      bd = items %>%
        count(.data[[friends]]) %>%
        mutate(bin = size_groups(n), k=nbins)
      
      list(items= inner_join(select(bd,-.data$n), items, by=friends), 
           bin_prop = rowsum(bd$n,bd$bin)) 
      
    }
  }
}

bin_items2 = function(items, nbins, balance=NULL, friends=NULL)
{
  
  if(is.null(balance) && is.null(friends))
  {
    r = nrow(items) %% nbins
    m = rep(floor(nrow(items)/nbins), nbins)
    m[1:r] = m[1:r] + 1L
    items$bin = rep(1:nbins, m)
    list(items=items,bin_prop=matrix(m,nrow=1))  
  }
  
  if(is.null(friends)) friends = 'item_id'
  
  np = 1 
  if(is.null(balance))
    np = np + sum(sapply(balance, function(n) n_distinct(items[[n]])))
  
  nit = nrow(items)
  bd = count(items, .data[[friends]]) %>% arrange(desc(n))
  
  for(cn in balance)
  {
    bd = items %>%
      group_by_at(c(friends,cn)) %>%
      tally() %>%
      pivot_wider(id_cols=friends, names_from=cn, values_from='n', 
                  values_fill=0L, names_prefix = paste0(cn,'_')) %>%
      inner_join(bd,., by=friends)
  }
  iprop = t(as.matrix(bd[,-1])) + 0
  
  bin_opt = rowSums(iprop) / nbins
  
  perm = 1:nbins
  res = bdivide2(bd$n[perm], iprop[,perm], bin_opt, as.integer(nbins))
  for(i in 1:8)
  {
    p = perm_bdivide2(bd$n, iprop, bin_opt, as.integer(nbins))
    r = bdivide2(bd$n[p], iprop[,p], bin_opt, as.integer(nbins))
    if(r$cost<res$cost)
    {
      res=r; perm=p
    }
  }
  
  
  bins = select(bd,.data[[friends]])
  bins$bin = res$item_bin[order(perm)]
  list(items=inner_join(items,bins,by=friends), bin_properties = res$bin_properties)
}

