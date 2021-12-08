
check_constraint = function(indices, cn, values, constraints, tolerance)
{
  cs = constraints[[cn]]
  tol = tolerance[[cn]]
  nit = sum(values$n_items[indices])
  if(cn == 'n_items')
  {
    between(nit, cs[1]-tol, cs[2])
  } else
  {
    m = cs * nit
    between(sum(values[[cn]][indices]), floor(m-tol), ceiling(m+tol))
  }
}


increase_tolerance = function(cn, values, constraints, tolerance)
{
  if(cn == 'n_items') stop("impossible")
  small = 1e-8
  
  cs = constraints[[cn]]
  tol = tolerance[[cn]]
  
  if(cn == 'n_items')
    return(tol+1L)
  
  nit = constraints$n_items[1]:constraints$n_items[2]
  tol + min((nit * cs - tol) %% 1, 1 - (nit * cs + tol) %% 1 ) + small
}

block_design = function(values, constraints, relax=TRUE)
{
  nblocks = nrow(values)
  cc = list()
  tolerance = list()
  for(s in seq_along(constraints))
  {
    cn = names(constraints)[s]
   
    m = matrix(FALSE,nblocks,nblocks)
    tolerance[[cn]] = 0L

    repeat
    {
      m[]=FALSE
      for(i in 1:(nblocks-1)) 
        for(j in (i+1):nblocks)
          m[i,j] = check_constraint(c(i,j), cn, values, constraints, tolerance)

      m[lower.tri(m)] = t(m)[lower.tri(m)]
      if(!relax || !is.null(full_cycle(m)))
        break
      tolerance[[cn]] = increase_tolerance(cn, values, constraints, tolerance)
    }
    cc[[cn]] = m
  }
  
  f = Reduce('&',cc )
  cl = full_cycle(f)
  
  # hole in one
  if(!is.null(cl) || !relax)
    return(list(solution = cl, tolerance = tolerance, ccm = cc))
  
  # the irreducible set of constraints is guaranteed larger than one
  # I would like to find the smallest one
  # this is an NP complete problem, have to find all combinations of constraints
  IIS = list()
  for(i in 2:length(constraints))
  {
    apply(combn(2:length(constraints),i),2, function(s)
    {
      if(!is.null(full_cycle(Reduce('&',cc[-s]))))
        IIS <<- append(IIS, list(names(constraints)[s]))
    })
    if(length(IIS)>0)
      break
  }
  IIS = tibble(cn = unlist(IIS[[1]])) %>%
    count(.data$cn) %>%
    mutate(tol = unlist(tolerance[.data$cn])) %>%
    arrange(.data$tol,desc(.data$n)) %>%
    pull(cn)

  repeat
  {
    for(cn in IIS)
    {
      tolerance[[cn]] = increase_tolerance(cn, values, constraints, tolerance)
      m[] = FALSE
      for(i in 1:(nblocks-1)) 
        for(j in (i+1):nblocks)
          m[i,j] = check_constraint(c(i,j), cn, values, constraints, tolerance)
      
      m[lower.tri(m)] = m[upper.tri(m)]
      cc[[cn]] = m 
      f = Reduce('&',cc )
      cl = full_cycle(f)
      if(!is.null(cl)) break
    }
    if(!is.null(cl)) break
  }
  list(solution = cl, tolerance = tolerance, ccm = cc)
}

full_cycle = function(m)
{
  if(any(rowSums(m)<2))
    return(NULL)
  
  N = nrow(m)
  nds = 1:N
  
  fcc = function(seen)
  {
    if(length(seen)==N)
    {
      if(m[seen[1],seen[N]])
        return(seen)
    } else
      for(i in nds)
        if(! i %in% seen)
          if(m[i,seen[1]])
          {
            a = fcc(c(i,seen))
            if(!is.null(a))
              return(a)
          }
    NULL
  }
  fcc(1)
}