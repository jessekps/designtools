
#' Simple pretest designs
#'
#' @param items data.frame with columns item_id and optionally weight and other columns. 
#' See details 
#' @param max_nit maximum number of items in a booklet; if items contains a column
#' `weight` then this is the maximum sum of item weight (e.g. expected time spent on an item)
#' @param balance categorical variables in items that will be balanced across booklets
#' @param friends name of an id column in items that specifies friend groups (e.g. reading texts)
#' 
#' @return data.frame
#' 
#' @details 
#' This function returns a simple chained block design which minimizes the number of booklet versions
#' 
#' to do: items$weight does not yet do anything
#'
pretest_design = function(items, max_nit, balance=NULL, friends=NULL)
{
  overlap = 2L
  
  items = arrange(items,.data$item_id)
  items = mutate_if(items,is.numeric,as.character)
  nit = nrow(items)

  nblocks = ceiling(2*nit/max_nit)
  
  nit_block = rep(ceiling(nit/nblocks),nblocks)

  if(sum(nit_block)>nit)
  {
    e = sum(nit_block)-nit
    nit_block[1:e] = nit_block[1:e]-1L
    nit_block = rev(nit_block)
  }
  
  
  if(is.null(friends))
  {
    items$item_quality = runif(nrow(items))
  } else
  {
    items = items %>%
      group_by(.data[[friends]]) %>%
      mutate(item_quality = n()) %>%
      ungroup()
  }

  args = list(pars=tibble(item_id=items$item_id,item_score=1L,beta=0), 
              optimize='quality',item_properties=items,
              friends=friends)

  adapt_args = function(nit_block, args, tol=0, tol_nit=NULL)
  {
    args$nit = nit_block
    if(!is.null(tol_nit))
      args$nit = args$nit + tol_nit 
    nit_block = max(args$nit)
    for(cn in balance)
    {
      b = list()
      p = args$item_properties[[cn]]
      for(catg in unique(p))
      {
        w = as.integer(p == catg)
        m = mean(w) * nit_block
        # let's say minimal tolerance
        b[[catg]] = c(floor(m-tol),ceiling(m+tol))
      }
      args[[cn]] = b
    }

    args
  }
  out=list()
  
  for(i in 1:(nblocks-1))
  {
    bd = do.call(test_design,adapt_args(nit_block[i], args))
    if(is.null(bd))
    {
      # very crude, we need an infeasibility indicator
      # next week
      tol=0
      # another stopgap measure
      tol_nit = 0
      if(!is.null(friends))
      {
        # see if the nit is possible
        tst_nit = test_design(pars=args$pars,optimize='length', nit=c(nit_block[i],1e8),
                              friends=friends, item_properties=args$item_properties)
        mx_nit = nrow(tst_nit)
        if(mx_nit > nit_block[i])
          tol_nit = c(-1,mx_nit-nit_block[i] )
      } else
      {
        tol = 0.2
      }

      while(NROW(bd)==0)
      {
        bd = do.call(test_design,adapt_args(nit_block[i], args,tol=tol,tol_nit=tol_nit))
        tol = tol +.1
      }
    }
    args$item_properties = anti_join(args$item_properties, bd, by='item_id')
    args$pars = anti_join(args$pars, bd, by='item_id')
    out[[i]] = bd
  }
  out[[nblocks]] = args$item_properties

  # nu samenvoegen
  d = bind_rows(out,.id='block') %>% select(-.data$item_quality)
  
  bd = count(d,.data$block,name='n_items')
  
  constraints = list(n_items = c(2L * min(bd$n_items), max_nit))
  
#  if(min(bd$n_items) < max(bd$n_items) && mean(bd$n_items == min(bd$n_items)) < 0.5)
#    constraints$n_items[1] = constraints$n_items[1]+1L
  
  for(cn in balance)
  {
    bd = d %>%
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
  d$block = as.integer(d$block)
  r = rep(res$solution, each=2)
  r = c(r[-1],r[1])
  
  tibble( booklet_id = rep(1:nblocks, each=2),
          block = bd$block[r],
          bn = c(1L,rep(2:nblocks,each=2),1L)
          ) %>%
      inner_join(d,by='block') %>%
    select(-.data$block) %>%
    rename(block='bn')
}



#' Linear and simple MST test design through LP
#'
#' 
#' @param routing_module data.frame with a column item_id giving the selection of items for the routing module
#' @param nmod number of modules in second stage
#' @param pars item parameters
#' @param item_properties data.frame with columns item_id, and one or more of
#' item_weight (for test length optimization, defaults to 1),
#' item_quality (for quality optimization),
#' other discrete properties (e.g. for test matrix constraints, item groupings)
#' @param optimize minimal test length (optionally weighted), maximum information (on `information_points`),
#' or maximum sum of `items$item_quality`
#' @param populaton_density density function for the population, the default is standard normal. 
#' the function should accept a vector of theta and return a vector of densities. 
#' To use for example a normal distribtution with a different mean, specify `population density = function(x) dnorm(x, mean=1)`
#' @param sem_constraint list/data.frame with columns theta and sem, maximum allowed SEM on these points
#' @param mean_pvalue vector of length two, minimum and maximum average expected pvalue.
#' Of course this depends on population_density
#' @param nit number of items, single integer or integer range
#' @param information_points only needed when `optimize = 'information'`, 
#' will optimize the mean information on these points
#' @param friends name of an id column in items that specifies friend groups (e.g. reading texts)
#' @param ... content constraints, the argument names must refer to columns in item_properties, see details
#'
#' @return test_design returns a tibble with the selected items. mst_design returns a list with elements design and routing
#'
#' @details
#' Content constraints refer to columns in item_properties that must contain discrete values. For example suppose
#' that `item_properties` contains a column named `bloom_tax` which contains the values 'know', 'apply' and 'evaluate'.
#' You can add a constraint specifying exactly how many of each should be in the test, e.g. 
#' `bloom_tax = list(know=10, apply=15, evaluate=5)`. You may also specify ranges, e.g.
#' `bloom_tax = list(know=c(10,15), apply=15, evaluate=(5,10))`. Lastly, you can specify these as ratio's. In this case
#' only ranges are allowed. For example: `bloom_tax = list(know=c(.2,.3), apply=c(.5,.6), evaluate=c(.2,.3))`.
#'  
#' MST designs are currently limited to two stages. The first stage must be supplied by the user. Routing is
#' optimized for equal size groups for the second stage.
#'
#'
#' @examples
#' library(dplyr)
#'
#' set.seed(123)
#' pars = tibble(item_id=sprintf('item%03i',1:300), item_score=1,beta=runif(300,-2,2))
#' items = tibble(item_id=sprintf('item%03i',1:300), category = sample(c('A','B','C'),300,TRUE))
#' 
#' ds = test_design(pars, item_properties=items,optimize='information', mean_pvalue = c(.6,.7),
#'                  category=list(A=40,B=15,C=20))
#'
#' items2 = tibble(item_id=sprintf('item%03i',1:300),
#'                 category = sample(c('A','B','C'),300,TRUE),
#'                 text_id = sample(1:100,300,TRUE))
#'
#' ds2 = test_design(pars, item_properties=items2, optimize='length', mean_pvalue = c(.6,.8),
#'                  category=list(A=c(0.5,0.55),B=c(0.15,0.20),C=c(0.25,0.30)),
#'                  sem_constraint = list(theta=c(0,1),sem=c(.5,.5)),
#'                  friends='text_id')
#'
#' ds2 %>% count(category) %>% mutate(p=n/sum(n))
#'
test_design = function(pars,
                        optimize = c('length','information','quality'),
                        item_properties=NULL,
                        population_density = dnorm,
                        sem_constraint = NULL,
                        mean_pvalue = NULL,
                        information_points = 0,
                        nit=NULL,
                        friends=NULL,
                        ...)
{
  # input sanitation
  optimize = match.arg(optimize)
  if(!inherits(pars,'data.frame'))
    pars = coef(pars)
  items = if(is.null(item_properties)) select(pars,.data$item_id) else item_properties
  items = arrange(items, .data$item_id)
  pars = arrange(pars, .data$item_id, .data$item_score)
  if(!setequal(items$item_id,pars$item_id))
    stop('item_properties and parameters must specify the exact same set of items')
  # content constraints
  cnt = list(...)
  
  if(optimize == 'length')
  {
    lpm = LPM(if('item_weight' %in% colnames(items)) items$item_weight else rep(1L,nrow(items)),'min')
  } else if(optimize == "information")
  {
    I=sapply(split(pars,pars$item_id), function(r) info(r$item_score, r$beta, information_points))
    if(!is.null(dim(I)))
      I = colSums(I)
    
    lpm = LPM(I,'max')
  } else
  {
    lpm = LPM(items$item_quality,'max')
  }
  
  
  if(!is.null(mean_pvalue))
  {
    # quadpoints will be slightly better but NC is adequate
    s = seq(-6,6,0.1)
    w = population_density(s)
    w = w/sum(w)
    pval = sapply(split(pars,pars$item_id), 
                  function(r){sum(w*Es(r$item_score, r$beta, s)/max(r$item_score))} )
    
    if(!all(pval>mean_pvalue[1]))
      lpm$add_constraint(pval-mean_pvalue[1], '>=', 0)
    
    if(!all(pval<mean_pvalue[2]))
      lpm$add_constraint(pval-mean_pvalue[2], '<=', 0)
  }
  if(!is.null(sem_constraint))
  {
    for(i in seq_along(sem_constraint$theta))
    {
      lpm$add_constraint(
        sapply(split(pars,pars$item_id), function(r) info(r$item_score, r$beta, sem_constraint$theta[i])),
      '>=', 1/sem_constraint$sem[i]**2)
    }
  }
    
  
  if(length(nit)==1)
  {
    lpm$add_constraint(1L, '=', as.integer(nit)) 
  } else if(length(nit)==2)
  {
    if(nit[1]>0)
      lpm$add_constraint(1L, '>=', as.integer(nit[1])) 
    
    if(nit[2] < nrow(items))
      lpm$add_constraint(1L, '<=', as.integer(nit[2])) 
  }
  # content categories
  for(cn in names(cnt))
  {
    categories = items[[cn]]
    if(is.null(categories))
      stop('Column',cn,'not found in item_properties')
    
    for(cat in names(cnt[[cn]]))
    {
      if(length(cnt[[cn]][[cat]])==1)
      {
        if(cnt[[cn]][[cat]]<1)
          stop("content category ratio's have to be specified as a range")
        
        lpm$add_constraint(as.integer(categories==cat), '=', as.integer(cnt[[cn]][[cat]])) 
      }
      else if(length(cnt[[cn]][[cat]])==2) #range
      {
        if(cnt[[cn]][[cat]][2] < cnt[[cn]][[cat]][1] )
          stop("content category contains an invalid range")
        
        if(cnt[[cn]][[cat]][1] < 1) #ratios
        {
          
          lpm$add_constraint(as.double(categories==cat) - cnt[[cn]][[cat]][1], '>=', 0) 
          lpm$add_constraint(as.double(categories==cat) - cnt[[cn]][[cat]][2], '<=', 0) 
        } else # absolute numbers
        {
          lpm$add_constraint(as.integer(categories==cat), '>=', cnt[[cn]][[cat]][1]) 
          lpm$add_constraint(as.integer(categories==cat), '<=', cnt[[cn]][[cat]][2]) 
        }
      }
    }
  }
  if(!is.null(friends))
  {
    if(!friends %in% colnames(items))
      stop('Column',friends,'not found in item_properties')
    
    l = lpm$get_model(group=items[[friends]])
    if(lpsolve(l))
    {
      res = as.logical(get.variables(l))
      g = sort(unique(items[[friends]]))[res]
      res = filter(items, .data[[friends]] %in% g)
    } else
    {
      return(NULL)
    }
    
  } else
  {
    l = lpm$get_model()
    if(lpsolve(l))
    {
      res = items[as.logical(get.variables(l)),]
    } else
    {
      return(NULL)
    }
    
  }
  #delete.lp(l)
  res
}

#' @rdname test_design
#' 
mst_design = function(pars,
                       routing_module,
                       nmod=3L,
                       optimize = c('length','information','quality'),
                       item_properties=NULL,
                       population_density = dnorm,
                       sem_constraint = NULL,
                       mean_pvalue = NULL,
                       information_points = 0,
                       nit=NULL,
                       friends=NULL,
                       ...
)
{
  # input sanitation
  optimize = match.arg(optimize)
  if(!inherits(pars,'data.frame'))
    pars = coef(pars)
  items = if(is.null(item_properties)) select(pars,.data$item_id) else item_properties
  items = arrange(items, .data$item_id)
  pars = arrange(pars, .data$item_id, .data$item_score)
  if(!setequal(items$item_id,pars$item_id))
    stop('item_properties and parameters must specify the exact same set of items')
  # content constraints
  cnt = list(...)
  
  
  rt_nit = nrow(routing_module)
  
  rt_pars = semi_join(pars, routing_module, by='item_id')
  
  # choose routing points to divide population in equal groups
  s = seq(-6,6,0.1)
  w = population_density(s)
  w = w/sum(w)
  P = w %*% Ps_test(rt_pars,s)
  rt_cut = cut(cumsum(P),nmod,labels=FALSE)
  
  # let's say we need minimal sem on a number of points and a test matrix
  # test matrix with % of (unknown nbr of) items is something we still have to do
  rt_items = semi_join(items,routing_module,by='item_id')
  items = anti_join(items,routing_module,by='item_id')
  pars = semi_join(pars,items, by='item_id')
  
  N = nrow(items)
  
  
  if(optimize == 'length')
  {
    lpm = LPM(if('item_weight' %in% colnames(items)) items$item_weight else 1L,'min', nbooklets = nmod, nitems=N)
  } else if(optimize == "information")
  {
    ps = Ps_test(rt_pars,information_points)
    if(!is.matrix(ps))
      ps=matrix(ps,nrow=1)
    
    p_mod = t(apply(ps,1,rowsum,group=rt_cut))
    I = sapply(split(pars,pars$item_id), function(r) info(r$item_score, r$beta, information_points))
    if(!is.matrix(I))
      I = matrix(I,nrow=1)
    
    I = colSums(t(sapply(1:nrow(I), function(i){
      as.double(sapply(p_mod[i,], function(p){ I[i,] * p}))
    })))
    
    lpm = LPM(I, 'max', nbooklets=nmod,nitems=N)
  } else
  {
    lpm = LPM(items$item_quality,'max',nbooklets=nmod,nitems=N)
  }
  
  if(!is.null(sem_constraint))
  {
    ps = Ps_test(rt_pars,sem_constraint$theta)
    p_mod = t(apply(ps,1,rowsum,group=rt_cut))
    
    I_rt = sapply(sem_constraint$theta, function(tht){
      sum(sapply(split(rt_pars,rt_pars$item_id), function(r) info(r$item_score, r$beta, tht)))
    })
    
    I = sapply(split(pars,pars$item_id), function(r) info(r$item_score, r$beta, sem_constraint$theta))
    for(i in 1:nrow(I))
    {
      lpm$add_constraint(as.double(sapply(p_mod[i,], function(p){ I[i,] * p})),
                         '>=', 1/sem_constraint$sem[i]**2 - I_rt)
    }
  }
  # length constraints
  if(!is.null(nit))
  {
    for(i in 1:nmod)
    {
      if(length(nit)==1)
        lpm$add_constraint(1L,'=', nit - rt_nit, bk=i)
      else
      {
        lpm$add_constraint(1L,'>=', nit[1] - rt_nit, bk=i)
        lpm$add_constraint(1L,'<=', nit[2] - rt_nit, bk=i)
      }
    }
  }
  
  if(!is.null(mean_pvalue))
  {
    s = seq(-6,6,0.1)
    w = population_density(s)
    w = w/sum(w)
    ps = Ps_test(rt_pars,s)
    p_mod = t(apply(ps,1,rowsum,group=rt_cut))
    rt_scores = 0:(ncol(ps) - 1)
    
    for(i in 1:nmod)
    {
      wm = w * p_mod[,i]
      wm = wm/sum(wm)
      rts = rt_scores[rt_cut==i]
      
      pval = sapply(split(pars,pars$item_id), 
                    function(r){sum(wm*Es(r$item_score, r$beta, s)/max(r$item_score))} )
      
      
      rt_pval = sum(wm * rowSums((ps[,rts+1]/rowSums(ps[,rts+1])) %*% diag(rts)))
      
      lpm$add_constraint(pval-mean_pvalue[1], '>=', nrow(rt_items)*mean_pvalue[1] - rt_pval, bk=i)
      lpm$add_constraint(pval-mean_pvalue[2], '<=', nrow(rt_items)*mean_pvalue[2] - rt_pval, bk=i)
    }
  }
  
  # content categories
  for(cn in names(cnt))
  {
    categories = items[[cn]]
    if(is.null(categories))
      stop('Column',cn,'not found in item_properties')
    for(cat in names(cnt[[cn]]))
    {
      stf = sum(rt_items[[cn]]==cat)
      ctf = as.integer(categories==cat)
      if(length(cnt[[cn]][[cat]])==1)
      {
        if(cnt[[cn]][[cat]]<1)
          stop("content category ratio's have to be specified as a range")
        
        for(i in 1:nmod)
          lpm$add_constraint(ctf, '=', cnt[[cn]][[cat]] - stf,bk=i)
      }
      else if(length(cnt[[cn]][[cat]])==2) #range
      {
        if(cnt[[cn]][[cat]][2] < cnt[[cn]][[cat]][1] )
          stop("content category contains an invalid range")
        
        if(cnt[[cn]][[cat]][1] < 1) #ratios
        {
          for(i in 1:nmod)
          {
            lpm$add_constraint(ctf - cnt[[cn]][[cat]][1], '>=', cnt[[cn]][[cat]][1] * nrow(rt_items) - stf,bk=i)
            lpm$add_constraint(ctf - cnt[[cn]][[cat]][2], '<=', cnt[[cn]][[cat]][2] * nrow(rt_items) - stf,bk=i)
          }
        } else
        {
          for(i in 1:nmod)
          {
            lpm$add_constraint(ctf, '>=', cnt[[cn]][[cat]][1] - stf,bk=i)
            lpm$add_constraint(ctf, '<=', cnt[[cn]][[cat]][2] - stf,bk=i)
          }
        }
      }
    }
  }
  
  if(!is.null(friends))
  {
    if(!friends %in% colnames(items))
      stop('Column',friends,'not found in item_properties')
    
    l = lpm$get_model(group = items[[friends]])
    status = solve(l)
    if(!status %in% 0:1)
    {
      message('no solution found')
      return(NULL)
    }
    res = as.logical(get.variables(l))
    f = sort(unique(items[[friends]]))
    res = apply(matrix(res,byrow=TRUE,nrow=nmod),1, function(g)
    {
      tibble(item_id = items$item_id[items[[friends]] %in% f[g]],
             stage=2L)  
    }) %>%
      bind_rows(.id='module') %>%
      mutate(module=as.integer(.data$module))
    
  } else
  {
    l = lpm$get_model()
    status = solve(l)
    if(!status %in% 0:1)
    {
      message('no solution found')
      return(NULL)
    }
    res = tibble(item_id=rep(items$item_id,nmod),
                 module = rep(1:nmod,each=nrow(items)),
                 stage=2L,booklet=.data$module)[as.logical(get.variables(l)),]
  }
  #delete.lp(l)
  
  res = bind_rows(res, tibble(item_id = rep(rt_items$item_id,nmod),
                              module=0L,stage=1L,
                              booklet=rep(1:nmod,each=nrow(rt_items)))) 
  
  if(!is.null(item_properties))
    res = inner_join(res,item_properties,by='item_id')
  
  #dmst_rt = tibble(module_id='M0',exit_min)
  
  out = list(design=res,routing=tibble(module_score=1:length(rt_cut)-1L, next_module=rt_cut))
  
  # dexter mst formats
  rr = out$routing %>% 
    group_by(.data$next_module) %>% 
    summarise(exit_min=min(.data$module_score),exit_max=max(.data$module_score)) %>%
    rename(booklet_id=.data$next_module) %>%
    mutate(module_nbr=1L, module_id='M0') %>%
    bind_rows(tibble( module_nbr=2L,booklet_id=1:nmod,
                      module_id = paste0('M',1:nmod),
                      exit_min=NA_integer_, exit_max=NA_integer_)) %>%
    mutate(booklet_id=paste0('bk',.data$booklet_id))
  
  td = distinct(out$design, .data$module, .data$item_id) %>%
    mutate(module_id=paste0('M',.data$module)) %>%
    select(.data$module_id,.data$item_id)
  
  out$dexterMST = list(routing_rules=rr,test_design=td)
  out
}

