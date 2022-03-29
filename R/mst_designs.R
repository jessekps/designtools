

valid_pval = function(pval) all(pval >= 0) && all(pval <= 1) && length(pval)==2 && pval[1] < pval[2]

pval_mst_input = function(pval, nmod, name = deparse(substitute(pval)))
{
  if(is.list(pval))
  {
    if(length(pval != nmod)) 
      stop(sprintf('if %s is a list it should have length %i, the same as the number of modules', name,nmod),.call=FALSE)
  } else
  {
    pval = lapply(1:nmod, function(i) pval)
  }
  
  if(!all(sapply(pval, valid_pval)))
    stop('invalid pvalue range in path_mean_pvalue') 
  pval
}

#' Simple 2-stage MST test design through LP
#'
#' Determinse a design for the secondary modules and routing give a pretest an specifications.
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
#' @param population_density density function for the population, the default is standard normal. 
#' the function should accept a vector of theta and return a vector of densities. 
#' To use for example a normal distribtution with a different mean, specify `population density = function(x) dnorm(x, mean=1)`
#' @param sem_constraint list/data.frame with columns theta and sem, maximum allowed SEM on these points
#' @param path_mean_pvalue Either a vector of length two, minimum and maximum average expected pvalue, or a list of length \code{nmod} with minimum and maximum pvalue per path (including the touring module).
#' @param module_mean_pvalue Either a vector of length two, minimum and maximum average expected pvalue, or a list of length \code{nmod} with minimum and maximum pvalue per module in the second stage.
#' @param nit number of items, single integer or integer range
#' @param information_points only needed when `optimize = 'information'`, 
#' will optimize the mean information on these points
#' @param friends name of an id column in items that specifies friend groups (e.g. reading texts)
#' @param ... content constraints, the argument names must refer to columns in item_properties, see details
#'
#' @return A list with elements design and routing
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
mst_design = function(pars,
                      routing_module,
                      nmod=3L,
                      optimize = c('length','information','quality'),
                      item_properties=NULL,
                      population_density = dnorm,
                      sem_constraint = NULL,
                      path_mean_pvalue = NULL,
                      module_mean_pvalue = NULL,
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
  
  
  
  if(!is.null(path_mean_pvalue)) path_mean_pvalue = pval_mst_input(path_mean_pvalue, nmod)
  if(!is.null(module_mean_pvalue)) module_mean_pvalue = pval_mst_input(module_mean_pvalue, nmod)
  
  if(!is.null(path_mean_pvalue) || !is.null(module_mean_pvalue))
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
      
      if(!is.null(path_mean_pvalue))
      {
        lpm$add_constraint(pval-path_mean_pvalue[[i]][1], '>=', nrow(rt_items)*path_mean_pvalue[[i]][1] - rt_pval, bk=i)
        lpm$add_constraint(pval-path_mean_pvalue[[i]][2], '<=', nrow(rt_items)*path_mean_pvalue[[i]][2] - rt_pval, bk=i)
      }
      if(!is.null(module_mean_pvalue))
      {
        lpm$add_constraint(pval-module_mean_pvalue[[i]][1], '>=', nrow(rt_items)*module_mean_pvalue[[i]][1], bk=i)
        lpm$add_constraint(pval-module_mean_pvalue[[i]][2], '<=', nrow(rt_items)*module_mean_pvalue[[i]][2], bk=i)
      }
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

