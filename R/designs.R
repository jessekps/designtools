



#' Linear  test design through LP
#'
#' @param items data.frame with column item_id and any other columns used in the constraints or as 'friends'
#' @param ... constraints on test content or psychometric properties, see details
#' @param objective function to optimize, see details
#' @param pars item parameters, either a data.frame or a fit object from dexter::fit_enorm, dexterMST::fit_enorm_MST
#' @param friends name of a column in items that specifies friend groups (e.g. a text_id for reading texts)
#' @param population_density density function, eg `dnorm` or `function(x) dnorm(x,mean=1,1)`
#' @param nmod numer of modules
#' @param routing_module data.frame with column item_id, specifying the routing module
#'
#' @return A tibble with the selected items.
#'
#'
#' @examples
#' library(dplyr)
#'
#' set.seed(123)
#' pars = tibble(item_id=sprintf('item%03i',1:300), item_score=1,beta=runif(300,-2,2))
#' items = tibble(item_id=sprintf('item%03i',1:300), category = sample(c('A','B','C'),300,TRUE))
#' 
#' ds = test_design(items, 
#'                  item_count(category=='A') == 40,
#'                  item_count(category=='B') == 15,
#'                  item_count(category=='C') == 20,
#'                  difficulty_constraint(0.6, 0.7),
#'                  objective = maximize_information(theta=0),
#'                  pars=pars)
#' 
#' items2 = tibble(item_id = sprintf('item%03i',1:300),
#'                  category = sample(c('A','B','C'),300,TRUE),
#'                  text_id = sample(1:100,300,TRUE))
#' 
#' ds2 = test_design(items2, 
#'                   difficulty_constraint(0.6, 0.8),
#'                   item_fraction(category=='A') >= 0.5,
#'                   item_fraction(category=='A') <= 0.55,
#'                   item_fraction(category=='B') >= 0.15,
#'                   item_fraction(category=='B') <= 0.20,
#'                   item_fraction(category=='C') >= 0.25,
#'                   item_fraction(category=='C') <= 0.30,
#'                   irt_constraint(theta=c(0,1), max_sem=0.5),
#'                   pars=pars,
#'                   friends='text_id')
#' 
#' ds2 %>% count(category) %>% mutate(p=n/sum(n))
#' 
test_design = function(items, ...,
                        objective = random_test(),
                        pars = NULL,
                        population_density = dnorm,
                        friends = NULL)
{
  env = caller_env()
  constraints = list(...)
  items = arrange(items, .data$item_id)
  
  ### sanity checks
  if(!is.null(names(constraints)))
    stop("constraints should not be named")
  
  if(n_distinct(items$item_id) != nrow(items))
    stop('item_id in items is not unique')
  
  if(is.null(pars))
  {
    if(needs_pars(objective))
      stop('parameters must be specified with this objective function')
    
    for(ctr in constraints)
      if(needs_pars(ctr))
        stop(sprintf('parameters must be specified to use %s', class(ctr)),call.=FALSE)
    # need some sort of repr or as.character method for constraints
  } else
  {
    if(!inherits(pars,'data.frame'))
      pars = coef(pars)
    pars = semi_join(pars,items,by='item_id') %>%
      arrange(.data$item_id)
    if(length(setdiff(items$item_id,pars$item_id)) > 0)
      stop('parameters not found for some items')
  }
  if(!inherits(objective,'objective'))
    stop('Objective not valid, see ?objective_items')
  
  
  ### objective function
  objective = to_lp(objective, items=items,env=env,pars=pars)
  lpm = LPM(objective$xt, objective$direction)
  
  
  ### constraints
  for(ctr in constraints)
  {
    ctr_lp = to_lp(ctr, items=items, env=env, population_density=population_density, pars=pars)
    
    for(ct in ctr_lp)
      lpm$add_constraint(ct$xt, ct$type, ct$rhs)
  }
  
  ### friends
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
mst_design = function(items, routing_module,
                       ...,
                       objective = random_test(),
                       nmod=3L,
                       pars = NULL,
                       population_density = dnorm,
                       friends = NULL)
{
  env = caller_env()
  constraints = list(...)
  items = arrange(items, .data$item_id)
  
  ### sanity checks
  if(!is.null(names(constraints)))
    stop("constraints should not be named")
  
  if(n_distinct(items$item_id) != nrow(items))
    stop('item_id in items is not unique')
  
  if(is.null(pars))
  {
    if(needs_pars(objective))
      stop('parameters must be specified with this objective function')
    
    for(ctr in constraints)
      if(needs_pars(ctr))
        stop(sprintf('parameters must be specified to use %s', class(ctr)),call.=FALSE)
    # need some sort of repr or as.character method for constraints
  } else
  {
    if(!inherits(pars,'data.frame'))
      pars = coef(pars)
    if(length(setdiff(items$item_id,pars$item_id)) > 0)
      stop('parameters not found for some items')
  }
  if(!inherits(objective,'objective'))
    stop('Objective not valid, see ?objective_items')
  
  
  ## mst sanity
  nmod = as.integer(nmod)
  if(nmod<2) stop("nmod must be 2 or more")
  
  if(length(setdiff(routing_module$item_id, items$item_id)) > 0)
    stop('One or more items in the routing module are not found in "items"')
  
  routing_module = select(routing_module,.data$item_id) %>% 
    inner_join(items, by='item_id') %>%
    arrange(.data$item_id)
  
  items = anti_join(items, routing_module, by='item_id') %>%
    arrange(.data$item_id)
  
  if(nmod > (nrow(routing_module)+1)/2)
    stop('the routing module is too short to support this many modules in stage 2')
  
  if(nmod > nrow(items))
    stop('not enough items to support this many modules')
  
  
  # choose routing points to divide population in equal groups
  rt_pars = semi_join(pars, routing_module, by='item_id')
  s = seq(-6,6,0.1)
  w = population_density(s)
  w = w/sum(w)
  P = w %*% Ps_test(rt_pars,s)
  rt_cut = cut(cumsum(P),nmod,labels=FALSE)
  
  mst = list(rt_pars=rt_pars, rt_cut=rt_cut,items=routing_module)
  
  pars = semi_join(pars, items,by='item_id') %>%
    arrange(.data$item_id, .data$item_score)
  
  # objective function
  objective = to_lp(objective, items=items,env=env,pars=pars, mst=mst)
  lpm = LPM(objective$xt, objective$direction,nbooklets=nmod,nitems = nrow(items))
  
  ### constraints
  for(ctr in constraints)
  {
    ctr_lp = to_lp(ctr, items=items, env=env, population_density=population_density, pars=pars,mst=mst)
    
    for(ct in ctr_lp)
      do.call(lpm$add_constraint,ct)
  }
  
  
  if(!is.null(friends))
  {
    if(!friends %in% colnames(items))
      stop('Column',friends,'not found in item_properties')
    
    res = lpm$solve_separate(group = items[[friends]])

    if(!res$success)
    {
      message('no solution found')
      return(NULL)
    }
    res =res$result
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
    res = lpm$solve_separate()
    if(!res$success)
    {
      message('no solution found')
      return(NULL)
    }

    res = tibble(item_id=rep(items$item_id,nmod),
                 module = rep(1:nmod,each=nrow(items)),
                 stage=2L,booklet=.data$module)[res$result,]
  }
  #delete.lp(l)
  routing_module$module=0L
  routing_module$stage=1L
  routing_module = inner_join(routing_module, tibble(module=0L,booklet=1:nmod),by='module')
  res = inner_join(res,items,by='item_id')
  res = union_all(routing_module,res) %>%
    arrange(.data$booklet,.data$module,.data$item_id)
  
  
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


report_mst = function(design, routing,pars,population_density,theta=NULL)
{
  nmod = max(design$module)
  cols = setdiff(colnames(design),c('item_id','module','booklet','stage'))
  smr = function(items, grp)
  {
    if(grp == 'module')
      items = filter(items,.data$module != 0 | .data$booklet == 1)
    res = lapply(cols, function(cn){
      if(is.numeric(items[[cn]]))
      {
        items %>%
          group_by(.data[[grp]]) %>%
          summarise(mean=mean(.data[[cn]]))
      } else
      {
        items %>%
          add_count(.data[[grp]],name='nn') %>%
          group_by(.data[[grp]], .data[[cn]]) %>%
          summarise(n = n(), fraction = .data$n/.data$nn[1])
      }
    })
    names(res) = cols
    res
  }
  out = list(modules=smr(design,'module'), booklets = smr(design,'booklet'))
  
  items = filter(design, .data$module>0) %>%
    distinct(.data$item_id,.keep_all=TRUE) %>%
    arrange(.data$item_id)
  
  rt_items = filter(design, .data$module==0 & .data$booklet==1) %>%
    arrange(.data$item_id)
  
  if(!inherits(pars,'data.frame')) pars=coef(pars)
  
  rt_pars = semi_join(pars, rt_items, by='item_id') %>% arrange(.data$item_id)
  pars = semi_join(pars, items, by='item_id') %>% arrange(.data$item_id)
  mst = list(rt_items = rt_items, rt_cut=routing$next_module, rt_pars=rt_pars)
  
  included = tibble(item_id=rep(items$item_id,3), module=rep(1:nmod,each=nrow(items))) %>%
    left_join(design,by=c('module','item_id')) %>%
    mutate(included = !is.na(.data$booklet)) 
  
  if(!is.null(theta))
  {
    ii = sapply(theta,function(tht)
    {
      sum(to_lp(maximize_information(theta=tht),items=items,pars=pars,mst=mst)$xt[included$included]) +
        sum(to_lp(maximize_information(theta=tht),items=mst$rt_items,pars=mst$rt_pars)$xt)
    })
    out$information = tibble(theta=theta, information=ii, sem=1/sqrt(ii))
  }
  
  df = to_lp(difficulty_constraint(0.5), items=items, pars=pars, population_density=population_density,mst=mst)
  pval = bind_rows(lapply(1:nmod, function(i) tibble(module=i,item_id=items$item_id,pvalue=df[[i]]$xt+0.5))) 

  
  pval = included |>
    filter(included) |>
    select(item_id,module) |>
    mutate(booklet_id=module) |>
    inner_join(pval,by=c('item_id','module'))
  
  df = to_lp(difficulty_constraint(0.5), items=rt_items, pars=rt_pars, population_density=population_density,mst=mst)
  pval2 = bind_rows(lapply(1:nmod, function(i) tibble(booklet_id=i,module=0,item_id=rt_items$item_id,pvalue=df[[i]]$xt+0.5))) 
  out$pvalue = union_all(pval,pval2)
  out
}

report_linear = function(items,pars,population_density,theta=NULL)
{
  cols = setdiff(colnames(items),c('item_id','module','booklet','stage'))
  if(!inherits(pars,'data.frame'))
    pars=coef(pars)

  res = lapply(cols, function(cn)
    {
      if(is.numeric(items[[cn]]))
      {
        items %>%
          summarise(mean=mean(.data[[cn]]))
      } else
      {
        items %>%
          group_by(.data[[cn]]) %>%
          summarise(n = n(), fraction = n()/nrow(items))
      }
    })
  names(res) = cols
  out = list(content = res)

  items = arrange(items,.data$item_id)
  pars = semi_join(pars,items,by='item_id') %>% arrange(.data$item_id,.data$item_score)

  if(!is.null(theta))
  {
    ii = sapply(theta,function(tht) sum(to_lp(maximize_information(theta=tht),items=items,pars=pars)$xt)) 
    out$information = tibble(theta=theta, information=ii, sem=1/sqrt(ii))
  }
  
  xt = to_lp(difficulty_constraint(0.5), items=items, pars=pars, population_density=population_density)[[1]]$xt + 0.5

  out$pvalue = tibble(item_id=items$item_id,pvalue=xt)
  out
}


