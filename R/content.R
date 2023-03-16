
# to do: try https://github.com/dirkschumacher/rcbc as a back end, we could try to handle infeasibility
# (or manually try to do adjust the optimalization function but give it a try)

setClassUnion("call_or_NULL_or_name", c("call", "NULL","name"))
setClassUnion("integer_or_NULL", c("integer", "NULL"))

to_lp = function(x,...) UseMethod("to_lp")
needs_pars = function(f)
{
  m = getS3method('to_lp', class(f)[1])
  'pars' %in% names(formals(m))
}

# content constraints -----------------------------------------------------


content_spec = setClass("cspec",slots=c(qpredicate='call_or_NULL_or_name',qsubset='call_or_NULL_or_name', paths="integer_or_NULL",
                                 type='character'))

spec_paths = function(paths)
{
  if(is.null(paths)) return(NULL)
  if(is.integer(paths)) return(paths)
  if(is.double(paths) && all(paths %% 1 == 0)) return(as.integer(paths))
  stop("paths must be NULL or an integer vector")
}


#' Content specification
#' 
#' Specify content for test designs. These functions are meant to be used with \code{\link{test_design}} or \code{\link{mst_design}} 
#' 
#' @param predicate predicate
#' @param subset subset
#' @param paths paths to which the specification should apply. NULL means all paths. Paths are indicated by consecutive integers
#' where 1 stands for the path to which persons with the lowest score on the routing tests are routed. Argument is ignored for linear tests
#' 
#' 
#' 
item_fraction = function(predicate, subset=NULL, paths=NULL)
{
  qp = eval(substitute(quote(predicate)))
  qs = eval(substitute(quote(subset)))
  
  content_spec(qpredicate=qp,qsubset=qs, paths=spec_paths(paths),type='fraction')
}

#' @rdname item_fraction
item_sum = function(predicate,paths=NULL)
{
  qp = eval(substitute(quote(predicate)))

  content_spec(qpredicate=qp, qsubset=NULL, type='sum', paths=spec_paths(paths))
}

#' @rdname item_fraction
item_count = function(predicate=NULL,paths=NULL)
{
  qp = eval(substitute(quote(predicate)))
  if(is.null(qp))
    qp = quote(rep(TRUE,length(item_id)))
  
  content_spec(qpredicate=qp, qsubset=NULL, type='count', paths=spec_paths(paths))
}


cspec_compare = function(cspec, op = c('>','>=','<','<=','=='), e2)
{
  out = list(qpredicate = cspec@qpredicate, qsubset=cspec@qsubset, paths=cspec@paths, type=cspec@type,
       op = match.arg(op), e2=e2)
  class(out) = c('content_specification',class(out))
  out
}

cspec_compare2 = function(cspec, op = c('>','>=','<','<=','=='), cspec2)
{
  if(cspec@type != cspec2@type)
    stop('canot compare different types of content specification')
  
  out = list(e1=list(qpredicate = cspec@qpredicate, qsubset=cspec@qsubset, paths=cspec@paths, type=cspec@type,
                    op = match.arg(op),e2=0),
             e2=list(qpredicate = cspec2@qpredicate, qsubset=cspec2@qsubset, paths=cspec2@paths, type=cspec2@type,
                     op = '==',e2=0))
  
  class(out$e1) = c('content_specification',class(out$e1))
  class(out$e2) = c('content_specification',class(out$e2))
  class(out) = c('content_specification2',class(out))
  out
}


#' Comparison methods for content specifications
#' 
#' @param e1 content specification or a scalar number
#' @param e2 content specification or a scalar number
#' 
#' @rdname  content-specification-methods
setMethod('>=',c(e1="cspec",e2="numeric"), function(e1,e2) cspec_compare(e1,'>=',e2))
#' @rdname  content-specification-methods
setMethod('>', c(e1="cspec",e2="numeric"), function(e1,e2) cspec_compare(e1,'>',e2))
#' @rdname  content-specification-methods
setMethod('<=',c(e1="cspec",e2="numeric"), function(e1,e2) cspec_compare(e1,'<=',e2))
#' @rdname  content-specification-methods
setMethod('<', c(e1="cspec",e2="numeric"), function(e1,e2) cspec_compare(e1,'<',e2))
#' @rdname  content-specification-methods
setMethod('==', c(e1="cspec",e2="numeric"), function(e1,e2) cspec_compare(e1,'==',e2))
#' @rdname  content-specification-methods
setMethod('>=',c(e1="numeric",e2="cspec"), function(e1,e2) cspec_compare(e2,'<=',e1))
#' @rdname  content-specification-methods
setMethod('>', c(e1="numeric",e2="cspec"), function(e1,e2) cspec_compare(e2,'<',e1))
#' @rdname  content-specification-methods
setMethod('<=',c(e1="numeric",e2="cspec"), function(e1,e2) cspec_compare(e2,'>=',e1))
#' @rdname  content-specification-methods
setMethod('<', c(e1="numeric",e2="cspec"), function(e1,e2) cspec_compare(e2,'>',e1))
#' @rdname  content-specification-methods
setMethod('==', c(e1="numeric",e2="cspec"), function(e1,e2) cspec_compare(e2,'==',e1))

#' @rdname  content-specification-methods
setMethod('>=',c(e1="cspec",e2="cspec"), function(e1,e2) cspec_compare2(e2,'<=',e1))
#' @rdname  content-specification-methods
setMethod('>', c(e1="cspec",e2="cspec"), function(e1,e2) cspec_compare2(e2,'<',e1))
#' @rdname  content-specification-methods
setMethod('<=',c(e1="cspec",e2="cspec"), function(e1,e2) cspec_compare2(e2,'>=',e1))
#' @rdname  content-specification-methods
setMethod('<', c(e1="cspec",e2="cspec"), function(e1,e2) cspec_compare2(e2,'>',e1))
#' @rdname  content-specification-methods
setMethod('==', c(e1="cspec",e2="cspec"), function(e1,e2) cspec_compare2(e2,'==',e1))


to_lp.content_specification2 = function(f,...,items,env, mst=NULL)
{
  a = to_lp(f[[1]],...,items=items,env=env, mst=mst)
  b = to_lp(f[[2]],...,items=items,env=env, mst=mst)
  for(i in seq_along(a))
    a[[i]]$xt = a[[i]]$xt - b[[i]]$xt
  a
}


to_lp.content_specification = function(f,...,items,env, mst=NULL)
{

  p = eval_tidy(f$qpredicate, data=items,env=env)
  s = eval_tidy(f$qsubset, data=items,env=env)
  e2 = f$e2
  op = f$op
  
  if(length(e2)!=1) stop(paste0('item_',f$type,'can only be compared with an element of length 1'))
  if(f$type == 'fraction')
  {
    if(e2<0 || e2 > 1) stop('A fraction must be between 0 and 1')
  } 
  if(length(p) != nrow(items))
    stop('predicate must evaluate to a single value per row of items')
  
  if(f$type %in% c('fraction','count'))
  {
    if(!is.logical(p))
      stop('Predicate must evaluate to a logical vector')
    if(!is.null(s) && !is.logical(s))
      stop('Subset must evaluate to a logical vector')
    p = as.integer(p)
    s = as.integer(s)
  }
  
  if(any(is.na(p)) || any(is.na(s)))
    stop('Predicate or subset may not contain NA values')
  
  if(f$type == 'sum')
  {
    p = as.numeric(p)
    e2 = as.numeric(e2)
    #if(any(p<0))
    #  stop('predicate may not result in negative values')
  }
  
  if(f$type == 'fraction')
  {
    xt = as.double(p) - e2
    if(!is.null(s))
      xt[!s] = 0
    
    rhs=0
  } else
  {
    xt = p
    rhs = e2
  }
  
  if(op == '>=')
  {
    op = '>'
    rhs = rhs + 1e-15
  } else if(op == '<=')
  {
    op = '<'
    rhs = rhs - 1e-15
  }
  if(op == '==')
    op = '='
  
  
  if(!is.null(mst))
  {
    list(list(xt=xt,type=op,rhs=rhs - sum(to_lp(f,items=mst$items,env=env)[[1]]$xt), all_bk=TRUE))
  } else
  {
    list(list(xt=xt,type=op,rhs=rhs))
  }
}


# IRT constraints ---------------------------------------------------------

#' IRT constraints on test designs
#' 
#' Typically used when the length of the test is minimized. You must specify either min_information or max_sem but not both.
#' 
#' @param theta points on the ability scale to evaluate information
#' @param min_information the minimum information required
#' @param max_sem the maximum allowed Standard Error of Measurement
#' 
#' 
#' 
irt_constraint = function(theta, min_information=NULL, max_sem=NULL)
{
  if(is.null(min_information))
  {
    if(is.null(max_sem))
      stop('irt_constraint requires either "min_information" or "max_sem" to be specified')
    
    if(anyNA(max_sem))
      stop('NA values in max_sem')
    
    if(any(max_sem <= 0))
      stop('max_sem needs to be larger than 0')
    
    min_information = 1/max_sem**2
  } else
  {
    if(anyNA(min_information))
      stop('NA values in min_information')
    if(any(min_information <=0))
      stop('min_information needs to be larger than 0')
  }
  if(anyNA(theta))
    stop('NA values in theta')
  
  if(length(min_information)==1)
    min_information = rep(min_information, length(theta))
  if(length(theta)!=length(min_information))
    stop('vector theta needs to be of the same length as min_information or max_sem')
  
  out = list(theta=theta, min_information=min_information)
  class(out) = c('constraint_information',class(out))
  out
}

to_lp.constraint_information = function(f,...,items, pars,mst=NULL)
{
  if(!is.null(mst))
  {
    #deduct information gained from routing test
    f$min_information = f$min_information - sapply(f$theta, function(tht){
      sum(sapply(split(mst$rt_pars,mst$rt_pars$item_id), function(r) info(r$item_score, r$beta, tht)))
    })
  }
  
  mapply(f$theta,f$min_information,SIMPLIFY=FALSE, FUN=function(tht,i)
  {
    xt = to_lp(maximize_information(tht),items=items,pars=pars,mst=mst)$xt
    list(xt=xt,type='>',rhs=i)
  })
}




# classical constraints ---------------------------------------------------

# sanity check functions
null_or_in_range = function(val,mn=0,mx=1,name=deparse(substitute(val)))
{
  if(!is.null(val))
  {
    if(length(val) != 1)
      stop(sprintf('"%s" must have length 1',name))
    if(!between(val,mn,mx))
      stop(sprintf("%s needs to be between %f and %f",name,mn,mx))
  }
}
null_or_is_range = function(a,b,name_a=deparse(substitute(a)), name_b=deparse(substitute(b)))
{
  if(!is.null(a) && !is.null(b) && a>=b)
    stop(sprintf('"%s" and "%s" do not form a valid range',name_a,name_b))
}

#' Constraints on test or item difficulty
#' 
#' These contraints use the item parameters and population distribution to limit the difficulty of the test or individual items
#' to a specified range
#'
#' @param avg_pval_min minimum allowed average pvalue (proportion correct) on the test
#' @param avg_pval_max maximum allowed average pvalue on the test
#' @param pval_min minimum allowed pvalue for any item in the test
#' @param pval_max maximum allowed pvalue for any item in the test
#' @param modules integer vector of modules to apply to, with the 'easiest' module specified by 1. If NULL, applies to all modules
#' 
#' @details
#' You must specify at least one argument.
#'
#'
difficulty_constraint = function(avg_pval_min=NULL, avg_pval_max=NULL, pval_min=NULL, pval_max=NULL, modules=NULL)
{
  null_or_in_range(avg_pval_min)
  null_or_in_range(avg_pval_max)
  null_or_in_range(pval_min)
  null_or_in_range(pval_max)
  null_or_is_range(avg_pval_min,avg_pval_max)
  null_or_is_range(pval_min,pval_max)

  out = list(avg_pval_min=avg_pval_min, avg_pval_max=avg_pval_max, pval_min=pval_min, pval_max=pval_max)
  if(all(is.null(unlist(out))))
    stop('Difficulty constraint must have at least one argument specified')
  if(!is.null(modules))
    out$modules = as.integer(modules)
  class(out) = c('constraint_difficulty',class(out))
  out
}

to_lp.constraint_difficulty = function(f,...,items, pars, population_density, mst=NULL)
{
  s = seq(-6,6,0.1)
  w = population_density(s)
  w = w/sum(w)
  
  pvc = function(s,w)
  {
    out = list()
  
    pval = pars %>%
      group_by(.data$item_id) %>%
      arrange(.data$item_score) %>%
      do({
        r=.
        tibble(pvalue = sum(w*Es(r$item_score, r$beta, s)/max(r$item_score)))
      }) %>%
      ungroup() %>%
      semi_join(items,by='item_id') %>%
      arrange(.data$item_id) %>%
      pull(.data$pvalue)
    
    nit = nrow(items)
    
    if(!is.null(f$pval_min) || !is.null(f$pval_max))
    {
      xt = integer(nit)
      if(!is.null(f$pval_min))
        xt[pval<f$pval_min] = 1L
      if(!is.null(f$pval_max))
        xt[pval>f$pval_max] = 1L
      if(any(xt==1L))
        out$item_pval = list(xt=xt, type='=', rhs=0L)
    }
    if(!is.null(f$avg_pval_min))
    {
      out$pval_min = list(xt=pval-f$avg_pval_min, type = '>', rhs=0)
    }
    if(!is.null(f$avg_pval_max))
    {
      out$pval_max = list(xt=pval-f$avg_pval_max, type = '<', rhs=0)
    }
    out
  }

  if(is.null(mst))
  {
    return(pvc(s,w))
  } else
  {
    nmod = max(mst$rt_cut)
    if(is.null(f$modules))
      modules = 1:nmod
    else modules = f$modules
    if(any(modules<1 | modules>nmod))
      stop("invalid modules specifier in difficulty constraint")
    
    ps = Ps_test(mst$rt_pars,s)
    p_mod = t(apply(ps,1,rowsum,group=mst$rt_cut))
    rt_scores = 0:(ncol(ps) - 1)
    
    out = list()
    n=1
    
    for(i in modules)
    {
      wm = w * p_mod[,i]
      wm = wm/sum(wm)
      # per path would just mean a correction on the means
      for(ctr in pvc(s,wm))
      {
        ctr$bk = i
        out[[n]] = ctr
        n = n+1
      }
    }
  }
  out
}


# objective function ------------------------------------------------------

#' Property of the test to optimize
#' 
#' Minimize_items or maximize_items minimize of maximize some property of the items in the test. 
#' Random_test creates a random selection of items.
#' Maximize_information maximizes the information on specified points of the ability scale.
#' 
#' 
#' @param property unquoted column name or statement to be evaluated in the context of items, e.g. item_max_score, time_required. If NULL
#' optimizes the number of items in the test.
#' @param theta point(s) on the ability scale on which the information should be maximized
#' @param weight relative weights attached to the information on theta
#' 
#' 
#' 
minimize_items = function(property=NULL)
{
  qp = substitute(eval(quote(property)))
  out = list(qproperty = qp, direction='min')
  class(out) = c('objective_item','objective',class(out))
  out
}

#' @rdname minimize_items
maximize_items = function(property=NULL)
{
  qp = substitute(eval(quote(property)))
  out = list(qproperty = qp, direction='max')
  class(out) = c('objective_item','objective',class(out))
  out
}

#' @rdname minimize_items
random_test = function()
{
  out=list()
  class(out) = c('objective_random','objective')
  out
}

#' @rdname minimize_items
maximize_information = function(theta=0,weight=1)
{
  if(any(weight <= 0))
    stop('weights must be positive')
  
  if(length(theta)>1 && length(weight)==1)
    weight = rep(1,length(theta))
  
  if(length(theta) != length(weight))
    stop('Weight and theta do not have the same length')
  
  out = list(theta=theta,weight=weight)
  class(out) = c('objective_information','objective',class(out))
  out
}


to_lp.objective_item = function(f, ..., items, env)
{
  if(is.null(f$qproperty))
    xt = rep(1,nrow(items))
  else
    xt = as.numeric(eval_tidy(f$qproperty, data=items, env=env))
  
  list(xt=xt, direction = f$direction)
}

to_lp.objective_random = function(f, ..., items, mst=NULL)
{
  n = nrow(items)
  if(!is.null(mst)) n = n * max(mst$rt_cut)
  list(xt = runif(n), direction='max')
}

# do param and item check in main function
to_lp.objective_information = function(f,...,items, pars, mst=NULL)
{
  ii = function(theta, weight)
  {
    pars %>%
      group_by(.data$item_id) %>%
      arrange(.data$item_score) %>%
      do({
        r=.
        tibble(i = sum(info(r$item_score, r$beta, theta) * weight))
      }) %>%
      ungroup() %>%
      semi_join(items,by='item_id') %>%
      arrange(.data$item_id) %>%
      pull(.data$i)
  }  
  if(is.null(mst))
    return(list(xt=ii(f$theta,f$weight),direction='max'))
  
  ps = Ps_test(mst$rt_pars,f$theta)
  if(!is.matrix(ps))
    ps=matrix(ps,nrow=1)
  
  #rows: theta, cols: modules, normalized per row
  p_mod = t(apply(ps,1,rowsum,group=mst$rt_cut))
  
  xt = as.double(apply(p_mod,2, function(w) ii(f$theta, f$weight * w))) 

  list(xt=xt, direction='max')
}
