


lpsolve = function(lp) solve(lp) %in% c(0,1,9,12)

scale_objective = function(val)
{
  if(any(val<0))
    stop('negative values in objective are not allowed')
  val/max(val)
}

LPM = setRefClass('LPM',
	fields = list(objf = 'list', cm = 'list', tp = 'character', rh = 'list', nit='integer',
	              nbk = 'integer', booklet='integer', cname='character'),
	methods = list(
	  initialize = function(weights, type=c('min','max'), nbooklets=1L, nitems=length(weights))
	  {
	    type=match.arg(type)
	    nit <<- as.integer(nitems)
	    nbk <<- as.integer(nbooklets)
	    
	    if(length(weights)==1)
	      weights = rep(weights, nit*nbk)
	    else if(length(weights) == nit)
	      weights = rep(weights,nbk)
	    
	    objf <<- list(w=weights, type=type)
	  },
	  add_constraint = function(xt, type = c("<=", "=", ">="), rhs, bk=NA_integer_, all_bk=FALSE,
	                            name=NULL)
	  {
	    if(length(xt)==1)
	      xt = rep.int(xt,nit)
	    if(is.null(name))
	      name = sprintf('constr%03i',length(cm) + 1)
	    if(all_bk)
	    {
	      for(i in 1:nbk)
	      {
	        add_constraint(xt,type=type,rhs=rhs,bk=i,name=name)
	      }
	    } else
	    {
	      cm <<- append(cm,list(xt))
	      tp <<- c(tp,match.arg(type))
	      rh <<- append(rh,rhs)
	      booklet <<- c(booklet,bk) 
	      cname <<- c(cname,name)
	    }
	  },
	  
	  get_constraint_names = function() cname,
	  
	  # single_booklet=NULL gives the combined model
	  # presently, we have no cross booklet constraints so it is always save
	  # and almost always faster and sometimes better to solve per booklet
	  get_model = function(group=NULL, constraints=NULL, single_booklet=NULL) 
	  {
	    if(is.null(constraints))
	      constraints = cname
	    
	    constraints = unique(constraints)
	    
	    if(!is.null(group))
	    {
	      M = n_distinct(group)
        group2 = sprintf("%05i %s", rep(1:nbk,each=nit), rep(group,nbk))
	      g = function(x) rowsum(x,if(length(x)==nit) group else group2)
	    } else
	    {
	      M = nit
	      g = identity
	    }
	    if(is.null(single_booklet) || nbk == 1)
	    {
	      l = make.lp(0, M*nbk)
	      set.type(l, 1:(M*nbk), 'binary')
	      set.objfn(l, scale_objective(g(objf$w)))
	    } 
	    else
	    {
	      stopifnot(single_booklet %in% 1:nbk)
	      l = make.lp(0, M)
	      set.type(l, 1:M, 'binary')
	      set.objfn(l, scale_objective(g(objf$w)[1:M + (single_booklet-1L)*M]))
	    }
	    
	    lp.control(l, sense=objf$type)
	    
	    for(i in seq_along(rh))
	    {
	      if(cname[i] %in% constraints)
	      {
  	      if(nbk == 1L || is.na(booklet[i]))
  	        add.constraint(l, g(cm[[i]]), tp[i], rh[[i]])
  	      else if(is.null(single_booklet))
  	        add.constraint(l, g(cm[[i]]), tp[i], rh[[i]], indices = (1L+(booklet[i]-1L)*M):(booklet[i]*M))
	        else if(single_booklet == booklet[i])
	        {
	          add.constraint(l, g(cm[[i]]), tp[i], rh[[i]])
	        } 
	      }
	     }
	    
	    l
	  },
	  solve_separate = function(..., fail=c('abort','investigate','force'), optimal=TRUE, timeout=120)
	  {
	    slv = lapply(1:nbk,function(i){
	      solve_model(..., single_booklet=i, fail=fail, optimal=optimal,timeout=timeout)
	    })
	    if(all(sapply(slv,'[[','success')))
	    {
	      return(list(success=TRUE,result = unlist(lapply(slv,'[[','result'))))
	    } 
	    list(success=FALSE,result=NULL)
	  },
	  solve_model = function(..., fail=c('abort','investigate','force'), optimal=TRUE, timeout=120)
	  {
	    res=NULL
	    fail = match.arg(fail)
	    lpm = get_model(...)
	    if(!optimal)
	      lp.control(lpm, break.at.first=TRUE)
	    
	    if(timeout>0)
	      lp.control(lpm,timeout=timeout)
	    
	    status = solve(lpm)
	    
	    if(status %in% c(0,1,9,11,12))
	    {
	      res = list(success = TRUE, result = as.logical(get.variables(lpm)),status=status)
	    }
	    else if(fail == 'abort')
	    {
        res = list(success=FALSE, result=NULL,status=status)
	    }
	    else
	    {
 	    }
	    #delete.lp(lp)
	    res
	  }
	))
