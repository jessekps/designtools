
lpsolve = function(lp) solve(lp) %in% c(0,1,9,12)

# in case of infeasibility we can replace the objective function by a penalized constraint
# if it turns out positive/negative then this constraint (psbl in combination with another)
# is infeasible

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
	  
	  get_model = function(group=NULL, item_use = NULL, constraints=NULL)
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
	    l = make.lp(length(tp), M*nbk)
	    set.type(l, 1:(M*nbk), 'binary')
	    set.objfn(l, g(objf$w))
	    lp.control(l, sense=objf$type)
	    
	    for(i in seq_along(rh))
	    {
	      if(cname[i] %in% constraints)
	      {
  	      if(nbk == 1L || is.na(booklet[i]))
  	        add.constraint(l, g(cm[[i]]), tp[i], rh[[i]])
  	      else
  	        add.constraint(l, g(cm[[i]]), tp[i], rh[[i]], indices = (1L+(booklet[i]-1L)*M):(booklet[i]*M))
	      }
	     }
	    
	    if(!is.null(item_use) && nbk>1)
	    {
	      w = rep(1L,nbk)
	      for(i in 1:M)
	      {
	        add.constraint(l,w,'=', item_use, indices=seq(i, nbk*M, M))
	      }
	    }
	    
	    l
	  },
	  solve_model = function(..., fail=c('abort','investigate','force'), optimal=TRUE)
	  {
	    res=NULL
	    fail = match.arg(fail)
	    lpm = get_model(...)
	    if(!optimal)
	      lp.control(lpm, break.at.first=TRUE)
	    
	    status = solve(lpm)
	    
	    if(status %in% c(0,1,9,12))
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
