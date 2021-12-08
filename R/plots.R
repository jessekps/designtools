
#' IRT plot of mst design effects
#' 
#' @param md mst design, output of function \code{\link{mst_design}}
#' @param pars item parameters
#' @param populaton_density density function for the population, the default is standard normal. 
#' @param type information or SEM
#' @param infolines for the information plot vertical lines at these points
#' 
#' @returns ggplot2 like object to be plotted
#' 
plot_mst_design = function(md,pars, populaton_density=dnorm,type=c('sem','info'), 
                    infolines=NULL)
{
  type=match.arg(type)
  rt_pars = md$design %>%
    filter(.data$module==0) %>%
    distinct(.data$item_id) %>%
    inner_join(pars,by='item_id')
  
  grd = seq(-4,4,.02)
  w = populaton_density(grd)
  w = w/sum(w)
  ps = Ps_test(rt_pars,grd)
  p_mod = t(apply(ps,1,rowsum,group=md$routing$next_module))
  
  
  dat = lapply(split(md$design,md$design$booklet), function(ds){
    i = ds$booklet[1]
    wm = w * p_mod[,i]
    wm = wm/sum(wm)
    info = info_test(inner_join(ds,pars,by='item_id'),grd)
    
    tibble(theta=grd,booklet=i, 
           info=info,
           sem=1/sqrt(info),
           w=wm)               
  }) %>%
    bind_rows() %>%
    mutate(booklet=factor(.data$booklet))
  
  dat2 = dat %>%
    group_by(.data$theta) %>%
    summarise(avg_sem = sum(.data$w*.data$sem)/sum(.data$w),
              avg_info = sum(.data$w*.data$info)/sum(.data$w))
  
  if(type=='info'){
    p = ggplot(dat,aes(x=.data$theta,y=.data$info)) +
      geom_line(aes(colour=.data$booklet)) +
      geom_line(aes(y=.data$avg_info),data=dat2,linetype='dashed') +
      ylab('information') +
      theme_bw()
    
    if(!is.null(infolines))
    {
      ymx = sapply(infolines,function(tht) dat2$avg_info[which.min(abs(dat2$theta-tht))])
      p = p + geom_segment(data=tibble(theta=infolines,ymx=ymx),
                           aes(xend=.data$theta,yend=.data$ymx),y=0,linetype='dotdash')
    }
  } else
  {
    p = ggplot(dat,aes(x=.data$theta,y=.data$sem)) +
      geom_line(aes(colour=.data$booklet)) +
      geom_line(aes(y=.data$avg_sem),data=dat2,linetype='dashed') +
      ylab('sem') +
      theme_bw()
  }
  
  m = ggplot(dat,aes(x=.data$theta,y=.data$w)) +
    geom_line(aes(colour=.data$booklet),show.legend=F) +
    geom_area(aes(fill=.data$booklet),show.legend=FALSE,position='identity',alpha=0.2) +
    theme_void() 
  
  res=list(p=p,m=m)
  class(res) = c("ds_marginal_plot",class(res))
  res
}

print.ds_marginal_plot = function(x,...)
{
  ggarrange(x$m,x$p,heights=c(2,7),ncol=1)
  invisible(NULL)
}

plot.ds_marginal_plot = function(x,...)
{
  print(x)
}

`+.ds_marginal_plot` = function(mgg,add)
{
  if(inherits(add,'labels'))
  {
    addm = add
    addm[!names(addm) %in% c('title','subtitle')] = NULL
    if(length(addm)>0)
      mgg$m = mgg$m + addm
    
    add[names(add) %in% c('title','subtitle')] = NULL
    if(length(add)==0)
      return(mgg)
  }
  
  mgg$p = mgg$p + add
  mgg
}