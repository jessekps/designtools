# all the necessary psychometrics
# P, E and information on scores c(0,item_score) for theta
Ps = function(item_score, beta, theta)
{
  a = item_score - lag(item_score,default=0)
  t(sapply(theta, function(tht){
    e = exp(cumsum(a*(tht-beta)))
    c(1,e)/(1+sum(e))
  }))
}

Es = function(item_score, beta, theta)
{
  Ps(item_score, beta, theta) %*% c(0L, item_score)
}

info = function(item_score, beta, theta)
{
  P = Ps(item_score, beta, theta)
  E = P %*% c(0L, item_score)
  rowSums(sapply(c(0L,item_score),'-',E)**2 * P)
}

info_test = function(pars,theta)
{
  rowSums(sapply(split(pars,pars$item_id), function(p) info(p$item_score,p$beta,theta)))
}

Ps_test = function(pars, theta)
{
  i = split(pars,pars$item_id)
  maxs = sum(sapply(i, function(r) max(r$item_score)))
  pt = array(0,c(length(theta),maxs+1,2))
  pt[,1,2] = 1
  
  ms = 0L
  j=0L
  for(r in i)
  {
    p = Ps(r$item_score,r$beta,theta)
    score = c(0,r$item_score)
    for(xi in 1:ncol(p))
    {
      for(st in 0:ms)
      {
        pt[,st+score[xi]+1,1+j] = pt[,st+score[xi]+1,1+j] + pt[,st+1,2-j] * p[,xi]
      }
    }
    ms = ms + max(r$item_score)
    j = 1L - j
    pt[,,1+j]=0
  }
  pt[,,2-j]
}

