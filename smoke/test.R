
library(dplyr)

items = tibble(item_id=sprintf('itm%03i',1:200), 
               text_id=sample(c(letters,LETTERS),200,replace=TRUE),
               category = sample(c('know','apply','evaluate'),200,TRUE,prob=c(3,2,1)),
               difficulty = sample(1:3,200,TRUE,prob=c(1,2,1)))




d = pretest_design(items, max_nit=30, balance = c('category','difficulty'))


d = pretest_design(items, max_nit=30, 
                    balance = c('category','difficulty'),
                    friends = 'text_id')


d = pretest_design(items, max_nit=30)
d

rowSums(d)

# with this example presumably solvable without relaxation 1 in 10 times

library(igraph)

g = graph_from_adjacency_matrix(d,mode='undirected')

plot(g)





# for(i in 1:max_iter)
# {
#   mn = which.min(n)
#   mx = which.max(n)
#   rng = n[mx] - n[mn]
#   swp = NULL
#   for(i1 in seq_along(s[[mn]]))
#   {
#     for(i2 in seq_along(s[[mx]]))
#     {
#       prop = abs((n[mx]-s[[mx]][i2] + s[[mn]][i1]) - (n[mn]-s[[mn]][i1] + s[[mx]][i2]))
#       if(  prop < rng )
#       {
#         #swap
#         swp = c(i1,i2)
#         rng = prop
#       }
#     }
#   }
#   if(is.null(swp)) break
#   tmp = s[[mn]][swp[1]]
#   s[[mn]][swp[1]] = s[[mx]][swp[2]]
#   s[[mx]][swp[2]] = tmp
#   n = sapply(s, sum)
# }
# if(max(n)-min(n) < qal)
# {
#   qal = max(n) - min(n)
#   res = s
# }
# if(qal <= crit) break
# }
