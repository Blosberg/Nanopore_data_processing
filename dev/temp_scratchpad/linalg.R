#  General linear algebra functions 


#------------------------------------
proj <- function( u, v ) #---- PROJECTS v ONTO u 
{
  temp_u = dim( as.matrix(u) )
  temp_v = dim( as.matrix(v) )
  
  if( (temp_u[2] != 1 || temp_v[2] != 1) || temp_u[1] != temp_v[1] ){print("yup")}
  
  #  { stop("input vectors incompatible with projection") }

result =  (sum(u*v) / sum(u*u)) *u
    
return(result)
}

  # outliers     = which( colSums( all_olapreads_normedkvecs * all_olapreads_normedkvecs ) >1 )
  #  outlier_mean = rowMeans( all_olapreads_normedkvecs [,outliers])
  #  omu = outlier_mean*(1/sqrt(sum(outlier_mean*outlier_mean) ) )
  #  first_projection <- omu %*%   all_olapreads_normedkvecs
  #--- unit vector in the direction of the outlier means.
  # ----------------------------------------------------
  # primary display vector is the x-axis of the 2-D plot
  # primary_disp_vec = outlier_mean/(sqrt(sum(outlier_mean*outlier_mean))); outlier_mean
  # proj ( omu, all_olapreads_normedkvecs[,1] ) #---- PROJECTS v ONTO u 


#------------------------------------

gen_orthonorm_vecs <- function(vec_in) 
{ 
  temp=dim(vec_in)
  
  if(temp[1] <= 1 || temp[2] != 1)
  {
    stop(paste("vector input of dimension",temp[1]," x ", temp[2],"; incompatible with orthogonal generation") )
  }
  N=temp[1]
  
  result     = matrix(0,N,N)
  result[,1] = vec_in/(norm(vec_in,type="F"))
  
  for(i in c(2:N))
  {
    temp_a   = matrix(0,N,1)
    temp_b   = matrix(0,N,1)
    
    temp_a[1:(i-1)] = 1
    temp_a[i]       = 2
    
    for( j in c(1:(i-1) ) )
    {
      temp_a = temp_a - proj(result[, j, drop = FALSE],temp_a) #-- subtract projections onto each other dimension
    }
    
    result[,i] = temp_a/(norm(temp_a,type="F"))
  }
  
  return(result)
}

