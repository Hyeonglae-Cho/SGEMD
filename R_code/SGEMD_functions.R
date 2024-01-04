# Necessary functions for executing numerical experiments and visualization
# Install packages
install.packages('igraph')
install.packages('gasper')
install.packages('ggplot2')
install.packages('Matrix')
install.packages('plotly')
install.packages('writexl')
install.packages('EbayesThresh')
install.packages('tidyverse')

library(gasper)
library(igraph)
library(ggplot2)
library(Matrix)
library(plotly)
library(writexl)
library(EbayesThresh)
library(tidyverse)


# Obtain knn vertices from a vertex

knn_vertices = function(graph, k) {
  dist = c()
  coord = graph$xy
  n = nrow(coord)
  idx = sample(x=1:n, size=1)
  for (j in 1:n) {
    dist = c(dist, sqrt((graph$xy[idx,'x'] - graph$xy[j,'x'])^2 + (graph$xy[idx,'y'] - graph$xy[j,'y'])^2))
  }
  nearest_vertices = order(dist)[1:(k+1)]
  return(nearest_vertices)
}


# Obtain a matrix representing distances between the vertices

distance_matrix = function(dimension, adjacency_matrix, coordinates) {
  n = dimension
  for (i in 1:n) {
    for (j in 1:n) {
      if (adjacency_matrix[i,j] != 0) {
        adjacency_matrix[i,j] = sqrt((coordinates[i,'x'] - coordinates[j,'x'])^2 + (coordinates[i,'y'] - coordinates[j,'y'])^2)
      } 
    }
  }
  return(adjacency_matrix)
}


# Obtain weighted adjacency matrix by connecting the vertices withing certain distance 'delta'

weighted_ad_matrix = function(dimension, delta, coordinates) {
  n = dimension
  ad_mat = matrix(nrow=n, ncol=n)
  for (i in 1:n) {
    for (j in 1:n) {
      dist_ij = sqrt((coordinates[i,'x'] - coordinates[j,'x'])^2 + (coordinates[i,'y'] - coordinates[j,'y'])^2)
      if (dist_ij <= delta) {
        ad_mat[i,j] = exp((-(dist_ij^2)) / (2*(delta^2)))
      } else {
        ad_mat[i,j] = 0
      }
    }
  }
  diag(ad_mat) = 0
  return(ad_mat)
}


# Obtain weighted adjacency matrix by connecting nearest 'k' vertices

weighted_ad_matrix_k = function(dimension, k, coordinates) {
  n = dimension
  ad_mat = matrix(nrow=n, ncol=n)
  dst_mat = matrix(nrow=n, ncol=n)
  for (i in 1:n) {
    dst_i = c()
    for (j in 1:n) {
      dst_ij = sqrt((coordinates[i,'x'] - coordinates[j,'x'])^2 + (coordinates[i,'y'] - coordinates[j,'y'])^2)
      dst_i = c(dst_i, dst_ij)
    }
    idx_connected_vertices = order(dst_i)[1:k+1]
    dst_i_threshold = rep(0,n) 
    dst_i_threshold[idx_connected_vertices] = dst_i[idx_connected_vertices]
    dst_mat[i,] = dst_i_threshold
  }
  delta = max(dst_mat)
  ad_mat = exp((-(dst_mat^2)) / (2*(delta^2)))
  ad_mat[dst_mat == 0] = 0
  diag(ad_mat) = 0
  
  for (k in 1:n) {
    for (l in union(which(ad_mat[k,] != 0), which(ad_mat[,k] != 0))) {
      if (ad_mat[k,l] != 0 & ad_mat[l,k] == 0) {
        ad_mat[l,k] = ad_mat[k,l]
      } else if (ad_mat[k,l] == 0 & ad_mat[l,k] != 0) {
        ad_mat[k,l] = ad_mat[l,k]
      }
    }
  }
  
  return(ad_mat)
}


# Obtain weighted adjacency matrix from graph

weighted_ad_matrix_from_graph = function(graph, delta) {
  ad_mat = as.matrix(graph, matrix.type=c('adjacency'))
  n = dim(ad_mat)
  for (i in 1:n) {
    for (j in 1:n) {
      dist_ij = sqrt((graph$xy[i,'x'] - graph$xy[j,'x'])^2 + (graph$xy[i,'y'] - graph$xy[j,'y'])^2)
      if (dist_ij <= delta) {
        ad_mat[i,j] = exp((-(dist_ij^2)) / (2*(delta^2)))
      } else {
        ad_mat[i,j] = 0
      }
    }
  }
  diag(ad_mat) = 0
  return(ad_mat)
}


# Obtain extrema (maxima and minima) of a signal on a graph

get_extrema = function(graph, signal) {
  maxima_list = c()
  minima_list = c()
  admat = graph$sA
  
  for (i in 1:dim(admat)[1]) {
    vertex = i
    connected_vertices = which(admat[vertex,] != 0)
    
    if (length(connected_vertices) == sum(signal[vertex] >= signal[connected_vertices])) {
      maxima_list = c(maxima_list, vertex)
    } else if (length(connected_vertices) == sum(signal[vertex] <= signal[connected_vertices])) {
      minima_list = c(minima_list, vertex)
    }
    
  }
  extrema = list(maxima_list, minima_list)
  names(extrema) = c('maxima_list', 'minima_list')
  return(extrema)
}


# Graph interpolation and graph empirical mode decomposition is proposed in proposed in Tremblay, N., Borgnat, P., and Flandrin, P. (2014). Graph empirical mode decomposition. In 2014 22nd European Signal Processing Conference (EUSIPCO), pages 2350–2354.
# Graph interpolation 

graph_interpolation = function(graph, signal, vertices) {
  adjacency_matrix = graph$sA
  Laplacian = diag(rowSums(adjacency_matrix)) - adjacency_matrix
  n = dim(Laplacian)[1]
  s.known = vertices
  s.unknown = setdiff(1:n, s.known)
  
  ordered_Laplacian = Laplacian[c(s.known, s.unknown), c(s.known, s.unknown)]
  ordered_signal = signal[c(s.known, s.unknown)]
  
  RT = ordered_Laplacian[(length(s.known)+1):n,1:length(s.known)]
  Lu = ordered_Laplacian[(length(s.known)+1):n,(length(s.known)+1):n]
  sb = ordered_signal[1:length(s.known)]
  
  su = solve(Lu, -RT %*% sb)
  
  signal[s.unknown] = su
  
  return(signal)
}


# Graph empirical mode decomposition

GEMD = function(graph, signal, K) {
  elements = list()
  m = signal
  
  for (i in 1:K) {
    s0 = m
    s1 = m
    mu = m
    
    j = 1
    while (1) {
      extrema = get_extrema(graph, s1)
      maxima = extrema$maxima_list
      minima = extrema$minima_list
      
      inter.max = graph_interpolation(graph, s1, maxima)
      inter.min = graph_interpolation(graph, s1, minima)
      
      mu = (inter.max + inter.min) / 2
      
      if (mu %*% mu < (s1 %*% s1)/1000) {
        break
      }
      
      s1 = s1 - mu
      
      j = j + 1
      
    }
    elements = c(elements, list(s1))
    m = s0-s1
  }
  elements = c(elements, list(m))
  return(elements)
}


# Obtain average signal of knn vertices

knn_average = function(graph, idx, signal, Laplacian_matrix, k) {
  dist = c()
  L = Laplacian_matrix
  
  while(1) {
    if (sum(L[idx,] != 0) - 1 >= k) {
      break
    } else {L = L %*% L}
    
  }
  
  for (j in which(L[idx,] != 0)) {
    dist = c(dist, sqrt((graph$xy[idx,'x'] - graph$xy[j,'x'])^2 + (graph$xy[idx,'y'] - graph$xy[j,'y'])^2))
  }
  
  nearest_vertices = order(dist)[1:(k+1)]
  nearest_vertices = nearest_vertices[! nearest_vertices %in% idx]
  average = mean(signal[nearest_vertices])
  
  return(average)
}


# Obtain average signal of adjacent vertices

neighbor_average = function(graph, idx, signal) {
  A = graph$sA
  
  neighbor_vertices = which(A[idx, ] != 0)
  constants = A[idx, neighbor_vertices]
  weights = constants / sum(constants)
  
  weighted_average = sum(weights * signal[neighbor_vertices])
  
  return(weighted_average)
}


# Graph denoising

GFT.smoothing.ebayesthresh = function(graph, signal) {
  adjacency_matrix = graph$sA
  L = diag(rowSums(adjacency_matrix)) - adjacency_matrix
  eigenvectors = eigen(L)$vectors
  
  GFT = t(eigenvectors) %*% signal # Calculate Graph Fourier Transforms
  GFT.Thresh = ebayesthresh(GFT, verbose = TRUE, threshrule = 'soft') # Get Empirical Bayes Thresholding from GFT
  
  GFT.inverse.smoothing = eigenvectors %*% GFT.Thresh$muhat # Get smoothed signal through Inverse Graph Fourier Transform with thresholded graph Fourier coefficients
  
  return(GFT.inverse.smoothing)
}


# Statistical graph empirical mode decomposition

GEMD.refl.GFT.ebayesthresh = function(graph, signal, K, method, delta, c, connection = 'dst') {
  
  #Set Reflected Graph
  ad_mat_tmp = graph$sA
  xy_tmp = graph$xy
  n_vertices = dim(ad_mat_tmp)[1]
  
  left_end = min(xy_tmp[,1]) ; right_end = max(xy_tmp[,1]) ; lower_end = min(xy_tmp[,2]) ; upper_end = max(xy_tmp[,2])
  
  left_reflection_idx = as.numeric(rownames(xy_tmp[xy_tmp[,1] > left_end & xy_tmp[,1] <= left_end + c*delta, ])) ; n_left = length(left_reflection_idx)
  right_reflection_idx = as.numeric(rownames(xy_tmp[xy_tmp[,1] < right_end & xy_tmp[,1] >= right_end - c*delta, ])) ; n_right = length(right_reflection_idx)
  lower_reflection_idx = as.numeric(rownames(xy_tmp[xy_tmp[,2] > lower_end & xy_tmp[,2] <= lower_end + c*delta, ])) ; n_lower = length(lower_reflection_idx)
  upper_reflection_idx = as.numeric(rownames(xy_tmp[xy_tmp[,2] < upper_end & xy_tmp[,2] >= upper_end - c*delta, ])) ; n_upper = length(upper_reflection_idx)
  
  left_upper_reflection_idx = intersect(left_reflection_idx, upper_reflection_idx) ; n_left_upper = length(left_upper_reflection_idx)
  right_upper_reflection_idx = intersect(right_reflection_idx, upper_reflection_idx) ; n_right_upper = length(right_upper_reflection_idx)
  left_lower_reflection_idx = intersect(left_reflection_idx, lower_reflection_idx) ; n_left_lower = length(left_lower_reflection_idx)
  right_lower_reflection_idx = intersect(right_reflection_idx, lower_reflection_idx) ; n_right_lower = length(right_lower_reflection_idx)
  
  reflection_idx_side = c(left_reflection_idx, right_reflection_idx, lower_reflection_idx, upper_reflection_idx)
  reflection_idx_apex = c(left_upper_reflection_idx, right_upper_reflection_idx, left_lower_reflection_idx, right_lower_reflection_idx)
  
  n_reflection = length(c(reflection_idx_side, reflection_idx_apex))
  n_reflection_side = length(reflection_idx_side)
  n_reflection_apex = length(reflection_idx_apex)
  

  left_reflection = xy_tmp[left_reflection_idx,]
  left_reflection$x = left_reflection$x * (-1) + 2 * left_end
  cal_left_refl = union(which(max(left_reflection$x)-delta < left_reflection$x & left_reflection$x <= max(left_reflection$x)),
                        union(which(max(left_reflection$y)-delta < left_reflection$y & left_reflection$y <= max(left_reflection$y)),
                              which(min(left_reflection$y) <= left_reflection$y & left_reflection$y < min(left_reflection$y)+delta)))
  
  right_reflection = xy_tmp[right_reflection_idx,]
  right_reflection$x = right_reflection$x * (-1) + 2 * right_end
  cal_right_refl = union(which(min(right_reflection$x) <= right_reflection$x & right_reflection$x < min(right_reflection$x)+delta),
                         union(which(max(right_reflection$y)-delta < right_reflection$y & right_reflection$y <= max(right_reflection$y)),
                               which(min(right_reflection$y) <= right_reflection$y & right_reflection$y < min(right_reflection$y)+delta)))
  
  lower_reflection = xy_tmp[lower_reflection_idx,]
  lower_reflection$y = lower_reflection$y * (-1) + 2 * lower_end
  cal_lower_refl = union(which(max(lower_reflection$y)-delta < lower_reflection$y & lower_reflection$y <= max(lower_reflection$y)),
                         union(which(max(lower_reflection$x)-delta < lower_reflection$x & lower_reflection$x <= max(lower_reflection$x)),
                               which(min(lower_reflection$x) <= lower_reflection$x & lower_reflection$x < min(lower_reflection$x)+delta)))
  
  upper_reflection = xy_tmp[upper_reflection_idx,]
  upper_reflection$y = upper_reflection$y * (-1) + 2 * upper_end
  cal_upper_refl = union(which(min(upper_reflection$y) <= upper_reflection$y & upper_reflection$y < min(upper_reflection$y)+delta),
                         union(which(max(upper_reflection$x)-delta < upper_reflection$x & upper_reflection$x <= max(upper_reflection$x)),
                               which(min(upper_reflection$x) <= upper_reflection$x & upper_reflection$x < min(upper_reflection$x)+delta)))
  
  left_upper_reflection = xy_tmp[left_upper_reflection_idx,]
  left_upper_reflection$x = left_upper_reflection$x * (-1) + 2 * left_end ; left_upper_reflection$y = left_upper_reflection$y * (-1) + 2 * upper_end
  cal_left_upper_refl = union(which(min(left_upper_reflection$y) <= left_upper_reflection$y & left_upper_reflection$y < min(left_upper_reflection$y)+delta),
                              which(max(left_upper_reflection$x)-delta < left_upper_reflection$x & left_upper_reflection$x <= max(left_upper_reflection$x)))
  
  right_upper_reflection = xy_tmp[right_upper_reflection_idx,]
  right_upper_reflection$x = right_upper_reflection$x * (-1) + 2 * right_end ; right_upper_reflection$y = right_upper_reflection$y * (-1) + 2 * upper_end
  cal_right_upper_refl = union(which(min(right_upper_reflection$y) <= right_upper_reflection$y & right_upper_reflection$y < min(right_upper_reflection$y)+delta),
                               which(min(right_upper_reflection$x) <= right_upper_reflection$x & right_upper_reflection$x < min(right_upper_reflection$x)+delta))
  
  left_lower_reflection = xy_tmp[left_lower_reflection_idx,]
  left_lower_reflection$x = left_lower_reflection$x * (-1) + 2 * left_end ; left_lower_reflection$y = left_lower_reflection$y * (-1) + 2 * lower_end
  cal_left_lower_refl = union(which(max(left_lower_reflection$y)-delta < left_lower_reflection$y & left_lower_reflection$y <= max(left_lower_reflection$y)),
                              which(max(left_lower_reflection$x)-delta < left_lower_reflection$x & left_lower_reflection$x <= max(left_lower_reflection$x)))
  
  right_lower_reflection = xy_tmp[right_lower_reflection_idx,]
  right_lower_reflection$x = right_lower_reflection$x * (-1) + 2 * right_end ; right_lower_reflection$y = right_lower_reflection$y * (-1) + 2 * lower_end
  cal_right_lower_refl = union(which(max(right_lower_reflection$y)-delta < right_lower_reflection$y & right_lower_reflection$y <= max(right_lower_reflection$y)),
                               which(min(right_lower_reflection$x) <= right_lower_reflection$x & right_lower_reflection$x < min(right_lower_reflection$x)+delta))
  
  cal_center = union(which((min(xy_tmp$x) <= xy_tmp$x & xy_tmp$x < min(xy_tmp$x)+delta)),
                     union(which((max(xy_tmp$x)-delta < xy_tmp$x & xy_tmp$x <= max(xy_tmp$x))),
                           union(which((min(xy_tmp$y) <= xy_tmp$y & xy_tmp$y < min(xy_tmp$y)+delta)),
                                 which((max(xy_tmp$y)-delta < xy_tmp$y & xy_tmp$y <= max(xy_tmp$y))))))
  
  n_center = length(cal_center)
  
  cal_idx = c(cal_center, n_vertices+cal_left_refl, n_vertices+n_left+cal_right_refl, n_vertices+n_left+n_right+cal_lower_refl, n_vertices+n_left+n_right+n_lower+cal_upper_refl,
              n_vertices+n_left+n_right+n_lower+n_upper+cal_left_upper_refl, n_vertices+n_left+n_right+n_lower+n_upper+n_left_upper+cal_right_upper_refl,
              n_vertices+n_left+n_right+n_lower+n_upper+n_left_upper+n_right_upper+cal_left_lower_refl, n_vertices+n_left+n_right+n_lower+n_upper+n_left_upper+n_right_upper+n_left_lower+cal_right_lower_refl)
  
  
  reflection_coord_side = rbind(left_reflection, right_reflection, lower_reflection, upper_reflection)
  reflection_coord_apex = rbind(left_upper_reflection, right_upper_reflection, left_lower_reflection, right_lower_reflection)
  
  
  graph.refl.coord = rbind(xy_tmp, reflection_coord_side, reflection_coord_apex)
  
  ad_matrix.reflection = cbind(ad_mat_tmp, Matrix(0, nrow = n_vertices, ncol = n_reflection, sparse = TRUE))
  ad_matrix.reflection = rbind(ad_matrix.reflection, Matrix(0, nrow = n_reflection, ncol = n_vertices + n_reflection))
  n_ad_matrix = dim(ad_matrix.reflection)[1]
  
  re_ad_matrix = ad_matrix.reflection[(n_vertices+1):n_ad_matrix,(n_vertices+1):n_ad_matrix]
  
  re_ad_matrix[1:n_left,1:n_left] = ad_mat_tmp[left_reflection_idx,left_reflection_idx]
  re_ad_matrix[(n_left+1):(n_left+n_right),(n_left+1):(n_left+n_right)] = ad_mat_tmp[right_reflection_idx,right_reflection_idx]
  re_ad_matrix[(n_left+n_right+1):(n_left+n_right+n_lower),(n_left+n_right+1):(n_left+n_right+n_lower)] = ad_mat_tmp[lower_reflection_idx,lower_reflection_idx]
  re_ad_matrix[(n_left+n_right+n_lower+1):(n_left+n_right+n_lower+n_upper),(n_left+n_right+n_lower+1):(n_left+n_right+n_lower+n_upper)] = ad_mat_tmp[upper_reflection_idx,upper_reflection_idx]
  
  try({
  re_ad_matrix[(n_reflection-n_right_lower-n_left_lower-n_right_upper-n_left_upper+1):(n_reflection-n_right_lower-n_left_lower-n_right_upper),
               (n_reflection-n_right_lower-n_left_lower-n_right_upper-n_left_upper+1):(n_reflection-n_right_lower-n_left_lower-n_right_upper)] = ad_mat_tmp[left_upper_reflection_idx,left_upper_reflection_idx]
  re_ad_matrix[(n_reflection-n_right_lower-n_left_lower-n_right_upper+1):(n_reflection-n_right_lower-n_left_lower),
               (n_reflection-n_right_lower-n_left_lower-n_right_upper+1):(n_reflection-n_right_lower-n_left_lower)] = ad_mat_tmp[right_upper_reflection_idx,right_upper_reflection_idx]
  re_ad_matrix[(n_reflection-n_right_lower-n_left_lower+1):(n_reflection-n_right_lower),(n_reflection-n_right_lower-n_left_lower+1):(n_reflection-n_right_lower)] = ad_mat_tmp[left_lower_reflection_idx,left_lower_reflection_idx]
  re_ad_matrix[(n_reflection-n_right_lower+1):n_reflection,(n_reflection-n_right_lower+1):n_reflection] = ad_mat_tmp[right_lower_reflection_idx,right_lower_reflection_idx]})
  
  ad_matrix.reflection[(n_vertices+1):n_ad_matrix,(n_vertices+1):n_ad_matrix] = re_ad_matrix
  
  if (connection == 'dst') {
    dst_max = delta
  } else if (connection == 'neighbor') {
    dst_mat = matrix(nrow=n_vertices, ncol=n_vertices)
    for (i in 1:n_vertices) {
      dst_i = c()
      for (j in 1:n_vertices) {
        dst_ij = sqrt((xy_tmp[i,'x'] - xy_tmp[j,'x'])^2 + (xy_tmp[i,'y'] - xy_tmp[j,'y'])^2)
        dst_i = c(dst_i, dst_ij)
      }
      dst_mat[i,] = dst_i
    }
    dst_max = max(dst_mat[as.logical(ad_matrix.reflection[1:n_vertices, 1:n_vertices] != 0)])
  }

  
  # Connectivity between center and the other areas
  for (i in cal_center) {
    
    for (j in cal_idx[(n_center+1):length(cal_idx)]) {
      dist_ij = sqrt((graph.refl.coord[i,'x'] - graph.refl.coord[j,'x'])^2 + (graph.refl.coord[i,'y'] - graph.refl.coord[j,'y'])^2)
      if (dist_ij <= delta & dist_ij != 0) {
        ad_matrix.reflection[i,j] = exp((-(dist_ij^2)) / (2*(dst_max^2)))
      } else {
        ad_matrix.reflection[i,j] = 0
      }
    }
  }
  
  # Connectivity between left upper and left, upper
  for (i in n_vertices+n_left+n_right+n_lower+n_upper+cal_left_upper_refl) {
    
    for (j in n_vertices+cal_left_refl) {
      dist_ij = sqrt((graph.refl.coord[i,'x'] - graph.refl.coord[j,'x'])^2 + (graph.refl.coord[i,'y'] - graph.refl.coord[j,'y'])^2)
      if (dist_ij <= delta & dist_ij != 0) {
        ad_matrix.reflection[i,j] = exp((-(dist_ij^2)) / (2*(dst_max^2)))
      } else {
        ad_matrix.reflection[i,j] = 0
      }
    }
    
    for (j in n_vertices+n_left+n_right+n_lower+cal_upper_refl) {
      dist_ij = sqrt((graph.refl.coord[i,'x'] - graph.refl.coord[j,'x'])^2 + (graph.refl.coord[i,'y'] - graph.refl.coord[j,'y'])^2)
      if (dist_ij <= delta & dist_ij != 0) {
        ad_matrix.reflection[i,j] = exp((-(dist_ij^2)) / (2*(dst_max^2)))
      } else {
        ad_matrix.reflection[i,j] = 0
      }
    }
  }
  
  # Connectivity between right upper and right, upper
  for (i in n_vertices+n_left+n_right+n_lower+n_upper+n_left_upper+cal_right_upper_refl) {
    
    for (j in n_vertices+n_left+cal_right_refl) {
      dist_ij = sqrt((graph.refl.coord[i,'x'] - graph.refl.coord[j,'x'])^2 + (graph.refl.coord[i,'y'] - graph.refl.coord[j,'y'])^2)
      if (dist_ij <= delta & dist_ij != 0) {
        ad_matrix.reflection[i,j] = exp((-(dist_ij^2)) / (2*(dst_max^2)))
      } else {
        ad_matrix.reflection[i,j] = 0
      }
    }
    
    for (j in n_vertices+n_left+n_right+n_lower+cal_upper_refl) {
      dist_ij = sqrt((graph.refl.coord[i,'x'] - graph.refl.coord[j,'x'])^2 + (graph.refl.coord[i,'y'] - graph.refl.coord[j,'y'])^2)
      if (dist_ij <= delta & dist_ij != 0) {
        ad_matrix.reflection[i,j] = exp((-(dist_ij^2)) / (2*(dst_max^2)))
      } else {
        ad_matrix.reflection[i,j] = 0
      }
    }
  }
  
  # Connectivity between left lower and left, lower
  for (i in n_vertices+n_left+n_right+n_lower+n_upper+n_left_upper+n_right_upper+cal_left_lower_refl) {
    
    for (j in n_vertices+cal_left_refl) {
      dist_ij = sqrt((graph.refl.coord[i,'x'] - graph.refl.coord[j,'x'])^2 + (graph.refl.coord[i,'y'] - graph.refl.coord[j,'y'])^2)
      if (dist_ij <= delta & dist_ij != 0) {
        ad_matrix.reflection[i,j] = exp((-(dist_ij^2)) / (2*(dst_max^2)))
      } else {
        ad_matrix.reflection[i,j] = 0
      }
    }
    
    for (j in n_vertices+n_left+n_right+cal_lower_refl) {
      dist_ij = sqrt((graph.refl.coord[i,'x'] - graph.refl.coord[j,'x'])^2 + (graph.refl.coord[i,'y'] - graph.refl.coord[j,'y'])^2)
      if (dist_ij <= delta & dist_ij != 0) {
        ad_matrix.reflection[i,j] = exp((-(dist_ij^2)) / (2*(dst_max^2)))
      } else {
        ad_matrix.reflection[i,j] = 0
      }
    }
  }
  
  # Connectivity between left lower and left, lower
  for (i in n_vertices+n_left+n_right+n_lower+n_upper+n_left_upper+n_right_upper+n_left_lower+cal_right_lower_refl) {
    
    for (j in n_vertices+n_left+cal_right_refl) {
      dist_ij = sqrt((graph.refl.coord[i,'x'] - graph.refl.coord[j,'x'])^2 + (graph.refl.coord[i,'y'] - graph.refl.coord[j,'y'])^2)
      if (dist_ij <= delta & dist_ij != 0) {
        ad_matrix.reflection[i,j] = exp((-(dist_ij^2)) / (2*(dst_max^2)))
      } else {
        ad_matrix.reflection[i,j] = 0
      }
    }
    
    for (j in n_vertices+n_left+n_right+cal_lower_refl) {
      dist_ij = sqrt((graph.refl.coord[i,'x'] - graph.refl.coord[j,'x'])^2 + (graph.refl.coord[i,'y'] - graph.refl.coord[j,'y'])^2)
      if (dist_ij <= delta & dist_ij != 0) {
        ad_matrix.reflection[i,j] = exp((-(dist_ij^2)) / (2*(dst_max^2)))
      } else {
        ad_matrix.reflection[i,j] = 0
      }
    }
  }
  
  
  for (k in 1:dim(ad_matrix.reflection)[1]) {
    for (l in union(which(ad_matrix.reflection[k,] != 0), which(ad_matrix.reflection[,k] != 0))) {
      if (ad_matrix.reflection[k,l] != 0 & ad_matrix.reflection[l,k] == 0) {
        ad_matrix.reflection[l,k] = ad_matrix.reflection[k,l]
      } else if (ad_matrix.reflection[k,l] == 0 & ad_matrix.reflection[l,k] != 0) {
        ad_matrix.reflection[k,l] = ad_matrix.reflection[l,k]
      }
    }
  }
  
  graph.reflection = graph_from_adjacency_matrix(ad_matrix.reflection, mode=c('undirected'), weighted=TRUE)
  graph_attr(graph.reflection, 'xy') = graph.refl.coord
  graph_attr(graph.reflection, 'sA') = Matrix(ad_matrix.reflection, sparse=TRUE)
  
  # Discarding some disconnected parts, obtain connected graph 
  graph.reflection.connectivity = components(graph.reflection)
  graph.reflection.main.index = which(graph.reflection.connectivity$membership == which.max(graph.reflection.connectivity$csize))
  
  graph.refl.coord.connected = graph.refl.coord[graph.reflection.main.index,]
  ad_matrix.reflection.connected = ad_matrix.reflection[graph.reflection.main.index,graph.reflection.main.index]
  
  graph.reflection.connected = graph_from_adjacency_matrix(ad_matrix.reflection.connected, mode=c('undirected'), weighted=TRUE)
  graph_attr(graph.reflection.connected, 'xy') = graph.refl.coord.connected
  graph_attr(graph.reflection.connected, 'sA') = Matrix(ad_matrix.reflection.connected, sparse=TRUE)
  
  graph.refl = graph.reflection.connected
  
  elements = list()
  m = signal
  
  # Method : KnnAvg || Same
  if (method == 'KnnAvg') {
    for (i in 1:K) {
      s0 = m
      s1 = m
      mu = m
      
      j = 1
      while (1) {
        print(paste('j =', j))
        
        get_refl_signal_KnnAvg = function(graph, signal) {
          
          tmp = c()
          for (idx in left_reflection_idx) {
            tmp = c(tmp, neighbor_average(graph, idx, signal))
          }
          signal.left_reflection = tmp
          
          tmp = c()
          for (idx in right_reflection_idx) {
            tmp = c(tmp, neighbor_average(graph, idx, signal))
          }
          signal.right_reflection = tmp
          
          tmp = c()
          for (idx in lower_reflection_idx) {
            tmp = c(tmp, neighbor_average(graph, idx, signal))
          }
          signal.lower_reflection = tmp
          
          tmp = c()
          for (idx in upper_reflection_idx) {
            tmp = c(tmp, neighbor_average(graph, idx, signal))
          }
          signal.upper_reflection = tmp
          
          tmp = c()
          for (idx in left_upper_reflection_idx) {
            tmp = c(tmp, neighbor_average(graph, idx, signal))
          }
          signal.left_upper_reflection = tmp
          
          tmp = c()
          for (idx in right_upper_reflection_idx) {
            tmp = c(tmp, neighbor_average(graph, idx, signal))
          }
          signal.right_upper_reflection = tmp
          
          tmp = c()
          for (idx in left_lower_reflection_idx) {
            tmp = c(tmp, neighbor_average(graph, idx, signal))
          }
          signal.left_lower_reflection = tmp
          
          tmp = c()
          for (idx in right_lower_reflection_idx) {
            tmp = c(tmp, neighbor_average(graph, idx, signal))
          }
          signal.right_lower_reflection = tmp
          
          signal.reflection = c(signal.left_reflection, signal.right_reflection, signal.lower_reflection, signal.upper_reflection, 
                                signal.left_upper_reflection, signal.right_upper_reflection, signal.left_lower_reflection, signal.right_lower_reflection)
          
          graph.refl.sig = c(signal, signal.reflection)
          graph.refl.sig.connected = graph.refl.sig[graph.reflection.main.index]
          
          
          rst = list(graph.refl.sig.connected)
          names(rst) = c('graph.refl.signal')
          
          return(rst)
        }
        
        s1.refl = get_refl_signal_KnnAvg(graph, s1)$graph.refl.signal
        
        extrema = get_extrema(graph.refl, s1.refl)
        maxima = extrema$maxima_list
        minima = extrema$minima_list
        
        inter.max.refl = graph_interpolation(graph.refl, s1.refl, maxima)
        inter.min.refl = graph_interpolation(graph.refl, s1.refl, minima)
        
        inter.max = inter.max.refl[1:dim(graph$sA)[1]]
        inter.min = inter.min.refl[1:dim(graph$sA)[1]]
        
        upper.envelope = GFT.smoothing.ebayesthresh(graph, inter.max)
        lower.envelope = GFT.smoothing.ebayesthresh(graph, inter.min)
        
        mu = (upper.envelope + lower.envelope) / 2
        
        if (t(mu) %*% mu < (t(s1) %*% s1)/1000) {
          break
        }
        
        s1 = s1 - mu
        
        j = j + 1
        
      }
      elements = c(elements, list(s1))
      m = s0-s1
    }
    elements = c(elements, list(m))
  } else if (method == 'Same') {
    for (i in 1:K) {
      s0 = m
      s1 = m
      mu = m
      
      j = 1
      while (1) {
        print(paste('j =', j))
        
        get_refl_signal_same = function(graph, signal) {
          
          signal.left_reflection = signal[left_reflection_idx]
          signal.right_reflection = signal[right_reflection_idx]
          signal.lower_reflection = signal[lower_reflection_idx]
          signal.upper_reflection = signal[upper_reflection_idx]
          signal.left_upper_reflection = signal[left_upper_reflection_idx]
          signal.right_upper_reflection = signal[right_upper_reflection_idx]
          signal.left_lower_reflection = signal[left_lower_reflection_idx]
          signal.right_lower_reflection = signal[right_lower_reflection_idx]
          
          signal.reflection = c(signal.left_reflection, signal.right_reflection, signal.lower_reflection, signal.upper_reflection, 
                                signal.left_upper_reflection, signal.right_upper_reflection, signal.left_lower_reflection, signal.right_lower_reflection)
          
          graph.refl.sig = c(signal, signal.reflection)
          graph.refl.sig.connected = graph.refl.sig[graph.reflection.main.index]
          
          rst = list(graph.refl.sig.connected)
          names(rst) = c('graph.refl.signal')
          
          return(rst)
        }
        
        s1.refl = get_refl_signal_same(graph, s1)$graph.refl.signal
        
        extrema = get_extrema(graph.refl, s1.refl)
        maxima = extrema$maxima_list
        minima = extrema$minima_list
        
        inter.max.refl = graph_interpolation(graph.refl, s1.refl, maxima)
        inter.min.refl = graph_interpolation(graph.refl, s1.refl, minima)
        
        inter.max = inter.max.refl[1:dim(graph$sA)[1]]
        inter.min = inter.min.refl[1:dim(graph$sA)[1]]
        
        upper.envelope = GFT.smoothing.ebayesthresh(graph, inter.max)
        lower.envelope = GFT.smoothing.ebayesthresh(graph, inter.min)
        
        mu = (upper.envelope + lower.envelope) / 2
        
        if (t(mu) %*% mu < (t(s1) %*% s1)/1000) {
          break
        }
        
        s1 = s1 - mu
        
        j = j + 1
        
      }
      elements = c(elements, list(s1))
      m = s0-s1
    }
    elements = c(elements, list(m))
  }
  
  return(elements)
}


# Graph Fourier Decomposition

GFD = function(graph, signal, K) {
  decomposition_rst = list()
  adjacency_matrix = graph$sA
  dimension = dim(adjacency_matrix)[1]
  
  Laplacian_matrix = diag(rowSums(adjacency_matrix)) - adjacency_matrix
  eigen.Laplacian = eigen(Laplacian_matrix)
  eigenvalues = eigen.Laplacian$values
  eigenvectors = eigen.Laplacian$vectors
  
  GFT.signal = t(eigenvectors) %*% signal
  sqrt_periodogram = abs(GFT.signal)
  idx.significant_component = order(sqrt_periodogram)[(dimension-K+1):dimension]
  
  component_sum = 0
  
  for (i in idx.significant_component[order(idx.significant_component)]) {
    component = GFT.signal[i] * eigenvectors[,i]
    component_sum = component_sum + component
    decomposition_rst = c(decomposition_rst, list(component))
  }
  residue = signal - component_sum
  
  decomposition_rst = c(decomposition_rst, list(residue))
  
  return(decomposition_rst)
}