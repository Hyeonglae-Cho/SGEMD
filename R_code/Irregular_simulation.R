# Define graph

x = rep(0:29, 30) 
y = rep(0:29, rep(30,30))

xy = data.frame(x,y)

set.seed(315)
xy$x = (xy$x + runif(dim(xy)[1], -0.5,0.5) + 0.5) / 30
xy$y = (xy$y + runif(dim(xy)[1], -0.5,0.5) + 0.5) / 30

g2.wei_ad_mat = weighted_ad_matrix_k(900, 7, xy) # Obtain weighted adjacency matrix of irregular graph g, by connecting 7 nearest vertices for each vertex

g2 = graph_from_adjacency_matrix(g2.wei_ad_mat, mode=c('undirected'), weighted=TRUE) # Obtain irregular graph from adjacency matrix
graph_attr(g2, 'xy') = xy 
graph_attr(g2, 'sA') = Matrix(g2.wei_ad_mat, sparse=TRUE) 

nv.2 = dim(g2$sA)[1] # Number of vertices

# Define graph signal

g2.sg1 = 1 * sin(2*pi*1*g2$xy$x-0.5*pi) * sin(2*pi*1*g2$xy$y-0.5*pi) # Low-frequency component
g2.sg2 = 1 * sin(2*pi*3*g2$xy$x-0.5*pi) * sin(2*pi*3*g2$xy$y-0.5*pi) # High-frequency component

g2.MSE.SNR = list()
empirical_SNR = list()

for (SNR in c(3, 5, 7)) {
  Var.noise = var(g2.sg1 + g2.sg2) / SNR
  
  set.seed(777)
  indices = sample(1:1000, 30)
  
  g2.MSE = list()
  SNR.simulation = c()
  
  for (i in indices) {
    # Define signal
    
    set.seed(i)
    g2.noise = rnorm(nv.2, 0, sqrt(Var.noise)) 
    SNR.simulation = c(SNR.simulation, var(g2.sg1 + g2.sg2) / var(g2.noise))
    
    g2.sig1 = g2.sg1 + g2.sg2 + g2.noise
    
    # Graph EMD
    
    g2.GEMD = GEMD(g2, g2.sig1, 3)
    
    g2.MSE.GEMD.denoised = sum((g2.GEMD[[1]] - g2.noise)^2) / nv.2
    g2.MSE.GEMD.1st = sum((g2.GEMD[[2]] - g2.sg2)^2) / nv.2
    g2.MSE.GEMD.2nd = sum((g2.GEMD[[3]] - g2.sg1)^2) / nv.2
    g2.MSE.GEMD.res = sum((g2.GEMD[[4]] - rep(0,nv.2))^2) / nv.2
    
    g2.MSE.GEMD = c(g2.MSE.GEMD.denoised,
                    g2.MSE.GEMD.1st,
                    g2.MSE.GEMD.2nd,
                    g2.MSE.GEMD.res)

    # Statistical graph EMD
    
    g2.GEMD.refl_KnnAvg.GFT = GEMD.refl.GFT.ebayesthresh(g2, g2.sig1, 3, 'KnnAvg', 0.05, 10, connection = 'neighbor')
    
    g2.MSE.GEMD.refl_KnnAvg.GFT.denoised = sum((g2.GEMD.refl_KnnAvg.GFT[[1]] - g2.noise)^2) / nv.2
    g2.MSE.GEMD.refl_KnnAvg.GFT.1st = sum((g2.GEMD.refl_KnnAvg.GFT[[2]] - g2.sg2)^2) / nv.2
    g2.MSE.GEMD.refl_KnnAvg.GFT.2nd = sum((g2.GEMD.refl_KnnAvg.GFT[[3]] - g2.sg1)^2) / nv.2
    g2.MSE.GEMD.refl_KnnAvg.GFT.res = sum((g2.GEMD.refl_KnnAvg.GFT[[4]] - rep(0,nv.2))^2) / nv.2
    
    g2.MSE.GEMD.refl_KnnAvg.GFT = c(g2.MSE.GEMD.refl_KnnAvg.GFT.denoised,
                                    g2.MSE.GEMD.refl_KnnAvg.GFT.1st,
                                    g2.MSE.GEMD.refl_KnnAvg.GFT.2nd,
                                    g2.MSE.GEMD.refl_KnnAvg.GFT.res)
    
    # Graph Fourier decomposition
    
    g2.GFD = GFD(g2, g2.sig1, 2)
    
    g2.MSE.GFD.denoised = sum((g2.GFD[[3]] - g2.noise)^2) / nv.2
    g2.MSE.GFD.1st = sum((g2.GFD[[1]] - g2.sg2)^2) / nv.2
    g2.MSE.GFD.2nd = sum((g2.GFD[[2]] - g2.sg1)^2) / nv.2
    
    g2.MSE.GFD = c(g2.MSE.GFD.denoised, 
                   g2.MSE.GFD.1st,
                   g2.MSE.GFD.2nd)
    
    g2.MSE = c(g2.MSE, list(list(g2.MSE.GEMD[c(1,2,3)], g2.MSE.GEMD.refl_KnnAvg.GFT[c(1,2,3)], g2.MSE.GFD)))
  }
  g2.MSE.SNR = c(g2.MSE.SNR, list(g2.MSE))
  empirical_SNR = c(empirical_SNR, list(SNR.simulation))
}

# Calculate average and standard error of Mean Squared Errors (MSEs) from 30 replications

get_df = function(lst, methods) {
  
  tmp = data.frame()
  
  for (i in 1:30) {
    tmp[1:3, i] = lst[[i]][[methods]]
  }
  
  return(tmp)
}

rows = c('denoised_part', '1st_IMF', '2nd_IMF')

g2.MSE.simulation.mean.1 = data.frame(rows)
g2.MSE.simulation.sd.1 = data.frame(rows)

g2.MSE.simulation.mean.2 = data.frame(rows)
g2.MSE.simulation.sd.2 = data.frame(rows)

g2.MSE.simulation.mean.3 = data.frame(rows)
g2.MSE.simulation.sd.3 = data.frame(rows)

g2.MSE.simulation.mean = list(g2.MSE.simulation.mean.1, g2.MSE.simulation.mean.2, g2.MSE.simulation.mean.3)
g2.MSE.simulation.sd = list(g2.MSE.simulation.sd.1, g2.MSE.simulation.sd.2, g2.MSE.simulation.sd.3)

for (k in 1:3) {
  
  for (i in 1:3) {
    df.tmp = get_df(g2.MSE.SNR[[k]], i)
    
    for (j in 1:3){
      avg = mean(as.numeric(df.tmp[j,])) ; sd = sqrt(var(as.numeric(df.tmp[j,])))
      
      g2.MSE.simulation.mean[[k]][j, i+1] = avg
      g2.MSE.simulation.sd[[k]][j, i+1] = sd
    }
    
  }
}

for (k in 1:3) {
  colnames(g2.MSE.simulation.mean[[k]]) = c('rows', 'GEMD', 'SGEMD', 'GFD')
  colnames(g2.MSE.simulation.sd[[k]]) = c('rows', 'GEMD', 'SGEMD', 'GFD')
}

g2.MSE.final.1 = data.frame(rows)
g2.MSE.final.2 = data.frame(rows)
g2.MSE.final.3 = data.frame(rows)
g2.MSE.final = list(g2.MSE.final.1, g2.MSE.final.2, g2.MSE.final.3)

for (k in 1:3) {
  
  for (i in 2:4) {
    g2.MSE.final[[k]][, i] = paste(round(g2.MSE.simulation.mean[[k]][, i], digits=3), ' (', round(g2.MSE.simulation.sd[[k]][, i], digits=3), ')', sep='')
  }
  
}

for (k in 1:3) {
  colnames(g2.MSE.final[[k]]) = c('rows', 'GEMD', 'SGEMD', 'GFD')
}

g2.MSE.final # Average and standard error of the MSEs from 30 replications for three level of SNRs (3, 5, and 7)
