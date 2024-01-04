# Define graph

x = rep(0:29, 30) 
y = rep(0:29, rep(30,30))

xy = data.frame(x,y)

g1.wei_ad_mat = weighted_ad_matrix(900, 1, xy) # Obtain weighted adjacency matrix of 30 x 30 lattice graph g, based on distance between vertices

g1 = graph_from_adjacency_matrix(g1.wei_ad_mat, mode=c('undirected'), weighted=TRUE) # Obtain 30 x 30 lattice graph from adjacency matrix
graph_attr(g1, 'xy') = xy
graph_attr(g1, 'sA') = Matrix(g1.wei_ad_mat, sparse=TRUE)

nv.1 = dim(g1$sA)[1] # Number of vertices

# Define graph signal

g1.sg1 = 1*sin(2*pi*2*seq(0, 1.25, length=30) + (-0.5*pi)) # Low-frequency component
g1.sg1 = g1.sg1 * seq(from=1, to=0.25, length=30)
g1.sg1 = rep(g1.sg1, 30)
g1.sg2 = 1*sin(2*pi*5*seq(0, 1.25, length=30) + (-1.25*pi)) # High-frequency component
g1.sg2 = g1.sg2 * seq(from=1, to=0.25, length=30)
g1.sg2 = rep(g1.sg2, 30)

set.seed(777)
indices = sample(1:1000, 30)

g1.MSE.SNR = list()
empirical_SNR = list()

for (SNR in c(3,5,7)) {
  Var.noise = var(g1.sg1 + g1.sg2) / SNR
  
  set.seed(777)
  indices = sample(1:1000, 30)
  
  g1.MSE = list()
  SNR.simulation = c()
  
  for (i in indices) {
    # Define signal
    
    set.seed(i)
    g1.noise = 0.9 * rnorm(900, 0, sqrt(Var.noise/100.81)) + 0.1 * rnorm(900, 0, sqrt(Var.noise*10000/100.81)) # Gaussian mixture noise where sigma1 << sigma2 (=100*sigma1)
    SNR.simulation = c(SNR.simulation, var(g1.sg1 + g1.sg2) / var(g1.noise))
    
    g1.sig2 = g1.sg1 + g1.sg2 + g1.noise
    
    # Graph EMD
    
    g1.GEMD = GEMD(g1, g1.sig2, 3)
    
    g1.MSE.GEMD.denoised = sum((g1.GEMD[[1]] - g1.noise)^2) / nv.1
    g1.MSE.GEMD.1st = sum((g1.GEMD[[2]] - g1.sg2)^2) / nv.1
    g1.MSE.GEMD.2nd = sum((g1.GEMD[[3]] - g1.sg1)^2) / nv.1
    g1.MSE.GEMD.res = sum((g1.GEMD[[4]] - rep(0,nv.1))^2) / nv.1
    
    g1.MSE.GEMD = c(g1.MSE.GEMD.denoised,
                    g1.MSE.GEMD.1st,
                    g1.MSE.GEMD.2nd,
                    g1.MSE.GEMD.res)
    
    # Statistical graph EMD
    
    g1.GEMD.refl_KnnAvg.GFT = GEMD.refl.GFT.ebayesthresh(g1, g1.sig2, 3, 'KnnAvg', 1, 6)
    
    g1.MSE.GEMD.refl_KnnAvg.GFT.denoised = sum((g1.GEMD.refl_KnnAvg.GFT[[1]] - g1.noise)^2) / nv.1
    g1.MSE.GEMD.refl_KnnAvg.GFT.1st = sum((g1.GEMD.refl_KnnAvg.GFT[[2]] - g1.sg2)^2) / nv.1
    g1.MSE.GEMD.refl_KnnAvg.GFT.2nd = sum((g1.GEMD.refl_KnnAvg.GFT[[3]] - g1.sg1)^2) / nv.1
    g1.MSE.GEMD.refl_KnnAvg.GFT.res = sum((g1.GEMD.refl_KnnAvg.GFT[[4]] - rep(0,nv.1))^2) / nv.1
    
    g1.MSE.GEMD.refl_KnnAvg.GFT = c(g1.MSE.GEMD.refl_KnnAvg.GFT.denoised,
                                    g1.MSE.GEMD.refl_KnnAvg.GFT.1st,
                                    g1.MSE.GEMD.refl_KnnAvg.GFT.2nd,
                                    g1.MSE.GEMD.refl_KnnAvg.GFT.res)
    
    # Graph Fourier decomposition
    
    g1.GFD = GFD(g1, g1.sig2, 2)
    
    g1.MSE.GFD.denoised = sum((g1.GFD[[3]] - g1.noise)^2) / nv.1
    g1.MSE.GFD.1st = sum((g1.GFD[[1]] - g1.sg2)^2) / nv.1
    g1.MSE.GFD.2nd = sum((g1.GFD[[2]] - g1.sg1)^2) / nv.1
    
    g1.MSE.GFD = c(g1.MSE.GFD.denoised, 
                   g1.MSE.GFD.1st,
                   g1.MSE.GFD.2nd)
    
    
    g1.MSE = c(g1.MSE, list(list(g1.MSE.GEMD[c(1,2,3)], g1.MSE.GEMD.refl_KnnAvg.GFT[c(1,2,3)], g1.MSE.GFD)))
  }
  g1.MSE.SNR = c(g1.MSE.SNR, list(g1.MSE))
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
g1.MSE.simulation.mean.1 = data.frame(rows)
g1.MSE.simulation.sd.1 = data.frame(rows)

g1.MSE.simulation.mean.2 = data.frame(rows)
g1.MSE.simulation.sd.2 = data.frame(rows)

g1.MSE.simulation.mean.3 = data.frame(rows)
g1.MSE.simulation.sd.3 = data.frame(rows)

g1.MSE.simulation.mean = list(g1.MSE.simulation.mean.1, g1.MSE.simulation.mean.2, g1.MSE.simulation.mean.3)
g1.MSE.simulation.sd = list(g1.MSE.simulation.sd.1, g1.MSE.simulation.sd.2, g1.MSE.simulation.sd.3)

for (k in 1:3) {
  
  for (i in 1:3) {
    df.tmp = get_df(g1.MSE.SNR[[k]], i)
    
    for (j in 1:3){
      avg = mean(as.numeric(df.tmp[j,])) ; sd = sqrt(var(as.numeric(df.tmp[j,])))
      
      g1.MSE.simulation.mean[[k]][j, i+1] = avg
      g1.MSE.simulation.sd[[k]][j, i+1] = sd
    }
    
  }
}

for (k in 1:3) {
  colnames(g1.MSE.simulation.mean[[k]]) = c('rows', 'GEMD', 'SGEMD', 'GFD')
  colnames(g1.MSE.simulation.sd[[k]]) = c('rows', 'GEMD', 'SGEMD', 'GFD')
}

g1.MSE.final.1 = data.frame(rows)
g1.MSE.final.2 = data.frame(rows)
g1.MSE.final.3 = data.frame(rows)
g1.MSE.final = list(g1.MSE.final.1, g1.MSE.final.2, g1.MSE.final.3)

for (k in 1:3) {
  
  for (i in 2:4) {
    g1.MSE.final[[k]][, i] = paste(round(g1.MSE.simulation.mean[[k]][, i], digits=3), ' (', round(g1.MSE.simulation.sd[[k]][, i], digits=3), ')', sep='')
  }
  
}

for (k in 1:3) {
  colnames(g1.MSE.final[[k]]) = c('rows', 'GEMD', 'SGEMD', 'GFD')
}

g1.MSE.final # Average and standard error of the MSEs from 30 replications for three level of SNRs (3, 5, and 7)
