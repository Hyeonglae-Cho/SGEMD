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

SNR = 3
Var.noise = var(g1.sg1 + g1.sg2) / SNR

set.seed(510)
g1.noise = 0.9 * rnorm(900, 0, sqrt(Var.noise/100.81)) + 0.1 * rnorm(900, 0, sqrt(Var.noise*10000/100.81)) # Gaussian mixture noise where sigma1 << sigma2 (=100*sigma1)

var(g1.sg1 + g1.sg2) / var(g1.noise)

g1.sig2 = g1.sg1 + g1.sg2 + g1.noise

# Plot graph signal

plot_signal1(g1, g1.sig2, g1.sig2, size=6, custom_colours=custom_colours1)

# Graph EMD

start.time = Sys.time()
g1.GEMD = GEMD(g1, g1.sig2, 3)
end.time = Sys.time()
req.time_GEMD = end.time - start.time # Required time for GEMD

plot_signal1(g1, g1.sig2, g1.GEMD[[1]], size=6, custom_colours=custom_colours1) # Extracted noise
plot_signal1(g1, g1.sig2, g1.noise, size=6, custom_colours=custom_colours1) # True noise

plot_signal1(g1, g1.sig2, g1.GEMD[[2]], size=6, custom_colours=custom_colours1) # 1st IMF
plot_signal1(g1, g1.sig2, g1.sg2, size=6, custom_colours=custom_colours1) # High-frequency component

plot_signal1(g1, g1.sig2, g1.GEMD[[3]], size=6, custom_colours=custom_colours1) # 2nd IMF
plot_signal1(g1, g1.sig2, g1.sg1, size=6, custom_colours=custom_colours1) # Low-frequency component

plot_signal1(g1, g1.sig2, g1.GEMD[[4]], size=6, custom_colours=custom_colours1) # Residue
plot_signal1(g1, g1.sig2, rep(0,nv.1), size=6, custom_colours=custom_colours1) # True residue

# Calculate MSE of each component from GEMD

g1.MSE.GEMD.1st = sum((g1.GEMD[[1]] - g1.noise)^2) / nv.1
g1.MSE.GEMD.2nd = sum((g1.GEMD[[2]] - g1.sg2)^2) / nv.1
g1.MSE.GEMD.3rd = sum((g1.GEMD[[3]] - g1.sg1)^2) / nv.1
g1.MSE.GEMD.res = sum((g1.GEMD[[4]] - rep(0,nv.1))^2) / nv.1

g1.MSE.GEMD = c(g1.MSE.GEMD.1st,
                g1.MSE.GEMD.2nd,
                g1.MSE.GEMD.3rd,
                g1.MSE.GEMD.res)

# Statistical graph EMD

start.time = Sys.time()
g1.GEMD.refl_KnnAvg.GFT = GEMD.refl.GFT.ebayesthresh(g1, g1.sig2, 3, 'KnnAvg', 1, 6)
end.time = Sys.time()
req.time_SGEMD = end.time - start.time # Required time for SGEMD

plot_signal1(g1, g1.sig2, g1.GEMD.refl_KnnAvg.GFT[[1]], size=6, custom_colours=custom_colours1) # Extracted noise
plot_signal1(g1, g1.sig2, g1.noise, size=6, custom_colours=custom_colours1) # True noise

plot_signal1(g1, g1.sig2, g1.GEMD.refl_KnnAvg.GFT[[2]], size=6, custom_colours=custom_colours1) # 1st IMF
plot_signal1(g1, g1.sig2, g1.sg2, size=6, custom_colours=custom_colours1) # High-frequency component

plot_signal1(g1, g1.sig2, g1.GEMD.refl_KnnAvg.GFT[[3]], size=6, custom_colours=custom_colours1) # 2nd IMF
plot_signal1(g1, g1.sig2, g1.sg1, size=6, custom_colours=custom_colours1) # Low-frequency component

plot_signal1(g1, g1.sig2, g1.GEMD.refl_KnnAvg.GFT[[4]], size=6, custom_colours=custom_colours1) # Residue
plot_signal1(g1, g1.sig2, rep(0,nv.1), size=6, custom_colours=custom_colours1) # True residue

# Calculate MSE of each component from SGEMD

g1.MSE.GEMD.refl_KnnAvg.GFT.1st = sum((g1.GEMD.refl_KnnAvg.GFT[[1]] - g1.noise)^2) / nv.1
g1.MSE.GEMD.refl_KnnAvg.GFT.2nd = sum((g1.GEMD.refl_KnnAvg.GFT[[2]] - g1.sg2)^2) / nv.1
g1.MSE.GEMD.refl_KnnAvg.GFT.3rd = sum((g1.GEMD.refl_KnnAvg.GFT[[3]] - g1.sg1)^2) / nv.1
g1.MSE.GEMD.refl_KnnAvg.GFT.res = sum((g1.GEMD.refl_KnnAvg.GFT[[4]] - rep(0,nv.1))^2) / nv.1

g1.MSE.GEMD.refl_KnnAvg.GFT = c(g1.MSE.GEMD.refl_KnnAvg.GFT.1st,
                                g1.MSE.GEMD.refl_KnnAvg.GFT.2nd,
                                g1.MSE.GEMD.refl_KnnAvg.GFT.3rd,
                                g1.MSE.GEMD.refl_KnnAvg.GFT.res)

# Graph Fourier decomposition

g1.GFD = GFD(g1, g1.sig2, 2)

plot_signal1(g1, g1.sig2, g1.GFD[[3]], size = 6, custom_colours = custom_colours1) # Extracted Noise (remaining part after extracting high- and low-frequency components)
plot_signal1(g1, g1.sig2, g1.noise, size = 6, custom_colours = custom_colours1) # True noise

plot_signal1(g1, g1.sig2, g1.GFD[[1]], size = 6, custom_colours = custom_colours1) # High-frequency component
plot_signal1(g1, g1.sig2, g1.sg2, size = 6, custom_colours = custom_colours1) # True high-frequency component

plot_signal1(g1, g1.sig2, g1.GFD[[2]], size = 6, custom_colours = custom_colours1) # Low-frequency component
plot_signal1(g1, g1.sig2, g1.sg1, size = 6, custom_colours = custom_colours1) # True low-frequency component

# Calculate MSE of each component from GFD

g1.MSE.GFD.1st = sum((g1.GFD[[3]] - g1.noise)^2) / nv.1
g1.MSE.GFD.2nd = sum((g1.GFD[[1]] - g1.sg2)^2) / nv.1
g1.MSE.GFD.3rd = sum((g1.GFD[[2]] - g1.sg1)^2) / nv.1

g1.MSE.GFD = c(g1.MSE.GFD.1st, 
               g1.MSE.GFD.2nd, 
               g1.MSE.GFD.3rd)


# Obtain and visualize denoised signals for each method

true.denoised = g1.sg2 + g1.sg1 # True
g1.GEMD.denoised = g1.GEMD[[2]] + g1.GEMD[[3]] + g1.GEMD[[4]] # Graph EMD
g1.GEMD.GFT.KnnAvg.denoised = g1.GEMD.refl_KnnAvg.GFT[[2]] + g1.GEMD.refl_KnnAvg.GFT[[3]] + g1.GEMD.refl_KnnAvg.GFT[[4]] #Statistical graph EMD
g1.GFD.denoised = g1.GFD[[1]] + g1.GFD[[2]] #GFD

plot_signal1(g1, g1.sig2, true.denoised, size=6, custom_colours=custom_colours1) # True
plot_signal1(g1, g1.sig2, g1.GEMD.denoised, size=6, custom_colours=custom_colours1) # Graph EMD
plot_signal1(g1, g1.sig2, g1.GEMD.GFT.KnnAvg.denoised, size = 6, custom_colours = custom_colours1) #Statistical graph EMD
plot_signal1(g1, g1.sig2, g1.GFD.denoised, size = 6, custom_colours = custom_colours1) #GFD

# MSE for this case

rows = c('denoised_part', '1st_IMF', '2nd_IMF')
g1.MSE = data.frame(rows, round(g1.MSE.GEMD[c(1,2,3)], 3), round(g1.MSE.GEMD.refl_KnnAvg.GFT[c(1,2,3)], 3), round(g1.MSE.GFD, 3))
colnames(g1.MSE) = c('rows', 'GEMD', 'SGEMD', 'GFD')

g1.MSE

plot_signal1(g1, g1.noise, g1.noise, size=6, custom_colours=custom_colours1) # Visualizing true noise with noise scale

g1.noise.thresh = rep(0, 900)
g1.noise.thresh[which((g1.noise <= summary(g1.noise)[2] - 1.5 * (summary(g1.noise)[5] - summary(g1.noise)[2]) | g1.noise >= summary(g1.noise)[5] + 1.5 * (summary(g1.noise)[5] - summary(g1.noise)[2])))] = g1.noise[which((g1.noise <= summary(g1.noise)[2] - 1.5 * (summary(g1.noise)[5] - summary(g1.noise)[2]) | g1.noise >= summary(g1.noise)[5] + 1.5 * (summary(g1.noise)[5] - summary(g1.noise)[2])))]
g1.noise.thresh
plot_signal1(g1, g1.noise, g1.noise.thresh, size=6, custom_colours=custom_colours1) # Visualizing outliers
