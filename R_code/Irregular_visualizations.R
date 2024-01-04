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

SNR = 5
Var.noise = var(g2.sg2 + g2.sg1) / SNR

set.seed(452)
g2.noise = rnorm(nv.2, 0, sqrt(Var.noise)) 

g2.sig1 = g2.sg1 + g2.sg2 + g2.noise

# Plot graph signal

plot_signal1(g2, g2.sig1, g2.sig1, size=5, custom_colours=custom_colours1)

# Graph EMD

g2.GEMD = GEMD(g2, g2.sig1, 3)

plot_signal1(g2, g2.sig1, g2.GEMD[[1]], size=5, custom_colours=custom_colours1) # Extracted noise
plot_signal1(g2, g2.sig1, g2.noise, size=5, custom_colours=custom_colours1) # True noise

plot_signal1(g2, g2.sig1, g2.GEMD[[2]], size=5, custom_colours=custom_colours1) # 1st IMF
plot_signal1(g2, g2.sig1, g2.sg2, size=5, custom_colours=custom_colours1) # High-frequency component

plot_signal1(g2, g2.sig1, g2.GEMD[[3]], size=5, custom_colours=custom_colours1) # 2nd IMF
plot_signal1(g2, g2.sig1, g2.sg1, size=5, custom_colours=custom_colours1) # Low-frequency component

plot_signal1(g2, g2.sig1, g2.GEMD[[4]], size=5, custom_colours=custom_colours1) # Residue
plot_signal1(g2, g2.sig1, rep(0,nv.2), size=5, custom_colours=custom_colours1) # True residue

# Calculate MSE of each component from GEMD

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

plot_signal1(g2, g2.sig1, g2.GEMD.refl_KnnAvg.GFT[[1]], size=5, custom_colours=custom_colours1) # Extracted noise
plot_signal1(g2, g2.sig1, g2.noise, size=5, custom_colours=custom_colours1) # True noise

plot_signal1(g2, g2.sig1, g2.GEMD.refl_KnnAvg.GFT[[2]], size=5, custom_colours=custom_colours1) # 1st IMF
plot_signal1(g2, g2.sig1, g2.sg2, size=5, custom_colours=custom_colours1) # High-frequency component

plot_signal1(g2, g2.sig1, g2.GEMD.refl_KnnAvg.GFT[[3]], size=5, custom_colours=custom_colours1) # 2nd IMF
plot_signal1(g2, g2.sig1, g2.sg1, size=5, custom_colours=custom_colours1) # Low-frequency component

plot_signal1(g2, g2.sig1, g2.GEMD.refl_KnnAvg.GFT[[4]], size=5, custom_colours=custom_colours1) # Residue
plot_signal1(g2, g2.sig1, rep(0,nv.2), size=5, custom_colours=custom_colours1) # True residue

# Calculate MSE of each component from SGEMD

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

plot_signal1(g2, g2.sig1, g2.GFD[[3]], size=5, custom_colours=custom_colours1) # Extracted Noise (remaining part after extracting high and low frequency components)
plot_signal1(g2, g2.sig1, g2.noise, size=5, custom_colours=custom_colours1) # True noise

plot_signal1(g2, g2.sig1, g2.GFD[[1]], size=5, custom_colours=custom_colours1) # High-frequency component
plot_signal1(g2, g2.sig1, g2.sg2, size=5, custom_colours=custom_colours1) # True high-frequency component

plot_signal1(g2, g2.sig1, g2.GFD[[2]], size=5, custom_colours=custom_colours1) # Low-frequency component
plot_signal1(g2, g2.sig1, g2.sg1, size=5, custom_colours=custom_colours1) # True low-frequency component

# Calculate MSE of each component from GFD

g2.MSE.GFD.denoised = sum((g2.GFD[[3]] - g2.noise)^2) / nv.2
g2.MSE.GFD.1st = sum((g2.GFD[[1]] - g2.sg2)^2) / nv.2
g2.MSE.GFD.2nd = sum((g2.GFD[[2]] - g2.sg1)^2) / nv.2

g2.GFD.MSE = c(g2.MSE.GFD.denoised, 
               g2.MSE.GFD.1st, 
               g2.MSE.GFD.2nd)


# Obtain and visualize denoised signals for each method

true.denoised = g2.sg2 + g2.sg1 # True
g2.GEMD.denoised = g2.GEMD[[2]] + g2.GEMD[[3]] + g2.GEMD[[4]] # Graph EMD
g2.GEMD.GFT.KnnAvg.denoised = g2.GEMD.refl_KnnAvg.GFT[[2]] + g2.GEMD.refl_KnnAvg.GFT[[3]] + g2.GEMD.refl_KnnAvg.GFT[[4]] #Statistical graph EMD
g2.GFD.denoised = g2.GFD[[1]] + g2.GFD[[2]] # GFD

plot_signal1(g2, g2.sig1, true.denoised, size=5, custom_colours=custom_colours1) # True
plot_signal1(g2, g2.sig1, g2.GEMD.denoised, size=5, custom_colours=custom_colours1) # Graph EMD
plot_signal1(g2, g2.sig1, g2.GEMD.GFT.KnnAvg.denoised, size=5, custom_colours=custom_colours1) # Statistical graph EMD
plot_signal1(g2, g2.sig1, g2.GFD.denoised, size=5, custom_colours=custom_colours1) # GFD

# MSE for this case

rows = c('denoised_part', '1st_IMF', '2nd_IMF')
g2.MSE = data.frame(rows, g2.MSE.GEMD[c(1,2,3)], g2.MSE.GEMD.refl_KnnAvg.GFT[c(1,2,3)], g2.GFD.MSE)
colnames(g2.MSE) = c('rows', 'GEMD', 'SGEMD', 'GFD')

g2.MSE
