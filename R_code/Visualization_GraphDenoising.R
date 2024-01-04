# Define graph

x = rep(0:29, 30)
y = rep(0:29, rep(30,30))

xy = data.frame(x,y)

g1.wei_ad_mat = weighted_ad_matrix(900, 1, xy) #get weighted adjacency matrix of graph g, based on distance between vertices

g1 = graph_from_adjacency_matrix(g1.wei_ad_mat, mode=c('undirected'), weighted=TRUE) #get graph from adjacency matrix
graph_attr(g1, 'xy') = xy
graph_attr(g1, 'sA') = Matrix(g1.wei_ad_mat, sparse=TRUE)

nv.1 = dim(g1$sA)[1]

# Define Signal

g1.sg1 = 1*sin(2*pi*2*seq(0, 1.25, length=30) + (-0.5*pi)) # Low-frequency component
g1.sg1 = rep(g1.sg1, 30)
g1.sg2 = 1*sin(2*pi*5*seq(0, 1.25, length=30) + (-1.25*pi)) # High-frequency component
g1.sg2 = rep(g1.sg2, 30)

SNR = 7
Var.noise = var(g1.sg1 + g1.sg2) / SNR

set.seed(68)
g1.noise = rnorm(900, 0, sqrt(Var.noise)) 

g1.sig1 = g1.sg1 + g1.sg2 + g1.noise

# Compare graph interpolation and graph denoising
maxima.1 = get_extrema(g1, g1.sig1)$maxima_list
minima.1 = get_extrema(g1, g1.sig1)$minima_list

upper_envelope = graph_interpolation(g1, g1.sig1, maxima.1)
lower_envelope = graph_interpolation(g1, g1.sig1, minima.1)

# Visualization of empirical Bayes thresholding

adjacency_matrix = g1$sA
L = diag(rowSums(adjacency_matrix)) - adjacency_matrix
eigenvectors = eigen(L)$vectors

GFT = t(eigenvectors) %*% upper_envelope # Calculate graph Fourier transforms
GFT.Thresh = ebayesthresh(GFT, verbose = TRUE, threshrule = 'soft') # Obtain graph Fourier coefficients thresholded by empirical Bayes thresholding

plot(GFT, type = 'h', tck = FALSE, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
plot(GFT.Thresh$muhat, type = 'h', tck = FALSE, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')

# 3-dimensional visualization of graph interpolation

tmp.1 = data.frame(x, y, upper_envelope, lower_envelope, idx.extrema.1)
plot_signal1(g1, idx.extrema.1, idx.extrema.1, size=4, custom_colours=custom_colours1)

fig.1.upper = plot_ly(z = ~matrix(upper_envelope, nrow=30, byrow=TRUE), type='surface', showscale=FALSE, showlegend=FALSE) %>% layout(scene = list(xaxis = list(title = '', showticklabels=FALSE),
                                                                                                                                                   yaxis = list(title = '', showticklabels=FALSE),
                                                                                                                                                   zaxis = list(title = '', showticklabels=FALSE)),
                                                                                                                                      legend = list())
fig.1.upper = add_markers(fig.1.upper, type='scatter', x = tmp.1$x[maxima.1], y = tmp.1$y[maxima.1], z = g1.sig1[maxima.1])
fig.1.upper

fig.1.lower = plot_ly(z = ~matrix(lower_envelope, nrow=30, byrow=TRUE), type='surface', showscale=FALSE, showlegend=FALSE) %>% layout(scene = list(xaxis = list(title = '', showticklabels=FALSE),
                                                                                                                                                   yaxis = list(title = '', showticklabels=FALSE),
                                                                                                                                                   zaxis = list(title = '', showticklabels=FALSE)),
                                                                                                                                      legend = list())
fig.1.lower = add_markers(fig.1.lower, type='scatter', x = tmp.1$x[minima.1], y = tmp.1$y[minima.1], z = g1.sig1[minima.1])
fig.1.lower


# 3-dimensional visualization of graph denoising

upper_envelope_smoothing = GFT.smoothing.ebayesthresh(g1, upper_envelope)
lower_envelope_smoothing = GFT.smoothing.ebayesthresh(g1, lower_envelope)

tmp.2 = data.frame(x, y, upper_envelope_smoothing, lower_envelope_smoothing, idx.extrema.1)
plot_signal1(g1, idx.extrema.1, idx.extrema.1, size=4, custom_colours=custom_colours1)

fig.2.upper = plot_ly(z = ~matrix(upper_envelope_smoothing, nrow=30, byrow=TRUE), type='surface', showscale=FALSE, showlegend=FALSE) %>% layout(scene = list(xaxis = list(title = '', showticklabels=FALSE),
                                                                                                                                                             yaxis = list(title = '', showticklabels=FALSE),
                                                                                                                                                             zaxis = list(title = '', showticklabels=FALSE)),
                                                                                                                                                legend = list())
fig.2.upper = add_markers(fig.2.upper, type='scatter', x = tmp.2$x[maxima.1], y = tmp.2$y[maxima.1], z = g1.sig1[maxima.1])
fig.2.upper

fig.2.lower = plot_ly(z = ~matrix(lower_envelope_smoothing, nrow=30, byrow=TRUE), type='surface', showscale=FALSE, showlegend=FALSE) %>% layout(scene = list(xaxis = list(title = '', showticklabels=FALSE),
                                                                                                                                                             yaxis = list(title = '', showticklabels=FALSE),
                                                                                                                                                             zaxis = list(title = '', showticklabels=FALSE)),
                                                                                                                                                legend = list())
fig.2.lower = add_markers(fig.2.lower, type='scatter', x = tmp.2$x[minima.1], y = tmp.2$y[minima.1], z = g1.sig1[minima.1])
fig.2.lower


