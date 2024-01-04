# Function 'plot_signal1' is a modified version of the 'plot_signal' function from the 'gasper' package and is licensed under LGPL 2.1.
# For details, see: https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html
# Changes:
# - Modified the function name from 'plot_signal' to 'plot_signal1'.
# - Added new parameters 'original_signal', 'custom_colours' and renamed parameter 'f' to 'signal'.
# - Replaced the original color scale 'scale_colour_distiller' with 'scale_color_gradient'.
# - Adjusted the theme settings.

plot_signal1 <- function(z, original_signal, signal, size=0.75, custom_colours, limits=range(signal)) {
  boundary = range(original_signal)
  b = seq(boundary[1], boundary[2], length=length(custom_colours))
  
  if(is(z$sA, 'sparseMatrix')){
    z$sA <- summary(z$sA)
  }
  x <- z$xy[, 1]
  y <- z$xy[, 2]
  ind_i <- z$sA[, 1]
  ind_j <- z$sA[, 2]
  y1 <- x[ind_j]
  y2 <- y[ind_j]
  df1 <- data.frame(x = x, y = y)
  df2 <- data.frame(x = x[ind_i],
                    y = y[ind_i],
                    y1 = y1, y2 = y2)
  p2 <- ggplot(df1, aes(x, y)) +
    geom_segment(aes(x = x, y = y,
                     xend = y1, yend = y2),
                 color = "gray", data = df2) +
    geom_point(size = size, aes(colour = signal)) +
    scale_color_gradientn(limits = c(boundary[1], boundary[2]), colours = custom_colours, 
                          breaks = b, labels = format(b)) +
    theme(#legend.position = "bottom",
      legend.title=element_blank(),
      legend.text=element_text(size=8),
      legend.key.size = unit(0.8,"line"),
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.title.y=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank(),
      legend.margin = margin(0.0,0.0,0.0,0.0),
      plot.margin = unit(c(0.01,0.01,0.01,0.01), "cm"))
  print(p2)
}


# The 'plot_signal1' below is the original version of the 'plot_signal1' function and is licensed under LGPL 2.1.
#' Plot a Signal on Top of a Given Graph
#'
#' Visualize a signal over a graph.
#'
#'@details
#' This function allows visualization of a graph signal \code{f} superimposed on the structure of a graph defined by \code{z}. It offers an intuitive way to analyze  the behavior of graph signals in the vertex domain.
#'
#'@note If node coordinates \code{xy} are not provided, they will be calculated using spectral methods \code{\link{spectral_coords}}. For large graphs, this can be computationally intensive and may take significant time. Use with caution for large graphs if node coordinates are not supplied.
#'
#' @export plot_signal
#' @importFrom methods is
#' @importFrom Matrix summary
#' @param z A list containing graph data. This list must have the following components:
#'          \itemize{
#'            \item{sA}  An adjacency matrix or a sparse Matrix representation of the graph.
#'            \item{xy}  A matrix or dataframe containing the x and y coordinates of each node in the graph.
#'          }
#' @param f Signal to plot.
#' @param size Numeric. Dot size for nodes. Default is 0.75.
#' @param limits Set colormap limits.
#' @examples
#' f <- rnorm(length(grid1$xy[,1]))
#' plot_signal(grid1, f)
#' @seealso \code{\link{plot_graph}}, \code{\link{spectral_coords}}

plot_signal <- function(z, f, size=0.75, limits=range(f)) {
  if(!"xy" %in% names(z)){
    z$xy <- spectral_coords(z$sA)
  }
  if(is(z$sA, 'sparseMatrix')){
    z$sA <- summary(z$sA)
  }
  x <- z$xy[, 1]
  y <- z$xy[, 2]
  ind_i <- z$sA[, 1]
  ind_j <- z$sA[, 2]
  y1 <- x[ind_j]
  y2 <- y[ind_j]
  df1 <- data.frame(x = x, y = y)
  df2 <- data.frame(x = x[ind_i],
                    y = y[ind_i],
                    y1 = y1, y2 = y2)
  p2 <- ggplot(df1, aes(x, y)) +
    geom_segment(aes(x = x, y = y,
                     xend = y1, yend = y2),
                 color = "gray", data = df2) +
    geom_point(size = size, aes(colour = f)) +
    scale_colour_distiller(palette = "Spectral",
                           limits = limits) +
    theme_void() +
    theme(#legend.position = "bottom",
      legend.text=element_text(size=8),
      legend.key.size = unit(0.8,"line"),
      
      legend.margin = margin(0.0,0.0,0.0,0.0),
      plot.margin = unit(c(0.01,0.01,0.01,0.01), "cm"))
  print(p2)
}
