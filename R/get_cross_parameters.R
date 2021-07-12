
#############
get_cross_params <- function(z_mean){
  pt_cross_data = NULL
  xparams = peak_trove(apply(z_mean, 2, mean))
  pt_cross_data = rbind(pt_cross_data, unlist(xparams))
  yparams = peak_trove(apply(z_mean, 1, mean))
  pt_cross_data = rbind(pt_cross_data, unlist(yparams))
  diag1=diag(z_mean)
  diag1params=peak_trove(diag1)
  pt_cross_data = rbind(pt_cross_data, unlist(diag1params))
  flipped_matrix <- z_mean[,c(ncol(z_mean):1)]
  diag2=diag(flipped_matrix)
  diag2params=peak_trove(diag2)
  pt_cross_data = rbind(pt_cross_data, unlist(diag2params))
  return(as.data.frame(pt_cross_data))
}

#########################

threed_cross_data <- function(pt_cross_data){
  threed_cross_data = NULL
  threed_cross_data = rbind(threed_cross_data, xdataxcoord = c(0, pt_cross_data$depths1[[1]], pt_cross_data$depths2[[1]]))
  threed_cross_data = rbind(threed_cross_data, xdataycoord = c(0, -(pt_cross_data$radius1[[1]]), pt_cross_data$radius2[[1]]))
  threed_cross_data = rbind(threed_cross_data, xdatazcoord = c(0,0,0))
  threed_cross_data = rbind(threed_cross_data, ydataxcoord = c(0, pt_cross_data$depths1[[2]], pt_cross_data$depths2[[2]]))
  threed_cross_data = rbind(threed_cross_data, ydataycoord = c(0,0,0))
  threed_cross_data = rbind(threed_cross_data, ydatazcoord = c(0, -(pt_cross_data$radius1[[2]]), pt_cross_data$radius2[[2]]))
  angle = pi/4
  threed_cross_data = rbind(threed_cross_data, posdiagxcoord = c(0, pt_cross_data$depths1[[3]], pt_cross_data$depths2[[3]]))
 #threed_cross_data = rbind(threed_cross_data, posdiagycoord = c(0, ((cos(angle))*pt_cross_data$radius1[[3]])*-1, (cos(angle))*pt_cross_data$radius2[[3]])
 # threed_cross_data = rbind(threed_cross_data, posdiagzcoord = c(0, (cos(angle))*pt_cross_data$radius1[[3]], ((cos(angle))*pt_cross_data$radius2[[3]])))
 # threed_cross_data = rbind(threed_cross_data, negdiagxcoord = c(0, pt_cross_data$depths1[[4]], pt_cross_data$depths2[[4]]))
 # threed_cross_data = rbind(threed_cross_data, negdiagycoord = c(0, (cos(angle))*pt_cross_data$radius1[[4]], -((cos(angle))*pt_cross_data$radius2[[4]]))
 # threed_cross_data = rbind(threed_cross_data, negdiagzcoord = c(0, (cos(angle))*pt_cross_data$radius1[[4]], -((cos(angle))*pt_cross_data$radius2[[4]]))
  
}

threed_cross_plot <- function(threed_cross_data){
  plot_3d_cross <- data.frame(threed_cross_data$xdataxcoord, threed_cross_data$xdataycoord, threed_cross_data$xdataxcoord)
  colnames(plot_3d_cross) <- c("x", "y", "x")
  ydata = data.frame(threed_cross_data$ydataxcoord, threed_cross_data$ydataycoord, threed_cross_data$ydatazcoord)
  colnames(ydata) <- c("x","y", "z")
  plot_3d_cross = rbind(plot_3d_cross, ydata)
  posdiag = data.frame(threed_cross_data$posdiagxcoord, threed_cross_data$posdiagycoord, threed_cross_data$posdiagzcoord)
  colnames(posdiag) <- c("x","y", "z")
  plot_3d_cross = rbind(plot_3d_cross, posdiag)
  negdiag = data.frame(threed_cross_data$negdiagxcoord, threed_cross_data$negdiagycoord, threed_cross_data$negdiagzcoord)
  colnames(negdiag) <- c("x", "y", "z")
  plot_3d_cross = rbind(plot_3d_cross, negdiag)
  pdf("./3dcross.pdf")
  scatterplot3d(plot_3d_cross$x, plot_3d_cross$y, plot_3d_cross$z)
  dev.off
}


cross_fovea_shoulder_plot <- function(threed_cross_data){
  plot_3d_cross <- data.frame(threed_cross_data$xdataxcoord, threed_cross_data$xdataycoord, threed_cross_data$xdataxcoord)
  colnames(plot_3d_cross) <- c("x", "y", "x")
  ydata = data.frame(threed_cross_data$ydataxcoord, threed_cross_data$ydataycoord, threed_cross_data$ydatazcoord)
  colnames(ydata) <- c("x","y", "z")
  plot_3d_cross = rbind(plot_3d_cross, ydata)
  posdiag = data.frame(threed_cross_data$posdiagxcoord, threed_cross_data$posdiagycoord, threed_cross_data$posdiagzcoord)
  colnames(posdiag) <- c("x","y", "z")
  plot_3d_cross = rbind(plot_3d_cross, posdiag)
  negdiag = data.frame(threed_cross_data$negdiagxcoord, threed_cross_data$negdiagycoord, threed_cross_data$negdiagzcoord)
  colnames(negdiag) <- c("x", "y", "z")
  plot_3d_cross = rbind(plot_3d_cross, negdiag)
  pdf("./cross_fovea_shoulder.pdf")
  plot(plot_3d_cross$y, plot_3d_cross$z)
  dev.off()
}

#################################################################################
#z_mean = get_z_mean(image_data$img)
#image(z_mean)
#z_mean=make_square(z_mean)
#image(z_mean)
