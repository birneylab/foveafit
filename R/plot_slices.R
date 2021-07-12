#Two functions, one creates a table that extracts parameters for every single row and the other creates a table that extracts parameters for every column and creates an output file 

walk_through_slice_rows <- function(zmean){
  pt_dataxaxis = NULL
  for(x in 1:nrow(z_mean)) {
    p = peak_trove(z_mean[x,])
    pt_dataxaxis = rbind(pt_dataxaxis, unlist(p))
  }
  return(as.data.frame(pt_dataxaxis))
}


walk_through_slice_cols <- function(zmean){
  pt_datayaxis = NULL
  for(x in 1:ncol(z_mean)) {
    p = peak_trove(z_mean[,x])
    pt_datayaxis = rbind(pt_datayaxis, unlist(p))
  }
  return(as.data.frame(pt_datayaxis))
}

plot_slices_data <- function(pt_dataxaxis){
  plotslice_data=NULL
  foveaslice=NULL
  if(nrow(pt_dataxaxis)%%2 == 0){
    foveaslice = nrow(pt_dataxaxis)/2
    zaxis=c((-(foveaslice-1)):foveaslice)
  } else{
    foveaslice = (nrow(pt_dataxaxis)%/%2)+1
    zaxis=c((-(foveaslice-1):(foveaslice-1)))
  }
  plotslice_data = cbind(plotslice_data, zaxis)
  pt_dataxaxis$foveaxcoord <- pt_dataxaxis$position2 - pt_dataxaxis$position2[[foveaslice]]
  plotslice_data <- cbind(plotslice_data, foveaxcoord = pt_dataxaxis$foveaxcoord)
  pt_dataxaxis$foveaycoord <- pt_dataxaxis$fovealoc - pt_dataxaxis$fovealoc[[foveaslice]]
  plotslice_data <- cbind(plotslice_data, foveaycoord = pt_dataxaxis$foveaycoord)
  pt_dataxaxis$lhxcoord <- pt_dataxaxis$foveaxcoord - pt_dataxaxis$radius1
  plotslice_data <- cbind(plotslice_data, lhxcoord = pt_dataxaxis$lhxcoord)
  pt_dataxaxis$lhycoord <- pt_dataxaxis$foveaycoord + pt_dataxaxis$depths1
  plotslice_data <- cbind(plotslice_data, lhycoord = pt_dataxaxis$lhycoord)
  pt_dataxaxis$rhxcoord <- pt_dataxaxis$foveaxcoord + pt_dataxaxis$radius2
  plotslice_data <- cbind(plotslice_data, rhxcoord = pt_dataxaxis$rhxcoord)
  pt_dataxaxis$rhycoord <- pt_dataxaxis$foveaycoord + pt_dataxaxis$depths2
  plotslice_data <- cbind(plotslice_data, rhycoord = pt_dataxaxis$rhycoord)
  plotslice_data <- as.data.frame(plotslice_data)
}

plot_all_slice_points <- function(plot_slices_data){
  plot_all_data = data.frame(plot_slices_data$zaxis, plot_slices_data$foveaxcoord, plot_slices_data$foveaycoord)
  colnames(plot_all_data) <- c("z", "x", "y")
  lhpoint = data.frame(plot_slices_data$zaxis, plot_slices_data$lhxcoord, plot_slices_data$lhycoord)
  colnames(lhpoint) <- c("z", "x", "y")
  plot_all_data = rbind(plot_all_data, lhpoint)
  rhpoint = data.frame(plot_slices_data$zaxis, plot_slices_data$rhxcoord, plot_slices_data$rhycoord)
  colnames(rhpoint) <- c("z", "x", "y")
  plot_all_data = rbind(plot_all_data, rhpoint)
  library(scatterplot3d)
  pdf("./all_points_slices_plot.pdf")
  scatterplot3d(plot_all_data$x, plot_all_data$z, plot_all_data$y)
  dev.off()
}  

plot_one_point <- function(plot_slices_data, data_point=""){
  if(data_point == "trove"){
    pdf("./trove_slices_plot.pdf")
    scatterplot3d(plot_slices_data$foveaxcoord, plot_slices_data$zaxis, plot_slices_data$foveaycoord)
    dev.off()
  } else if(data_point == "left shoulder"){
    pdf("./leftshoulder_slices_plot.pdf")
    scatterplot3d(plot_slices_data$lhxcoord, plot_slices_data$zaxis, plot_slices_data$lhycoord)
    dev.off
  } else if(data_point == "right shoulder"){
    pdf("./rightshoulder_slices_plot.pdf")
    scatterplot3d(plot_slices_data$rhxcoord, plot_slices_data$zaxis, plot_slices_data$rhycoord)
    dev.off
  }
}

plot_fovea_shoulder_ring <- function(plot_slices_data){
  plot_all_data = data.frame(plot_slices_data$zaxis, plot_slices_data$foveaxcoord, plot_slices_data$foveaycoord)
  colnames(plot_all_data) <- c("z", "x", "y")
  lhpoint = data.frame(plot_slices_data$zaxis, plot_slices_data$lhxcoord, plot_slices_data$lhycoord)
  colnames(lhpoint) <- c("z", "x", "y")
  plot_all_data = rbind(plot_all_data, lhpoint)
  rhpoint = data.frame(plot_slices_data$zaxis, plot_slices_data$rhxcoord, plot_slices_data$rhycoord)
  colnames(rhpoint) <- c("z", "x", "y")
  plot_all_data = rbind(plot_all_data, rhpoint)
  pdf("./fovea_shoulder_ring.pdf")
  plot(plot_all_data$x, plot_all_data$z)
  dev.off()
}


#z_mean = get_z_mean(image_data$img)
#pt_dataxaxis = walk_through_slice_rows(zmean)
#pt_datayaxis = walk_through_slice_cols(zmean)

