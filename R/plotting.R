# plotting function creating a pdf of old fitting method across the foveal dip slices
# @ visulisation
create_fit_across_foveal_dip_pdf <- function(slices, outname, cut = 0.5) {
	pdf(outname, onefile=T)
	for(slice in slices$l_shoulder:slices$r_shoulder) {
		fovea_slice = t(image_data$img[,,slice])
		y_fovea_params = fit_y_fovea(fovea_slice, x_fovea_params)
		dip = fit_sholders(fovea_slice, x_fovea_params, y_fovea_params)
		plot_fitted_slice(fovea_slice, x_fovea_params, y_fovea_params, dip)
	}
	dev.off()
}

# plotting function creating a pdf of the predicted max dip slice across a directory of dicom files
# @ visulisation
create_max_foveal_dip_pdf <- function(direct, outname, smooth_factor = 0.2) {
	pdf(outname, onefile=T)
	files = dir(direct)
	for(x in 1:length(files)) {
		image_file = paste(direct, "/", files[x], sep="")
		image_data = dicomInfo(image_file)
		fpos = get_fovea_position(image_data$img)
		slices = get_fovea_slices(image_data$img, fpos[1])
		slice = slices$midpoint
		fovea_slice = t(image_data$img[,,slice])
		x_fovea_params = fit_x_fovea(fovea_slice)
		y_fovea_params = fit_y_fovea(fovea_slice, x_fovea_params)
		dip = fit_sholders(fovea_slice, x_fovea_params, y_fovea_params)
		plot_fitted_slice(fovea_slice, x_fovea_params, y_fovea_params, dip)
	}
	dev.off()
}

# plotting function for a fitted slice - old - not robust
# @ visulisation 
plot_fitted_slice <- function(fovea_slice, x_fovea_params, y_fovea_params, dip) {
	p_ind = dim(fovea_slice)
	m1 = x_fovea_params$midpoint/p_ind[1]
	x1 = x_fovea_params$l_shoulder/p_ind[1]
	x2 = x_fovea_params$r_shoulder/p_ind[1]
	
	# this for dynamic x points
	#x1 = y_fovea_params$left_position/p_ind[1]
	#x2 = y_fovea_params$right_position/p_ind[1]
	
	y1 = y_fovea_params$dip_pos/p_ind[2]
	y2 = y_fovea_params$shoulder_r_pos/p_ind[2]
	ly1 = y_fovea_params$left_top_sholder/p_ind[2]
	ry1 = y_fovea_params$right_top_sholder/p_ind[2]
	ly2 = y1
	ry2 = y1
	coor = c(ly1, ry1)
	if(which.min(coor)==1) {
		ly2 = ly2 - (ry1-ly1)
		ry2 = ry2 + (ry1-ly1)
	} else {
		ly2 = ly2 + (ry1-ly1)
		ry2 = ry2 - (ry1-ly1)
	}

	image(fovea_slice, breaks=c(0,50,200), col=c("black", "green"))
	#segments(x1,ly2,x2,ry2, col="blue")
	segments(x1,y1,x2,y1, col="blue")
	segments(x1,ly1,x2,ry1, col="blue")
	segments(x1,y1,x1,ly1, col="blue")
	segments(x2,y1,x2,ry1, col="blue")
	segments(m1,y1,m1,mean(c(ly1, ry1)), col="red")

	co = NULL
	if(length(dip$y[dip$y<mean(c(ly1, ry1))])/length(dip$y)>0.2) {
		co = "red"
		matplot(dip$x, dip$y, lwd=2, type='l', col=co, add=T)
	}
	if(length(dip$y[dip$y<mean(c(ly1, ry1))])/length(dip$y)>0.5) {
		co = "blue"
		matplot(dip$x, dip$y, lwd=2, type='l', col=co, add=T)
	}
}
