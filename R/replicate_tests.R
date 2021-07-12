# test code for fovea modelling functions - old
# @ unnessecery 
test <- function() {
	data_folder = "/Users/tomas/projects/2016/eye_imaging/30_replicates"
	file_dirs = paste(data_folder, dir(data_folder), "OCT", sep="/")
	pdf("fovea_reps_better.pdf")
	par(mfrow=c(1,2))
	params = NULL
	for(x in 1:length(file_dirs)) {
		files = paste(file_dirs[x], dir(file_dirs[x]), sep="/")
		image_data = dicomInfo(files[1])
		xy_positions = get_fovea_slice_and_position(image_data)
		image_data = dicomInfo(files[2])
		xy_positions = get_fovea_slice_and_position(image_data)

	}
	dev.off()
}

# fovea modelling functions - old
# @ unnessecery
get_params <- function(image_data) {
	fpos = get_fovea_position(image_data)
	slices = get_fovea_slices(image_data, fpos)
	slice = slices$midpoint
	fovea_slice = t(image_data$img[,,slice])
	x_fovea_params = fit_x_fovea(fovea_slice)
	y_fovea_params = fit_y_fovea(fovea_slice, x_fovea_params)
	dip = fit_sholders(fovea_slice, x_fovea_params, y_fovea_params)
	plot_fitted_slice(fovea_slice, x_fovea_params, y_fovea_params, dip)
return(c(unlist(x_fovea_params), unlist(y_fovea_params)))
}

# test code for fovea modelling functions - espliod fitting
# @ nice example of current fitting methods
test.elipse <- function() {
	data_folder = "/Users/tomas/projects/2016/eye_imaging/30_replicates"
	file_dirs = paste(data_folder, dir(data_folder), "OCT", sep="/")
	pdf("fovea_reps_elipsoids.pdf", width=10, height=5)
	par(mfrow=c(1,2))
	params = NULL
	for(x in 1:length(file_dirs)) {
		files = paste(file_dirs[x], dir(file_dirs[x]), sep="/")
		dev_register(files[1])
		dev_register(files[2])
	}
	dev.off()
}


# random code from development
# @ unnessecery 
replicate_example <- function() {
	setwd("~/projects/2016/eye_imaging/30_replicates/")
	data_folder = "/Users/tomas/projects/2016/eye_imaging/30_replicates"
	file_dirs = paste(data_folder, dir(data_folder), "OCT", sep="/")

	par(mfrow=c(1,2))
	params = NULL
	for(x in 1:length(file_dirs)) {
		files = paste(file_dirs[x], dir(file_dirs[x]), sep="/")
		image_data = dicomInfo(files[1])
		#get_fovea_slice_and_positions(image_data)
		p = get_params(image_data)
		k = unlist(strsplit(files[1], "/"))
		params = rbind(params, c(k[8], k[10], p))

		image_data = dicomInfo(files[2])
		#get_fovea_slice_and_positions(image_data)
		p = get_params(image_data)
		k = unlist(strsplit(files[1], "/"))
		params = rbind(params, c(k[8], k[10], p))
		print(x)

	}

	params = as.data.frame(params)
	p_ind = c(dim(image_data$img)[2], dim(image_data$img)[1])
	m1 = as.numeric(as.character(params$midpoint))/p_ind[1]
	x1 = as.numeric(as.character(params$l_shoulder))/p_ind[1]
	x2 = as.numeric(as.character(params$r_shoulder))/p_ind[1]
		
	# this for dynamic x points
	#x1 = y_fovea_params$left_position/p_ind[1]
	#x2 = y_fovea_params$right_position/p_ind[1]
		
	y1 = as.numeric(as.character(params$dip_pos))/p_ind[2]
	y2 = as.numeric(as.character(params$shoulder_r_pos))/p_ind[2]
	ly1 = as.numeric(as.character(params$left_top_sholder))/p_ind[2]
	ry1 = as.numeric(as.character(params$right_top_sholder))/p_ind[2]

	lens = x2-x1
	lens=cbind(lens[seq(1, length(x1), by=2)], lens[seq(2, length(x1), by=2)])

	heis = ((ly1 + ry1)/2)-y1
	heis=cbind(heis[seq(1, length(heis), by=2)], heis[seq(2, length(heis), by=2)])

}