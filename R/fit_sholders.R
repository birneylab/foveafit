# fitting across the x axis of a oct slice
# @ fovea modelling
fit_x_fovea <- function(fovea_slice, fac=0.2) {
	x_dist = apply(fovea_slice, 1, mean)
	x_dist = x_dist-median(x_dist)
	l = lowess(x_dist, f=fac)
	up_down = check_up_and_down(l)
	fovea_params = up_down
	if(length(up_down)>1) {
		fovea_params = get_max_peak_trove(up_down, l)
	}
return(fovea_params)
}

# fitting across the y axis of a oct slice
# @ fovea modelling
fit_y_fovea <- function(fovea_slice, x_fovea_params, fac=0.1) {
    x_dist = apply(fovea_slice[x_fovea_params$l_shoulder:x_fovea_params$r_shoulder,], 2, mean)
    x_dist = x_dist - median(x_dist)
    x_dist[which.min(x_dist):1] = -1
    params = NULL
    params$dip_pos = which.max(x_dist)
    x_dist[params$dip_pos:1] = -1
    l = lowess(x_dist, f = fac)
    params$shoulder_r_pos = min(which(l$y < 0 & l$x > params$dip_pos))
    left = fovea_slice[x_fovea_params$l_shoulder:x_fovea_params$midpoint, params$dip_pos:dim(fovea_slice)[2]]
    lef = apply(left, 2, mean)
    lef = lef - median(lef)
    mp = check_up_and_down(lowess(cumsum(lef), f = 0.01))[1]
    x1_adjust = check_up_and_down(lowess(cumsum(left[, mp]), f = 0.01))
    params$left_position = x_fovea_params$l_shoulder
    if (!is.null(x1_adjust)) {
        params$left_position = params$left_position + max(x1_adjust)
    }
    x_fovea_params$l_shoulder = params$left_position
    params$left_top_sholder = params$dip_pos + mp
    right = fovea_slice[x_fovea_params$midpoint:x_fovea_params$r_shoulder, params$dip_pos:dim(fovea_slice)[2]]
    rig = apply(right, 2, mean)
    rig = rig - median(rig)
    mp = check_up_and_down(lowess(cumsum(rig), f = 0.01))[1]
    x2_adjust = check_up_and_down(lowess(cumsum(right[, mp]), f = 0.01))
    params$right_position = x_fovea_params$r_shoulder
    if (!is.null(x2_adjust)) {
        params$right_position = params$right_position - min(x2_adjust)
    }
    x_fovea_params$r_shoulder = params$right_position
    params$right_top_sholder = params$dip_pos + mp
return(params)
}

# function to fit the dip and sholders of the fovea - old
# @ fovea modelling
fit_sholders <- function(fovea_slice, x_fovea_params, y_fovea_params) {
	#max_y = max(c(y_fovea_params$left_top_sholder, y_fovea_params$right_top_sholder))
	#zoom_fovea = fovea_slice[x_fovea_params$l_shoulder:x_fovea_params$r_shoulder,y_fovea_params$dip_pos:max_y]
	zoom_fovea = fovea_slice[x_fovea_params$l_shoulder:x_fovea_params$r_shoulder,y_fovea_params$dip_pos:y_fovea_params$shoulder_r_pos]
	y_dist = apply(zoom_fovea, 1, mean)
	l = lowess(y_dist, f=0.1)
	#image(zoom_fovea)
	#matplot(l$x/dim(zoom_fovea)[1], l$y/dim(zoom_fovea)[2], type='l', col="blue", add=T)

	#image(fovea_slice)
	scales = dim(zoom_fovea)/dim(fovea_slice)
	x_add = x_fovea_params$l_shoulder/dim(fovea_slice)[1]
	y_add = y_fovea_params$dip_pos/ dim(fovea_slice)[2]
	#matplot(((l$x/dim(zoom_fovea)[1])*scales[1])+x_add, ((l$y/dim(zoom_fovea)[2])*scales[2])+y_add, type='l', col="blue", add=T)
return(list("x" = ((l$x/dim(zoom_fovea)[1])*scales[1])+x_add, "y" = ((l$y/dim(zoom_fovea)[2])*scales[2])+y_add))
}

# function thinking about dip volume - looks strange
# @ fovea modelling
calc_dip_volume <- function(dip) {
	positions = check_up_and_down(dip)
	pin_pos = which(dip$x%in%positions)
	indices = which.min(dip$y[dip$x%in% positions])
	indices = c(pin_pos[indices-1], pin_pos[indices], pin_pos[indices+1])
	x_pos = dip$x[indices]
	y_pos = dip$y[indices]
return(list("xs"=x_pos, "ys"=y_pos))
}

# function thinking about extracting dip params 
# @ parameter extraction	
get_max_dip_params <- function(image_data) {
	slices = get_fovea_slices(image_data$img, get_fovea_position(image_data$img))
	x_fovea_params = fit_x_fovea(t(image_data$img[,,slices$midpoint]))
	return(	
		calc_dip_volume(
			fit_sholders(t(image_data$img[,,slices$midpoint])
				, x_fovea_params
				, fit_y_fovea(t(image_data$img[,,slices$midpoint]), x_fovea_params)
			)
		)
	)
}

# example of uses some of these functions
# @ examples
shoulder_example <- function() {
	root_dir = "/Users/tomas/projects/2016/eye_imaging/30_replicates/"
	image_dirs = dir(root_dir)
	x = 1

	image_d = paste(root_dir, "/", image_dirs[x], "/OCT/", sep="")
	rep_images = dir(image_d)

	rep_params = NULL
	for(y in 1:length(rep_images)) {
		rep_file = paste(image_d, "/", rep_images[y], sep="")
		image_data = dicomInfo(rep_file)
		

		# find fovea position
		fpos = get_fovea_position(image_data$img)
		slices = get_fovea_slices(image_data$img, fpos[1])
		slice = slices$midpoint

		fovea_slice = t(image_data$img[,,slice])
		x_fovea_params = fit_x_fovea(fovea_slice)
		y_fovea_params = fit_y_fovea(fovea_slice, x_fovea_params)

		dip = fit_sholders(fovea_slice, x_fovea_params, y_fovea_params)

		plot_fitted_slice(fovea_slice, x_fovea_params, y_fovea_params, dip)
		pars = get_max_dip_params(image_data)
		rep_params = rbind(rep_params, unlist(pars))
	}
}