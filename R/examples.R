fovea_extract_example <- function() {
	
  load(system.file("data/example1", package="foveafit"))

	# find fovea position
	fpos = get_fovea_position(image_data)
	slices = get_fovea_slices(image_data, fpos$midpoint)
	slice = slices$midpoint

	fovea_slice = t(image_data$img[,,slice])
	x_fovea_params = fit_x_fovea(fovea_slice)
	y_fovea_params = fit_y_fovea(fovea_slice, x_fovea_params)

	dip = fit_sholders(fovea_slice, x_fovea_params, y_fovea_params)
	plot_fitted_slice(fovea_slice, x_fovea_params, y_fovea_params, dip)
}


model_dip <- function() {
	load(system.file("data/example1", package="foveafit"))
	fpos = get_fovea_position(image_data)
	slices = get_fovea_slices(image_data, fpos)
	x = 1
	z_dips = list()
	p_dips = list()
	for(slice in slices$l_shoulder:slices$r_shoulder) {
		fovea_slice = t(image_data$img[,,slice])
		x_fovea_params = fit_x_fovea(fovea_slice)
		if(length(x_fovea_params)>2) {
			y_fovea_params = fit_y_fovea(fovea_slice, x_fovea_params)
			if(length(y_fovea_params)>2 ) {
				z_dips[[x]] = fit_sholders(fovea_slice, x_fovea_params, y_fovea_params)
				ps = calc_dip_volume(z_dips[[x]])
				p_dips[[x]] = ps
				plot(z_dips[[x]])
				abline(v=ps$xs)
			}
		}
		x = x+1
	}
invisible(list("zs"=z_dips, "ps"=p_dips))
}

extract_peak_troves <- function(data, direction=1) {
	len = nrow(data)
	if(direction==2) {
		len = ncol(data)
	}
	pt_data = NULL
	for(x in 1:len) {
		if(direction==2) {
  		p = peak_trove(data[,x])
  	} else {
  		p = peak_trove(data[x,])
  	}
  	pt_data = rbind(pt_data, unlist(p))
	}
return(pt_data)
}


make_square <- function(data) {
	ndata = NULL
	ds = dim(data)
	if(ds[1]>ds[2]) {
		for(x in 1:ds[1]) {
			d = approx(data[x,], n=ds[1])
			ndata = rbind(ndata, d$y)
		}
	} else {
		for(x in 1:ds[2]) {
			d = approx(data[,x], n=ds[2])
			ndata = cbind(ndata, d$y)
		}
	}
return(ndata)
} 

