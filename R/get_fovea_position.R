# calculate the x or y means across the z plane
# @ auxcillary 
get_z_mean <- function(img, rcol=2) {
	r = NULL
	for(x in 1:dim(img)[3]) {
		a = apply(img[,,x], rcol, mean)
		r = rbind(r, a)
	}
return(r)
}

# calculate the x or y means across the all plane
# @ auxcillary 
get_xyz_mean <- function(img, rcol=2) {
	x1<- y1 <- z1 <- NULL
	for(x in 1:dim(img)[1]) {
		a = apply(img[x,,], rcol, mean)
		x1 = rbind(x1, a)
	}
	for(x in 1:dim(img)[2]) {
		a = apply(img[,x,], rcol, mean)
		y1 = rbind(y1, a)
	}
	for(x in 1:dim(img)[3]) {
		a = apply(img[,,x], rcol, mean)
		z1 = rbind(z1, a)
	}
return(list("x"=x1, "y"=y1, "z"=z1))
}

# function attempting to get the slices that conatin the fovea - old
# @ fovea modelling 
get_fovea_slices <- function(image_data, fpos, fac=0.2) {
	img = image_data$img
	r = NULL
	for(x in 1:dim(img)[1]) {
		a = apply(img[x,(fpos$midpoint-5):(fpos$midpoint+5),], 2, mean)
		r = rbind(r, a)
	}
	l = lowess(apply(r, 2, mean), f=fac)
	up_down = check_up_and_down(l)
	while(length(up_down)<3) {
		fac = fac-0.01
		l = lowess(apply(r, 2, mean), f=fac)
		up_down = check_up_and_down(l)
	}
	spos = get_max_peak_trove(up_down, l)
	 #image(r)
	 #abline(h=unlist(spos)/128)
return(spos)
}

# function attempting to get the position of maximum dip - old but good
# @ fovea modelling 
get_fovea_position <- function(image_data, fac=0.2) {
	img = image_data$img
	r = NULL
	for(x in 1:dim(img)[3]) {
		a = apply(img[,,x], 2, mean)
		r = rbind(r, a)
	}
	l = lowess(apply(r, 2, mean), f=fac)
	up_down = check_up_and_down(l)
	while(length(up_down)<3) {
		fac = fac-0.01
		l = lowess(apply(r, 2, mean), f=fac)
		up_down = check_up_and_down(l)
	}
	fpos = get_max_peak_trove(up_down, l)
	#image(r)
	#abline(v=unlist(fpos)/dim(r)[2])
return(fpos)
}

# same of "check_up_and_down" but better name - TODO: replace with this better name
# @ fovea modelling 
switch_points <- function(l) {
	check = 1
	up_down = NULL
	for(x in 2:length(l$x)) {
		if(check==1 & l$y[x]<l$y[x-1]) {
			check = 2
			up_down = c(up_down, l$x[x])
		}
		if(check==2 & l$y[x]>l$y[x-1]) {
			check = 1
			up_down = c(up_down, l$x[x])
		}
	}
return(up_down)
}

# exampe of using some of these functions to extract slices and positions
# quite good approach but maybe not robust enough - needs further development
# @ fovea modelling 
get_fovea_slice_and_positions <- function(image_data, fac=0.2, band_f=10) {
	img = image_data$img
	r = sapply(1:dim(img)[3], function(x) rbind(apply(img[,,x], 2, mean)))
	l = lowess(apply(r, 2, mean), f=fac)
	up_down = up_down_reduce(l, r, fac)
	zpos = get_max_peak_trove(up_down, l)$midpoint		
	r = apply(img[,,zpos], 1, mean)
	sp_ep = get_peak_positions(r, fac)
	sp = min(sp_ep$sp,sp_ep$ep)-100; ep = max(sp_ep$sp,sp_ep$ep)+100;

	r = apply(img[sp:ep,,zpos], 2, mean)
	ys = lowess(r, f=fac)
	yy = up_down_reduce(ys, r, fac)
	ypos = which(ys$y==min(ys$y[yy]))

	yst = yy[which.min(abs(yy-ypos))-1]
	yso = yy[which.min(abs(yy-ypos))+1]
	clipped_img = img[sp:ep,yst:yso,zpos]
	
	r = apply(clipped_img, 1, mean)
	l = lowess(diff(r), f=0.1)
	peaks = which(diff(sign(diff(l$y)))==-2)+1
	troves = which(diff(sign(diff(-l$y)))==-2)+1
	swp = c(peaks, troves) 
	swp = swp[order(abs(l$y[swp]), decreasing=T)][1:5]

	image(clipped_img)
	abline(h=ypos/dim(img)[2])
	abline(v=swp/dim(clipped_img)[1])
	plot(l); abline(v=swp);

return(list("zslice"  = zpos, "ypos"=ypos, "xpoints"=swp))
}


