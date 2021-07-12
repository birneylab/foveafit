# function to find the maximum foveal dip slice - not robust enough
# @ fovea slice 
get_slice <- function(image_data) {
		img = image_data$img
		r = sapply(1:dim(img)[3], function(x) rbind(apply(img[,,x], 2, mean)))
		l = lowess(apply(r, 2, mean), f=fac)
		up_down = up_down_reduce(l, r, fac)
		zpos = get_max_peak_trove(up_down, l)$midpoint
return(zpos)
}

# robust function to find the maximum foveal dip slice 
# @ fovea slice 
get_max_fovea_slice_robust <- function(image_data, fac=0.15) {
	bulls_eye_x = sapply(1:dim(image_data$img)[3], function(i) apply(image_data$img[,,i],1, mean))
	aband = apply(bulls_eye_x, 2, mean)
	l1 = lowess(aband, f=fac)
	center = get_fovea_slice_robust(l1)
return( list("center" = center, "fit"=l1, "aband"=aband) )
}

# helper function to solve formulas 
# @ MATLAB backslash function 
mldivide <- function (A, B, pinv = TRUE) {
    stopifnot(is.numeric(A) || is.complex(A), is.numeric(B) ||  is.complex(B))
    if (is.vector(A)) 
        A <- as.matrix(A)
    if (is.vector(B)) 
        B <- as.matrix(B)
    if (nrow(A) != nrow(B)) 
        stop("Matrices 'A' and 'B' must have the same number of rows.")
    if (pinv) {
        pinv(t(A) %*% A) %*% t(A) %*% B
    }
    else {
        qr.solve(A, B)
    }
}

# function to calculate the curvature of based on two vectors
# @ parameter extraction 
calculate_curvature <- function(x, y) {
	mx = mean(x); my = mean(y)
 	X = x - mx; Y = y - my; # Get differences from means
 	dx2 = mean(X^2); dy2 = mean(Y^2); # Get variances
 	t = mldivide(cbind(X,Y), (X^2-dx2+Y^2-dy2)/2); # Solve least mean squares problem
 	a0 = t[1]; b0 = t[2]; # t is the 2 x 1 solution array [a0;b0]
 	r = sqrt(dx2+dy2+a0^2+b0^2); # Calculate the radius
 	a = a0 + mx; b = b0 + my; # Locate the circles center
	t <- seq(0,2*pi,length=100) 
	coords <- t(rbind( a+sin(t)*r, b+cos(t)*r)) 
	plot(y, ylim=c(min(coords[,2]), max(coords[,2])))
	matplot(coords[,1], coords[,2], col="red", pch=20, add=T)
 	curv = 1/r;  #Get the curvature
}

# calculate curvature of each side of the oct slice - assume a dip is present
# @ parameter extraction
find_dip_curves <- function(iband, fac=0.1) {
	mid = get_fovea_slice_robust(lowess(iband, f=fac))
	c1 = calculate_curvature(1:mid, iband[1:mid])
	c2 = calculate_curvature(mid:length(iband), iband[mid:length(iband)])
return( list("c1"=c1, "c2"=c2) )
}

# calculate curvature of fovea bands
# @ parameter extraction
extract_fovea_band_parameters <- function(iband, fac=0.1) {
	# fovea mid point
	mid = get_fovea_slice_robust(lowess(iband, f=fac))
	# curvature parameters
	c = calculate_curvature(1:length(iband), iband)
	c1 = calculate_curvature(1:mid, iband[1:mid])
	c2 = calculate_curvature(mid:length(iband), iband[mid:length(iband)])
	# slope parameters

	ps = 1:length(iband)
	s = coef(lm(iband~ps))[2]
	s1 = 0#lm(iband[1:mid]~1:mid)$slope
	s2 = 0#lm(iband[mid:length(iband)]~mid:length(iband))$slope
	eturn( list("s"=s, "s1"=s1, "s2"=s2, "c"=c, "c1"=c1, "c2"=c2) )
}

get.ellipse <- function( fit, n=360 )  {
  # Calculate points on an ellipse described by 
  # the fit argument as returned by fit.ellipse 
  # 
  # n is the number of points to render 
  
  tt <- seq(0, 2*pi, length=n) 
  sa <- sin(fit$angle) 
  ca <- cos(fit$angle) 
  ct <- cos(tt) 
  st <- sin(tt) 
  
  x <- fit$center[1] + fit$maj * ct * ca - fit$min * st * sa 
  y <- fit$center[2] + fit$maj * ct * sa + fit$min * st * ca 
  
  cbind(x=x, y=y) 
}


# noise measure robust to ouliers - IQR of differences between consecutive points
# @ auxcillary
dLRs <- function (x) {
	return(IQR(diff(na.omit(x)))/(4 * qnorm((1 + 0.5)/2)/sqrt(2)))
}

# fits the top and bottom major bands of a oct slice
# @ surface modelling
fit_major_bands <-function(image_data, zpos, s=5, f=2) {
		img = image_data$img
		outer_band_positions = NULL
		for(y in 1:dim(img)[2]) {
			r = runmed(img[,y,zpos], s)
			r = r-median(r)
			ps = which(abs(r) > median(r)+(sd(r)*f))
			outer_band_positions = rbind(outer_band_positions, c(min(ps), max(ps)))
		}
return(outer_band_positions)
}

# fits the top and bottom major bands of a oct slice
# @ surface modelling
fit_major_bands_bayes <-function(image_data, zpos, s=5, f=2) {
		img = image_data$img
		outer_band_positions = NULL
		for(y in 1:dim(img)[2]) {
			r = runmed(img[,y,zpos], s)
			r = r-median(r)
			ps = which(abs(r) > median(r)+(sd(r)*f))
			outer_band_positions = rbind(outer_band_positions, c(min(ps), max(ps)))
		}
return(outer_band_positions)
}


# interpolate outlier values uses smoothed values
# @ surface modelling
interpolate_outliers <- function(data, bf=3, sf=51) {
	ofac = dLRs(data)*bf
	sdata = runmed(data, sf)
	w = which(abs(sdata-ofac) > data | abs(sdata+ofac) < data)
	data[w] = sdata[w]
return(data)
}

# models the surface of the top and bottom major bands acorss oct scan slices
# @ surface modelling
model_fovea_surface <- function(image_data) {
	model_top_surface = NULL
	model_bottom_surface = NULL
	pb <- txtProgressBar(min = 0, max = dim(image_data$img)[3], style = 3);
	for(x in 1:dim(image_data$img)[3]) {
		bands = fit_major_bands(image_data, x)
		ibands = interpolate_outliers(bands)
		model_bottom_surface = rbind(model_bottom_surface, ibands[,1])
		model_top_surface = rbind(model_top_surface, ibands[,2])
		setTxtProgressBar(pb, x)
	} 
	close(pb);
return( list("top_band"= model_top_surface, "bottom_band"=model_bottom_surface) )
}

# plotting function for a fitted slice 
# @ visulisation
eye <- function(image_data, zpos, bands, name) {
	img = image_data$img
	image(t(img[,,zpos]), breaks=c(0,50,200), col=c("black", "green"), main=name)
	matplot(1:dim(img)[2]/dim(img)[2], bands[,1]/dim(img)[1], col="red", type="l",add=T)
	matplot(1:dim(img)[2]/dim(img)[2], bands[,2]/dim(img)[1], col="red", type="l",add=T)
}

# plotting function for surface model
# visulisation - uses "plot.ly" web service
surface_plot <- function(model_surface, file_name="temp") {
	require(plotly)
	Sys.setenv("plotly_username"="tomas81")
	Sys.setenv("plotly_api_key"="qmqhcq94nv")
	p <- plot_ly(z = ~model_surface) %>% add_surface()
	chart_link = plotly_POST(p, filename=file_name)
	chart_link
}

# example of running some modelling and parameter extraction and plotting
# @ example
run_example <- function(image_file, fac=0.15) {
	image_data = dicomInfo(image_file)
	zpos = get_max_fovea_slice_robust(image_data, fac)
	bands = fit_major_bands(image_data, zpos$center)
	ibands = interpolate_outliers(bands)
	cs = find_dip_curves(ibands[,2])
	eye(image_data, zpos$center, ibands, paste(unlist(strsplit(image_file, "/"))[c(8,10)], collapse="-"))
	plot(zpos$aband)
	matplot(zpos$fit$x, zpos$fit$y, pch=20, col="green", add=T)
	abline(v=zpos$center, col="green")
	legend("topleft", c(paste("zpos:", zpos$center)), col=c("black"))
return(cs)
}

# example of current best modelling and fovea position finding approach
# @ example
run_generate_model_and_find_fovea_dip_poisition <- function(image_file = '~/projects/2016/eye_imaging/30_replicates/1846403/OCT/300812.dcm') {
	image_data = dicomInfo(image_file)
	smodel= model_fovea_surface(image_data)
	center = get_max_fovea_slice_robust(image_data, 0.15)$center
	y_center = get_fovea_slice_robust(lowess(smodel$top_band[center,], f=0.15), lwb=1, upb=ncol(smodel$top_band))
	z_center = get_fovea_slice_robust(lowess(smodel$top_band[,y_center], f=0.15), lwb=1, upb=nrow(smodel$top_band))
return (list("model"=smodel, "zcenter"=z_center, "ycenter"=y_center))
}

 eye <- function (image_data, zpos, bands, name) {
    img = image_data$img
    image(t(img[, , zpos]), breaks = c(0, 50, 200), col = c("black", 
        "green"), main = name)
    for(x in 1:ncol(bands)) {
    	matplot(1:dim(img)[2]/dim(img)[2], bands[, x]/dim(img)[1], 
        col = "red", type = "l", add = T)
  	}
    #matplot(1:dim(img)[2]/dim(img)[2], bands[, 2]/dim(img)[1], 
    #    col = "red", type = "l", add = T)
}

fit_four_bands <- function(image_file, f = 2, s = 5, debug=FALSE) {
	model = run_generate_model_and_find_fovea_dip_poisition(image_file)
		image_data = dicomInfo(image_file)
		img = image_data$img
		zpos = model$zcenter
		outer_band_positions = NULL
		for(y in 1:dim(img)[2]) {
			d = img[,y,zpos]
			r = runmed(img[,y,zpos], s)
			r = r-median(r)
			ps = which(abs(r) > median(r)+(sd(r)*f))
			nps = which(r < -sd(abs(d))*f)
			i1 = min(nps)
			if(i1==Inf) {
				i1=1
			}
			d2 = d[i1:max(ps)]
			pe = step4(d2)
			index = which(diff(sign(diff(pe))) == -2)+1
			inds_14 = index[order(pe[index], decreasing=T)[c(1:2)]]
	
			pr = c(min(ps), max(ps), min(nps), min(nps)+inds_14)

			outer_band_positions = rbind(outer_band_positions, pr[order(pr)])
		}
		for(y in 1:ncol(outer_band_positions)) {
			outer_band_positions[is.na(outer_band_positions[,y]),y] = median(outer_band_positions[,y], na.rm=T)
			outer_band_positions[outer_band_positions[,y]==Inf,y]= median(outer_band_positions[,y], na.rm=T)
			
		}
	if(debug) {
		eye(image_data, zpos, interpolate_outliers(outer_band_positions), "")
	}
return(list("image_data"=image_data,"zpos"=zpos, "outer_band_positions"=outer_band_positions))
}

create_tessalation <- function(m, fac=10) {
	
	inds = NULL
	s1 = seq(1, nrow(m), by=fac)
	for(x in 1:length(s1)) {
		inds = rbind(inds, (s1))
	}
	inds = cbind(as.vector(inds), as.vector(t(inds)))
	dxy1 <- deldir(inds[,1], inds[,2])
	delsgs <- dxy1$delsgs/nrow(m)
	dirsgs = dxy1$dirsgs/nrow(m)
        x1 <- delsgs[, 1]
        y1 <- delsgs[, 2]
        x2 <- delsgs[, 3]
        y2 <- delsgs[, 4]
        u1 <- dirsgs[, 1]
        v1 <- dirsgs[, 2]
        u2 <- dirsgs[, 3]
        v2 <- dirsgs[, 4]
   image((m), col=topo.colors(10))     
   segments(u1, v1, u2, v2, col = "blue", lty = 1, lwd=0.2)
	 segments(x1, y1, x2, y2, col = "blue", lty = 1, lwd=0.2)
}

# random functions used during development
# @ development
random_functions <- function() {

	setwd("~/projects/2016/eye_imaging/30_replicates/")
	data_folder = "/Users/tomas/projects/2016/eye_imaging/30_replicates"
	file_dirs = paste(data_folder, dir(data_folder), "OCT", sep="/")

	#pdf("test.pdf", width=800, height=800)
	ccs = NULL
	par(mfrow=c(2,2))
	l1 = list()
	l2 = list()
	for(x in 1:length(file_dirs)) {
		files = paste(file_dirs[x], dir(file_dirs[x]), sep="/")
		t1 = test(files[1], 0.15)
		t2 = test(files[2], 0.15)
		ccs = rbind(ccs, cbind(t1$c, t1$s, t2$c, t2$s))
		print(x)
	}

	mat = matrix(nrow=30, ncol=30)
	for(x in 1:30) {
		x1 = l1[[x]][,2]
		for( y in 1:30) {
			mat[x,y] = cor(l2[[y]][,2], x1)
		}
	}

	par(mfrow=c(1,2))
	params = NULL
	for(x in 1:length(file_dirs)) {
		files = paste(file_dirs[x], dir(file_dirs[x]), sep="/")
		image_data = dicomInfo(files[1])
		run(image_data,  paste(unlist(strsplit(files[1], "/"))[c(8,10)], collapse="-"))
		image_data = dicomInfo(files[2])
		run(image_data,  paste(unlist(strsplit(files[2], "/"))[c(8,10)], collapse="-"))
	}

	# max - mean which is use - how to combine?
	run <- function(image_file, fac=0.2) {
		image_data = dicomInfo(image_file)
		#zpos = get_zpos(image_data, fac)
		zpos = get_max_fovea_slice_robust(image_data, fac)
		bands = fit_major_bands(image_data, zpos$center)
		ibands = interpolate_outliers(bands)
		eye(image_data, zpos$center, ibands, paste(unlist(strsplit(image_file, "/"))[c(8,10)], collapse="-"))

		bulls_eye_x = sapply(1:dim(image_data$img)[3], function(i) apply(image_data$img[,,i],1, mean))
		aband = apply(bulls_eye_x, 2, mean)
		aband_max = apply(bulls_eye_x, 2, max)
		aband_min = apply(bulls_eye_x, 2, min)
		l1 = lowess(aband, f=fac)
		center1 = get_fovea_slice_robust(l1)
		l2 = lowess(aband_min, f=fac)
		center2 = get_fovea_slice_robust(l2)
		l3 = lowess(aband_max, f=fac)
		center3 = get_fovea_slice_robust(l3)
		image(image_data$img[,,center1])
		plot(aband)
		matplot(l1$x, l1$y, pch=20, col="green", add=T)
		abline(v=center1, col="green")
		abline(v=center2, col="blue")
		abline(v=center3, col="red")
		print(center)
	}
}


