create_tessalation <- function(m, fac=10, main="", brks = 10) {
	
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
   image((m), col=topo.colors(brks), main=main)     
  # segments(u1, v1, u2, v2, col = "red", lty = 1, lwd=0.2)
	# segments(x1, y1, x2, y2, col = "red", lty = 1, lwd=0.2)
}

square_tes_extract <- function(m, fac=0.9) {
	values = NULL
	m = m-median(m)
	s1 = seq(1, dim(m)[1], length=dim(m)[1]*fac)
	s2 = seq(1, dim(m)[2], length=dim(m)[2]*fac)
	for(x in 2:length(s1)) {
		m1 = m[s1[x-1]:s1[x],]
		s = sapply(2:length(s2), function(i) mean(m1[,s2[i-1]:s2[i]]))
		values = rbind(values, s)
	}
	mm = m
	for(x in 2:length(s1)) {
		for(i in 2:length(s2)) { 
			mm[s1[x-1]:s1[x],s2[i-1]:s2[i]]= values[x-1,i-1]
		}
	}
	#image(values)
return(mm)
}

tri_tes_extract <- function(m, fac=10, debug=FALSE) {
	triangles_from_square <- function(m1) {
			tri1 = NULL; tri2 = NULL
			for(y in 1:nrow(m1)) {
				tri1 = c(tri1, m1[y, 1:y])
				tri2 = c(tri2, m1[y,1:(nrow(m1)-y)])
			}
		return(list("t1" =tri1, "t2"=tri1))
	}
	values = NULL
	s1 = seq(1, dim(m)[1], by=fac)
	s2 = seq(1, dim(m)[2], by=fac)
	for(x in 2:length(s1)) {
		m1 = m[s1[x-1]:s1[x],]
		s = sapply(2:length(s2), function(i) unlist(lapply(triangles_from_square(m1[,s2[i-1]:s2[i]]), mean)))
		values = rbind(values, s)
	}
	if(debug) {
		image((values), col=topo.colors(1000))
		abline(v=seq(0,1, length=length(s1)))
		abline(h=seq(0,1, length=length(s2)))
		#abline(v=s1/dim(m)[1], col="lightgrey")
		#abline(h=s2/dim(m)[2], col="lightgrey")
		for(x in 2:length(s1)) {
			for(y in 2:length(s2)) {
				#segments(s1[x-1]/dim(m)[1], s2[y-1]/dim(m)[2], s1[x]/dim(m)[1], s2[y]/dim(m)[2])
			}
		} 
	}
return(values)
}

align_model <- function(top_band_surface, ns=c(75, 306)) {
	ndata = NULL
	ds = dim(top_band_surface)
	st = ns/2
	ndata = top_band_surface[((ds[1]/2)-st[1]): ((ds[1]/2)+st[1]), ((ds[2]/2)-st[2]): ((ds[2]/2)+st[2])]
return(ndata)
} 

center_align<- function(model, f=0.6) {
	model$ycenter = fit_deriv(model$model$top_band[model$zcenter,])$midpoint
	model$zcenter = fit_deriv(model$model$top_band[,model$ycenter])$midpoint
	centers =(dim(model$model$top_band)/2)-c(model$zcenter, model$ycenter)
	rangs = list()
	for(x in 1:length(centers)) {
	 	if(centers[x]<0) {
	 		r = abs(centers[x]*2):dim(model$model$top_band)[x]
	 	} else {
	 		r = 1:(dim(model$model$top_band)[x]-centers[x]*2)
	 }
	 rangs[[x]] = r
	}
	ns=dim(model$model$top_band) * f
	centered_band = model$model$top_band[rangs[[1]], rangs[[2]]]
	centered_band = align_model(centered_band, ns)
return(centered_band)
} 

example <- function() {
	image_file = "/Users/tomas/projects/2016/eye_imaging/30_replicates/1851037/OCT/405448.dcm"
	model = run_generate_model_and_find_fovea_dip_poisition(image_file)


	par(mfrow=c(2,4))
	s = seq(0.9, 0.3, length=7)
	image(model$model$top_band, main="100%")
	for(x in 1:length(s)) {
		sc = s[x]*100
		amodel = center_align(model, s[x])
		image(amodel, main=paste(sc, "%", sep=""))
	}

	tes_f = 5
	v1 = square_tes_extract(amodel, tes_f)





	data_folder = "/Users/tomas/projects/2016/eye_imaging/30_replicates"
	file_dirs = paste(data_folder, dir(data_folder), "OCT", sep="/")
	image(t(m))
	cens = NULL
	for(x in 1:length(file_dirs)) {
		files = paste(file_dirs[x], dir(file_dirs[x]), sep="/")
		image_file = files[1]
		model = run_generate_model_and_find_fovea_dip_poisition(image_file)
		amodel = center_align(model, 0.7)
		image(amodel, main=paste(sc, "%", sep=""))

		#model$ycenter = fit_deriv(model$model$top_band[model$zcenter,])$midpoint
		#model$zcenter = fit_deriv(model$model$top_band[,model$ycenter])$midpoint
		#centers =(dim(model$model$top_band)/2)-c(model$zcenter, model$ycenter)
		#mod = model$model$top_band
		#m = matrix(nrow=dim(mod)[1], ncol=dim(mod)[2], 1)
		#abline(v=model$ycenter/dim(mod)[2], lty="dashed")
		#abline(h=model$zcenter/dim(mod)[1], lty="dashed")
		#matplot(model$zcenter/dim(mod)[1], model$ycenter/dim(mod)[2], pch="+", add=T)
		#cens = rbind(cens, c(model$zcenter, model$ycenter))
		print(x)
	}

}




