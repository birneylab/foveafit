
fit_deriv <- function(data, debug=F, f=3) {
	time = 1:length(data)
	div1 = get_derivative(data)
	p.t = c(time[which(diff(sign(diff(div1)))==+2)+1]
					, time[which(diff(sign(diff(div1)))==-2)+1])
	p.t = p.t[order(p.t)]
	div2 = get_derivative(div1)
	p.t2 = c(time[which(diff(sign(diff(div2)))==+2)+1]
					, time[which(diff(sign(diff(div2)))==-2)+1])
	p.t2 = p.t2[order(p.t2)]

	pt = p.t2[p.t2>length(data)/f &  p.t2<length(data)-(length(data)/f)]
	center = pt[which.max(div2[pt])]
	pt = p.t[p.t<center]; max.s1 = pt[which.min(abs(center-pt))]
	pt = p.t[p.t>center]; max.s2 = pt[which.min(abs(center-pt))]

	pt = p.t2[p.t2<max.s1]; max.l1 = pt[which.min(abs(center-pt[pt<max.s1]))]
	pt = p.t2[p.t2>max.s2]; max.l2 = pt[which.min(abs(center-pt[pt>max.s2]))]

	if(debug) {	
		plot(time, data)
		matplot(time, rescale(div1, range(data)), add=T, type="l", col="blue", lty="dashed", lwd=0.5)
		matplot(time, rescale(div2, range(data)), add=T, type="l", col="green", lty="dashed", lwd=0.5)
		abline(v=max.s1, col="blue")
		abline(v=max.s2, col="blue")
		abline(v=center, col="red")
		abline(v=max.l1, col="green")
		abline(v=max.l2, col="green")
	}

return(list("midpoint"=center
						, "l_div1"=max.s1
						, "r_div1"=max.s2
						, "l_div2"=max.l1
						, "r_div2"=max.l2))
}

extract_circle_params <- function(m) {
	center = floor(dim(m)/2)
	radius = dim(m)/2
	coors = data.frame(center + radius * sin(seq(0, 2*pi, length = 360))
					, center + radius * cos(seq(0, 2*pi, length = 360)))
	colnames(coors) = c("x", "y") 
	s.ss = list()
	i.ss = list()
	for( ang in 1:360) {
		ang2 = abs(ang+179)
		if(ang2>360) {
			ang2 = ang2-360
		}
		inds1 = unlist(abs(floor(coors[ang,])))
		inds2 = unlist(abs(floor(coors[ang2,])))
		inds1[inds1[1]==0][1] = 1; inds1[inds1[2]==0][2] = 1;
		inds2[inds2[1]==0][1] = 1; inds2[inds2[2]==0][2] = 1; 
		s1 = seq(inds1[1], inds2[1], length=dim(m)[1])
		s2 = seq(inds1[2], inds2[2], length=dim(m)[1])
		d = sapply(1:dim(m)[1], function(x) m[s1[x], s2[x]])
		s.ss[[ang]] = d
		i.ss[[ang]] = list()
		i.ss[[ang]]$i1 = inds1
		i.ss[[ang]]$i2 = inds2
	}
return(list("coors"=coors, "s.ss"=s.ss, "i.ss"=i.ss, "n"=dim(m)[1]))
}

convert_coords <- function(v360, divs) {
	f.xy1 = NULL
	f.xy2 = NULL
	for(ang  in 1:length(v360$i.ss)) {		
		ang2 = abs(ang+179)
		if(ang2>360) {
			ang2 = ang2-360
		}
		inds1 = unlist(abs(floor(v360$coors[ang,])))
		inds2 = unlist(abs(floor(v360$coors[ang2,])))
		inds1[inds1[1]==0][1] = 1; inds1[inds1[2]==0][2] = 1;
		inds2[inds2[1]==0][1] = 1; inds2[inds2[2]==0][2] = 1; 
		s1 = seq(inds1[1], inds2[1], length=v360$n)
		s2 = seq(inds1[2], inds2[2], length=v360$n)
		f.xy1 =rbind(f.xy1, s1[divs[ang,]])
		f.xy2 =rbind(f.xy2, s2[divs[ang,]])
	}
return(list("f.xy1"=f.xy1, "f.xy2"=f.xy2))
}

center_and_square <- function(model) {
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
	centered_band = model$model$top_band[rangs[[1]], rangs[[2]]]
	centered_band = make_square(centered_band)
return(centered_band)
} 

fit_centered_surface <- function(centered_band, debug=FALSE) {
	v360 = extract_circle_params(centered_band)
	divs = do.call(rbind, lapply(lapply(v360$s.ss, fit_deriv, debug), unlist))
	xy.s = convert_coords(v360, divs)
	max_slope = cbind(c(xy.s$f.xy1[,2], xy.s$f.xy1[,3]), c(xy.s$f.xy2[,2], xy.s$f.xy2[,3]))
	plateu = cbind(c(xy.s$f.xy1[,4], xy.s$f.xy1[,5]), c(xy.s$f.xy2[,4], xy.s$f.xy2[,5]))
return(list("circle_params"=v360, "max_slope"=max_slope, "plateu"=plateu, "divs"=divs))
}

registar <- function(image_file, debug=FALSE) {
	model = run_generate_model_and_find_fovea_dip_poisition(image_file)
	if(abs(0.5 - c(model$zcenter)/dim(model$model$top_band)[1])>0.15) {
		model$zcenter = dim(model$model$top_band)[1]/2
	}
	if(abs(0.5 - c(model$ycenter)/dim(model$model$top_band)[2])>0.15) {
		model$ycenter = dim(model$model$top_band)[2]/2
	}
	centered_band = center_and_square(model)
	fitted_center = fit_centered_surface(centered_band, debug)

	e1 = get.ellipse(fit.ellipse(fitted_center$max_slope[,1], fitted_center$max_slope[,2]))
	e2 = get.ellipse(fit.ellipse(fitted_center$plateu[,1], fitted_center$plateu[,2]))

	if(debug) {
		image(t(centered_band-median(centered_band)), col=topo.colors(10))
		lines(e1/dim(centered_band)[1], lwd=2, col="blue")
		lines(e2/dim(centered_band)[1], lwd=2, col="blue")
		matplot(fitted_center$max_slope[,1]/dim(centered_band), fitted_center$max_slope[,2]/dim(centered_band), col="red", add=T, pch="+", cex=0.5)
		matplot(fitted_center$plateu[,1]/dim(centered_band), fitted_center$plateu[,2]/dim(centered_band), col="red", add=T, pch="+", cex=0.5)
	}
return(list("e1"=e1, "e2"=e2, "model"=model$model, "ycenter"  = model$ycenter, "zcenter" = model$zcenter
			, "centered_band"=centered_band, "fitted_center"=fitted_center))
}

fovea_overlay <- function(data_folder = "/Users/tomas/projects/2016/eye_imaging/30_replicates") {

	file_dirs = paste(data_folder, dir(data_folder), "OCT", sep="/")
	es = list()
	average_model = NULL
	for(x in 1:length(file_dirs)) {
	#for(x in 1:10) {	
		files = paste(file_dirs[x], dir(file_dirs[x]), sep="/")
		e1 = registar(files[1], TRUE)
		e2 = registar(files[2], TRUE)
		es[[x]] = list()
		es[[x]]$rep1 = e1
		es[[x]]$rep2 = e2
		if(is.null(average_model)) {
			average_model = e1$model-median(e1$model)
			average_model = average_model + e2$model-median(e2$model)
		} else {
			average_model = average_model + e1$model-median(e1$model)
			average_model = average_model + e2$model-median(e2$model)
		}	
	}
	
	image(t(average_model%/%10), col=topo.colors(50))
	par(mfrow=c(1,2))
	for(z in 1:30) {
		centers =(dim(es[[z]]$rep1$model)/2)-c(es[[z]]$rep1$zcenter, es[[z]]$rep1$ycenter)+( max(dim(es[[z]]$rep1$model)-es[[z]]$rep1$ds))
		#image(t(es[[z]]$rep1$model), col=topo.colors(10))
		#lines(((es[[z]]$rep1$e[,1])/512)+((es[[z]]$rep1$zcenter/dim(es[[z]]$rep1$model)[1])/2), ((es[[z]]$rep1$e[,2])/512)+((es[[z]]$rep1$ycenter/dim(es[[z]]$rep1$model)[2])/2), col="red")

		lines((es[[z]]$rep1$e[,1]+centers[1])/512+((es[[z]]$rep1$ycenter/dim(es[[z]]$rep1$model)[2])-0.5), (es[[z]]$rep1$e[,2]+centers[2])/512+((es[[z]]$rep1$zcenter/dim(es[[z]]$rep1$model)[1])-0.5), col="red")

		centers =(dim(es[[z]]$rep2$model)/2)-c(es[[z]]$rep2$zcenter, es[[z]]$rep2$ycenter)+( max(dim(es[[z]]$rep2$model)-es[[z]]$rep2$ds))
		#image(t(es[[z]]$rep2$model), col=topo.colors(10))
		#lines((es[[z]]$rep2$e[,1]+centers[1])/512, (es[[z]]$rep2$e[,2]+centers[2])/512, col="red")

		lines((es[[z]]$rep2$e[,1]+centers[1])/512+((es[[z]]$rep2$ycenter/dim(es[[z]]$rep2$model)[2])-0.5), (es[[z]]$rep2$e[,2]+centers[2])/512+((es[[z]]$rep2$zcenter/dim(es[[z]]$rep2$model)[1])-0.5), col="red")

	}
}



#### currently unused:



get_circle_vectors <- function(centered_band) {
 	diam <- unique(dim(centered_band))
 	center <- radius <- diam/2
  coors = data.frame(center + radius * cos(seq(0, 2*pi, length = 360))
   										,  center + radius * sin(seq(0, 2*pi, length = 360))) 
  grand = diam/180
  x.xs = NULL
  for (x in 1:length(coors[,1])) {
   	x0 = floor(grand*x)
   	x1 = coors[x,1]
   	y0 = dim(centered_band)[1]-coors[x,2]
   	y1 = coors[x,2]
   	x.xs = rbind(x.xs, cbind(x0, x1, y0, y1))
   }
   x.xs = x.xs[x.xs[,1]<=diam,]
   s.ss = list()
   i.ss = list()
   for(x in 1:length(x.xs[,1])) {
   	se1 = round(seq(x.xs[x,1], x.xs[x,2], length=diam))
   	se2 = round(seq(x.xs[x,3], x.xs[x,4], length=diam))
   	se1[se1<1] = 1; se1[se1>diam] = diam;
		se2[se2<1] = 1; se2[se2>diam] = diam;
		ss = vector()
		 for(y in 1:length(se1)) {
				ss[y] = centered_band[se2[y],se1[y]]
			}
		i.ss[[x]] = list()
		i.ss[[x]]$se1 = se1
		i.ss[[x]]$se2 = se2 	
		s.ss[[x]] = ss	
	}
return(list("x.xs"=x.xs, "s.ss"=s.ss, "i.ss"=i.ss))
}

convert_x_y <- function(p.ps, n) {
	xy.s = NULL
	xy.ss = NULL
	for(a in 1:dim(p.ps)[1]) {
		cx <- cy <- (n/2)
		r = mean(p.ps[,1]) - p.ps[a,2]
		x = cx + r * cos(a)
		y = cy + r * sin(a)
		xy.s = rbind(xy.s, cbind(x, y))
		a = (180-a)+1
		r = p.ps[a,3] - mean(p.ps[,1])
		cx <- cy <- (n/2)
		x = cx + r * cos(a)
		y = cy + r * sin(a)
		xy.s = rbind(xy.s, cbind(x, y))
		r = abs (p.ps[a,1]-(diam/2))
		x = cx + r * cos(a)
		y = cy + r * sin(a)
		xy.ss = rbind(xy.ss, cbind(x, y))
	}
return(xy.s)
}

convert_x_y2 <- function(p.ps, n) {
	xy.s = NULL
	for(a in 1:dim(p.ps)[1]) {
		cx <- cy <- (n/2)
		r = mean(p.ps[,1]) - p.ps[a,2]
		x = cx + r * cos(a)
		y = cy + r * sin(a)
		xy.s = rbind(xy.s, cbind(x, y))
		#a = (180-a)+1
		r = p.ps[a,3] - mean(p.ps[,1])
		cx <- cy <- (n/2)
		x = cx + r * cos(a)
		y = cy + r * sin(a)
		xy.s = rbind(xy.s, cbind(x, y))
	}
return(xy.s)
}



