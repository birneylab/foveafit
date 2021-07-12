# @@ functions for doing things with peaks and troves - need documentating at somepoint

peak_trove <- function(dy, num=3, sdf = 6.4, lincr= 1000) {
	dx = 1:length(dy)
	spl <- smooth.spline(dx, dy, df = sdf)
	xx <- unique(sort(c(seq(0, 30, length = lincr), kn <- unique(dx))))
	pp = predict(spl, xx, deriv=1)
	ps = data.frame(pp$x, abs(pp$y))
	pt = ps[order(ps[,2]),1][1:num]
	pt = pt[order(pt)]
	depths = dy[pt[order(pt)]]-min( dy[pt[order(pt)]])
	depths = depths[depths!=0]
	a1 = (depths[1]*(pt[2]-pt[1]))/2
	a2 = (depths[2]*(pt[3]-pt[2]))/2 
	fovealoc = dy[pt[order(pt)[2]]]
	return( list("position" = pt, "diameter"= pt[length(pt)]-pt[1], "depths"= depths, "areas"= c(a1, a2), "radius1" = pt[2]-pt[1], "radius2" = pt[3]-pt[2], "fovealoc"=fovealoc )) 
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

get_peak_positions <- function(r, fac=0.2) {
	l = lowess(r, f=fac-0.1)
	t1 = sapply(2:length(l$y), function(x) l$y[x-1] > l$y[x])
	peakp =which(sapply(2:length(t1), function(x) t1[x+1]==TRUE & t1[x]==FALSE))
	trovep =which(sapply(2:length(t1), function(x) t1[x+1]==FALSE & t1[x]==TRUE))
	sp = which(l$y==max(l$y[peakp]))
	peakp = peakp[-(which.max(l$y[peakp]))]
	ep = which(l$y==max(l$y[peakp]))
	tp = min(trovep[trovep>max(sp, sp)])
return(list("sp"=min(sp, ep), "ep"=tp))
}

get_peak_troves <- function(l, fac=0.2) {
	t1 = sapply(2:length(l), function(x) l[x-1] > l[x])
	peakp =which(sapply(2:length(t1), function(x) t1[x+1]==TRUE & t1[x]==FALSE))
	trovep =which(sapply(2:length(t1), function(x) t1[x+1]==FALSE & t1[x]==TRUE))
return(list("peaks"=peakp, "troves"=trovep))
}

make_square <- function(data) {
	ndata = NULL
	ds = dim(data)
	if(ds[1]%%2!=0) {
		ds[1]=ds[1]-1
	}
	if(ds[2]%%2!=0) {
		ds[2]=ds[2]-1
	}
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

get_max_peak_trove <- function(up_down, l) {
	ds = vector()
	for(x in 2:(length(up_down)-1)) {
		if(l$y[up_down[x-1]]> l$y[up_down[x]] & l$y[up_down[x+1]]> l$y[up_down[x]])
			ds[x] = abs(l$y[up_down[x]] - l$y[up_down[x-1]]) + abs(l$y[up_down[x]] - l$y[up_down[x+1]])
	}
	ds[1] = 0; ds = c(ds, 0)
	#pos = which.min(abs(diff(l$y[up_down])))+1
return(list("midpoint"=up_down[which.max(ds)]
						, "l_shoulder"=up_down[which.max(ds)-1]
						, "r_shoulder"=up_down[which.max(ds)+1]))
}
