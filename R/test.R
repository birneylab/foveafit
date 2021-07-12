create_tessalation <- function(m, fac=10) {
	inds = NULL
	s1 = seq(1, nrow(m), by=10)
	for(x in 1:length(s1)) {
		inds = rbind(inds, (s1))
	}
	inds = cbind(as.vector(inds), as.vector(t(inds)))
	dxy1 <- deldir(inds[,1], inds[,2])
	plot(inds/nrow(m))
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
        
        segments(u1, v1, u2, v2, col = "blue", lty = 1)
	 segments(x1, y1, x2, y2, col = "red", lty = 1)
}

#ind = 1
#summary(es[[ind]]$rep1$smodel)
#dim(es[[ind]]$rep1$smodel)
#dim(es[[ind]]$rep1$model)


#par(mfrow=c(1,2))
#image(t(es[[ind]]$rep1$model))
#image(t(es[[ind]]$rep1$smodel))


eye <- function (image_data, zpos, bands, name) 
{
    img = image_data$img
    image(t(img[, , zpos]), breaks = c(0, 50, 200), col = c("black", 
        "green"), main = name)
    for(x in 1:ncol(bands)) {
    matplot(1:dim(img)[2]/dim(img)[2], bands[, x]/dim(img)[1], 
        col = "red", type = "l", add = T)
    }
}


#m = mods$smodel




#r = apply(img[,,zpos], 1, mean)
#r = r-median(r)
#ps = which(abs(r) > median(r)+(sd(r)*f))
#m_bands = r[ps]
#ps = vector()
#for(x in 1:4) {
#	ps[x] = which.max(step4(m_bands))
#		plot(m_bands)
#		abline(v=ps[x])
#	m_bands = m_bands[ps[x]:length(m_bands)]
#}


#				img = image_data$img
#		outer_band_positions = NULL
#		for(y in 1:dim(img)[2]) {
#			r = runmed(img[,y,zpos], s)
#			r = r-median(r)
#			ps = which(abs(r) > median(r)+(sd(r)*f))
#			outer_band_positions = rbind(outer_band_positions, c(min(ps), max(ps)))
#		}
