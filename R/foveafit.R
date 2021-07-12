#library(foveafit)

oct <- setRefClass("oct",
      fields = list( imagefile = "character"
      						, imagedata = "list"
      						, model = "list"
      						, fitted_model = "list")
)

oct$methods (
		read = function() {
			.self$imagedata = dicomInfo(.self$imagefile)
		},
		generate_model = function(f = 0.15) {
			smodel= model_fovea_surface(.self$imagedata)
			center = get_max_fovea_slice_robust(.self$imagedata, f)$center
			y_center = get_fovea_slice_robust(lowess(smodel$top_band[center,], f=0.15), lwb=1, upb=ncol(smodel$top_band))
			z_center = get_fovea_slice_robust(lowess(smodel$top_band[,y_center], f=0.15), lwb=1, upb=nrow(smodel$top_band))
			.self$model = list("model"=smodel, "zcenter"=z_center, "ycenter"=y_center)
		},
		register = function() {
			centered_band = center_and_square(.self$model)
			fitted_center = fit_centered_surface(centered_band, FALSE)
			e1 = get.ellipse(fit.ellipse(fitted_center$max_slope[,1], fitted_center$max_slope[,2]))
			e2 = get.ellipse(fit.ellipse(fitted_center$plateu[,1], fitted_center$plateu[,2]))
			.self$fitted_model = list("e1"=e1, "e2"=e2, "centered_band"=centered_band, "fitted_center"=fitted_center)
		},
		example = function() {
			o = oct(imagefile="/Users/tomas/projects/2016/eye_imaging/30_replicates/1331287/OCT/306418.dcm")
			o$read()
			o$generate_model()
			o$register()
		}

)


setClass("oct", slots=list( imagefile = "character"
      						, imagedata = "list"
      						, model = "list"
      						, fitted_model = "list")
) 


setGeneric("oct", function(x, ...) standardGeneric("oct"))

setMethod("oct", "character", function(x, ...) {
	imagedata = dicomInfo(x)
	smodel= model_fovea_surface(imagedata)
	center = get_max_fovea_slice_robust(imagedata, f=0.15)$center
	y_center = get_fovea_slice_robust(lowess(smodel$top_band[center,], f=0.15), lwb=1, upb=ncol(smodel$top_band))
	z_center = get_fovea_slice_robust(lowess(smodel$top_band[,y_center], f=0.15), lwb=1, upb=nrow(smodel$top_band))
  model = list("model"=smodel, "zcenter"=z_center, "ycenter"=y_center)
  centered_band = center_and_square(model)
	fitted_center = fit_centered_surface(centered_band, FALSE)
	e1 = get.ellipse(fit.ellipse(fitted_center$max_slope[,1], fitted_center$max_slope[,2]))
	e2 = get.ellipse(fit.ellipse(fitted_center$plateu[,1], fitted_center$plateu[,2]))
	fitted_model = list("e1"=e1, "e2"=e2, "centered_band"=centered_band, "fitted_center"=fitted_center)
	new("oct", imagefile=x, imagedata=imagedata, model=model, fitted_model=fitted_model)
})


plot.oct <- function(oct) {
		plot(0)
		#image(t(oct$fitted_model$centered_band-median(oct$fitted_model$centered_band)), col=topo.colors(10))
		#lines(e1/dim(oct$fitted_model$centered_band)[1], lwd=2, col="blue")
		#lines(e2/dim(oct$fitted_model$centered_band)[1], lwd=2, col="blue")
		#matplot(oct$fitted_model$fitted_center$max_slope[,1]/dim(oct$fitted_model$centered_band), oct$fitted_model$fitted_center$max_slope[,2]/dim(oct$fitted_model$centered_band), col="red", add=T, pch="+", cex=0.5)
		#matplot(oct$fitted_model$fitted_center$plateu[,1]/dim(oct$fitted_model$centered_band), oct$fitted_model$fitted_center$plateu[,2]/dim(oct$fitted_model$centered_band), col="red", add=T, pch="+", cex=0.5)
}





#setMethod("read", signature(x="oct", y="character"), function(x, y, ...) {
#	x$imagedata = dicomInfo(y)
#})




#setGeneric("oct", function(x, ...) standardGeneric("oct"))


#setClass("oct", slots=list(imagefile="character", scandata="matrix")) 
#setGeneric("oct", function(x, ...) { 
	#standardGeneric("oct")
#	.self$imagefile = x
#})

#setMethod("oct", "character", function(x, ...) {
#	.self$imagefile = x
#	.self$scandata = dicomInfo(x)
#})

#s <- new("oct",imagefile="")
