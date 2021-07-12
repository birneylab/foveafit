# plotting function creating a pdf of the predicted max dip slice across a directory of dicom files
# @ visulisation
create_max_foveal_dip_pdf <- function(direct, outname, smooth_factor = 0.2) {
	pdf(outname, onefile=T)
	files = dir(direct)
	for(x in 1:length(files)) {
		image_file = paste(direct, "/", files[x], sep="")
		image_data = dicomInfo(image_file)
		fpos = get_fovea_position(image_data$img)
		slices = get_fovea_slices(image_data$img, fpos[1])
		slice = slices$midpoint
		fovea_slice = t(image_data$img[,,slice])
		x_fovea_params = fit_x_fovea(fovea_slice)
		y_fovea_params = fit_y_fovea(fovea_slice, x_fovea_params)
		dip = fit_sholders(fovea_slice, x_fovea_params, y_fovea_params)
		plot_fitted_slice(fovea_slice, x_fovea_params, y_fovea_params, dip)
	}
	dev.off()
}
