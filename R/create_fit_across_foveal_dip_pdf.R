# plotting function creating a pdf of old fitting method across the foveal dip slices
# @ visulisation
create_fit_across_foveal_dip_pdf <- function(slices, outname, cut = 0.5) {
	pdf(outname, onefile=T)
	for(slice in slices$l_shoulder:slices$r_shoulder) {
		fovea_slice = t(image_data$img[,,slice])
		y_fovea_params = fit_y_fovea(fovea_slice, x_fovea_params)
		dip = fit_sholders(fovea_slice, x_fovea_params, y_fovea_params)
		plot_fitted_slice(fovea_slice, x_fovea_params, y_fovea_params, dip)
	}
	dev.off()
}
