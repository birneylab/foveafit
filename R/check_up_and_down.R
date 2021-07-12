# function to find the switch point positions 
# where values change from accending to deceeding and visa versa 
# @ curve modelling
check_up_and_down <- function(l) {
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

# iteratively reduces smoothness is less than cst switch points exist
# @ curve modelling
up_down_reduce <- function(l, r, cst=3, fac=0.2) {
	up_down = check_up_and_down(l)
	while(length(up_down)<cst) {
		fac = fac-0.01
		l = lowess(apply(r, 2, mean), f=fac)
		up_down = check_up_and_down(l)
	}
return(up_down)
}

# function to find the switch point positions 
# also return which type of switch acc decc 
# @ curve modelling
check_up_and_down_which <- function(l) {
	check = 1
	up_down = NULL
	for(x in 2:length(l$x)) {
		if(check==1 & l$y[x]<l$y[x-1]) {
			check = 2
			up_down = rbind(up_down, cbind(check, l$x[x]-1))
		}
		if(check==2 & l$y[x]>l$y[x-1]) {
			check = 1
			up_down = rbind(up_down, cbind(check, l$x[x]-1))
		}
	}
return(up_down)
}

# robust implementation for finding the most likely max foveal dip slice
# @ fovea slice
get_fovea_slice_robust <- function(l, lwb=42, upb=85) {
	uds = check_up_and_down_which(l)
	uds = uds[uds[,1]==1 & uds[,2]>lwb & uds[,2]<upb,2]
	center = NULL
	if(length(uds)[1]>0) {
		center = uds[which.min(abs((uds/length(l$x))-0.5))]
	}
return(ifelse(length(center)>0, center, length(l$x)/2))
}

