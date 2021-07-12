# main changepoint detection algorithm
# single changepoint model
step4 <- function(d) {
    n = length(d)
    dbar = mean(d)
    dsbar = mean(d*d)
    fac = dsbar-(dbar^2)
    summ = sum(d)
    summup = cumsum(d)
    y = vector()
    for (m in 1:(n-1)) {
        pos=m+1
        mscale = 4*(pos)*(n-pos)
        Q = summup[m]-(summ-summup[m])
        U = -(dbar*(n-2*pos) + Q)^2/mscale + fac
        y[m] = (-(n/2-1)*log(n*U/2) - 0.5*log((pos*(n-pos))))
    }
return (c(NaN,y))
}

example_step4 <- function() {
    d = c(6,5,8,5,5,6,4,6,8,7,8,6,9,5,7,7,9,4,10,11,8,5,5,7,6,5,8,12,8,
        11,7,6,4,4,8,12,19,9,9,4,10,9,12,7,6,7,6,11,9,7,12,5,10,5,5,11,5,10,11,
        5,6,7,6,5,6,15,4,6,5,5,7,5,9,11,4,11,7,6,3,2,9,11,21,11,6,16,9,8,5,8,6,
        10,10,11,11,6,9,9,13,13,10,7,14,11,12,8,8,10,11,7,11,5,9,12,6,11,8,8,
        18,10,13,7,8,11,9,10,6,12,9,11,12,6,4,13,11,6,7,12,11,14,8,7,8,7,19,
        13,8,14,12,15,15,12,16,15,28,10,15,16,9,19,19,14,8)
    years = 1851:2013
    step_like = step4(d)
    par(mfrow=c(2,1))
    plot(years, d)
    plot(years, step_like)
}
