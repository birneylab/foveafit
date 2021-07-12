triangle <- function() {
image_file = "/Users/currant/OCTeyes/OCT_30_replicates_20161024/1444558/OCT/679115.dcm"
model = run_generate_model_and_find_fovea_dip_poisition(image_file)


center_f = 0.6
amodel = center_align(model, center_f)

image(amodel)

#####makes a grid######
posrange <- seq(from = 1, to = -1.732, by = -0.1)
for (n in posrange){
  abline(a=n, b='1.732051')
}
negrange <- seq(from = 0, to = 2.732, by = 0.1)
for (n in negrange){
  abline(a=n, b='-1.732051')
}
horizontal1 <- seq(from = 0, to = 1, by=0.05)
for (n in horizontal1){
  abline(a=n, b=0)
}
for (n in 0.5:77){
  abline(a=n, b=0)
}

#Find veritces####
#positive equation: y = 1.732051x + n
#negative equation: y = -1.732051x + n
#positive equation intercept range  = -1.732:1
#negative equation intercept range = 0:2.732
posrange <- seq(from = 1, to = -1.732, by = -0.1)
posrange
posrange <- sort(posrange)
posrange
negrange <- seq(from = 0, to = 2.732, by = 0.1)
negrange 


vertices1 <- NULL
for (n in posrange){
  for (j in negrange){
    #if (n != j){
      coords1 <- NULL
      xval1 = (j - n)/3.464102
      coords1 <- append(coords1, xval1)
      yval1 = (1.732051*((j - n)/3.46102)) + n 
      coords1 <- append(coords1, yval1)
   # }
    vertices1 <- rbind(vertices1, coords1)
  }
  
}

vertices2 <- NULL
for (k in horizontal1){
  for (n in posrange){
    coords2 <- NULL
    xval2 <- (k - n)/1.732051
    coords2 <- append(coords2, xval2)
    yval2 <- (1.732051*((k - n)/1.732051)) + n
    coords2 <- append(coords2, yval2)
    vertices1 <- rbind(vertices1, coords2)
    #vertices2 <- rbind(vertices2, coords2)
  }
}

vertices3 <- NULL
for (k in horizontal1){
  for (j in negrange){
    coords3 <- NULL
    xval3 <- (k - j)/-1.732051
    coords3 <- append(coords3, xval3)
    yval3 <- (-1.732051*((k-j)/-1.732051)) + j
    coords3 <- append(coords3, yval3)
    #vertices3 <- rbind(vertices3, coords3)
    vertices1 <- rbind(vertices1, coords3)
  }
}
#######

alltrianglevertices <- NULL
for (n in 2:length(posrange)){
    for (j in 1:j){
      vertexofeach <- NULL
      vertexax <- (negrange[j] - posrange[n])/3.464102
      vertexofeach <- append(vertexofeach, vertexax)
      vertexay <- (1.732051*((negrange[j] - posrange[n])/3.464102)) + posrange[n]
      vertexofeach <- append(vertexofeach, vertexay)
      vertexbx <- (negrange[j+1] - posrange[n])/3.464102
      vertexofeach <- append(vertexofeach, vertexbx)
      vertexby <- (1.732051*((negrange[j+1] - n)/3.464102)) + posrange[n]
      vertexofeach <- append(vertexofeach, vertexby)
      vertexcx <- (negrange[j+1] - posrange[n-1])/3.464102
      vertexofeach <- append(vertexofeach, vertexcx)
      vertexcy <- (1.732051*((negrange[j+1] - posrange[n-1])/3.464102)) + posrange[n]
      vertexofeach <- append(vertexofeach, vertexcy) 
      alltrianglevertices <- rbind(alltrianglevertices, vertexofeach)
    }
}

alltrianglevertices_new <- NULL
for (n in 2:length(posrange)){
  for (j in 1:j){
    vertexofeach1 <- NULL
    vertexax <- (negrange[j] - posrange[n])/3.464102
    vertexofeach1 <- append(vertexofeach1, vertexax)
    vertexay <- (1.732051*((negrange[j] - posrange[n])/3.464102)) + posrange[n]
    vertexofeach1 <- append(vertexofeach1, vertexay)
    alltrianglevertices_new <- rbind(alltrianglevertices_new, vertexofeach1)
    vertexofeach2 <- NULL
    vertexbx <- (negrange[j+1] - posrange[n])/3.464102
    vertexofeach2 <- append(vertexofeach2, vertexbx)
    vertexby <- (1.732051*((negrange[j+1] - n)/3.464102)) + posrange[n]
    vertexofeach2 <- append(vertexofeach2, vertexby)
    alltrianglevertices_new <- rbind(alltrianglevertices_new, vertexofeach2)
    vertexofeach3 <- NULL
    vertexcx <- (negrange[j+1] - posrange[n-1])/3.464102
    vertexofeach3 <- append(vertexofeach3, vertexcx)
    vertexcy <- (1.732051*((negrange[j+1] - posrange[n-1])/3.464102)) + posrange[n]
    vertexofeach3 <- append(vertexofeach3, vertexcy) 
    alltrianglevertices_new <- rbind(alltrianglevertices_new, vertexofeach3)
  }
}















atriangle <- NULL

xandy <- alltrianglevertices[1, 1:2]
atriangle <- rbind(atriangle, xandy)

xandy1 <- alltrianglevertices[1, 3:4]
atriangle <- rbind(atriangle, xandy1)

xandy2 <- alltrianglevertices[1, 5:6]
atriangle <- rbind(atriangle,xandy2)




n = -1.7
j = 0
coords1 <- NULL
xval = (j - n)/3.464102
coords1 <- append(coords1, xval)
yval = (1.732051*((j - n)/3.46102)) + n 
coords1 <- append(coords1, yval)
yval2 = (-1.732051*((j - n)/3.46102)) + j
coords1 <- append(coords1, yval2)

coords1
########




















tes_f = 5
v1 = square_tes_extract(amodel, tes_f)

x <- -20:20
y <- x*1.732051


plot(y, type = 'l')
for (n in -534:77){
  abline(a=n, b='1.732051')
}
for (n in 0:610.4717){
  abline(a=n, b='-1.732051')
}
for (n in 0:77.5){
  abline(a=n, b=0)
}
for (n in 0.5:77){
  abline(a=n, b=0)
}
}
