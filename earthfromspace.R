# Earth seen from space: Moon and ISS
# www.overfitting.net
# https://www.overfitting.net/2023/01/proyeccion-lineal-de-escenas-3d-sobre.html


library(ggmap)  # map_data()
library(data.table)


# BITMAP DRAWING FUNCTIONS

NewBitmap = function(dimx, dimy, val=0) {
    # Crea bitmap de dimensiones dimx y dimy
    return(array(val,c(dimx,dimy)))
}

DrawCircle = function(img, x0, y0, r, inc=T, val=1, fill=F, thick=1) {
    # Dibuja círculo de centro (x0,y0) y radio r
    # Por defecto método no destructivo, con valor=1 y sin relleno
    # Puede elegirse el grosor si no se rellena
    img=DrawEllip(img, x0, y0, r, r, inc, val, fill, thick)
    
    return(img)
}

DrawEllip = function(img, x0, y0, a, b, inc=T, val=1, fill=F, thick=1) {
    # Dibuja elipse de centro (x0,y0) y radios a y b
    # Por defecto método no destructivo, con valor=1 y sin relleno
    # Puede elegirse el grosor si no se rellena
    # Aquí no redondeamos para tener más precisión en la división
    if (fill) {
        indices=which( ((row(img)-x0)/a)^2 + ((col(img)-y0)/b)^2 < 1 )
    } else {
        indices=which( ((row(img)-x0)/(a+thick/2))^2 + ((col(img)-y0)/(b+thick/2))^2 <  1 &
                           ((row(img)-x0)/(a-thick/2))^2 + ((col(img)-y0)/(b-thick/2))^2 >= 1 )
    }
    if (inc) img[indices]=img[indices]+val
    else img[indices]=val
    
    return(img)
}

SaveBitmap = function(img, name, trunc=TRUE, gamma=1) {
    # Guarda bitmap en formato PNG
    # Solo si trunc=FALSE y la imagen excede de 1 se reescala a 1
    library(png)
    img[img<0]=0
    if (trunc) img[img>1]=1
    if (tolower(substr(name, nchar(name)-3, nchar(name))) != ".png") name=paste0(name,".png")
    writePNG(t(img[,ncol(img):1] / max(max(img),1))^(1/gamma), name)
}


# COORDINATES CONVERSION

polar2x = function(R, phi, theta) {R*cos(theta)*sin(phi)}
polar2y = function(R, phi, theta) {R*sin(theta)}
polar2z = function(R, phi, theta) {-R*cos(theta)*cos(phi)}


# DRAW EARTH

# Physical parameters
Rearth=6371.23  # Earth average radius (km)
Rmoon=1737.4  # Moon average radius (km)
dearthmoon=385000  # centre to centre Earth to Moon distance (km)
dobserver.iss=408  # ISS average altitude (km)
dobserver.moon=dearthmoon-Rearth-Rmoon  # observation point to Earth surface distance (km)


# Obtain and process coordinates
DT=data.table(map_data("world"))  # long/lat pairs for all countries
DT=DT[, .(num=.N), by=.(long, lat)]  # summarize to deduplicate points

# deg to rad conversion
DT$phi=DT$long*pi/180  # longitude
DT$theta=DT$lat*pi/180  # latitude

# polar to XYZ coordinates conversion
# NOTE: x and y don't depend on d, but z depends on d so z is nested
DT$x=polar2x(Rearth, DT$phi, DT$theta)
DT$y=polar2y(Rearth, DT$phi, DT$theta)

NIMG=100  # number of images
IMAGESIZE=800  # 512  # animation dimensions in pixels
TH=0.9  # allow border
t=1

for (d in dobserver.iss+seq(0, 1 , length.out=NIMG)^5*(dobserver.moon-dobserver.iss)) {
    DT$z=polar2z(Rearth, DT$phi, DT$theta) + Rearth+d  # Earth along Z axis
    # Distance from each map point to observation point (0,0,0)
    DT$dist=(DT$x^2+DT$y^2+DT$z^2)^0.5
    
    distmax=( (d+Rearth)^2 - Rearth^2 )^0.5  # max distance to visible points
    DTplot=DT[DT$dist<=distmax]  # keep only visible points
    
    # Draw globe map
    img=NewBitmap(IMAGESIZE, IMAGESIZE)
    NCOLDIV2=ncol(img)/2
    NROWDIV2=nrow(img)/2
    
    # Calculate focal length to fit Earth in the final image and FOV
    thetamax=acos(Rearth/(d+Rearth))
    f=min(NCOLDIV2,NROWDIV2)*
        (d+Rearth-Rearth*cos(thetamax))/(Rearth*sin(thetamax))*TH
    FOV=(pi-2*thetamax)*180/pi  # FOV in deg
    # https://www.scantips.com/lights/fieldofview.html
    # ISS FOV:  140º diagonal -> 7.87mm FF
    # Moon FOV: 1.9º diagonal -> 1300mm FF
    
    print(paste0("Calculating ", t, "/", NIMG, " with distance d=", d,
                 "km and FOV=", FOV, "º"))
    
    img=DrawCircle(img, NCOLDIV2, NROWDIV2, min(NCOLDIV2,NROWDIV2)*TH,
                   fill=TRUE, val=0.25)
    DTplot$factor=f/DTplot$z
    DTplot$xp=DTplot$x*DTplot$factor + NCOLDIV2
    DTplot$yp=DTplot$y*DTplot$factor + NROWDIV2
    img[round(cbind(DTplot$xp, DTplot$yp))]=1  # draw points7
    
    SaveBitmap(img, paste0("img", ifelse(t<10, "00", ifelse(t<100, "0", "")), t))
    t=t+1
}
