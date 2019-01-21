from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import math
from math import pi, cos



num_bins = 100
max_vel = 500000.
central_wavelength1 = 656.28 #nano metres
central_wavelength2 = 486.1
npixels = 250


MAX_LAMDA = 750.
MIN_LAMDA = 380.
HUE_MIN = 0.
HUE_MAX = 240.
c = 299792458.

def hue(wavelength):
    if(wavelength < MIN_LAMDA or wavelength > MAX_LAMDA):
        print("wavelenght out of bounds")
        exit()
    return int(HUE_MAX-((wavelength-MIN_LAMDA)*(HUE_MAX-HUE_MIN)/(MAX_LAMDA-MIN_LAMDA)))

def hsi_2_rgb(intensity, hue, saturation):
    if(hue > HUE_MAX):
        hue -= HUE_MAX
    h = hue * pi / 180.
    x = intensity * (1 - saturation)
    R=0.
    G=0.
    B=0.
    if(h < 2.*pi/3):
        y =  intensity * (1.+(saturation * cos(h)) / cos((pi / 3.) - h))
        z = (3. * intensity) - (x + y)
        R=y
        G=z
        B=x
    elif(h < (4./3.)*pi):
        h -= 2.*pi/3.
        y =  intensity * (1.+((saturation * cos(h)) / cos((pi / 3.) - h)))
        z = (3. * intensity) - (x + y)
        R=x
        G=y
        B=z
    elif(h < 2.*pi):
        h -= 4.*pi/3.
        y =  intensity * (1.+(saturation * cos(h)) / cos((pi / 3.) - h))
        z = (3. * intensity) - (x + y)
        R=z
        G=x
        B=y

    return(int(255*R/3.),int(255*G/3.),int(255*B/3.))

def vel_2_wavelength(velocity, lamda): #lamda is the central wavelength
    centre_freq = c / lamda
    return c/(centre_freq*((c + velocity) / (c - velocity)))

def normalize_point(x,min,max,new_min,new_max):
    return ((x-min)*(new_max-new_min))/(max-min)

def add_colour(c1,c2,weight):#c1 and c2 are 3 element turples

    """if(c1[0] == 0 and c1[1] == 0 and c1[2] == 0):
        weight = 1.
    if(c2[0] == 0 and c2[1] == 0 and c2[2] == 0):
        weight = 1."""
    R = (c1[0] + weight*c2[0])
    if(R > 255):
        R = 255
    G = (c1[1] + weight*c2[1])
    if(G > 255):
        G = 255
    B = (c1[2] + weight*c2[2])
    if(B > 255):
        B = 255

    return (int(R),int(G),int(B))

def create_image(filename,wavelength,numframes,output):
    with fits.open(filename) as data_file:
        data = data_file[0].data

    data_size = np.shape(data)
    max_int = np.amax(data)
    count = 0
    im = Image.new("RGB", (data_size[1], data_size[2]))
    pix = im.load()
    max_int = np.amax(data)


    for image_num in range(0,data_size[0],int(data_size[0]/numframes)):

        velocity = (image_num * (2*max_vel / num_bins)) - max_vel
        lamda = vel_2_wavelength(velocity,wavelength,)
        print("processing image {} at wavelength {}".format(image_num,lamda))

        for x in range(data_size[1]):
            for y in range(data_size[2]):
                norm_int =  normalize_point(data[image_num,x,y],0.,max_int,0.,1.,)
                r,g,b = hsi_2_rgb(norm_int,hue(lamda),1)
                pix[x,y] = add_colour(pix[x,y],(r,g,b),norm_int)

    im.rotate(90).save(output, "PNG")
    return im
##########################################################################

im1 = Image.new("RGB", (npixels,npixels))
im2 = Image.new("RGB", (npixels,npixels))


im1 = create_image("halpha_055_06563.fits",central_wavelength1,10,"ttauri1.png")
im2 = create_image("halpha_055_06563.fits",central_wavelength2,10,"ttauri2.png")

pix1 = im1.load()
pix2 = im2.load()

for x in range(npixels):
    for y in range(npixels):
        pix1[x,y] = add_colour(pix1[x,y],pix2[x,y],0.5)

im1.rotate(90).save("combined_50frame.png", "PNG")
