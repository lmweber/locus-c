#!/bin/bash

##########################################################
# Code to downsample image using ImageMagick (if required)
##########################################################

# show details about image file, including size (width x height) in pixels
identify filename.jpg

# downsample image (where '-resize 5000' is the new image width in pixels)
convert filename.jpg -resize 5000 filename_new.jpg

# check
identify filename_new.jpg

