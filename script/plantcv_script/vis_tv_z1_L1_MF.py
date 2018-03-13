#!/usr/bin/env python
import sys, traceback
import cv2
import numpy as np
import argparse
import string
import plantcv as pcv

### Parse command-line arguments
def options():
  parser = argparse.ArgumentParser(description="Imaging processing with opencv")
  parser.add_argument("-i", "--image", help="Input image file.", required=True)
  parser.add_argument("-m", "--roi", help="Input region of interest file.", required=False)
  parser.add_argument("-o", "--outdir", help="Output directory for image files.", required=True)
  parser.add_argument("-r","--result", help="result file.", required= False )
  parser.add_argument("-D", "--debug", help="Turn on debug, prints intermediate images.", action="store_true")
  args = parser.parse_args()
  return args

### Main pipeline
def main():
  # Get options
  args = options()
  path_mask = '/home/mfeldman/tester/mask/mask_brass_tv_z1_L1.png'
  
  # Read image
  img, path, filename = pcv.readimage(args.image)
  brass_mask = cv2.imread(path_mask)
  
  # Pipeline step
  device = 0

  # Mask pesky brass piece
  device, brass_mask1 = pcv.rgb2gray_hsv(brass_mask, 'v', device, args.debug)
  device, brass_thresh = pcv.binary_threshold(brass_mask1, 0, 255, 'light', device, args.debug)
  device, brass_inv=pcv.invert(brass_thresh, device, args.debug)
  device, masked_image = pcv.apply_mask(img, brass_inv, 'white', device, args.debug)
  
  # We can do a pretty good job of identifying the plant from the s channel
  device, s = pcv.rgb2gray_hsv(masked_image, 's', device, args.debug)
  s_thresh = cv2.inRange(s, 100, 190)
  
  # Lets blur the result a bit to get rid of unwanted noise
  s_blur = cv2.medianBlur(s_thresh,5)
  
  # Lets mask out that cart from the s channel using the b channel
  device, b = pcv.rgb2gray_lab(masked_image, 'b', device, args.debug)
  b_top = cv2.inRange(b, 0, 115)
  b_top_inv = cv2.bitwise_not(b_top)
  s_masked = cv2.bitwise_and(s_blur, b_top_inv)

  # The a channel is good too
  device, a = pcv.rgb2gray_lab(masked_image, 'a', device, args.debug)
  a_thresh = cv2.inRange(a, 100, 120)
  a_blur = cv2.medianBlur(a_thresh,5)
  
  # Now lets set of a series of filters to remove unwanted background
  plant_shape = cv2.bitwise_and(a_blur, s_masked)
  
  # Lets remove all the crap on the sides of the image
  plant_shape[:,:330] = 0
  plant_shape[:,2100:] = 0
  plant_shape[:200,:] = 0
  
  # Now remove all remaining small points using erosion with a 3 x 3 kernel
  kernel = np.ones((3,3),np.uint8)
  erosion = cv2.erode(plant_shape ,kernel,iterations = 1)
  
  # Now dilate to fill in small holes
  kernel = np.ones((3,3),np.uint8)
  dilation = cv2.dilate(erosion ,kernel,iterations = 1)
  
  # Apply mask to the background image
  device, masked = pcv.apply_mask(masked_image, dilation, 'white', device, args.debug)
  
  # Identify objects
  device, id_objects, obj_hierarchy = pcv.find_objects(masked, dilation, device, args.debug)
  
  # Get ROI contours
  device, roi, roi_hierarchy = pcv.define_roi(masked_image, 'circle', device, None, 'default', args.debug, True, x_adj=0, y_adj=0, w_adj=0, h_adj=-1200)
  
  # ROI
  device,roi_objects, hierarchy3, kept_mask, obj_area = pcv.roi_objects(masked_image,'partial',roi, roi_hierarchy, id_objects,obj_hierarchy,device, args.debug)
  
  # Get object contour and masked object
  device, obj, mask = pcv.object_composition(img, roi_objects, hierarchy3, device, args.debug)
  
  ############## Landmarks    ################
  
  device, points = pcv.acute_vertex(obj, 40, 40, 40, img, device, args.debug)
  boundary_line = 'NA'
  # Use acute fxn to estimate tips
  device, points_r, centroid_r, bline_r = pcv.scale_features(obj, mask, points, boundary_line, device, args.debug)
    # Get number of points
  tips = len(points_r)
  # Use turgor_proxy fxn to get distances 
  device, vert_ave_c, hori_ave_c, euc_ave_c, ang_ave_c, vert_ave_b, hori_ave_b, euc_ave_b, ang_ave_b = pcv.turgor_proxy(points_r, centroid_r, bline_r, device, args.debug)
  # Get pseudomarkers along the y-axis
  device, left, right, center_h = pcv.y_axis_pseudolandmarks(obj, mask, img, device, args.debug)
  # Re-scale the points
  device, left_r, left_cr, left_br = pcv.scale_features(obj, mask, left, boundary_line, device, args.debug)
  device, right_r, right_cr, right_br = pcv.scale_features(obj, mask, right, boundary_line, device, args.debug)
  device, center_hr, center_hcr, center_hbr = pcv.scale_features(obj, mask, center_h, boundary_line, device, args.debug)
  
  # Get pseudomarkers along the x-axis
  device, top, bottom, center_v = pcv.x_axis_pseudolandmarks(obj, mask, img, device, args.debug)
  
  # Re-scale the points
  device, top_r, top_cr, top_br = pcv.scale_features(obj, mask, top, boundary_line, device, args.debug)
  device, bottom_r, bottom_cr, bottom_br = pcv.scale_features(obj, mask, bottom, boundary_line, device, args.debug)
  device, center_vr, center_vcr, center_vbr = pcv.scale_features(obj, mask, center_v, boundary_line, device, args.debug)
  
  ## Need to convert the points into a list of tuples format to match the scaled points
  points = points.reshape(len(points),2)
  points = points.tolist()
  temp_out = []
  for p in points:
    p = tuple(p)
    temp_out.append(p)
  points = temp_out
  left = left.reshape(20,2)
  left = left.tolist()
  temp_out = []
  for l in left:
    l = tuple(l)
    temp_out.append(l)
  left = temp_out
  right = right.reshape(20,2)
  right = right.tolist()
  temp_out = []
  for r in right:
    r = tuple(r)
    temp_out.append(r)
  right = temp_out
  center_h = center_h.reshape(20,2)
  center_h = center_h.tolist()
  temp_out = []
  for ch in center_h:
    ch = tuple(ch)
    temp_out.append(ch)
  center_h = temp_out
  ## Need to convert the points into a list of tuples format to match the scaled points
  top = top.reshape(20,2)
  top = top.tolist()
  temp_out = []
  for t in top:
    t = tuple(t)
    temp_out.append(t)
  top = temp_out
  bottom = bottom.reshape(20,2)
  bottom = bottom.tolist()
  temp_out = []
  for b in bottom:
    b = tuple(b)
    temp_out.append(b)
  bottom = temp_out
  center_v = center_v.reshape(20,2)
  center_v = center_v.tolist()
  temp_out = []
  for cvr in center_v:
    cvr = tuple(cvr)
    temp_out.append(cvr)
  center_v = temp_out
  
  #Store Landmark Data
  landmark_header=(
    'HEADER_LANDMARK',
    'tip_points',
    'tip_points_r',
    'centroid_r',
    'baseline_r',
    'tip_number',
    'vert_ave_c',
    'hori_ave_c',
    'euc_ave_c',
    'ang_ave_c',
    'vert_ave_b',
    'hori_ave_b',
    'euc_ave_b',
    'ang_ave_b',
    'left_lmk',
    'right_lmk',
    'center_h_lmk',
    'left_lmk_r',
    'right_lmk_r',
    'center_h_lmk_r',
    'top_lmk',
    'bottom_lmk',
    'center_v_lmk',
    'top_lmk_r',
    'bottom_lmk_r',
    'center_v_lmk_r'
    )

  landmark_data = (
    'LANDMARK_DATA',
    points,
    points_r,
    centroid_r,
    bline_r,
    tips,
    vert_ave_c,
    hori_ave_c,
    euc_ave_c,
    ang_ave_c,
    vert_ave_b,
    hori_ave_b,
    euc_ave_b,
    ang_ave_b,
    left,
    right,
    center_h,
    left_r,
    right_r,
    center_hr,
    top,
    bottom,
    center_v,
    top_r,
    bottom_r,
    center_vr
    )
    
  
  
 ############## VIS Analysis ################
  
  outfile=False
  #if args.writeimg==True:
  #  outfile=args.outdir+"/"+filename
  
  # Find shape properties, output shape image (optional)
  device, shape_header,shape_data,shape_img = pcv.analyze_object(img, args.image, obj, mask, device,args.debug,outfile)
  
  # Shape properties relative to user boundary line (optional)
  device, boundary_header,boundary_data, boundary_img1= pcv.analyze_bound(img, args.image,obj, mask, 330, device,args.debug,outfile)
  
  # Determine color properties: Histograms, Color Slices and Pseudocolored Images, output color analyzed images (optional)
  device, color_header,color_data,color_img= pcv.analyze_color(img, args.image, mask, 256, device, args.debug,None,'v','img',300,outfile)
  
  # Output shape and color data

  result=open(args.result,"a")
  result.write('\t'.join(map(str,shape_header)))
  result.write("\n")
  result.write('\t'.join(map(str,shape_data)))
  result.write("\n")
  for row in shape_img:
    result.write('\t'.join(map(str,row)))
    result.write("\n")
  result.write('\t'.join(map(str,color_header)))
  result.write("\n")
  result.write('\t'.join(map(str,color_data)))
  result.write("\n")
  result.write('\t'.join(map(str,boundary_header)))
  result.write("\n")
  result.write('\t'.join(map(str,boundary_data)))
  result.write("\n")
  result.write('\t'.join(map(str,boundary_img1)))
  result.write("\n")
  for row in color_img:
    result.write('\t'.join(map(str,row)))
    result.write("\n")
  result.write('\t'.join(map(str,landmark_header)))
  result.write("\n")
  result.write('\t'.join(map(str,landmark_data)))
  result.write("\n")
  result.close()
  
if __name__ == '__main__':
  main()