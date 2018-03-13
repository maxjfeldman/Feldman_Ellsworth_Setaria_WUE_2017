#!/usr/bin/env python

import sys, traceback
import cv2
import os
import re
import numpy as np
import argparse
import string
import plantcv as pcv


def options():
    parser = argparse.ArgumentParser(description="Imaging processing with opencv")
    parser.add_argument("-i", "--image", help="Input image file.", required=True)
    parser.add_argument("-o", "--outdir", help="Output directory for image files.", required=False)
    parser.add_argument("-r","--result", help="result file.", required= False )
    parser.add_argument("-r2","--coresult", help="result file.", required= False )
    parser.add_argument("-w","--writeimg", help="write out images.", default=False, action="store_true")
    parser.add_argument("-D", "--debug", help="Turn on debug, prints intermediate images.", action="store_true")
    args = parser.parse_args()
    return args

### Main pipeline
def main():
  # Get options
  args = options()
  
  # Read image
  img, path, filename = pcv.readimage(args.image)
    
  # Pipeline step
  device = 0

  # Convert RGB to HSV and extract the Saturation channel
  device, s = pcv.rgb2gray_hsv(img, 's', device, args.debug)
  
  # Threshold the Saturation image
  device, s_thresh = pcv.binary_threshold(s, 36, 255, 'light', device, args.debug)
  
  # Median Filter
  device, s_mblur = pcv.median_blur(s_thresh, 0, device, args.debug)
  device, s_cnt = pcv.median_blur(s_thresh, 0, device, args.debug)
  
  # Fill small objects
  #device, s_fill = pcv.fill(s_mblur, s_cnt, 0, device, args.debug)
  
  # Convert RGB to LAB and extract the Blue channel
  device, b = pcv.rgb2gray_lab(img, 'b', device, args.debug)
  
  # Threshold the blue image
  device, b_thresh = pcv.binary_threshold(b, 137, 255, 'light', device, args.debug)
  device, b_cnt = pcv.binary_threshold(b, 137, 255, 'light', device, args.debug)
  
  # Fill small objects
  #device, b_fill = pcv.fill(b_thresh, b_cnt, 10, device, args.debug)
  
  # Join the thresholded saturation and blue-yellow images
  device, bs = pcv.logical_and(s_mblur, b_cnt, device, args.debug)
  
  # Apply Mask (for vis images, mask_color=white)
  device, masked = pcv.apply_mask(img, bs, 'white', device, args.debug)
  
  # Convert RGB to LAB and extract the Green-Magenta and Blue-Yellow channels
  device, masked_a = pcv.rgb2gray_lab(masked, 'a', device, args.debug)
  device, masked_b = pcv.rgb2gray_lab(masked, 'b', device, args.debug)
  
  # Threshold the green-magenta and blue images
  device, maskeda_thresh = pcv.binary_threshold(masked_a, 127, 255, 'dark', device, args.debug)
  device, maskedb_thresh = pcv.binary_threshold(masked_b, 128, 255, 'light', device, args.debug)
  
  # Join the thresholded saturation and blue-yellow images (OR)
  device, ab = pcv.logical_or(maskeda_thresh, maskedb_thresh, device, args.debug)
  device, ab_cnt = pcv.logical_or(maskeda_thresh, maskedb_thresh, device, args.debug)
  
  # Fill small noise
  device, ab_fill1 = pcv.fill(ab, ab_cnt, 2, device, args.debug)
  
  # Dilate to join small objects with larger ones
  device, ab_cnt1=pcv.dilate(ab_fill1, 3, 2, device, args.debug)
  device, ab_cnt2=pcv.dilate(ab_fill1, 3, 2, device, args.debug)
  
  # Fill dilated image mask
  device, ab_cnt3=pcv.fill(ab_cnt2,ab_cnt1,150,device,args.debug)
  img2 = np.copy(img)
  device, masked2 = pcv.apply_mask(img2, ab_cnt3, 'white', device, args.debug)
  
  # Convert RGB to LAB and extract the Green-Magenta and Blue-Yellow channels
  device, masked2_a = pcv.rgb2gray_lab(masked2, 'a', device, args.debug)
  device, masked2_b = pcv.rgb2gray_lab(masked2, 'b', device, args.debug)
  
  # Threshold the green-magenta and blue images
  device, masked2a_thresh = pcv.binary_threshold(masked2_a, 127, 255, 'dark', device, args.debug)
  device, masked2b_thresh = pcv.binary_threshold(masked2_b, 128, 255, 'light', device, args.debug)
  device, ab_fill = pcv.logical_or(masked2a_thresh, masked2b_thresh, device, args.debug)
  
  # Identify objects
  device, id_objects,obj_hierarchy = pcv.find_objects(masked2, ab_fill, device, args.debug)
  
  # Define ROI
  device, roi1, roi_hierarchy= pcv.define_roi(masked2,'rectangle', device, None, 'default', args.debug,True, 650, 10,-575,-415)
  
  # Decide which objects to keep
  device,roi_objects, hierarchy3, kept_mask, obj_area = pcv.roi_objects(img,'partial',roi1,roi_hierarchy,id_objects,obj_hierarchy,device, args.debug)
  
  # Object combine kept objects
  device, obj, mask = pcv.object_composition(img, roi_objects, hierarchy3, device, args.debug)
  
  ############## Landmarks    ################
  device, points = pcv.acute_vertex(obj, 40, 40, 40, img, device, args.debug)  
  boundary_line = 405
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
    #outfile=args.outdir+"/"+filename
  
  # Find shape properties, output shape image (optional)
  device, shape_header,shape_data,shape_img = pcv.analyze_object(img, args.image, obj, mask, device,args.debug,outfile)
  
  # Shape properties relative to user boundary line (optional)
  device, boundary_header,boundary_data, boundary_img1= pcv.analyze_bound(img, args.image,obj, mask, 405, device,args.debug,outfile)
  
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

