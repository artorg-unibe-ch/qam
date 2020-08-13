# qam
Library to compute 3D surface-distances between max. 3 binary segmentation images.



You can find more about ablation treatments for liver cancer in the following open access book chapter ["Stereotactic Image-Guidance for Ablation of Malignant Liver Tumors"](https://www.intechopen.com/online-first/stereotactic-image-guidance-for-ablation-of-malignant-liver-tumors) from Liver Pathology.


![image](https://user-images.githubusercontent.com/20581812/88671582-c38bdc80-d0e6-11ea-9dc8-325ca8a878c1.png)


We used an adapted Dice Soerensen coefficient to evaluate the ablation completness in treating liver tumors. The quantitative ablation margin (QAM) calculation is illustrated  in the next figure. The output is an array of 3D surface distances that can also be visualized as traffic-light colored histogram (see last step in the pipeline).
![image](https://user-images.githubusercontent.com/20581812/88669690-80306e80-d0e4-11ea-8524-494217bd58f5.png)

To save the contours and the distance map for overlaying on the (segmentation images), the affine transform must be the same as the source image and crop must be set to False.
Example:

    surface_distance = compute_distances(mask_gt=tumor_np, mask_pred=ablation_np,
                                         exclusion_zone=liver_np if has_liver_segmented else None,
                                         spacing_mm=spacing, connectivity=1, crop=False)

    if surface_distance['distances_gt_to_pred'].size > 0:
         borders_gt = nib.Nifti1Image(surface_distance["borders_gt"].astype(np.uint8), affine=tumor.affine)
         nib.save(borders_gt, "borders_gt.nii.gz")
The "borders_gt.nii.gz" contour (purple) can be overlayed over the tumor segmentation (white circle): 

![image](https://user-images.githubusercontent.com/20581812/90117807-a8140900-dd57-11ea-8d9b-200fbd3146c5.png)
