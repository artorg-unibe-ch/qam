# qam
Library to compute 3D surface-distances between max. 3 binary segmentation images.



You can find more about ablation treatments for liver cancer in the following open access book chapter ["Stereotactic Image-Guidance for Ablation of Malignant Liver Tumors"](https://www.intechopen.com/online-first/stereotactic-image-guidance-for-ablation-of-malignant-liver-tumors) from Liver Pathology.


![image](https://user-images.githubusercontent.com/20581812/88671582-c38bdc80-d0e6-11ea-9dc8-325ca8a878c1.png)


We used an adapted Dice Soerensen coefficient to evaluate the ablation completness in treating liver tumors. The quantitative ablation margin (QAM) calculation is illustrated  in the next figure. The output is an array of 3D surface distances that can also be visualized as traffic-light colored histogram (see last step in the pipeline).
![image](https://user-images.githubusercontent.com/20581812/88669690-80306e80-d0e4-11ea-8524-494217bd58f5.png)
