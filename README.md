# 3DMM_edges
This is a Matlab implementation of an algorithm for fully automatically fitting a 3D Morphable Model to a single image using landmarks and edge features.

Please note, this is in development and we are in the process of uploading the required files. The current version should run and the provided demo.m will execute the basic algorithm. Some of the additional scripts and functions related to using other morphable models are not yet in the repository.

## References

If you use this code in your research, you should cite [the following paper](http://arxiv.org/abs/1602.01125):

[DOI](http://dx.doi.org/10.1007/978-3-319-54427-4_28) A. Bas, W.A.P. Smith, T. Bolkart and S. Wuhrer. "Fitting a 3D Morphable Model to Edges: A Comparison Between Hard and Soft Correspondences". In Proc. ACCV Workshop on Facial Informatics, LNCS vol. 10117, pp. 377-391, 2016.

and (for the landmark detector):

X. Zhu and D. Ramanan. "Face detection, pose estimation and landmark localization in the wild" in Proceedings of IEEE Conference on Computer Vision and Pattern Recognition (CVPR), 2012.

## Running the code

The easiest way to run the code is to run the demo.m script. This loads the model, loads an image, runs the landmark detector and then calls the edge fitting code.

To recreate the quantitative experiment in our paper, run the evaluation.m script. Again, you will need to set the base path to the Basel Face Model. The script first renders 9 views of each of the 10 out of sample faces. It then fits the model to each image, computes errors and stores these in a matrix.

## Dependencies

In order to use this code, you need to provide your own 3D Morphable Model. One such model (and the one we used while developing the code) is the [Basel Face Model](http://faces.cs.unibas.ch/bfm/?nav=1-0&id=basel_face_model). This model is freely available upon signing a license agreement. If you use the Basel Face Model, then all you need to do is set the base path to your model in the demo.m file:

```matlab
BFMbasedir = '...'; % Set this to your Basel Face Model base directory
```

If you wish to use a different morphable model, this should be fine but you will need to follow these steps:

1. Your model must provide four variables:
  * shapePC is a 3n by k matrix where n is the number of model vertices and k the number of principal components
  * shapeMU is a 3n by 1 vector containing the vertices of the mean shape
  * shapeEV is a k by 1 vector containing the sorted standard deviations of each principal component (note: standard deviations not variances as the BFM name would imply)
  * tl is an f by 3 matrix containing the face list for the model
2. You need to precompute two structures that allow fast lookup of edges adjacent to vertices and faces. You should save the two structures since they are fixed for a given triangulation. To do so, follow these two steps:
 * First compute the edge/vertex list by doing:
 ```matlab
 TR = triangulation(tl,ones(k,1),ones(k,1),ones(k,1));
 Ev = TR.edges;
 clear TR;
 ```
 * Second, use the provided script to compute the edge/face list by doing:
 ```matlab
 Ef = meshFaceEdges(tl,Ev);
 ```
3. You need to provide the morphable model indices that correspond to the output of the landmark detector (see below). This is done for the Zhu and Ramanan detector and Basel model in the function ZR2BFM.m function.

The code also requires a landmark detector. We provide the Zhu and Ramanan detector (see license below) but it would be easy to modify the code to use another detector. You may need to compile the mex files for your platform (we include precompiled mex files for win64 and Mac OS X).

## Third party licenses

This repository ships with a copy of the [Zhu and Ramanan facial feature detector](https://www.ics.uci.edu/~xzhu/face/), which was released under the following license:

> Copyright (C) 2012 Xiangxin Zhu, Deva Ramanan
> 
> Permission is hereby granted, free of charge, to any person obtaining
> a copy of this software and associated documentation files (the
> "Software"), to deal in the Software without restriction, including
> without limitation the rights to use, copy, modify, merge, publish,
> distribute, sublicense, and/or sell copies of the Software, and to
> permit persons to whom the Software is furnished to do so, subject to
> the following conditions:
>
> The above copyright notice and this permission notice shall be
> included in all copies or substantial portions of the Software.
>
> THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
> EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
> MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
> NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
> LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
> OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
> WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
