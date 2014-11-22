NNREGCG
=======

##A "dynamic" 3 layered neural net regression code using conjugate gradient algorithm for training
======================================================================================

Caution: This code was used for an engineering design application during my PhD work nearly 15 years ago. I abandoned this work after a year and switched to other topics. Hence I neither cleaned this code up or optimized it for performance. Luckily the code does work but do use it at your own risk.

#Features of this code:
=====================

1. Sum of squared errors cost function
2. A Polak-Ribiere CG training algorithm (with Powell restart strategy) http://en.wikipedia.org/wiki/Nonlinear_conjugate_gradient_method
3. Reinitialization of Input-Hidden and Hidden-Output Weight to random values in case training error target is not met or error decay is very small
4. Hidden neurons automatically increased after certain weight weight re-initializations (currently hard coded at 50)
5. Training exits after 10,000 iterations (hard coded)


##How to use:
===========

Compile nnopt.C using g++. Only supported platform is Linux, Mac or Cygwin (for Windows) right now since a word count utility is run internally to count the number of data points. Below is an example run for fitting y=x^2 (sample points are given in file 'sq')

```
$ g++ nnopt.C
$ ./a.out
Input the file name
sq
Input the no. of input, output and determinant
1 1 1
Input the Min. & Max. Scaling Values resp.
0.1 0.9
Input the Noise Factor
0.2
Input the tolerance for line searches
0.001
Input the step size for bounding phase
0.5
Input the max. allowed training error
1e-12
Input the tolerance for adding noise
1e-9
1 210.25
Maximum no. of Iterations reached

Least global error obtained 0.000135455
1 pattern error 7.53815e-11
2 pattern error 2.15721e-09
3 pattern error 9.17356e-09
4 pattern error 7.71574e-09
5 pattern error 6.22023e-10
6 pattern error 4.80327e-08
7 pattern error 1.80104e-07
8 pattern error 3.10334e-07
9 pattern error 3.00632e-07
10 pattern error 1.33142e-07
11 pattern error 4.95689e-10
12 pattern error 1.70654e-07
13 pattern error 6.94364e-07
14 pattern error 1.23439e-06
15 pattern error 1.27919e-06
16 pattern error 6.79519e-07
17 pattern error 3.96919e-08
18 pattern error 4.08953e-07
19 pattern error 2.21664e-06
20 pattern error 4.28525e-06
21 pattern error 4.33653e-06
22 pattern error 1.67802e-06
23 pattern error 9.98668e-08
24 pattern error 5.79084e-06
25 pattern error 1.69038e-05
26 pattern error 1.6058e-05
27 pattern error 7.89808e-08
28 pattern error 7.85083e-05
------------------------------------------------
Max pattern error 7.85083e-05 pattern no.28
Min pattern error 7.53815e-11 pattern no.1
Input file name to store input-hidden weights
wih1.dat
Input file name to store hidden-output weights
who1.dat
Input file name to store hidden bias values
whb.dat
Input file name to store network arch.
nnarch.dat
Input file name to store network parms
nnparms.dat
Wanna predict? <type 'y' for yes>
y
Type the input
5
Predicted : 24.7972
Wanna predict? <type 'y' for yes>
y
Type the input
10
Predicted : 99.4493
Wanna predict? <type 'y' for yes>
y
Type the input
12
Predicted : 144.117
Wanna predict? <type 'y' for yes>
n
```

##Explanation of Input Arguments to be supplied:
============================================

1. Input the file name (needs to have data in x1,x2,x3,..,y format where y is the response variable - whitespace is the delimiter supported)
sq
2. Input the no. of input, output and determinant (determinant means the number of hidden neurons)
1 1 1
3. Input the Min. & Max. Scaling Values resp. (response variable is uniformly scaled based on min, max values)
0.1 0.9
4. Input the Noise Factor (this is currently not used - I had originally coded a feature to randomly perturb the optimizer to get out of local optima. Only random re-initialization is carried out currently)
0.2
5. Input the tolerance for line searches (this is exit criteria for line search - how small the step-size gets during the line search)
0.001
6. Input the step size for bounding phase (this is the step-size is used for bounding a minimum before line search is used)
0.5
7. Input the max. allowed training error (tolerance on 0.5* sum of squared error)
1e-12
8. Input the tolerance for adding noise (error decay tolerance that triggers a random initialization)
1e-9


##Improvements to this code
==========================

1. Ability to read in weight & hidden bias files and NN architecture for prediction (scale parameters not stored currently)
2. generalizing this for classification problems, 
3. support for categorical variables, 
4. Adding a Stochastic Gradient Descent variant for large scale data, 
5. Better interface (programmatic as supposed to interactive)
