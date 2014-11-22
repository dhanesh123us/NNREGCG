NNREGCG
=======

A "dynamic" 3 layered neural net regression code using conjugate gradient algorithm for training
======================================================================================

Features of this code:

1. least squares cost function for error
2. a Polak-Ribiere CG training algorithm (with Powell restart strategy)
3. Ability to add some noise to the optimizer to prevent it from getting stuck in local optima
4. Weight reinitializations to random weights in case training error target is not met
5. Hidden neurons automatically increased after certain number of weight re-initializations

http://en.wikipedia.org/wiki/Nonlinear_conjugate_gradient_method

How to use:

Compile nnopt.C using g++. Only supported platform is Unix right now since a word count utility is run internally to count the number of data points.

```
Dhaneshs-MacBook-Air-9:NNREGCG dhanesh.padmanabhan$ g++ nnopt.C
Dhaneshs-MacBook-Air-9:NNREGCG dhanesh.padmanabhan$ ./a.out
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
Dhaneshs-MacBook-Air-9:NNREGCG dhanesh.padmanabhan$ 
```


I guess lot of improvements can be done to this code: 
1. generalizing this for classification problems, 
2. support for categorical variables, 
3. Adding a Stochastic Gradient Descent variant for large scale data, 
4. Better interface (programmatic as supposed to interactive)
5. Ability to read in weight & hidden bias files and NN architecture for prediction (scale parameters not stored currently)

