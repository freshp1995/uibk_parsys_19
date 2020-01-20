use Random;
use Time;
use LinearAlgebra;

/*
  Here we set up several config consts that represent, in order:
    - ``N``, the dimension for the square array ``A``
    - ``K``, the second dimension for arrays ``X`` and ``B``
    - ``seed``, a seed for random number generation
*/
config const N = 200;//2552;
config const seed = 100;
config const seed2 = 105;

// variables for the time taken by the program.
var t:Timer;

// Create the arrays ``A``, ``X``, and ``B``. Fill ``A`` and ``X`` with random
var A, X, B : [1..N, 1..N] real;  // first matrix
fillRandom(A, seed);

// Matrix multiply ``A*X``, store result in ``B``
proc PMmultiply(N: int, K: int): real {

  t.start();
  forall i in 1..N do {
    forall j in 1..K do {
      forall k in 1..N do {
         B[i,j] += A[i,k] * X[j,k];
      }
    }
  }
  t.stop();
  var elapsed = t.elapsed();
  t.clear();

  return elapsed;
}


// loop to print results of 10 executions of Mmultiply proc. 
var sum: real = 0;
for f in 1..10 do {
    sum += PMmultiply(N, N);
    writeln("Step ", f, "of", 10);
}



writeln(sum / 10, "Elapsed Time");