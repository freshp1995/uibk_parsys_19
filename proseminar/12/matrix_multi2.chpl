use LinearAlgebra;
use Random;
use Time;


config const N = 1000;//2552;
config const seed = 100;
config const seed2 = 105;

var A, X: [0..#N, 0..#N] real;

var t:Timer;

fillRandom(A, seed);
fillRandom(X, seed2);

proc PMmultiply(): real {
    t.start();
    var MM = dot(A, X.T);
    t.stop();
    var elapsed = t.elapsed();
    t.clear();

    return elapsed;
}

var sum: real = 0;
for f in 1..10 do {
    sum += PMmultiply();
    writeln("Step ", f, "of", 10);
}

writeln(sum / 10, "Elapsed Time");