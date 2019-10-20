# Exercise 1

## Write a sequential application ```pi_seq``` in C or C++ that computes π for a given number of samples (command line argument). Test your application for various, large sample sizes to verify the correctness of your implementation.

Execute make file and ```pi_seq```.

Results are close to the value of π.

## Consider a parallelization strategy using MPI. Which communication pattern(s) would you choose and why?

Since we create random values, we only have to prallelize this part. We don't have to send any initialize data to other ranks. For gathering results, we used  ```MPI_Reduce``` with the ```MPI_SUM``` method.
