# EPFL IC, Mini-Project 1: Fingerprints Comparison

## Function getConnectedPixels:

A general idea of algorithm was proposed in the pdf file to implement for this function.
The first thing we realised when writing this method is that it has an undefinite number of iterations in a while loop.
Furthermore, its time complexity is O(n^2) because it loops over the whole image each time.

So we tried to write a more time-efficient algorithm.
At each iteration the process executes the following steps:

-   Compute the eight neighbours of the current position on the following two arrays: image and connectedPixels
-   If there is at least one black neighbour in image not in connectedPixels
    -   Check if its in distance, if false go to next unconnected pixel in the neighbours or if none left go to else statement
    -   Set the current position to that neighbour
    -   Add the previous position to a stack of old positions, the base one being the minutia.
-   Else:
    -   Pop the element of the stack last added
    -   Set the position to the top of the stack
-   While loop condition: If stack is empty exit.

In this implementation, we can see that we pass over each connected pixel twice, the first time when
added to the connection network, the second time when backtracking to find new ones.
So this algorithm doesn't depend much on the size of the image as the distance is a constant, so there
will be constantly many possible connected pixels.

## Return type in tests for function match

We changed some testing functions provided in the main class. We made them return boolean for testCompareFingerprints.
And in the other looping testing functions, they return int type which is the number of errors in total over the loop.
This makes the process of computing the correctness of the algorithm easier.
