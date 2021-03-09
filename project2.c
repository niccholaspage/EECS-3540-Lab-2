//==============================================================================================
// project2.c - Multithreaded hybrid sort program
//
// This program implements a hybrid sorting algorithm, utilizing
// both the quick sort and shell sort algorithms on an array, as well as
// POSIX threads to simultaneously sort multiple segments of an array
// at once. This can greatly improve the speed of sorting arrays with
// many elements depending on the commandline arguments passed into
// the program. The program takes multiple arguments, including the size
// of the array it will be sorting, the threshold of a segment before switching
// to shell sort, the seed for the random number generator that shuffles the array,
// whether to use multithreading or not, how many pieces to split the array into,
// and how many threads to run concurrently while sorting the pieces.
//
// Author:     Nicholas Nassar, University of Toledo
// Class:      EECS 3540-001 Operating Systems and Systems Programming, Spring 2021
// Instructor: Dr. Thomas
// Date:       March 8, 2021
// Copyright:  Copyright 2021 by Nicholas Nassar. All rights reserved.

// We define  _GNU_SOURCE so we can
// use the nonportable pthread_tryjoin_np
// method.
#define _GNU_SOURCE

// We love macros! Since macros let me write
// function-like calls that get inlined by the
// preprocessor, it is an incredible feature.

// We first define a macro that lets us pass
// in a range pointer and get its size,
// which is simply its high index minus
// its low index plus one.
#define getRangeSize(range) (range)->high - (range)->low + 1

// Next, we define a macro that lets us pass
// a start and end time, and returns the difference
// in seconds between the two times. We utilize the
// CLOCKS_PER_SEC definition instead of hardcoding
// its value to be more proper.
#define getTimeBetweenClocks(startTime, endTime) ((double)(endTime - startTime)) / CLOCKS_PER_SEC

// Finally, we define a macro that lets us pass
// a start and end time of day, and returns the
// difference in seconds between the two times.
// We simply get the difference between the seconds
// of the two times, then add it to the difference in
// microseconds after converting the microseconds into
// seconds.
#define getTimeBetweenTimeOfDays(startTimeOfDay, endTimeOfDay) (endTimeOfDay.tv_sec - startTimeOfDay.tv_sec) + (endTimeOfDay.tv_usec - startTimeOfDay.tv_usec) / 1e6

// We also define a constant representing the seed
// value when no seed is provided by the user. An
// important note is that this value DOES NOT actually
// get seeded with srand - it is just used so that
// our final output can show that we did not provide
// a seed, which we will output as 00.
#define NO_SEED -2

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <pthread.h>
#include <unistd.h>
#include <string.h>
#include <sys/time.h>

// A structure used to pass a range which
// includes a low and high index of our
// array to our worker threads. It is also
// to represent each piece of the array we
// create and sort.
typedef struct
{
    int low;
    int high;
} range;

// We setup our global variables. Since our sorting methods
// need access to the array and array size to sort segments
// of the array, we make them global. The threshold also needs
// to be global so that our hybridSort can properly call shell
// sort on a certain threshold.
int *array;
int arraySize;
int threshold;

// Method prototypes so the main method can be the first
// method with a body in my source code. For descriptions
// on each method, see the implementations of each method
// below.
void swap(int *a, int *b);
void shellSort(int low, int high);
int partition(int p, int r);
void shellSort(int low, int high);
void hybridSort(int low, int high);
void *runner(void *param);
bool isSorted();

int main(int argc, char *argv[])
{
    // The entry point of the program.

    // We start by getting the clock value
    // and storing it in our program start time,
    // so we can calculate the entire time it took
    // to run the program. We also get the current
    // time of day, so that we can calculate our
    // wall clock time as well.
    clock_t programStartTime = clock();
    struct timeval programStartTimeOfDay;
    gettimeofday(&programStartTimeOfDay, NULL);

    // We now handle the arguments the user can
    // pass into our program.
    if (argc < 3) // If the user doesn't pass in the size and threshold,
    {
        // we print out a line detailing the arguments for the
        // program, then immediately return.
        printf("./project2 SIZE THRESHOLD [SEED [MULTITHREAD [PIECES [THREADS]]]]\n");

        return -1;
    }

    // The first argument the user passes in (ignoring the executable path)
    // is the size they want the array to be, so we parse the argument as an
    // int and set our array size to it. atoi is kind of strange since it counts
    // certain strings with letters somewhere throughout the string (ex. 123ad)
    // as valid, and would return 123. Since we aren't worried about invalid
    // parameters for this project, we just don't worry about this.
    arraySize = atoi(argv[1]);

    if (arraySize <= 0) // If the array size is less or equal 0,
    {
        // we tell the user that the size must be positive
        // and then immediately return. We can't have an array
        // with a negative size, and an array with zero elements
        // doesn't need to be sorted.
        printf("The size must be a positive integer.\n");

        return -1;
    }

    // The second argument is the threshold at which our hybrid
    // sort switches to utilizing shell sort instead of quick
    // sort.
    threshold = atoi(argv[2]);

    if (threshold < 0) // If our threshold is less than zero,
    {
        // we tell the user our threshold must be zero or a positive
        // integer and return. Since our threshold is checked against
        // an array size and we can't have a negative array size, a
        // negative threshold doesn't make sense.
        printf("The threshold must be zero or a positive integer.\n");

        return -1;
    }

    // We setup the defaults for our
    // parameters.
    bool isMultithreaded = true; // We multithread by default
    int numberOfPieces = 10;     // We default to 10 pieces,
    int maxThreads = 4;          // and four max threads.
    int seed = NO_SEED;          // We default to no seed.

    if (argc >= 4) // If we have four or more arguments,
    {
        // then the user has provided a seed, so we
        // parse it as an integer.
        seed = atoi(argv[3]);

        if (seed < -1) // If our seed is less than -1,
        {
            // we error out and tell the user their seed must be greater than
            // or equal to one and return. srand takes unsigned integers, so
            // we need to make sure we seed srand with an unsigned int. -1 is
            // a special case, since we will seed srand with a call to clock().
            printf("You must specify a seed greater than or equal to -1.\n");

            return -1;
        }

        if (seed == -1) // If the seed is -1,
        {
            srand(clock()); //we seed our generator with the current clock value.
        }
        else
        {
            srand(seed); // otherwise, we seed our generator with the seed argument provided.
        }

        if (argc >= 5) // If we have five or more arguments,
        {
            char *multiThreadedArg = argv[4]; // we get the multithreaded argument.

            // If our multithreaded argument is n or N, we should not multithread.
            // Otherwise, we do multithread.
            if (strcmp(multiThreadedArg, "n") == 0 || strcmp(multiThreadedArg, "N") == 0)
            {
                isMultithreaded = false;
            }

            if (argc >= 6) // If we have six or more arguments,
            {
                // we get the number of peices by parsing the
                // argument as an integer.
                numberOfPieces = atoi(argv[5]);

                if (numberOfPieces < 1) // If we have less than one piece,
                {
                    // we print a message stating that the user must
                    // specify a positive number of pieces and return.
                    // You cannot split an array in a negative number
                    // of pieces, or zero pieces.
                    printf("You must specify a positive number of pieces.\n");

                    return -1;
                }

                if (argc >= 7) // If we have seven or more arguments,
                {
                    // we get the max number of threads by
                    // parsing the argument as an integer.
                    maxThreads = atoi(argv[6]);

                    // Again, we must have a positive number of threads.
                    // We can't spawn a negative number of threads, and
                    // with zero threads, we could not sort anything.
                    if (maxThreads < 1)
                    {
                        printf("You must specify a positive number of threads.\n");

                        return -1;
                    }
                }

                // If our max threads is greater than our number of pieces,
                if (maxThreads > numberOfPieces)
                {
                    // we print an error and return immediately. Because we
                    // don't have enough pieces to schedule work for each
                    // thread, we can't continue.
                    printf("THREADS cannot be bigger than PIECES.\n");

                    return -1;
                }
            }
        }
    }

    // To time our program, we declare
    // variables for start and end times,
    // as well as start and end time of
    // days, so we can calculate CPU time
    // as well as time of day.
    clock_t startTime, endTime;
    struct timeval startTimeOfDay, endTimeOfDay;

    // We also create variables to hold our
    // time calculations for array creation time,
    // array initialization/fill time, array shuffling
    // time, array partition time, and CPU time needed
    // to sort.
    double createTime, initTime, shuffleTime, partitionTime, sortCpuTime;

    // We first mark the start time, allocate memory
    // for our array of integers, mark the end time,
    // and calculate the time it took to perform the
    // allocation. To determine the correct number of
    // bytes to allocate, we multiply the given array
    // size by the size of bytes needed for a single int.
    startTime = clock();
    array = (int *)malloc(arraySize * sizeof(int));
    endTime = clock();
    createTime = getTimeBetweenClocks(startTime, endTime);

    // Next, we need to fill our array with elements.
    // We mark the start time, fill our array with
    // elements, mark the end time, then calculate the
    // time it took to fill the array.
    startTime = clock();
    for (int i = 0; i < arraySize; i++) // For each item in the array,
    {
        array[i] = i; // we set the item at index i of our array to the int i.
    }
    endTime = clock();
    initTime = getTimeBetweenClocks(startTime, endTime);

    // Finally, we shuffle the elements in our array.
    // We mark the start time, shuffle the elements of
    // our array, mark the end time, then calculate the
    // time it took to shuffle the array.
    startTime = clock();

    // We loop through every item in the array except
    // the last one.
    for (int i = 0; i < arraySize - 1; i++)
    {
        // To shuffle this particular element of the array,
        // we first get the index we will swap it with by
        // calculating a random number and adding i, then
        // swapping the elements at the two indices.
        int j = i + rand() / (RAND_MAX / (arraySize - i) + 1);
        swap(&array[j], &array[i]);
    }

    endTime = clock();
    shuffleTime = getTimeBetweenClocks(startTime, endTime);

    if (isMultithreaded) // If multithreading is on,
    {
        // we mark the start time so we can
        // calculate the time it takes to split
        // the array into multiple pieces.
        startTime = clock();

        // We now need to split our array
        // into multiple pieces so that
        // we can start threads to perform
        // sorting on these pieces of the
        // array.

        // We start with a current number of
        // pieces of 1, because we will start
        // with one piece, which begins at index
        // 0 of the array and ends at the last index
        // in the array.
        int currentNumberOfPieces = 1;

        // We setup our pieces array, which will
        // have the number of pieces specified
        // in the commandline arguments, or
        // 10 by default.
        range pieces[numberOfPieces];

        // We then grab the initial piece of the array
        // and set its low index to zero and its
        // right index to the last element index in the
        // array, as we described above.
        range *initialPiece = &pieces[0];
        initialPiece->low = 0;
        initialPiece->high = arraySize - 1;

        // While our current number of pieces is less than
        // the number of pieces we need to split the array
        // into,
        while (currentNumberOfPieces < numberOfPieces)
        {
            // we find the biggest piece in the array.
            // We need to find it because we will be
            // splitting it into two pieces.
            range *biggestPiece = &pieces[0];                  // We start with the first piece
            int biggestPieceSize = getRangeSize(biggestPiece); // and retrieve its size.

            // We then loop through all the pieces we have,
            // and set biggestPiece and biggestPieceSize to
            // the biggest piece and biggest piece size in
            // the array.
            for (int i = 1; i < currentNumberOfPieces; i++)
            {
                range *piece = &pieces[i];

                int pieceSize = getRangeSize(piece);

                if (biggestPieceSize < pieceSize)
                {
                    biggestPiece = piece;
                    biggestPieceSize = pieceSize;
                }
            }

            // We then call our partition method, passing in
            // the left and right index of our biggest piece
            // so that we can find the pivot index. Elements to
            // the left of the pivot index will be in one piece,
            // and elements to the right of it wil be in another.
            int pivotIndex = partition(biggestPiece->low, biggestPiece->high);

            // Now we need to create our new pieces. We can place
            // the indexes for the left piece into our biggest piece,
            // and create a new piece for the right piece, so we
            // start by grabbing our new piece.
            range *newPiece = &pieces[currentNumberOfPieces];

            // We then set its left index to the index of the
            // element right after our pivot index. We also set
            // its right index to the right index of the biggest
            // piece.
            newPiece->low = pivotIndex + 1;
            newPiece->high = biggestPiece->high;

            // We then shrink our biggest piece by setting
            // its right index to our pivot index minus one,
            // so that it compromises the values to to the left
            // of the pivot index.
            biggestPiece->high = pivotIndex - 1;

            // We then increment our current number of pieces,
            // since we have added a new piece.
            currentNumberOfPieces++;
        }

        // We then mark our end time and
        // calculate and set the time it
        // took to partition our array.
        endTime = clock();
        partitionTime = getTimeBetweenClocks(startTime, endTime);

        // Now we have all of our pieces, so we need
        // to sort them in descending order by size.
        // We do this so we can have our threads sort
        // the pieces from biggest to smallest, as
        // bigger pieces will take longer than smaller
        // ones to sort. I did this utilizing an
        // insertion sort.
        //
        // We start at the second element in the array.
        for (int i = 1; i < currentNumberOfPieces; i++)
        {
            range key = pieces[i]; // We copy the piece to our key,
            int j = i - 1;         // and look at the element to the left of the element at j.

            // While our j value is still in the bounds of the array and
            // the size of the key piece is greater than the size of the
            // piece at j,
            while (j >= 0 && getRangeSize(&key) > getRangeSize(&pieces[j]))
            {
                // we set the item at j + 1 to the item at j
                // and we decrement j.
                pieces[j + 1] = pieces[j];
                --j;
            }

            // Finally, we set the index at j + 1
            // to our key piece.
            pieces[j + 1] = key;
        }

        // Now that we have sorted our pieces from largest to smallest,
        // we can start timing how long it takes to sort the entire
        // array.
        startTime = clock();
        gettimeofday(&startTimeOfDay, NULL);

        pthread_t threads[maxThreads];
        pthread_attr_t threadAttributes[maxThreads];

        int nextPieceToSchedule = 0;

        // We start our first set of threads which will sort
        // the first few biggest pieces. The amount of pieces
        // we will sort is equivalent to the amount of max
        // threads we have.
        for (int i = 0; i < maxThreads; i++)
        {
            // We get a pointer to the piece
            // we are about to run a thread for.
            range *piece = &pieces[nextPieceToSchedule];

            // We then print the low index, high index, and range
            // size.
            printf("(%9d, %9d, %9d)\n", piece->low, piece->high, getRangeSize(piece));

            // We initialize the thread attributes
            // for the thread we are about to start.
            pthread_attr_init(&threadAttributes[i]);

            // We then create the thread, passing in our
            // the thread and thread attributes, as well as
            // our runner function which the thread will
            // run when it is started. We also pass our piece
            // pointer in so that the thread knows what segment
            // of the array to sort.
            pthread_create(&threads[i], &threadAttributes[i], runner, piece);

            // We increment our next piece to schedule
            // since we have successfully scheduled a piece.
            nextPieceToSchedule++;
        }

        // We now need to check for threads
        // that have finished sorting their
        // segment of the array so we can
        // start threads in their place to
        // sort our remaining pieces.
        int threadIndex = 0; // We start with checking the first thread.

        // While we still have pieces to schedule, we
        // check for a thread that has finished execution.
        while (nextPieceToSchedule < currentNumberOfPieces)
        {
            // We get the thread at our thread index.
            pthread_t *thread = &threads[threadIndex];

            // We check if the thread has terminated.
            // If pthread_tryjoin_np returns zero,
            // then it has succeeded in joining to
            // our thread, so we can start sorting our
            // next piece.
            if (pthread_tryjoin_np(*thread, NULL) == 0)
            {
                // We get a pointer to the piece
                // we are about to run a thread for.
                range *piece = &pieces[nextPieceToSchedule];

                // We then print the low index, high index, and range
                // size.
                printf("(%9d, %9d, %9d)\n", piece->low, piece->high, getRangeSize(piece));

                // We initialize the thread attributes
                // for the thread we are about to start.
                pthread_attr_init(&threadAttributes[threadIndex]);

                // We then create the thread, passing in our
                // the thread and thread attributes, as well as
                // our runner function which the thread will
                // run when it is started. We also pass our piece
                // pointer in so that the thread knows what segment
                // of the array to sort.
                pthread_create(&threads[threadIndex], &threadAttributes[threadIndex], runner, piece);

                // We increment our next piece to schedule
                // since we have successfully scheduled a piece.
                nextPieceToSchedule++;

                // We set our thread index back to
                // zero, since we will start checking
                // threads from the beginning again.
                threadIndex = 0;
            }
            else // If we weren't able to join the thread,
            {
                // we increment our thread index
                // so that we can check the next
                // one.
                threadIndex++;

                // If our thread index is greater than
                // or equal to our max threads, we have
                // reached the end of our threads array
                // without finding an available thread.
                // To avoid heavy polling on the CPU,
                if (threadIndex >= maxThreads)
                {
                    // we sleep for 50 ms. We also
                    // set our thread index back to 0
                    // so we can start over from the
                    // beginning after our slumber.
                    threadIndex = 0;
                    usleep(50000);
                }
            }
        }

        // At this point, we have started threads
        // to sort all of our pieces, so now we
        // just have to go through all of our
        // threads and ensure they have finished.
        for (int i = 0; i < maxThreads; i++)
        {
            // To do this, we simply call pthread_join
            // and pass in each thread. This will cause
            // our program to wait until the thread
            // terminates.
            pthread_join(threads[i], NULL);
        }
    }
    else // without multithreading,
    {
        // we mark the start time so we can
        // calculate the time it takes to sort
        // the array.
        startTime = clock();

        // We do the same with time of day
        // as well.
        gettimeofday(&startTimeOfDay, NULL);

        // Finally, we just run our hybrid sort,
        // passing in index zero as our low index
        // and arraySize - 1 as our high index.
        // This will cause hybridSort to sort
        // the entire array.
        hybridSort(0, arraySize - 1);
    }

    // At this point, we have finished sorting the array,
    // so let's get the end time so we can calculate
    // the time taken to sort the array.
    endTime = clock();
    gettimeofday(&endTimeOfDay, NULL);
    sortCpuTime = getTimeBetweenClocks(startTime, endTime);
    double sortWallClockTime = getTimeBetweenTimeOfDays(startTimeOfDay, endTimeOfDay);

    // We have also finished sorting in general,
    // so we get the end time and end time of day so
    // we can calculate the total program time.
    clock_t programEndTime = clock();
    struct timeval programEndTimeOfDay;
    gettimeofday(&programEndTimeOfDay, NULL);
    double allCpuTime = getTimeBetweenClocks(programStartTime, programEndTime);
    double allWallClockTime = getTimeBetweenTimeOfDays(programStartTimeOfDay, programEndTimeOfDay);

    // We check if the array is sorted and print a message
    // if it is not. Theoretically, this should never happen,
    // but it is good to know that it did.
    if (!isSorted())
    {
        printf("This is very bad - the array is not sorted!\n");
    }

    // Finally, we free the memory for our array.
    // We need to be responsible with our memory!
    free(array);

    // Print out our labels and final statistics.
    // For some values like seed and multithreading,
    // we check them just to determine special cases,
    // like not providing a seed or turning off
    // multithreading.
    printf("     SIZE    THRESHOLD SD PC T CREATE   INIT  SHUFFLE   PART  SrtWall Srt CPU ALLWall ALL CPU\n");
    printf("   --------- --------- -- -- - ------ ------- ------- ------- ------- ------- ------- -------\n");
    printf("F: %9d ", arraySize);
    printf("%9d ", threshold);
    if (seed == NO_SEED)
    {
        printf("00 ");
    }
    else
    {
        printf("%2d ", seed);
    }

    if (!isMultithreaded) // If we have multithreading off,
    {
        // we set the number of pieces to one
        // and our max thread count to one for
        // printing purposes, since single threading
        // only uses a single piece and thread.
        numberOfPieces = 1;
        maxThreads = 1;
    }

    printf("%2d ", numberOfPieces);
    printf("%1d ", maxThreads);
    printf("%6.3f %7.3f %7.3f ", createTime, initTime, shuffleTime);
    printf("%7.3f ", partitionTime);
    printf("%7.3f %7.3f ", sortWallClockTime, sortCpuTime);
    printf("%7.3f %7.3f\n", allWallClockTime, allCpuTime);

    return 0; // Finally, we return zero and exit.
}

void swap(int *a, int *b)
{
    // A simple method that swaps the values at two
    // int pointers. Having this simply cuts down on
    // repetition in the code above. I have benchmarked
    // both versions of the program, one with this method
    // and one without, and had practically identical
    // execution times averaged over multiple runs. I imagine
    // gcc is being pretty smart and is automatically inling
    // this function wherever it is called. Technically,
    // this could have been written as a macro to ensure
    // inlining.
    int temp = *a;
    *a = *b;
    *b = temp;
}

int partition(int p, int r)
{
    // This method partitions a segment of our
    // array with a left index of p and a right
    // index of r, sorting a single element and
    // returning a pivot index. This pivot index
    // can be used to split the segment into
    // two smaller segments.
    int i = p;        // Start at i on the left side
    int j = r + 1;    // and one past r on the right side,
    int x = array[p]; // and set x to the value on the left.
    do
    {
        // We then increase i while the value at i is less than
        // x, and we decrease j while the value at j is greater
        // than x.
        do
            i++;
        while (array[i] < x);
        do
            j--;
        while (array[j] > x);
        if (i < j) // If i < j, we swap the values at i and j,
        {
            swap(&array[i], &array[j]);
        }
        else // otherwise we break out of our do...while.
        {
            break;
        }
    } while (true);

    // Finally, we swap th values at indices p and j,
    // and return j, our pivot index.
    swap(&array[p], &array[j]);
    return j;
}

void shellSort(int low, int high)
{
    // This method performs a shell sort
    // using Hibbard's sequence on a segment
    // of our global array ranging from the low
    // index to the high index.
    int k = 1;                 // We set k equal to 1
    int size = high - low + 1; // and calculate the size of our segment.

    while (k <= size) // While k <= size,
    {
        k *= 2; // we multiply it by two,
    }
    k = (k / 2) - 1; // and finally cut it in half and subtract one.

    // We then loop while our k > 0 and ensure
    // we enter the loop at least once.
    do
    {
        for (int i = low; i < (low + size - k); i++)
        {
            // We start out our low value and loop until we have
            // have reached our low value plus the segment size minus
            // k.
            for (int j = i; j >= low; j -= k)
            {
                // We then loop starting at i again, going until we have reached
                // the low index of our segment. If the array value of j is less
                // than or equal to the value at j + k, we break. Otherwise,
                // we swap the two elements.
                if (array[j] <= array[j + k])
                {
                    break;
                }
                else
                {
                    swap(&array[j], &array[j + k]);
                }
            }
        }

        // We divide k by 2, then potentially
        // enter the loop again.
        k /= 2;
    } while (k > 0);
}

void hybridSort(int low, int high)
{
    // This method performs our hybrid sort on
    // a certain segment of the array, using quick
    // sort in most cases, but switching to shell sort
    // if the size of the segment is less than or equal
    // to the threshold.
    if (low < high) // If our low index is less than our high index,
    {
        // we get the size of the segment.
        int size = high - low + 1;

        if (size < 2) // If our segment has a size of one,
        {
            // a segment with one element is already sorted! Sweet.
            return;
        }

        if (size == 2) // If our segment has a size of two,
        {
            // we compare the two elements, and
            // swap them if the value at the low index
            // is bigger than the value at the high index.
            if (array[low] > array[high])
            {
                swap(&array[low], &array[high]);
            }
        }
        else if (size <= threshold) // If our size is <= our threshold,
        {
            // we just use shell sort to sort the segment.
            shellSort(low, high);
        }
        else // Otherwise, we use quicksort,
        {
            // determine our pivot value using partition,
            // then recursively call our hybridSort method
            // to sort the two smaller segments of this particular
            // segment.
            int q = partition(low, high);
            hybridSort(low, q - 1);
            hybridSort(q + 1, high);
        }
    }
}

void *runner(void *parameter)
{
    // This method is simply a runner function
    // that will take the range parameter and
    // perform a hybrid sort on the segment it
    // represents, passing the low index and high
    // index from the range. After it is done,
    // it exits the thread with an exit code
    // of zero.
    range *rangeParameter = (range *)parameter;
    hybridSort(rangeParameter->low, rangeParameter->high);
    pthread_exit(0);
}

bool isSorted()
{
    // This method returns whether our array is
    // sorted. Since our array will always
    // have unique elements ranging from 0
    // to arraySize - 1, we can simply check
    // if the value at each index in the array
    // equals the index itself.
    for (int i = 0; i < arraySize; i++)
    {
        if (array[i] != i)
        {
            return false;
        }
    }

    return true;
}
