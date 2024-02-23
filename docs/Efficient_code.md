# Good practices for computational efficiency

The complexity of the SPH algorithm requires the code implementation to display good computational performance in order for the program to be able to handle cases involving large numbers of particles. In this project we attempted to apply several standard good practices in order to avoid unecessary computations or accesses to memory locations, and therefore to achieve a reasonable and useful execution time. It is good practice to apply such optimisations and techniques without compromising the memory usage, readability and future maintenance of the code, and always to balance the trade-off between these three.

## Optimisation flags

Optimisation flags are compiler directives (introduced in the Makefile) that instruct the compiler to apply various optimisations when translating high-level source code (in this case C++ code) into machine code (the binary executable). These flags are essential tools for improving the performance and efficiency of a program. Common optimisation flags, such as -O1, -O2, or -O3 in GCC and Clang, enable different levels of optimisation, ranging from basic to more aggressive transformations. These optimisations include inlining functions, constant folding, dead code elimination, loop unrolling, and many others. In practice, judicious use of optimisation flags, coupled with profiling tools, allow for fine-tuning of the code for optimal performance while maintaining a balance between speed and other considerations.

## Properly arranged `for` loops

When using dynamically allocated containers such as arrays or std::vectors, these structures allocate a consecutive block of memory on the heap. In the case of std::vector, although the elements are physically contiguous, there might be cases in which the memory that was allocated is not large enough and the vector needs to be resized. This resize will copy all the elements from the old location in the memory to the new larger space, which might lead to a performance penalty. If the size of the container is known at the size of initialisation, the program can reserve the exact space that it needs, preventing any unnecessary copying of elements.

When accessing elements within these containers, the program checks if the data are present in the cache. If not, it needs to fetch the data from memory. Accessing data from the cache is significantly faster than fetching it from memory, providing a performance improvement.

It is worth noting that the hardware cache management involves retrieving data in larger chunks known as cache lines, optimising for spatial locality. This means that when data are retrieved from memory, neighbouring data that might be used soon are also fetched to the cache, contributing to improved overall performance.

To leverage this functionality, the nested `for` loops can be organized to facilitate faster access to the next entry of the array. When the outer loop iterates over the `i` elements and the nested loop iterates over the `j` elements (i * elements + j), this arrangement increases the chances of accessing elements that are already in the cache line. This optimisation makes the program cache-friendly, resulting in greater performance. In simpler terms, by structuring the loops to access nearby elements in memory, the program becomes more efficient and runs faster.

```cpp
  void Fluid::calculateDensity() {
  double phi, normalisedDistance, normalisedDistanceSqr;

  // find φ
  for (size_t i = 0; i < nbParticles; i++) {
    density[i] = 0.0;

    for (size_t j = 0; j < nbParticles; j++) {
      normalisedDistance = distanceQ[i * nbParticles + j] =
          std::abs(distance[i * nbParticles + j] * hInverse);

      normalisedDistanceSqr = (1.0 - normalisedDistance * normalisedDistance);
      if (normalisedDistance < 1.0) {
        phi = fourPih2 * normalisedDistanceSqr * normalisedDistanceSqr *
              normalisedDistanceSqr;

      } else {
        phi = 0.0;
      }

      density[i] += mass * phi;
    }
  }
}
```

## Pre-calculated values

The algorithm provides the equations which describe the calculations to be performed during one iteration. However, when we transfer the algorithm in a code, some of the equations have constant parts (either throughout the whole simulation or throughout a timestep), and these do not need to be re-calculated every time the code reaches the corresponding algorithmic step. It is the developer's responsibility to identify these calculations and pre-calculate all such values. A demonstration of this concept in our SPH code is the calculation of the pressure force. A somewhat naive implementation of this step (on branch `v2` in the repository of the project) would be:

```cpp
  double SphSolver::calculatePressureForce(Fluid &data,
                                         std::function<double(int)> getPosition,
                                         int particleIndex) {
  double sum = 0.0;  // Initializing the summation
  double radiusOfInfluence = data.getRadInfl();
  double thirtyPih3 =
      (-30.0 / (M_PI * radiusOfInfluence * radiusOfInfluence *
                radiusOfInfluence));  // Precalculated value used to avoid
                                      // multiple divisions and multiplications

  for (int j = 0; j < numberOfParticles; j++) {
    if (particleIndex != j) {
      if (data.getDistanceQ(particleIndex * numberOfParticles + j) < 1.0) {
        sum +=
            (data.getMass() / data.getDensity(j)) *
            ((data.getPressure(particleIndex) + data.getPressure(j)) / 2.0) *
            (thirtyPih3 * (getPosition(particleIndex) - getPosition(j))) *
            (((1.0 - data.getDistanceQ(particleIndex * numberOfParticles + j)) *
              (1.0 -
               data.getDistanceQ(particleIndex * numberOfParticles + j))) /
             data.getDistanceQ(particleIndex * numberOfParticles + j));
      }
    }
  }
  return -sum;
}
```

In the above piece of code there is an attempt to pre-calculate some values by defining the `double thirtyPih3`. However, not only are there a lot more things to be improved, but this specific optimisation could be achieved in a better way. More specifically, this value is constant throughout the entire simulation, since none of its components change at runtime. It is therefore better practice to define it as a data member of the `SphSolver` class and calculate it only once during the initialization process.

```cpp
  // Assign values to pre-calculated values
  inline void setPrecalculatedValues(double radiusOfInfluence) {
    thirtyPih3 = -30.0 / (M_PI * std::pow(radiusOfInfluence, 3.0));

    fourtyPih4 = 40.0 / (M_PI * std::pow(radiusOfInfluence, 4.0));
  }
```

In addition to that, we can identify several values which can be stored as local variables within the scope of the function and be used in the `j` loop, instead of invoking the getter functions of the `Fluid` class for every `j`. As we saw earlier, there is a possible overhead associated with calling a function and accessing elements in an array or vector, and this practice helps us alleviate it. Finally, since the iterations are performed on the `j` index, there is a constant value in the index of `data.getDistanceQ(particleIndex * numberOfParticles + j)` which is the `particleIndex * numberOfParticles` and this can also be pre-calculated and stored as a function variable. Finally, the above function has taken the following form:

```cpp
  double SphSolver::calculatePressureForce(Fluid &data,
                                         std::function<double(int)> getPosition,
                                         int particleIndex) {
  double sum = 0.0;  // Initializing the summation
  double normalisedDistance;
  double position = getPosition(particleIndex);
  double pressure = data.getPressure(particleIndex);
  double mass = data.getMass();
  int index = particleIndex * numberOfParticles;

  for (int j = 0; j < numberOfParticles; j++) {
    if (particleIndex != j) {
      normalisedDistance =
          data.getDistanceQ(index + j);  // Store this variable to avoid
                                         // multiple calls of the get function
      if (normalisedDistance < 1.0) {
        sum += (mass / data.getDensity(j)) *
               ((pressure + data.getPressure(j)) / 2.0) *
               (thirtyPih3 * (position - getPosition(j))) *
               (((1.0 - normalisedDistance) * (1.0 - normalisedDistance)) /
                normalisedDistance);
      }
    }
  }
  return -sum;
}
```

One may notice that the repeated operation `1.0 - normalisedDistance` was not replaced. However, this is part of the trade-off, balancing readability with performance. This would require creating a new variable with a name such as `oneMinusNormDistance` and setting a value for this variable at every iteration, in order to avoid a single subtraction. The computational saving is likely negligible (if not none) and the code readability would have been compromised. By looking at the updated version of the pressure force calculation function, although this is not always the case, one may notice that apart from accelerating the code, its readability in this case has improved drastically.

## Inlined functions

Inlining functions in C++ is a compiler optimisation technique where the compiler replaces the function call with the actual body of the function at the call site. This eliminates the overhead of the function call itself, leading to potentially faster execution. Inlining is particularly useful for small, frequently called functions, as it reduces the function call overhead and can result in more efficient code. Additionally, inlining can facilitate further optimisations by providing the compiler with opportunities for constant folding, dead code elimination, and other performance enhancements. Prudent use of inline functions, particularly for small, performance-critical code sections, can contribute to more efficient and streamlined C++ programs. In our code, we made use of the `inline` keyword for all the getter and setter functions which are small and qualify for this usage.

```cpp
  // Function to get the density felt by a single particle
  inline double getDensity(int index) { return density[index]; }
```
