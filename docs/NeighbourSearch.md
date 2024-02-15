# Neighbour Searching Algorithm

In SPH, each particle is only influenced by the particles that lie within a maximum distance equal to the radius of influence (`h`) parameter. These particles are usually called the particle's neighbours. 

## Big O Notation

Big O notation is a mathematical notation used to describe the asymptotic behavior of a function or algorithm. It's commonly used in computer science to analyze the performance of algorithms in terms of their time complexity.

Big O notation provides a way to classify algorithms based on their efficiency and scalability. It helps us understand how the time or space requirements of an algorithm grow as the size of the input increases. In Big O notation, we express the upper bound or worst-case scenario of an algorithm's time complexity in terms of a mathematical function.

For example, an algorithm with a time complexity of O(n) means that its execution time grows linearly with the size of the input. An algorithm with O(n^2) time complexity means that its execution time grows quadratically with the size of the input, and so on.

## Brute Force Approach

A brute force algorithm is a straightforward approach to solving a problem that exhaustively tries all possible solutions. For certain problems, this can lead to exponential time complexity, making brute force algorithms impractical for large inputs.

In out initial implementation, we followed a "brute force" approach, where, for each particle, we would iterate over all particles, calculate their distance, compare it to the radius of influence, and finally, if the particle was a neighbour, we we would proceed to the appropriate calculations. This approach resulted in a Big O notation complexity of O(n^2), denoting an inefficient algorithm that should be replaced.

## Cells Neighbour Searching

A cells neighbor searching algorithm is an optimization technique used to efficiently find neighboring elements within a spatial grid or partitioned space. Instead of comparing each element with every other element, this approach divides the space into smaller cells or buckets and only considers elements within the same or neighboring cells.

For example, in the context of particle simulation or computational fluid dynamics, a cells neighbor searching algorithm can significantly reduce the number of pairwise comparisons needed to calculate interactions between particles. By organizing particles into spatial cells and only considering interactions between particles within the same or adjacent cells, the algorithm achieves a time complexity of O(n), where n is the number of particles, making it much more scalable and efficient than brute force approaches for large datasets.

In the current implementation, we follow a neighbour searching algorithm that is based on partitioning the domain into square grid cells. These cells are of side equal to the radius of influence. The benefit of this design decision is that we know in advance that a particle's neighbours may only lie in the same cell as the particle, or in one of the neighbouring cells. Therefore, we do not need to iterate over all partices and check their distances before every calculation.

```
void SphSolver::createGrid(Fluid &data) {
  numberOfParticles = data.getNumberOfParticles();
  neighbourParticles.resize(numberOfParticles);

  double radiusOfInfluence = data.getRadInfl();
  int cellsRows =
      static_cast<int>(std::ceil((topWall - bottomWall) / radiusOfInfluence));
  int cellsCols =
      static_cast<int>(std::ceil((rightWall - leftWall) / radiusOfInfluence));
  numberOfCells = cellsRows * cellsCols;
  cells.resize(numberOfCells);
  neighbourCells.resize(numberOfCells);

  assignNeighbourCells(cellsRows, cellsCols);
}
```

The first step of the algorithm is to turn the domain to a grid, by splitting it into cells. This logic is implemented in the `SphSolver::createGrid` function. It should be noted that all neighbour searching functions have been implemented as member functions of the `SphSolver` class, since the algorithm is a part of the SPH system. In this function, we calculate the number cell rows and columns that form the grid, by deviding the corresponding domain dimension (length and width accordingly) by the cell sidelength (radius of influence). Since the result of this division is not guaranteed to be an integer, we use `ceil` to get the next integer number, and `static_cast<int>` to change the `double` to and `int`, which is finally assigned as the cell rows/columns. The reason we use `ceil` to round up to the next integer, instead of just rounding down by simply casting the `double` to an `int`, is to be sure that the whole domain is considered. This way, even though parts of some cells may lie outside the domain, it is guaranteed that the algorithm will not miss any particles. Finally, we resize the two containers that will be used to store each cell's particles (`cells`) and neighbouring cells (`neighbourCells`). Each of these containers is a vector of vectors. The outer vector is used to identify cells. In the case of the `cells` container, each inner vector will hold the indices of the particles that lie within the cell, and in the case of the `neighbourCells`, the inner vector will hold the indices of the current cell's neighbouring cells. This way, the algorithm will be able to search for neighbour particles both in a cell, as well as in its neighbouring cells.

```
void SphSolver::assignNeighbourCells(int cellsRows, int cellsCols) {
  // Flags to check if the cell is on the edge or in the middle
  bool top = false;
  bool left = false;
  bool right = false;
  bool bottom = false;
  bool middle = false;

  for (int i = 0; i < numberOfCells; i++) {
    // Cell has a bottom neighbour
    if (i >= cellsCols) {
      neighbourCells[i].push_back(i - cellsCols);
      bottom = true;
    }
    ...
    // If the cell is not on the edge, add the diagonal neighbours
    if (bottom && top && left && right) {
      neighbourCells[i].push_back(i - 1 - cellsCols);
      neighbourCells[i].push_back(i - 1 + cellsCols);
      neighbourCells[i].push_back(i + 1 - cellsCols);
      neighbourCells[i].push_back(i + 1 + cellsCols);
      middle = true;
    }
    // If the cell is on the edge, add only specific neighbours
    if (!middle) {
      // Add bottom-left diagonal neighbour
      if (bottom && left) {
        neighbourCells[i].push_back(i - 1 - cellsCols);
      }
    ...
  }
}
```
The next step of the algorithm is to assign to each cell all its neighbouring cells. This logic is implemented in the `SphSolver::assignNeighbourCells` function. Depending to its position in the grid (middle, edge, corner) a cell could have 8, 5, or 3 neighbouring cells accordingly. In this function, while we iterate over all cells, we use a number of flags (`top`, `left`, `right`, `bottom`, `middle`) in combination to a number of rules, to find eacg cell's neighbouring cells. For instance, all cell indices `>=` than the number of cell columns, denote that the cell iterator has passed the first (lowest) row of the grid, a fact that means that the current cell has a "bottom neighbour", and therefore, in the `neighbourCells` container, we add the index of the "bottom cell" (`push_back(i - cellsCols)`) to the current cell's inner vector (`neighbourCells[i]`), and set `bottom` to `true`. Following similar rules, we identify each cells position and assign the right neigbhouring cells to it, including the diagonal neighours. Finally, we reset the flags and move to the next cell.

Up to this point, the algorithm's steps included logic that does not change during the time integration procedure, and thefore, the above mentioned functions are only used during initialisation (called at the end of `initialise` in `SPH-main.cpp`). However, the next steps are steps that need to be followed in each iteration, as the particles' positions change.

```
void SphSolver::placeParticlesInCells(Fluid &data) {
  for (int i = 0; i < numberOfCells; i++) {
    cells[i].clear();
  }

  double radiusOfInfluence = data.getRadInfl();
  int cellsCols =
      static_cast<int>(std::ceil((rightWall - leftWall) / radiusOfInfluence));
  for (int i = 0; i < numberOfParticles; i++) {
    double positionX = data.getPositionX(i);
    double positionY = data.getPositionY(i);
    int j = static_cast<int>(positionX / radiusOfInfluence) +
            static_cast<int>(positionY / radiusOfInfluence) * cellsCols;
    cells[j].push_back(i);
  }
}
```

The first thing that the algorithm should do in each iteration is to place the existing particles to their corresponding cells, according to their position in the domain. This functionality is included in the `SphSolver::placeParticlesInCells` function. First, the `cells` vector is cleared to remove any particles placement from a previous iteration. Then, the function iterates over all particles, gets their coordinates, and based on them, decided on the index of the cell that each particle should be placed at. Finally, the index of the particle is added to the appropriate inner vector of the `cells` container (`cells[j].push_back[i]`).

```
void SphSolver::neighbourParticlesSearch(Fluid &data) {
  for (int i = 0; i < numberOfParticles; i++) {
    neighbourParticles[i].clear();
  }

  placeParticlesInCells(data);

  // For each cell, for each particle in the cell, find neighbour particles in
  // the cell
  for (int i = 0; i < numberOfCells; i++) {
    for (int j = 0; j < cells[i].size(); j++) {
      for (int k = 0; k < cells[i].size(); k++) {
        if (cells[i][j] != cells[i][k]) {
          double distance = sqrt(pow(data.getPositionX(cells[i][j]) -
                                         data.getPositionX(cells[i][k]),
                                     2) +
                                 pow(data.getPositionY(cells[i][j]) -
                                         data.getPositionY(cells[i][k]),
                                     2));

          if (distance <= data.getRadInfl()) {
            neighbourParticles[cells[i][j]].push_back({cells[i][k], distance});
          }
        }
      }
    }
  }
  ...
}
```

After all particles have been placed in cells, the actual neighbours searching can occur. We know that a particle's neighbours could only lie in the same cell or in a neighbouring cell. In the `SphSolver::neighbourParticlesSearch` function, for each cell, for each particle in the cell, we iterate over all other particles in the same cell, calculate the pair's distance, and if it is lower or equal to the radius of influence, we count the second particle as a neighbour of the first particle. In order to store each particle's neighbours, as well as their distance to the particle, we utilise a container of type `vector<vector<pair<int, double>>>`. In this container, the outer vector accounts for each particle and its size should be equal to the `numberOfParticles`. The inner vector is used to push back "neighbour-distance pairs". For each particle's neighbour, we store a pair of the neighbour's index (`int`) and the distance between the particle and its neighbour (`double`). This container provides a way to store and access each particle's neighbours and their corresponding distances during our calculations. Finally, the whole process is repeated for each neighbouring cell of the current cell, so that we are sure that we have identified all neighbours for each particle.

```
double SphSolver::calcViscousForce(Fluid &data,
                                   std::function<double(int)> getVelocity,
                                   int particleIndex) {
  double sum = 0.0;  // Initializing the summation
  double velocity = getVelocity(particleIndex);
  double mass = data.getMass();
  double radiusOfInfluence = data.getRadInfl();

  for (int j = 0; j < neighbourParticles[particleIndex].size(); j++) {
    if (particleIndex != neighbourParticles[particleIndex][j].first) {
      sum +=
          (mass / data.getDensity(neighbourParticles[particleIndex][j].first)) *
          (velocity - getVelocity(neighbourParticles[particleIndex][j].first)) *
          (fourtyPih4 * (1.0 - neighbourParticles[particleIndex][j].second /
                                   radiusOfInfluence));
    }
  }

  return -data.getViscosity() * sum;
}
```

The neighbour searching algorithm provides a significantly more efficient way to identify the particles that affect each particle during the SPH calculations. Even though the code may seem more complex, with up to four nested `for` loops, most of these loops iterate over a low number of neighbour particles. Additionally, the initial `numberOfParticles` range loop that was used during each iteration's calculations, as well as the distance checks, have now been replaced by a simple and faster iteration over each particle's neighbours (`for (int j = 0; j < neighbourParticles[particleIndex].size(); j++)`), resulting in significantly lower execution times compared to the brute force approach, especially as the number of particles grows.