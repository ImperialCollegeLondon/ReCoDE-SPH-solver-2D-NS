# STL usage

Throughout the previous version of the code base there has been extensive use of C-style arrays and raw pointers (v0-v2 of the project's repository). This could be harmless in the case that the code base stays the same and the developers can be sure about their memory management and ownership semantics. But, as the code base expands, it will be more difficult to track every object, moving from function to function.

Modern C++ uses the STL which is a library that comprises a lot of algorithms as well as standard containers, most of which have been optimised in the best possible way (though some of [them](https://stackoverflow.com/a/17797560) still cause problems).

In our attempt to make the code base comply with Modern C++ practises, we have changed all of the dynamically allocated arrays with `std::vector` and all of the raw pointers with `std::unique_ptr`.

## std::vector 

The vector container behind the scenes is simply an array, minus the painful memory management that would be needed in case we have used raw arrays. It can create an overhead since the developer can expand the container by using, i.e. `push_back()` method, but there is a catch. The `std::vector` allocates a specified space in memory at initialization. If the developer "pushes" more objects in it than memory exists, it would need to allocate another space, bigger this time, and copy all (!!) the elements from one container to the other. As the reader can imagine, this can create a serious overhead in really big vectors. But, the authors of STL have thought about this case and they have included a pretty handy method called `reserve`. With `reserve`, the developer can allocate exactly how much memory the object will need, given that the size is a known number. This could relax the constraints about performance hits. 

Another useful thing about vectors is that they use contiguous memory, meaning each element <i>should</i> be "side to side" inside memory. That means that a simple access should give an average of $O(1)$ in terms of time. Since the code base includes a lot of vector accesses, we assesed that this container is the best option for our cases.

## std::unique_ptr

Smart pointers was an addition to the STL in order to avoid manual memory allocation for an object and then (and most crucially) manual deallocation of its resources. Smart pointers are part of an idiom called [RAII](https://en.cppreference.com/w/cpp/language/raii) which binds the lifetime of a resource that must be acquired before use (like memory on the heap) to the lifetime of an object. What this means is that the smart pointer is only valid once it's in scope. Once it gets out of scope, all of its resources are properly free'd. This way the developer can be sure about proper memory management. There are several smart pointers (namely `std::shared_ptr`, `std::unique_ptr` and `std::weak_ptr`) but we are going to talk about `std::unique_ptr` since it's the one used in the code base and it's by far the most used one.

We used `std::unique_ptr` wherever there was a need for a pointer in the code base. There's nothing wrong with using raw pointer, but raw pointers shouldn't handle ownership, meaning a raw pointer is good enough only when it doesn't own anything (i.e passing it to a function that does some calculations).

One thing that stands out for `std::unique_ptr` is that it cannot be copied. That means that the ownership of the managed resource needs to be transferred to a new object, if the developer wants to move the pointer. The original object will no longer own the resource and the original `std::unique_ptr` will contain a `nullptr`. To create a `std::unique_ptr` it is good practise to use 
`std::make_unique` as can be seen in the snippet below:

```C++
auto fluidPtr = std::make_unique<Fluid>(nbParticles);
```

After this the user can use `fluidPtr` as he would use a raw pointer, with the difference that he doesn't (and he shouldn't) delete the pointer at the end of the scope.