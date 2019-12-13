# HPC239
1. Please first put data file into the repo folder. Our data can be downloaded from https://drive.google.com/drive/folders/1iHUwqaflpWB6_3leGFGzuSdD3Vljv063?usp=sharing
2. To run the basic Kernel K-Means, 

    'g++ -std=c++11 -c kmeans.cpp'
  
    'g++ -o kmeans kmeans.o'
  
    './kmeans'
    
   The number of data points can be edited in Line 17.
3. To run OpenMP parallel version,

    'g++ -std=c++11 kmeans_openmp.cpp -fopenmp -o kmeans_openmp'
    
    './kmeans_openmp'
    
   The number of threads and data points can be edited in Line 17 and 18.
4. To run MPI parallel version,

    'mpic++ -std=c++11 mpi_kernel_kmeans.cpp -o mpi_kernel_kmeans'
    
    'mpirun -np 6 mpi_kernel_kmeans' where 6 is the number of cores used for parallel computing
    
    The number of cores and data points can be edited in Line 378 and 379.
