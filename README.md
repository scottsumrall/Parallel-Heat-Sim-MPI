# Parallelized Heat Transfer Simulation

This program produces a simulation of heat transfer within an arbitrary grid based on a user-specified transfer rate. The transfer is represented graphically using the GFX graphics library. Performance has been optimized to perform at scale with use of the Message Passing Interface (MPI) to distribute the calculation load across several different processors.

## Getting Started

These instructions will give you a copy of the project up and running on
your local machine for development and testing purposes.

### Prerequisites

If utilizing remote linux server: 
- [Xming X Server (Windows)](https://sourceforge.net/projects/xming/)
- [XQuartz (macOS)](https://www.xquartz.org)

Set up MPI on the server: 
- [OpenMPI](https://rantahar.github.io/introduction-to-mpi/setup.html)

### Compilation

A step by step series of examples that tell you how to get a development
environment running

Deploy heat_parallel.c to server and compile

    mpicc heat_parallel.c -o heat_parallel

Run using MPI parallelized over 'X' threads

    mpirun -np X heat_parallel
    
Distributed hosts can be specified by adding a '--hostfile' tag

    mpirun -np X --hostfile <filename> heat_parallel

If on Windows OS, use the following command in Powershell for graphical output

    $env:Display="localhost:0"

## License

This project is licensed under the [CC0 1.0 Universal](LICENSE.md)
Creative Commons License - see the [LICENSE.md](LICENSE.md) file for
details
