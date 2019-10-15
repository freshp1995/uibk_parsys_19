# Exercise 2

Group members: **Patrick Lanzinger, Raphael Gruber**

## Download and build the OSU Micro-Benchmarks 


## After building, submit SGE jobs that run the `osu_latency` and `osu_bw` executables.


## Create a table and figures that illustrate the measured data and study them. What effects can you observe?

![Image](images/bandwith1_default.png) 
![Image](images/latency1_default.png) 

The bandwith is linear and the latency states same until 8 size then growth linear as well


## Modify your experiment such that the 2 MPI ranks are placed on

### different cores of the same socket

Bandwith
![Image](images/bandwith1_diff_core.png) 
Lantency
![Image](images/latency1_diff_core.png) 

### different sockets of the same node

![Image](images/bandwith1_diff_socket.png) 
![Image](images/latency1_diff_socket.png) 

### different nodes

![Image](images/bandwith1_diff_nodes.png) 
![Image](images/latency1_diff_nodes.png) 

## Ammend your table and figures to include these additional measurements. What effects can you observe? How can you verify rank placement without looking at performance?

### different cores of the same socket

![Image](images/bandwith1_diff_core.png) 
![Image](images/latency1_diff_core.png) 

![Image](images/bandwith2_diff_core.png) 
![Image](images/latency2_diff_core.png) 

![Image](images/bandwith3_diff_core.png) 
![Image](images/latency3_diff_core.png) 

### different sockets of the same node

![Image](images/bandwith1_diff_socket.png) 
![Image](images/latency1_diff_socket.png) 

![Image](images/bandwith1_diff_socket.png) 
![Image](images/latency1_diff_socket.png) 

![Image](images/bandwith1_diff_socket.png) 
![Image](images/latency1_diff_socket.png) 

### different nodes

![Image](images/bandwith1_diff_nodes.png) 
![Image](images/latency1_diff_nodes.png) 

![Image](images/bandwith1_diff_nodes.png) 
![Image](images/latency1_diff_nodes.png) 

![Image](images/bandwith1_diff_nodes.png) 
![Image](images/latency1_diff_nodes.png) 


## How stable are the measurements when running the experiments multiple times?
