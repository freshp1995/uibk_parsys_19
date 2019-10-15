# Exercise 2

Group members: **Patrick Lanzinger, Raphael Gruber**

## Download and build the OSU Micro-Benchmarks 


## After building, submit SGE jobs that run the `osu_latency` and `osu_bw` executables.


## Create a table and figures that illustrate the measured data and study them. What effects can you observe?

Bandwith

![Image](images/bandwith1_default.png) 

Lantency

![Image](images/latency1_default.png) 

The bandwith is linear and the latency states same until 8 size then growth linear as well


## Modify your experiment such that the 2 MPI ranks are placed on

### different cores of the same socket

Bandwith

![Image](images/bandwith1_diff_core.png) 

Lantency

![Image](images/latency1_diff_core.png) 

### different sockets of the same node

Bandwith

![Image](images/bandwith1_diff_socket.png) 

Lantency

![Image](images/latency1_diff_socket.png) 

### different nodes

Bandwith

![Image](images/bandwith1_diff_nodes.png) 

Lantency

![Image](images/latency1_diff_nodes.png) 

## Ammend your table and figures to include these additional measurements. What effects can you observe? How can you verify rank placement without looking at performance?

The latency is less in by using different sockets of the same node compare to the other two. The bandwith is higher by using different cores of the same socket. In general all growth lineraly.

## How stable are the measurements when running the experiments multiple times?


### different cores of the same socket

Bandwith 1

![Image](images/bandwith1_diff_core.png) 

Latency 1

![Image](images/latency1_diff_core.png) 

Bandwith 2

![Image](images/bandwith2_diff_core.png) 

Latency 2

![Image](images/latency2_diff_core.png) 

Bandwith 3

![Image](images/bandwith3_diff_core.png) 

Latency 3

![Image](images/latency3_diff_core.png) 

### different sockets of the same node

Bandwith 1

![Image](images/bandwith1_diff_socket.png) 

Lantency 1

![Image](images/latency1_diff_socket.png) 

Bandwith 2

![Image](images/bandwith2_diff_socket.png) 

Lantency 2

![Image](images/latency2_diff_socket.png) 

Bandwith 3

![Image](images/bandwith3_diff_socket.png) 

Lantency 3

![Image](images/latency3_diff_socket.png) 

### different nodes

Bandwith 1

![Image](images/bandwith1_diff_nodes.png) 

Lantency 1

![Image](images/latency1_diff_nodes.png) 

Bandwith 2

![Image](images/bandwith2_diff_nodes.png) 

Lantency 2

![Image](images/latency2_diff_nodes.png) 

Bandwith 3

![Image](images/bandwith3_diff_nodes.png) 

Lantency 3

![Image](images/latency3_diff_nodes.png) 

We get similar results. The latency sometimes varies slightly in low sizes, but after that it always grows linearly
