# Exercise 1

Group members: **Patrick Lanzinger, Raphael Gruber**

## Study how submit jobs in SGE, how to check their state and how to cancel them.

To submit the following command is used with **qsub**:

```
qsub [ -q std.q ] [options] job_script [ job_script_arguments ...]
```

To check the state of a running job you can use **qstat**. It has some options like -u, -j, -f to print out the jobs of a user, the information about a specific job or to print out all queues and all jobs.

The only way to stop a job that already runs is to use **qdel**.

## In your opinion, what are the 5 most important parameters available when submitting a job and why? What are possible settings of these parameters, and what effect do they have?

### Queuename

With -q you can set the name of the queue in which the program will run. There exist three possible parameters std.q, short.q, bigmem.q. These three queues have all different properties. The first one std.q is the standard queue it runs for a maximum of 240 hours. The short.q queue is for small jobs, you only get a limited number of CPU slots and the maximum runtime is only 10 hours. The last one bigmem.q is meant for jobs with high memory requirements.

### Email address

With the option -M you can set an email address. This will then be used to notify you about the state of a job.

### Notifications

With the option -m you have a number of options that set for what event a notification will be sent (b = begin, e = end, a = abort, s = suspend, n = no mail). 

### CWD

This option (-cwd) is very important if you want to control the directory from which the job is run. When omitted the **$home** directory is used and this is usually not a good idea.

### Syntax checking

There is also an option to check if the syntax of the job is correct (-w v).


## How do you run your program in parallel? What environment setup is required? 

If you don't add the option -pe, SGE just assumes the job is sequential. The -pe option has two parameters the first one is the environment and the second one is the number of slots. The cluster has three environments openmpi-Xperhost, openmpi-fillup, openmp. The most important thing to note here is that this only reserves CPU-cores but you yourself are responsible that the job actually uses them.

To run the test script: 

```
qsub job.script
```
