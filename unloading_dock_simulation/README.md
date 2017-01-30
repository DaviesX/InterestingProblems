#Unloading dock simulation

You will write a simulation of a train unloading dock. Trains arrive at the station as a Poisson process on average once every 10 hours. Each train takes between 3.5 and 4.5 hours, uniformly at random, to unload. If the loading dock is busy, then trains wait in a first-come, first-served queue outside the loading dock for the currently unloading train to finish. Negligible time passes between the departure of one train, and the entry of the next train (if any) into the loading dock---unless the entering train has no crew (see below).
        
There are a number of complications. Each train has a crew that, by union regulations, cannot work more than 12 hours at a time. When a train arrives at the station, the crew’s remaining work time is uniformly distributed at random between 6 and 11 hours. When a crew abandons their train at the end of their shift, they are said to have “hogged out”. A train whose crew has hogged out cannot be moved, and so if a hogged-out train is at the front of the queue and the train in front finishes unloading, it cannot be moved into the loading dock until a replacement crew arrives (crews from other trains cannot be used). Furthermore, a train that is already in the loading dock cannot be unloaded in the absence of its crew, so once the crew hogs out, unloading must stop temporarily until the next crew arrives for that train. This means that the unloading dock can be idle even if there is a train in it, and even if it is empty and a (hogged-out) train is at the front of the queue. Once a train’s crew has hogged out, the arrival of a replacement crew takes between 2.5 and 3.5 hours, uniformly at random. However, due to union regulations, the new crew’s 12-hour clock starts ticking as soon as they are called in for replacement (i.e., at the instant the previous crew hogged out); i.e., their 2.5-3.5 hour travel time counts as part of their 12-hour shift. You will simulate for 7200 hours, and output the following statistics at the end:
1. Total number of trains served.
2. Average and maximum over trains of the time-in-system.
3. The percentage of time the loading dock spent busy, idle, and hogged-out (does this add to
100%? Why or why not?)
4. Time average and maximum number of trains in the queue.
5. A histogram of the number of trains that hogged out 0, 1, 2, etc times.

###Bonus work: 
 If all other parameters are fixed, estimate the long-term maximum train arrival rate that this system can accommodate before overloading. Perform the estimate in two ways: first, a completely paper-and-pencil analysis; second, using your simulation. Explain your rationale in both cases; especially, in the case of the simulation, explain what criteria you used to decide if the system was “overloaded”.
