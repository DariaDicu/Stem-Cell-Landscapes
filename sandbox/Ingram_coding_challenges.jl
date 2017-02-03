File containing summary of tasks and coding tidbits that I aim to work on.

1. Proof of concept for combining ODEs and simulation code.
...Take a simple ODE system, conduct simulations based on Ivan's latin
...hypercube, and track trajectories through quadrants. Sum up each quadrant
...as a trajectory passes through. Should this be just the same as the
...gradient of a probability landscape? May reveal certain points that a cell
...would always pass through, that you wouldn't realise otherwise.

Note: this is a great way to form the link between Ivan's and Daria's code.

2. Stochastic DE over time.
...Here, a trajectory will eventually jump out and travel towards another
...stationary state. I should track when this jump occurs and plot the
...trajectory it takes. Take elements of 1.

3. Explore and tweak pitchfork bifurcation.
...Take an ODE system that produces this bifurcation and change the parameter
...value such that the landscape changes. See how trajectories behave with
...this value being tweaked? Simulate in between stationary states? Take
...elements of 1.

Note: as the parameter changes, the relative intensities AND positions of the
stationary states could change. Not much work has gone on this. With simple
systems, it's expected that position won't change much, so can still visualise
paths between stationary states.

4. Path of least resistance.
...Perhaps the simpest way of visualising the reprogramming routes. Get matrix
values for heights of surrounding contours, and simulate the different paths
that could be taken between stationary states. Look up this algorithm online.
