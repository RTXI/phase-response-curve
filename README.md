###Phase Response Curve

**Requirements:** None  
**Limitations:** None  

![PRC GUI](phase-response-curve.png)

<!--start-->
This module applies an alpha-shaped conductance to the cell at a fixed delay after 10 interspike intervals (ISI). It computes an intrinsic period P0 by averaging the most recent 5 of 10 ISIs. Thus, the value of P0 can change over time. The period after the perturbed period is designated as P1 and the period following P1 is designated as P2. The first and second order PRCs are computed as PRC1=(P1-P0)/P0 and PRC2=(P2-P0)/P0. By this convention, a delay in the spike yield a positive value and an advance yields a negative value. The phase of the perturbation stimulus is computed based on the most recent measured P0. To help you choose the fixed delay for the perturbations, use the “Measure intrinsic P0″ at the top of the module to begin computing a running average of measured ISIs. After you start running the module, the ISI textbox will simply display the current measured ISI. To run a PRC protocol by stepping through the specified fixed delays, click the “Measure PRC” button. Stop either mode of data collection by clicking the “Pause” button. The source code is heavily commented and the program logic is based on indexing the number of spikes that have occurred in each stimulus cycle. It is easy to change how the intrinsic period is computed and whether it is held constant during the entire experiment or recomputed. Similarly, you can easily change the sign convention for computing the phase shifts and extend the code to compute other phase response measures.
<!--end-->

This screenshot was made using the Connor Stevens model to generate spikes. A complete demo is available for re-creating this screenshot. 

![PRC and Scope](PRC-scope.png)

This module has a second plotting mode that creates a ts-tr plot of the stimulus times and response times for the perturbed period.

![PRC (ts-tr)](PRC-tstr.png)

#### Input
1. input(0) - Vm : Membrane voltage of real cell (V)

#### Output
1. output(0) - Command : Command perturbation current
2. output(1) - Cell spike state : Real cell spike state

#### Parameters
1. Min Delay (ms) - Minimum stimulus delay (ms)
2. Max Delay (ms) - Maximum stimlus delay (ms)
3. Step Size (ms) - Step size between minimum and maximum delay (ms)
4. Gmax (nS) - Maximum synaptic conductance for stimulus
5. Time Constant tau (ms) - Time constant for alpha-shaped conductance
6. Esyn (mV) - Reversal potential for stimulus
7. Repeat - Number of times to run cycle
8. Threshold (mV) - Threshold (mV) at which to detect a spike
9. Min Interval (s) - Minimum interval(s) that must pass between spikes

#### States
1. ISI (ms) - ISI (ms)
2. Period Number - Period Number
3. Time (s) - Time
