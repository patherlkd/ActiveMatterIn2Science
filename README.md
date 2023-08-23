# ActiveMatterIn2Science

Click here for direct access: https://colab.research.google.com/github/patherlkd/ActiveMatterIn2Science/blob/main/VicsekModelCode-binder.ipynb

###### Lead 
Dr. Luke K. Davis (@patherlkd) www.drlukekdavis.com

![Alt text](./sample.png?raw=true "What the simulation looks like")

## How to use this repository

### Running the interactive simulation

1. Click the google collab button at the top of this README.

2. Click "Runtime" then "Run all", which will run the whole notebook ready for you to start playing with.

3. Scroll down to the last cell and run it as per instructions. Click the ``Run Interact`` button to run a simulation.

Note: The simulations are not instantaneous, your are running full-fledged active matter simulation(!), though should not take too long (~1-3 minutes for moderate parameters)

4. Play with the sliders:
  
  
![Alt text](./sliders.png?raw=true "What the simulation looks like")

### Guiding questions/pointers to help you explore the active matter simulations

The main thing is for you to explore the system, code, and setup to get a feel of running computational simulations. So, follow you curiosity. If you want some nudges of what to look at then see below.

This questions/pointers typically get more difficult as you progress through them.

+ For each slider, investigate what they do [check at least three values] and how they affect the runtime of the simulation.
+ How can I get all the arrows aligned? There might be more than one way to do this, so write them down!
+ How can I get all the arrows to be disordered?
+ Is the transition between disordered and all-aligned gradual or sharp?
+ How does the number of particles affect what you see?
+ How does the initial condition of the particles change things?
+ What could be the problem of only looking at short simulation "Iterations"?
+ From the code cells, extract and explain the lines of code implementing periodic boundary conditions? Research why are they used?
+ From the code cells, extract and explain the lines of code implementing the equations of motion. Write down the equations and show that discretising them allows for their numerical solution.
+ How could you edit the code to implement (isotropic) repulsion between the particles?
+ Edit the code so that the self-propulsion can be changed using a slider.
+ Edit the code to compute (and then plot) the average alignment of particles as a function of time.
+ Edit the code to include obstacles for the particles.
+ Edit the code to incorporate non-reciprocal interactions.
