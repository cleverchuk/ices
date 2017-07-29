## Phase Field Modelling nano-particle sintering
### Introduction
My project involved simulating, visualizing and obtaining parameters used to study nano-particle sintering using the phase field model.
Sintering is a process of compacting two or more particle powder into single solid at high temperatures, but below the particles’ melting points. This process is used in powder metallurgy to produce high strength and durable materials. To gain understand into how this process takes place in a nano-scale, a phase field modeling approach is used to study the evolution of the microstructures during sintering. Phase field is a set of values of a single variable called order parameter that represent the entire state of the microstructure during sintering[1]. Phase field has proven to be a good model in  modelling sintering because it enables the computational simulation of microstructure evolution without making restrictive assumptions. This feature of the phase field model  makes it an attractive approach to study nano-scale sintering. Computational simulation of nano-scale sintering is presented.

#### Goal
The motivation of this project is to understand the temporal evolution of nano-scale sintering.
The project objective is to successfully perform a 3D simulation and visualization of nano-scale particle sintering using the Phase Field Model.
#### Approach
**Mass Transport Phenomena**
* Grain boundary diffusion
* Surface diffusion
* Volume diffusion

**Methods**
* Phase field model
* Finite difference method

### Result
![GIF](https://github.com/CleverChuk/ICES/blob/master/images/simulation.gif)

### Summary
**Accomplishments**
1. Temporal growth of nano particle sintering was successfully simulated yielding insights into how particle microstructures evolve during sintering.
2. Most important observation is change in density of individual particles as it sinters but the overall density of the system remains relatively constant as it evolves.

**Acknowledgement**
*Department of Mechanical Engineering
*Institute for Computational Engineering and Sciences
*Dr M. Cullinan
*Georgina Obehi Dibua

### References
[1] Qin, R. S., and H. K. Bhadeshia. "Phase Field Method." Review. n.d.: n. pag. Web. 20 May 2017.

# Code Blocks
#### Simulation
The simulation was written by Georgina Obehi Dibua. I added some routines to make it work the purpose of parameter value search. The important.

Modifiers | Methods and Description
--- | --- 
**double** | **inRange**(double low, double high, double seed)
**void** | **resetParams**()

Markdown | Less
--- | --- | ---
*Still* | `renders`
1 | 2 

#### Visualization


#### Parameter Search




