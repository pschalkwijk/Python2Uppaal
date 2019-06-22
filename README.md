# Thesis
Repository for my master thesis on automating optimal scheduler synthesis using traffic abstractions based on Timed Automata.

*(Scheduling of Event-Triggered Networked Control Systems using Timed Game Automata, Dieky Adzkiya and Manuel Mazo, Jr.)*

In the *matlab* folder, an edited version of C. Hop's code to generate abstractions for traffic models can be found.
The main.m file contains all parameters to describe a system and find the abstraction.

Changes w.r.t. C. Hop's code:

- Implemented a bisection search to find the lower and upper bounds on triggering time per section
- The bounds and reachability analysis are done in a function and using a parallel loop to speed up the process.
- The result is saved to a .mat file automatically

In the *python* folder, a python implementation of timed automata and the abstractions is found.
The abstractions can be created using ETCTime by G.Gleizer or by loading a .mat file as saved by the matlab code

To recreate the network and control loop as in the paper mentioned above, the ControlLoop and Network class are implemented
These inherit from the TGA class, and can export as an Uppaal XML diagram.
To do this an adapted (partial) version of pyuppaal is included in the ta, which these classes use.

An example is included in test.py, which will be updated as the code progresses.
Currently it:
 - Creates two Matlab traffic abstractions from .mat files
 - Creates two control loops
 - Creates a network
 - Exports the three combined in a NTA to a Uppaal XML diagram
 - Creates a query to find a strategy that avoids 'losing'
 - Execute verifyta to check this query and write the strategy to ./strat/{filename}
 
 To use these files you need to have some python packages installed. I've used conda to generate an environment, and 
 in conda/conda-thesis.yml you can find the exported environment. 
 
NOTE:
- I'm using uppaal stratego 4.1.20-5 for linux, and finding a strategy can (and probably will) take a lot of time and memory.
- Uppaal Stratego can be found at http://people.cs.aau.dk/~marius/stratego/download.html
- The matlab code is not optimised in any other way than adding a BSS instead of increasing along the line to find the bounds and adding some paralellisation. 
