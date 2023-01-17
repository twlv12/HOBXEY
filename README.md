![HOBXEY](https://user-images.githubusercontent.com/79507121/212982079-f6b68cf5-3565-4e88-b70d-0c229eadcd36.png)
# HOBXEY (Hybrid Organism Behavior eXperimentation sYstem.)

This program uses genetic algorithms to simulate the evolution of organisms in a virtual world. Creatures are controlled by a neural network defined by a genome, a string of binary digits. This program uses the pygame library to visualize the world.
## Requirements

    Python 3
    Pickle
    Pygame
    Numpy

## Structure

The program is structured as follows:

    1. The program starts by defining global variables and functions that are used throughout the program. These include the types of neurons, the properties of the simulation (such as population size and mutation chance), and basic functions (such as mutation and decoding of the genome).
    2. The program then defines the classes used in the simulation. These include the Connection and Neuron classes, which are used to create the neural network of the creatures.
    3. The program then creates the World, Creature, and Brain classes, which define the simulation and its components.
    4. The program then enters the main loop, which handles the simulation's logic and visualization. The loop runs for a set number of steps (or generations) and updates the creatures' neural networks and movements based on their genomes.
    5. The program saves and loads the simulation's data to and from a file.

## Theory

The program uses a genetic algorithm to simulate the evolution of creatures in a virtual world. The creatures are controlled by a neural network, which is defined by a genome, a string of binary digits. The creatures' genomes are randomly generated at the start of the simulation, and are mutated with a certain probability at each step. The creatures' fitness is determined by their ability to move in the virtual world. The creatures with the highest fitness have a higher probability of reproducing and passing on their genome to the next generation. Over time, the creatures' genomes will evolve to optimize their movement in the virtual world, leading to the emergence of intelligent behavior.
## Usage

To run the program, simply execute the file main.py with python. The properties of the simulation can be adjusted in the global variables section at the beginning of the file. Hotkeys are available with the display enabled only. 

    R -> Reset sim.
    D -> While hovering over a creature, decode and print brain structure.
    S -> Save to file. If specific file set, use that one.
##Customization

The program can be easily modified to change the properties of the world, the behavior of the creatures and the parameters of the genetic algorithm. You can also change the behavior of the creature by changing the types of neurons and the logic of their function. To add a new neuron type, add it with its name to the corresponding type dictionary. Create a function inside the Neuron class to define the neuron's behaviour. Make sure the function name is exactly as it is in the dictionary.
