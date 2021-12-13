# FSO system simulation in MATLAB

Writer: Mojtaba Mansour Abadi

This simulation code uses the models and simulation MATLAB code from my PhD thesis.
Description: This code is simulating a free-space optical communication system. The simulation parameters are defined in 'GlobalParameters.m' file.
To perform the simulation run 'FSO_System.m'. The simulation can be done for various parameter sweeps. In the 'Main Loop' section you can follow the steps from generating random bits to extracting the bits from the received signal and calculating BER.
Due to an unknown memory error, instead of using function in the code, I have used a mechanism similar to MACRO definition in C/C++ language. Therefore, whenever needed, a separate file containing the codes to perform a task is called with no argument/return values. This way all the variables defined in this file will be available to the code in the separate file.
After each run, depending on the channel condition picked in 'GlobalParameters.m' file, a graph is created with a file name reflecting the channel condition.

For more information about models refer to my thesis [A hybrid free space optics/radio frequency antenna - design and evaluation](http://nrl.northumbria.ac.uk/36012/) or the book [Optical Wireless Communications: System and Channel Modelling with MATLAB](https://uk.mathworks.com/academia/books/optical-wireless-communications-ghassemlooy.html).

A brief document (.\Documents\FSO Simulation.pdf) is also included to explain the equations and mathematics used in the simulation. 

![Screenshot](Screenshot.jpg)

### Note: to run the code, you need to have all folders except "documents" on your local storage.
