# Code to simulate the MMT model
Code to simulate the one-dimensional Majda-McLaughlin-Tabak (MMT) model used in 
-   Dù, R.S., Bühler, O. _The impact of frequency bandwidth on a one-dimensional model for dispersive wave turbulence_. Submitted to Journal of Nonlinear Science.

In short, the code solve the MMT model pseudospectrally using an Integrating Factor method with fourth-order-Runge–Kutta timestepping (IF-RK4). For more details, please refer to our paper.

To start, run [``outer_driver.m``](https://github.com/Empyreal092/MMT_Public/blob/main/outer_driver.m). This outer wrapper specifies some the key tunable parameteres of the MMT model, then it calls the main solver: [``MMT_solver.m``](https://github.com/Empyreal092/MMT_Public/blob/main/MMT_solver.m).

The main solver then calls functions and subscripts in the folder [``SubRoutine``](https://github.com/Empyreal092/MMT_Public/blob/main/SubRoutine). Here we list some key functions:

1. [``DNL_MMT.m``](https://github.com/Empyreal092/MMT_Public/blob/main/SubRoutine/DNL_MMT.m): This function evaluate the nonlinear part of the MMT model. It is key for the pseudospectral solver.
2. [``RuntimePlotting_ForDiss.m``](https://github.com/Empyreal092/MMT_Public/blob/main/SubRoutine/RuntimePlotting_ForDiss.m): Subscript for plotting some information (e.g.: solution in real space, action spectrum, and integral quatities) about the MMT system as it evolves in time. Turn on runtime plotting by setting the ``doruntimeplots`` parameter in [``outer_driver.m``](https://github.com/Empyreal092/MMT_Public/blob/main/outer_driver.m). 
