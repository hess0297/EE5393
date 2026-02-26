# Problem 1
## A
Code was initially written with ChatGPT. The user then edited the code manually and with the help of ChatGPT.
### choose2(n)
Calculates the cominatorial factor for second-order reactions involving two identical reactants. This calculates how many distinct pairs of a species exists for use in propensity calculation.
### one_trial()
Runs a singular trial of N_STEPS of the set of reactions firing. For loop calculates the propensities for each reaction and then chooses a random number between 0 and the sum of propensities. It then goes through a running sum of the propensities to determine which reaction fires. This is done for N_STEPS iterations. The outcomes are then checked.
### main()
Runs one_trial() for TRIALS and prints the ratio of times each outcome was hit for the amount of TRIALS.

## B
Code was initially written with ChatGPT. The user then edited the code manually and with the help of ChatGPT.
### choose2()
See A
### step()
Calculates propensities of each reaction. Checks if any reactions can fire. Chooses a random number between 0 and the sum of propensities and then determines which reaction fires.
### run_one()
Starts from initial conditions. Perform step() for STEPS. For each iteration of step() the molecule counts are calculated and the current state is updated. The final state is returned.
### main()
Generates a random seed. Final states for TRIALS are set to zero. Updates final states by running run_one() for the amount of trials. The mean and variance for the final states for all trials are calculated and printed

# Problem 2
Code was initially written with ChatGPT. The user then edited the code manually and with the help of ChatGPT.
### nCk()
Computes the number of ways to choose distinct sets of molecules with the total number of molecules and the size of each set
### parse_stoich()
Reads the reaction file and creates a dictionary for each reactant with the reactants and products.
### load_reactions()
Uses the output from parse_stoich() to create a dictionary of reactions with reactants, products, and rates for each reaction
### load_initial_counts()
Reads the input file and determines the intial count for each specie.
### propensity()
Computes the gillespie propensity for one reaction
### classify()
Determines whether a system has reached a terminal fate
### run_one()
Gathers a copy of the initial molecule counts and sets the MOI value. Checks if a terminal state has already been reached. A for loop is created to run until MAX_STEPS or MAX_TIME has been reached. For each step the propensities of each reaction is calculated and creates a sum. It breaks if the sum is 0 and no reactions can fire. It then determines the time until the next reaction using Gillespie's theorem and chooses what reaction fires using similar principle with a random number between 0 and the sum of propensities as used in Problem 1. Stoichiometry is then applied to determine the state after that reaction. Finally, it is checked if a terminal fate has been reached. If the time or step limits has been reached then the "neither" is returned as no terminal fate was reached. 
### main()
Creates a random seed and determines the file path. The reactions and intial molecule counts are then read from the file. For each MOI value, TRIALS_PER_MOI trials are ran and it is determined if a terminal fate has was reached. Then for each MOI value the ratio of each terminal fate is calculated and printed.

# Problem 3
## A
A stoichiometric simulation was created using ChatGPT following the same structure as used in Problem 2 with the reaction network outlined in EE5393_HW1_3A.md. The code was initially created with ChatGPT and further changes were made manually and with the help of ChatGPT. It was found that in such a simulation to ensure accurate computations with the chemical reaction networks. Reaction rates were tuned with the help of ChatGPT to ensure proper outcomes.
## B
Stoichiometric and continuous simulations could not accurately simulate the chemical reaction network outlined in EE5393_HW1_3A.md. A deterministic simulation was created to mathematically prove this CRN using ChatGPT. An explanation of the chemical reaction network is also provided in EE5393_HW1_3A.md.