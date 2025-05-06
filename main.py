from Conformation import Conformation
from Protein import Protein
from Population import Population

# Global variables
global_fittest_ptr = None
isTerminated = False
minimum_energy = -9
max_evaluations = 100000

def calculation(pop: Population):
    global global_fittest_ptr, isTerminated
    global_fittest_ptr = pop.get_fittest()
    
    # Continue until the fittest's fitness reaches threshold or max evaluations are exceeded.
    while pop.get_fittest().get_fitness() > minimum_energy and Conformation.energyEvalSteps < max_evaluations:
        pop.crossover()
        
        # Check if an improvement has been made.
        if pop.get_fittest().get_fitness() < global_fittest_ptr.get_fitness():
            global_fittest_ptr = pop.get_fittest()
            print(global_fittest_ptr.getStatusString())
            global_fittest_ptr.printAsciiPicture()
    
    isTerminated = True

def main():
    # Seed the random number generator for reproducibility
    #random.seed(42)
    
    # Initialize a Protein object with an example sequence
    
    prot = Protein("BBWWBWWBWWBWWBWWBWWBWWBB")

    population_size = 10000
    mutation_probability = 0.05
    crossover_probability = 0.85

    
    # Create the Population object
    pop = Population(population_size, prot, mutation_probability, crossover_probability)
    
    # Run the calculation loop
    calculation(pop)
    
    # Output the final fittest individual
    fittest = pop.get_fittest()
    print("Final Fittest Individual:", fittest.getStatusString())
    fittest.printAsciiPicture()
    enc = fittest.getConformationString()
    print("Optimal encoding:", enc)
    print("Fitness:", fittest.get_fitness())

if __name__ == "__main__":
    main()



