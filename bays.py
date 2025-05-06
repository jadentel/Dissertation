from bayes_opt import BayesianOptimization
from Conformation import Conformation
from Protein import Protein
from Population import Population

# Global variables as before.
global_fittest_ptr = None
isTerminated = False
switch_enable_graphics = False  
switch_minen = -52              
switch_max_evaluations = 1000000 

def calculation(pop: Population):
    global global_fittest_ptr, isTerminated
    global_fittest_ptr = pop.get_fittest()
    
    # Continue until the fittest's fitness reaches threshold or max evaluations are exceeded.
    while pop.get_fittest().get_fitness() > switch_minen and Conformation.energyEvalSteps < switch_max_evaluations:
        pop.crossover()
        
        if pop.get_fittest().get_fitness() < global_fittest_ptr.get_fitness():
            global_fittest_ptr = pop.get_fittest()
            if not switch_enable_graphics:
                print(global_fittest_ptr.getStatusString())
                global_fittest_ptr.printAsciiPicture()
    
    isTerminated = True

def run_ga(population_size: float, mutation_probability: float, crossover_probability: float) -> float:
    """
    Run the genetic algorithm with given hyperparameters and return the negative final fitness.
    We return the negative fitness because the GA minimizes fitness and BayesianOptimization maximizes.
    """
    # Convert population_size to integer for the GA.
    pop_size = int(population_size)
    # Use the same protein sequence as before.
    prot = Protein("BBBBWWWWBBBBBBBBBBBBWWWWWWBBBBBBBBBBBBWWWBBBBBBBBBBBBWWWBBBBBBBBBBBBWWWBWWBBWWBBWWBWB")
    
    # Reset evaluation counter for fairness.
    Conformation.energyEvalSteps = 0
    # Create the Population instance with the given hyperparameters.
    pop = Population(pop_size, prot, mutation_probability, crossover_probability)
    
    # Run the GA.
    calculation(pop)
    
    fittest = pop.get_fittest()
    final_fitness = fittest.get_fitness()
    print(f"Run complete with params: pop_size={pop_size}, mut_prob={mutation_probability}, cross_prob={crossover_probability} => final fitness: {final_fitness}")
    
    # Return negative fitness so that a lower fitness (better) gives a higher objective.
    return -final_fitness

def bayesian_optimization():
    # Define the parameter bounds.
    pbounds = {
        'population_size': (500, 2000),
        'mutation_probability': (0.01, 0.4),
        'crossover_probability': (0.4, 0.9)
    }
    
    optimizer = BayesianOptimization(
        f=run_ga,
        pbounds=pbounds,
        random_state=42,
        verbose=2
    )
    
    # Run with a few initial random points and then iterations.
    optimizer.maximize(init_points=5, n_iter=10)
    return optimizer.max

def bayesian_GA():
    print("Running Bayesian Optimization for Hyperparameter Tuning")
    best = bayesian_optimization()
    best_params = best['params']
    best_score = best['target']
    print("Best hyperparameters found:", best_params)
    print("Best GA negative fitness with optimized parameters:", best_score)
    
def main():
    bayesian_GA()
    
if __name__ == "__main__":
    main()
