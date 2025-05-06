import random
from typing import Set

from Conformation import Conformation, Protein

class Population:
    def __init__(self, size, prot, mutProb, crossProb):
        # Conformation of parents used during crossover
        self.parent1 = None
        self.parent2 = None
        self.mutProb = mutProb
        self.crossProb = crossProb
        self.protein = prot
        self.size = size

        # This set is passed to Conformation for collision checking
        self.collisionSet: Set = set()
        # To track already encountered conformations to ensure theyre unique
        self.setOfConformations = set()
        # List to hold all individuals (Conformation objects)
        self.individuals = []

        print("Generate Population:")
        i = 0
        # Keep generating until the population is filled
        while i < self.size:
            # Create a temporary conformation(Conformation's constructor generates a random valid conformation when passed protein and collision set)
            temp = Conformation(self.protein, self.collisionSet)
            temp.calculate_fitness()
            # Only insert if fitness is not zero and is unique
            if temp.get_fitness() != 0 and temp.getConformationString() not in self.setOfConformations:
                self.individuals.append(temp)
                self.setOfConformations.add(temp.getConformationString())
                print(f"{i}.", end="", flush=True)
                i += 1

        # Initialize the fittest individual as the first one and then update
        self.theFittest = self.individuals[0]
        self.set_fittest()
        print()

    def is_insertable(self, candidate):
        # Check if a conformation is new to the population
        confStr = candidate.getConformationString()
        if confStr not in self.setOfConformations:
            self.setOfConformations.add(confStr)
            return True
        else:
            return False

    def set_fittest(self):
        # Identifies the individual with the lowest energy (best fitness)
        for indiv in self.individuals:
            if indiv.get_fitness() < self.theFittest.get_fitness():
                self.theFittest = indiv

    def get_fittest(self):
        return self.theFittest


    # New tournament selection method
    def tournament_select(self, tournament_size: int = 3) -> Conformation:
        # Randomly select a subset of individuals
        participants = random.sample(self.individuals, tournament_size)
        # For protein folding, a lower fitness (more negative energy) is better.
        best = min(participants, key=lambda indiv: indiv.get_fitness())
        return best

    # Updated crossover method to use tournament selection
    def crossover(self):
        # Select parents using tournament selection instead of roulette wheel selection.
        self.parent1 = self.tournament_select()
        if self.crossProb < random.random():
            return  # skip crossover based on probability

        self.parent2 = self.tournament_select()
        if self.crossProb < random.random():
            return  # skip crossover based on probability

        # Create two children via crossover (recombination).
        child1 = Conformation.crossover(self.parent1, self.parent2, self.collisionSet)
        child2 = Conformation.crossover(self.parent2, self.parent1, self.collisionSet)

        # Mutate child1, recalc validity and fitness.
        child1.mutate(self.mutProb)
        child1.calculate_validity()
        if child1.isValid():
            child1.calculate_fitness()
            if self.is_insertable(child1):
                # Replace the less fit parent with child1 if fitter.
                if child1.get_fitness() < self.parent1.get_fitness():
                    self.parent1.__dict__.update(child1.__dict__)
                elif child1.get_fitness() < self.parent2.get_fitness():
                    self.parent2.__dict__.update(child1.__dict__)

        # Mutate child2, recalc validity and fitness.
        child2.mutate(self.mutProb)
        child2.calculate_validity()
        if child2.isValid():
            child2.calculate_fitness()
            if self.is_insertable(child2):
                if child2.get_fitness() < self.parent2.get_fitness():
                    self.parent2.__dict__.update(child2.__dict__)
                elif child2.get_fitness() < self.parent1.get_fitness():
                    self.parent1.__dict__.update(child2.__dict__)

        # Update the fittest individual if any parent improved.
        if self.parent1.get_fitness() < self.theFittest.get_fitness():
            self.theFittest = self.parent1
        if self.parent2.get_fitness() < self.theFittest.get_fitness():
            self.theFittest = self.parent2