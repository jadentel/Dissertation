
import os, sys

project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if project_root not in sys.path:
    sys.path.insert(0, project_root)


import unittest
import random
import copy
import main

from termcolor import colored
from Conformation import Conformation
from Protein import Protein
from Population import Population
from main import calculation
from visual_utils import MutationVisualizer

# Direction constants.
FORWARD = 0
LEFT = -1
RIGHT = 1
SEQUENCE = "BBWWBWWBWWBWWBWWBWWBWWBB"
OPTIMAL_ENCODING = "LRRLRRLRRLFFLRRLRRLRRL"
OPTIMAL_FITNESS = -9
MUT_PROB = 0.19
CROSS_PROB = 0.77
POP_SIZE = 500


class TestProtein(unittest.TestCase):
    def test_protein_properties(self):
        p = Protein(SEQUENCE)
        length = p.getLength()
        first = p.getNth(0)
        last = p.getNth(-1)
        print(f"\n[Protein] length={length}, first={first}, last={last}")
        self.assertEqual(length, len(SEQUENCE))
        self.assertEqual(first, "B")
        self.assertEqual(last, "B")


class TestConformation(unittest.TestCase):
    def setUp(self):
        random.seed(42)
        self.prot = Protein(SEQUENCE)
        self.conf = Conformation(self.prot, set())

    def test_random_conformation_validity(self):
        self.conf.generate_random_conformation(valid=True)
        self.assertTrue(self.conf.isValid())

    def test_encoding_length(self):
        self.assertEqual(len(self.conf.get_encoding()), len(SEQUENCE) - 2)

    def test_get_conformation_string_format(self):
        self.conf.generate_random_conformation(valid=True)
        s = self.conf.getConformationString()
        self.assertEqual(len(s), len(SEQUENCE) - 2)
        for ch in s:
            self.assertIn(ch, "FLR")

    def test_calc_validity_overlap_detection(self):
        prot = Protein("BBBBB")
        conf = Conformation(prot, set())
        conf.encoding = [RIGHT, RIGHT, RIGHT]
        conf.calculate_validity()
        self.assertFalse(conf.isValid())

    def test_calc_fitness(self):
        self.conf.generate_random_conformation(valid=True)
        self.conf.calculate_fitness()
        self.assertIsInstance(self.conf.get_fitness(), int)

    def test_absolute_positions(self):
        self.conf.generate_random_conformation(valid=True)
        self.conf.calculate_absolute_position()
        self.assertEqual(len(self.conf.absPositions), len(SEQUENCE))

    def test_abs_at(self):
        self.conf.generate_random_conformation(valid=True)
        self.conf.calculate_absolute_position()
        for i in range(self.conf.getLength()):
            pos = self.conf.getAbsAt(i)
            self.assertIsInstance(pos, tuple)
            self.assertEqual(len(pos), 2)

    def test_generation_methods(self):
        gen0 = self.conf.get_generation()
        self.conf.olden()
        self.assertEqual(self.conf.get_generation(), gen0 + 1)

    def test_mutations_visualizer(self):
        original = copy.deepcopy(self.conf)
        self.conf.generate_random_conformation(valid=True)
        self.conf.mutate_directed(1.0)
        self.conf.calculate_validity(); self.conf.calculate_fitness()
        MutationVisualizer.print_conformation_pair(original, self.conf, "Directed Mutation")

    def test_corner_flip_mutation(self):
        # Test corner flip mutation and print
        original = copy.deepcopy(self.conf)
        self.conf.generate_random_conformation(valid=True)
        self.conf.calculate_validity(); self.conf.calculate_fitness()
        self.conf.mutate_corner_flip(1.0)
        self.conf.calculate_validity(); self.conf.calculate_fitness()
        MutationVisualizer.print_conformation_pair(original, self.conf, "Corner Flip Mutation")

    def test_crankshaft_mutation(self):
        # Test crankshaft mutation and print
        original = copy.deepcopy(self.conf)
        self.conf.generate_random_conformation(valid=True)
        self.conf.calculate_validity(); self.conf.calculate_fitness()
        self.conf.mutate_crankshaft(1.0)
        self.conf.calculate_validity(); self.conf.calculate_fitness()
        MutationVisualizer.print_conformation_pair(original, self.conf, "Crankshaft Mutation")


class TestPopulation(unittest.TestCase):
    def setUp(self):
        random.seed(42)
        self.prot = Protein(SEQUENCE)
        self.pop = Population(POP_SIZE, self.prot, MUT_PROB, CROSS_PROB)

    def test_population_initialization(self):
        count = len(self.pop.individuals)
        print(f"\n[Population] size={count}")
        self.assertEqual(count, POP_SIZE)

    def test_fittest_selection(self):
        fittest = self.pop.get_fittest()
        fitnesses = [indiv.get_fitness() for indiv in self.pop.individuals]
        print(f"\n[Population] fittest fitness={fittest.get_fitness()}")
        self.assertTrue(all(fittest.get_fitness() <= f for f in fitnesses))

    def test_insertable_check(self):
        clone = copy.deepcopy(self.pop.individuals[0])
        ok = self.pop.is_insertable(clone)
        print(f"\n[Population] isInsertable(same) = {ok}")
        self.assertFalse(ok)

    def test_tournament_selection(self):
        sel = self.pop.tournament_select()
        print(f"\n[Population] tournamentSelect -> fitness={sel.get_fitness()}")
        self.assertIn(sel, self.pop.individuals)

    def test_crossover_function(self):
        pre = self.pop.get_fittest().get_fitness()
        self.pop.crossover()
        post = self.pop.get_fittest().get_fitness()
        print(f"\n[Population] crossover: before={pre}, after={post}")
        self.assertLessEqual(post, pre)

    def test_crossover_visualization(self):
        # Select two real parents
        random.seed(42)
        p1 = copy.deepcopy(self.pop.tournament_select())
        p2 = copy.deepcopy(self.pop.tournament_select())
        # Guarantee distinct
        if p1.get_encoding() == p2.get_encoding():
            p2.mutate_directed(1.0)
            p2.calculate_validity(); p2.calculate_fitness()

        # Perform crossover
        child = Conformation.crossover(p1, p2, self.pop.collisionSet)
        child.calculate_validity(); child.calculate_fitness()

        # Determine crossover index
        enc1 = p1.get_encoding()
        encC = child.get_encoding()
        cx = next((i for i,(a,b) in enumerate(zip(encC,enc1)) if a!=b), len(encC))

        print(f"\n--- Crossover at index {cx} ---")

        # Generic printer
        def print_conf(conf, color, label):
            print(f"\n{label} conformation ({label} color):")
            conf.calculate_absolute_position()
            xs = [x for x,y in conf.absPositions]
            ys = [y for x,y in conf.absPositions]
            lowX, highX = min(xs), max(xs)
            lowY, highY = min(ys), max(ys)
            grid = [[' ']*((highX-lowX)*2+1) for _ in range((highY-lowY)*2+1)]
            last = conf.absPositions[0]
            for idx,(x,y) in enumerate(conf.absPositions):
                nx, ny = (x-lowX)*2, (y-lowY)*2
                char = 'H' if conf.protein.getNth(idx)=='B' else 'P'
                grid[ny][nx] = colored(char, color)
                if idx>0:
                    lx,ly = last
                    lx2, ly2 = (lx-lowX)*2, (ly-lowY)*2
                    conn = colored('-', color); vert = colored('|', color)
                    if lx2 < nx: grid[ny][nx-1] = conn
                    if lx2 > nx: grid[ny][nx+1] = conn
                    if ly2 < ny: grid[ny-1][nx] = vert
                    if ly2 > ny: grid[ny+1][nx] = vert
                last = (x,y)
            for row in grid: print(''.join(row))

        # Child printer
        def print_child(conf, cx):
            print("\nChild conformation (red=Parent1, blue=Parent2):")
            conf.calculate_absolute_position()
            xs = [x for x,y in conf.absPositions]
            ys = [y for x,y in conf.absPositions]
            lowX, highX = min(xs), max(xs)
            lowY, highY = min(ys), max(ys)
            grid = [[' ']*((highX-lowX)*2+1) for _ in range((highY-lowY)*2+1)]
            last = conf.absPositions[0]
            for idx,(x,y) in enumerate(conf.absPositions):
                nx, ny = (x-lowX)*2, (y-lowY)*2
                char = 'H' if conf.protein.getNth(idx)=='B' else 'P'
                color = 'red' if idx-2 < cx else 'blue'
                grid[ny][nx] = colored(char, color)
                if idx>0:
                    lx,ly = last
                    lx2, ly2 = (lx-lowX)*2, (ly-lowY)*2
                    conn = colored('-', color); vert = colored('|', color)
                    if lx2 < nx: grid[ny][nx-1] = conn
                    if lx2 > nx: grid[ny][nx+1] = conn
                    if ly2 < ny: grid[ny-1][nx] = vert
                    if ly2 > ny: grid[ny+1][nx] = vert
                last = (x,y)
            for row in grid: print(''.join(row))

        print_conf(p1, 'blue', 'Parent1')
        print_conf(p2, 'red', 'Parent2')
        print_child(child, cx)

        self.assertLessEqual(child.get_fitness(), min(p1.get_fitness(), p2.get_fitness()))


class TestCalculation(unittest.TestCase):
    def test_final_result_reasonable(self):
        random.seed(42)
        main.switch_minen = OPTIMAL_FITNESS
        main.switch_max_evaluations = 10000
        prot = Protein(SEQUENCE)
        pop = Population(POP_SIZE, prot, MUT_PROB, CROSS_PROB)
        calculation(pop)
        best = pop.get_fittest()
        fit = best.get_fitness()
        enc = best.getConformationString()
        print(f"\n[Calculation] final fitness={fit}, encoding={enc}")
        self.assertLessEqual(fit, OPTIMAL_FITNESS)

if __name__ == "__main__":
    unittest.main()
