import random
from typing import List, Tuple
import Protein

# Direction constants.
FORWARD = 0
LEFT = -1
RIGHT = 1



class Conformation:
    # Tracks total number of energy evaluations
    energyEvalSteps = 0

    def __init__(self, protein = None, setOfPoints = None):
        self.protein = protein
        self.setOfPoints = setOfPoints
        if protein is not None:
            self.length = protein.getLength()
        else:
            self.length = 0

        # The encoding holds directions for positions 2..n, so there are length-2 entries.
        if self.length >= 2:
            self.encoding = [0] * (self.length - 2)
        else:
            self.encoding = []

        # Absolute positions for each residue (x, y).
        self.absPositions = [(0, 0)] * self.length
        self.fitness = 0
        self.generation = 0
        self.validState = False
        
        # Generate initial conformation if both protein and collision set are provided
        if protein is not None and setOfPoints is not None:
            self.generate_random_conformation(valid=True)

    @classmethod
    def crossover(cls, p1, p2, setOfPoints = None):
        # One-point crossover: combine encodings from two parents
        new_conf = cls()
        new_conf.setOfPoints = setOfPoints
        new_conf.protein = p1.protein
        new_conf.length = p1.length
        new_conf.generation = (p1.generation + p2.generation) // 2 + 1
        new_conf.encoding = [0] * (new_conf.length - 2)
        new_conf.absPositions = [(0, 0)] * new_conf.length

        if new_conf.length - 2 > 0:
            randI = random.randint(0, new_conf.length - 3)
        else:
            randI = 0

        # First part from p1, second from p2.
        for i in range(0, randI):
            new_conf.encoding[i] = p1.encoding[i]
        for i in range(randI, new_conf.length - 2):
            new_conf.encoding[i] = p2.encoding[i]

        if new_conf.setOfPoints is not None:
            new_conf.calculate_validity()
        return new_conf


    def get_encoding(self):
        return self.encoding

    def getProtein(self):
        return self.protein

    def isValid(self):
        return self.validState

    def getLength(self):
        return self.length

    def calculate_fitness(self):
        # Calculate energy (fitness) based on number of hydrophobic (B) contacts
        Conformation.energyEvalSteps += 1
        fitness = 0
        # For each amino acid, count hydrophobic contacts.
        for i in range(self.length):
            if self.protein.getNth(i) == 'B':
                ori = self.absPositions[i]
                for j in range(i + 2, self.length):
                    if self.protein.getNth(j) == 'B':
                        dest = self.absPositions[j]
                        distX = abs(ori[0] - dest[0])
                        distY = abs(ori[1] - dest[1])
                        if (distX == 1 and distY == 0) or (distX == 0 and distY == 1):
                            fitness -= 1
        self.fitness = fitness

    def get_fitness(self):
        return self.fitness

    def getConformationString(self):
        result = ""
        for d in self.encoding:
            if d == FORWARD:
                result += "F"
            elif d == LEFT:
                result += "L"
            elif d == RIGHT:
                result += "R"
            else:
                result += "?"
        return result

    def generate_random_conformation(self, valid: bool = False):
        # Initialize encoding with random directions.
        for i in range(self.length - 2):
            self.encoding[i] = random.randint(-1, 1)
        if valid:
            self.calculate_validity()
            while not self.validState:
                self.mutate(0.1)
                self.calculate_validity()

    def calculate_validity(self):
        self.calculate_absolute_position()
        seen = set()
        self.validState = True
        for pos in self.absPositions:
            if pos in seen:
                self.validState = False
                return
            seen.add(pos)
        if self.setOfPoints is not None:
            self.setOfPoints.clear()
            self.setOfPoints.update(seen)

  
    # randomly chooses one of three mutations.
    def mutate(self, probability):
        op_choice = random.random()
        if op_choice < 0.33:
            self.mutate_directed(probability)
        elif op_choice < 0.66:
            self.mutate_corner_flip(probability)
        else:
            self.mutate_crankshaft(probability)

    # 1. Directed local perturbation mutation.
    def mutate_directed(self, probability):
        for i in range(self.length - 2):
            if self.randomFloat() <= probability:
                current = self.encoding[i]
                # Choose a new value different from the current one.
                alternatives = [d for d in [-1, 0, 1] if d != current]
                self.encoding[i] = random.choice(alternatives)

    # 2. Corner flip mutation: flip a corner by inverting the turning move.
    def mutate_corner_flip(self, probability):
        self.calculate_absolute_position()
        # Consider residues 1 through length-2 as potential corners.
        for i in range(1, self.length - 1):
            if self.randomFloat() <= probability and self.isCorner(i):
                # The responsible encoding is at index i-1.
                if self.encoding[i-1] == LEFT:
                    self.encoding[i-1] = RIGHT
                elif self.encoding[i-1] == RIGHT:
                    self.encoding[i-1] = LEFT

    # 3. Crankshaft mutation: rotate a short segment if the endpoints form a square.
    def mutate_crankshaft(self, probability):
        if self.length < 4:
            return
        if self.randomFloat() <= probability:
            self.calculate_absolute_position()
            # Choose a random segment (indices i and i+1 will be rotated).
            i = random.randint(1, self.length - 3)
            p0 = self.absPositions[i-1]
            p3 = self.absPositions[i+2]
            # Check if p0 and p3 are diagonal neighbors (i.e. form a square).
            if abs(p0[0] - p3[0]) == 1 and abs(p0[1] - p3[1]) == 1:
                # The two candidates to complete the square.
                candidate1 = (p0[0], p3[1])
                candidate2 = (p3[0], p0[1])
                # For this example, choose candidate1 for residue i and candidate2 for residue i+1.
                new_p1 = candidate1
                new_p2 = candidate2
                # Compute the move from p0 to new_p1.
                move1 = (new_p1[0] - p0[0], new_p1[1] - p0[1])
                if i-1 == 0:
                    prev_abs_dir = 0  # Starting direction: up.
                else:
                    prev_move = (self.absPositions[i-1][0] - self.absPositions[i-2][0],
                                 self.absPositions[i-1][1] - self.absPositions[i-2][1])
                    prev_abs_dir = self.absolute_direction(prev_move)
                new_move1_abs = self.absolute_direction(move1)
                new_rel1 = self.relative_move(prev_abs_dir, new_move1_abs)
                # Compute the move from new_p1 to new_p2.
                move2 = (new_p2[0] - new_p1[0], new_p2[1] - new_p1[1])
                new_move2_abs = self.absolute_direction(move2)
                new_rel2 = self.relative_move(new_move1_abs, new_move2_abs)
                # Update the encoding for the affected segment.
                self.encoding[i-1] = new_rel1
                self.encoding[i] = new_rel2

    # Helper: determine if residue at index i is a corner.
    def isCorner(self, i):
        self.calculate_absolute_position()
        if i <= 0 or i >= self.length - 1:
            return False
        v1 = (self.absPositions[i][0] - self.absPositions[i-1][0],
              self.absPositions[i][1] - self.absPositions[i-1][1])
        v2 = (self.absPositions[i+1][0] - self.absPositions[i][0],
              self.absPositions[i+1][1] - self.absPositions[i][1])
        return v1 != v2

    # Compute absolute positions from the relative encoding.
    def calculate_absolute_position(self):
        self.absPositions = [None] * self.length
        x, y = 0, 0
        self.absPositions[0] = (x, y)
        y = 1
        self.absPositions[1] = (x, y)
        lastDirection = 0  # Starting "up" direction.
        for i in range(self.length - 2):
            d = self.encoding[i]
            if lastDirection == 0:  # up
                if d == FORWARD:
                    y += 1
                elif d == RIGHT:
                    x += 1
                    lastDirection = 2
                elif d == LEFT:
                    x -= 1
                    lastDirection = 3
            elif lastDirection == 1:  # down
                if d == FORWARD:
                    y -= 1
                elif d == RIGHT:
                    x -= 1
                    lastDirection = 3
                elif d == LEFT:
                    x += 1
                    lastDirection = 2
            elif lastDirection == 2:  # right
                if d == FORWARD:
                    x += 1
                elif d == RIGHT:
                    y -= 1
                    lastDirection = 1
                elif d == LEFT:
                    y += 1
                    lastDirection = 0
            elif lastDirection == 3:  # left
                if d == FORWARD:
                    x -= 1
                elif d == RIGHT:
                    y += 1
                    lastDirection = 0
                elif d == LEFT:
                    y -= 1
                    lastDirection = 1
            self.absPositions[i + 2] = (x, y)

    # Helper: convert a move vector (dx, dy) to an absolute direction code.
    def absolute_direction(self, move):
        if move == (0, 1):
            return 0  # up
        elif move == (1, 0):
            return 1  # right
        elif move == (0, -1):
            return 2  # down
        elif move == (-1, 0):
            return 3  # left
        else:
            raise ValueError("Invalid move vector: " + str(move))

    # Helper: given a previous absolute direction and a new absolute direction,
    # return the relative move (-1 for left turn, 0 for forward, 1 for right turn).
    def relative_move(self, prev_abs_dir, new_abs_dir):
        diff = (new_abs_dir - prev_abs_dir) % 4
        if diff == 0:
            return FORWARD
        elif diff == 1:
            return RIGHT
        elif diff == 3:
            return LEFT
        else:
            # A 180° turn should not occur in a valid conformation.
            return FORWARD

    def get_generation(self):
        return self.generation

    def olden(self):
        self.generation += 1

    def getStatusString(self):
        return f"Fitness: {self.fitness}   Generation: {self.generation}"

    def printAsciiPicture(self):
        xs = [pos[0] for pos in self.absPositions]
        ys = [pos[1] for pos in self.absPositions]
        lowestX = min(xs)
        highestX = max(xs)
        lowestY = min(ys)
        highestY = max(ys)
        width = (highestX - lowestX) * 2 + 1
        height = (highestY - lowestY) * 2 + 1
        grid = [[' ' for _ in range(width)] for _ in range(height)]
        lastX, lastY = self.absPositions[0]
        for idx, (x, y) in enumerate(self.absPositions):
            normX = (x - lowestX) * 2
            normY = (y - lowestY) * 2
            acid = self.protein.getNth(idx)
            if acid == 'B':
                grid[normY][normX] = 'H'
            elif acid == 'W':
                grid[normY][normX] = 'P'
            else:
                grid[normY][normX] = acid
            if idx > 0:
                if lastX > x:
                    if normX + 1 < width:
                        grid[normY][normX + 1] = '-'
                elif lastX < x:
                    if normX - 1 >= 0:
                        grid[normY][normX - 1] = '-'
                if lastY < y:
                    if normY - 1 >= 0:
                        grid[normY - 1][normX] = '|'
                elif lastY > y:
                    if normY + 1 < height:
                        grid[normY + 1][normX] = '|'
            lastX, lastY = x, y
        for row in grid:
            print(''.join(row))



    def randomFloat(self):
        return random.random()

    def getAbsAt(self, i):
        return self.absPositions[i]
