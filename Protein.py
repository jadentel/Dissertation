class Protein:
    def __init__(self, sequence: str):
        self.sequence = sequence

    def getNth(self, i: int) -> str:
        return self.sequence[i]

    def getLength(self) -> int:
        return len(self.sequence)