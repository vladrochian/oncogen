# Object for storing a single substitution mutation on a protein
class Mutation:
    def __init__(self, position: int, original: str, actual: str):
        self.position = position
        self.original = original
        self.actual = actual

    def to_string(self):
        return self.original + str(self.position) + self.actual
