# Object for storing a single substitution mutation on a protein
class Mutation:
    def __init__(self, position: int, original: str, actual: str):
        self.position = position
        self.original = original
        self.actual = actual

    @staticmethod
    def parse(code: str):
        if len(code) < 3:
            raise Exception('Invalid code')

        original = code[0]
        actual = code[-1]
        try:
            position = int(code[1:-1])
        except ValueError:
            raise Exception('Invalid code')

        return Mutation(position, original, actual)

    def to_string(self):
        return self.original + str(self.position) + self.actual
