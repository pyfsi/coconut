

class Variable:
    def __init__(self, name, type):
        self.__name = name
        if type in ['scalar', 'vector']:
            self.__type = type
        else:
            raise ValueError(f'unknown type "{type}"')

    @property
    def name(self):
        return self.__name

    @property
    def type(self):
        return self.__type

    def __eq__(self, other):
        return (self.name == other.name and self.type == other.type)

    def __hash__(self):
        return hash((self.name, self.type))

    def __repr__(self):
        return f'Variable "{self.name}" of type {self.type}'





# create some Variable objects
a = Variable('pressure', 'vector')
b = Variable('pressure', 'scalar')
c = Variable('pressure', 'vector')

# test their equality
print(a == b)
print(a == c)

# test hashing
for var in [a, b, c]:
    print(hash(var))

# try to use Variable as key in dict
pres = Variable('pressure', 'vector')
d = {}

d['a'] = 123
d[2] = 'sth'
d[pres] = [1, 2, 3]

print(type(d))
print(d)
print(d[pres])