import cmath
import math
import operator
import random

"""
    Attempt to see if a lattice flow with values in the Eisenstein integers can be decomposed
     in a "well-behaved" way (with positive inner-product maybe).
"""


class Eisen:
    """ Class for Eisenstein integers """

    # Store Eisenstein integers as coordinates with respect to the ordered basis (1, zeta_6)
    _reCoord = 0;
    _z6Coord = 0;

    def __init__(self, a = 0, b = 0):
        """ Initialize Eisenstein integer of the form a + b*zeta_6
            Default values for a and b are 0
        """
        if not isinstance(a, int):
            raise TypeError("Eisen(a,b): expected argument a of type {}, received {}: '{}'".format(type(a), a))
        if not isinstance(b, int):
            raise TypeError("Eisen(a,b): expected argument b of type {}, received {}: '{}'".format(type(b), b))
        self._reCoord = a
        self._z6Coord = b


    def __str__(self):
        """ Returns a string corresponding to the Eisenstein integer
            The string has form a + b*zeta_6
        """
        if self._z6Coord == 0:
            return str(self._reCoord)
        elif self._reCoord == 0:
            return '{0}*zeta_6'.format(str(self._z6Coord))
        else:
            return '{0} + {1}*zeta_6'.format(str(self._reCoord), str(self._z6Coord))


    def __complex__(self):
        """ Returns the Eisenstain integer as a complex number """
        return complex(self._reCoord + 0.5 * self._z6Coord, math.sqrt(3) / 2 * self._z6Coord)


    def __eq__(self, other):
        """ Equality operator
            Will work with integers and complex numbers
        """
        if isinstance(other, self.__class__):
            return (self._reCoord == other._reCoord) and (self._z6Coord == other._z6Coord)
        elif isinstance(other, int):
            return (self._reCoord == other) and (self._z6Coord == 0)
        elif isinstance(other, complex):
            return complex(self) == other
        else:
            return False


    def __ne__(self, other):
        """ Negative equality operator """
        return not (self == other)


    """ Reverse equality/inequality opeartor """
    __req__ = __eq__
    __rne__ = __ne__


    def __neg__(self):
        """ Negation Operator """
        return Eisen(-self._reCoord, -self._z6Coord)


    def __add__(self, other):
        """ Addition operator
            Will accept integers as other argument
        """
        if isinstance(other, self.__class__):
            return Eisen(self._reCoord + other._reCoord, self._z6Coord + other._z6Coord)
        elif isinstance(other, int):
            return Eisen(self._reCoord + other, self._z6Coord)
        else:
            raise TypeError("unsupported operand type(s) for +: '{}' and '{}'".format(self.__class__, type(other)))


    # Define "reverse addition" operator (for things like 2 + Eisen(1,0))
    __radd__ = __add__


    def __sub__(self, other):
        """ Subtraction operator
            Will accept integers as other argument
        """
        if isinstance(other, self.__class__):
            return Eisen(self._reCoord - other._reCoord, self._z6Coord - other._z6Coord)
        elif isinstance(other, int):
            return Eisen(self._reCoord - other, self._z6Coord)
        else:
            raise TypeError("unsupported operand type(s) for -: '{}' and '{}'".format(self.__class__, type(other)))


    def __rsub__(self, other):
        """ Reverse subtraction operator
            Will accept integers as other argument
        """
        if isinstance(other, self.__class__):
            return Eisen(other - self._reCoord._reCoord, -self._z6Coord - other._z6Coord)
        elif isinstance(other, int):
            return Eisen(other - self._reCoord, -self._z6Coord)
        else:
            raise TypeError("unsupported operand type(s) for -: '{}' and '{}'".format(type(other), self.__class__))


    def __mul__(self, other):
        """ Multiplication operator
            Will accept integers as other argument
        """
        if isinstance(other, self.__class__):
            a = self._reCoord
            b = self._z6Coord
            c = other._reCoord
            d = other._z6Coord
            return Eisen(a*c - b*d, a*d + c*b + b*d)
        elif isinstance(other, int):
            return Eisen(self._reCoord * other, self._z6Coord * other)
        else:
            try:
                return other.__rmul__(self)
            except:
                raise TypeError("unsupported operand type(s) for *: '{}' and '{}'".format(self.__class__, type(other)))


    """ Reverse multiplication operator """
    __rmul__ = __mul__


    def __iadd__(self, other):
        """ Extended addition operator """
        if isinstance(other, self.__class__):
            self._reCoord += other._reCoord
            self._z6Coord += other._z6Coord
        elif isinstance(other, int):
            self._reCoord += other
        else:
            raise TypeError("unsupported operand type(s) for +=: '{}' and '{}'".format(self.__class__, type(other)))
        return self


    def __isub__(self, other):
        """ Extended subtraction operator """
        if isinstance(other, self.__class__):
            self._reCoord -= other._reCoord
            self._z6Coord -= other._z6Coord
        elif isinstance(other, int):
            self._reCoord -= other
        else:
            raise TypeError("unsupported operand type(s) for -=: '{}' and '{}'".format(self.__class__, type(other)))
        return self


    def __imul__(self, other):
        """ Extended multiplication operator """
        if isinstance(other, self.__class__):
            a = self._reCoord
            b = self._z6Coord
            c = other._reCoord
            d = other._z6Coord
            self._reCoord = a*c - b*d
            self._z6Coord = a*d + c*b + b*d
        elif isinstance(other, int):
            self._reCoord *= other
            self._z6Coord *= other
        else:
            raise TypeError("unsupported operand type(s) for *=: '{}' and '{}'".format(self.__class__, type(other)))
        return self


    def re(self):
        """ Returns the real part of the Eisenstein integer """
        return self._reCoord + 0.5 * self._z6Coord


    def conjugate(self):
        """ Returns the complex conjugate of Eisenstein integer """
        return Eisen(self._reCoord + self._z6Coord, -self._z6Coord)


    def distance_to_origin(self):
        """ Returns the distance from the origin as the length of a shortest path
            in the lattice of Eisenstein integers

            If the coordinates are both positive or both negative we add them and take
             make the result positive, otherwise we take the largest of the two coordintes
        """

        return max(abs(self._reCoord + self._z6Coord), abs(self._reCoord), abs(self._z6Coord))



class EisenFlow:
    """ Class for an eisenstein integer valued flow """

    # Flow is stored as a list of of _numElts Eisenstein integers
    #  (ground set is assumed to be {1,...,_numElts})
    _numElts = 0
    _flow = []


    def __init__(self, n, init_list = None):
        """ Create a new flow on ground set {1,...,n},
            Initialize flow with the given list if provided
        """
        self._numElts = n
        if init_list is None:
            self._flow = [Eisen() for i in range(n)]
        else:
            assert(len(init_list) == n)
            self._flow = list(init_list)


    def __len__(self):
        """ Return the size of the ground set as the "length" of the flow """
        return self._numElts


    def __str__(self):
        """ Returns a string corresponding to the flow represented as a row vector"""
        return '[{}]'.format(', '.join(map(str, self._flow)))


    def __add__(self, other):
        """ Addition operator """
        if isinstance(other, self.__class__):
            assert(self._numElts == other._numElts)
            newFlow = EisenFlow(self._numElts)
            newFlow._flow = map(operator.add, self._flow, other._flow)
            return newFlow
        else:
            raise TypeError("unsupported operand type(s) for +: '{}' and '{}'".format(self.__class__, type(other)))


    def __sub__(self, other):
        """ Subtraction operator """
        if isinstance(other, self.__class__):
            assert(self._numElts == other._numElts)
            newFlow = EisenFlow(self._numElts)
            newFlow._flow = map(operator.sub, self._flow, other._flow)
            return newFlow
        else:
            raise TypeError("unsupported operand type(s) for +: '{}' and '{}'".format(self.__class__, type(other)))


    def __mul__(self, other):
        """ Multiplcation operator used for scalar multiplcation """
        # No type checking, let the multiplication for Eisenstein integers do that
        newFlow = EisenFlow(self._numElts)
        newFlow._flow = [x * other for x in self._flow]
        return newFlow

    """ Reverse multiplication operator """
    __rmul__ = __mul__


    def __iadd__(self, other):
        """ Extended addition operator """
        if isinstance(other, self.__class__):
            assert(self._numElts == other._numElts)
            for i in range(self._numElts):
                self._flow[i] += other._flow[i]
            return self
        else:
            raise TypeError("unsupported operand type(s) for +=: '{}' and '{}'".format(self.__class__, type(other)))


    def __isub__(self, other):
        """ Extended subtraction operator """
        if isinstance(other, self.__class__):
            assert(self._numElts == other._numElts)
            for i in range(self._numElts):
                self._flow[i] -= other._flow[i]
            return self
        else:
            raise TypeError("unsupported operand type(s) for -=: '{}' and '{}'".format(self.__class__, type(other)))


    def __imul__(self, other):
    """ Extended multiplcation operator used for scalar multiplcation """
        self._flow = [x * other for x in self._flow]
        return self


    def iprod(self, other):
        """ Standard complex inner product """
        assert(self._numElts == other._numElts)
        result = Eisen(0,0)
        for i in range(self._numElts):
            result += self._flow[i] * other._flow[i].conjugate()
        return result


    def value(self, elt):
        """ Returns the value of the flow on element elt """
        assert(0 <= elt < self._numElts)
        return self._flow[elt]

    def is_zero(self):
        return all(x == 0 for x in self._flow)



def eisen_pow(base, exponent):
    """ Returns the power of the Eisenstein integer base, raised to the given integer exponent """
    result = Eisen(1,0)
    for i in range(exponent):
        result *= base
    return result


def flow_entry_product(f1, f2):
    assert(len(f1) == len(f2))
    return [(f1.value(i) * f2.value(i).conjugate()).re() for i in range(len(f1))]




# List of 6-th roots of unity
zeta6 = [eisen_pow(Eisen(0,1), i) for i in range(6)]


# Conisder a matroid on 9 elements with the following circuits:
#   [1, 1, 0, 1, 0, 0]
#   [0,-1,-1, 0,-1, 0]
#   [1, 1, 0, 0, 1, 1]
#   [1, 0,-1, 0, 0, 1]
#   [0, 0, 0,-1, 1, 1]
#   [1, 0,-1, 1,-1, 0]
#   [0,-1,-1,-1, 0, 1]

# Previously we had this matroid, but (1,3), (2,5), and (8,9) were seies pairs
# contracting those down to single elements we obtain the matroid above
#   1,3,6,4         [1, 0, 1, 1, 0, 1, 0, 0, 0]
#   2,5,7,4         [0, 1, 0, 1, 1, 0, 1, 0, 0]
#   1,3,8,9,7,4     [1, 0, 1, 1, 0, 0, 1, 1, 1]
#   1,3,8,9,-5,-2   [1,-1, 1, 0,-1, 0, 0, 1, 1]
#   -6,8,9,7        [0 ,0 ,0 ,0 ,0,-1, 1, 1, 1]
#   1,3,6,-7,-5,-2  [1,-1, 1, 0,-1, 1,-1, 0, 0]
#   -4,-6,8,9,-5,-2 [0,-1, 0,-1,-1,-1, 0, 1, 1]
#
# Corresponding simple flows
C = [0 for i in range(7)]
C[0] = EisenFlow(6, list(map(Eisen, [1, 1, 0, 1, 0, 0])))
C[1] = EisenFlow(6, list(map(Eisen, [0,-1,-1, 0,-1, 0])))
C[2] = EisenFlow(6, list(map(Eisen, [1, 1, 0, 0, 1, 1])))
C[3] = EisenFlow(6, list(map(Eisen, [1, 0,-1, 0, 0, 1])))
C[4] = EisenFlow(6, list(map(Eisen, [0, 0, 0,-1, 1, 1])))
C[5] = EisenFlow(6, list(map(Eisen, [1, 0,-1, 1,-1, 0])))
C[6] = EisenFlow(6, list(map(Eisen, [0,-1,-1,-1, 0, 1])))


def makeRandomFlow(complexity):
    f = EisenFlow(6)
    for i in range(complexity):
        circ = random.randint(0,6)
        phse = random.randint(0,5)
        # print(str(zeta6[phse]), str(C[circ]))
        f += zeta6[phse] * C[circ]
    return f

def agree(a, b):
    return (a * b.conjugate()).re() >= 0
    # p1 = cmath.phase(complex(a))
    # p2 = cmath.phase(complex(b))
    # return abs(p1 - p2) - 0.00001 <= math.pi / 3


def conforms(sflow, flow):
    assert(len(sflow) == len(flow))
    for i in range(len(sflow)):
        if sflow.value(i) != 0 and not (flow.value(i) != 0 and agree(sflow.value(i), flow.value(i))):
            return False
    if all((flow - sflow).value(i).distance_to_origin() >= flow.value(i).distance_to_origin() for i in range(len(flow))):
        return False
    return True


def exists_decomposition(flow, timeout, verbose = False):
    numIters = 0
    while not flow.is_zero():
        numIters += 1
        if numIters > timeout:
            print("timeout!")
            return False
        found = False
        for circ in C:
            for phase in zeta6:
                if conforms(phase * circ, flow):
                    found = True
                    flow -= phase * circ
                    if verbose:
                        print(str(phase), str(circ))
                        print("now", str(flow))
                    break
            if found:
                break
        if not found:
            return False
    return True


# print("flow", str(f))
# exists_decomposition(f, 6, verbose = True)

for i in range(1000):
    f = makeRandomFlow(4)
    # print("flow", str(f))
    if not exists_decomposition(f, 10):
        print("flow", str(f))
        exists_decomposition(f, 10, True)
        break


def test(e):
    e += 1
    return e

a = Eisen(0,1)
print(a)
print(test(a + 0))
print(a)

# for i in range(2):
#     f = makeRandomFlow(7)
#     for c in C:
#         g = c - f
#         if f.is_zero() or g.is_zero():
#             continue
#         print(str(f.iprod(g)), complex(f.iprod(g)))
#         if f.iprod(g).re() >= 0:
#             print("AHA! ---------------")
#             print(f)
#             print(g)
#             print(c)