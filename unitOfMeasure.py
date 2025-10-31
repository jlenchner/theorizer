""""The base UofM class for managing units of measure.  See also the BaseUnit class
    in the file baseUnit.py. the currently defined BaseUnits are m,kg,s,mol,A,cd,K
    (for meter, kilogram, second, mole, Ampere, candela and Kelvin, respectively).
    SUofMs that are not BaseUnits are products of BaseUnits and their reciprocals.
    In addition to supporting explicitly constructed OofMs, this class also handles
    the inferring or guessing of UofMs based on the supplied name of a Variable,
    Derivative or Constant. Changing the set of inferred UofMs can be accomplished
    by altering the parallel arrays: ALL_MEASURED_QUANTITIES and ALL_UNITS.
"""

# Author: Jonathan Lenchner (lenchner@us.ibm.com)
#
# License: BSD 3-Clause

from sympy import *
from baseUnit import *

class UofM:  
    m,kg,s,mol,A,cd,K = base_units('m,kg,s,mol,A,cd,K')
    BASE_MEASURED_QUANTITY = ["l", "m", "t", "n", "i", "I", "T"]
    BASE_UNITS = [1,m,kg,s,mol,A,cd,K]
    BASE_UNITS_TO_DESC_DICT = {1: "unity",
                               m: "meter",
                              kg: "kilogram",
                              s: "second",
                              mol: "mole",
                              A: "Ampere",
                              cd: "candela",
                              K: "Kelvin"}
    ALL_MEASURED_QUANTITIES = ["l", "d", "x", "y", "z", "m", "t", "n", "i", "I", "T", "F", "v", "a", "W","p", "P",
                               "c","G", "h","h-bar", "E", "k", "f", "e", "pi"]
    ALL_UNITS = [m,m,m,m,m,kg,s,mol,A,cd,K, kg*m/(s*s), m/s, m/(s*s), kg*m*m/(s*s), kg*m/s, kg/(m*s*s),
                                m/s, m**3/(kg*s*s), kg*m/s, kg*m/s, kg*m*m/(s*s), kg*m*m/(s*s*K), 1/s, A*s, 1]

    ALL_MEASURED_QUANTITIES += ["theta", "sinTheta", "cosTheta", "eTheta"]
    ALL_UNITS += [1, 1, 1, 1]
    
    def __init__(self, units):
        """
            Ordinarily BaseUnits are created at the time of construction of a Variable, Derivative
            or Constant. However it is also possible to create derived UofMs and utilize the derived
            UofMs. For example one could define a Newton via N = UofM(UofM.kg*UofM.m/UofM.s),
            and then, rather than using F = Variable('F'), or F = Variable('F', UofM(UofM.kg*UofM.m/UofM.s)),
            one could accomplish the same thing with F = Variable ('F', N) and reuse N in other
            Variable/Deriviatif/Constant or UofM definitions.
        """
        self._units = units #this can be a Mul, Pow or BaseUnit
        if units != 1:  #if not dimensionless
            for sym in self._units.free_symbols:
                if sym not in UofM.BASE_UNITS:
                    UofM.BASE_UNITS.append(sym)

    def getUnits(self):
        return self._units.free_symbols

    @classmethod
    def SetBaseUnits(cls, symbols):
        UofM.BASE_UNITS = symbols

    @classmethod
    def GuessUofM(cls, measuredQtyDesc):
        if measuredQtyDesc in UofM.ALL_MEASURED_QUANTITIES:
            which = UofM.ALL_MEASURED_QUANTITIES.index(measuredQtyDesc)
            return UofM(UofM.ALL_UNITS[which])
        else:
            for i in range(len(UofM.ALL_MEASURED_QUANTITIES)):
                if measuredQtyDesc.startswith(UofM.ALL_MEASURED_QUANTITIES[i]):
                    return UofM(UofM.ALL_UNITS[i])

        return None

    def __eq__(self, UofM):
        return self._units == UofM._units

    def __hash__(self):
        return hash(self._units)

    def __str__(self):
        #return str(Poly(self._units))
        if isinstance(self, BaseUnit):
            return self.name
        if isinstance(self._units, BaseUnit):
            return self._units.name
        elif isinstance(self._units, Mul): #can have more than 2 args!!
            res = ""
            for i in range(len(self._units.args)):
                if i == 0:
                    res = str(self._units.args[i])
                else:
                    res += '*' + str(self._units.args[i])
            return res
        elif isinstance(self._units, Pow):
            return str(self._units.args[0]) + '**' + str(self._units.args[1])
        else:
            return str(self._units)





