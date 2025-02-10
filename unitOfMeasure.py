from sympy import *
from baseUnit import *

class UofM:  #Need methods to get UofM for a term and to determine if an equation is dimensionally consistent (perhaps not in this class0
   #How to deal with a N (Newton) = kg*m/(sec*sec)?  UofMs need to be polys not symbols!!
    m,kg,s,mol,A,cd,K = base_units('m,kg,s,mol,A,cd,K')
    BASE_MEASURED_QUANTITY = ["l", "m", "t", "n", "i", "I", "T"]
    BASE_UNITS = [m,kg,s,mol,A,cd,K]
    BASE_UNITS_TO_DESC_DICT = {m: "meter",
                              kg: "kilogram",
                              s: "second",
                              mol: "mole",
                              A: "Ampere",
                              cd: "candela",
                              K: "Kelvin"}
    ALL_MEASURED_QUANTITIES = ["l", "d", "x", "y", "z", "m", "t", "n", "i", "I", "T", "F", "v", "a", "W","p", "P",
                               "c","G", "h","h-bar", "E", "k", "e", "pi"]
    ALL_UNITS = [m,m,m,m,m,kg,s,mol,A,cd,K, kg*m/(s*s), m/s, m/(s*s), kg*m*m/(s*s), kg*m/s, kg/(m*s*s),
                                m/s, m**3/(kg*s*s), kg*m/s, kg*m/s, kg*m*m/(s*s), kg*m*m/(s*s*K), A*s, m/m]


    def __init__(self, units):
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
        return str(Poly(self._units))











