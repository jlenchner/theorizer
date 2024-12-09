from sympy import *

class UofM:  #Need methods to get UofM for a term and to determine if an equation is dimensionally consistent (perhaps not in this class0
   #How to deal with a N (Newton) = kg*m/(sec*sec)?  UofMs need to be polys not symbols!!
    m,kg,s,mol,A,cd,K = symbols('m,kg,s,mol,A,cd,K')
    BASE_MEASURED_QUANTITY = ["l", "m", "t", "n", "i", "I", "T"]
    BASE_UNITS = [m,kg,s,mol,A,cd,K]
    BASE_UNITS_TO_DESC_DICT = {m: "meter",
                              kg: "kilogram",
                              s: "second",
                              mol: "mole",
                              A: "Ampere",
                              cd: "candela",
                              K: "Kelvin"}
    ALL_MEASURED_QUANTITY = ["l", "d", "x", "y", "z", "m", "t", "n", "i", "I", "T", "F"]
    ALL_UNITS = [Mul(m),Mul(m),Mul(m),Mul(m),Mul(m),Mul(kg),Mul(s),Mul(mol),Mul(A),Mul(cd),Mul(K), Mul(kg*m/(s*s))]


    def __init__(self, units):
        self._units = Mul(units)
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
        if measuredQtyDesc in UofM.ALL_MEASURED_QUANTITY:
            which = UofM.ALL_MEASURED_QUANTITY.index(measuredQtyDesc)
            return UofM.ALL_UNITS[which]
        else:
            for i in range(len(UofM.ALL_MEASURED_QUANTITY)):
                if measuredQtyDesc.startswith(UofM.ALL_MEASURED_QUANTITY[i]):
                    return UofM.ALL_UNITS[i]

        return None


