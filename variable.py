from unitOfMeasure import *

class Variable:
    VAR_TO_UofM_DICT = dict()
    DEFAULT_MAPPING = dict()

    def __init__(self, sym, u_of_m = None, isConstant = False, derivativeOrder = 0, depVariable = None, indepVariable = None):
        self._sym = sym
        if derivativeOrder == 0:
            ableToGuessDerivativeInfo = self.guessDerivativeInfo()
        if u_of_m is None:
            self.u_of_m = UofM.GuessUofM(measuredQtyDesc=sym.name)
        else:
            self._u_of_m = u_of_m
        Variable.VAR_TO_UofM_DICT.update({sym: u_of_m})
        self._isConstant = isConstant
        if derivativeOrder > 0 or not ableToGuessDerivativeInfo:
            self._derivativeOrder = derivativeOrder
            self._depVariable = depVariable
            self._indepVariable = indepVariable


    def guessDerivativeInfo(self):
        ableToGuess = False
        if self._sym.name.startswith('d'):
            d_count = self._sym.name.count('d')
            two_count = self._sym.name.count('2')
            three_count = self._sym.name.count('3')
            if d_count == 2:
                if two_count == 0 and three_count == 0:
                    l1 = self._sym.name.find('d')
                    l2 = self._sym.name.find('d', l1+1)
                    self._derivativeOrder = 1
                    self._depVariable = Variable(Symbol(self._sym.name[l1+1:l2]))
                    self._indepVariable = Variable(Symbol(self._sym.name[l2:]))
                    ableToGuess = True
                elif two_count == 2:
                    l1 = self._sym.name.find('d2')
                    l2 = self._sym.name.find('d', l1 + 1)
                    l3 = self._sym.name.find('2', l2 + 1)
                    self._derivativeOrder = 2
                    self._depVariable = Variable(Symbol(self._sym.name[l1 + 2:l2]))
                    self._indepVariable = Variable(Symbol(self._sym.name[l2:l3]))
                    ableToGuess = True
                elif three_count == 2:
                    l1 = self._sym.name.find('d3')
                    l2 = self._sym.name.find('d', l1 + 1)
                    l3 = self._sym.name.find('3', l2 + 1)
                    self._derivativeOrder = 3
                    self._depVariable = Variable(Symbol(self._sym.name[l1 + 2:l2]))
                    self._indepVariable = Variable(Symbol(self._sym.name[l2:l3]))
                    ableToGuess = True

        return ableToGuess

    @classmethod
    def SetAll(cls, var_to_UofM_dict):  #perhaps do this with parallel arrays
        Variable.VAR_TO_UofM_DICT = var_to_UofM_dict

    @classmethod
    def GetAll(cls):
        allVars = []
        keys = Variable.VAR_TO_UofM_DICT.keys()
        for key in keys:
            var = Variable(key, Variable.VAR_TO_UofM_DICT.get(key))
            allVars.append(var)

        return allVars

    def isConstant(self):
        return self._isConstant

    def isDerivate(self):
        return self._derivativeOrder > 0

    def getDerivativeOrder(self):
        return self._derivativeOrder

    def getIndepVariable(self):
        return self._indepVariable
