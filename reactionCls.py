class rxn(object):
    """docstring for rxn."""
    def __init__(self, name, rateConst, initConc):
        self.name = name
        self.rateConst = rateConst
        self.initConc = initConc
    def rxnEq(self):
        r = self.rateConst * self.initConc
        return(r)

class difeq(object):
    """docstring for difeq."""
    def __init__(self, name,rxns):
        self.name = name
        self.rxns = rxns
    # eq creates the dif eq by adding the elements of list rxns, which contains objects created using rxn
    def eq(self):
        sol = 0
        for i in range(0, len(self.rxns)):
            sol = sol + self.rxns[i]
        return(sol)


'''
        if(self.type == 'MM'):
            self.rateConst = [input('Please enter reaction constants')]
        elif(self.type == "MA"):
            self.rateConst = input('Please enter reaction constants')
        elif(self.type == 'transcription'):
            self.type = input('please enter whatever')
        else:
            print("Please enter MM or MA")
'''
const = [0.1,0.2,0.2,0.3,0.4]
concentration = [2,3,4,3,1,7]

reactions = []
reactions.append(rxn('first', const[0], concentration[0]))
reactions[0] = reactions[0].rxnEq()
concentration[0] = difeq('stuff', reactions[0]).eq()
print(concentration[0])

f = difeq('bla', concentration)
print(f.eq())

def sol()

'''
f = rxn('bla' ,  const[1], concentration[0])
f.rxnEq()
print(f.rateConst)


'''
