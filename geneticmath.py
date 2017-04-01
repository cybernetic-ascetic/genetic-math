from random import *

#the interpretations of the genes, stored as a key-value map
geneInterpretations = {"0000": 0, "0001": 1, "0010": 2, "0011": 3,
                       "0100": 4, "0101": 5, "0110": 6, "0111": 7,
                       "1000": 8, "1001": 9, "1010": "+", "1011": "-",
                       "1100": "*", "1101": "/", "1110": "_", "1111": "_"}

class Organism:
    def __init__(self, chromosome, target):
        self.target = target
        self.chromosome = chromosome

        self.update()

        #probability of a gene mutating
        self.mutationProb = 0.5

    def __str__(self):
        expression = ""
        for symbol in self.expression:
            expression += str(symbol)

        rexpression = ""
        for symbol in self.rawExpression:
            rexpression += str(symbol)
            
        return "Gens: " + self.chromosome + " Expr: " + expression + " Rexp: " + rexpression + " Val: " + str(self.value) + " Fitn: " + str(self.fitness)

    def update(self):
        self.rawExpression = self.makeRawExpression()
        self.expression = self.makeExpression()
        self.value = self.evalExpression(self.expression[:])
        self.fitness = self.fitnessFunction(self.value, self.target)

    def evalExpression(self, expression):
        if len(expression) == 0:
            return 0

        if len(expression) == 1:
            return expression[0]

        value = expression.pop()
        operator = expression.pop()

        if operator == "+":
            value = self.evalExpression(expression) + value
        if operator == "-":
            value = self.evalExpression(expression) - value
        if operator == "*":
            value = self.evalExpression(expression) * value
        if operator == "/":
            value = self.evalExpression(expression) / value
        if operator == "^":
            value = self.evalExpression(expression) ** value

        return value

    def fitnessFunction(self, value, target):
        if value == target:
            return 1000000000
        else:
            return abs(1 / (target - value))

    def makeExpression(self):
        chromosome = self.chromosome[:]
        expression = []
        genes = []

        #split the chromosome into groups of 4 gbits, stored as a list
        for k in range(0, len(chromosome), 4):
            genes.append(chromosome[k:k+4])

        #converts each gene into a gene interpretation, stored as a string
        for gene in genes:
            expression.append(geneInterpretations[gene])
        
        prevSymbol = None
        newExpression = []

        #goes through the expression and excludes characters that are illegal
        for symbol in expression:
            if prevSymbol is None:
                if type(symbol) is int:
                    newExpression.append(symbol)
                    prevSymbol = symbol
            if type(prevSymbol) is int:
                if type(symbol) is str:
                    newExpression.append(symbol)
                    prevSymbol = symbol
            if type(prevSymbol) is str and prevSymbol != "_" and prevSymbol != "/":
                if type(symbol) is int:
                    newExpression.append(symbol)
                    prevSymbol = symbol
            if type(prevSymbol) is str and prevSymbol == "/":
                if type(symbol) is int and symbol != 0:
                    newExpression.append(symbol)
                    prevSymbol = symbol

        if len(newExpression) > 0:
            if type(newExpression[-1]) is str:
                newExpression.pop()

        return newExpression

    #changes the last n genes to the string of genetic material received from another organism
    def mate(self, geneticMaterial):
        segmentDifference = len(self.chromosome) - len(geneticMaterial)
        initSegment = self.chromosome[0: segmentDifference]

        self.chromosome = initSegment + geneticMaterial
        self.update()

    #every gene in every chromosome flips bits with probability of mutationProb
    def mutate(self):
        newChromosome = ""
        
        for gene in self.chromosome:
            randNum = random()

            if randNum <= self.mutationProb:
                if gene == "0":
                    newChromosome += "1"
                else:
                    newChromosome += "0"

            else:
                newChromosome += gene

        self.chromosome = newChromosome
        self.update()

    #creates the unmodified arithmetic expression from the gene interpretation
    #stored as a string
    def makeRawExpression(self):
        chromosome = self.chromosome[:]
        expression = []
        genes = []

        #split the chromosome into groups of 4 gbits, stored as a list
        for k in range(0, len(chromosome), 4):
            genes.append(chromosome[k:k+4])

        #converts each gene into a gene interpretation, stored as a string
        for gene in genes:
            expression.append(geneInterpretations[gene])

        return expression
            

class Population:
    def __init__(self, populationSize, chromosomeLength, target):
        self.population = []
        self.chromosomeLength = chromosomeLength
        self.target = target
        self.generation = 1
        self.mateProb = 1

        for k in range(populationSize):
            organism = Organism(self.randomChromosome(chromosomeLength * 4), self.target)
            self.population.append(organism)

    def computeGeneration(self):
        bestOrganism = self.maxFitness()
        solutionOrganism = self.checkAnswers()

        if not solutionOrganism is None:
            print("Solution found at Generation", self.generation, "!!!")
            print(solutionOrganism)
            return False

        print("Generation", self.generation, " ", bestOrganism)

        self.mutatePopulation()

        randNum = random()

        if randNum <= self.mateProb:
            self.selectMatingPair()

        self.generation += 1

        return True

    #looks to see if any organism has a solution, and returns the first one that does
    #if no organism has a solution, return null
    def checkAnswers(self):
        for organism in self.population:
            if organism.fitness == 1000000000:
                return organism

        return None

    #induces combination on the selected organisms
    def matePair(self, boy, girl):
        mateLength = randint(1, self.chromosomeLength - 1)
        fatherGenes = boy.chromosome[mateLength:]
        motherGenes = girl.chromosome[mateLength:]

        boy.mate(motherGenes)
        girl.mate(fatherGenes)

    #returns the organism with the highest fitness
    def maxFitness(self):
        bestOrganism = self.population[0]
        
        for organism in self.population:
            if organism.fitness > bestOrganism.fitness:
                bestOrganism = organism

        return bestOrganism

    #runs the mutate method on every organism
    def mutatePopulation(self):
        for organism in self.population:
            organism.mutate()

    #creates a random chromosome at a specified length, stored as a string
    def randomChromosome(self, length):
        chromosome = ""

        for k in range(length):
            chromosome += str(randint(0,1))

        return chromosome

    #selects two organisms from the population to mate
    def selectMatingPair(self):
        selectionWeights = []
        fitnessSum = self.totalFitness()

        #computes the normalized weighting of fitness, stored as a list
        for organism in self.population:
            selectionWeights.append(organism.fitness / fitnessSum)

        selectionList = []
        sumOfPreviousWeights = 0

        for weight in selectionWeights:
            selectionList.append(weight + sumOfPreviousWeights)
            sumOfPreviousWeights += weight

        #randomly choose first organism from the selection list
        randomNum = random()

        for selection in range(len(selectionList)):
            if randomNum <= selectionList[selection]:
                father = self.population[selection]

        #randomly choose second organism from the selection list
        randomNum = random()

        for selection in range(len(selectionList)):
            if randomNum <= selectionList[selection]:
                mother = self.population[selection]

        self.matePair(father, mother)

    #sums the fitness of every organism in the population, stored as a float
    def totalFitness(self):
        result = 0
        for organism in self.population:
            result += organism.fitness

        return result

target = 1000
pop = Population(50, 15, target)
print("Target is", target)

runSimulation = True

while(runSimulation):
    runSimulation = pop.computeGeneration()

