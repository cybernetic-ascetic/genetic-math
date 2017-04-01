from random import *

def product(_list):
    result = 1
    for number in _list:
        result *= number
    return result

def summ(_list):
    result = 0
    for number in _list:
        result += number
    return result

class Organism:
    def __init__(self, chromosome="random"):
        if chromosome == "random":
            self.chromosome = self.randomChromosome()
        else:
            self.chromosome = chromosome
        
        self.fitness = self.fitnessFunction()

    def __str__(self):
        return str(self.chromosome) + " " + str(self.fitness)

    def cardNotation(self):
        sumPile = ""
        mulPile = ""

        for gene in range(len(self.chromosome)):
            if self.chromosome[gene] == 0:
                sumPile += str(gene + 1) + "+"
            else:
                mulPile += str(gene + 1) + "*"

        return sumPile[:-1] + " = 36 " + mulPile[:-1] + " = 360"

    def fitnessFunction(self):
        addGenes = []
        mulGenes = []
        
        for gene in range(len(self.chromosome)):
            if self.chromosome[gene] == 0:
                addGenes.append(gene + 1)
            else:
                mulGenes.append(gene + 1)

        addValue = summ(addGenes)
        mulValue = product(mulGenes)

        return (abs(36 - addValue) + abs(360 - mulValue))

    def mutate(self):
        mutatedGene = randint(0,9)

        if self.chromosome[mutatedGene] == 0:
            self.chromosome[mutatedGene] = 1
        else:
            self.chromosome[mutatedGene] = 0

        self.update()

    def randomChromosome(self):
        result = []

        for k in range(10):
            result.append(randint(0,1))

        return result

    def update(self):
        self.fitness = self.fitnessFunction()

class Population:
    def __init__(self, size):
        self.population = []
        self.size = size
        self.generation = 1

        for k in range(size):
            organism = Organism()
            self.population.append(organism)

    def __str__(self):
        result = "Generation " + str(self.generation) + "\n"

        for organism in self.population:
            result += str(organism) + "\n"

        return result

    def advanceGeneration(self):
        bestOrganism = self.bestOrganism()
        print(self.generation, bestOrganism)

        for organism in range(self.size):
            newOrganism = Organism(bestOrganism.chromosome[:])
            newOrganism.mutate()
            self.population[organism] = newOrganism

        self.generation += 1

    def bestOrganism(self):
        result = self.population[0]

        for organism in self.population:
            if organism.fitness < result.fitness:
                result = organism

        return result

    def averageFitness(self):
        result = 0

        for organism in self.population:
            result += organism.fitness

        return result / size

def main():
    x = Population(3)
    contSim = True

    while(contSim):
        if x.bestOrganism().fitness == 0:
            contSim = False
            print(x.bestOrganism())
            print(x.bestOrganism().cardNotation())
        else:
            x.advanceGeneration()

main()
