from random import (random, uniform, choice, randint)
import csv
from math import (floor, ceil)
import matplotlib.pyplot as plt


#global data variable
data = []

#Chromosome class which contains parameters and fitness for each chromosome
class Chromosome:
    __slots__ = ['parameters', 'fitness']
    def __init__(self, parameters):
        self.parameters = parameters
        self.fitness = fitness(parameters)


#returns the mutated gene
#uniform mutation
def mutate(gene, minimumParameter, maximumParameter):
    #select random mutation index
    mutationIndex = randint(0, len(gene.parameters)-1)
    #mutate current index (assign new random value from interval [minimumParameter, maximumParameter])
    gene.parameters[mutationIndex] = uniform(minimumParameter, maximumParameter)
    #update fitness
    gene.fitness = fitness(gene.parameters)
    return gene


#returns two children chromosomes given two parents
#single-point technique
def mate(parent1, parent2):
    #select random point
    point = randint(0, len(parent1.parameters)-1)
    #add first part till point from the first parent and second part from the second parent
    child1 = parent1.parameters[:point]
    child1.extend(parent2.parameters[point:])
    #add first part till point from the second parent and second part from the first parent
    child2 = parent2.parameters[:point]
    child2.extend(parent1.parameters[point:])
    #return new chromosomes made from children parameters
    return Chromosome(child1), Chromosome(child2)


#return the fittest chromosome out of randomly selected
def tournamentSelection(population, tournamentSize):
    #select a random chromosome from population
    fittest = choice(population)
    #for the tournament size
    for i in range(tournamentSize):
        #select a random chromosome from population
        selected = choice(population)
        #if it is fitter(fitness value is less) than previously selected, let it survive(save it as the fittest)
        if selected.fitness < fittest.fitness: fittest = selected #chem menshe, tem silnee(loss->0)
    #return the fittest chromosome after tournament selection
    return fittest


#returns a new population sorted by fitness
def evolve(population, crossover, elitism, mutation, tournamentSize, minimumParameter, maximumParameter):
    #given elitism, copy the last(best fitness) population*elitism number of chromosomes to the new population
    index = floor(len(population)*(1.0-elitism))
    new_population = population[index:]
    while index > 0:
        #with probability of crossover and if index is not equal to 1 (because two elements are copied)
        if random() <= crossover and index != 1:
            #select two parents using tournament selection algorithm
            (parent1, parent2) = (tournamentSelection(population, tournamentSize), tournamentSelection(population, tournamentSize))
            #mate them using single-point technique
            children = mate(parent1, parent2)
            #for each child(out of two)
            for child in children:
                #either mutate and add them to new population
                if random() <= mutation:
                    new_population.append(mutate(child, minimumParameter, maximumParameter))
                #or just add them
                else:
                    new_population.append(child)
            #two chromosomes were added to the new population
            index -= 2
        else:
            #either mutate chromosome from old population add it to new population
            if random() <= mutation:
                new_population.append(mutate(population[index], minimumParameter, maximumParameter))
            #or just add it
            else:
                new_population.append(population[index])
            #only one chromosome was added to the new population
            index -= 1
    #return new population sorted by fitness
    return list(sorted(new_population, key = lambda f : f.fitness))


#measures fitness function given theta, which is basically mean squared error
def fitness(theta):
    sumLoss = 0
    for datum in data:
        sumLoss += (datum[0] - sum([a*b for a, b in zip(datum[1:], theta[1:])]) - theta[0])**2
    return (sumLoss/len(data))


#returns a random chromosome
def generateRandomChromosome(numberOfFeatures, minimumParameter, maximumParameter):
    gene = []
    #append as many random numbers in interval [minimumParameter, maximumParameter] as number of features + 1 for the bias term
    for parameter in range(numberOfFeatures+1):
        gene.append(uniform(minimumParameter, maximumParameter))
    return Chromosome(gene)


#generates random population and returns the list of population sorted by fitness
def generateInitialPopulation(populationSize, numberOfFeatures, minimumParameter, maximumParameter):
    population = []
    for chromosome in range(populationSize):
        population.append(generateRandomChromosome(numberOfFeatures, minimumParameter, maximumParameter))
    return list(sorted(population, key = lambda f : f.fitness))


#returns data given file name and delimiter
#data should be in this format: y feature feature feature ...
def readData(filename, delimit):
    dataFile = open(filename, 'r')
    dataReader = csv.reader(dataFile, delimiter=delimit)
    d = []
    for row in dataReader:
        d.append([float(i) for i in row])
    return d


def main():
    #read data to global variable so that it can be accessed from other functions
    global data
    data = readData("data.csv", ',')
    #number of features in data
    numberOfFeatures = len(data[0])-1
    #genetic algorithm parameters
    populationSize = 1000
    maximumNumberOfGenerations = 100
    tournamentSize = 3
    crossover = 0.8
    elitism = 0.1
    mutation = 0.3
    #specify the range of parameters of linear regression
    minimumParameter = 0
    maximumParameter = 100
    #uncomment this if you want to see a plot
    """
    x = []
    y = []
    for i in data:
        y.append(i[0])
        x.append(i[1])
    plt.ion()
    plt.scatter(x, y, color = 'green')
    """
    #generate initial population
    population = generateInitialPopulation(populationSize, numberOfFeatures, minimumParameter, maximumParameter)
    #for each generation
    for i in range(maximumNumberOfGenerations):
        #and this
        """
        yy = []
        for ybl in x:
            yy.append(population[0].parameters[0] + population[0].parameters[1] * ybl)
        plt.plot(x, yy, color='red', linewidth=2)
        plt.scatter(x, y, color = 'green')
        plt.draw()
        plt.pause(0.05)
        """
        #print current best fitness(loss)
        print("Generation %d | Cost: %f" % (i, population[0].fitness))
        #evolve the current population
        population = evolve(population, crossover, elitism, mutation, tournamentSize, minimumParameter, maximumParameter)
    #print the best parameters after all generations
    print('Parameters:')
    print(population[0].parameters[1:])
    print('Bias term')
    print(population[0].parameters[0])


if __name__ == '__main__':
    main()
