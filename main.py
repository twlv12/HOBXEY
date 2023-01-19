import random
import pickle
import time as t
import gc
import pygame
import numpy as np


gc.enable()
#this is a comment


inputNeuronTypes = {
    "0001" : "RandomInput",
    "0010" : "BlockageUp",
    "0011" : "BlockageDown",
    "0100" : "XGradient",
    "0101" : "YGradient",
    "0110" : "ChangeX",
    "0111" : "ChangeY",
    "1000" : "Density",
    "1001" : "Pheromone",
    "1010" : "BlockageLeft",
    "1011" : "BlockageRight",
    "1100" : "ConstantPos",
    "1101" : "ConstantNeg",
    }
internalNeuronTypes = {
    "0001" : "Invert",
    "0010" : "Multiply"
    }
outputNeuronTypes = {
    "0001" : "MoveX",
    "0010" : "MoveY",
    "0011" : "MoveRandom",
    "0100" : "MoveToPheromone"
    }
neuronTypeIndicators = {
    "00" : internalNeuronTypes,
    "01" : inputNeuronTypes,
    "10" : outputNeuronTypes
    }



#--------Properties-----------#
genomeLength = 6                #How many genes in a genome
worldSize = (128, 128)          #World size in cells
population = 1000                #How many creatures
windowSize = (1280, 1280)       #Window size
stepsPerGeneration = 200        #How many steps in each generation
mutationChance = 5           #Chance in percent for a mutation to occur on each creature
customGenomes = []              #Will be overwritten by save file
quickMode = False               #About 3x faster, minimal console and display
saveFile = ''    #File to save to
loadFile = 'autosave.pickle'    #File to load from
#--------Properties-----------#



#--------Basic Functions-----------#

""" Within 20x20 square at the center of the world:
if(x > (world.sizeX//2 - 10)) and (x < (world.sizeX//2 + 10)) and (y > (world.sizeY//2 - 10)) and (y < (world.sizeY//2 + 10)):
    return 1
else: return 0
"""
   
""" Within 5 cells of the boundaries:
if 4 <= x <= world.sizeX-4 and 4 <= y <= world.sizeY-4:
    return 0
else: return 1
"""
        

def isValid(x,y,world):
    if(x > (world.sizeX//2 - 20)) and (x < (world.sizeX//2 + 20)) and (y > (world.sizeY//2 - 20)) and (y < (world.sizeY//2 + 20)):
        return 1
    else: return 0


    

def remapWeight(x):
    new_value = np.interp(x, [1,15], [-1,1])
    return new_value
    
    
    

def mutate(genome):
    genes = [genome[i:i+20] for i in range(0, len(genome), 20)]
    index = random.randint(0, len(genes)-1)
    genes.pop(index)
    mutated_gene = createGene() + bin(random.getrandbits(8))[2:].zfill(8)
    genes.insert(index, mutated_gene)
    return "".join(genes)



 
def decodeGenome(input):
    genomeLen = len(input)
    genes = []
    struct = ""
    
    for i in range(0, genomeLen, 20):
        genes.append(input[i:i+20])
    print(f"""
        Current Genome -> {input} ({genomeLen//20} Genes)\n
    Breakdown:""")
    
    for gene in genes:
        try:
            st = gene[:2]; sid = gene[2:6]
            tt = gene[6:8]; tid = gene[8:12]
            weight = gene[12:]
            source = neuronTypeIndicators[st][sid]
            target = neuronTypeIndicators[tt][tid]
            
            struct += f"       {source} -> {target} : {weight}\n"
        except Exception:
            source = "Invalid"
            target = "Invalid"
            struct += f"       Invalid Gene: {gene}\n"
            
    return source, target, struct




def saveCreatures(world, file):
    genomes_data = [creatures.genome for creatures in world.creatures]
    with open(file, "wb") as f:
        pickle.dump(genomes_data, f)
def loadCreatures():
    try:
        with open(loadFile, "rb") as f:
            genomes = pickle.load(f)
        return genomes
    except Exception:
        print("No load file found.")
        return None
        


        
def createGene():
        sourceType = random.choice(["00", "01"])
        if sourceType == "00":
            sourceID = random.choice(list(internalNeuronTypes.keys()))
        elif sourceType == "01":
            sourceID = random.choice(list(inputNeuronTypes.keys()))
        source = sourceType + sourceID
        
        targetType = random.choice(["00", "10"])
        if targetType == "00":
            targetID = random.choice(list(internalNeuronTypes.keys()))
        elif targetType == "10":
            targetID = random.choice(list(outputNeuronTypes.keys()))
        target = targetType + targetID
        
        if sourceID == targetID:
            return createGene()
        return source + target




def createGenome(length):  
    genome = ""
    for i in range(0, length):
        gene = createGene() + bin(random.getrandbits(8))[2:].zfill(8)
        genome += gene
    return genome




def selectCreatures(world):
    selected = []
    mutations = []
    i = 0
    while len(world.creatures) > 0:
        creature = world.creatures.pop()
        i += 1
        if creature.canReproduce:
            x = int((1 / mutationChance) * 100) - 1
            if random.randint(0, x) == 1:
                mutatedGenome = (mutate(creature.genome))
                selected.append(mutatedGenome); mutations.append(mutatedGenome)
            else:
                selected.append(creature.genome)
        creature.kill()

    print(f"\nProcessed: {i}")
    print(f"Creatures Left: {len(world.creatures)}\n")
    
    if showDisplay: world.draw()
    world.initializePopulation(population, genomes = selected)
    
    intelligence = round(len(selected) / population, 3)
    print(f"Killed {population-len(selected)} : Intelligence {intelligence} : Mutations {len(mutations)}\n\n")
    return intelligence
#--------Basic Functions-----------#

     

#--------Classes-----------#
class Connection:
    def __init__(self, source, target, weight):
        self.source = source
        self.target = target
        self.weight = remapWeight(int(weight, 2))
        
    def update(self):
        self.target.inputValue += float(self.source.outputValue) * self.weight
        

        
        
        
class Neuron:
    def __init__(self, type, brain):
        self.inputValue = 0.0
        self.outputValue = 0.0
        self.type = type
        self.creature = brain.creature
        self.world = self.creature.world
        self.connectionsOut = []
        self.connectionsIn = []
        
    def connect(self, target, weight):
        self.connectionsOut.append(Connection(self, target, weight))
        target.connectionsIn.append(Connection(target, self, weight))
        
    def update(self):
        self.outputValue = getattr(self, self.type)(self.inputValue)
        self.inputValue = 0
        for connection in self.connectionsOut:
            connection.update()
            
    def RandomInput(self, input):
        return random.random()

    def BlockageUp(self, input):
        try:
            if self.creature.getCell(self.creature.x, self.creature.y+1).creature: return 1
            else: return 0
        except: return 1

    def BlockageDown(self,input):
        try:
            if self.creature.getCell(self.creature.x, self.creature.y-1).creature: return 1
            else: return 0
        except: return 1

    def BlockageLeft(self,input):
        try:
            if self.creature.getCell(self.creature.x-1, self.creature.y).creature: return 1
            else: return 0
        except: return 1

    def BlockageRight(self,input):
        try:
            if self.creature.getCell(self.creature.x+1, self.creature.y).creature: return 1
            else: return 0
        except: return 1

    def XGradient(self, input):
        gradient = 2 * (self.creature.x / worldSize[0] - 0.5)
        return gradient
       
    def YGradient(self, input):
        gradient = 2 * (self.creature.y / worldSize[1] - 0.5)
        return gradient

    def ChangeX(self, input):
        try:
            if self.creature.x > self.creature.lastX: return 1
            elif self.creature.x < self.creature.lastX: return -1
        except: return 0

    def ChangeY(self, input):
        try:
            if self.creature.y > self.creature.lastY: return 1
            elif self.creature.y < self.creature.lastY: return -1
        except: return 1

    def Density(self, input):
        return self.creature.getPopDensity()
       
    def ConstantPos(self, input):
        return 1

    def ConstantNeg(self, input):
        return -1

    def Invert(self, input):
        return -input

    def Multiply(self,input):
        try:
            product = self.connectionsIn[0].outputValue
            for connection in self.connectionsIn[1:]:
                product *= connection.outputValue
            return product
        except Exception: return 0

    def MoveX(self, input):
        if input >= 0.2:
            if self.creature.x + 1 < self.world.sizeX:
                self.creature.moveToCell(self.world.getCell(self.creature.x + 1, self.creature.y))
        elif input <= -0.2:
            if self.creature.x - 1 >= 0:
                self.creature.moveToCell(self.world.getCell(self.creature.x - 1, self.creature.y))
        
    def MoveY(self, input):
        if input >= 0.2:
            if self.creature.y + 1 < self.world.sizeY:
                self.creature.moveToCell(self.world.getCell(self.creature.x, self.creature.y + 1))
        elif input <= -0.2:
            if self.creature.y - 1 >= 0:
                self.creature.moveToCell(self.world.getCell(self.creature.x, self.creature.y - 1))
        
    def MoveRandom(self, input):
        if input <= -0.2 or input >= 0.2:
            cell = random.choice(self.creature.getAdjacentCells())
            self.creature.moveToCell(cell)
        
    def Pheromone(self, input):
        return self.creature.getPheromoneDensity()[0]

    def MoveToPheromone(self, input):
        if input >= 0.4 or input <= -0.4:
            self.creature.moveToCell(self.creature.getPheromoneDensity()[1])

        
        
        
            
class Brain:
    def __init__(self, genome, creature):
        self.genome = genome
        self.creature = creature
        self.neurons = []
        self.genes = [genome[i:i+20] for i in range(0, len(genome), 20)]
            
        i = 0
        for gene in self.genes:
            st = gene[:2]; sid = gene[2:6]
            tt = gene[6:8]; tid = gene[8:12]
            weight = gene[12:]
            
            sourceType = neuronTypeIndicators[st][sid]
            targetType = neuronTypeIndicators[tt][tid]
            
            if not any(neuron.type == sourceType for neuron in self.neurons):
                sourceNeuron = Neuron(sourceType, self)
                self.neurons.append(sourceNeuron)
            else:
                sourceNeuron = [neuron for neuron in self.neurons if neuron.type == sourceType][0]
            if not any(neuron.type == targetType for neuron in self.neurons):
                targetNeuron = Neuron(targetType, self)
                self.neurons.append(targetNeuron)
            else:
                targetNeuron = [neuron for neuron in self.neurons if neuron.type == targetType][0]
                
            sourceNeuron.connect(targetNeuron, weight)
            i += 1
        
    def update(self):
        for neuron in self.neurons:
            neuron.update()



            
             
class Cell:
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.creature = None
        self.pheromone = 0
    def __str__(self):
        return f"Cell({self.x},{self.y})"
            
        
            
            

class Creature:
    def __init__(self, world, genome = None):
        if genome == None:
            self.genome = createGenome(genomeLength)
        else:
            self.genome = genome
        self.world = world
        self.brain = Brain(self.genome, self)

        self.cell = random.choice(world.freeCells)
        self.density = 0
        self.x = self.cell.x
        self.y = self.cell.y
        self.canReproduce = False
        self.cell.creature = self
        world.freeCells.remove(self.cell)
        
    def update(self):
        self.brain.update()
        self.cell.pheromone += 0.01
        self.canReproduce = isValid(self.x,self.y,self.world)

    def getPopDensity(self):
        cells = self.getAdjacentCells()
        nearbyCreatures = []
        for cell in cells:
            if cell.creature:
                nearbyCreatures.append(cell.creature)
        return len(nearbyCreatures) / 8
    
    def getPheromoneDensity(self):
        adjacent_cells = self.getAdjacentCells()
        pheromones = [cell.pheromone for cell in adjacent_cells]
        if pheromones:
            avg_pheromone = sum(pheromones) / len(pheromones)
            max_pheromone_cell = max(adjacent_cells, key=lambda cell: cell.pheromone)
            return avg_pheromone, max_pheromone_cell
        else:
            return None, None

    def getAdjacentCells(self):
        adjacent_cells = []
        for i in range(self.x-1, self.x+2):
            for j in range(self.y-1, self.y+2):
                if i >= 0 and i < self.world.sizeX and j >= 0 and j < self.world.sizeY:
                    if i != self.x or j != self.y:
                        adjacent_cells.append(self.world.getCell(i, j))
        return adjacent_cells
        
    def moveToCell(self, cell):
        if not cell.creature:
            lastCell = self.cell
            newCell = cell
            self.cell = newCell
            self.x = newCell.x; self.y = newCell.y
            newCell.creature = self
            lastCell.creature = None
            self.world.freeCells.append(lastCell)
            self.world.freeCells.remove(newCell)
            self.world.movements += 1

    def getColour(self):
        random.seed(self.genome)
        r = random.randint(0, 255)
        g = random.randint(0, 255)
        b = random.randint(0, 255)
        return (r, g, b)

    def kill(self):
        self.cell.creature = None
        self.world.freeCells.append(self.cell)
        del self
    
    def toDict(self):
        return {"genome": self.genome, 
                "x": self.x, 
                "y": self.y, 
                "cell": self.cell,
                "brain": self.brain}
    
    
        
class World:
    def __init__(self, size, population, genomes = None):
        self.sizeX = size[0]
        self.sizeY = size[1]
        self.cells = []
        self.creatures = []
        self.movements = 0
        self.cells = [[Cell(i,j) for j in range(self.sizeY)] for i in range(self.sizeX)]
        self.freeCells = [cell for row in self.cells for cell in row]

        if showDisplay:
            self.initializeDisplay()
        self.initializePopulation(population, genomes)
                
    def getCell(self, x, y):
       if x >= 0 and x < self.sizeX and y >= 0 and y < self.sizeY:
           return self.cells[x][y]

    def initializeDisplay(self):
        pygame.init()
        self.screen = pygame.display.set_mode((windowSize[0], windowSize[1]))
        pygame.display.set_caption("Genetic Drift Simulation")
        self.clock = pygame.time.Clock()
        self.cellSize = (windowSize[0] / self.sizeX, windowSize[1] / self.sizeY)

    def initializePopulation(self, population, genomes):
        self.population = population
        if genomes == None:
            for i in range(self.population):
                self.creatures.append(Creature(self))
        else:
            for i in range(self.population):
                genome = random.choice(genomes)
                self.creatures.append(Creature(self, genome))

    def update(self):
        for creature in self.creatures:
            creature.update()
        if showDisplay:
            self.draw()
        
    def draw(self):
        self.screen.fill((0, 0, 0))
        for row in self.cells:
            for cell in row:
                if cell.creature:
                    pygame.draw.rect(self.screen, cell.creature.getColour(), (cell.x * self.cellSize[0], cell.y * self.cellSize[1], self.cellSize[0], self.cellSize[1]))
                else:
                    blue = cell.pheromone * 50
                    if blue > 50: blue = 50
                    green = isValid(cell.x,cell.y,self) * 100
                    pygame.draw.rect(self.screen, (0,green,blue), (cell.x * self.cellSize[0], cell.y * self.cellSize[1], self.cellSize[0], self.cellSize[1]))
        pygame.display.flip()

#--------Classes-----------#            


usingLoadFile = False; usingSaveFile = False; showDisplay = True; usingCustomGenomes = False
if len(customGenomes) > 0: usingCustomGenomes = True
if len(saveFile) > 0: usingSaveFile = True
if len(loadFile) > 0: usingLoadFile = True
if quickMode: showDisplay = False


#--------Loop-----------#
def main():
    print("Initializing World...")
    if usingLoadFile: world = World(worldSize, population, loadCreatures())
    elif usingCustomGenomes: world = World(worldSize, population, customGenomes)
    else: world = World(worldSize, population)
    print("World Initialized.")
    print("------------------\nStarting Master Update Loop...\n")
    
    counter = 1; gen = 1; iq = "N/A"
    running = True
    lastSave = 0
    while running:
        if counter > stepsPerGeneration:
            print(f"Gen: {gen}")
            counter = 1
            gen += 1
            iq = selectCreatures(world)
            for row in world.cells:
                for cell in row:
                    cell.pheromone = 0
            saveCreatures(world, "autosave.pickle")
            print("Autosaved.")

        #if gen > 10 and iq < 0.5:
        #    main()
            
        world.movements = 0
        world.update()
        
        if showDisplay:
            pos = pygame.mouse.get_pos()
            cell = world.getCell(int(pos[0] / world.cellSize[0]), int(pos[1] / world.cellSize[1]))
            if cell.creature: genome = cell.creature.genome
            else: genome = "N/A"

            for event in pygame.event.get():
                if event.type == pygame.QUIT:
                    running = False

                if event.type == pygame.KEYDOWN:
                    if event.key == pygame.K_r:
                        main()
                    if event.key == pygame.K_d:
                        print(decodeGenome(genome)[2])
                        t.sleep(3)
                    if event.key == pygame.K_s:
                        if lastSave is not gen:
                            lastSave = gen
                            if usingSaveFile:
                                saveCreatures(world, saveFile)
                                print("Saved to file: " + saveFile)
                            else:
                                current_time_struct = t.gmtime(t.time())
                                time = t.strftime("%Y-%m-%d_%H-%M", current_time_struct)
                                saveFile = f"usersave{time}.pickle"
                                saveCreatures(world, saveFile)
                                print("Saved to file: " + saveFile)
                        else:
                            print("Please wait before saving.")

        if not quickMode:
            if showDisplay: print(f"Step {str(counter).zfill(3)} | Gen {gen} | Win {str(iq).zfill(3)} | Moves {str(world.movements).zfill(3)} | Cell ({str(cell.x).zfill(3)},{str(cell.y).zfill(3)}) | Ph {str(round(cell.pheromone, 3)).zfill(4)} | Pop {population}")
            else: print(f"Step {str(counter).zfill(3)} | Gen {gen} | Win {str(iq).zfill(3)} | Moves {str(world.movements).zfill(3)} | Pop {population}")
        counter += 1
        gc.collect()
#--------Loop-----------#

print("Starting...")
main()