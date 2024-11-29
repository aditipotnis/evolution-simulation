import numpy as np
import random
import math
import matplotlib.pyplot as plt

def read_vcf(file_name):
    initial_population = []
    
    with open(file_name, 'r') as vcf_file:
        for line in vcf_file:
            if line.startswith('#'):
                continue  
            
            cols = line.strip().split('\t')
            individuals = cols[9:]  # VCF genotypes usually start from column 9 onwards
            population = []
            for ind in individuals:
                maternal, paternal = map(int, ind.split('|'))
                population.append(([maternal], [paternal]))  
            initial_population.append(population)
    
    # Transpose the population so each entry is a diploid individual
    transposed_population = [[list(map(int, [ind[0][0] for ind in individual])),
                              list(map(int, [ind[1][0] for ind in individual]))] 
                              for individual in zip(*initial_population)]
    return transposed_population

def select_parent_neutral(population):
    return random.choice(population)

def fitness(individual):
    return 1 # Fitness neutral for all

def reproduce(parent1, parent2, L):
    # Simulate recombination and create child
    cross_pos_maternal = random.randint(0, L-1)
    if random.random() < 0.5:
        child_maternal = parent1[0][:cross_pos_maternal] + parent1[1][cross_pos_maternal:]
    else:
        child_maternal = parent1[1][:cross_pos_maternal] + parent1[0][cross_pos_maternal:]

    cross_pos_paternal = random.randint(0, L-1)
    if random.random() < 0.5:
        child_paternal = parent2[0][:cross_pos_paternal] + parent2[1][cross_pos_paternal:]
    else:
        child_paternal = parent2[1][:cross_pos_paternal] + parent2[0][cross_pos_paternal:]

    return (child_maternal, child_paternal)

def simulate_evolution(initial_population, generations, fitness_function):
    population = initial_population
    N = len(population)
    
    for generation in range(generations):
        fitness_values = [fitness_function(individual) for individual in population]
        new_population = []
    
        for _ in range(N):
            parent1 = select_parent_neutral(population)
            parent2 = select_parent_neutral(population)
            while parent2 == parent1:
                parent2 = select_parent_neutral(population)
            
            child = reproduce(parent1, parent2, len(parent1[0]))
            new_population.append(child)
    
        population = new_population
    
    return population

def calculate_extinct_alleles(initial_population, final_population):
    L = len(initial_population[0][0])
    N = len(final_population)
    extinct_count = 0

    for i in range(L):
        # Check if the alternate allele was present in the initial population
        alt_present_initially = False
        for individual in initial_population:
            alleles = [
                individual[0][i],  # maternal chromosome
                individual[1][i]   # paternal chromosome
            ]
            if 1 in alleles:
                alt_present_initially = True
                break

        if not alt_present_initially:
            continue  # Skip SNPs where the alternate allele was not present

        # Check if the alternate allele is extinct in the final population
        alt_present_finally = False
        for individual in final_population:
            alleles = [
                individual[0][i],  # maternal chromosome
                individual[1][i]   # paternal chromosome
            ]
            if 1 in alleles:
                alt_present_finally = True
                break

        if not alt_present_finally:
            extinct_count += 1

    return extinct_count

def allele_frequencies(population, L):
    N = len(population)
    freqs = [0] * L
    
    for i in range(L):
        allele_count = 0
        for individual in population:
            alleles = [
                individual[0][i],  # maternal chromosome
                individual[1][i]   # paternal chromosome
            ]
            allele_count += sum(alleles)  # Count the number of alternate alleles (1's)
        freqs[i] = allele_count / (2 * N)  # Divide by total number of alleles (2 per individual)
    
    return freqs

def calculate_extinction_and_fixation(population, L):
    N = len(population)
    extinct_count = 0
    fixation_count = 0

    for i in range(L):
        allele_count = 0
        for individual in population:
            alleles = [
                individual[0][i],  # maternal chromosome
                individual[1][i]   # paternal chromosome
            ]
            allele_count += sum(alleles)  # Count the number of alternate alleles (1's)

        if allele_count == 0:
            extinct_count += 1  # All alleles are reference (0)
        elif allele_count == 2 * N:
            fixation_count += 1  # All alleles are alternate (1)

    extinction_prob = extinct_count / L
    fixation_prob = fixation_count / L
    return extinction_prob, fixation_prob



'-----------------------------------------------------------------------------------------------------------------------------------------------------'
'''
# Q1 PART C

# Parameters for the simulation
N = 10000  # Population size
L = 1  # Number of SNPs
generations = 1  # Simulate for only 1 generation

# Read the initial population from VCF file
initial_population = read_vcf('/Users/aditipotnis/Downloads/initial_population.vcf')
final_population = simulate_evolution(initial_population, generations, fitness)

# Run the simulation
extinction_count = calculate_extinct_alleles(initial_population, final_population)

# Calculate the extinction probability
total_snp_count = L * N  # Each SNP can be in one of N individuals
extinction_probability = extinction_count / total_snp_count

# Output the extinction probability and comparison to theoretical value
theoretical_probability = 1 / math.e
print(f"Estimated extinction probability: {extinction_probability}")
print(f"Theoretical extinction probability: {theoretical_probability}")
print(f"Extinction count: {extinction_count}, Total SNPs: {total_snp_count}")
'''

'----------------------------------------------------------------------------------------------------------------------------------------------------'
'''
# Q1 PART D

# Parameters for the simulation
N = 10000  # Population size
L = 100  # Number of SNPs
generations = 20  # Simulate for only 1 generation

# Read the initial population from VCF file
initial_population = read_vcf('/Users/aditipotnis/Downloads/initial_population.vcf')

# Track allele frequencies over generations for the first 100 SNPs
allele_frequency_over_time = []

# Store initial population allele frequencies
initial_frequencies = allele_frequencies(initial_population, L)
allele_frequency_over_time.append(initial_frequencies)

# Simulate for 20 generations
for gen in range(generations):
    population = simulate_evolution(initial_population, 1, fitness)
    frequencies = allele_frequencies(population, L)
    allele_frequency_over_time.append(frequencies)
    initial_population = population  # Update the population for the next generation

# Convert list of frequencies to a numpy array for plotting
allele_frequency_over_time = np.array(allele_frequency_over_time)

# Plot the allele frequencies for the first 100 SNPs over 20 generations
plt.figure(figsize=(10, 6))
for snp_index in range(L):
    plt.plot(range(generations + 1), allele_frequency_over_time[:, snp_index], label=f"SNP {snp_index+1}", alpha=0.5)

plt.title('Alternate Allele Frequency Over 20 Generations for First 100 SNPs')
plt.xlabel('Generation')
plt.ylabel('Alternate Allele Frequency')
plt.grid(True)
plt.show()
'''
'----------------------------------------------------------------------------------------------------------------------------------------------------'
'''
# Q1 PART E

# Parameters for the simulation
N = 10000  # Population size
L = 10000  # Total number of SNPs
generations = 1000  # Simulate for 1000 generations

# Read the initial population from VCF file
initial_population = read_vcf('/Users/aditipotnis/Downloads/initial_population.vcf')

# Track extinction and fixation probabilities over generations
extinction_probs = []
fixation_probs = []

# Simulate for 1000 generations and calculate extinction and fixation probabilities
for gen in range(generations + 1):
    extinction_prob, fixation_prob = calculate_extinction_and_fixation(initial_population, L)
    extinction_probs.append(extinction_prob)
    fixation_probs.append(fixation_prob)

    if gen < generations:  # Continue simulating for the next generation
        population = simulate_evolution(initial_population, 1, fitness)
        initial_population = population  # Update the population for the next generation

# Plot the extinction and fixation probabilities over generations
plt.figure(figsize=(10, 6))
plt.plot(range(generations + 1), extinction_probs, label='Extinction Probability')
plt.plot(range(generations + 1), fixation_probs, label='Fixation Probability')
plt.title('Allele Extinction and Fixation Probability Over 1000 Generations')
plt.xlabel('Generation')
plt.ylabel('Probability')
plt.legend()
plt.grid(True)
plt.show()
'''
'----------------------------------------------------------------------------------------------------------------------------------------------------'
'''
# Q1 PART F

# Function to calculate fitness based on the genotype at SNP42
def fitness(individual, snp_index=42):
    maternal = individual[0][snp_index]
    paternal = individual[1][snp_index]
    
    # Heterozygous (0|1 or 1|0) at SNP42 -> fitness = 1.5
    if (maternal == 1 and paternal == 0) or (maternal == 0 and paternal == 1):
        return 1.5
    
    # Homozygous for alternate allele (1|1) at SNP42 -> fitness = 2
    if maternal == 1 and paternal == 1:
        return 2
    
    # All other genotypes -> fitness = 1
    return 1

def select_parent_weighted(population, fitness_values):
    # Select a parent with probability proportional to its fitness
    total_fitness = sum(fitness_values)
    selection_probs = [f / total_fitness for f in fitness_values]
    return random.choices(population, weights=selection_probs, k=1)[0]

# Simulate evolution with fitness function incorporating SNP42 special rules
def simulate_evolution_with_fitness(initial_population, generations, snp_index=42):
    population = initial_population
    N = len(population)
    
    for generation in range(generations):
        # Calculate fitness for each individual based on SNP42
        fitness_values = [fitness(individual, snp_index) for individual in population]
        new_population = []
    
        for _ in range(N):
            # Select parents with weighted probabilities based on fitness
            parent1 = select_parent_weighted(population, fitness_values)
            parent2 = select_parent_weighted(population, fitness_values)
            while parent2 == parent1:
                parent2 = select_parent_weighted(population, fitness_values)
            
            child = reproduce(parent1, parent2, len(parent1[0]))
            new_population.append(child)
    
        population = new_population
    
    return population

# Function to calculate the extinction of alternate allele at SNP42
def calculate_extinction_at_snp42(population, snp_index=42):
    N = len(population)
    allele_count = 0
    for individual in population:
        alleles = [
            individual[0][snp_index],  # maternal chromosome
            individual[1][snp_index]   # paternal chromosome
        ]
        allele_count += sum(alleles)  # Count the number of alternate alleles (1's)
    
    # If all alleles are reference alleles (allele_count == 0), then it's extinct
    if allele_count == 0:
        return True
    return False

# Run the simulation repeatedly and calculate extinction probability for SNP42
def run_extinction_simulation(initial_population, generations, repetitions, snp_index=42):
    extinction_count = 0
    
    for _ in range(repetitions):
        final_population = simulate_evolution_with_fitness(initial_population, generations, snp_index)
        if calculate_extinction_at_snp42(final_population, snp_index):
            extinction_count += 1
    
    extinction_probability = extinction_count / repetitions
    return extinction_probability

# Parameters
N = 100  # Population size
L = 100  # Total number of SNPs
generations = 1  # Simulate for only 1 generation
snp_index = 42  # The SNP with special fitness rules
repetitions = 1000  # Number of repeated simulations to estimate extinction probability

# Read the initial population from VCF file
initial_population = read_vcf('/Users/aditipotnis/Downloads/initial_population.vcf')

# Run the simulation multiple times and calculate extinction probability
extinction_probability = run_extinction_simulation(initial_population, generations, repetitions, snp_index)

# Output the estimated extinction probability
print(f"Estimated extinction probability at SNP42 after 1 generation: {extinction_probability}")

# Output : Estimated extinction probability at SNP42 after 1 generation: 0.214 
'''

'----------------------------------------------------------------------------------------------------------------------------------------------------'
'''
# Q1 PART G

# Fitness function based on genotype at SNP42
def fitness_snp42(individual):
    maternal, paternal = individual
    
    # Heterozygous (0|1 or 1|0) at SNP42 -> fitness = 1.5
    if (maternal == 1 and paternal == 0) or (maternal == 0 and paternal == 1):
        return 1.5
    
    # Homozygous for alternate allele (1|1) at SNP42 -> fitness = 2
    if maternal == 1 and paternal == 1:
        return 2
    
    # All other genotypes -> fitness = 1
    return 1

# Select a parent with probability proportional to its fitness
def select_parent_weighted(population, fitness_values):
    total_fitness = sum(fitness_values)
    selection_probs = [f / total_fitness for f in fitness_values]
    return random.choices(population, weights=selection_probs, k=1)[0]

# Reproduce a child based on parents' alleles at SNP42
def reproduce_snp42(parent1, parent2):
    # Randomly select one allele from each parent
    child_maternal = random.choice(parent1)
    child_paternal = random.choice(parent2)
    return (child_maternal, child_paternal)

# Simulate evolution for SNP42 only
def simulate_evolution_snp42(initial_population, generations):
    population = initial_population
    N = len(population)
    
    for generation in range(generations):
        # Calculate fitness for each individual based on SNP42
        fitness_values = [fitness_snp42(individual) for individual in population]
        new_population = []
    
        for _ in range(N):
            # Select parents with weighted probabilities based on fitness
            parent1 = select_parent_weighted(population, fitness_values)
            parent2 = select_parent_weighted(population, fitness_values)
            child = reproduce_snp42(parent1, parent2)
            new_population.append(child)
    
        population = new_population
    
    return population

# Function to check if alternate allele (1) is extinct at SNP42
def calculate_extinction_snp42(population):
    for individual in population:
        if 1 in individual:
            return False
    return True

# Function to check if alternate allele (1) has fixed (everyone has 1|1 genotype)
def calculate_fixation_snp42(population):
    for individual in population:
        if individual != (1, 1):
            return False
    return True

# Run the simulation multiple times to estimate extinction and fixation probabilities
def run_extinction_and_fixation_simulation_snp42(initial_population, generations, repetitions):
    extinction_count = 0
    fixation_count = 0
    
    for _ in range(repetitions):
        final_population = simulate_evolution_snp42(initial_population, generations)
        
        # Check for extinction
        if calculate_extinction_snp42(final_population):
            extinction_count += 1
        # Check for fixation
        elif calculate_fixation_snp42(final_population):
            fixation_count += 1
    
    extinction_probability = extinction_count / repetitions
    fixation_probability = fixation_count / repetitions
    return extinction_probability, fixation_probability


# Parameters for part (g)
generations = 100  # Simulate for 100 generations
repetitions = 1000  # Number of repeated simulations to estimate extinction and fixation probabilities

# Initial population for SNP42 (you can load from a file if needed)
# Example population structure: list of (maternal, paternal) tuples
# This will need to be initialized properly from your VCF file
initial_population = [(0, 0)] * 95 + [(1, 0)] * 5  # Example population with 5 individuals having the alternate allele

# Run the simulation multiple times and calculate extinction and fixation probabilities
extinction_probability, fixation_probability = run_extinction_and_fixation_simulation_snp42(initial_population, generations, repetitions)

# Output the estimated probabilities
print(f"Estimated extinction probability at SNP42 after {generations} generations: {extinction_probability}")
print(f"Estimated fixation probability at SNP42 after {generations} generations: {fixation_probability}")
'''
'----------------------------------------------------------------------------------------------------------------------------------------------------'
'''
# Q1 PART H

# Deleterious fitness function based on genotype at SNP42
def fitness_deleterious_snp42(individual):
    maternal, paternal = individual
    
    # Heterozygous (0|1 or 1|0) at SNP42 -> fitness = 0.9
    if (maternal == 1 and paternal == 0) or (maternal == 0 and paternal == 1):
        return 0.9
    
    # Homozygous for alternate allele (1|1) at SNP42 -> fitness = 0.8
    if maternal == 1 and paternal == 1:
        return 0.8
    
    # All other genotypes -> fitness = 1 (homozygous reference)
    return 1

# Select a parent with probability proportional to its fitness
def select_parent_weighted(population, fitness_values):
    total_fitness = sum(fitness_values)
    selection_probs = [f / total_fitness for f in fitness_values]
    return random.choices(population, weights=selection_probs, k=1)[0]

# Reproduce a child based on parents' alleles at SNP42
def reproduce_snp42(parent1, parent2):
    # Randomly select one allele from each parent
    child_maternal = random.choice(parent1)
    child_paternal = random.choice(parent2)
    return (child_maternal, child_paternal)

# Simulate evolution for SNP42 only with deleterious fitness
def simulate_evolution_snp42_deleterious(initial_population, generations):
    population = initial_population
    N = len(population)
    
    for generation in range(generations):
        # Calculate fitness for each individual based on SNP42 deleterious effect
        fitness_values = [fitness_deleterious_snp42(individual) for individual in population]
        new_population = []
    
        for _ in range(N):
            # Select parents with weighted probabilities based on fitness
            parent1 = select_parent_weighted(population, fitness_values)
            parent2 = select_parent_weighted(population, fitness_values)
            child = reproduce_snp42(parent1, parent2)
            new_population.append(child)
    
        population = new_population
    
    return population

# Function to check if alternate allele (1) is extinct at SNP42
def calculate_extinction_snp42(population):
    for individual in population:
        if 1 in individual:
            return False
    return True

# Function to check if alternate allele (1) has fixed (everyone has 1|1 genotype)
def calculate_fixation_snp42(population):
    for individual in population:
        if individual != (1, 1):
            return False
    return True

# Run the simulation multiple times and calculate extinction and fixation probabilities for deleterious allele
def run_extinction_and_fixation_simulation_snp42_deleterious(initial_population, generations, repetitions):
    extinction_count = 0
    fixation_count = 0
    
    for _ in range(repetitions):
        final_population = simulate_evolution_snp42_deleterious(initial_population, generations)
        
        # Check for extinction
        if calculate_extinction_snp42(final_population):
            extinction_count += 1
        # Check for fixation
        elif calculate_fixation_snp42(final_population):
            fixation_count += 1
    
    extinction_probability = extinction_count / repetitions
    fixation_probability = fixation_count / repetitions
    return extinction_probability, fixation_probability

# Parameters for deleterious scenario (part g with deleterious fitness)
generations = 100  # Simulate for 100 generations
repetitions = 1000  # Number of repeated simulations to estimate extinction and fixation probabilities

# Initial population for SNP42 (you can load from a file if needed)
# Example population structure: list of (maternal, paternal) tuples
# This will need to be initialized properly from your VCF file
initial_population = [(0, 0)] * 95 + [(1, 0)] * 5  # Example population with 5 individuals having the alternate allele

# Run the simulation for deleterious allele at SNP42
extinction_probability_deleterious, fixation_probability_deleterious = run_extinction_and_fixation_simulation_snp42_deleterious(
    initial_population, generations, repetitions)

# Output the estimated probabilities for the deleterious scenario
print(f"Estimated extinction probability at SNP42 (deleterious) after {generations} generations: {extinction_probability_deleterious}")
print(f"Estimated fixation probability at SNP42 (deleterious) after {generations} generations: {fixation_probability_deleterious}")
'''