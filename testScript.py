import msprime
import core

from importlib import reload

reload(core)

# Times are provided in years, so we convert into generations.
generation_time = 29
T_HOM = 575e3 / generation_time
T_ARC = 415e3 / generation_time
T_INT = 50e3 / generation_time

# We need to work out the starting population sizes based on
# the growth rates provided for these two populations
N_HOM = 18500
N_ARC = 7100
N_ND = 3400
N_DN = 2500
N_HUM = 23700


speed_of_growth = 0.002

demography = msprime.Demography()

demography.add_population(
    name="HOM",
    description="Hominin",
    initial_size=N_HOM,
    growth_rate=speed_of_growth,
)
demography.add_population(
    name="ARC",
    description="Archaic",
    initial_size=N_ARC,
    growth_rate=speed_of_growth,
)
demography.add_population(
    name="ND",
    description="Neanderthal",
    initial_size=N_ND,
    growth_rate=speed_of_growth,
)
demography.add_population(
    name="DN",
    description="Denisovan",
    initial_size=N_DN,
    growth_rate=speed_of_growth,
)
demography.add_population(
    name="HUM",
    description="Modern human",
    initial_size=N_HUM,
    growth_rate=speed_of_growth,
)

demography.add_population(
    name="Af",
    description="Modern human which stayed in Africe",
    initial_size=N_HUM,
    growth_rate=speed_of_growth,
)

demography.add_population(
    name="EuAs",
    description="Modern human which left Africa for Eurasia",
    initial_size=N_HUM,
    growth_rate=speed_of_growth,
)

demography.add_population(
    name="EuAs2",
    description="Modern human admixed with Neanderthal",
    initial_size=N_HUM,
    growth_rate=speed_of_growth,
)

demography.add_admixture(time=T_INT, derived="EuAs2", ancestral=["EuAs", "ND"], proportions=[0.5, 0.5])
demography.add_population_split(time=3700, derived=["Af", "EuAs"], ancestral="HUM")
demography.add_population_split(time=T_ARC, derived=["ND", "DN"], ancestral="ARC")
demography.add_population_split(time=T_HOM, derived=["ARC", "HUM"], ancestral="HOM")

demography.sort_events()

# Recombination rate is set to the average value for humans. 
# record_migrations=true  is better when we want to do some ancestry analysis afterward, it offers more tools.
ts = msprime.sim_ancestry({"EuAs2": 1 , "Af": 100}, demography=demography, ploidy = 1,sequence_length=4000000,recombination_rate=2.5e-9,record_migrations=True, random_seed=5432) 

 # We add mutations, rate is set to the average value for humans, but in real life it actually varies depending on the position on the chromosome.
ts = msprime.sim_mutations(ts, rate=1.25e-8, random_seed=4321)    


tracts = core.clean_tracts(core.get_migrating_tracts(ts))

seq = core.createSeqObs(ts,10000,0,1,400,'Af')

print("Positions (given as intervals) in the genome inherited from a Neanderthal individual: \n")
print(tracts)
print("\n Number of mutations that are present in the individual from Eurasia, but absent from individuals from Africa. Genome is cut in segments of size 10000. \n")
print(seq)

