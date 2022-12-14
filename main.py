import msprime

def main():
    # Times are provided in years, so we convert into generations.
    generation_time = 29
    T_HOM = 575e3 / generation_time
    T_ARC = 415e3 / generation_time

    # We need to work out the starting population sizes based on
    # the growth rates provided for these two populations
    N_HOM = 18500
    N_ARC = 7100
    N_ND = 3400
    N_DN = 2500
    N_HUM = 23700

    speed_of_growth = 0.75

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

    demography.add_population_split(time=T_HOM, derived=["ARC", "HUM"], ancestral="HOM")
    demography.add_population_split(time=T_ARC, derived=["ND", "DN"], ancestral="ARC")
    demography.sort_events()

    # Recombination rate is set to the average value for humans. 
    # record_migrations=true  is better when we want to do some ancestry analysis afterward, it offers more tools.
    ts = msprime.sim_ancestry({"HUM": 4}, demography=demography, ploidy = 1, sequence_length=10000000,recombination_rate=2.5e-9,record_migrations=True, random_seed=543212) 
    
     # We add mutations, rate is set to the average value for humans, but in real life it actually varies depending on the position on the chromosome.
    ts = msprime.sim_mutations(ts, rate=1.25e-8, random_seed=4321)    

    
if __name__ == "__main__":
    main()
