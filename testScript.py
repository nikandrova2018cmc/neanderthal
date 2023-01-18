import msprime
import core
import HMMS
from comparaison import *
from importlib import reload



reload(core)

def start():

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

    demography.add_admixture(time=T_INT, derived="EuAs2", ancestral=["EuAs", "ND"], proportions=[0.95, 0.05])
    demography.add_population_split(time=3700, derived=["Af", "EuAs"], ancestral="HUM")
    demography.add_population_split(time=T_ARC, derived=["ND", "DN"], ancestral="ARC")
    demography.add_population_split(time=T_HOM, derived=["ARC", "HUM"], ancestral="HOM")

    demography.sort_events()

    # Recombination rate is set to the average value for humans.
    # record_migrations=true  is better when we want to do some ancestry analysis afterward, it offers more tools.
    def test(X):
        ts = msprime.sim_ancestry({"EuAs2": 1 , "Af": X}, demography=demography, ploidy = 1,sequence_length=100000000,recombination_rate=2.5e-9,record_migrations=True, random_seed=5432)
        return ts
    x = []
    y = []
    X = 1

    while X <= 400:

        ts = test(X)

         # We add mutations, rate is set to the average value for humans, but in real life it actually varies depending on the position on the chromosome.
        ts = msprime.sim_mutations(ts, rate=1.25e-8, random_seed=4321)


        tracts = core.clean_tracts(core.get_migrating_tracts(ts))


        cut = 10000
        seq = core.createSeqObs(ts,cut,0,1,1000,'Af')

        final_array = [list(map(int,array)) for array in np.ndarray.tolist(np.around(tracts/cut, decimals=0, out=None))]


        print("Positions (given as intervals) in the genome inherited from a Neanderthal individual: \n",\
              final_array,end="\n")
        #print("\n Number of mutations that are present in the individual from Eurasia, but absent from individuals from Africa. Genome is cut in segments of size 10000. \n")
        #print(seq)


        states = [0,1]
        S = HMMS.initS(0.95)
        A = HMMS.initA(T_INT,2.5e-9,cut,0.05)
        B = HMMS.initB(1.25e-8,cut,3700,T_HOM)
        res = HMMS.viterbi( seq, S, A, B)
        tractsHMM = core.get_HMM_tracts(res)
        try:
            compar = comparaison(final_array,tractsHMM[1])
            Y = compar[1]
        except:
            break

        print("Calculated data",tractsHMM[1],sep="\n")
        print("total genomes from Nean:", compar[0], "%", sep=" ")
        x.append(X)
        y.append(Y)

        X += 10

    print(x,y,sep="  ")

    graph(x,y)

