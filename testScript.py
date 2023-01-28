import msprime
import core
import HMMS
from importlib import reload

reload(core)


class Neanderthal:
    def __init__(self,s=5,f=5,l=5):
        self.s = s
        self.f = f
        self.l = l
        # Times are provided in years, so we convert into generations.
        self.generation_time = 29
        self.T_HOM = 575e3 / self.generation_time
        self.T_ARC = 415e3 / self.generation_time
        self.T_INT = 50e3 / self.generation_time

        # We need to work out the starting population sizes based on
        # the growth rates provided for these two populations
        self.N_HOM = 18500
        self.N_ARC = 7100
        self.N_ND = 3400
        self.N_DN = 2500
        self.N_HUM = 23700

        self.speed_of_growth = 0.002



    def start(self):
        demography = msprime.Demography()

        demography.add_population(
            name="HOM",
            description="Hominin",
            initial_size=self.N_HOM,
            growth_rate=self.speed_of_growth,
        )
        demography.add_population(
            name="ARC",
            description="Archaic",
            initial_size=self.N_ARC,
            growth_rate=self.speed_of_growth,
        )
        demography.add_population(
            name="ND",
            description="Neanderthal",
            initial_size=self.N_ND,
            growth_rate=self.speed_of_growth,
        )
        demography.add_population(
            name="DN",
            description="Denisovan",
            initial_size=self.N_DN,
            growth_rate=self.speed_of_growth,
        )
        demography.add_population(
            name="HUM",
            description="Modern human",
            initial_size=self.N_HUM,
            growth_rate=self.speed_of_growth,
        )

        demography.add_population(
            name="Af",
            description="Modern human which stayed in Africe",
            initial_size=self.N_HUM,
            growth_rate=self.speed_of_growth,
        )

        demography.add_population(
            name="EuAs",
            description="Modern human which left Africa for Eurasia",
            initial_size=self.N_HUM,
            growth_rate=self.speed_of_growth,
        )

        demography.add_population(
            name="EuAs2",
            description="Modern human admixed with Neanderthal",
            initial_size=self.N_HUM,
            growth_rate=self.speed_of_growth,
        )

        demography.add_admixture(time=self.T_INT, derived="EuAs2", ancestral=["EuAs", "ND"], proportions=[0.95, 0.05])
        demography.add_population_split(time=3700, derived=["Af", "EuAs"], ancestral="HUM")
        demography.add_population_split(time=self.T_ARC, derived=["ND", "DN"], ancestral="ARC")
        demography.add_population_split(time=self.T_HOM, derived=["ARC", "HUM"], ancestral="HOM")

        demography.sort_events()

        # Recombination rate is set to the average value for humans.
        # record_migrations=true  is better when we want to do some ancestry analysis afterward, it offers more tools.
        return demography

    def slow(self,X):
        ts = msprime.sim_ancestry({"EuAs2": 1, "Af": X}, demography=self.start(), ploidy=1,
                                  sequence_length=100000000, recombination_rate=2.5e-9, record_migrations=True,
                                  random_seed=5432)
        ts = msprime.sim_mutations(ts, rate=1.25e-8, random_seed=4321)
        # We add mutations, rate is set to the average value for humans, but in real life it actually varies depending on the position on the chromosome.
        tracts = core.clean_tracts(core.get_migrating_tracts(ts))
        cut = 10000
        seq = core.createSeqObs(ts, cut, 0, 1, 1000, 'Af')  # 1000
        final_array = [list(map(int, array)) for array in
                       np.ndarray.tolist(np.around(tracts / cut, decimals=0, out=None))]

        print("Positions (given as intervals) in the genome inherited from a Neanderthal individual: \n", \
              final_array, end="\n")

        states = [0, 1]
        S = HMMS.initS(0.95)
        A = HMMS.initA(self.T_INT, 2.5e-9, cut, 0.05)
        B = HMMS.initB(1.25e-8, cut, 3700, self.T_HOM)
        res = HMMS.viterbi(seq, S, A, B)
        tractsHMM = core.get_HMM_tracts(res)
        compar = self.comparaison(final_array, tractsHMM[1])
        Y = compar
        return Y
    def fast(self,Y):
        ts = msprime.sim_ancestry({"EuAs2": 1, "Af": 100}, demography=self.start(), ploidy=1,
                                  sequence_length=100000000, recombination_rate=2.5e-9, record_migrations=True,
                                  random_seed=5432)
        ts = msprime.sim_mutations(ts, rate=1.25e-8, random_seed=4321)
        # We add mutations, rate is set to the average value for humans, but in real life it actually varies depending on the position on the chromosome.
        tracts = core.clean_tracts(core.get_migrating_tracts(ts))
        cut = 10000
        seq = core.createSeqObs(ts, cut, 0, 1, Y, 'Af')  # 1000
        final_array = [list(map(int, array)) for array in
                       np.ndarray.tolist(np.around(tracts / cut, decimals=0, out=None))]

        print("Positions (given as intervals) in the genome inherited from a Neanderthal individual: \n", \
              final_array, end="\n")

        states = [0, 1]
        S = HMMS.initS(0.95)
        A = HMMS.initA(self.T_INT, 2.5e-9, cut, 0.05)
        B = HMMS.initB(1.25e-8, cut, 3700, self.T_HOM)
        res = HMMS.viterbi(seq, S, A, B)
        tractsHMM = core.get_HMM_tracts(res)
        compar = self.comparaison(final_array, tractsHMM[1])
        Y = compar
        return Y
    def length_dna(self,len):
        ts = msprime.sim_ancestry({"EuAs2": 1, "Af": 100}, demography=self.start(), ploidy=1,
                                  sequence_length=len, recombination_rate=2.5e-9, record_migrations=True,
                                  random_seed=5432)
        ts = msprime.sim_mutations(ts, rate=1.25e-8, random_seed=4321)
        # We add mutations, rate is set to the average value for humans, but in real life it actually varies depending on the position on the chromosome.
        tracts = core.clean_tracts(core.get_migrating_tracts(ts))
        cut = 10000
        seq = core.createSeqObs(ts, cut, 0, 1, 1000, 'Af')  # 1000
        final_array = [list(map(int, array)) for array in
                       np.ndarray.tolist(np.around(tracts / cut, decimals=0, out=None))]

        print("Positions (given as intervals) in the genome inherited from a Neanderthal individual: \n", \
              final_array, end="\n")

        states = [0, 1]
        S = HMMS.initS(0.95)
        A = HMMS.initA(self.T_INT, 2.5e-9, cut, 0.05)
        B = HMMS.initB(1.25e-8, cut, 3700, self.T_HOM)
        res = HMMS.viterbi(seq, S, A, B)
        tractsHMM = core.get_HMM_tracts(res)
        compar = self.comparaison(final_array, tractsHMM[1])
        Y = compar
        return Y

    def comparaison(self,reality, calculated):
        total_gen = sum(list(map(sum, reality)))
        errors = []

        precision_for_dectected = []

        not_inside = 0
        not_inside1 = 0
        for calc in calculated:
            flag = True
            for real in reality:
                start_c, end_c, start_r, end_r = calc[0], calc[1], real[0], real[1]
                maximum_margin = int(end_c - start_c) * 0.75
                if abs(start_c - start_r) < maximum_margin and abs(end_c - end_r) < maximum_margin:
                    error = abs(start_c - start_r) + abs(end_c - end_r)
                    errors.append(error)

                    precision_inside = 1 - error / (end_r - start_r)
                    precision_for_dectected.append((end_r - start_r, precision_inside))

                    flag = False
                    break
                error = abs(start_c - end_c)
                errors.append(error)
            if flag:
                not_inside += 1
            not_inside1 += 1
        rate_I_NI = not_inside / not_inside1

        rate_more = (len(calculated) - (not_inside1 - not_inside)) / len(reality)

        error_rate = (sum(errors) / total_gen) * 100
        return error_rate, rate_I_NI, rate_more, precision_for_dectected

    def graph(self,X, Y):
        x, y = np.array(X), np.array(Y)
        X_Y_Spline = make_interp_spline(x, y)

        X_ = np.linspace(x.min(), x.max(), 500)
        Y_ = X_Y_Spline(X_)

        plt.plot(X_, Y_)
        plt.title("Relation between the number of African genomes with the incorrect intervals ratio")
        plt.xlabel("Number of Afrian genomes")
        plt.ylabel("Incorrect intervals ratio")
        plt.show()
