import starsim as ss
import stisim as sti

class YouthPrev(ss.Analyzer):
    """Track gonorrhea prevalence among 15-24 year olds."""

    def init_results(self):
        super().init_results()
        # define_results registers a result with the sim's results system.
        # After sim.run(), it appears in sim.results.youthprev.ng_prev_15_24
        # and in sim.to_df() as 'youthprev.ng_prev_15_24'.
        self.define_results(
            ss.Result('ng_prev_15_24', dtype=float, scale=False, label='NG prev (15-24)'),
        )

    def step(self):
        ppl = self.sim.people
        youth = (ppl.age >= 15) & (ppl.age < 25)       # Boolean mask for 15-24 year olds
        n_youth = youth.count()                          # How many are in this group
        if n_youth > 0:
            infected = self.sim.diseases.ng.infected     # Boolean: who is infected
            prev = (infected & youth).count() / n_youth  # Prevalence = infected / total
            self.results['ng_prev_15_24'][self.ti] = prev


youth_prev = YouthPrev()
ng = sti.Gonorrhea(init_prev=ss.bernoulli(p=0.1))  # Higher init_prev so the demo has visible infections among youth
sim = sti.Sim(diseases=ng, analyzers=youth_prev, n_agents=2000, start=2010, stop=2030)
sim.run(verbose=0)

print(sim.results.youthprev.ng_prev_15_24)
sim.results.youthprev.ng_prev_15_24.plot()


#How it works:

#init_results is called during sim.init(). Call self.define_results(...) to register your results. Each ss.Result gets a pre-allocated array matching the simulation’s number of timesteps. The result is stored under sim.results.<module_name>.<result_name> – in this case, sim.results.youthprev.ng_prev_15_24.

#step is called every timestep. Access the sim via self.sim, compute your quantity, and store it in self.results[name][self.ti] where self.ti is the current timestep index.

#After sim.run(), the result is available everywhere: sim.results, sim.to_df(), sim.plot(), and as a calibration target (see Calibration tutorial).