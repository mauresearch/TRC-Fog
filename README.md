# Energy-Aware Workflow Scheduling in Fog Computing with Coupled Thermal and Reliability Modeling, Joint Replication and Predictive SLA Control


MATLAB implementation of the TRC framework for scheduling scientific workflows on fog nodes whose speed and reliability depend on temperature. The code plans the placement and standby replication of workflow tasks with the improved multi-objective walrus optimizer (IMOWO), evaluates every schedule with Monte-Carlo ground-truth simulators, and includes the predictive SLA controller (PSC) used at run time.

The repository accompanies the paper

> Energy-Aware Workflow Scheduling in Fog Computing with Coupled Thermal and Reliability Modeling, Joint Replication and Predictive SLA Control (under review).

## Requirements

- MATLAB or GNU Octave 8.4 or later
- No additional toolboxes

## Repository structure

```
framework/                  scheduling framework and experiments
  Main_TRC_Experiments.m    runs one experiment and stores the raw results
  Summarize_TRC.m           builds the CSV tables, paired tests and figures
  experiment_config.m       design of every experiment
  trc_models.m              planning model of every compared method
  trc_evaluate_schedule.m   thermal-reliability coupled planning model
  trc_operators.m           thermal repair, FTRS and DRFA operators
  trc_runtime_plan.m        plan information for the PSC controller
  simulate_TRC_truth.m      ground truth 1 (steady-state thermal model)
  simulate_TRC_truth2.m     ground truth 2 (transient thermal model, step-wise DTM, Weibull failures)
  ARSENAL_forecast.m        load forecast of the fog nodes
  CARS_forecast.m           criticality ranking of the tasks
  psb_schedule.m            PSB baseline
  tdrm_schedule.m           TDRM baseline
  eft_schedule.m            EFT baseline
  generate_workflow_DAG.m   Montage, Epigenomics, CyberShake, Inspiral, SIPHT
  generate_fog_nodes.m      heterogeneous fog environment
optimizers/
  IMOWO.m                   improved multi-objective walrus optimizer
  MOWO.m, NSGA2.m, SPEA2.m, MO_GWO.m, MO_PSO.m
  MOWO_published/           original MOWO code used by MOWO.m
```

## Quick start

1. Open `framework/Main_TRC_Experiments.m`.
2. Set `Experiment` to one of the names in the table below. Set `QuickTest = true` for a short check run.
3. Run the script. The raw results are written to `framework/Results_TRC/<Experiment>/<Experiment>_<RunTag>.mat`.
4. Open `framework/Summarize_TRC.m`, set the same `Experiment` and run it. It writes `<Experiment>_AllRuns.csv`, `<Experiment>_Summary.csv` and `<Experiment>_Stats.csv` and plots the main metrics.

`HV_RefPoint_Sensitivity.m` recomputes the hypervolume of the `E2` fronts for several reference points.

An experiment can be split over several MATLAB sessions by giving each session a different `RunTag` and a subset of `WorkloadsToRun` or `MethodsToRun`. The random seeds depend only on the workflow, scenario and run, so the split does not change the results. `Summarize_TRC.m` merges all files of an experiment.

## Experiments

| Name | Content | Methods |
|------|---------|---------|
| `E2` | Optimizers on the TRC problem | IMOWO, MOWO, NSGA-II, MO-GWO, MO-PSO, SPEA2 |
| `E3` | Comparison with existing schedulers, ground truth 1 | TRC, TRC-EqE, EFT, Thermal-only, Failure-only, PSB, TDRM |
| `E3T` | Same comparison, ground truth 2 | as `E3` |
| `E3W` | Ground truth 2 with Weibull shapes 1.0 and 3.0 | as `E3` |
| `E4` | Ablation of the planning model | TRC, TRC-noThrottle, TRC-noArrhenius, TRC-noReplica, TRC-noSleep, TRC-noForecast, Sequential |
| `E5` | Ablation of the optimizer | IMOWO, IMOWO-noOps, -noLHS, -noOBL, -noCDLeader, -noKNN, -noLevy, MOWO-core |
| `E6` | Sensitivity analysis (p0, Ea, ambient, kappa, deadline, scale, reservation cost) | TRC, EFT, Thermal-only, Failure-only |
| `E6p` | Low base failure probabilities (p0 = 0.0001, 0.001) | as `E6` |
| `E8` | Sensitivity to the CARS weights (equal weights and one dominant indicator) | TRC |
| `E7` | Run-time SLA control | TRC, TRC+Local, TRC+PSC-U, TRC+PSC, EFT, EFT+PSC |
| `E7b` | Extensions of the controller (PSC-R, temperature update after a move, @T) | TRC, TRC+PSC, TRC+PSC-R, TRC+PSC@T, TRC+PSC-R@T, TRC+Local@T, EFT, EFT+PSC@T |

With the default settings (population 100, 100 iterations, 200 Monte-Carlo samples), `E3` takes several hours on a desktop computer.

## Main parameters

All parameters are set at the top of `Main_TRC_Experiments.m`. Their values correspond to Table 6 of the paper.

| Parameter | Value | Meaning |
|-----------|-------|---------|
| `P.p0`, `P.Tref`, `P.Ea` | 0.02, 60 °C, 0.7 eV | base failure probability, reference temperature, activation energy |
| `P.throttleSlope`, `P.minSpeed` | 0.03, 0.4 | speed loss per °C above the limit, lowest speed |
| `P.thermalIters` | 3 | fixed-point iterations between schedule and temperature (Q) |
| `P.sleepFrac`, `P.standbyFrac`, `P.replicaCostFrac` | 0.10, 0.05, 0.50 | sleep power, standby power, reservation cost of a replica |
| `P.epsMin`, `P.epsMax`, `P.pscGain`, `P.pscRiskTol` | 0.05, 0.30, 0.05, 1.0 | PSC tolerances, minimum gain of a move, risk tolerance of PSC-R |
| `P.numMC` | 200 | Monte-Carlo samples per schedule |
| `P.etaError`, `P.EaError`, `P.bgNoise` | 0.10, 0.15, 0.05 | model error of the simulator |
| `P.detectDelay`, `P.replicaDelay`, `P.wakeDelay` | 2 s, 1 s, 10 s | recovery delays |
| `DeadlineFactor` | 1.25 | deadline relative to the EFT makespan |
| `PopSize`, `MaxIter` | 100, 100 | population size and iterations of all optimizers |

## Compared methods

- **TRC**: the proposed framework, compromise schedule of the IMOWO front.
- **TRC-EqE**: the schedule of the same front with the lowest expected critical loss whose planned energy does not exceed that of PSB.
- **EFT**: earliest-finish-time list scheduling.
- **Thermal-only** and **Failure-only**: IMOWO with a thermal-aware model without replicas, and with a temperature-independent failure model with replicas.
- **PSB**: primary–secondary backup scheduling, implemented from S. She et al., Future Generation Computer Systems 168 (2025) 107738.
- **TDRM**: temperature- and deadline-aware resource management, implemented from M. Niknafs, P. Eles, Z. Peng, Journal of Systems Architecture 167 (2025) 103506.

## Third-party code

The folder `optimizers/MOWO_published` contains the published code of the multi-objective walrus optimizer and routines of the walrus optimizer (M. Han et al., Expert Systems with Applications, 2023) and of NSGA-II, with their original headers. `non_domination_sort_mod.m` is a vectorised version of the original sorting routine that returns identical results.

## Citation

If you use this code, please cite the paper above. The full reference will be added after publication.
