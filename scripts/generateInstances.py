import argparse
import json
import pathlib
import warnings

def main(args):

    """
    :param args: all arguments given to the python script
    :return: Nothing, generate all configuration file for solve or generate instances
    """

    # Initialize parameters from command-line arguments

    """
    GENERATION INSTANCE
    """
    # characteristics of instance
    noIdenticalJobs = False if args.no_identical_jobs == 0 else True
    TF = args.tf
    RDD = args.rdd
    NList = args.n_list
    fracOfN = args.frac_of_n
    M = [tuple(map(int, pair.split(','))) for pair in args.m]
    Vmax = args.vmax
    V0 = args.v0
    # seed and nb instance to generate
    nbInstanceToGenerate = args.nb_instance_to_generate
    seed = args.seed
    # paths
    pathSaveInstance = pathlib.Path(str(pathlib.Path(__file__).parent.resolve()) + "/.."+ args.path_save_instance ).resolve()
    configFileNameGenerate = args.config_file_name_generate

    """
    GENERAL SOLVER PARAM
    """
    verbose = args.verbose
    outputResult = pathlib.Path(str(pathlib.Path(__file__).parent.resolve()) + "/.."+ args.output_result).resolve()
    configFileNameSolve = args.config_file_name_solve
    methods = args.methods
    timeLimit = int(args.timeLimit)

    """
    BRANCH AND BOUND PARAM
    """
    memo = False if args.memorization == 0 else True
    useHeuristic = args.use_heuristic_gen_col
    maxNbCallHeuristic = args.maxNbCallHeuristic
    nbMinStateDP = args.nbMinStateDP
    strategy = args.strategy

    """
    BEAM SEARCH PARAM
    """
    beamSize = args.beamSize
    recovering = args.recovering

    """
    Generate instances config 
    """
    # list of names for each instance
    listNameInstances = []

    # JSON
    configGenerate = {}

    generate = {
        "seed": seed
    }

    # the list of configuration of each instance
    listConfigInstance = []
    # define the config files for generate instance
    for (mMax, M0) in M:
        for N in NList:
            for n in [int(N * x) for x in fracOfN]:
                for tf in TF:
                    for rdd in RDD:
                        # Keep the list of instances for make the solve config file
                        instancePath = "_n_" + str(n) + "_N_" + str(N) + "_tf_" + str(tf) + "_rdd_" + str(
                            rdd) + "_mMax_" + str(mMax) + "_m0_" + str(M0) + "_Vmax_" + str(Vmax) + "_V0_" + str(V0)
                        listNameInstances.append(instancePath)

                        instance = {
                            "basePath": str(pathSaveInstance)+ "/",
                            "numberInstance": nbInstanceToGenerate}
                        param = {
                                "noIdenticalJobs": noIdenticalJobs ,
                                "n": n,
                                "N": N,
                                "M_max": mMax,
                                "M_0": M0,
                                "V_max": Vmax,
                                "V_0": V0,
                                "pi": {
                                    "inf": 1,
                                    "sup": 100
                                },
                                "di": {
                                    "TF": tf,
                                    "RDD": rdd
                                },
                                "wi": {
                                    "inf": 1,
                                    "sup": 10
                                }
                                }
                        instance["paramInstance"] = param

                        listConfigInstance.append(instance)

    generate["instances"] = listConfigInstance
    configGenerate["generate"] = generate
    if configFileNameGenerate is not None:
        configPath = pathlib.Path(str(pathlib.Path(str(pathlib.Path(__file__).parent.resolve()))) + f"/../instances/configs/{configFileNameGenerate}.json").resolve()
        configPath.parent.mkdir(parents=True, exist_ok=True)
        with open(configPath, 'w') as file:
            file.write(json.dumps(configGenerate))
        print(f"Generated config files:\n- {configPath}")

    """
    Generate solver config file 
    """
    if configFileNameSolve is not None:
        # Add methods based on user selection
        for method in methods:
            for x in listNameInstances:
                for i in range(nbInstanceToGenerate):
                    # Solve instances config
                    configSolve = {"solve": {"verbose": verbose, "output": str(outputResult), "methods": []}}
                    if method == "MIP":
                        configSolve["solve"]["methods"].append({
                            "name": "MIP",
                            "verbose": verbose,
                            "timeLimits": timeLimit,
                            "relaxation": False,
                            "instances": [
                                {"path": f"{pathSaveInstance}/instance{i}{x}.txt"}
                            ]
                        })
                    elif method == "LP":
                        configSolve["solve"]["methods"].append({
                            "name": "LP",
                            "verbose": verbose,
                            "timeLimits": timeLimit,
                            "relaxation": True,
                            "instances": [
                                {"path": f"{pathSaveInstance}/instance{i}{x}.txt"}
                            ]
                        })
                    elif method == "CG":
                        configSolve["solve"]["methods"].append({
                            "name": "CG",
                            "verbose": verbose,
                            "debug": True,
                            "gen_columns": 0,
                            "nbMinStateDP": nbMinStateDP,
                            "thresholdSetCol": 0.8,
                            "nbTimeNotUsed": 20,
                            "timeLimits": timeLimit,
                            "instances": [
                                {"path": f"{pathSaveInstance}/instance{i}{x}.txt"}
                            ]
                        })
                    elif method == "BaB_CG":
                        configSolve["solve"]["methods"].append({
                            "name": "BaB",
                            "verbose": verbose,
                            "strategy": str(strategy),
                            "scheme": "location",
                            "memorization":memo,
                            "LB_parameters": {
                                "name": "CG",
                                "verbose": verbose,
                                "debug": False,
                                "gen_columns": useHeuristic,
                                "maxNbCallHeuristic": maxNbCallHeuristic,
                                "nbMinStateDP": nbMinStateDP,
                                "thresholdSetCol": 0.8,
                                "nbTimeNotUsed": 20,
                                "timeLimits": 300
                            },
                            "timeLimits": timeLimit,
                            "instances": [
                                {"path": f"{pathSaveInstance}/instance{i}{x}.txt"}
                            ]
                        })
                    elif method == "BaB_MIP":
                        configSolve["solve"]["methods"].append({
                            "name": "BaB",
                            "verbose": verbose,
                            "strategy": str(strategy),
                            "scheme": "location",
                            "memorization": memo,
                            "LB_parameters": {
                                "name": "CG",
                                "verbose": 1,
                                "debug": False,
                                "gen_columns": useHeuristic,
                                "maxNbCallHeuristic": maxNbCallHeuristic,
                                "nbMinStateDP": nbMinStateDP,
                                "thresholdSetCol": 0.8,
                                "nbTimeNotUsed": 20,
                                "timeLimits": 300
                            },
                            "timeLimits": timeLimit,
                            "instances": [
                                {"path": f"{pathSaveInstance}/instance{i}{x}.txt"}
                            ]
                        })
                    elif method == "BeamSearch":
                        configSolve["solve"]["methods"].append({
                            "name": "BeamSearch",
                            "verbose": verbose,
                            "beamSize": beamSize,
                            "timeLimits": timeLimit,
                            "recovering": recovering==1,
                            "instances": [
                                {"path": f"{pathSaveInstance}/instance{i}{x}.txt"}
                            ]
                        })
                    elif method == "LocalSearch":
                        configSolve["solve"]["methods"].append({
                            "name": "LocalSearch",
                            "verbose": verbose,
                            "timeLimits": timeLimit,
                            "predictor": False,
                            "instances": [
                                {"path": f"{pathSaveInstance}/instance{i}{x}.txt"}
                            ]
                        })


                    solveConfigPath = pathlib.Path(str(pathlib.Path(str(pathlib.Path(__file__).parent.resolve()))) + f"/../instances/configs/{configFileNameSolve}i{i}{x}_memo_{memo}.json").resolve()
                    solveConfigPath.parent.mkdir(parents=True, exist_ok=True)
                    with open(solveConfigPath, "w") as file:
                        file.write(json.dumps(configSolve))

        dirOfConfiFiles = pathlib.Path(str(pathlib.Path(str(pathlib.Path(__file__).parent.resolve()))) + f"/../instances/configs/{configFileNameSolve}").resolve()
        if len(methods) > 0:
            print(f"- {dirOfConfiFiles}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate and solve scheduling instances.")

    # Add arguments
    parser.add_argument("--no-identical-jobs",type=int, default=0, help="Disable identical jobs, 0 for disable 1 for enable. By default, identical jobs are generated." )
    parser.add_argument("--memorization",type=int, default=1, help="Enable memorization, 0 for disable 1 for enable. By default, memorization is activated." )
    parser.add_argument("--tf", type=float, nargs='+', default=[0.2, 0.4, 0.6, 0.8, 1.0], help="TF values.")
    parser.add_argument("--rdd", type=float, nargs='+', default=[0.2, 0.4, 0.6, 0.8, 1.0], help="RDD values.")
    parser.add_argument("--n-list", type=int, nargs='+', default=[20], help="List of N values.")
    parser.add_argument("--frac-of-n", type=float, nargs='+', default=[0.5], help="Fraction of N jobs.")
    parser.add_argument("--m", type=str, nargs='+', default=["1,1", "2,2"], help="List of M values as 'mMax,M0'.")
    parser.add_argument("--vmax", type=float, default=2.0, help="Maximum velocity (Vmax).")
    parser.add_argument("--timeLimit", type=float, default=300, help="Maximum velocity (Vmax).")
    parser.add_argument("--v0", type=float, default=1.0, help="Initial velocity (V0).")
    parser.add_argument("--nb-instance-to-generate", type=int, default=5, help="Number of instances to generate.")
    parser.add_argument("--seed", type=int, default=1, help="Random seed.")
    parser.add_argument("--path-save-instance", type=str, default= "/instances/default/", help="Path to save instances. By default it's in directory ./instances/")
    parser.add_argument("--config-file-name-generate", type=str, default=None, help="Name of the generate config file.")
    parser.add_argument("--verbose", type=int, default=0, help="Verbosity level.")
    parser.add_argument("--output-result", type=str, default= "/instances/default/", help="Path to save results. By default it's in directory ./instances/default/")
    parser.add_argument("--config-file-name-solve", type=str, default=None, help="Name of the solve config file. If not set, then we don't generate solver configuration file")
    parser.add_argument("--methods", type=str, nargs='+', choices=["MIP", "LP", "CG", "BaB_CG","BaB_MIP","LocalSearch","BeamSearch"], default=[], help="Methods to use for solving.")
    parser.add_argument('--strategy',type=str,choices=['depth-first','breadth-first', 'best-first'],default='depth-first',help="The strategy to use for the BaB method. Options: 'depth-first','breadth-first', 'best-first'. Defaults to 'depth-first'. Only applicable if method is 'BaB'.")
    parser.add_argument('--use-heuristic-gen-col',type=int, default=0, help="Use an heuristic to generate columns in the Column Generation. By default, we don't use it." )
    parser.add_argument('--maxNbCallHeuristic',type=int, default=1, help="The maximum number of calls where the heuristic failed before use only dynamic programming . By default is 1" )
    parser.add_argument('--nbMinStateDP',type=int, default=1, help="The number of column that will be generated. By default is 1" )
    parser.add_argument('--beamSize',type=int, default=1, help="The beam size of the beamSearch algorithm. By default is 1" )
    parser.add_argument('--recovering',type=int, default=1, help="Diseable the recovering, 0 for disable 1 for enable. By default is 1" )
    args, leftovers = parser.parse_known_args()
    main(args)
