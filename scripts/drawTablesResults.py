#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
from itertools import repeat
from functools import reduce
import os

# Configuration parameters
N_VALUES = [40,50,60,70,80]
M_VALUES = [2, 4]
MAX_NB_CALL_HEURISTIC = 1  # max number of failing of heuristic we tolerate
NB_MIN_COLUMN = 3  # the minimal number of column we generate with DP (generate several with backtracking)
METHODS = ['MIP', 'BaB']
STRATEGIES = ['depth-first', np.nan]
LOWER_BOUNDS = ['Columns_Generation', np.nan]
MEMORIZATION = True
HEURISTIC = 0
OUTPUT_FILE = "temp/test.tex"
BASE_PATH = "/home/schau/Documents/These/Code/Bilevel-Scheduling-I4.0/instances/N{}/{}"

if not os.path.exists(os.path.dirname(OUTPUT_FILE)):
    os.makedirs(os.path.dirname(OUTPUT_FILE))

def load_data(n_values, m_values):
    """Load and concatenate all CSV files into a single DataFrame."""
    dfs = []
    for n in n_values:
        for m in m_values:
            # Load MIP and BaB files for each (N, m) combination
            mip_path = BASE_PATH.format(n, f"{m}M_resultsMIP.csv")
            bab_path = BASE_PATH.format(n, f"{m}M_resultsBaB.csv")
            dfs.append(pd.read_csv(mip_path, sep='\t'))
            dfs.append(pd.read_csv(bab_path, sep='\t'))
    return pd.concat(dfs, ignore_index=True)


def calculate_gap(row):
    """Calculate optimality gap percentage."""
    if row['isOptimal'] == 1:
        return 0.0
    elif row['bestKnowObj'] > 0:
        return ((row['Objective'] - row['bestKnowObj']) / row['bestKnowObj']) * 100.0
    return None


def NbOPT(x, y):
    return f"{int(x)}/{int(y)}"


def filter_dataframe(df,nb_machines, ):
    """
    Filter dataframe based on specified criteria.

    Args:
        is_optimal: Whether to include only instances solved optimally by all methods
        nb_machines: Number of machines to filter by
    """
    cond = (df['Method'].isin(METHODS) | (df["Method"].isna()))
    cond &= df['$m$'] == nb_machines
    cond &= (df['Memorization'] == MEMORIZATION) | (df["Memorization"].isna())
    cond &= (df['GenerateCol'] == HEURISTIC) | (df["GenerateCol"].isna())
    cond &= (df["Strategy"].isin(STRATEGIES) | (df["Strategy"].isna()))
    cond &= (df["LowerBound"].isin(LOWER_BOUNDS) | (df["LowerBound"].isna()))
    cond &= (df["maxNbCallHeuristics"] == MAX_NB_CALL_HEURISTIC) | (df["maxNbCallHeuristics"].isna())
    cond &= (df["NBMinColum"] == NB_MIN_COLUMN) | (df["NBMinColum"].isna())
    filtered = df[cond]
    return filtered


def count_instances(df):
    """
    Count the number of unique instances based on method type.

    For MIP: Group by (N, n, $m$)
    For non-MIP: Group by (N, n, $m$, LowerBound, Strategy, GenerateCol, maxNbCallHeuristics, NBMinColum)

    Returns:
        int: The uniform count of instances if all groups agree

    Raises:
        ValueError: If different groups have different instance counts
    """
    # Define base columns that always apply
    base_cols = ['N', 'n', '$m$']

    # Check if DataFrame contains any MIP methods
    has_mip = 'MIP' in df['Method'].values
    counts = None
    if has_mip:
        # For MIP: Count using only base columns
        counts = df.groupby(base_cols)['InstanceName'].nunique()
    else:
        # For non-MIP: Use extended set of columns
        extended_cols = base_cols + [
            'LowerBound', 'Strategy', 'GenerateCol',
            'maxNbCallHeuristics', 'NBMinColum'
        ]
        counts = df.groupby(extended_cols)['InstanceName'].nunique()

    if counts is None or not counts.nunique() == 1:
        raise ValueError(f"Inconsistent instance counts found:\n{counts}")

    # Return the uniform count
    return counts.iloc[0]


def generate_table(df, metric, only_OPT_instance=False):
    """
    Generate a statistics table for the specified metric and condition.

    Args:
        metric: Column name to analyze ('Time', 'NBNodes', 'GAP')
        condition: 'optimal' or 'non_optimal'
        only_OPT_instance: compare only optimal solved instances
    """
    groupByIndex = ['N', 'n',"$m$"]  # group index of the table by N and n
    first_col = []
    if (only_OPT_instance):
        first_col += ["\\#OPT"]
        first_col += reduce(lambda x, y: x + y, [list(repeat(x, 3)) for x in METHODS])
    else:
        first_col += reduce(lambda x, y: x + y, [list(repeat(x, 4)) for x in METHODS])
    second_col = []
    if (only_OPT_instance):
        second_col += ["\\#OPT"]
        second_col += reduce(lambda x, y: x + y, list(
            repeat([f"${metric}_{{min}}$", f"${metric}_{{avg}}$", f"${metric}_{{max}}$"], len(METHODS))))
    else:
        second_col += reduce(lambda x, y: x + y, list(
            repeat(["\\#OPT", f"${metric}_{{min}}$", f"${metric}_{{avg}}$", f"${metric}_{{max}}$"], len(METHODS))))

    arrays = [np.array(first_col), np.array(second_col)]
    df_table = []
    # for each set of machines
    for nbMachine in M_VALUES:
        # sort the dataframe to keep the data with the right number of machines and all optimal solved instances
        df_all_opt = filter_dataframe(df, nbMachine)
        if (only_OPT_instance):
            df_all_opt= df_all_opt.groupby('InstanceName').filter(lambda group: group['isOptimal'].all())
        else:
            df_all_opt= df_all_opt.groupby('InstanceName').filter(lambda group: not group['isOptimal'].all())
        if (only_OPT_instance and df_all_opt.groupby('InstanceName').filter(lambda group: group['Objective'].nunique()==1).shape != df_all_opt.shape):
            raise ValueError(f"Instance indicate to solve optimaly have not same result with all methods")
        # We use the index of df_all_opt, i.e., (N,n) and we set the desired column with 'arrays'
        df_table_machine = pd.DataFrame(index=df_all_opt.groupby(groupByIndex).describe().index, columns=arrays)
        # Change the value of the first column, instead of an index, we set the number of machines
        df_table_machine.iloc[:, 0] = [nbMachine for _ in range(df_table_machine.shape[0])]
        # define filter to select only instances with the right maxNbCallHeuristic and nbMinColumn
        filterCondition = (df_all_opt["maxNbCallHeuristics"] == MAX_NB_CALL_HEURISTIC) & (
                df_all_opt["NBMinColum"] == NB_MIN_COLUMN)

        if (only_OPT_instance):
            # compute the number of instance to solve for a set of parameters
            nbInstanceToSolve = count_instances(df)
            # use the number of optimal instance solve by the MIP, to get the number of optimal instance for both method
            df_table_machine.iloc[:, 0] = df_all_opt[df_all_opt["Method"]=="MIP"].groupby(groupByIndex).describe().sort_index(level=groupByIndex).loc[:, (metric, 'count')].apply(lambda x: NbOPT(x, nbInstanceToSolve))

        for i, method in enumerate(METHODS):
            if (method != "MIP"):
                filter = (df_all_opt["Method"] == method) & filterCondition
            else:
                filter = (df_all_opt["Method"] == method)
            # for each method, compute the dataframe of statistics
            df_Method_stats = df_all_opt[filter].groupby(groupByIndex).describe()

            # if we compare only optimal instance, then we compute only one time the number of optimal solved instance
            if (only_OPT_instance):
                df_table_machine.iloc[:, (i * 3) + 1] = df_Method_stats.loc[:, (metric, 'min')]
                df_table_machine.iloc[:, (i * 3) + 2] = df_Method_stats.loc[:, (metric, 'mean')]
                df_table_machine.iloc[:, (i * 3) + 3] = df_Method_stats.loc[:, (metric, 'max')]
            else:
                indexSelection = df_all_opt[groupByIndex].drop_duplicates().set_index(groupByIndex).sort_index().index
                nbOPT = df_all_opt[filter & (df_all_opt["isOptimal"] == True)].groupby(groupByIndex).count().sort_index(level=groupByIndex).iloc[:, 0].reindex(indexSelection, fill_value=0)
                nbInstanceToSolve= df_all_opt[filter & (df_all_opt["Method"] == method)].groupby(groupByIndex).count().sort_index(level=groupByIndex).iloc[:, 0].reindex(indexSelection, fill_value=0)
                df_table_machine.iloc[:, (i * 4) + 0] = [f"{int(x)}/{int(y)}" for x, y in zip(nbOPT, nbInstanceToSolve)]
                df_table_machine.iloc[:, (i * 4) + 1] = df_Method_stats.loc[:, (metric, 'min')]
                df_table_machine.iloc[:, (i * 4) + 2] = df_Method_stats.loc[:, (metric, 'mean')]
                df_table_machine.iloc[:, (i * 4) + 3] = df_Method_stats.loc[:, (metric, 'max')]

        df_table.append(df_table_machine)

    return pd.concat(df_table)


def write_latex_table(f, df, caption=None):
    """Write DataFrame as LaTeX table to file."""
    f.write(r"\begin{table}[ht!]" + "\n")
    if caption:
        f.write(r"\caption{" + caption + "}\n")
    f.write(df.to_latex(index=True, float_format=lambda x: str.format("{:.2f}", x), multirow=True, multicolumn=True,
                        multicolumn_format="c"))
    f.write(r"\end{table}" + "\n\n")


def main():
    # Load and preprocess data
    df = load_data(N_VALUES, M_VALUES)
    min_objective = df.groupby('InstanceName')['Objective'].min()
    df['bestKnowObj'] = df['InstanceName'].map(min_objective)
    df['GAP'] = df.apply(calculate_gap, axis=1)
    df["$m$"] = df["m_Max"] + df["m_0"]
    df.rename(columns={
        "m_Max": "$N_{max}$",
        "m_0": "$N_0$",
        "V_max": "$V_{max}$",
        "V_0": "$V_0$"
    }, inplace=True)

    # Generate LaTeX report
    with open(OUTPUT_FILE, 'w') as f:
        f.write(r"""\documentclass[11pt]{article}
        \usepackage{array}
        \usepackage{pdflscape}
        \usepackage{booktabs}
        \usepackage{multirow}
        \begin{document}
        \begin{landscape}
        """)
        # Time tables
        for condition in ['non-optimal','optimal']:
            df_table = generate_table(
                df, 'Time', (condition == 'optimal')
            )
            write_latex_table(f, df_table,
                              f"Time Statistics all {condition} for {MAX_NB_CALL_HEURISTIC} maxNbCall and {NB_MIN_COLUMN} nb min col)")

        # Node tables
        for condition in ['non-optimal', 'optimal']:
            df_table = generate_table(
                df, 'NBNodes', (condition == 'optimal')
            )
            write_latex_table(f, df_table,
                              f"NBNodes Statistics all {condition} for {MAX_NB_CALL_HEURISTIC} maxNbCall and {NB_MIN_COLUMN} nb min col)")

        # Gap table (only for non-optimal)
        for condition in ['non-optimal', 'optimal']:
            df_table = generate_table(
                df, 'GAP', (condition == 'optimal')
            )
            write_latex_table(f, df_table,
                              f"GAP Statistics all {condition} for {MAX_NB_CALL_HEURISTIC} maxNbCall and {NB_MIN_COLUMN} nb min col)")

        f.write(r"""\end{landscape}
                \end{document}""")


if __name__ == "__main__":
    main()