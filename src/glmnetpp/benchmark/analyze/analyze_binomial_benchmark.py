import io
import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from subprocess import check_output

TESTNAME = 'binomial_benchmark'
n = 1000
ps = 2**(np.linspace(5, 11, 7)).astype(int)

def plot(df, fig_dir):
    grouped = df.groupby(['sp', 'kopt', 'two_class'])
    nrows=2
    ncols = int(np.ceil(grouped.ngroups/nrows))

    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(12,6), sharey=True)

    for (key, ax) in zip(grouped.groups.keys(), axes.flatten()):
        print(grouped.get_group(key))
        grouped.get_group(key).plot(ax=ax, x='p', y='relative')
        ax.set_title('kopt={kopt}, do_2={do_two_class}, sp={sp}'.format(sp=key[0], do_two_class=key[2], kopt=key[1]))
        ax.legend()
        ax.set_ylabel('Time (s)')

    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, TESTNAME + '_fig.png'))

# Run benchmark
# bench_dir     directory to glmnetpp benchmark program (e.g. build/release/benchmark)
# data_dir      directory to store our timing data (e.g. docs/data)
# ref_dir       directory to reference (glmnet) program for comparison (e.g. benchmark/reference)
# data_scr_dir  directory to scripts that generate data (e.g. benchmark/data/script)
# gen           boolean whether to generate data or not
def run(bench_dir, data_dir, ref_dir, data_scr_dir, gen=False):
    df = pd.DataFrame()

    # save current working directory
    cur_path = os.getcwd()

    # change directory to glmnetpp benchmark location
    os.chdir(bench_dir)

    bench_path = os.path.join('.', TESTNAME)
    print('Benchmark path: {p}'.format(p=bench_path))

    # run our benchmark and get output
    args = (bench_path, "--benchmark_format=csv")
    data = io.StringIO(check_output(args).decode("utf-8"))
    os.chdir(cur_path)

    df_bench = pd.read_csv(data, sep=',')
    n = df_bench['n'][0]    # assume n is constant throughout
    df = df_bench[['p', 'kopt', 'two_class', 'sp', 'real_time']]
    df_glmnet = df[df_bench['glmnetpp'] == 0]
    df_glmnetpp = df[df_bench['glmnetpp'] == 1]
    df_glmnet.set_index(['p', 'kopt', 'two_class', 'sp'], inplace=True)
    df_glmnetpp.set_index(['p', 'kopt', 'two_class', 'sp'], inplace=True)
    df = pd.concat([df_glmnet, df_glmnetpp], axis=1)
    df.columns = ['glmnet', 'glmnetpp']
    df *= 1e-9
    df.reset_index(inplace=True)

    df['relative'] = df['glmnet'] / df['glmnetpp']

    # save absolute time
    data_path = os.path.join(data_dir, TESTNAME + ".csv")
    df.to_csv(data_path)

    return df

if __name__ == '__main__':
    import path_names as pn
    df = pd.read_csv(os.path.join(pn.data_dir, TESTNAME + ".csv"))
    plot(df, pn.fig_dir)
