import io
import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from subprocess import check_output

TESTNAME = 'binomial_two_class_benchmark'
n = 1000
ps = 2**(np.linspace(5, 11, 7)).astype(int)

def plot(df, fig_dir):
    ax = df.plot(x='p', y='relative')
    ax.set_title('N=1000')
    ax.legend()
    ax.set_ylabel('Relative Time')

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
    df = df_bench[['p', 'real_time']]
    df_glmnet = df[df_bench['glmnetpp'] == 0]
    df_glmnetpp = df[df_bench['glmnetpp'] == 1]
    df_glmnet.set_index(['p'], inplace=True)
    df_glmnetpp.set_index(['p'], inplace=True)
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
