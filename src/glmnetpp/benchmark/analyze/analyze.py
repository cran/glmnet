# A simple script to run benchmark programs and visualize the data.
# Assumes that the benchmark program has already been built.

import argparse
import matplotlib.pyplot as plt
import path_names
import analyze_binomial_benchmark as abb
import analyze_gaussian_benchmark as agb
import analyze_binomial_two_class_benchmark as btcb

parser = argparse.ArgumentParser(description='Collects data and produces plots of benchmark programs.')
parser.add_argument('bench_names', nargs='*',
                    help='list of benchmark program names to analyze.')
parser.add_argument('-a', action='store_const', const=True,
                    help='analyze all benchmark programs in build/release/benchmark.')
args = parser.parse_args()

if len(args.bench_names) == 0 and not args.a:
    raise RuntimeError(
        'At least one benchmark name must be specified if -a is not specified.')

# Dictionary of bench name to module name
bench_to_module = {
    abb.TESTNAME : abb,
    agb.TESTNAME : agb,
    btcb.TESTNAME : btcb
}

mods = [bench_to_module[bench_name] for bench_name in args.bench_names]
for mod in mods:
    mod.plot(mod.run(path_names.bench_dir,
                     path_names.data_dir,
                     path_names.ref_dir,
                     path_names.data_scr_dir),
             path_names.fig_dir)
