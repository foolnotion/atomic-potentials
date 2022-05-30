# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: Copyright 2019-2022 Heal Research

from colorama import Back as bg
from colorama import Fore as fg
from colorama import Style as st
from colorama import init
import argparse
import coloredlogs
import itertools
import json
import logging
import math
import os
import pandas as pd
import numpy as np
import subprocess
import sys
import re

import optuna

init(autoreset=True)

parser = argparse.ArgumentParser()
parser.add_argument('--bin', help='Path to algorithm executable', type=str)
parser.add_argument('--energy', help='Path to the file containing the energy', type=str)
parser.add_argument('--coordinates', help='Path to the file containing the atomic coordinates', type=str)
parser.add_argument('--evaluations', help='Evaluation budget', type=str)
parser.add_argument('--train', help='Training range', type=str)
parser.add_argument('--target', help='Regression target', type=str)
parser.add_argument('--test', help='Test range', type=str)
parser.add_argument('--trials', help='The number of trials for optuna', type=int)
parser.add_argument('--prefix', help='Prefix to add to output filenames', type=str)
parser.add_argument('--out', help='Location where the produced result files should be saved', type=str)

args = parser.parse_args()

bin_path     = args.bin
energy_path  = args.energy
coords_path  = args.coordinates
base_path    = os.path.dirname(energy_path)
trials       = args.trials
prefix       = args.prefix
target       = args.target
results_path = args.out

#coloredlogs.install(stream=sys.stdout, level=logging.INFO)
logger = logging.getLogger("atomic")

rng = np.random.default_rng()

def is_float(v):
    try:
        float(v)
        return True
    except ValueError:
        return False


def run(params):
    cmd = ['stdbuf', '-oL',
            bin_path,
            '--shuffle', '--symbolic'
            ]

    for p in params:
        cmd += [ f'--{p}', str(params[p]) ]

    print(' '.join(cmd))
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True);
    for l in iter(p.stdout.readline, ''):
        yield l
    p.stdout.close()
    p.wait()


params = {
        'dataset' : energy_path,
        'coordinates' : coords_path,
        'cluster-size' : '32',
        'inputs' : 'r,q',
        'target' : target,
        'train' : args.train,
        'test' : args.test,
        'iterations' : '0',
        'evaluations' : str(int(1e7)),
        'generations' : str(1000),
        'offspring-generator' : 'basic',
        'reinserter' : 'keep-best',
        'maxdepth' : str(20),
        'seed' : str(rng.bit_generator.random_raw()),
        'threads' : '16',
        }


def objective(trial):
    params['population-size'] = trial.suggest_int("population-size", 100, 1000, 100) 
    params['pool-size'] = trial.suggest_int("pool-size", 100, 1000, 100)
    #params['iterations'] = str(trial.suggest_int("iterations", 0, 10, 1))
    params['maxlength'] = trial.suggest_int("maxlength", 10, 30, 5)
    params['enable-symbols'] = trial.suggest_categorical("enable-symbols", ['add,sub,mul,div', 'add,sub,mul,div,pow', 'add,sub,mul,aq', 'add,sub,mul,aq,pow'])
    params['epsilon'] = trial.suggest_categorical("epsilon", [ 1e-5, 1e-4, 1e-3 ])

    #seed = rng.bit_generator.random_raw()
    lines = []
    for l in run(params):
        print(l, end='')
        lines.append(l)

    model = lines[-1]
    stats = lines[-2]
    vals = re.split('\s+', stats) 
    return float(vals[2])

study = optuna.create_study(direction='maximize')
study.optimize(objective, n_trials=trials)

params = params | study.best_params
print('params:', params)

output_header = [
        'run',
        'iter',
        'r2_tr',
        'r2_te',
        'mae_tr',
        'mae_te',
        'nmse_tr',
        'nmse_te',
        'avg_fit',
        'avg_len',
        'eval_cnt',
        'res_eval',
        'jac_eval',
        'seed',
        'elapsed',
        'model'
        ]

data = []
for i in range(50):
    logger.info(f'run {i} with params {params}')
    if i > 0:
        params['seed'] = rng.bit_generator.random_raw()

    lines = []
    for l in run(params):
        print(l, end='')
        lines.append(l)

    stats = [i] + [float(x) for x in re.split('\s+', lines[-2]) if is_float(x)]
    df = pd.DataFrame({output_header[i] : stats[i] for i in range(len(stats))}, index=[i])
    df['model'] = lines[-1].strip()
    data.append(df)

df = pd.concat(data)
df.sort_values(by=['mae_te'], inplace=True)
df.to_csv(f'{rng.bit_generator.random_raw()}.csv', index=None)

print('best parameters')
print(params | study.best_params)
