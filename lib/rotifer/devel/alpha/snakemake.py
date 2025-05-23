from rotifer.db import ncbi
import pandas as pd
from rotifer.core import config as CoreConfig
from rotifer.devel.beta import sequence as rdbs
from rotifer import config
import shutil
import tempfile
import subprocess
from pathlib import Path
from functools import wraps
import pickle



def run_snakemake_in_tmp(snake_src_dir, result_file_relpath, snakemake_args=None, input_files=None):
    """
    Run a Snakemake pipeline in a temporary folder and copy result file back.

    Parameters:
    - snake_src_dir: Path to Snakemake project
    - result_file_relpath: Relative path to output file (from Snakemake root)
    - snakemake_args: List of command-line arguments to Snakemake
    - input_files: List of (src, rel_dest) tuples to copy into the temp dir
                   Example: [("/path/to/input1.txt", "data/input1.txt")]
    """
    snake_src_dir = Path(snake_src_dir).resolve()
    result_file_relpath = Path(result_file_relpath)
    snakemake_args = snakemake_args or []
    input_files = input_files or []

    with tempfile.TemporaryDirectory(dir=f"{config['cache']}") as tmpdirname:
        tmpdir = Path(tmpdirname)
        print(f"[INFO] Temporary Snakemake directory: {tmpdir}")

        # Copy Snakemake pipeline to temp folder
        tmp_snakemake_dir = tmpdir / snake_src_dir.name
        shutil.copytree(snake_src_dir, tmp_snakemake_dir)

        # Copy input files into appropriate relative locations
        for src, rel_dest in input_files:
            src_path = Path(src).resolve()
            dest_path = tmp_snakemake_dir / rel_dest
            dest_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src_path, dest_path)
            print(f"[INFO] Copied input: {src_path} → {dest_path}")

        # Build and run the command
        cmd = ["snakemake", "--snakefile", "Snakefile", *snakemake_args]
        print(f"[INFO] Running: {' '.join(cmd)}")
        subprocess.run(cmd, cwd=tmp_snakemake_dir, check=True)

        # Copy result file back to current working directory
        src_result = tmp_snakemake_dir / result_file_relpath
        #dest_result = Path.cwd() / result_file_relpath.name
        #shutil.copy2(src_result, dest_result)
        print(src_result)
        if src_result.suffix =='.pkl':
            with open(src_result, "rb") as f:
                obj = pickle.load(f)
        elif src_result.suffix =='.tsv':
            obj = pd.read_csv(src_result,sep="\t", header=None)

        return obj



## The Decorator factory to re uses the more broader snakemake function:
from functools import wraps

def snakemake_pipeline(snake_src_dir, result_file_relpath):
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            result = func(*args, **kwargs)

            if not isinstance(result, dict):
                raise ValueError("Decorated function must return a dict with keys like 'snakemake_args', 'input_files', 'cores'")

            args_list = result.get("snakemake_args", [])
            input_files = result.get("input_files", [])
            cores = result.get("cores", None)

            if cores is not None:
                args_list = ["--cores", str(cores)] + args_list
            elif not any(arg.startswith("--cores") or arg.startswith("-c") for arg in args_list):
                args_list = ["--cores", "1"] + args_list  # default fallback

            return run_snakemake_in_tmp(
                snake_src_dir=snake_src_dir,
                result_file_relpath=result_file_relpath,
                snakemake_args=args_list,
                input_files=input_files
            )
        return wrapper
    return decorator

@snakemake_pipeline(
    snake_src_dir=f"{CoreConfig['baseDataDirectory']}/snakemake/neighborhood",
    result_file_relpath="results/all_neighborhood.pkl"
)
def run_neighbor(input_file, cores=2):
    # Automatically detect if the input is a list, pandas dataframe or string to run the functions
    temp_input_path = None
    if isinstance(input_file, (list, pd.Series)):
        temp = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt')
        for item in input_file:
            temp.write(f"{str(item)}\n")
        temp.close()
        temp_input_path = temp.name
    else:
        temp_input_path = input_file
    return {
        "snakemake_args": ["--config", "input_file=data/input.txt"],
        "input_files": [(temp_input_path, "data/input.txt")],
        "cores": cores
    }

@snakemake_pipeline(
    snake_src_dir=f"{CoreConfig['baseDataDirectory']}/snakemake/neighborhood",
    result_file_relpath="results/all_neighborhood.pkl"
)
def run_neighbor_farm(input_file,
                      cores=48,
                      batch_size=100,
                      jobs=100):
    batch_size = str(batch_size)
    jobs = str(jobs)
    # Automatically detect if the input is a list, pandas dataframe or string to run the functions
    temp_input_path = None
    if isinstance(input_file, (list, pd.Series)):
        temp = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt')
        for item in input_file:
            temp.write(f"{str(item)}\n")
        temp.close()
        temp_input_path = temp.name
    else:
        temp_input_path = input_file
    return {
        "snakemake_args": ["--config", "input_file=data/input.txt",
                           f"batch_size={batch_size}",
                           f"jobs={jobs}",
                           "--profile","profiles/farm"],
        "input_files": [(temp_input_path, "data/input.txt")],
        "cores": cores
    }

@snakemake_pipeline(
    snake_src_dir=f"{CoreConfig['baseDataDirectory']}/snakemake/rps",
    result_file_relpath="final.arch.tsv"
)
def run_dommain_annotation_farm(input_file,
                      cores=48,
                      batch_size=40,
                      profiledb=True,
                      pfam=True,
                      hmm_db_profiledb="/netmnt/vast01/cbb01/proteinworld/data/rpsdb/allprofiles",
                      hmm_db_pfam="/netmnt/vast01/cbb01/proteinworld/data/rpsdb/pwld_new_pfam",      
                      jobs=100):
    batch_size = str(batch_size)
    jobs = str(jobs)
    # Automatically detect if the input is a list, pandas dataframe or string to run the functions
    temp_input_path = None
    if isinstance(input_file, rdbs.sequence):
        temp = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt')
        input_file.to_file(temp.name)
        temp.close()
        temp_input_path = temp.name
    else:
        temp_input_path = input_file
    return {
        "snakemake_args": ["--config", "input_fasta=data/fasta.fa",
                           f"run_profiledb={profiledb}",
                           f"hmmdb_profiledb={hmm_db_profiledb}",
                           f"run_pfam={pfam}",
                           f"hmmdb_pfam={hmm_db_pfam}",
                           f"batch_size={batch_size}",
                           f"jobs={jobs}",
                           "--profile","profiles/farm"],
        "input_files": [(temp_input_path, "data/fasta.fa")],
        "cores": cores
    }
