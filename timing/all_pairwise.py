#!/usr/bin/env python
import subprocess
import argparse
import sys
import os

try:
    assert sys.version_info[0] == 3 and sys.version_info[1] >= 6
except AssertionError:
    raise ImportError("Error: Python 3.6+ required.")


SMD = {"original": " -E ", "ertl_improved": " -I ",
       "improved": " -I ", "ertl_mle": "", "ertl_jmle": " -J "}


def sm2str(sketchmethod):
    return SMD[sketchmethod.lower()]


def bashcall(cstr):
    return subprocess.check_call(cstr, shell=True, executable="/bin/bash")


BBITS = 16  # To make it evenly divisible for convenience in size comparisons.


def make_bd_ss64(hl2, bbits=BBITS):
    nbytes = (1 << hl2) >> 3
    nbytes //= BBITS
    nbytes //= 4  # I don't really know why, but this makes the sketch size equal what is expected.
    return nbytes


def submit_bindashonly(threads, ksz, sketchsz, outsizes, outdists,
                       fname_paths, logfile, sketchmethod, executable,
                       time_path, skip_already):
    bindashlog = logfile + ".bindashlog"
    ns64 = make_bd_ss64(sketchsz)  # 2**9 == 64 * 8, number of bytes
    if skip_already and os.path.isfile(bindashlog):
        sys.stderr.write("Now recalculating bindash sketch")
    else:
        cstr = f"nohup {time_path} -v bindash sketch --minhashtype=2 --nthreads={threads} --bbits={BBITS} --sketchsize64={ns64} --kmerlen={ksz} --outfname=bindash.{ksz}.{sketchsz}.bdsh --listfname={fname_paths} &> {bindashlog}"
        bashcall(cstr)
    bindashout = logfile + ".bindashout"
    if skip_already and os.path.isfile(bindashout):
        sys.stderr.write("Now recalculating bindash sketch")
    else:
        cstr = f"nohup {time_path} -v bindash dist --nthreads={threads} bindash.{ksz}.{sketchsz}.bdsh 1> {bindashout} 2>> {bindashlog}"
        bashcall(cstr)


def submit_mashonly(threads, ksz, sketchsz, outsizes, outdists, fname_paths,
                    logfile, sketchmethod, executable, time_path,
                    skip_already):
    mashlog = logfile + ".mashlog"
    if skip_already and os.path.isfile(mashlog):
        sys.stderr.write("Not recalculating mash sketch")
    else:
        cstr = f"nohup {time_path} -v mash sketch -p {threads} -s {1<<(sketchsz - 3 + (ksz > 16))} -k {ksz} -o mash.{ksz}.{sketchsz}.msh -l {fname_paths} &> {mashlog}"
        bashcall(cstr)
    mashout = logfile + ".mashout"
    if skip_already and os.path.isfile(mashout):
        sys.stderr.write("Not recalculating mash sketch")
    else:
        cstr = f"nohup {time_path} -v mash triangle -p {threads} mash.{ksz}.{sketchsz}.msh 1> {mashout} 2>> {mashlog}"
        bashcall(cstr)


def submit(threads, ksz, sketchsz, outsizes, outdists, fname_paths, logfile, sketchmethod, executable="dashing", time_path="/usr/bin/time", skip_already=True):
    assert os.path.isfile(fname_paths)
    sketchstr = sm2str(sketchmethod)
    cstr = f"nohup {time_path} -v {executable} sketch {sketchstr} -k{ksz} -p{threads} -S{sketchsz} -F{fname_paths} &> {logfile}"
    bashcall(cstr)
    cstr = f"nohup {time_path} -v {executable} dist {sketchstr} -fbWk{ksz} -p{threads} -S{sketchsz} -o{outsizes} -O{outdists} -F{fname_paths} &> {logfile}.tmp"
    bashcall(cstr)
    bashcall(f"cat {logfile}.tmp >> {logfile} && rm {logfile}.tmp")


def submit_dashing_nocache(threads, ksz, sketchsz, outsizes, outdists, fname_paths, logfile, sketchmethod, executable="dashing", time_path="/usr/bin/time", skip_already=True):
    assert os.path.isfile(fname_paths)
    sketchstr = sm2str(sketchmethod)
    logfile += ".both.log"
    cstr = f"nohup {time_path} -v {executable} dist {sketchstr} -fbk{ksz} -p{threads} -S{sketchsz} -o{outsizes} -O{outdists} -F{fname_paths} &> {logfile}"
    bashcall(cstr)

def submit_mash_nocache(threads, ksz, sketchsz, outsizes, outdists, fname_paths, logfile, time_path="/usr/bin/time", skip_already=True):
    assert os.path.isfile(fname_paths)
    mashlog = logfile + ".both.mashlog"
    mashout = logfile + ".both.mashout"
    # mash sketch -p {threads} -s {1<<(sketchsz - 3 + (ksz > 16))} -k {ksz}
    cstr = f"nohup {time_path} -v mash triangle -k {ksz} -p {threads} -s {1<<(sketchsz - 3 + (ksz > 16))} -l {fname_paths} 1> {mashout} 2> {mashlog}"
    bashcall(cstr)


def main():
    a = argparse.ArgumentParser()
    a.add_argument("fname_paths", help="File with paths to all sketched genomes")
    a.add_argument("--sketch-sizes", help="Comma-delimited sketch sizes to use.", required=True)
    a.add_argument("--kmer-sizes", help="Comma-delimited kmer sizes to use.", required=True)
    a.add_argument("--threads", "-p", type=int, help="Number of threads to use", default=100)
    a.add_argument("--num-times", default=1, type=int)
    a.add_argument("--executable", default="dashing", type=str)
    a.add_argument("--time-path", default="/usr/bin/time", type=str)
    a.add_argument("--bdonly", action="store_true", help="Only perform bindash experiments.")
    args = a.parse_args()
    args.skip_already_done = True
    sizes = list(map(int, args.sketch_sizes.split(",")))
    ks = list(map(int, args.kmer_sizes.split(",")))
    fname_paths = args.fname_paths
    bn = "all_against_all_" + os.path.basename(fname_paths).split(".")[0]
    sketch_methods = ["ertl_mle", "original", "ertl_jmle", "improved"]
    assert sizes
    assert ks
    print(f"Running {len(sizes) * len(ks) * len(sketch_methods)} runs", file=sys.stderr)
    run_num = 0
    for size in sizes:
        for k in ks:
            for sm in sketch_methods:
                for i in range(1, args.num_times + 1):
                    logf = f"{bn}.{size}.{k}.{sm}.{i}.log"
                    baselogf = f"{bn}.{size}.{k}.{i}.log"
                    #if args.skip_already_done:
                    if os.path.isfile(logf) and os.path.isfile(baselogf + ".mashlog") and os.path.isfile(baselogf + ".bindashlog"):
                        print(f"Skipping already done experiment {run_num} at index {i}", file=sys.stderr)
                    else:
                        if all(not os.path.isfile(x) for x in (logf, baselogf + ".mashlog", baselogf + ".bindashlog")):
                            submit(args.threads, k, size, f"{bn}.{size}.{k}.{sm}.{i}.txt", f"{bn}.{size}.{k}.{sm}.{i}.dist.bin",
                                   args.fname_paths, logf, sm, executable=args.executable, time_path=args.time_path, skip_already=args.skip_already_done)
                        elif not os.path.isfile(baselogf + ".mashlog"):
                            submit_mashonly(args.threads, k, size, f"{bn}.{size}.{k}.{sm}.{i}.txt", f"{bn}.{size}.{k}.{sm}.{i}.dist.bin", args.fname_paths, baselogf, sm, args.executable, args.time_path, args.skip_already_done)
                        if not os.path.isfile(baselogf + ".bindashlog"):
                            submit_bindashonly(args.threads, k, size, f"{bn}.{size}.{k}.{sm}.{i}.txt", f"{bn}.{size}.{k}.{sm}.{i}.dist.bin", args.fname_paths, baselogf, sm, args.executable, args.time_path, args.skip_already_done)
                    #else:
                    #    submit(args.threads, k, size, f"{bn}.{size}.{k}.{sm}.{i}.txt", f"{bn}.{size}.{k}.{sm}.{i}.dist.bin", args.fname_paths, logf, sm, executable=args.executable, time_path=args.time_path, skip_already=args.skip_already_done)
                    #    submit_mashonly(args.threads, k, size, f"{bn}.{size}.{k}.{sm}.{i}.txt", f"{bn}.{size}.{k}.{sm}.{i}.dist.bin", args.fname_paths, baselogf, sm, args.executable, args.time_path, args.skip_already_done)
                    #    submit_bindashonly(args.threads, k, size, f"{bn}.{size}.{k}.{sm}.{i}.txt", f"{bn}.{size}.{k}.{sm}.{i}.dist.bin", args.fname_paths, baselogf, sm, args.executable, args.time_path, args.skip_already_done)
                    submit_dashing_nocache(args.threads, k, size, f"{bn}.{size}.{k}.{sm}.{i}.both.txt", f"{bn}.{size}.{k}.{sm}.{i}.both.dist.bin",
                                           args.fname_paths, baselogf, sm, executable=args.executable, time_path=args.time_path, skip_already=args.skip_already_done)
                    submit_mash_nocache(args.threads, k, size, f"{bn}.{size}.{k}.{sm}.{i}.txt", f"{bn}.{size}.{k}.{sm}.{i}.dist.bin",  args.fname_paths, baselogf, args.time_path)

                    run_num += 1
                    print(f"Successfully ran experiment {run_num} in batch {i}", file=sys.stderr)


if __name__ == "__main__":
    sys.exit(main())
