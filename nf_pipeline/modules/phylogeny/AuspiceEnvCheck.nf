#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process AuspiceEnvCheck {
    // Builds (on first use) and verifies the Nextstrain env of the Auspice export ONCE per --auspice run, at the start,
    // in parallel with PhyloPreflight and before any tree work (plan D21, user request 2026-09-30: "make sure that the
    // env can be built always"). PhyloGeography and PhyloAuspice take this process's report as an input, so no
    // geography or per-tree export task ever starts without a verified env, and a broken env is one early, actionable
    // run failure instead of a "no Auspice dataset" line per tree.
    // The env is envs/nextstrain.yaml (Nextflow builds it into conda.cacheDir, keyed on the file's content) or, with
    // --auspiceCondaEnv, an existing pre-built prefix (offline / HPC: nothing is built; the up-front parameter check
    // already rejected a path that is not a conda prefix). Verified here: the augur CLI is the env's own (under
    // CONDA_PREFIX), python3 runs from the env (sys.prefix), and treetime (the library augur imports; the activated
    // FluTyper env has its own treetime 0.11.4), geonamescache, Bio and pandas import from it. The versions are compared with the pins of envs/nextstrain.yaml: a
    // difference (possible with --auspiceCondaEnv) is a WARNING, not a failure.
    // Verified on Nextflow 26.04.6 (D21): a failure to CREATE the env (solver/download error, bad prefix) fails the
    // whole run before this script runs and bypasses errorStrategy (no retry, not even with 'ignore'); main.nf's
    // onComplete then prints auspiceEnvHelp(). The retry below covers this script's own checks (e.g. a transient
    // filesystem error); after it the run terminates, never 'ignore'.
    conda params.auspiceCondaEnv ?: "${moduleDir}/../../../envs/nextstrain.yaml"
    errorStrategy { task.attempt <= 1 ? 'retry' : 'terminate' }
    maxRetries 1

    input:
    path(env_spec, stageAs: "in_spec/*")   // envs/nextstrain.yaml: the pins the versions are compared with

    output:
    path("nextstrain_env.txt"), emit: versions   // AuspiceEnvCheck INFO/WARNING lines: prefix + versions (into every phylo.log)

    script:
    """
    fail() {
        cat >&2 <<'MSG'
AuspiceEnvCheck: ERROR the Nextstrain env is not usable (augur/treetime/geonamescache missing or broken).
${auspiceEnvHelp()}
MSG
        echo "AuspiceEnvCheck: details: \$1" >&2
        exit 1
    }
    [ -n "\${CONDA_PREFIX:-}" ] || fail "no conda env is active (CONDA_PREFIX is empty)"
    # augur: the CLI PhyloAuspice calls must be the env's own
    where=\$(command -v augur) || fail "augur is not on PATH"
    case "\$where" in "\$CONDA_PREFIX"/*) ;; *) fail "augur resolves to \$where, outside the env \$CONDA_PREFIX" ;; esac
    augur --version > augur.version 2>&1 || fail "augur --version: \$(tail -3 augur.version)"
    treetime --version > treetime.version 2>&1 || fail "treetime --version: \$(tail -3 treetime.version)"
    # python3 (what PhyloGeography/PhyloAuspice call): checked by the interpreter that actually RUNS (sys.prefix), not by
    # where the name resolves: conda swaps the env into the activated FluTyper env's PATH slot, which can sit behind
    # e.g. pyenv shims that then delegate to it (seen on the dev machine, D21)
    python3 - "${env_spec}" > nextstrain_env.txt 2> python.err <<'PYEOF' || fail "python3 check: \$(tail -3 python.err)"
import os, re, sys
prefix = os.path.realpath(os.environ["CONDA_PREFIX"])
if os.path.realpath(sys.prefix) != prefix:
    sys.exit(f"python3 runs from {sys.prefix}, not from the env {prefix}")
import geonamescache, Bio, pandas, treetime
from importlib.metadata import version
spec = open(sys.argv[1]).read()
def pin(name, sep="="):
    m = re.search(r"^\\s*-\\s*" + re.escape(name) + re.escape(sep) + r"([^\\s#]+)", spec, re.M)
    return m.group(1) if m else None
def first_version(path):
    m = re.search(r"(\\d+\\.\\d+(?:\\.\\d+)?)", open(path).read())
    return m.group(1) if m else "unknown"
found = [("augur", first_version("augur.version"), pin("augur")),
         ("treetime", treetime.version, pin("treetime")),    # the library augur imports
         ("python", sys.version.split()[0], pin("python")),
         ("geonamescache", version("geonamescache"), pin("geonamescache", "==")),
         ("biopython", Bio.__version__, pin("biopython")),
         ("pandas", pandas.__version__, pin("pandas"))]
print(f"AuspiceEnvCheck: INFO Nextstrain env verified at {prefix}: " + ", ".join(f"{n} {v}" for n, v, _ in found) + ".")
cli = first_version("treetime.version")
if cli != treetime.version:
    print(f"AuspiceEnvCheck: WARNING the treetime command on PATH is {cli}, the library augur uses is {treetime.version} "
          "(another env's treetime shadows the env's script; augur is unaffected).")
for n, v, p in found:
    if p and v != p:
        print(f"AuspiceEnvCheck: WARNING {n} {v} differs from the envs/nextstrain.yaml pin {p} (results may differ "
              "from the tested versions).")
PYEOF
    """
}

// Likely causes and fixes of an unusable/unbuildable Nextstrain env: printed by AuspiceEnvCheck when its checks fail,
// and by main.nf's onComplete when Nextflow could not create the env at all (that failure happens before any script).
def auspiceEnvHelp() {
    return """\
The Auspice export (--auspice) needs the Nextstrain env of envs/nextstrain.yaml (augur, treetime, geonamescache).
Likely causes and fixes:
  - No internet on the first build: the env is built once (conda packages + pip geonamescache) and then reused from
    conda.cacheDir. Build it on a machine with internet, or pre-build a prefix and pass it with --auspiceCondaEnv:
        mamba env create -p <prefix> -f envs/nextstrain.yaml
  - A stale mamba index ("... is not installable", "... does not exist"): run  mamba clean --index-cache  and rerun.
  - conda.cacheDir is not writable or is full (default ~/.flutyper/conda): set NXF_CONDA_CACHEDIR to a writable folder.
  - The build timed out (conda.createTimeout, 30 min in nf_pipeline/conf/conda.config): rerun (downloads are cached),
    or pre-build with --auspiceCondaEnv.
  - A killed build left a partial env: delete its env-<hash> folder in conda.cacheDir and rerun.
  - --auspiceCondaEnv points to a prefix without these tools: rebuild it with the command above.
Details: .nextflow.log (search for "conda" / "AuspiceEnvCheck")."""
}
