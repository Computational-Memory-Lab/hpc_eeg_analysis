# Stage 2 ASR OOM on Fir

## Summary

Subject `1001` Stage 2 job `28187150` on Fir reached:

1. ICA completion
2. ICLabel completion
3. artifact IC removal

It then failed in `STEP 12: ASR - Artifact Subspace Reconstruction` with SLURM state `OUT_OF_MEMORY`.

## Symptom

The Stage 2 log showed:

- `Memory too low, increasing it...`
- `Now cleaning data in 16 blocks`
- `slurmstepd: error: Detected 1 oom_kill event`

This happened after the previous channel-location and ICA stdout issues were already fixed.

## Root Cause

The clean_rawdata plugin's ASR code does not respect the SLURM cgroup memory limit by default.

Relevant behavior in the plugin on Fir:

- `clean_asr.m` defaults `maxmem = 64`
- `asr_process.m` then falls back to `hlp_memfree`
- `hlp_memfree` reads node-wide free physical memory via Java

On Fir, node-wide free memory can be much larger than the job's `--mem` allocation. That causes ASR to underestimate how aggressively it should split the work into blocks, which can exceed the job memory limit and trigger an OOM kill.

## Why Subject 1001 Failed but 1002 Succeeded

Subject `1001` entered ASR with a larger workload than successful subject `1002`.

Comparison from Stage 2 logs:

- `1001`: `201` channels, `1,832,280` samples, ASR clean reference `630` seconds
- `1002`: `174` channels, `1,688,754` samples, ASR clean reference `429` seconds

ASR memory use scales approximately with `channels^2 * samples`, so `1001` had a materially larger working set.

## Implemented Fix

Two protections were added:

1. Raise Stage 2 SLURM memory from `128G` to `192G`
2. Pass an explicit `maxmem` value to `clean_asr` from `hpc_set_to_interpol.m`

The explicit `maxmem` logic:

- reads `SLURM_MEM_PER_NODE` when available
- falls back to the Stage 2 default of `192G`
- caps ASR to `50%` of the SLURM limit

This keeps substantial memory headroom for MATLAB, EEGLAB, ICA products, and ASR intermediates while preventing ASR from using node-wide free-memory heuristics.

## Code Changes

- `hpc_set_to_interpol.m`
  - added `resolve_asr_maxmem_mb()`
  - added `parse_memory_limit_mb()`
  - changed the `clean_asr(...)` call to pass the explicit `maxmem` argument
  - logs the detected SLURM memory limit and ASR cap
- `hpc_set_to_interpol.slurm`
  - Stage 2 memory increased to `192G`

## Operational Note

If Stage 2 still OOMs after this fix, the next step should be to lower the ASR cap further, not to remove the cap:

- try `40%` of the SLURM limit instead of `50%`
- or raise Stage 2 memory again if runtime and queue policy allow it

The important rule is: do not let `clean_asr` infer memory from node-wide free RAM on Fir.
