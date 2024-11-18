# Multi-GPU Computing with SELF

SELF builds come with domain-decomposition enabled by default, which allows you to run in parallel using distributed memory parallelism with MPI. When you build SELF with GPU acceleration enabled, native support is provided for multi-GPU platforms. The key requirement is that you are working with a GPU Aware MPI installation (see [GPU Aware MPI with ROCm](https://gpuopen.com/learn/amd-lab-notes/amd-lab-notes-gpu-aware-mpi-readme/) and [GPU Aware MPI with CUDA](https://developer.nvidia.com/blog/introduction-cuda-aware-mpi/)).

## Key Concepts:
1. **Process Affinity**: Bind MPI processes to specific CPUs/cores and GPUs to optimize performance by reducing interconnect overhead.
2. **Sequential Rank Assignment**: Assign MPI ranks such that ranks for a node are packed together before moving to the next node.
3. **Mapping Policies**: Use `mpirun`'s mapping and binding options to control how MPI processes are distributed across nodes and resources.


When deploying SELF on multi-GPU platforms, each MPI rank is assigned to a single GPU. The GPU assignment algorithm is simple and is implemented in the `Init_DomainDecomposition` method defined in the `src/gpu/SELF_DomainDecomposition.f90` module. Each MPI rank queries HIP or CUDA for the number of GPU devices on your server. The device id assigned to the MPI process is set as the modulo of the MPI rank id and the number of devices.

If you are working on clusters of multi-GPU accelerated nodes, this implies that all servers have the same number of GPUs. Additionally, if you are explicitly setting your process affinity, you will want to assign MPI ranks to pack servers sequentially.


## `mpirun` Options for Sequential Affinity:
Use these options when launching your application with `mpirun`:

- `--map-by ppr:N:node`: Places `N` processes per node.
- `--bind-to core` or `--bind-to socket`: Binds processes to specific cores or sockets.
- `--rank-by slot` or `--rank-by node`: Determines how ranks are ordered within the mapping.

### Example Command:
```bash
mpirun --map-by ppr:1:node --bind-to core --rank-by node -np <num_procs> ./your_self_application
```

### Explanation of Options:
1. `--map-by ppr:1:node`: Places 1 MPI process per GPU (assuming one GPU per rank) and packs the ranks sequentially on each node.
2. `--bind-to core`: Binds each MPI process to a core, ensuring proper CPU affinity.
3. `--rank-by node`: Ensures that ranks are assigned sequentially across nodes.
4. `-np <num_procs>`: Specifies the total number of MPI processes.

### Example for Multi-GPU Nodes:
If each node has 4 GPUs and you want 4 MPI processes per node:
```bash
mpirun --map-by ppr:4:node --bind-to core --rank-by node -np <total_procs> ./your_self_application
```

This ensures:
- MPI ranks `0-3` are on Node 1 (bound to GPUs `0-3`).
- MPI ranks `4-7` are on Node 2 (bound to GPUs `0-3`), and so on.

To validate process affinity and report process bindings, both MPI and Slurm provide tools and options to display detailed information about how MPI ranks and processes are mapped to hardware resources (e.g., CPUs, GPUs).

### Validating process bindings with `mpirun`

Most MPI implementations have options to print detailed binding information.

#### OpenMPI
OpenMPI provides options to report process bindings:
- Add `--report-bindings` to your `mpirun` command:
  ```bash
  mpirun --report-bindings --map-by ppr:1:node --bind-to core --rank-by node -np <num_procs> ./your_application
  ```

- Example Output:
  ```
  Rank 0: bound to socket 0[core 0]
  Rank 1: bound to socket 0[core 1]
  Rank 2: bound to socket 0[core 2]
  Rank 3: bound to socket 0[core 3]
  ```

- Use `--display-map` to visualize the mapping of ranks across nodes:
  ```bash
  mpirun --display-map --map-by ppr:1:node --bind-to core ./your_application
  ```

- Example Output:
  ```
  Data for JOB [12345,1] offset 0
  Mapping policy: BYNODE, Binding policy: CORE
  ...
  Node: node01  Ranks: 0, 1, 2, 3
  Node: node02  Ranks: 4, 5, 6, 7
  ```

#### MPICH
MPICH can display binding information using the `MPICH_RANK_REORDER_DISPLAY` environment variable:
- Set the environment variable before running:
  ```bash
  export MPICH_RANK_REORDER_DISPLAY=1
  mpirun -np <num_procs> ./your_application
  ```
- Example Output:
  ```
  Rank 0 running on node01
  Rank 1 running on node01
  ...
  ```

## Launching with Slurm:
For Slurm-managed clusters, the equivalent command is:
```bash
srun --ntasks-per-node=4 --cpus-per-task=<cpus_per_mpi_process> --gpus-per-task=1 --distribution=block:block ./your_application
``` 

This approach ensures proper resource packing and sequential affinity per node.


### Validating process bindings with Slurm 

Slurm provides options to display detailed information about how tasks are distributed and bound.

#### Job Execution Information
- Add the `--cpu-bind` or `--gpu-bind` flag to `srun` to specify binding and display it:
  ```bash
  srun --cpu-bind=verbose --gpus-per-task=1 --ntasks-per-node=4 ./your_application
  ```

- Example Output:
  ```
  srun: CPU Bindings: rank 0 on cores 0-3
  srun: CPU Bindings: rank 1 on cores 4-7
  ```

#### Slurm Environment Variables
- Slurm provides environment variables during job execution, which you can print from within your application:
  ```bash
  getenv("SLURM_TASKS_PER_NODE")   ! Number of tasks per node
  getenv("SLURM_CPUS_ON_NODE")     ! Number of CPUs allocated
  getenv("SLURM_JOB_NODELIST")     ! List of allocated nodes
  getenv("SLURM_LOCALID")          ! Local rank ID
  getenv("SLURM_PROCID")           ! Global rank ID
  getenv("SLURM_CPU_BIND")         ! CPU binding
  ```

- Example Fortran snippet:
  ```fortran
  print *, "SLURM_PROCID:", getenv("SLURM_PROCID")
  print *, "SLURM_CPU_BIND:", getenv("SLURM_CPU_BIND")
  ```

#### Job Allocation Report
Run `scontrol show job <job_id>` to display task and resource binding for a running or completed job:
```bash
scontrol show job <job_id>
```

- Relevant fields in the output include:
  - `TaskAffinity`
  - `CpusPerTask`
  - `TRES` (e.g., GPUs, memory)
  - `Nodes` and `Nodelist`



### System-Specific Considerations:
1. **`mpirun` Implementation**: The specific MPI implementation (e.g., OpenMPI, MPICH) might have slightly different syntax or options. Check your implementation’s documentation.
2. **Resource Manager Integration**: If using a resource manager like Slurm, consider its process binding flags (e.g., `--distribution block:block` or `--ntasks-per-node`).
3. **NUMA domains** : When assing process affinity, you should also consider the latency between CPUs and GPUs on your system. On AMD platforms, you can use `rocm-bandwidth-test` to report on your system's topology; On Nvidia platforms, you can use `nvidia-smi -topo`. Ideally, MPI ranks should be assigned to the NUMA domain closest to their assigned GPU.

