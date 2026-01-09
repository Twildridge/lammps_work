Simulation runs comparison
20000 timesteps with stress calculations only
these run times start after the initial soft push step

nodes: 1
ntasks: 4
cpus_per_ntask: 5
run time: 15:18

nodes: 1
ntasks: 4
cpus_per_ntask: 10
run time: 20:02

nodes: 1
ntasks: 8
cpus_per_ntask: 5
run time: 09:19

nodes: 1
ntasks: 10
cpus_per_ntask: 4
run time: 08:37

nodes: 1
ntasks: 20
cpus_per_ntask: 2
run time: 07:42
bridges-2 run time: 10:14

nodes: 1
ntasks: 35
cpus_per_ntask: 1
run time: 06:27
run time: 11:43

nodes: 1
ntasks: 40
cpus_per_ntask: 1
run time: 06:24
run time: 10:36 (12/31)
bridges-2 run time: 10:15

nodes: 1
ntasks: 70
cpus_per_ntask: 1
bridges-2 run time: 07:17

nodes: 1
ntasks: 80
cpus_per_ntask: 1
bridges-2 run time: 07:04

nodes: 1
ntasks: 100
cpus_per_ntask: 1
bridges-2 run time: 05:38

nodes: 1
ntasks: 110
cpus_per_ntask: 1
bridges-2 run time: 04:49

nodes: 1
ntasks: 110
cpus_per_ntask: 1
bridges-2 run time: 05:06

nodes: 1
ntasks: 120
cpus_per_ntask: 1
bridges-2 run time: 04:44

nodes: 1
ntasks: 128
cpus_per_ntask: 1
bridges-2 run time: 08:10

nodes: 1
ntasks: 128
cpus_per_ntask: 1
bridges-2 run time: 08:22






PURE SOLVENT (9,000 atoms, 200000 timesteps):

120 CPUs, no GPU: 0:56
1 CPU, no GPU: 22:33

8 GPUs, 4 threads/GPU, 1 task/node: 2:56
8 GPUs, 4 threads/GPU, 8 tasks/node: 4:59
1 GPU, 4 threads/GPU, 1 task/node: 2:59
1 GPU, 2 threads/GPU, 1 task/node: 2:57
1 GPU, 1 threads/GPU, 5 task/node: >18 minutes

PURE SOLVENT (800,000 atoms, 200000 timesteps):

120 CPUs, no GPU, 5 nodes: 1:02:30

PURE SOLVENT (800,000 atoms, 20000 timesteps):

120 CPUs, no GPU: 30:52
120 CPUs, no GPU, 2 nodes: 16:19
120 CPUs, no GPU, 3 nodes: 11:16
120 CPUs, no GPU, 4 nodes: 08:38
120 CPUs, no GPU, 5 nodes: 07:22

40 CPUs (Pod), no GPU: 43:20
40 CPUs (Pod), no GPU, 2 nodes: 43:49
35 CPUs (Pod), no GPU, 2 nodes: 52:15
40 CPUs (Pod), no GPU, 3 nodes: 28:29

8 GPUs, 4 threads/GPU, 1 tasks/node (541): > 30 mins
1 GPU, 4 threads/GPU, 1 task/node (586): " "
2 GPUs, 4 threads/GPU, 1 task/node (852): " "
4 GPUs, 4 threads/GPU, 1 task/node (594): " "

NOTE: internode communication speed is 200 Gb/s on Bridges-2, but only 100 Gb/s on Pod!!

--gpus=l40s-48:4 gives error! (also l40s-48:8)



SLAB (547,000 mobile atoms, 20000 timesteps):

120 CPUs, no GPU, 1 node: 04:48 
120 CPUs, no GPU, 2 node: 03:03
120 CPUs, no GPU, 3 node: 02:18
120 CPUs, no GPU, 4 node: 01:55


SLAB (547,000 mobile atoms, 10 M timesteps):

120 CPUs, no GPU, 1 node: 10:48:36 




Still waiting on Pod: 
81279 (slab): timed out because was 10 M timesteps with only one node
82665 (slab): also timed out



Walled piston (4 nodes):

double wall: 34:06
single wall with skipped wall–wall interactions (732): 20:36
single wall with skipped wall/piston/support interactions (818): 18:57




