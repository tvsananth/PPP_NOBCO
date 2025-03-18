### Create nodes and edges instances given input parameters
# (x,y) coordinates are X=U(0,1), Y=U(0,1)
# 1/10 nodes are exit nodes, specially created on the boundary of unit square
# capacities c0 are N(50,20)
# free-flow-times t0 are distance-based
# other n-(1/10)n node have SR_demand N(100,10)
# output file instance nodes_n_p.csv and edges_n_p.csv
import sys
import pandas
import numpy
import math

def create_instance_files(n,p,instance_num,save_dir):
# Create non-exit node coords
  num_exit_entry_nodes = int(n / 10)
  x_coords = [numpy.random.uniform(0,1) for x in range(0,n-num_exit_entry_nodes)]
  y_coords = [numpy.random.uniform(0,1) for x in range(0,n-num_exit_entry_nodes)]
  demands = [max(int(numpy.random.normal(100,10)), 0) for x in range(0,n-num_exit_entry_nodes)]
  exits = [0 for x in range(0,n-num_exit_entry_nodes)]
  injureds = [1 for x in range(0,n-num_exit_entry_nodes)]

# Create exit node coords
# Bottom of the unit square
  for x in numpy.arange(0.0,1.0,1.0/num_exit_entry_nodes):
    x_coords.append(x)
  for y in [0] * num_exit_entry_nodes:
    y_coords.append(y)
    demands.append(0)
    exits.append(1)
    injureds.append(0)

# Top of the unit square
  for x in numpy.arange(0.0,1.0,1.0/num_exit_entry_nodes):
    x_coords.append(x)
  for y in [1] * num_exit_entry_nodes:
    y_coords.append(y)
    demands.append(0)
    exits.append(1)
    injureds.append(0)

# Create edges with probability p, scaled for distance
  froms = []
  tos = []
  t0s = []
  c0s = []
  lanes = []
  initial_SRs = []
  initial_FRs = []
  for i in range(len(x_coords)-1):
    for j in range(i+1,len(x_coords)):
      distance = math.sqrt((x_coords[j]-x_coords[i])**2 + (y_coords[j]-y_coords[i])**2)
      if distance > 0.8:
        distance_factor = -0.5
      elif distance > 0.6:
        distance_factor = -0.3
      elif distance > 0.4:
        distance_factor = -0.2
      elif distance > 0.2:
        distance_factor = -0.1
      else:
        distance_factor = 0

      should_create_edge = numpy.random.uniform(0,1)
      if should_create_edge < (p + distance_factor):
        froms.append(i)
        tos.append(j)
        t0s.append(int(distance*60))
        capacity = max(int(numpy.random.normal(50,20)),0)
        c0s.append(capacity)
        lanes.append(2)
        initial_SRs.append(-1)
        initial_FRs.append(-1)

        # and reverse edges
        froms.append(j)
        tos.append(i)
        t0s.append(int(distance*60))
        c0s.append(capacity)
        lanes.append(2)
        initial_SRs.append(-1)
        initial_FRs.append(-1)

# Create edges file
  edges_dict = {'from' : froms,
                'to' : tos,
                't0' : t0s,
                'c0' : c0s,
                'lanes' : lanes,
                'initial_SR' : initial_SRs,
                'initial_FR' : initial_FRs}

  edges_df = pandas.DataFrame(edges_dict)
  edges_df.to_csv(f"{save_dir}/instances/random_instances/edges_{n}_{p}_{instance_num}.csv",index=False)

# Create nodes file
  nodes_dict = {'x' : x_coords,
                'y' : y_coords,
                'SR_demand' : demands,
                'exits' : exits,
                'injured' : injureds}
  nodes_df = pandas.DataFrame(nodes_dict)
  nodes_df.to_csv(f"{save_dir}/instances/random_instances/nodes_{n}_{p}_{instance_num}.csv",index=False)

  print(n, p, instance_num, len(edges_dict['c0']))


save_dir = '../data/instances/random_instances'
"""
for instance_num in range(10):
  for n in [10]:
    for p in [0.75]:
      create_instance_files(n,p,instance_num,save_dir)
"""
