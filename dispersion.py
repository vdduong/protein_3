# dispersion
import math
import re

# distance matrix
# emb = embedding
# neighbors ?

def dispersion_score(u,v,neighbors):
  neighbors_u = neighbors[u]
  neighbors_v = neighbors[v]
  set_uv = set(neighbors_u).intersection(set(neighbors[v]))
  if len(set_uv) = 0.0:
    if u in neighbors_v: return 0.5
    else:
      return 0.
  else:
    score_ = 0.0
    for node_1 in set_uv:
        for node_2 in set_uv:
          if node_1!=node_2 and dispersion_score(node_1,node_2,neighbors) = 0:
            score_ +=1.
    return score_


                      
