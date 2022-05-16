
# Run SpatialDE

import pandas as pd
# Need to install NaiveDE and SpatialDE 
import NaiveDE
import SpatialDE
import numpy as np

# load the count data
counts = pd.read_csv('bc_pattern_zero_10_replicate_1_count.csv', index_col=0)

print(counts)


# load the location data
sample_info = pd.read_csv('bc_pattern_zero_10_replicate_1_loc.csv', index_col=0)

sample_info['total_counts'] = counts.sum(1)

print(sample_info)

import time

st = time.time()

norm_expr = NaiveDE.stabilize(counts.T).T
resid_expr = NaiveDE.regress_out(sample_info, norm_expr.T, 'np.log(total_counts)').T


X = sample_info[['x', 'y']]
results = SpatialDE.run(X, resid_expr)

end = time.time()
print(end-st)

results.to_csv('SpatialDE_result.csv')
