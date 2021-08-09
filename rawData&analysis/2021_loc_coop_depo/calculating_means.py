import numpy as np
import pandas as pd

af_range = np.arange(0.85, 0.96, 0.01)
res = []
for af in af_range:
    af_str = str(af).replace('.', '')[0: 3]
    folder_name = 'init_den=005/results/s3rr7af' + af_str + '/states/'
    para_res = []
    for i in range(0, 30):
        file_name = folder_name + 'lineage{}.csv'.format(i)
        df = pd.read_csv(file_name)
        para_res.append(df.mean().mean())
    print(para_res)
    res.append(para_res)
res_df = pd.DataFrame(res, index=af_range)
res_df.to_csv('init_den=005/means.csv')
