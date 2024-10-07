
for i in range(100):
    idx = str(i).zfill(2)
    file_root = 'grid_{idx}'
    with open(f'sim/{file_root}.py') as f:

        f.writelines([
            'from diss_func import dissolve\n',
            'import numpy as np\n',
            'from tqdm import tqdm\n',
            'from multiprocessing import Pool\n',
            '\n',
            "if __name__ == '__main__':\n",
            f'\tgrid=np.loadtxt("../grids/subgrid_{idx}.csv",delimiter=",")\n',
            '\twith Pool(32) as pool:\n',
            "\t\tresults = list(tqdm(pool.imap(dissolve,grid),total=len(grid), desc= 'Processing'))\n",
            '\tdata=np.hstack((grid,results))\n',
            "\tnp.savetxt(f'grid_sim_{idx}}.csv',data,delimiter=',',fmt='%1.6f')"
        ])

    with open(f'submit/{file_root}.sl') as f:

        f.writelines([
            
        ])