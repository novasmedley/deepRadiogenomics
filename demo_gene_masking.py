'''Take pathways and put through neural networks (individual predictions) to see resultant values'''

import time
from model_perturbation import *

print('gene masking demo...')
start_time = time.time()

# get args
p = argparse.ArgumentParser()
p.add_argument('--label', help='label name', choices=['f5', 'f6', 'f7', 'f9', 'f10', 'f14'])
p.add_argument('--geneset', help='gene set name', choices=['verhaak', 'random_verh', 'puchalski', 'hallmark'])
p.add_argument('--save', help='save predictions? 1 = yes 0 = no', type=int, default=0)
p.add_argument('--gs_name', help='a specific gene set in collection', type=str, default='None')
p.add_argument('--cpus', help='number of cpus', type=int, default=4)
args = vars(p.parse_args())  # convert to dictionary

# dirs
d = 'demo'
datadir = 'demo_data'
tab_dir = os.path.join(d, 'nn_gene_masking')
if not os.path.exists(tab_dir):
        os.makedirs(tab_dir)

# ae data, the input scaler file
ae_dir = 'ae_cv/autoencoder/neuralnets/'
ae_dir = os.path.join(d, ae_dir)
pretrains = {'f5': os.path.join(ae_dir, '200_100_50_0_0_tanh_decay_0_drop_0_opt_Nadam_loss_mae_bat_10_eph_2')}

# model data
mdirs = {'f5': 'nn_retrain'}
exps = {'f5': '200_100_50_0_0_tanh_decay_0_drop_0_opt_Nadam_loss_binary_crossentropy_bat_10_eph_2_ae_nl_3_freeze_0'}
pred_types = {'f5': 'binaryClass'}
y_colnames = {'f5': 'f5.proportion.enhancing'}  # original/input output name
data_type = {'f5': 'vasari'}

# for heatmaps if want to plot probability,
pos_class_names = {'f5': 'More than 1/3'}  # original class name, check
order_bys = {'f5': 'f5 < 1/3 probability'}  # model output name
plot_these_classes = {'f5': ['f5 < 1/3']}  # ordering of classes based on log print out

# gene set data
gs, gs_names, gs_sizes = get_all_gs(datadir=datadir)
gs_keys_list = None
if args['gs_name'] != 'None':
    tab_dir = os.path.join(tab_dir, args['gs_name'])  # get specific gene set in collection
    these_gs = [args['gs_name']]
else:
    these_gs = None  # get all gene sets in the collection
if not (os.path.exists(tab_dir)):
    os.mkdir(tab_dir)

# gene masking
g = args['geneset']
k = args['label']

pert = Perturbation(dir=d,
                    data_subdir=datadir,
                    model_subdir=mdirs[k],
                    ae=ae_dir,
                    pretrain=pretrains[k],
                    exp_name=exps[k],
                    label=k,
                    pos_class_name=pos_class_names[k],
                    pred_type=pred_types[k],
                    order_by=order_bys[k],
                    y_colname=y_colnames[k],
                    plot_these_classes=plot_these_classes[k],
                    data_type=data_type[k],
                    cpus=args['cpus'])

model = load_model(pert.model_fn)

get_all_perts(Pertrubation=pert,
              model=model,
              gs=gs,
              gs_names=these_gs,
              gs_type=g,
              tab_dir=tab_dir,
              save_preds=args['save'])

print(time.time() - start_time)
