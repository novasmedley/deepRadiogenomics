import time
from model_importance import *
from matplotlib import pyplot as plt

print('gene saliency demo...')
start_time = time.time()
plt.rcParams['figure.figsize'] = (18, 6)

# dirs
d = 'demo'
datadir = 'demo_data'
tab_dir = os.path.join(d, 'nn_gene_importance')  # output dir
if not os.path.exists(tab_dir):
    os.makedirs(tab_dir)

# ae data, the input scaler file
ae_dir = 'ae_cv/autoencoder/neuralnets/'
ae = os.path.join(d,ae_dir, '200_100_50_0_0_tanh_decay_0_drop_0_opt_Nadam_loss_mae_bat_10_eph_2')

# get model
mdir = 'nn_retrain'  # model dir
exp_name = '200_100_50_0_0_tanh_decay_0_drop_0_opt_Nadam_loss_binary_crossentropy_bat_10_eph_2_ae_nl_3_freeze_0'
label = 'f5'
f_info = get_model_data(dir=d,
                        datadir=datadir,
                        model_subdir=mdir,
                        exp_name=exp_name,
                        data_type='vasari',
                        label_name=label,
                        ae=ae)

ppl = [ get_importance(tab_dir=tab_dir,
                       fig_dir=tab_dir,
                       model_info=f_info,
                       label_node=0,
                       label=label,
                       in_min=f_info['in_min'],
                       in_max=f_info['in_max'],
                       x_sample=f_info.get('genes')[i, :, :],
                       save=False,
                       backprop_modifier='guided')
        for i in range(len(f_info['y_labels']))]

a = pd.concat(ppl, axis=1)
a.columns = f_info['y_labels'].index
a.loc['y_label'] = f_info['y_labels'].iloc[:, 0].values  # attach labels (a row)
fn = 'saliency_patientwise_' + label + '_'
a.to_csv(os.path.join(tab_dir, fn + '.csv'))  # save

print(time.time() - start_time)
