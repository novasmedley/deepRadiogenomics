from model_importance import *

# get args
p = argparse.ArgumentParser()
p.add_argument('--label', help='label name', choices=['f5', 'f6', 'f7', 'f9', 'f10', 'f14', 'subtype'])

args = vars(p.parse_args())  # convert to dictionary

# dirs
d = '/media/nova/Seagate Backup Plus Drive/deepRadiogenomics'
datadir = os.path.join(d, 'gbm_tcga', 'model_ready')
ae_dir = 'gbm_tcga/experiments/gbm_ae_19_04_15_retrain/autoencoder/neuralnets/'  # ae data, the input scaler file
ae_dir = os.path.join(d, ae_dir)
mdirs = {'f5': 'gbm_tcga/experiments/gbm_ae_19_04_15_vasari_retrain',
         'f6': 'gbm_tcga/experiments/gbm_ae_19_04_15_vasari_retrain',
         'f7': 'gbm_tcga/experiments/gbm_ae_19_04_15_vasari_retrain',
         'f9': 'gbm_tcga/experiments/gbm_ae_19_04_15_vasari_retrain',
         'f10': 'gbm_tcga/experiments/gbm_ae_19_04_15_vasari_retrain',
         'f14': 'gbm_tcga/experiments/gbm_ae_19_04_15_vasari_retrain',
         'subtype': 'gbm_tcga/experiments/gbm_veerhak_retrain_18_12_05'
         }  # model data
fig_dir = 'with_ae_paper_figs'
tab_dir = 'with_ae_paper_tabular'

# ae data, the input scaler file
ae_path = os.path.join(ae_dir, '4000_2000_1000_0_0_tanh_decay_0_drop_0_opt_Adadelta_loss_mae_bat_50_eph_500')
pretrains = {'f5': ae_path,
             'f6': ae_path,
             'f7': ae_path,
             'f9': ae_path,
             'f10': ae_path,
             'f14': ae_path,
             'subtype': None
             }
exps = {'f5': '4000_2000_1000_0_0_tanh_decay_0_drop_0.6_opt_Adadelta_loss_binary_crossentropy_bat_10_eph_500_ae_nl_3_freeze_0',
        'f6': '4000_2000_1000_0_0_tanh_decay_0_drop_0.0_opt_Adadelta_loss_binary_crossentropy_bat_10_eph_500_ae_nl_1_freeze_1',
        'f7': '4000_2000_1000_0_0_tanh_decay_0_drop_0.0_opt_Adadelta_loss_binary_crossentropy_bat_10_eph_500_ae_nl_1_freeze_1',
        'f9': '4000_2000_1000_0_0_tanh_decay_0_drop_0.6_opt_Adadelta_loss_binary_crossentropy_bat_10_eph_500_ae_nl_3_freeze_0',
        'f10': '4000_2000_1000_0_0_tanh_decay_0_drop_0.0_opt_Adadelta_loss_binary_crossentropy_bat_10_eph_500_ae_nl_2_freeze_1',
        'f14': '4000_2000_1000_0_0_tanh_decay_0_drop_0.0_opt_Adadelta_loss_binary_crossentropy_bat_10_eph_500_ae_nl_1_freeze_1',
        'subtype': '3000_1500_750_0_0_sigmoid_decay_0_drop_0.4_opt_Nadam_loss_categorical_crossentropy_bat_10_eph_200'
        }
data_type = {'f5': 'vasari',
             'f6': 'vasari',
             'f7': 'vasari',
             'f9': 'vasari',
             'f10': 'vasari',
             'f14': 'vasari',
             'subtype': 'verhaak'
             }

fig_dir = os.path.join(d, 'gbm_tcga', fig_dir, 'nn_gene_masking')
tab_dir = os.path.join(d, 'gbm_tcga', tab_dir, 'nn_gene_masking')

for i in [fig_dir, tab_dir]:
    if not os.path.exists(i):
        os.makedirs(i)

k = args['label']

m_info = get_model_data(dir=d,
                        datadir=datadir,
                        model_subdir=mdirs[k],
                        exp_name=exps[k],
                        data_type=data_type[k],
                        label_name=k,
                        ae=ae_path)
if k == 'subtype':
    num_classes = 4
    classes = m_info['le'].classes_
else:
    num_classes = 1
    classes = ['']

for ci in range(0, num_classes):
    labname = classes[ci]
    ppl = [get_importance(fig_dir=fig_dir,
                          tab_dir=tab_dir,
                          model_info=m_info,
                          label_node=0,
                          label=k,
                          in_min=m_info['in_min'],
                          in_max=m_info['in_max'],
                          x_sample=m_info.get('genes')[i, :, :],
                          save=False,
                          backprop_modifier='guided')
           for i in range(len(m_info['y_labels']))]

    a = pd.concat(ppl, axis=1)
    a.columns = m_info['y_labels'].index
    a.loc['y_label'] = m_info['y_labels'].iloc[:, 0].values  # attach labels (a row)

    fn = 'saliency_patientwise_' + k + '_' + labname
    a.to_csv(os.path.join(tab_dir, fn + '.csv'))  # save

