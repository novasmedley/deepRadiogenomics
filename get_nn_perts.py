'''Take pathways and put through neural networks (individual predictions) to see resultant values'''

from model_perturbation import *

# get args
p = argparse.ArgumentParser()
p.add_argument('--label', help='label name', choices=['f5', 'f6', 'f7', 'f9', 'f10', 'f14', 'subtype', 'all'])
p.add_argument('--geneset', help='gene set collection name', choices=['verhaak', 'random_verh',
                                                           'puchalski', 'hallmark',
                                                           'canonical_paths', 'chem_gene_perts',
                                                           'reactome', 'biocarta', 'kegg',
                                                           'chromosome', 'motif',
                                                           'computational', 'GO',
                                                           'onco_sig', 'immuno_sig',
                                                           'single_gene', 'test'])
p.add_argument('--cpus', help='number of cpus', type=int, default=4)
p.add_argument('--save', help='save predictions? 1 = yes 0 = no', type=int, default=0)
p.add_argument('--gs_name', help=' a specific gene set in collection', type=str, default='None')

args = vars(p.parse_args())  # convert to dictionary

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

fig_dir = os.path.join(d, 'gbm_tcga', fig_dir, 'nn_gene_masking')
tab_dir = os.path.join(d, 'gbm_tcga', tab_dir, 'nn_gene_masking')

gs, gs_names, gs_sizes = get_all_gs(datadir=datadir)  # gene set data
gs_keys_list = None
if args['gs_name'] != 'None':
    tab_dir = os.path.join(tab_dir, args['gs_name'])  # get specific gene set in collection
    these_gs = [args['gs_name']]
else:
    these_gs = None  # get all gene sets in the collection
if not (os.path.exists(tab_dir)):
    os.mkdir(tab_dir)

pos_class_names = {'f5': 'More than 1/3', # for heatmaps if want to plot probability
                   'f6': 'More than 1/3',
                   'f7': 'More than 1/3',
                   'f9': 'non-focal',
                   'f10': 'infiltrative',
                   'f14': 'More than 1/3',
                   'subtype': 'Mesenchymal'
                   }  # an original class name, check
pred_types = {'f5': 'binaryClass',
              'f6': 'binaryClass',
              'f7': 'binaryClass',
              'f9': 'binaryClass',
              'f10': 'binaryClass',
              'f14': 'binaryClass',
              'subtype': 'multiClass'
              }
order_bys = {'f5': 'f5 < 1/3 probability',  # for heatmaps if want to plot probability
             'f6': 'f6 < 1/3 probability',
             'f7': 'f7 < 1/3 probability',
             'f9': 'f9 non-focal probability',
             'f10': 'f10 infiltrative probability',
             'f14': 'f14 < 1/3 probability',
             'subtype': 'MES probability'
             }  # model output name
y_colnames = {'f5': 'f5.proportion.enhancing',
              'f6': 'f6.proportion.ncet',
              'f7': 'f7.proportion.necrosis',
              'f9': 'f9.multifocal.or.multicentric',
              'f10': 'f10.t1.flair.ratio',
              'f14': 'f14.proportion.edema',
              'subtype': 'subtype'
              }  # original/input output name
plot_these_classes = {'f5': ['f5 < 1/3'],  # for heatmaps if want to plot probability
                      'f6': ['f6 < 1/3'],
                      'f7': ['f7 < 1/3'],
                      'f9': ['f9 non-focal'],
                      'f10': ['f10 infiltrative'],
                      'f14': ['f14 < 1/3'],
                      'subtype': ['CL', 'MES', 'NL', 'PN']
                      }  # ordering of classes based on log print out
data_type = {'f5': 'vasari',
             'f6': 'vasari',
             'f7': 'vasari',
             'f9': 'vasari',
             'f10': 'vasari',
             'f14': 'vasari',
             'subtype': 'verhaak'
             }

g = args['geneset']

if args['label'] == 'all':
    h = list(exps.keys())
else:
    h = [args['label']]

for k in h:
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
    get_all_perts(Pertrubation=pert, model=model, gs=gs, gs_names=these_gs, gs_type=g, tab_dir=tab_dir,
                  save_preds=args['save'])
    del model
    K.clear_session()
