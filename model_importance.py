from keras import activations
from vis.visualization import visualize_activation, visualize_saliency
from vis.utils import utils

from neuralnet import *
from setup import *


def get_model_data(dir, datadir, model_subdir, exp_name, data_type, label_name, ae=None):
    """

    Args:
        dir: directory to look for data
        datadir: folder in dir with model data
        model_subdir: folder in dir with model experiment
        exp_name: folder in model_subdir with experiment
        data_type: output data type, see args() in neuralnet.py
        label_name:
        ae: path to folder containing autoencoder's gene scaler file

    Returns:

    """
    md = {}
    exp = label_name + '/neuralnets/' + exp_name
    model_fn = os.path.join(dir, model_subdir, exp, 'model_retrained.h5')
    if ae is None:
        scaleFN_in = os.path.join(dir, model_subdir, exp, 'geneTrainScalers.pkl')
    else:
        print('using ae scaler')
        scaleFN_in = os.path.join(ae, 'geneTrainScalers.pkl')
    label_encoder = os.path.join(dir, model_subdir, exp, 'labelencoder.pkl')

    md['nn'] = load_model(model_fn)

    # prep model input
    genes, y_labels, md['ids'] = get_gbm_data(data_dir=datadir,
                                              data_type=data_type,
                                              label_name=label_name)  # read in rg data
    md['gene_ids'] = list(genes.columns.values)
    num_cases, input_dim = genes.shape

    scaler = pickle.load(open(scaleFN_in, 'rb'))  # scale data using model's training scaler
    genes = scaler.transform(genes)
    md['in_min'] = genes.min()
    md['in_max'] = genes.max()
    genes = genes.reshape(len(genes), input_dim, 1)  # reshape input into 3-D shape
    for i in y_labels.iloc[:, 0].unique():
        print(list(y_labels.iloc[:, 0]).count(i), '\t', list(y_labels.iloc[:, 0]).count(i) / y_labels.shape[0], '\t', i)

    if os.path.isfile(label_encoder):
        md['le'] = pickle.load(open(label_encoder, 'rb'))
    md['genes'] = genes
    md['y_labels'] = y_labels
    md['num_cases'] = num_cases
    md['input_dim'] = input_dim
    return md


def get_importance(tab_dir, fig_dir, model_info, label_node, label, x_sample, softmax=False, save=True, fn=None,
                   in_min=0., in_max=1., act_max_weight=1, tv_weight=10., lp_norm_weight=10., max_iter=200,
                   backprop_modifier=None, grad_modifier='absolute'):
    """

    Args:
        tab_dir: dir to save tabular data
        fig_dir: dir to save fig data
        model_info: output of get_model_data()
        label_node:
        label: name of label for output files
        softmax: was softmax used?
        save: save info?
        in_min:
        in_max:
        act_max_weight:
        tv_weight:
        lp_norm_weight:
        max_iter:
        x_sample:
        backprop_modifier:
        grad_modifier:

        for others see visualize_activation() from kera-vis

    Returns: dataframe, rows are genes, col is importance

    """

    # model
    # classes = model_info['le'].classes_
    #
    # print('original class names: ', classes)  # check classes
    # print('transformed class names: ', model_info['le'].transform(classes))
    model = model_info.get('nn')
    layer_idx = utils.find_layer_idx(model, 'preds')  # layer of interest
    if softmax:
        model.layers[layer_idx].activation = activations.linear  # swap softmax with linear
        model = utils.apply_modifications(model)

    filter_idx = label_node  # moe of interest, output node we want to maximize.

    print('getting saliency')
    # print(x_sample.shape)
    sals = visualize_saliency(model=model,
                              layer_idx=layer_idx,
                              filter_indices=filter_idx,
                              seed_input=x_sample,
                              backprop_modifier=backprop_modifier,
                              grad_modifier=grad_modifier)

    sals = sals.reshape(len(model_info.get('gene_ids')))
    sals = pd.DataFrame(data=sals, index=model_info.get('gene_ids'), columns=['activations'])

    if save:
        sals.to_csv(os.path.join(tab_dir, fn + '.csv'))  # save activations
        a = sals.values
        plt.figure()
        plt.hist(a)
        plt.savefig(os.path.join(fig_dir, fn + '.tiff'))

    return sals

