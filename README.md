# GLM_example
Example code for how to use GLM to characterize dF/F associated with behavioral events. This method is particularly useful when behavioral events are closely spaced in time, such that the dF/F responses to behavioral events would overlap. For references, see Pinto L and Dan Y, Neuron 2015, and Parker NF, ..., Witten IB, Nature Neuroscience, 2016.

Begin the example by running 'startExample.m'

Figures showing the expected outputs:
- ‘raw-mockdata’ shows how the mock data look like. There are stimuli and choices, each of which induces dF/F changes that are described by kernels. The goal of this exercise is to recover the kernels, i.e., the ground truth for dF/F responses associated with stimuli and response.

- ‘trialavg-mockdata’ shows what happens if we do trial averaging. That is, calculate mean signals using dF/F segments aligned to stimuli or responses. Because responses always occur shortly after stimuli, the dF/F responses would overlap. Therefore this method leads to inaccuracies.

- ‘glmfit-beta’mockdata’ shows results from a GLM fit of the dF/F data. You see the regressors from GLM do not fully recapitulate the ground truth. In particular, because of the regularization, the regressors underestimate the actual kernels. However, compared with trial averaging, this method leads to estimates of the dF/F responses with much less distortion.
