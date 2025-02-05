# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0. To view a copy of
# this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# To request access to the AlphaFold 3 model parameters, follow the process set
# out at https://github.com/google-deepmind/alphafold3. You may only use these
# if received directly from Google. Use is subject to terms of use available at
# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

"""AlphaFold 3 structure prediction script.

AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0. To view a copy of
this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/

To request access to the AlphaFold 3 model parameters, follow the process set
out at https://github.com/google-deepmind/alphafold3. You may only use these
if received directly from Google. Use is subject to terms of use available at
https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md
"""

from collections.abc import Callable, Sequence
import csv
import dataclasses
import datetime
import functools
import multiprocessing
import os
import pathlib
import shutil
import string
import textwrap
import time
import typing
from typing import overload

from absl import app
from absl import flags
from alphafold3.common import folding_input
from alphafold3.common import resources
from alphafold3.constants import chemical_components
import alphafold3.cpp
from alphafold3.data import featurisation
from alphafold3.data import pipeline
from alphafold3.jax.attention import attention
from alphafold3.model import features
from alphafold3.model import model
from alphafold3.model import params
from alphafold3.model import post_processing
from alphafold3.model.components import utils
import haiku as hk
import jax
from jax import numpy as jnp
import numpy as np


_HOME_DIR = pathlib.Path(os.environ.get('HOME'))
_DEFAULT_MODEL_DIR = _HOME_DIR / 'models'
_DEFAULT_DB_DIR = _HOME_DIR / 'public_databases'


# Input and output paths.
_JSON_PATH = flags.DEFINE_string(
    'json_path',
    None,
    'Path to the input JSON file.',
)
_INPUT_DIR = flags.DEFINE_string(
    'input_dir',
    None,
    'Path to the directory containing input JSON files.',
)
_OUTPUT_DIR = flags.DEFINE_string(
    'output_dir',
    None,
    'Path to a directory where the results will be saved.',
)
MODEL_DIR = flags.DEFINE_string(
    'model_dir',
    _DEFAULT_MODEL_DIR.as_posix(),
    'Path to the model to use for inference.',
)

# Control which stages to run.
_RUN_DATA_PIPELINE = flags.DEFINE_bool(
    'run_data_pipeline',
    True,
    'Whether to run the data pipeline on the fold inputs.',
)
_RUN_INFERENCE = flags.DEFINE_bool(
    'run_inference',
    True,
    'Whether to run inference on the fold inputs.',
)

# Binary paths.
_JACKHMMER_BINARY_PATH = flags.DEFINE_string(
    'jackhmmer_binary_path',
    shutil.which('jackhmmer'),
    'Path to the Jackhmmer binary.',
)
_NHMMER_BINARY_PATH = flags.DEFINE_string(
    'nhmmer_binary_path',
    shutil.which('nhmmer'),
    'Path to the Nhmmer binary.',
)
_HMMALIGN_BINARY_PATH = flags.DEFINE_string(
    'hmmalign_binary_path',
    shutil.which('hmmalign'),
    'Path to the Hmmalign binary.',
)
_HMMSEARCH_BINARY_PATH = flags.DEFINE_string(
    'hmmsearch_binary_path',
    shutil.which('hmmsearch'),
    'Path to the Hmmsearch binary.',
)
_HMMBUILD_BINARY_PATH = flags.DEFINE_string(
    'hmmbuild_binary_path',
    shutil.which('hmmbuild'),
    'Path to the Hmmbuild binary.',
)

# Database paths.
DB_DIR = flags.DEFINE_multi_string(
    'db_dir',
    (_DEFAULT_DB_DIR.as_posix(),),
    'Path to the directory containing the databases. Can be specified multiple'
    ' times to search multiple directories in order.',
)

_SMALL_BFD_DATABASE_PATH = flags.DEFINE_string(
    'small_bfd_database_path',
    '${DB_DIR}/bfd-first_non_consensus_sequences.fasta',
    'Small BFD database path, used for protein MSA search.',
)
_MGNIFY_DATABASE_PATH = flags.DEFINE_string(
    'mgnify_database_path',
    '${DB_DIR}/mgy_clusters_2022_05.fa',
    'Mgnify database path, used for protein MSA search.',
)
_UNIPROT_CLUSTER_ANNOT_DATABASE_PATH = flags.DEFINE_string(
    'uniprot_cluster_annot_database_path',
    '${DB_DIR}/uniprot_all_2021_04.fa',
    'UniProt database path, used for protein paired MSA search.',
)
_UNIREF90_DATABASE_PATH = flags.DEFINE_string(
    'uniref90_database_path',
    '${DB_DIR}/uniref90_2022_05.fa',
    'UniRef90 database path, used for MSA search. The MSA obtained by '
    'searching it is used to construct the profile for template search.',
)
_NTRNA_DATABASE_PATH = flags.DEFINE_string(
    'ntrna_database_path',
    '${DB_DIR}/nt_rna_2023_02_23_clust_seq_id_90_cov_80_rep_seq.fasta',
    'NT-RNA database path, used for RNA MSA search.',
)
_RFAM_DATABASE_PATH = flags.DEFINE_string(
    'rfam_database_path',
    '${DB_DIR}/rfam_14_9_clust_seq_id_90_cov_80_rep_seq.fasta',
    'Rfam database path, used for RNA MSA search.',
)
_RNA_CENTRAL_DATABASE_PATH = flags.DEFINE_string(
    'rna_central_database_path',
    '${DB_DIR}/rnacentral_active_seq_id_90_cov_80_linclust.fasta',
    'RNAcentral database path, used for RNA MSA search.',
)
_PDB_DATABASE_PATH = flags.DEFINE_string(
    'pdb_database_path',
    '${DB_DIR}/mmcif_files',
    'PDB database directory with mmCIF files path, used for template search.',
)
_SEQRES_DATABASE_PATH = flags.DEFINE_string(
    'seqres_database_path',
    '${DB_DIR}/pdb_seqres_2022_09_28.fasta',
    'PDB sequence database path, used for template search.',
)

# Number of CPUs to use for MSA tools.
_JACKHMMER_N_CPU = flags.DEFINE_integer(
    'jackhmmer_n_cpu',
    min(multiprocessing.cpu_count(), 8),
    'Number of CPUs to use for Jackhmmer. Default to min(cpu_count, 8). Going'
    ' beyond 8 CPUs provides very little additional speedup.',
)
_NHMMER_N_CPU = flags.DEFINE_integer(
    'nhmmer_n_cpu',
    min(multiprocessing.cpu_count(), 8),
    'Number of CPUs to use for Nhmmer. Default to min(cpu_count, 8). Going'
    ' beyond 8 CPUs provides very little additional speedup.',
)

# Template search configuration.
_MAX_TEMPLATE_DATE = flags.DEFINE_string(
    'max_template_date',
    '2021-09-30',  # By default, use the date from the AlphaFold 3 paper.
    'Maximum template release date to consider. Format: YYYY-MM-DD. All '
    'templates released after this date will be ignored.',
)

_CONFORMER_MAX_ITERATIONS = flags.DEFINE_integer(
    'conformer_max_iterations',
    None,  # Default to RDKit default parameters value.
    'Optional override for maximum number of iterations to run for RDKit '
    'conformer search.',
)

# JAX inference performance tuning.
_JAX_COMPILATION_CACHE_DIR = flags.DEFINE_string(
    'jax_compilation_cache_dir',
    None,
    'Path to a directory for the JAX compilation cache.',
)
_GPU_DEVICE = flags.DEFINE_integer(
    'gpu_device',
    0,
    'Optional override for the GPU device to use for inference. Defaults to the'
    ' 1st GPU on the system. Useful on multi-GPU systems to pin each run to a'
    ' specific GPU.',
)
_BUCKETS = flags.DEFINE_list(
    'buckets',
    # pyformat: disable
    ['256', '512', '768', '1024', '1280', '1536', '2048', '2560', '3072',
     '3584', '4096', '4608', '5120'],
    # pyformat: enable
    'Strictly increasing order of token sizes for which to cache compilations.'
    ' For any input with more tokens than the largest bucket size, a new bucket'
    ' is created for exactly that number of tokens.',
)
_FLASH_ATTENTION_IMPLEMENTATION = flags.DEFINE_enum(
    'flash_attention_implementation',
    default='triton',
    enum_values=['triton', 'cudnn', 'xla'],
    help=(
        "Flash attention implementation to use. 'triton' and 'cudnn' uses a"
        ' Triton and cuDNN flash attention implementation, respectively. The'
        ' Triton kernel is fastest and has been tested more thoroughly. The'
        " Triton and cuDNN kernels require Ampere GPUs or later. 'xla' uses an"
        ' XLA attention implementation (no flash attention) and is portable'
        ' across GPU devices.'
    ),
)
_NUM_RECYCLES = flags.DEFINE_integer(
    'num_recycles',
    10,
    'Number of recycles to use during inference.',
    lower_bound=1,
)
_NUM_DIFFUSION_SAMPLES = flags.DEFINE_integer(
    'num_diffusion_samples',
    5,
    'Number of diffusion samples to generate.',
    lower_bound=1,
)
_NUM_SEEDS = flags.DEFINE_integer(
    'num_seeds',
    None,
    'Number of seeds to use for inference. If set, only a single seed must be'
    ' provided in the input JSON. AlphaFold 3 will then generate random seeds'
    ' in sequence, starting from the single seed specified in the input JSON.'
    ' The full input JSON produced by AlphaFold 3 will include the generated'
    ' random seeds. If not set, AlphaFold 3 will use the seeds as provided in'
    ' the input JSON.',
    lower_bound=1,
)

# Output controls.
_SAVE_EMBEDDINGS = flags.DEFINE_bool(
    'save_embeddings',
    False,
    'Whether to save the final trunk single and pair embeddings in the output.',
)

_SAVE_DISTOGRAM = flags.DEFINE_bool(
    'save_distogram',
    False,
    'Whether to save the final distogram in the output.',
)

def make_model_config(
    *,
    flash_attention_implementation: attention.Implementation = 'triton',
    num_diffusion_samples: int = 5,
    num_recycles: int = 10,
    return_embeddings: bool = False,
    return_distogram: bool = False,
) -> model.Model.Config:
  """Returns a model config with some defaults overridden."""
  config = model.Model.Config()
  config.global_config.flash_attention_implementation = (
      flash_attention_implementation
  )
  config.heads.diffusion.eval.num_samples = num_diffusion_samples
  config.num_recycles = num_recycles
  config.return_embeddings = return_embeddings
  config.return_distogram = return_distogram
  return config


class ModelRunner:
  """Helper class to run structure prediction stages."""

  def __init__(
      self,
      config: model.Model.Config,
      device: jax.Device,
      model_dir: pathlib.Path,
  ):
    self._model_config = config
    self._device = device
    self._model_dir = model_dir

  @functools.cached_property
  def model_params(self) -> hk.Params:
    """Loads model parameters from the model directory."""
    return params.get_model_haiku_params(model_dir=self._model_dir)

  @functools.cached_property
  def _model(
      self,
  ) -> Callable[[jnp.ndarray, features.BatchDict], model.ModelResult]:
    """Loads model parameters and returns a jitted model forward pass."""

    @hk.transform
    def forward_fn(batch):
      return model.Model(self._model_config)(batch)

    return functools.partial(
        jax.jit(forward_fn.apply, device=self._device), self.model_params
    )

  def run_inference(
      self, featurised_example: features.BatchDict, rng_key: jnp.ndarray
  ) -> model.ModelResult:
    """Computes a forward pass of the model on a featurised example."""
    featurised_example = jax.device_put(
        jax.tree_util.tree_map(
            jnp.asarray, utils.remove_invalidly_typed_feats(featurised_example)
        ),
        self._device,
    )

    result = self._model(rng_key, featurised_example)
    result = jax.tree.map(np.asarray, result)
    result = jax.tree.map(
        lambda x: x.astype(jnp.float32) if x.dtype == jnp.bfloat16 else x,
        result,
    )
    result = dict(result)
    identifier = self.model_params['__meta__']['__identifier__'].tobytes()
    result['__identifier__'] = identifier
    return result

  def extract_structures(
      self,
      batch: features.BatchDict,
      result: model.ModelResult,
      target_name: str,
  ) -> list[model.InferenceResult]:
    """Generates structures from model outputs."""
    return list(
        model.Model.get_inference_result(
            batch=batch, result=result, target_name=target_name
        )
    )

  def extract_embeddings(
      self,
      result: model.ModelResult,
  ) -> dict[str, np.ndarray] | None:
    """Extracts embeddings from model outputs."""
    embeddings = {}
    if 'single_embeddings' in result:
      embeddings['single_embeddings'] = result['single_embeddings']
    if 'pair_embeddings' in result:
      embeddings['pair_embeddings'] = result['pair_embeddings']
    return embeddings or None

  def extract_distogram(
      self,
      result: model.ModelResult,
      num_tokens: int,
  ) -> np.ndarray | None:
    """Extracts distogram from model outputs."""
    if 'distogram' not in result['distogram']:
      return None
    distogram = result['distogram']['distogram'][:num_tokens, :num_tokens, :]
    return distogram

@dataclasses.dataclass(frozen=True, slots=True, kw_only=True)
class ResultsForSeed:
  """Stores the inference results (diffusion samples) for a single seed.

  Attributes:
    seed: The seed used to generate the samples.
    inference_results: The inference results, one per sample.
    full_fold_input: The fold input that must also include the results of
      running the data pipeline - MSA and templates.
    embeddings: The final trunk single and pair embeddings, if requested.
  """

  seed: int
  inference_results: Sequence[model.InferenceResult]
  full_fold_input: folding_input.Input
  embeddings: dict[str, np.ndarray] | None = None
  distogram: np.ndarray | None = None


def predict_structure(
    fold_input: folding_input.Input,
    model_runner: ModelRunner,
    buckets: Sequence[int] | None = None,
    conformer_max_iterations: int | None = None,
) -> Sequence[ResultsForSeed]:
  """Runs the full inference pipeline to predict structures for each seed."""

  print(f'Featurising data with {len(fold_input.rng_seeds)} seed(s)...')
  featurisation_start_time = time.time()
  ccd = chemical_components.cached_ccd(user_ccd=fold_input.user_ccd)
  featurised_examples = featurisation.featurise_input(
      fold_input=fold_input,
      buckets=buckets,
      ccd=ccd,
      verbose=True,
      conformer_max_iterations=conformer_max_iterations,
  )
  print(
      f'Featurising data with {len(fold_input.rng_seeds)} seed(s) took'
      f' {time.time() - featurisation_start_time:.2f} seconds.'
  )
  print(
      'Running model inference and extracting output structure samples with'
      f' {len(fold_input.rng_seeds)} seed(s)...'
  )
  all_inference_start_time = time.time()
  all_inference_results = []
  for seed, example in zip(fold_input.rng_seeds, featurised_examples):
    print(f'Running model inference with seed {seed}...')
    inference_start_time = time.time()
    rng_key = jax.random.PRNGKey(seed)
    result = model_runner.run_inference(example, rng_key)
    print(
        f'Running model inference with seed {seed} took'
        f' {time.time() - inference_start_time:.2f} seconds.'
    )
    print(f'Extracting output structure samples with seed {seed}...')
    extract_structures = time.time()
    inference_results = model_runner.extract_structures(
        batch=example, result=result, target_name=fold_input.name
    )
    print(
        f'Extracting {len(inference_results)} output structure samples with'
        f' seed {seed} took {time.time() - extract_structures:.2f} seconds.'
    )
    
    embeddings = model_runner.extract_embeddings(result)
    distogram = model_runner.extract_distogram(result, example['seq_length'])
    
    all_inference_results.append(
        ResultsForSeed(
            seed=seed,
            inference_results=inference_results,
            full_fold_input=fold_input,
            embeddings=embeddings,
            distogram=distogram,
        )
    )
  print(
      'Running model inference and extracting output structures with'
      f' {len(fold_input.rng_seeds)} seed(s) took'
      f' {time.time() - all_inference_start_time:.2f} seconds.'
  )
  return all_inference_results


def write_fold_input_json(
    fold_input: folding_input.Input,
    output_dir: os.PathLike[str] | str,
) -> None:
  """Writes the input JSON to the output directory."""
  os.makedirs(output_dir, exist_ok=True)
  path = os.path.join(output_dir, f'{fold_input.sanitised_name()}_data.json')
  print(f'Writing model input JSON to {path}')
  with open(path, 'wt') as f:
    f.write(fold_input.to_json())


def write_outputs(
    all_inference_results: Sequence[ResultsForSeed],
    output_dir: os.PathLike[str] | str,
    job_name: str,
) -> None:
  """Writes outputs to the specified output directory."""
  ranking_scores = []
  max_ranking_score = None
  max_ranking_result = None

  output_terms = (
      pathlib.Path(alphafold3.cpp.__file__).parent / 'OUTPUT_TERMS_OF_USE.md'
  ).read_text()

  os.makedirs(output_dir, exist_ok=True)
  for results_for_seed in all_inference_results:
    seed = results_for_seed.seed
    for sample_idx, result in enumerate(results_for_seed.inference_results):
      sample_dir = os.path.join(output_dir, f'seed-{seed}_sample-{sample_idx}')
      os.makedirs(sample_dir, exist_ok=True)
      post_processing.write_output(
          inference_result=result, output_dir=sample_dir
      )
      ranking_score = float(result.metadata['ranking_score'])
      ranking_scores.append((seed, sample_idx, ranking_score))
      if max_ranking_score is None or ranking_score > max_ranking_score:
        max_ranking_score = ranking_score
        max_ranking_result = result

    if embeddings := results_for_seed.embeddings:
      embeddings_dir = os.path.join(output_dir, f'seed-{seed}_embeddings')
      os.makedirs(embeddings_dir, exist_ok=True)
      post_processing.write_embeddings(
          embeddings=embeddings, output_dir=embeddings_dir
      )

    if (distogram := results_for_seed.distogram) is not None:
      print(distogram.shape)
      post_processing.write_distogram(
          distogram=distogram, output_path=os.path.join(output_dir, f'distogram_seed-{seed}.npz')
      )

  if max_ranking_result is not None:  # True iff ranking_scores non-empty.
    post_processing.write_output(
        inference_result=max_ranking_result,
        output_dir=output_dir,
        # The output terms of use are the same for all seeds/samples.
        terms_of_use=output_terms,
        name=job_name,
    )
    # Save csv of ranking scores with seeds and sample indices, to allow easier
    # comparison of ranking scores across different runs.
    with open(os.path.join(output_dir, 'ranking_scores.csv'), 'wt') as f:
      writer = csv.writer(f)
      writer.writerow(['seed', 'sample', 'ranking_score'])
      writer.writerows(ranking_scores)


@overload
def process_fold_input(
    fold_input: folding_input.Input,
    data_pipeline_config: pipeline.DataPipelineConfig | None,
    model_runner: None,
    output_dir: os.PathLike[str] | str,
    buckets: Sequence[int] | None = None,
) -> folding_input.Input:
  ...


@overload
def process_fold_input(
    fold_input: folding_input.Input,
    data_pipeline_config: pipeline.DataPipelineConfig | None,
    model_runner: ModelRunner,
    output_dir: os.PathLike[str] | str,
    buckets: Sequence[int] | None = None,
) -> Sequence[ResultsForSeed]:
  ...


def replace_db_dir(path_with_db_dir: str, db_dirs: Sequence[str]) -> str:
  """Replaces the DB_DIR placeholder in a path with the given DB_DIR."""
  template = string.Template(path_with_db_dir)
  if 'DB_DIR' in template.get_identifiers():
    for db_dir in db_dirs:
      path = template.substitute(DB_DIR=db_dir)
      if os.path.exists(path):
        return path
    raise FileNotFoundError(
        f'{path_with_db_dir} with ${{DB_DIR}} not found in any of {db_dirs}.'
    )
  if not os.path.exists(path_with_db_dir):
    raise FileNotFoundError(f'{path_with_db_dir} does not exist.')
  return path_with_db_dir


def process_fold_input(
    fold_input: folding_input.Input,
    data_pipeline_config: pipeline.DataPipelineConfig | None,
    model_runner: ModelRunner | None,
    output_dir: os.PathLike[str] | str,
    buckets: Sequence[int] | None = None,
    conformer_max_iterations: int | None = None,
) -> folding_input.Input | Sequence[ResultsForSeed]:
  """Runs data pipeline and/or inference on a single fold input.

  Args:
    fold_input: Fold input to process.
    data_pipeline_config: Data pipeline config to use. If None, skip the data
      pipeline.
    model_runner: Model runner to use. If None, skip inference.
    output_dir: Output directory to write to.
    buckets: Bucket sizes to pad the data to, to avoid excessive re-compilation
      of the model. If None, calculate the appropriate bucket size from the
      number of tokens. If not None, must be a sequence of at least one integer,
      in strictly increasing order. Will raise an error if the number of tokens
      is more than the largest bucket size.
    conformer_max_iterations: Optional override for maximum number of iterations
      to run for RDKit conformer search.

  Returns:
    The processed fold input, or the inference results for each seed.

  Raises:
    ValueError: If the fold input has no chains.
  """
  print(f'\nRunning fold job {fold_input.name}...')

  if not fold_input.chains:
    raise ValueError('Fold input has no chains.')

  if os.path.exists(output_dir) and os.listdir(output_dir):
    new_output_dir = (
        f'{output_dir}_{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}'
    )
    print(
        f'Output will be written in {new_output_dir} since {output_dir} is'
        ' non-empty.'
    )
    output_dir = new_output_dir
  else:
    print(f'Output will be written in {output_dir}')

  if data_pipeline_config is None:
    print('Skipping data pipeline...')
  else:
    print('Running data pipeline...')
    fold_input = pipeline.DataPipeline(data_pipeline_config).process(fold_input)

  write_fold_input_json(fold_input, output_dir)
  if model_runner is None:
    print('Skipping model inference...')
    output = fold_input
  else:
    print(
        f'Predicting 3D structure for {fold_input.name} with'
        f' {len(fold_input.rng_seeds)} seed(s)...'
    )
    all_inference_results = predict_structure(
        fold_input=fold_input,
        model_runner=model_runner,
        buckets=buckets,
        conformer_max_iterations=conformer_max_iterations,
    )
    print(f'Writing outputs with {len(fold_input.rng_seeds)} seed(s)...')
    write_outputs(
        all_inference_results=all_inference_results,
        output_dir=output_dir,
        job_name=fold_input.sanitised_name(),
    )
    output = all_inference_results

  print(f'Fold job {fold_input.name} done.\n')
  return output


def main(_):
  if _JAX_COMPILATION_CACHE_DIR.value is not None:
    jax.config.update(
        'jax_compilation_cache_dir', _JAX_COMPILATION_CACHE_DIR.value
    )

  if _JSON_PATH.value is None == _INPUT_DIR.value is None:
    raise ValueError(
        'Exactly one of --json_path or --input_dir must be specified.'
    )

  if not _RUN_INFERENCE.value and not _RUN_DATA_PIPELINE.value:
    raise ValueError(
        'At least one of --run_inference or --run_data_pipeline must be'
        ' set to true.'
    )

  if _INPUT_DIR.value is not None:
    fold_inputs = folding_input.load_fold_inputs_from_dir(
        pathlib.Path(_INPUT_DIR.value)
    )
  elif _JSON_PATH.value is not None:
    fold_inputs = folding_input.load_fold_inputs_from_path(
        pathlib.Path(_JSON_PATH.value)
    )
  else:
    raise AssertionError(
        'Exactly one of --json_path or --input_dir must be specified.'
    )

  # Make sure we can create the output directory before running anything.
  try:
    os.makedirs(_OUTPUT_DIR.value, exist_ok=True)
  except OSError as e:
    print(f'Failed to create output directory {_OUTPUT_DIR.value}: {e}')
    raise

  if _RUN_INFERENCE.value:
    # Fail early on incompatible devices, but only if we're running inference.
    gpu_devices = jax.local_devices(backend='gpu')
    if gpu_devices:
      compute_capability = float(
          gpu_devices[_GPU_DEVICE.value].compute_capability
      )
      if compute_capability < 6.0:
        raise ValueError(
            'AlphaFold 3 requires at least GPU compute capability 6.0 (see'
            ' https://developer.nvidia.com/cuda-gpus).'
        )
      elif 7.0 <= compute_capability < 8.0:
        xla_flags = os.environ.get('XLA_FLAGS')
        required_flag = '--xla_disable_hlo_passes=custom-kernel-fusion-rewriter'
        if not xla_flags or required_flag not in xla_flags:
          raise ValueError(
              'For devices with GPU compute capability 7.x (see'
              ' https://developer.nvidia.com/cuda-gpus) the ENV XLA_FLAGS must'
              f' include "{required_flag}".'
          )
        if _FLASH_ATTENTION_IMPLEMENTATION.value != 'xla':
          raise ValueError(
              'For devices with GPU compute capability 7.x (see'
              ' https://developer.nvidia.com/cuda-gpus) the'
              ' --flash_attention_implementation must be set to "xla".'
          )

  notice = textwrap.wrap(
      'Running AlphaFold 3. Please note that standard AlphaFold 3 model'
      ' parameters are only available under terms of use provided at'
      ' https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md.'
      ' If you do not agree to these terms and are using AlphaFold 3 derived'
      ' model parameters, cancel execution of AlphaFold 3 inference with'
      ' CTRL-C, and do not use the model parameters.',
      break_long_words=False,
      break_on_hyphens=False,
      width=80,
  )
  print('\n' + '\n'.join(notice) + '\n')

  if _RUN_DATA_PIPELINE.value:
    expand_path = lambda x: replace_db_dir(x, DB_DIR.value)
    max_template_date = datetime.date.fromisoformat(_MAX_TEMPLATE_DATE.value)
    data_pipeline_config = pipeline.DataPipelineConfig(
        jackhmmer_binary_path=_JACKHMMER_BINARY_PATH.value,
        nhmmer_binary_path=_NHMMER_BINARY_PATH.value,
        hmmalign_binary_path=_HMMALIGN_BINARY_PATH.value,
        hmmsearch_binary_path=_HMMSEARCH_BINARY_PATH.value,
        hmmbuild_binary_path=_HMMBUILD_BINARY_PATH.value,
        small_bfd_database_path=expand_path(_SMALL_BFD_DATABASE_PATH.value),
        mgnify_database_path=expand_path(_MGNIFY_DATABASE_PATH.value),
        uniprot_cluster_annot_database_path=expand_path(
            _UNIPROT_CLUSTER_ANNOT_DATABASE_PATH.value
        ),
        uniref90_database_path=expand_path(_UNIREF90_DATABASE_PATH.value),
        ntrna_database_path=expand_path(_NTRNA_DATABASE_PATH.value),
        rfam_database_path=expand_path(_RFAM_DATABASE_PATH.value),
        rna_central_database_path=expand_path(_RNA_CENTRAL_DATABASE_PATH.value),
        pdb_database_path=expand_path(_PDB_DATABASE_PATH.value),
        seqres_database_path=expand_path(_SEQRES_DATABASE_PATH.value),
        jackhmmer_n_cpu=_JACKHMMER_N_CPU.value,
        nhmmer_n_cpu=_NHMMER_N_CPU.value,
        max_template_date=max_template_date,
    )
  else:
    data_pipeline_config = None

  if _RUN_INFERENCE.value:
    devices = jax.local_devices(backend='gpu')
    print(
        f'Found local devices: {devices}, using device {_GPU_DEVICE.value}:'
        f' {devices[_GPU_DEVICE.value]}'
    )

    print('Building model from scratch...')
    model_runner = ModelRunner(
        config=make_model_config(
            flash_attention_implementation=typing.cast(
                attention.Implementation, _FLASH_ATTENTION_IMPLEMENTATION.value
            ),
            num_diffusion_samples=_NUM_DIFFUSION_SAMPLES.value,
            num_recycles=_NUM_RECYCLES.value,
            return_embeddings=_SAVE_EMBEDDINGS.value,
            return_distogram=_SAVE_DISTOGRAM.value,
        ),
        device=devices[_GPU_DEVICE.value],
        model_dir=pathlib.Path(MODEL_DIR.value),
    )
    # Check we can load the model parameters before launching anything.
    print('Checking that model parameters can be loaded...')
    _ = model_runner.model_params
  else:
    model_runner = None

  num_fold_inputs = 0
  for fold_input in fold_inputs:
    if _NUM_SEEDS.value is not None:
      print(f'Expanding fold job {fold_input.name} to {_NUM_SEEDS.value} seeds')
      fold_input = fold_input.with_multiple_seeds(_NUM_SEEDS.value)
    process_fold_input(
        fold_input=fold_input,
        data_pipeline_config=data_pipeline_config,
        model_runner=model_runner,
        output_dir=os.path.join(_OUTPUT_DIR.value, fold_input.sanitised_name()),
        buckets=tuple(int(bucket) for bucket in _BUCKETS.value),
        conformer_max_iterations=_CONFORMER_MAX_ITERATIONS.value,
    )
    num_fold_inputs += 1

  print(f'Done running {num_fold_inputs} fold jobs.')


if __name__ == '__main__':
  flags.mark_flags_as_required(['output_dir'])
  app.run(main)
