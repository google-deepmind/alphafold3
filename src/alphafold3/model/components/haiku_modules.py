# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under the Apache License, Version 2.0
# (the "License"); you may not use this file except in compliance with the
# License. You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# To request access to the AlphaFold 3 model parameters, follow the process set
# out at https://github.com/google-deepmind/alphafold3. You may only use these
# if received directly from Google. Use is subject to terms of use available at
# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

"""Common Haiku modules."""

from collections.abc import Sequence
import contextlib
import numbers
from typing import TypeAlias

import haiku as hk
import jax
import jax.numpy as jnp
import numpy as np


PRECISION: TypeAlias = (
    None
    | str
    | jax.lax.Precision
    | tuple[str, str]
    | tuple[jax.lax.Precision, jax.lax.Precision]
)

# Useful for mocking in tests.
DEFAULT_PRECISION = None

# Constant from scipy.stats.truncnorm.std(a=-2, b=2, loc=0., scale=1.)
TRUNCATED_NORMAL_STDDEV_FACTOR = np.asarray(
    0.87962566103423978, dtype=np.float32
)


class LayerNorm(hk.LayerNorm):
  """LayerNorm module.

  Equivalent to hk.LayerNorm but with an extra 'upcast' option that casts
  (b)float16 inputs to float32 before computing the layer norm, and then casts
  the output back to the input type.

  The learnable parameter shapes are also different from Haiku: they are always
  vectors rather than possibly higher-rank tensors. This makes it easier
  to change the layout whilst keep the model weight-compatible.
  """

  def __init__(
      self,
      *,
      axis: int = -1,
      create_scale: bool = True,
      create_offset: bool = True,
      eps: float = 1e-5,
      scale_init: hk.initializers.Initializer | None = None,
      offset_init: hk.initializers.Initializer | None = None,
      use_fast_variance: bool = True,
      name: str,
      param_axis: int | None = None,
      upcast: bool = True,
  ):
    super().__init__(
        axis=axis,
        create_scale=False,
        create_offset=False,
        eps=eps,
        scale_init=None,
        offset_init=None,
        use_fast_variance=use_fast_variance,
        name=name,
        param_axis=param_axis,
    )
    self.upcast = upcast
    self._temp_create_scale = create_scale
    self._temp_create_offset = create_offset

  def __call__(self, x: jnp.ndarray) -> jnp.ndarray:
    dtype = x.dtype
    is_16bit = x.dtype in [jnp.bfloat16, jnp.float16]
    if self.upcast and is_16bit:
      x = x.astype(jnp.float32)

    param_axis = self.param_axis[0] if self.param_axis else -1
    param_shape = (x.shape[param_axis],)

    param_broadcast_shape = [1] * x.ndim
    param_broadcast_shape[param_axis] = x.shape[param_axis]
    scale = None
    offset = None
    if self._temp_create_scale:
      scale = hk.get_parameter(
          'scale', param_shape, x.dtype, init=self.scale_init
      )
      scale = scale.reshape(param_broadcast_shape)

    if self._temp_create_offset:
      offset = hk.get_parameter(
          'offset', param_shape, x.dtype, init=self.offset_init
      )
      offset = offset.reshape(param_broadcast_shape)

    out = super().__call__(x, scale=scale, offset=offset)

    if self.upcast and is_16bit:
      out = out.astype(dtype)

    return out


def haiku_linear_get_params(
    inputs: jax.Array | jax.ShapeDtypeStruct,
    *,
    num_output: int | Sequence[int],
    use_bias: bool = False,
    num_input_dims: int = 1,
    initializer: str = 'linear',
    bias_init: float = 0.0,
    transpose_weights: bool = False,
    name: str | None = None,
) -> tuple[jax.Array, jax.Array | None]:
  """Get parameters for linear layer.

  Parameters will be at least float32 or higher precision.

  Arguments:
    inputs: The input to the Linear layer. Can be either a JAX array or a
      jax.ShapeDtypeStruct.
    num_output: The number of output channels. Can be an integer or a sequence
      of integers.
    use_bias: Whether to create a bias array.
    num_input_dims: The number of dimensions to consider as channel dims in the
      input.
    initializer: The name of the weight initializer to use.
    bias_init: A float used to initialize the bias.
    transpose_weights: If True, will create a transposed version of the weights.
    name: The Haiku namespace to use for the weight and bias.

  Returns:
    A tuple[weight, bias] if use_bias otherwise tuple[weight, None].
  """

  if isinstance(num_output, numbers.Integral):
    output_shape = (num_output,)
  else:
    output_shape = tuple(num_output)

  if num_input_dims > 0:
    in_shape = inputs.shape[-num_input_dims:]
  elif num_input_dims == 0:
    in_shape = ()
  else:
    raise ValueError('num_input_dims must be >= 0.')

  weight_init = _get_initializer_scale(initializer, in_shape)
  with hk.name_scope(name) if name else contextlib.nullcontext():

    if transpose_weights:
      weight_shape = output_shape + in_shape

      weights = hk.get_parameter(
          'weights', shape=weight_shape, dtype=inputs.dtype, init=weight_init
      )
    else:
      weight_shape = in_shape + output_shape
      weights = hk.get_parameter(
          name='weights',
          shape=weight_shape,
          dtype=inputs.dtype,
          init=weight_init,
      )

    bias = None
    if use_bias:
      bias = hk.get_parameter(
          name='bias',
          shape=output_shape,
          dtype=inputs.dtype,
          init=hk.initializers.Constant(bias_init),
      )
  return weights, bias


class Linear(hk.Module):
  """Custom Linear Module.

  This differs from the standard Linear in a few ways:
    * It supports inputs of arbitrary rank
    * It allows to use ntk parametrization
    * Initializers are specified by strings
    * It allows to explicitly specify which dimension of the input will map to
      the tpu sublane/lane dimensions.
  """

  def __init__(
      self,
      num_output: int | Sequence[int],
      *,
      initializer: str = 'linear',
      num_input_dims: int = 1,
      use_bias: bool = False,
      bias_init: float = 0.0,
      precision: PRECISION = None,
      fast_scalar_mode: bool = True,
      transpose_weights: bool = False,
      name: str,
  ):
    """Constructs Linear Module.

    Args:
      num_output: number of output channels. Can be tuple when outputting
        multiple dimensions.
      initializer: What initializer to use, should be one of {'linear', 'relu',
        'zeros'}.
      num_input_dims: Number of dimensions from the end to project.
      use_bias: Whether to include trainable bias (False by default).
      bias_init: Value used to initialize bias.
      precision: What precision to use for matrix multiplication, defaults to
        None.
      fast_scalar_mode: Whether to use optimized path for num_input_dims = 0.
      transpose_weights: decides whether weights have shape [input, output] or
        [output, input], True means [output, input], this is helpful to avoid
        padding on the tensors holding the weights.
      name: name of module, used for name scopes.
    """
    super().__init__(name=name)
    if isinstance(num_output, numbers.Integral):
      self.output_shape = (num_output,)
    else:
      self.output_shape = tuple(num_output)
    self.initializer = initializer
    self.use_bias = use_bias
    self.bias_init = bias_init
    self.num_input_dims = num_input_dims
    self.num_output_dims = len(self.output_shape)
    self.precision = precision if precision is not None else DEFAULT_PRECISION
    self.fast_scalar_mode = fast_scalar_mode
    self.transpose_weights = transpose_weights

  def __call__(self, inputs: jnp.ndarray) -> jnp.ndarray:
    """Connects Module.

    Args:
      inputs: Tensor of shape [..., num_channel]

    Returns:
      output of shape [..., num_output]
    """

    num_input_dims = self.num_input_dims

    # Adds specialized path for scalar inputs in Linear layer,
    # this means the linear Layer does not use the matmul units on the tpu,
    # which is more efficient and gives compiler more flexibility over layout.
    if num_input_dims == 0 and self.fast_scalar_mode:
      weight_shape = self.output_shape
      if self.initializer == 'zeros':
        w_init = hk.initializers.Constant(0.0)
      else:
        distribution_stddev = jnp.array(1 / TRUNCATED_NORMAL_STDDEV_FACTOR)
        w_init = hk.initializers.TruncatedNormal(
            mean=0.0, stddev=distribution_stddev
        )

      weights = hk.get_parameter('weights', weight_shape, inputs.dtype, w_init)

      inputs = jnp.expand_dims(
          inputs, tuple(range(-1, -self.num_output_dims - 1, -1))
      )
      output = inputs * weights
    else:
      if self.num_input_dims > 0:
        in_shape = inputs.shape[-self.num_input_dims :]
      else:
        in_shape = ()

      weight_init = _get_initializer_scale(self.initializer, in_shape)

      in_letters = 'abcde'[: self.num_input_dims]
      out_letters = 'hijkl'[: self.num_output_dims]

      if self.transpose_weights:
        weight_shape = self.output_shape + in_shape
        weights = hk.get_parameter(
            'weights', weight_shape, inputs.dtype, weight_init
        )
        equation = (
            f'...{in_letters}, {out_letters}{in_letters}->...{out_letters}'
        )
      else:
        weight_shape = in_shape + self.output_shape
        weights = hk.get_parameter(
            'weights', weight_shape, inputs.dtype, weight_init
        )

        equation = (
            f'...{in_letters}, {in_letters}{out_letters}->...{out_letters}'
        )

      output = jnp.einsum(equation, inputs, weights, precision=self.precision)

    if self.use_bias:
      bias = hk.get_parameter(
          'bias',
          self.output_shape,
          inputs.dtype,
          hk.initializers.Constant(self.bias_init),
      )
      output += bias

    return output


def _get_initializer_scale(initializer_name, input_shape):
  """Get initializer for weights."""

  if initializer_name == 'zeros':
    w_init = hk.initializers.Constant(0.0)
  else:
    # fan-in scaling
    noise_scale = 1.0
    for channel_dim in input_shape:
      noise_scale /= channel_dim
    if initializer_name == 'relu':
      noise_scale *= 2

    stddev = np.sqrt(noise_scale)
    # Adjust stddev for truncation.
    stddev = stddev / TRUNCATED_NORMAL_STDDEV_FACTOR
    w_init = hk.initializers.TruncatedNormal(mean=0.0, stddev=stddev)

  return w_init
