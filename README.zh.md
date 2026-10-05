<p align="center">
  <a href="README.md">English</a> · <b>简体中文</b>
</p>

![header](docs/header.jpg)

# AlphaFold 3

本软件包提供了 AlphaFold 3 推理管线（Inference Pipeline）的实现。有关如何获取模型参数，请参见下文。您只有直接从 Google 处获取才可以使用 AlphaFold 3 模型参数。使用须遵守这些[使用条款](https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md)。

任何公开发表披露因使用本源代码、模型参数或由其生成的输出结果的研究成果，均应[引用](#引用本研究-citing-this-work)论文：
[Accurate structure prediction of biomolecular interactions with AlphaFold 3](https://doi.org/10.1038/s41586-024-07487-w)（基于 AlphaFold 3 的生物分子相互作用精确结构预测）。

有关该方法的详细技术说明，请参阅补充信息（Supplementary Information）。

AlphaFold 3 亦可通过 [alphafoldserver.com](https://alphafoldserver.com) 免费用于非商业用途，但支持的配体（ligands）和共价修饰（covalent modifications）集合相对更为受限。

如需用于商业用途，AlphaFold 3 可通过 [Google Cloud 上的 Gemini Enterprise Agent Platform](https://docs.cloud.google.com/gemini-enterprise-agent-platform/models/open-models/alphafold-3) 获取。

如有任何疑问，请联系 AlphaFold 团队：[alphafold@google.com](mailto:alphafold@google.com)。

## 获取模型参数 (Obtaining Model Parameters)

本仓库包含了运行 AlphaFold 3 推理所需的全部代码。您可以从以下地址下载 AlphaFold 3 模型参数：
https://storage.googleapis.com/alphafold3/af3.bin.zst 。使用须遵守这些[使用条款](https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md)。

用于预测带有数字水印结构的 SynthID Bio-structure 版本 AlphaFold 3 模型参数可通过以下地址获取：
https://storage.googleapis.com/alphafold3/af3_synthid.bin.zst 。详情请参阅 [google-deepmind/synthidbio](https://github.com/google-deepmind/synthidbio)。

## 安装与运行首次预测 (Installation and Running Your First Prediction)

请参阅[安装文档](docs/installation.md)。

完成 AlphaFold 3 的安装后，您可以使用如下名为 `fold_input.json` 的输入 JSON 文件测试您的环境：

```json
{
  "name": "2PV7",
  "sequences": [
    {
      "protein": {
        "id": ["A", "B"],
        "sequence": "GMRESYANENQFGFKTINSDIHKIVIVGGYGKLGGLFARYLRASGYPISILDREDWAVAESILANADVVIVSVPINLTLETIERLKPYLTENMLLADLTSVKREPLAKMLEVHTGAVLGLHPMFGADIASMAKQVVVRCDGRFPERYEWLLEQIQIWGAKIYQTNATEHDHNMTYIQALRHFSTFANGLHLSKQPINLANLLALSSPIYRLELAMIGRLFAQDAELYADIIMDKSENLAVIETLKQTYDEALTFFENNDRQGFIDAFHKVRDWFGDYSEQFLKESRQLLQQANDLKQG"
      }
    }
  ],
  "modelSeeds": [1],
  "dialect": "alphafold3",
  "version": 1
}
```

随后，您可以使用以下命令运行 AlphaFold 3：

```bash
docker run -it \
    --volume $HOME/af_input:/root/af_input \
    --volume $HOME/af_output:/root/af_output \
    --volume <MODEL_PARAMETERS_DIR>:/root/models \
    --volume <DATABASES_DIR>:/root/public_databases \
    --gpus all \
    alphafold3 \
    python run_alphafold.py \
    --json_path=/root/af_input/fold_input.json \
    --model_dir=/root/models \
    --output_dir=/root/af_output
```

您可以向 `run_alphafold.py` 命令传递多种参数标志，若要列出全部标志请运行 `python run_alphafold.py --help`。其中控制 AlphaFold 3 运行哪些组件的两个核心参数标志为：

*   `--run_data_pipeline`（默认为 `true`）：是否运行数据预处理管线（即遗传序列与模板搜索）。该部分仅使用 CPU，耗时较长，可以在没有 GPU 的机器上单独运行。
*   `--run_inference`（默认为 `true`）：是否运行模型推理。该部分需要 GPU 支持。

## AlphaFold 3 输入数据

请参阅[输入文档](docs/input.md)。

## AlphaFold 3 输出结果

请参阅[输出文档](docs/output.md)。

## 性能基准 (Performance)

请参阅[性能文档](docs/performance.md)。

## 已知问题 (Known Issues)

已知问题记录在[已知问题文档](docs/known_issues.md)中。

若遇到未在[已知问题文档](docs/known_issues.md)或 [Issue 追踪器](https://github.com/google-deepmind/alphafold3/issues)中列出的问题，请[提交新的 Issue](https://github.com/google-deepmind/alphafold3/issues/new/choose)。

<a id="citing-this-work"></a>
## 引用本研究 (Citing This Work)

任何公开发表披露因使用本源代码、模型参数或由其生成的输出结果的研究成果，均应引用：

```bibtex
@article{Abramson2024,
  author  = {Abramson, Josh and Adler, Jonas and Dunger, Jack and Evans, Richard and Green, Tim and Pritzel, Alexander and Ronneberger, Olaf and Willmore, Lindsay and Ballard, Andrew J. and Bambrick, Joshua and Bodenstein, Sebastian W. and Evans, David A. and Hung, Chia-Chun and O’Neill, Michael and Reiman, David and Tunyasuvunakool, Kathryn and Wu, Zachary and Žemgulytė, Akvilė and Arvaniti, Eirini and Beattie, Charles and Bertolli, Ottavia and Bridgland, Alex and Cherepanov, Alexey and Congreve, Miles and Cowen-Rivers, Alexander I. and Cowie, Andrew and Figurnov, Michael and Fuchs, Fabian B. and Gladman, Hannah and Jain, Rishub and Khan, Yousuf A. and Low, Caroline M. R. and Perlin, Kuba and Potapenko, Anna and Savy, Pascal and Singh, Sukhdeep and Stecula, Adrian and Thillaisundaram, Ashok and Tong, Catherine and Yakneen, Sergei and Zhong, Ellen D. and Zielinski, Michal and Žídek, Augustin and Bapst, Victor and Kohli, Pushmeet and Jaderberg, Max and Hassabis, Demis and Jumper, John M.},
  journal = {Nature},
  title   = {Accurate structure prediction of biomolecular interactions with AlphaFold 3},
  year    = {2024},
  volume  = {630},
  number  = {8016},
  pages   = {493–-500},
  doi     = {10.1038/s41586-024-07487-w}
}
```

<a id="acknowledgements"></a>
## 致谢 (Acknowledgements)

AlphaFold 3 的发布离不开以下各位成员的宝贵贡献：

Andrew Cowie, Bella Hansen, Charlie Beattie, Chris Jones, Grace Margand,
Jacob Kelly, James Spencer, Josh Abramson, Kathryn Tunyasuvunakool, Kuba Perlin,
Lindsay Willmore, Max Bileschi, Molly Beck, Oleg Kovalevskiy,
Sebastian Bodenstein, Sukhdeep Singh, Tim Green, Toby Sargeant, Uchechi Okereke,
Yotam Doron, 以及 Augustin Žídek（工程主管 / engineering lead）。

我们同样向在 Google 和 Isomorphic Labs 的合作伙伴致以衷心感谢。

AlphaFold 3 使用了以下独立的开源库和软件包：

*   [abseil-cpp](https://github.com/abseil/abseil-cpp) 和
    [abseil-py](https://github.com/abseil/abseil-py)
*   [Docker](https://www.docker.com)
*   [DSSP](https://github.com/PDB-REDO/dssp)
*   [HMMER Suite](https://github.com/EddyRivasLab/hmmer)
*   [Haiku](https://github.com/deepmind/dm-haiku)
*   [JAX](https://github.com/jax-ml/jax/)
*   [libcifpp](https://github.com/pdb-redo/libcifpp)
*   [NumPy](https://github.com/numpy/numpy)
*   [pybind11](https://github.com/pybind/pybind11) 和
    [pybind11_abseil](https://github.com/pybind/pybind11_abseil)
*   [RDKit](https://github.com/rdkit/rdkit)
*   [Tokamax](https://github.com/openxla/tokamax)
*   [tqdm](https://github.com/tqdm/tqdm)

我们感谢所有开源贡献者和维护者！

## 联系我们 (Get in Touch)

如果您有本概述中未涵盖的任何疑问，请通过 [alphafold@google.com](mailto:alphafold@google.com) 联系 AlphaFold 团队。

我们非常期待听到您的反馈，了解 AlphaFold 3 在您的科研工作中发挥的作用。欢迎通过 [alphafold@google.com](mailto:alphafold@google.com) 与我们分享您的故事。

有关 SynthID Bio-structure 的相关疑问，请联系 [synthidbio@google.com](mailto:synthidbio@google.com)。

## 许可证与免责声明 (Licence and Disclaimer)

这不是 Google 官方支持的产品。

版权所有 2024 DeepMind Technologies Limited。

### AlphaFold 3 源代码与模型参数

AlphaFold 3 源代码采用 Apache License 2.0 许可证授权（“许可证”）；除非遵守许可证规定，否则您不得使用其源代码。您可在以下网址获取许可证副本：
http://www.apache.org/licenses/LICENSE-2.0

AlphaFold 3 模型参数根据 [AlphaFold 3 模型参数使用条款](https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md)（“条款”）提供；除非遵守这些条款，否则您不得使用它们。您可在以下网址获取条款副本：
[https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md](https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md)。

除非适用法律另有要求，AlphaFold 3 及其输出均按“原样”（AS IS）基础分发，不附带任何明示或暗示的保证或条件。您对确定使用 AlphaFold 3 或使用或分发其源代码或输出的适当性承担全部责任，并承担与此类使用或分发以及您行使相关条款下的权利和义务相关的任何及所有风险。输出结果是具有不同置信度水平的预测值，应仔细加以解释。在依赖、发表、下载或以其他方式使用 AlphaFold 3 资产前，请审慎判断。

AlphaFold 3 及其输出仅用于理论建模。它们未经临床使用验证或批准，不适用于临床用途。您不应将 AlphaFold 3 或其输出用于临床目的，或依赖它们获取医疗或其他专业建议。与这些主题相关的任何内容仅供参考，不能替代合格专业人员的建议。有关条款下许可与限制的具体语言表述，请参阅相关条款。

### 第三方软件 (Third-party Software)

使用上述[致谢](#致谢-acknowledgements)部分中提到的第三方软件、库或代码可能受单独的条款和条件或许可证规定的约束。您对第三方软件、库或代码的使用须遵守任何此类条款，在使用前应检查确认自己能够遵守任何适用的限制或条款条件。

### 镜像与参考数据库 (Mirrored and Reference Databases)

以下数据库已经由：(1) Google DeepMind 进行了镜像；以及 (2) 部分包含在推理代码包中用于测试目的，可参照以下信息使用：

*   [BFD](https://bfd.mmseqs.com/)（已修改），由 Steinegger M. 和 Söding J. 提供，经 Google DeepMind 修改，在 [Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/deed.en) 许可证下可用。详情请参阅 [AlphaFold 蛋白质组论文](https://www.nature.com/articles/s41586-021-03828-1) 的方法（Methods）部分。
*   [PDB](https://wwpdb.org)（未修改），由 H.M. Berman 等人提供，免除所有版权限制，并在 [CC0 1.0 Universal (CC0 1.0) Public Domain Dedication](https://creativecommons.org/publicdomain/zero/1.0/) 下完全免费供非商业和商业用途使用。
*   [MGnify: v2022\_05](https://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/2022_05/README.txt)（未修改），由 Mitchell AL 等人提供，免除所有版权限制，并在 [CC0 1.0 Universal (CC0 1.0) Public Domain Dedication](https://creativecommons.org/publicdomain/zero/1.0/) 下完全免费供非商业和商业用途使用。
*   [UniProt: 2021\_04](https://www.uniprot.org/)（未修改），由 The UniProt Consortium 提供，在 [Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/deed.en) 许可证下可用。
*   [UniRef90: 2022\_05](https://www.uniprot.org/)（未修改），由 The UniProt Consortium 提供，在 [Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/deed.en) 许可证下可用。
*   [NT: 2023\_02\_23](https://www.ncbi.nlm.nih.gov/nucleotide/)（已修改），详情请参阅 [AlphaFold 3 论文](https://nature.com/articles/s41586-024-07487-w) 的补充信息。
*   [RFam: 14\_4](https://rfam.org/)（已修改），由 I. Kalvari 等人提供，免除所有版权限制，并在 [CC0 1.0 Universal (CC0 1.0) Public Domain Dedication](https://creativecommons.org/publicdomain/zero/1.0/) 下完全免费供非商业和商业用途使用。详情请参阅 [AlphaFold 3 论文](https://nature.com/articles/s41586-024-07487-w) 的补充信息。
*   [RNACentral: 21\_0](https://rnacentral.org/)（已修改），由 The RNAcentral Consortium 提供，免除所有版权限制，并在 [CC0 1.0 Universal (CC0 1.0) Public Domain Dedication](https://creativecommons.org/publicdomain/zero/1.0/) 下完全免费供非商业和商业用途使用。详情请参阅 [AlphaFold 3 论文](https://nature.com/articles/s41586-024-07487-w) 的补充信息。

---

> 💡 **文档维护说明**：本中文文档由社区志愿者（[@JasonYeYuhe](https://github.com/JasonYeYuhe)）翻译维护，最后同步更新于 2026年10月05日。如发现内容与官方英文原版存在差异或新特性滞后，欢迎提交 PR 共同完善！
