# AlphaFold 3 Security Audit Report

This report documents potential vulnerabilities, security data leaks, bugs, and other findings from an automated and manual security review of the AlphaFold 3 repository.

## 1. Executive Summary
- **Critical Vulnerabilities:** None found.
- **Medium/Low Vulnerabilities:** Supply chain risks in database downloads (missing checksums), unsafe deserialization in test suites.
- **Data Leaks:** Automated scans flagged several potential secrets, but they were confirmed to be false positive biological sequence data.
- **Bugs and Code Quality:** Multiple `flake8` warnings and assertions used in production code paths.

---

## 2. Potential Vulnerabilities & Security Findings

### 2.1 Insecure Deserialization (Unsafe `pickle` usage in Tests)
- **Severity:** Medium
- **Location:**
  - `run_alphafold_test.py` (lines 158, 366, 361)
- **Description:**
  The test suite uses standard Python `pickle.loads()` to deserialize test data from the `test_data/` directory. While AlphaFold includes a `safe_pickle.py` module to mitigate arbitrary code execution in production data loading (`constants/chemical_components.py`), the test files still rely on the standard `pickle` library.
  If an attacker modifies the `pkl` files in `test_data/`, they could achieve arbitrary code execution on the machine running the tests.
- **Recommendation:**
  Update the test suites to utilize `alphafold3.common.safe_pickle.load()` or avoid `pickle` entirely in favor of `json` or other safe serialization formats.

### 2.2 Supply Chain Risk: Missing Checksum Verification in Database Downloads
- **Severity:** Low/Medium
- **Location:**
  - `fetch_databases.sh`
- **Description:**
  The `fetch_databases.sh` script downloads large datasets via `wget` and unpacks them using `tar` and `zstd`. Although the script downloads from an HTTPS endpoint (`https://storage.googleapis.com/alphafold-databases/v3.0`), it lacks cryptographic checksum (e.g., SHA256) verification of the downloaded files before untarring.
  If the upstream bucket is compromised or the data is tampered with in transit (e.g., via MITM if HTTPS fails or is bypassed), an attacker could provide maliciously crafted `tar` or `fasta` files that could exploit vulnerabilities in the decompression tools or poison the data pipeline.
- **Recommendation:**
  Provide a signed manifest or a list of known good SHA256 checksums alongside the script, and compute the checksum of the downloaded `.tar.zst` and `.zst` files before unzipping them.

### 2.3 Potential Command Injection Risks in External Tool Wrapping
- **Severity:** Low
- **Location:**
  - `src/alphafold3/data/tools/hmmsearch.py`
  - `src/alphafold3/data/tools/jackhmmer.py`
- **Description:**
  AlphaFold interacts with multiple external bioinformatics binaries (e.g., `hmmsearch`, `jackhmmer`). It utilizes `subprocess_utils.run()` which passes commands as a list (e.g., `[binary_path, flag1, flag2]`). While passing arguments as lists prevents simple shell injection (because `shell=False` is the default in `subprocess.run`), any user-controlled inputs that dictate the `binary_path` or arbitrary flags could still lead to unintended command execution.
- **Recommendation:**
  Ensure that any paths or external configurations passed to these wrappers are strictly sanitized. Given the tool's intended use case (research environment), this is a low risk, but should be noted.

---

## 3. Data Leaks & Hardcoded Secrets
- **Tool Used:** `detect-secrets`
- **Findings:** The tool identified 24 instances of potential "AWS Access Keys" and "Artifactory Credentials" in the miniature database fasta files (`uniprot_all__subsampled_1000.fasta` and `uniref90__subsampled_1000.fasta`).
- **Conclusion:** **False Positives.** Manual review of the lines flagged (e.g., `cd3412c795527ebdd001833d1f157e9941862709`) shows that these are purely biological sequence representations or hashed IDs used in standard FASTA format, not actual secrets.
- **Recommendation:** Configure a `.secrets.baseline` file to ignore these specific test data directories.

---

## 4. Bugs and Code Quality Issues

### 4.1 Improper Use of `assert` in Production Code
- **Severity:** Low
- **Location:**
  - 88 instances found across the codebase (e.g., `src/alphafold3/model/network/template_modules.py`, `src/alphafold3/model/scoring/alignment.py`, `src/alphafold3/structure/structure.py`).
- **Description:**
  The codebase uses the `assert` statement heavily to check input shapes, dimensions, and logical invariants in the neural network and data pipeline modules. Python's `assert` statements are removed when the code is compiled to optimized byte code (using `python -O`). This means any critical validations relying on `assert` will be skipped in an optimized environment, potentially leading to hard-to-debug logic errors or out-of-bounds array access.
- **Recommendation:**
  Replace critical `assert` statements with `if not condition: raise ValueError(...)` to ensure the checks run regardless of optimization flags.

### 4.2 Flake8 Findings
- **Severity:** Low
- **Location:** Across the repository.
- **Description:**
  Running `flake8` yielded over 6,000 warnings. The vast majority correspond to minor style issues (e.g., line lengths, spacing) that do not impact security. However, highly complex or messy code can obscure logical flaws over time.
- **Recommendation:**
  Adopt an automated linter/formatter (like `black` or `ruff`) and enforce it in a CI/CD pipeline.
