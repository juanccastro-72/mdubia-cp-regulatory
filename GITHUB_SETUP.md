# Publishing this repository to GitHub

The repository is named **mdubia-cp-regulatory** under your account **juanccastro-72**.

1. Create an empty repo on GitHub (no README/licence, we provide them):
   https://github.com/new  →  owner `juanccastro-72`, name `mdubia-cp-regulatory`,
   Private for now (make Public on acceptance).

2. From this folder:

   ```bash
   cd mdubia-cp-regulatory
   git init -b main
   git add .
   git commit -m "Initial release: M. dubia chloroplast regulatory-element pipeline (v2.0.0)"
   git remote add origin https://github.com/juanccastro-72/mdubia-cp-regulatory.git
   git push -u origin main
   ```

3. Tag the release and (on acceptance) archive on Zenodo:

   ```bash
   git tag -a v2.0.0 -m "v2.0.0"
   git push origin v2.0.0
   ```
   Then enable the GitHub–Zenodo integration (https://zenodo.org/account/settings/github/)
   and create the release to mint a DOI; add that DOI to `CITATION.cff` and the paper.

Large regenerated outputs (genomes/, results/, figures/) are git-ignored; only the
source, data inputs, environment and tests are tracked.
