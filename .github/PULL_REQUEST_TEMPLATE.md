# Description
<!--- Briefly describe the changes included in this pull request and the paths to the test cases below
 !--- starting with 'Closes #...' if appropriate --->

### Closes #...

## Testing Results
<!-- Delete unused cases below. Usually NFTest is desired and sufficient -->
- Case 1
    - sample:    <!-- e.g. A-mini S2.T-1, A-mini S2.T-n1 -->
    - input csv: <!-- path/to/input.csv -->
    - config:    <!-- path/to/xxx.config -->
    - output:    <!-- path/to/output -->
- Case 2
    - sample:    <!-- e.g. A-mini S2.T-1, A-mini S2.T-n1 -->
    - input csv: <!-- path/to/input.csv -->
    - config:    <!-- path/to/xxx.config -->
    - output:    <!-- path/to/output -->
- NFTest
    - output:    <!-- path/to/output -->
    - log:       <!-- path/to/log -->
    - cases:     <!-- list of test cases run -->

# Checklist
<!--- Please read each of the following items and confirm by replacing the [ ] with a [X] --->

- [ ] I have read the [code review guidelines](https://uclahs-cds.atlassian.net/wiki/spaces/BOUTROSLAB/pages/3187646/Code+Review+Guidelines) and the [code review best practice on GitHub check-list](https://uclahs-cds.atlassian.net/wiki/spaces/BOUTROSLAB/pages/3189956/Code+Review+Best+Practice+on+GitHub+-+Check+List).

- [ ] I have reviewed the [Nextflow pipeline standards](https://uclahs-cds.atlassian.net/wiki/spaces/BOUTROSLAB/pages/3193890/Nextflow+pipeline+standardization).

- [ ] The name of the branch is meaningful and well formatted following the [standards](https://uclahs-cds.atlassian.net/wiki/spaces/BOUTROSLAB/pages/3189956/Code+Review+Best+Practice+on+GitHub+-+Check+List), using \[AD_username (or 5 letters of AD if AD is too long)]-\[brief_description_of_branch].

- [ ] I have set up or verified the branch protection rule following the [github standards](https://uclahs-cds.atlassian.net/wiki/spaces/BOUTROSLAB/pages/3190380/GitHub+Standards#GitHubStandards-Branchprotectionrule) before opening this pull request.

- [ ] I have added my name to the contributors listings in the ``manifest`` block in the `nextflow.config` as part of this pull request, am listed already, or do not wish to be listed. (*This acknowledgement is optional.*)

- [ ] I have added the changes included in this pull request to the `CHANGELOG.md` under the next release version or unreleased, and updated the date.

- [ ] I have updated the version number in the `metadata.yaml` and `manifest` block of the `nextflow.config` file following [semver](https://semver.org/), or the version number has already been updated. (*Leave it unchecked if you are unsure about new version number and discuss it with the infrastructure team in this PR.*)

- [ ] I have tested the pipeline on at least one A-mini sample.
