# Contributing to phenum

:+1::tada: Thank you for taking the time to contribute! :tada::+1:

The following are guidlines to help contribute to the phonon
enumeration python package (phenum). These are guidlines are open to
adaptation upon suggestion and at the contributors discretion.

#### Table of Contents

[Code of Conduct](https://github.com/wsmorgan/phonon-enumeration/blob/master/CODE_OF_CONDUCT.md)

[What should I know before I get started?](#what-should-i-know-before-i-get-started)

[How Can I Contribute](#how-can-i-contribute)
  * [Reporting Bugs](#reporting-bugs)
  * [Suggesting Enhancements](#suggesting-enhancements)
  * [Places to Start](#places-to-start)
  * [Testing](#testing)
  * [Pull Requests](#pull-requests)

[Styleguides](#styleguides)
  * [Git Commit Messages](#git-commit-messages)
  * [History Messages](#history-messages)
  * [Python Styleguide](#python-styleguide)
  * [Documentation Styleguide](#documentation-styleguide)

[Additional Notes](#additional-notes)
  * [Issue and Pull Request Labels](#issue-and-pull-request-labels)
  * [Attribution](#attribution)

## What should I know before I get startedd?

'Phenum' is a code distribution designed to list all the unique
arrangements of atoms on a lattice and the unique directions that
atoms on a lattice can be displaced. In order to contribute it would
be a good idea to familiarize yourself with some of the algorithms in
use. They can be found in the following papers:

[Initial enumeration paper](http://msg.byu.edu/papers/GLWHart_enumeration.pdf)
[Second enumeration paper](http://msg.byu.edu/papers/multi.pdf)
[Recursively stabilized enumeration paper](https://arxiv.org/abs/1701.02382)

## How can I Contribute

### Reporting Bugs

This section guides you through submitting a bug report for
'phenum'. Following these guidlines helps developers understand your
report and reproduce the issue.

> **Note:** If you find a **Closed** issue that seems like it is the
    same thing that you're experiencing, open a new issue and include
    a link to the original issue in the body of your new one.

#### How do I Submit A (Good) Bug Report?

Bugs are tracked as [GitHub
issues](https://guides.github.com/features/issues/). To report a bug
simply open a new issue following these guidlines.

Explain the problem and include additional details to help reproduce
the problem:

* **Use a clear and descriptive title** for the issue to identify the
    problem.
* **Include the input files (enum.in, lattice.in, and polya.*) as
    attatchements. If the files aren't available descibe, with as much
    detail as possible, the system being enumerated, i.e., the parent
    lattice type, atomic basis, number of atomic species,
    concentrations of atoms and arrows.....
* **Include the full error stack that was printed to the terminal in
    the bug report.

### Suggesting Enhancements

This section guides you through submitting an enhancement suggestion
for 'phenum', including completely new features and minor improvements to
existing functionality. 

#### How Do I Submit A (Good) Enhancement Suggestion?

Enhancement suggestions are tracked as [GitHub
issues](https://guides.github.com/features/issues/). To create a new
enhancement suggestion open an issue and provide the following
information:

* **Use a clear and descriptive title** for the issue to identify the
    suggestion.
* **Provide a step-by-step description of the suggested enhancement**
    in as many details as possible.
* **Describe the current behavior**, if applicable, and **explain
    which behavior you expected to see instead** and why.

* **Explain why this enhancement would be useful** to 'phenum' users.

### Your First Code Contribution

Unsure where to begin contributing to Atom? You can start by looking
through these `beginner` and `help-wanted` issues:

* [Beginner issues][beginner] - issues which should only require a few lines of code, and a test or two.
* [Help wanted issues][help-wanted] - issues which should be a bit more involved than `beginner` issues.

Both issue lists are sorted by total number of comments. While not
perfect, number of comments is a reasonable proxy for impact a given
change will have.

### Testing

Any new subroutine or functionality for 'phenum' must be unit tested
before a pull request will be accepted. New unit tests should:

* Be placed in the tests folder within a file titled test_'module
  name'.py where 'module name' is the module being tested.
* Any test input or output not in the test_*.py file should be placed
  in a folder with the same name as the module being tested.
* Unit tests must be designed to cover any new code written.

'Phenum' is unit tested being using tox in the root directory. Unit
tests on pull requests are conducted by TravisCI.

### Pull Requests

* Do not include issue numbers in the PR title
* Follow the [Python](#python-styleguide) styleguide.
* Document new code based on the [Documentation Styleguide](#documentation-styleguide).
* List relavent issue numbers in the PR body.
* PR's will only be accepted if they pass all unit tests, don't lower
  the codecoverage, quantified code finds no new issues in the PR, and
  the HISTORY.md and README.md have been updated to reflect changes.

## Styleguides

### Git Commit Messages

* Git commit messages should be short and reference the new version
  number, found in HISTORY.md, for this commit. All details of the
  commit changes should be recorded in HISTORY.md.

### History Messages

* Reference issues and pull requests liberally.
* Detail which modules were changed and how.

### Python Styleguide

'Phenum' follows the [Google Python Style](https://google.github.io/styleguide/pyguide.html).

### Documentation Styleguide

All new subroutines must be documented in order to keep the code maintainable and easy to use.

* Use [Google Style Python](http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html).

## Additional Notes

### Issue and Pull Request Labels

This section lists the labels used to help track and manage issues and pull requests.

[GitHub search](https://help.github.com/articles/searching-issues/)
makes it easy to use labels for finding groups of issues or pull
requests you're interested in. To help you find issues and pull
requests, each label is listed with search links for finding open
items with that label in `phenum` only and also across all Atom
repositories. We encourage you to read about [other search
filters](https://help.github.com/articles/searching-issues/) which
will help you write more focused queries.

Please open an issue if you have suggestions for new labels.

#### Type of Issue and Issue State

| Label name | `phenum` :mag_right: | Description |
| --- | --- | --- |
| `enhancement` | [search][search-label-enhancement]  | Feature requests. |
| `bug` | [search][search-label-bug] | Confirmed bugs or reports that are very likely to be bugs. |
| `question` | [search][search-label-question] | Questions more than bug reports or feature requests (e.g. how do I do X). |
| `feedback` | [search][search-label-feedback] | General feedback more than bug reports or feature requests. |
| `help-wanted` | [search][search-label-help-wanted] | The phenum development team would appreciate help from the community in resolving these issues. |
| `beginner` | [search][search-label-beginner] | Less complex issues which would be good first issues to work on for users who want to contribute to phenum. |
| `more-information-needed` | [search][search-label-more-information-needed] | More information needs to be collected about these problems or feature requests (e.g. steps to reproduce). |
| `needs-reproduction` | [search][search-label-needs-reproduction] | Likely bugs, but haven't been reliably reproduced. |
| `duplicate` | [search][search-label-duplicate] | Issues which are duplicates of other issues, i.e. they have been reported before. |
| `wontfix` | [search][search-label-wontfix] | The phenum development team has decided not to fix these issues for now, either because they're working as intended or for some other reason. |
| `invalid` | [search][search-label-invalid] | Issues which aren't valid (e.g. user errors). |

### Attribution

This document was adapted from [Contribute].
    
The Code of Conduct in this document is adapted from the [Contributor Covenant][homepage], version 1.4,
available at [http://contributor-covenant.org/version/1/4][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/4/
[Contribute]: https://github.com/atom/atom/blob/master/CONTRIBUTING.md
    
[search-label-enhancement]: https://github.com/wsmorgan/phonon-enumeration/issues?q=is%3Aissue+is%3Aopen+label%3Aenhancement
[search-label-bug]: https://github.com/wsmorgan/phonon-enumeration/issues?utf8=%E2%9C%93&q=is%3Aissue%20is%3Aopen%20label%3Abug
[search-label-question]: https://github.com/wsmorgan/phonon-enumeration/issues?utf8=%E2%9C%93&q=is%3Aissue%20is%3Aopen%20label%3Aquestion
[search-label-feedback]: https://github.com/wsmorgan/phonon-enumeration/issues?utf8=%E2%9C%93&q=is%3Aissue%20is%3Aopen%20label%3Afeedback
[search-label-help-wanted]: https://github.com/wsmorgan/phonon-enumeration/issues?q=is%3Aissue+is%3Aopen+label%3A%22help+wanted%22
[search-label-beginner]: https://github.com/wsmorgan/phonon-enumeration/issues?utf8=%E2%9C%93&q=is%3Aissue%20is%3Aopen%20label%3A%22beginner%22%20
[search-label-more-information-needed]: https://github.com/wsmorgan/phonon-enumeration/issues?q=is%3Aopen+is%3Aissue+label%3A%22more+information+needed%22
[search-label-needs-reproduction]: https://github.com/wsmorgan/phonon-enumeration/issues?q=is%3Aopen+is%3Aissue+label%3A%22needs+reproduction%22
[search-label-duplicate]: https://github.com/wsmorgan/phonon-enumeration/issues?q=is%3Aopen+is%3Aissue+label%3Aduplicate
[search-label-wontfix]: https://github.com/wsmorgan/phonon-enumeration/issues?q=is%3Aopen+is%3Aissue+label%3Awontfix
[search-label-invalid]: https://github.com/wsmorgan/phonon-enumeration/issues?q=is%3Aopen+is%3Aissue+label%3Ainvalid

[beginner]: https://github.com/wsmorgan/phonon-enumeration/issues?utf8=%E2%9C%93&q=is%3Aissue%20is%3Aopen%20label%3A%22beginner%22%20
[help-wanted]: https://github.com/wsmorgan/phonon-enumeration/issues?q=is%3Aissue+is%3Aopen+label%3A%22help+wanted%22
