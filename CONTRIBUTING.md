# Contributing to Plant-Cluster-Finder.py
Thank you for your interest in contributing to this project! We look forward to working with you, and are happy that you are here to help.

To facilitate easy integration with existing changes, this document outlines what we expect for any contribution. Failing to follow these instructions
may delay your contribution, so please double check your work!

## Make a Bug Report
Many users encounter bugs they don't have the bandwidth to fix. A good bug report will help us fix the bug quickly and completely.

A bug report should have
1. a title clearly describing the problem
2. how the bug was encountered
3. what you expected
4. any details about your environment that might help us find the bug

A contributor may ask for additional detail, or test conditions to ensure the bug is completely fixed before closing the issue.

## Make a Feature Request
We **love** new features, sometimes too much. There isn't enough time in a day to add all the features we want, 
so here's how to write a feature request to give it the best chance of getting implemented. Most importantly, be polite.

It must have
1. limited and clearly defined scope
2. a need that it fills
3. clear description of what the feature will add

Obviously the more you are willing to contribute to a new feature, the more likely it will get implemented. Which brings us to:
## Adding to the Codebase
Whether you're fixing a bug or implementing a feature, these are the guidelines that need to be followed to get your PR approved.

Your change should serve a single purpose, be tested and documented before your submit. 
### Writing Your Code
For style, Rye fmt is used. The linter `ruff` can be helpful while you're coding, or `black` when you're done.
However, it can be easiest to just run `rye fmt` on your code before committing it. 

Docstrings are written using the [Google styleguide](https://github.com/google/styleguide/blob/gh-pages/pyguide.md#38-comments-and-docstrings) and
rendered via Sphinx automatically.

Type annotations are expected where possible. [Mypy](http://www.mypy-lang.org/) is used for type checking in testing, and adding type annotations
can prevent type-related bugs from creeping in.
### Testing Your Code
Changes without tests will not be accepted. The test should check that your change handles expeceted input, doesn't break existing code, and 
either handles exceptions or throws meaningful exceptions on failure.

If you haven't written a test before, it can be helpful to write the test before writing your change. That way you think through exactly what 
your change is going to do.
## Building for Contribution
To make things easy, [Rye](https://rye-up.com/guide/) handles package management, building and virtual environment management. Install 
Rye on your dev machine to get started.

Once you've forked and cloned the repository, run `rye sync` to download and install all necessary packages.
Test using `rye run pytest` and build using `rye build`.
## Submitting a change
Quick checklist before writing a PR:
- [ ] The branch I'm submitting only has one change
- [ ] My code is clearly documented
- [ ] My code follows the style guide enforced by `rye fmt`
- [ ] My change is tested
- [ ] All tests pass

Create a new PR, if there is a bug/feature issue be sure to mention it in the PR.
1. Describe the change in detail. This shouldn't be a rehash of your documentation, but a description of what your change *does*.
2. How was it tested? Are there any edge cases you may have missed?

A contributor will review your code and initiate a test workflow. If there are requested changes, make them in the branch you submitted for a PR 
and post a message. Once the contributor approves the code and all tests pass, your change will be integrated into the next version.
