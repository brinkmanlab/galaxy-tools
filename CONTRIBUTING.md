Contributing
============

Contributions are welcome.

Submit contributions as a pull request on GitHub.

Commits to the master branch of this repository are automatically deployed to the [Galaxy Main Toolshed](https://toolshed.g2.bx.psu.edu/) using the included 
[GitHub Actions recipe](https://github.com/brinkmanlab/galaxy-tools/blob/master/.github/workflows/main-shed.yml)

Use .shed.template.yml as a template when submitting a new tool to the Galaxy main toolshed.

Standards
---------
- All contributions must follow IUC Standards: https://galaxy-iuc-standards.readthedocs.io/en/latest/index.html
- Tools must represent a free standing tool. If you have written a tool that is only applicable to a specific workflow then you need to reconsider your approach.
- All tools must be added to tool_conf.xml. If it is not to be included in the install, comment out the `<tool file="" />` as `<!-- tool file="" / -->`.
