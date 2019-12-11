.. image:: https://zenodo.org/badge/162751955.svg
   :target: https://zenodo.org/badge/latestdoi/162751955

============
galaxy-tools
============

Central repo for all Brinkman Lab maintained Galaxy Tools

**Do not edit README.rst**, edit README.header.rst and regenerate README.rst using ./generate_readme

See CONTRIBUTING.rst_ for information on contributing to this repo.

.. _CONTRIBUTING.rst: CONTRIBUTING.rst

Manual installation
-------------------
- Clone this repo into the Galaxy server
- Add the full `./tool_conf.xml` path to the Galaxy config `tool_config_file:` list
- Use Galaxies `Manage dependencies` admin panel to install all tool dependencies

Tools
-----
See tool xml files for more information

=======================================  =========================================================================================================================  =========================  ====
Name                                     Desc                                                                                                                       ID                         Path
=======================================  =========================================================================================================================  =========================  ====
SigiHMM                                  Score-based prediction of genomic islands in prokaryotic genomes using hidden Markov models                                sigihmm                    colombo/SigiHMM.xml
GFF/GTF Feature Merge                    Merge GFF features based on a variety of criteria                                                                          feature-merge              feature_merge/feature_merge.xml
IslandPath-DIMOB                         predict genomic islands in bacterial and archaeal genomes based on the presence of dinucleotide biases and mobility genes  islandpath-dimob           islandpath-dimob/IslandPath-DIMOB.xml
MASH                                     Fast genome and metagenome distance estimation using MinHash                                                               mash                       mash/mash.xml
ParSNP                                   Efficient microbial core genome alignment and SNP detection                                                                parsnp                     parsnp/ParSNP.xml
MCL                                      The Markov Cluster Algorithm, a cluster algorithm for graphs                                                               mcl                        mcl/mcl.xml
BioPython SeqIO Converter                Interconvert between the various sequence file formats that BioPython supports                                             biopython-convert          biopython-convert/biopython-convert.xml
BioPython Phylo Parse Newick Leaf Order  Convert a newick tree to an ordered list of its leaves                                                                     extract-tree-order         extract_tree_order/extract_tree_order.xml
Mauve Contig Mover                       Reorder a multi-contig dataset against a reference genome                                                                  mauve-contig-mover         mauve_contig_mover/mauve_contig_mover.xml
Mauve Contig Mover - Stitch              Concatenate multiple contigs, complementing reversed sequences and rewriting all feature coordinates                       mauve-contig-mover-stitch  mauve_contig_mover/mcm_stitch.xml
BioPython Make Unique ID                 Makes all record ids unique across all input data                                                                          make-unique-id             make_unique_id/make_unique_id.xml
Send Email                               Send mail using templates and include information about datasets                                                           sendmail                   sendmail/sendmail.xml
Coreutils sha256sum                      Generate or check SHA256 (256-bit) checksums                                                                               sha256sum                  hash/sha256sum.xml
AWK Script                               Transform, modify, or generate data                                                                                        awkscript                  awkscript/awkscript.xml
Inspect                                  Dump Cheetah tool template environment. Useful for writing advanced tools. DO NOT INSTALL ON A PUBLIC SERVER               inspect                    inspect/inspect.xml
=======================================  =========================================================================================================================  =========================  ====
