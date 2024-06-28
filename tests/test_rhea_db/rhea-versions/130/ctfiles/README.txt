================================================================================
Rhea CT File
================================================================================

This directory contains Rhea data in chemical table file (CT File) formats, a
family of text formats that describe molecules and chemical reactions, see
https://en.wikipedia.org/wiki/Chemical%5Ftable%5Ffile for more information.

- rhea-rxn.tar.gz and rxn/

  Archive file and directory with one RXN file for each unidirectional Rhea
  reaction. These files are named after the Rhea IDs (i.e. the RXN file for
  RHEA:XXXXX is XXXXX.rxn).

- rhea-rd.tar.gz and rd/

  Archive file and directory with one RD file for each unidirectional Rhea
  reaction. These files are named after the Rhea IDs (i.e. the RD file for
  RHEA:XXXXX is XXXXX.rd).
  
  rhea-relXXX.rd.gz
  An RD file with all unidirectional Rhea reactions of release XXX.

- rhea-mol.tar.gz and mol/

  Archive file and directory with one MOL file for each ChEBI compound used to
  describe Rhea reaction participants. These files are named after the ChEBI IDs
  (i.e. the MOL file for CHEBI:XXXXX is XXXXX.mol).

- rhea.sdf.gz
  An SDF file with all ChEBI compounds used to describe Rhea reaction
  participants.

Please avoid browsing the subdirectories because they contain a large number of files.
